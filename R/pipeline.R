# ==============================================================================
# pipeline.R  --  End-to-end haplotype block pipeline
#
# Steps:
#   1. Open genotype backend
#   2. Streaming MAF pre-screen
#   3. Chromosome subset
#   4. Load genotype matrix
#   4.5 Call-rate -> MAF (authoritative) -> Imputation (all C++)
#   4.6 Build post-imputation backend
#   5. Genome-wide LD block detection
#   5b. Optional Beagle phasing (phase = TRUE only)
#   6. Haplotype extraction
#   7. Haplotype diversity
#   8. Haplotype feature matrix
#   9. Write outputs
# ==============================================================================

#' @noRd
.maf_filter_backend <- function(be, maf_cut = 0.05, verbose = TRUE) {
  if (verbose) message("[MAF filter] Computing MAF for ", be$n_snps, " SNPs ...")
  chrs <- unique(be$snp_info$CHR); keep_ids <- character(0)
  for (chr in chrs) {
    idx <- which(be$snp_info$CHR == chr); geno_chr <- read_chunk(be, idx)
    keep_ids <- c(keep_ids, be$snp_info$SNP[idx][maf_filter_cpp(geno_chr, maf_cut)])
    rm(geno_chr); gc(FALSE)
  }
  if (verbose) message("[MAF filter] ", length(keep_ids), " / ", be$n_snps, " SNPs pass MAF >= ", maf_cut)
  be$snp_info[be$snp_info$SNP %in% keep_ids, , drop = FALSE]
}

#' @noRd
.make_backed_backend_impl <- function(geno_mat, snp_info, sample_ids,
                                      stem_name, fingerprint = NULL,
                                      use_bigmemory = FALSE,
                                      bigmemory_path = tempdir(),
                                      bigmemory_type = "char",
                                      verbose = TRUE) {
  if (!is.matrix(geno_mat)) geno_mat <- as.matrix(geno_mat)
  rownames(geno_mat) <- sample_ids; colnames(geno_mat) <- snp_info$SNP
  if (!isTRUE(use_bigmemory) || !requireNamespace("bigmemory", quietly = TRUE))
    return(read_geno(path = geno_mat, format = "matrix",
                     snp_info = snp_info, sample_ids = sample_ids, verbose = FALSE))
  bigmemory_path <- gsub("\\\\", "/", normalizePath(bigmemory_path, mustWork = FALSE))
  stem <- file.path(bigmemory_path, stem_name)
  bin_file <- paste0(stem, ".bin"); desc_file <- paste0(stem, ".desc")
  si_file <- paste0(stem, "_snpinfo.rds"); sid_file <- paste0(stem, "_sampleids.rds")
  fp_file <- paste0(stem, "_params.rds"); bin_noext <- sub("[.]bin$", "", bin_file)
  bin_ok <- file.exists(bin_file) || file.exists(bin_noext)
  all_ok <- bin_ok && file.exists(desc_file) && file.exists(si_file) && file.exists(sid_file)
  if (all_ok && !is.null(fingerprint) && file.exists(fp_file)) {
    saved_fp <- tryCatch(readRDS(fp_file), error = function(e) NULL)
    if (identical(saved_fp, fingerprint)) {
      if (verbose) message("[", stem_name, "] Reattaching existing backend.")
      return(read_geno_bigmemory(source = desc_file, snp_info = readRDS(si_file),
                                 backingfile = basename(stem), backingpath = bigmemory_path))
    }
  }
  for (f in c(bin_file, bin_noext, desc_file, si_file, sid_file, fp_file))
    if (file.exists(f)) file.remove(f)
  if (verbose) message("[", stem_name, "] Building file-backed backend.")
  be_out <- read_geno_bigmemory(source = geno_mat, snp_info = snp_info,
                                backingfile = basename(stem), backingpath = bigmemory_path,
                                type = bigmemory_type, verbose = verbose)
  saveRDS(snp_info, si_file); saveRDS(sample_ids, sid_file)
  if (!is.null(fingerprint)) saveRDS(fingerprint, fp_file)
  be_out
}

#' @noRd
.make_imputed_backend <- function(geno_mat, snp_info, sample_ids,
                                  use_bigmemory = FALSE, bigmemory_path = tempdir(),
                                  bigmemory_type = "char", maf_cut = 0.05,
                                  min_callrate = 0.0, impute_method = "none",
                                  verbose = TRUE) {
  .make_backed_backend_impl(geno_mat, snp_info, sample_ids,
                            stem_name = "ldxblocks_bm_imputed",
                            fingerprint = list(maf_cut=maf_cut, min_callrate=min_callrate,
                                               impute_method=impute_method, n_snps=nrow(snp_info)),
                            use_bigmemory=use_bigmemory, bigmemory_path=bigmemory_path,
                            bigmemory_type=bigmemory_type, verbose=verbose)
}


# -- Internal: write filtered/cleaned genotype data as VCF.gz for Beagle ------
#
# After MAF filtering, call-rate filtering, chr subsetting, and optional
# malformed-line removal, the in-memory geno_mat represents a different SNP
# set than geno_source on disk. Passing geno_source directly to Beagle would:
#   a) include malformed records that were removed by clean_malformed = TRUE
#   b) include SNPs that failed MAF or call-rate filters
#   c) include chromosomes excluded by the chr= argument
#   d) cause SNP-count mismatches in .read_and_cache_phased_vcf() alignment
#
# This function writes geno_mat (post-filtering, imputed 0/1/2) as a minimal
# valid VCF.gz that Beagle 5.x can read directly. The output contains
# exactly the SNPs in snp_info_filtered, sorted by CHR then POS.
#
# Dosage encoding:
#   0  -> 0/0  (homozygous REF)
#   1  -> 0/1  (heterozygous)
#   2  -> 1/1  (homozygous ALT)
#   NA -> ./.  (missing)
#
# Beagle 5.x accepts plain gzip-compressed VCF (no tabix index required for
# the input file). gzfile() produces standard gzip compatible with Beagle.
#
# Returns the path to the written .vcf.gz file.
#' @noRd
.write_cleaned_vcf <- function(geno_mat,
                               snp_info,
                               sample_ids,
                               out_file,
                               verbose = TRUE) {
  # snp_info must have: CHR, POS, SNP (used as ID), REF, ALT
  # REF/ALT fall back to "A"/"T" if absent (Beagle does not use them for phasing)
  has_ref <- "REF" %in% names(snp_info)
  has_alt <- "ALT" %in% names(snp_info)

  # Sort by CHR then POS (Beagle requires sorted input)
  ord   <- order(snp_info$CHR, as.integer(snp_info$POS))
  si    <- snp_info[ord, , drop = FALSE]
  gm    <- geno_mat[, ord, drop = FALSE]   # individuals x sorted SNPs

  ns  <- nrow(gm)   # n samples
  np  <- ncol(gm)   # n SNPs

  if (verbose)
    message("[write_cleaned_vcf] Writing ", np, " SNPs x ", ns,
            " samples -> ", basename(out_file))

  # Dosage to GT lookup (vectorised)
  gt_lookup <- c("0" = "0/0", "1" = "0/1", "2" = "1/1")

  con <- gzfile(out_file, "wb")
  on.exit(try(close(con), silent = TRUE), add = TRUE)

  # VCF header
  writeLines(c(
    "##fileformat=VCFv4.2",
    paste0("##source=LDxBlocks_write_cleaned_vcf"),
    paste0("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">"),
    paste(c("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER",
            "INFO", "FORMAT", sample_ids), collapse = "\t")
  ), con = con)

  # Write one VCF record per SNP
  # Encode all genotypes for this SNP as a character vector, then paste
  for (k in seq_len(np)) {
    dosage <- gm[, k]
    gt_vec <- ifelse(
      is.na(dosage),
      "./.",
      gt_lookup[as.character(as.integer(dosage))]
    )
    gt_vec[is.na(gt_vec)] <- "./."   # guard for out-of-range values

    ref_k <- if (has_ref && !is.na(si$REF[k]) && nzchar(si$REF[k]))
      as.character(si$REF[k]) else "A"
    alt_k <- if (has_alt && !is.na(si$ALT[k]) && nzchar(si$ALT[k]))
      as.character(si$ALT[k]) else "T"

    # Ensure REF != ALT (some datasets have REF=ALT after normalisation)
    if (ref_k == alt_k) alt_k <- if (ref_k == "A") "T" else "A"

    line <- paste(
      c(as.character(si$CHR[k]),
        as.integer(si$POS[k]),
        as.character(si$SNP[k]),
        ref_k, alt_k, ".", "PASS", ".", "GT",
        gt_vec),
      collapse = "\t"
    )
    writeLines(line, con = con)
  }

  close(con)
  if (verbose)
    message("[write_cleaned_vcf] Done: ", out_file)
  invisible(out_file)
}


# -- Internal: read phased VCF and cache hap1/hap2/dosage as bigmemory --------
#
# After phase_with_beagle() produces a phased VCF, this function:
#   1. Checks whether bigmemory-backed caches already exist and are fresh.
#   2. If so, reattaches them instantly (no re-reading of the VCF).
#   3. Otherwise, reads the phased VCF via read_phased_vcf(), aligns samples
#      and SNPs to the cleaned/imputed set, then writes three bigmemory backends:
#        ldxblocks_bm_phased_hap1  (individuals x SNPs, char, 0/1)
#        ldxblocks_bm_phased_hap2  (individuals x SNPs, char, 0/1)
#        ldxblocks_bm_phased_dos   (individuals x SNPs, char, 0/1/2)
#      plus snp_info, sample_ids, and a fingerprint .rds in bigmemory_path.
#
# Returns a list(hap1, hap2, dosage, snp_info, sample_ids, phased=TRUE,
#                phase_method="beagle") compatible with extract_haplotypes()
# isp=TRUE path.  When use_bigmemory=FALSE the matrices are returned as plain
# R matrices (same behaviour as before, but still aligned and cached in RAM).
#' @noRd
.read_and_cache_phased_vcf <- function(
    phased_vcf,          # path to phased VCF.gz from phase_with_beagle()
    snp_info_filtered,   # cleaned/imputed snp_info (target SNP set)
    sample_ids_target,   # cleaned/imputed sample order
    chr_filter,          # character vector or NULL (from pipeline chr arg)
    norm_key_fn,         # .norm_key closure from pipeline scope
    bigmemory_path,
    bigmemory_type = "char",
    use_bigmemory  = FALSE,
    verbose        = TRUE
) {
  .msg <- function(...) if (verbose) message("[phased cache] ", ...)

  # Fingerprint: phased VCF path + mtime + target SNP count + sample count.
  # Built FIRST so it can be hashed for the stem name.
  vcf_mtime <- tryCatch(file.info(phased_vcf)$mtime, error = function(e) NA)
  curr_fp <- list(
    phased_vcf       = normalizePath(phased_vcf, mustWork = FALSE),
    vcf_mtime        = as.character(vcf_mtime),
    n_snps_target    = nrow(snp_info_filtered),
    n_samples_target = length(sample_ids_target),
    bigmemory_type   = bigmemory_type
  )

  # Hash the fingerprint to produce a stable, collision-safe stem name.
  # Different VCF files / SNP sets / samples get different stems, so
  # stale backends from prior runs cannot silently be reattached.
  run_tag   <- substr(
    digest::digest(curr_fp, algo = "md5", serialize = TRUE),
    1L, 10L
  )
  stem_base <- file.path(
    normalizePath(bigmemory_path, mustWork = FALSE),
    paste0("ldxblocks_bm_phased_", run_tag)
  )
  stem_base <- gsub("\\\\", "/", stem_base)

  si_file   <- paste0(stem_base, "_snpinfo.rds")
  sid_file  <- paste0(stem_base, "_sampleids.rds")
  fp_file   <- paste0(stem_base, "_params.rds")

  stems <- c(hap1 = paste0(stem_base, "_hap1"),
             hap2 = paste0(stem_base, "_hap2"),
             dos  = paste0(stem_base, "_dos"))


  # Check whether all cache files exist and fingerprint matches
  all_desc_ok <- all(file.exists(paste0(stems, ".desc")))
  fp_ok       <- FALSE
  if (all_desc_ok && file.exists(fp_file)) {
    saved_fp <- tryCatch(readRDS(fp_file), error = function(e) NULL)
    fp_ok    <- identical(saved_fp, curr_fp)
  }

  # --- Reattach path: all cache files fresh ---
  if (isTRUE(use_bigmemory) && all_desc_ok && fp_ok &&
      file.exists(si_file) && file.exists(sid_file)) {
    .msg("Reattaching cached phased backends from ", bigmemory_path)
    si_c  <- readRDS(si_file)
    sid_c <- readRDS(sid_file)
    attach_bm <- function(nm) {
      bm  <- bigmemory::attach.big.matrix(paste0(stems[[nm]], ".desc"))
      out <- bm[]
      rownames(out) <- sid_c
      colnames(out) <- si_c$SNP
      out
    }
    return(list(
      hap1       = attach_bm("hap1"),
      hap2       = attach_bm("hap2"),
      dosage     = attach_bm("dos"),
      snp_info   = si_c,
      sample_ids = sid_c,
      phased     = TRUE,
      phase_method = "beagle"
    ))
  }

  # --- Build path: read VCF, align, save ---
  # Single authoritative cleanup: remove all stale backend files before rebuild.
  # Uses unlink(force=TRUE) which is silent on missing files and works on
  # Windows locked handles better than file.remove().
  cleanup_stem <- function(stem) {
    files <- c(
      paste0(stem, ".bin"),
      paste0(stem, ".desc"),
      paste0(stem, "_params.rds"),
      paste0(stem, "_sampleids.rds"),
      paste0(stem, "_snpinfo.rds")
    )
    for (f in files) {
      if (file.exists(f)) {
        unlink(f, force = TRUE)
        if (file.exists(f))
          stop("Could not remove stale backend file: ", f, call. = FALSE)
      }
    }
  }
  invisible(lapply(stems, cleanup_stem))
  for (f in c(si_file, sid_file, fp_file))
    unlink(f, force = TRUE)

  .msg("Reading phased VCF: ", basename(phased_vcf))
  gp <- read_phased_vcf(phased_vcf, min_maf = 0.0, verbose = verbose)

  # Reorder samples to match cleaned/imputed matrix
  sid_match <- match(sample_ids_target, gp$sample_ids)
  if (anyNA(sid_match))
    stop("Phased VCF sample set does not match the cleaned genotype sample set.",
         call. = FALSE)
  gp$hap1       <- gp$hap1[, sid_match, drop = FALSE]
  gp$hap2       <- gp$hap2[, sid_match, drop = FALSE]
  gp$dosage     <- gp$dosage[, sid_match, drop = FALSE]
  gp$sample_ids <- gp$sample_ids[sid_match]

  # Restrict to exact cleaned/imputed SNP set via composite key
  target_df       <- snp_info_filtered
  target_df$CHR   <- sub("^chr", "", as.character(target_df$CHR), ignore.case = TRUE)
  phased_df       <- gp$snp_info
  phased_df$CHR   <- sub("^chr", "", as.character(phased_df$CHR), ignore.case = TRUE)
  if (!is.null(chr_filter))
    phased_df <- phased_df[phased_df$CHR %in% chr_filter, , drop = FALSE]

  target_key <- norm_key_fn(target_df)
  phased_key <- norm_key_fn(phased_df)
  idx_match  <- match(target_key, phased_key)
  if (anyNA(idx_match))
    stop("Phased VCF missing ", sum(is.na(idx_match)),
         " SNP(s) from the cleaned/imputed marker set after CHR+POS+REF+ALT+SNP matching.",
         call. = FALSE)

  h1  <- t(gp$hap1[idx_match,  , drop = FALSE])   # individuals x SNPs
  h2  <- t(gp$hap2[idx_match,  , drop = FALSE])
  dos <- t(gp$dosage[idx_match, , drop = FALSE])
  si  <- phased_df[idx_match, , drop = FALSE]
  si$CHR <- sub("^chr", "", as.character(si$CHR), ignore.case = TRUE)
  rownames(h1) <- rownames(h2) <- rownames(dos) <- gp$sample_ids
  colnames(h1) <- colnames(h2) <- colnames(dos) <- si$SNP

  .msg("Phased matrix: ", nrow(h1), " ind x ", ncol(h1), " SNPs after alignment")

  # Save snp_info + sample_ids
  saveRDS(si,            si_file)
  saveRDS(gp$sample_ids, sid_file)

  # Save to bigmemory or return as plain matrices
  if (isTRUE(use_bigmemory) && requireNamespace("bigmemory", quietly = TRUE)) {
    .msg("Writing phased bigmemory backends (type = '", bigmemory_type, "') ...")
    make_bm <- function(mat, stem_nm) {
      .bin  <- file.path(bigmemory_path, paste0(basename(stems[[stem_nm]]), ".bin"))
      .desc <- file.path(bigmemory_path, paste0(basename(stems[[stem_nm]]), ".desc"))
      # Cleanup is handled by cleanup_stem() above; if files still exist here
      # something went wrong - fail fast rather than silently overwriting.
      if (file.exists(.bin) || file.exists(.desc))
        stop("Stale bigmemory backend still exists after cleanup: ",
             basename(.bin), call. = FALSE)
      bm <- bigmemory::as.big.matrix(
        matrix(as.integer(mat), nrow = nrow(mat), ncol = ncol(mat)),
        type        = bigmemory_type,
        backingfile = paste0(basename(stems[[stem_nm]]), ".bin"),
        backingpath = bigmemory_path,
        descriptorfile = paste0(basename(stems[[stem_nm]]), ".desc")
      )
      # Restore dimnames: bigmemory bm[] strips rownames and colnames.
      # Without this, sg (SNP IDs) ends up NULL in extract_haplotypes().
      out <- bm[]
      rownames(out) <- rownames(mat)
      colnames(out) <- colnames(mat)
      out
    }
    h1_out  <- make_bm(h1,  "hap1")
    h2_out  <- make_bm(h2,  "hap2")
    dos_out <- make_bm(dos, "dos")
    .msg("Phased backends written. Rerun to reattach without re-reading VCF.")
  } else {
    h1_out  <- h1
    h2_out  <- h2
    dos_out <- dos
  }

  saveRDS(curr_fp, fp_file)

  list(hap1 = h1_out, hap2 = h2_out, dosage = dos_out,
       snp_info = si, sample_ids = gp$sample_ids,
       phased = TRUE, phase_method = "beagle")
}


#' End-to-End Haplotype Block Pipeline
#'
#' @description
#' Single-call wrapper: one genotype file in, complete haplotype dataset out.
#' Handles format detection, MAF filtering, call-rate filtering, imputation,
#' LD block detection, optional Beagle phasing, haplotype extraction, diversity
#' analysis, feature matrix construction, and output writing.
#'
#' @section Phasing modes:
#' \describe{
#'   \item{\code{phase = FALSE} (default)}{
#'     Dosage-pattern haplotypes extracted directly from the imputed matrix.
#'     No external tools required. Fast. Suitable for genomic prediction.
#'     Each block entry is a multi-SNP dosage string - not a true gametic
#'     haplotype. Frequencies are individual-level pattern proportions.}
#'   \item{\code{phase = TRUE}}{
#'     Beagle 5.x called on the original input VCF after LD block detection,
#'     producing true statistically-inferred gametic haplotypes using
#'     population-LD across all markers. Haplotype strings become
#'     \code{"g1|g2"}. Frequencies are computed over \eqn{2N} gamete
#'     observations. Recommended for diversity analysis and biologically
#'     interpretable results.
#'
#'     \strong{Requirements:} \code{geno_source} must be VCF/VCF.gz.
#'     Place \code{beagle.jar} in \code{out_dir} or supply via
#'     \code{beagle_jar}. Download:
#'     \url{https://faculty.washington.edu/browning/beagle/beagle.html}}
#' }
#'
#' @section Why Beagle and why a cleaned VCF:
#' Beagle 5.x performs chromosome-level statistical phasing using population
#' LD across all markers simultaneously, producing true inferred gametic
#' haplotypes. This is the only supported phasing method in LDxBlocks.
#'
#' When \code{phase = TRUE}, the pipeline does \strong{not} pass
#' \code{geno_source} directly to Beagle. Instead it first writes a
#' cleaned VCF (\code{<geno_source_stem>_cleaned.vcf.gz} in \code{out_dir})
#' containing exactly the SNPs that survived MAF filtering, call-rate
#' filtering, chromosome subsetting, and malformed-record removal. This
#' guarantees that the phased VCF's SNP set is identical to the marker set
#' used for LD block detection, so \code{.read_and_cache_phased_vcf()} can
#' align phased gametes to blocks without SNP-count mismatches.
#'
#' @param geno_source  File path or \code{LDxBlocks_backend}. Supported
#'   formats: VCF (\code{.vcf}/\code{.vcf.gz}), HapMap (\code{.hmp.txt}),
#'   CSV, GDS, PLINK BED. \code{phase = TRUE} requires VCF/VCF.gz.
#' @param out_dir      Output directory. Default \code{"."}. When
#'   \code{phase = TRUE}, \code{beagle.jar} is expected here and the phased
#'   VCF is written here.
#' @param out_blocks   Path for LD block table CSV.
#' @param out_diversity Path for haplotype diversity table CSV.
#' @param out_hap_matrix Path for haplotype genotype matrix file.
#' @param hap_format   \code{"numeric"} (default) or \code{"character"}.
#' @param phase        If \code{TRUE}, phase with Beagle. Default \code{FALSE}.
#' @param beagle_jar   Path to \code{beagle.jar}.
#'   Default: \code{file.path(out_dir, "beagle.jar")}.
#' @param beagle_threads Beagle threads. Default \code{1L}.
#' @param java_path    Java executable. Default \code{"java"}.
#' @param beagle_java_mem_gb JVM heap in GB (\code{-Xmx}). Default \code{NULL}.
#' @param beagle_args  Extra Beagle argument string. Default \code{""}.
#' @param beagle_ref_panel Phased reference VCF path. Default \code{NULL}.
#' @param beagle_map_file  Genetic map file path. Default \code{NULL}.
#' @param beagle_chrom Restrict Beagle to one chromosome. Default \code{NULL}
#'   (inherits \code{chr} if set).
#' @param beagle_seed  Integer seed for reproducibility. Default \code{NULL}.
#' @param beagle_burnin Burn-in iterations. Default \code{NULL}.
#' @param beagle_iterations Phasing iterations. Default \code{NULL}.
#' @param beagle_window Window size. Default \code{NULL}.
#' @param beagle_overlap Window overlap. Default \code{NULL}.
#' @param maf_cut      Minimum MAF. Default \code{0.05}.
#' @param impute       \code{"mean_rounded"} (default), \code{"mode"}, or
#'   \code{"none"}.
#' @param min_callrate Minimum per-SNP call rate. Default \code{0.0}.
#' @param CLQcut       r\eqn{^2} threshold for CLQD. Default \code{0.5}.
#' @param method       \code{"r2"} (default) or \code{"rV2"}.
#' @param kin_method   \code{"chol"} (default) or \code{"eigen"}.
#' @param CLQmode      \code{"Density"} (default), \code{"Maximal"},
#'   \code{"Louvain"}, or \code{"Leiden"}.
#' @param clstgap      Max bp gap within clique. Default \code{40000L}.
#' @param split        Split cliques at largest gap. Default \code{FALSE}.
#' @param appendrare   Append rare SNPs to nearest block. Default \code{FALSE}.
#' @param singleton_as_block Return singletons as blocks. Default \code{FALSE}.
#' @param checkLargest Dense-core pre-pass. Default \code{FALSE}.
#' @param digits       Round r\eqn{^2} (\code{-1L} = off). Default \code{-1L}.
#' @param leng         Boundary scan half-window (SNPs). Default \code{200L}.
#' @param subSegmSize  Max SNPs per sub-segment. Default \code{1500L}.
#' @param n_threads    OpenMP threads. Default \code{1L}.
#' @param min_snps_chr Skip chromosomes below this SNP count. Default \code{10L}.
#' @param max_bp_distance Max bp for r\eqn{^2} (\code{0L} = all). Default \code{0L}.
#' @param min_snps_block Min SNPs per haplotype block. Default \code{3L}.
#' @param top_n        Max alleles per block (\code{NULL} = all above
#'   \code{min_freq}). Default \code{NULL}.
#' @param min_freq     Min haplotype allele frequency. Default \code{0.01}.
#' @param scale_hap_matrix Scale haplotype matrix columns. Default \code{FALSE}.
#' @param chr          Chromosomes to process (\code{NULL} = all). Default
#'   \code{NULL}.
#' @param clean_malformed Remove malformed VCF lines. Default \code{FALSE}.
#' @param use_bigmemory File-backed bigmemory store. Default \code{FALSE}.
#' @param bigmemory_path Directory for backing files. Default \code{tempdir()}.
#' @param bigmemory_type \code{"char"} (default), \code{"short"},
#'   or \code{"double"}.
#' @param verbose      Print progress. Default \code{TRUE}.
#'
#' @return Named list (invisibly): \code{blocks}, \code{diversity},
#'   \code{hap_matrix}, \code{hap_matrix_info}, \code{haplotypes},
#'   \code{geno_matrix}, \code{snp_info_filtered}, \code{phased_vcf},
#'   \code{phased_backend_desc}, \code{phase_method},
#'   \code{n_blocks}, \code{n_hap_columns}.
#'
#' @examples
#' \donttest{
#' geno_file <- system.file("extdata", "example_genotypes_numeric.csv",
#'                           package = "LDxBlocks")
#' res <- run_ldx_pipeline(
#'   geno_source = geno_file, out_dir = tempdir(),
#'   out_blocks = tempfile(fileext=".csv"),
#'   out_diversity = tempfile(fileext=".csv"),
#'   out_hap_matrix = tempfile(fileext=".csv"),
#'   phase = FALSE, maf_cut = 0.05, CLQcut = 0.5,
#'   leng = 10L, subSegmSize = 80L, verbose = FALSE
#' )
#' \dontrun{
#' # With Beagle phasing (place beagle.jar in out_dir first):
#' res2 <- run_ldx_pipeline(
#'   geno_source = "data.vcf.gz", out_dir = "results/",
#'   out_blocks = "results/blocks.csv",
#'   out_diversity = "results/diversity.csv",
#'   out_hap_matrix = "results/hap_matrix.csv",
#'   phase = TRUE, beagle_threads = 4L,
#'   beagle_java_mem_gb = 8, beagle_seed = 42L
#' )
#' }}
#'
#' @seealso \code{\link{run_Big_LD_all_chr}}, \code{\link{extract_haplotypes}},
#'   \code{\link{build_haplotype_feature_matrix}},
#'   \code{\link{run_haplotype_prediction}}
#' @export
run_ldx_pipeline <- function(
    geno_source,
    out_dir        = ".",
    out_blocks,
    out_diversity,
    out_hap_matrix,
    hap_format     = c("numeric", "character"),

    # -- Phasing --------------------------------------------------------------
    phase              = FALSE,
    beagle_jar         = NULL,
    beagle_threads     = 1L,
    java_path          = "java",
    beagle_java_mem_gb = NULL,
    beagle_args        = "",
    beagle_ref_panel   = NULL,
    beagle_map_file    = NULL,
    beagle_chrom       = NULL,
    beagle_seed        = NULL,
    beagle_burnin      = NULL,
    beagle_iterations  = NULL,
    beagle_window      = NULL,
    beagle_overlap     = NULL,

    # -- Filtering & imputation -----------------------------------------------
    maf_cut        = 0.05,
    impute         = c("mean_rounded", "mode", "none"),
    min_callrate   = 0.0,

    # -- LD block detection ---------------------------------------------------
    CLQcut             = 0.5,
    method             = c("r2", "rV2"),
    kin_method         = "chol",
    CLQmode            = c("Density", "Maximal", "Louvain", "Leiden"),
    leng               = 200L,
    subSegmSize        = 1500L,
    clstgap            = 40000L,
    split              = FALSE,
    appendrare         = FALSE,
    singleton_as_block = FALSE,
    checkLargest       = FALSE,
    digits             = -1L,
    n_threads          = 1L,
    min_snps_chr       = 10L,
    max_bp_distance    = 0L,

    # -- Haplotype extraction -------------------------------------------------
    min_snps_block   = 3L,
    top_n            = NULL,
    min_freq         = 0.01,
    scale_hap_matrix = FALSE,

    # -- General --------------------------------------------------------------
    chr             = NULL,
    clean_malformed = FALSE,
    use_bigmemory   = FALSE,
    bigmemory_path  = tempdir(),
    bigmemory_type  = "char",
    verbose         = TRUE
) {
  hap_format     <- match.arg(hap_format)
  method         <- match.arg(method)
  CLQmode        <- match.arg(CLQmode)
  impute         <- match.arg(impute)
  bigmemory_type <- match.arg(bigmemory_type, choices = c("char", "short", "double"))

  .log <- function(...) {
    if (verbose)
      message(sprintf("[%s] %s", format(Sys.time(), "%H:%M:%S"), paste0(...)))
  }

  # SNP matching key: primary = CHR+POS+SNP, secondary = REF+ALT.
  # Beagle can normalise alleles (e.g. strand-flip, trim), which breaks
  # exact CHR+POS+REF+ALT+SNP matching. We match primarily by CHR+POS+SNP
  # and use REF/ALT only as a tie-breaker / sanity check, making alignment
  # robust to Beagle allele normalisation.
  .norm_key <- function(df) {
    req <- c("CHR", "POS", "SNP")
    miss <- setdiff(req, names(df))
    if (length(miss))
      stop("Missing required columns for alignment: ", paste(miss, collapse = ","), call. = FALSE)
    paste(
      .norm_chr_hap(df$CHR),
      as.character(df$POS),
      as.character(df$SNP),
      sep = "||"
    )
  }

  .make_backed_backend <- function(
    geno_mat,
    snp_info,
    sample_ids,
    stem_name,
    fingerprint = NULL,
    verbose = TRUE
  ) {
    .make_backed_backend_impl(
      geno_mat = geno_mat, snp_info = snp_info, sample_ids = sample_ids,
      stem_name = stem_name, fingerprint = fingerprint,
      use_bigmemory = use_bigmemory, bigmemory_path = bigmemory_path,
      bigmemory_type = bigmemory_type, verbose = verbose
    )
  }

  # -- Validate phasing requirements early ------------------------------------
  if (isTRUE(phase)) {
    if (inherits(geno_source, "LDxBlocks_backend"))
      stop(
        "phase = TRUE requires geno_source to be a VCF/VCF.gz file path. ",
        "Pre-built backends cannot be passed to Beagle.",
        call. = FALSE
      )

    if (!is.character(geno_source) ||
        length(geno_source) != 1L ||
        !grepl("\\.vcf(\\.gz)?$", geno_source, ignore.case = TRUE))
      stop(
        "phase = TRUE requires geno_source to be a VCF or VCF.gz file path.",
        call. = FALSE
      )

    if (!file.exists(geno_source))
      stop("geno_source not found: ", geno_source, call. = FALSE)

    if (is.null(beagle_jar))
      beagle_jar <- file.path(out_dir, "beagle.jar")
    if (!file.exists(beagle_jar))
      stop("beagle.jar not found at: ", beagle_jar,
           "\nPlace beagle.jar in out_dir ('", out_dir, "') or supply via beagle_jar.\n",
           "Download: https://faculty.washington.edu/browning/beagle/beagle.html",
           call. = FALSE)
  }

  # -- Step 1: Open backend ---------------------------------------------------
  if (inherits(geno_source, "LDxBlocks_backend")) {
    be <- geno_source
    .log("Using pre-built backend: ", be$type, " | ",
         be$n_samples, " ind | ", be$n_snps, " SNPs")
  } else {
    if (!is.character(geno_source) || length(geno_source) != 1L)
      stop("geno_source must be a file path or an LDxBlocks_backend.", call. = FALSE)

    if (isTRUE(use_bigmemory)) {
      if (!requireNamespace("bigmemory", quietly = TRUE))
        stop("use_bigmemory = TRUE requires the 'bigmemory' package.", call. = FALSE)

      bigmemory_path <- normalizePath(bigmemory_path, mustWork = FALSE)
      bigmemory_path <- gsub("\\\\", "/", bigmemory_path, fixed = TRUE)

      bm_stem      <- file.path(bigmemory_path, "ldxblocks_bm")
      bm_bin_file  <- paste0(bm_stem, ".bin")
      bm_desc_file <- paste0(bm_stem, ".desc")
      bm_si_file   <- paste0(bm_stem, "_snpinfo.rds")

      bm_bin_noext <- sub("\\.bin$", "", bm_bin_file)
      bm_bin_ok    <- file.exists(bm_bin_file) || file.exists(bm_bin_noext)
      bm_all_ok    <- bm_bin_ok && file.exists(bm_desc_file) && file.exists(bm_si_file)
      bm_partial   <- (bm_bin_ok || file.exists(bm_desc_file) || file.exists(bm_si_file)) && !bm_all_ok

      if (bm_partial) {
        for (f in c(bm_bin_file, bm_bin_noext, bm_desc_file, bm_si_file))
          if (file.exists(f)) {
            .log("[bigmemory] Removing stale: ", basename(f))
            file.remove(f)
          }
        bm_all_ok <- FALSE
      }

      if (bm_all_ok) {
        .log("[bigmemory] Reattaching: ", basename(bm_bin_file))
        si_bm <- readRDS(bm_si_file)
        be <- read_geno_bigmemory(
          source      = bm_desc_file,
          snp_info    = si_bm,
          backingfile = basename(bm_stem),
          backingpath = bigmemory_path
        )
      } else {
        .log("[bigmemory] Building file-backed matrix (type = '", bigmemory_type, "') ...")
        be_tmp <- read_geno(geno_source, clean_malformed = clean_malformed, verbose = verbose)
        be <- read_geno_bigmemory(
          source      = be_tmp,
          snp_info    = be_tmp$snp_info,
          backingfile = basename(bm_stem),
          backingpath = bigmemory_path,
          type        = bigmemory_type,
          verbose     = verbose
        )
        close_backend(be_tmp)
        saveRDS(be$snp_info, bm_si_file)
        saveRDS(be$sample_ids, paste0(bm_stem, "_sampleids.rds"))
      }
      on.exit(close_backend(be), add = TRUE)
    } else {
      .log("Opening: ", basename(geno_source))
      be <- read_geno(geno_source, clean_malformed = clean_malformed, verbose = verbose)
      on.exit(close_backend(be), add = TRUE)
    }

    .log("Backend: ", be$type, " | ", be$n_samples, " ind | ", be$n_snps, " SNPs")
  }

  # -- Step 2: Streaming MAF pre-screen ---------------------------------------
  orig_snp_info     <- be$snp_info
  snp_info_filtered <- .maf_filter_backend(be, maf_cut = maf_cut, verbose = verbose)

  if (nrow(snp_info_filtered) == 0L)
    stop("No SNPs remain after MAF pre-screen. Lower maf_cut.", call. = FALSE)

  be$snp_info <- snp_info_filtered
  be$n_snps   <- nrow(snp_info_filtered)

  # -- Step 3: Chromosome subset ----------------------------------------------
  if (!is.null(chr)) {
    chr <- .norm_chr_hap(as.character(chr))
    be$snp_info$CHR <- .norm_chr_hap(be$snp_info$CHR)
    be$snp_info <- be$snp_info[be$snp_info$CHR %in% chr, , drop = FALSE]
    be$n_snps   <- nrow(be$snp_info)
    .log("Chr filter: ", be$n_snps, " SNPs on chr ", paste(chr, collapse = ", "))
  }

  # -- Step 4: Load genotype matrix -------------------------------------------
  .log("Loading filtered genotype matrix ...")
  full_idx <- match(be$snp_info$SNP, orig_snp_info$SNP)
  geno_mat <- read_chunk(be, full_idx)
  rownames(geno_mat) <- be$sample_ids
  colnames(geno_mat) <- be$snp_info$SNP
  .log("Genotype matrix: ", nrow(geno_mat), " x ", ncol(geno_mat))

  # -- Step 4.5a: Call-rate filter --------------------------------------------
  if (min_callrate > 0) {
    cr_res <- impute_and_filter_cpp(
      geno         = matrix(as.integer(geno_mat), nrow = nrow(geno_mat)),
      min_callrate = min_callrate,
      method       = 0L
    )

    if (cr_res$n_filtered > 0L) {
      keep_cr           <- as.logical(cr_res$keep)
      geno_mat          <- geno_mat[, keep_cr, drop = FALSE]
      snp_info_filtered <- snp_info_filtered[keep_cr, , drop = FALSE]
      be$snp_info       <- snp_info_filtered
      be$n_snps         <- nrow(snp_info_filtered)

      .log(sprintf(
        "Call-rate filter: removed %d SNPs (< %.0f%%) | Remaining: %d",
        cr_res$n_filtered, min_callrate * 100, ncol(geno_mat)
      ))
    } else {
      .log(sprintf("Call-rate filter: all %d SNPs passed.", ncol(geno_mat)))
    }
  }

  # -- Step 4.5b: MAF filter (authoritative) ----------------------------------
  keep_maf      <- maf_filter_cpp(geno_mat, maf_cut = maf_cut)
  n_maf_removed <- sum(!keep_maf)

  if (n_maf_removed > 0L) {
    geno_mat          <- geno_mat[, keep_maf, drop = FALSE]
    snp_info_filtered <- snp_info_filtered[keep_maf, , drop = FALSE]
    be$snp_info       <- snp_info_filtered
    be$n_snps         <- nrow(snp_info_filtered)

    .log(sprintf(
      "MAF filter (>= %.2f): removed %d SNPs | Remaining: %d",
      maf_cut, n_maf_removed, ncol(geno_mat)
    ))
  }

  if (ncol(geno_mat) == 0L)
    stop("No SNPs remain after MAF filtering.", call. = FALSE)

  # -- Step 4.5c: Imputation --------------------------------------------------
  n_na <- sum(is.na(geno_mat))

  if (impute == "none" && n_na > 0L)
    stop(
      "impute = 'none' but ", n_na, " missing genotypes remain. ",
      "Set impute = 'mean_rounded' or 'mode'.",
      call. = FALSE
    )

  if (impute != "none" && n_na > 0L) {
    .log(sprintf("Imputation (%s): %d missing values ...", impute, n_na))
    imp_res <- impute_and_filter_cpp(
      geno         = matrix(as.integer(geno_mat), nrow = nrow(geno_mat)),
      min_callrate = 0.0,
      method       = if (impute == "mode") 1L else 0L
    )
    geno_mat <- imp_res$geno_imputed
    storage.mode(geno_mat) <- "numeric"
    rownames(geno_mat) <- be$sample_ids
    colnames(geno_mat) <- snp_info_filtered$SNP

    .log(sprintf(
      "Imputation: filled %d values | Remaining NA: %d",
      imp_res$n_imputed, sum(is.na(geno_mat))
    ))
  }

  # -- Step 4.6: Build cleaned/imputed backend --------------------------------
  .log("Building imputed backend (", nrow(geno_mat), " x ", ncol(geno_mat), ") ...")
  be_imputed <- .make_imputed_backend(
    geno_mat       = geno_mat,
    snp_info       = snp_info_filtered,
    sample_ids     = be$sample_ids,
    use_bigmemory  = use_bigmemory,
    bigmemory_path = bigmemory_path,
    bigmemory_type = bigmemory_type,
    maf_cut        = maf_cut,
    min_callrate   = min_callrate,
    impute_method  = impute,
    verbose        = verbose
  )
  on.exit(close_backend(be_imputed), add = TRUE)

  # -- Step 5: LD block detection ---------------------------------------------
  .log("Running genome-wide LD block detection ...")
  blocks <- run_Big_LD_all_chr(
    geno_matrix        = be_imputed,
    snp_info           = be_imputed$snp_info,
    method             = method,
    kin_method         = kin_method,
    CLQcut             = CLQcut,
    CLQmode            = CLQmode,
    leng               = leng,
    subSegmSize        = subSegmSize,
    clstgap            = clstgap,
    split              = split,
    MAFcut             = maf_cut,
    appendrare         = appendrare,
    singleton_as_block = singleton_as_block,
    checkLargest       = checkLargest,
    digits             = digits,
    n_threads          = n_threads,
    min_snps_chr       = min_snps_chr,
    max_bp_distance    = max_bp_distance,
    verbose            = verbose
  )

  if (is.null(blocks) || nrow(blocks) == 0L)
    stop("No LD blocks detected. Try lowering CLQcut or adjusting parameters.", call. = FALSE)

  .log("Detected ", nrow(blocks), " LD blocks")
  data.table::fwrite(blocks, file = out_blocks, sep = ",", quote = FALSE, na = "NA")
  .log("Block table written: ", out_blocks)

  # -- Step 5b: Optional Beagle phasing ---------------------------------------
  # Phasing occurs AFTER block detection (unphased dosage matrix is correct for
  # LD) and BEFORE haplotype extraction (which needs gametic strings).
  #
  # .read_and_cache_phased_vcf() handles:
  #   - Running phase_with_beagle() to produce the phased VCF
  #   - Reading via read_phased_vcf()
  #   - Aligning samples and SNPs to the cleaned/imputed set
  #   - Saving hap1, hap2, and dosage as bigmemory-backed matrices in
  #     bigmemory_path (ldxblocks_bm_phased_hap1, _hap2, _dos)
  #   - Fingerprint-based restart: reattaches from disk without re-reading
  #     the phased VCF on subsequent runs
  phased_vcf_path     <- NULL
  phased_backend_desc <- NULL
  phase_method        <- "unphased"
  geno_for_hap        <- be_imputed
  out_snp_info        <- be_imputed$snp_info

  if (isTRUE(phase)) {
    # beagle_chrom: use explicit value or inherit from chr filter
    phase_chr_use <- beagle_chrom
    if (is.null(phase_chr_use) && !is.null(chr))
      phase_chr_use <- paste(chr, collapse = ",")

    out_prefix <- file.path(
      out_dir,
      sub("\\.vcf(\\.gz)?$", "_phased", basename(geno_source), ignore.case = TRUE)
    )

    # Write a filtered VCF from geno_mat (post-MAF/callrate/chr/clean_malformed
    # filtering) for use as Beagle input.  This is ALWAYS done, not only when
    # clean_malformed = TRUE, because:
    #   - geno_source may contain SNPs that failed MAF or call-rate filters
    #   - geno_source may contain chromosomes excluded by the chr= argument
    #   - geno_source may contain malformed records if clean_malformed = TRUE
    # Feeding geno_source directly to Beagle would produce a phased VCF with a
    # different SNP set than snp_info_filtered, causing .read_and_cache_phased_vcf()
    # to report missing SNPs and fail alignment.
    #
    # The cleaned VCF contains exactly the SNPs in snp_info_filtered, sorted
    # by CHR and POS, encoded as 0/0 | 0/1 | 1/1 | ./. GT fields.
    cleaned_vcf_for_beagle <- file.path(
      out_dir,
      sub("\\.vcf(\\.gz)?$", "_cleaned.vcf.gz",
          basename(geno_source), ignore.case = TRUE)
    )
    .log("Writing cleaned VCF for Beagle input: ", basename(cleaned_vcf_for_beagle))
    .write_cleaned_vcf(
      geno_mat   = geno_mat,
      snp_info   = snp_info_filtered,
      sample_ids = be$sample_ids,
      out_file   = cleaned_vcf_for_beagle,
      verbose    = verbose
    )

    .log("Phasing with Beagle -> ", out_prefix, ".vcf.gz ...")
    phased_vcf_path <- phase_with_beagle(
      input_vcf   = cleaned_vcf_for_beagle,
      out_prefix  = out_prefix,
      beagle_jar  = beagle_jar,
      java_path   = java_path,
      java_mem_gb = beagle_java_mem_gb,
      nthreads    = as.integer(beagle_threads),
      ref_panel   = beagle_ref_panel,
      map_file    = beagle_map_file,
      chrom       = phase_chr_use,
      seed        = beagle_seed,
      burnin      = beagle_burnin,
      iterations  = beagle_iterations,
      window      = beagle_window,
      overlap     = beagle_overlap,
      beagle_args = beagle_args,
      verbose     = verbose
    )

    # Read phased VCF, align, and cache hap1/hap2/dosage to bigmemory.
    # On restart with the same bigmemory_path and matching fingerprint,
    # the VCF is NOT re-read - backends are reattached from disk instantly.
    geno_for_hap <- .read_and_cache_phased_vcf(
      phased_vcf         = phased_vcf_path,
      snp_info_filtered  = snp_info_filtered,
      sample_ids_target  = be$sample_ids,
      chr_filter         = chr,
      norm_key_fn        = .norm_key,
      bigmemory_path     = bigmemory_path,
      bigmemory_type     = bigmemory_type,
      use_bigmemory      = use_bigmemory,
      verbose            = verbose
    )

    out_snp_info <- geno_for_hap$snp_info

    if (isTRUE(use_bigmemory))
      phased_backend_desc <- file.path(
        normalizePath(bigmemory_path, mustWork = FALSE),
        "ldxblocks_bm_phased_dos.desc"
      )

    phase_method <- "beagle"
    .log("Phased data ready: ", nrow(geno_for_hap$snp_info), " SNPs x ",
         length(geno_for_hap$sample_ids), " individuals")
  }

  # -- Step 6: Haplotype extraction -------------------------------------------
  .log("Extracting haplotypes (min_snps = ", min_snps_block, ") ...")
  haplotypes <- extract_haplotypes(
    geno     = geno_for_hap,
    snp_info = out_snp_info,
    blocks   = blocks,
    chr      = NULL,
    min_snps = min_snps_block
  )

  n_hap_blocks <- length(haplotypes)
  .log("Haplotypes extracted: ", n_hap_blocks, " blocks | phased: ", isTRUE(phase))

  hap_bi <- attr(haplotypes, "block_info")
  if (!is.null(hap_bi) && nrow(hap_bi) > 0L) {
    bad_sing <- hap_bi[
      !is.na(hap_bi$start_bp) & !is.na(hap_bi$end_bp) & !is.na(hap_bi$n_snps) &
        hap_bi$start_bp == hap_bi$end_bp & hap_bi$n_snps > 1L, , drop = FALSE
    ]
    if (nrow(bad_sing) > 0L)
      stop("[LDxBlocks] Singleton-coordinate mismatch in ", nrow(bad_sing),
           " block(s). Run devtools::install() to recompile.", call. = FALSE)
    .log("QC: singleton-mismatch check passed.")
  }

  if (n_hap_blocks == 0L)
    stop("No haplotype blocks produced. Check min_snps_block vs block sizes.", call. = FALSE)

  # -- Step 7: Haplotype diversity --------------------------------------------
  .log("Computing haplotype diversity ...")
  diversity <- compute_haplotype_diversity(haplotypes)
  write_haplotype_diversity(diversity, out_diversity, append_summary = TRUE, verbose = FALSE)
  .log("Diversity table written: ", out_diversity)

  # -- Step 8: Build haplotype feature matrix ---------------------------------
  .log("Building haplotype feature matrix (top_n = ", top_n, ", min_freq = ", min_freq, ") ...")
  feat_out <- build_haplotype_feature_matrix(
    haplotypes     = haplotypes,
    top_n          = top_n,
    min_freq       = min_freq,
    scale_features = scale_hap_matrix
  )

  hap_matrix <- feat_out$matrix
  hap_info   <- feat_out$info
  n_hap_cols <- ncol(hap_matrix)
  .log("Feature matrix: ", nrow(hap_matrix), " ind x ", n_hap_cols, " haplotype allele columns")

  # -- Step 9: Write haplotype matrix -----------------------------------------
  .log("Writing haplotype matrix (format = ", hap_format, ") ...")

  if (identical(hap_format, "numeric")) {
    write_haplotype_numeric(
      hap_matrix = hap_matrix,
      out_file   = out_hap_matrix,
      haplotypes = haplotypes,
      snp_info   = out_snp_info,
      hap_info   = hap_info,
      sep        = "	",
      verbose    = verbose
    )
  } else {
    write_haplotype_character(
      haplotypes = haplotypes,
      snp_info   = out_snp_info,
      out_file   = out_hap_matrix,
      min_freq   = min_freq,
      top_n      = top_n,
      verbose    = verbose
    )
  }

  .log("Haplotype matrix written: ", out_hap_matrix)

  # -- Done -------------------------------------------------------------------
  .validate_hap_output(hap_matrix, hap_info, haplotypes)

  .log("Pipeline complete.")
  .log("  Phase method:        ", phase_method)
  .log("  Blocks:              ", nrow(blocks))
  .log("  Haplotype blocks:    ", n_hap_blocks)
  .log("  Haplotype columns:   ", n_hap_cols)
  .log("  Individuals:         ", nrow(hap_matrix))

  invisible(list(
    blocks              = blocks,
    diversity           = diversity,
    hap_matrix          = hap_matrix,
    hap_matrix_info     = hap_info,
    haplotypes          = haplotypes,
    geno_matrix         = geno_mat,
    snp_info_filtered   = snp_info_filtered,
    phased_vcf          = phased_vcf_path,
    phased_backend_desc = phased_backend_desc,
    phase_method        = phase_method,
    n_blocks            = nrow(blocks),
    n_hap_columns       = n_hap_cols
  ))
}
