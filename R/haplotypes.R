# ==============================================================================
# haplotypes.R - Comprehensive haplotype analysis for LDxBlocks
# Sections: 1-Phasing  2-Extraction  3-Diversity  4-QTL  5-Matrix  6-Writers
# ==============================================================================
.norm_chr_hap <- function(x) sub("^chr","",as.character(x),ignore.case=TRUE)
`%||%` <- function(x,y) if(is.null(x)) y else x

#' Read Pre-Phased VCF
#' @param vcf_file Path to phased VCF or VCF.gz with 0|1 GT fields.
#' @param min_maf Minimum MAF. Default 0.0.
#' @param verbose Logical. Default TRUE.
#' @return List: hap1, hap2 (SNPs x individuals, 0/1), dosage (0/1/2), snp_info, sample_ids.
#' @export
read_phased_vcf <- function(vcf_file, min_maf = 0.0, verbose = TRUE) {
  if (!is.numeric(min_maf) || length(min_maf) != 1L ||
      is.na(min_maf) || min_maf < 0 || min_maf > 0.5)
    stop("min_maf must be a single numeric in [0, 0.5].", call. = FALSE)
  if (!file.exists(vcf_file))
    stop("VCF not found: ", vcf_file, call. = FALSE)
  if (verbose) message("[read_phased_vcf] Reading: ", basename(vcf_file))

  # -- Step 1: Extract header lines to get sample IDs -------------------------
  # Read only until the #CHROM line - at most ~100 meta lines for any VCF.
  con <- gzfile(vcf_file, "r")
  on.exit(try(close(con), silent = TRUE), add = TRUE)
  header_line <- NULL
  repeat {
    ln <- readLines(con, n = 1L, warn = FALSE)
    if (!length(ln)) break
    if (startsWith(ln, "#CHROM")) { header_line <- ln; break }
  }
  close(con)
  if (is.null(header_line)) stop("No #CHROM header found.", call. = FALSE)

  cols <- strsplit(header_line, "\t", fixed = TRUE)[[1L]]
  sids <- cols[10:length(cols)]
  ns   <- length(sids)

  # -- Step 2: Read all data lines with data.table::fread() -------------------
  # Fast path: fread(cmd=) pipes through zcat + grep at C level (~50x speedup).
  # This requires zcat and grep on PATH. Check availability BEFORE calling fread
  # so that we never emit the shell-command-not-found WARNING that fread/shell()
  # produces on Windows when the command is missing. tryCatch only catches errors;
  # a shell warning leaks out even when the error handler runs the fallback.
  is_gz <- grepl("\\.gz$", vcf_file, ignore.case = TRUE)
  .use_shell <- nzchar(Sys.which(if (is_gz) "zcat" else "grep"))

  dt <- NULL

  if (.use_shell) {
    cmd <- if (is_gz) {
      paste0("zcat ", shQuote(vcf_file), " | grep -v '^#'")
    } else {
      paste0("grep -v '^#' ", shQuote(vcf_file))
    }
    dt <- tryCatch(
      data.table::fread(
        cmd          = cmd,
        header       = FALSE,
        sep          = "\t",
        colClasses   = "character",
        showProgress = FALSE,
        data.table   = FALSE
      ),
      error = function(e) {
        if (verbose)
          message("[read_phased_vcf] Shell fast-path failed; using R fallback reader")
        NULL
      }
    )
  }

  if (is.null(dt)) {
    if (verbose && !.use_shell)
      message("[read_phased_vcf] Shell tools not on PATH; using R fallback reader")
    pipe_con <- if (is_gz) gzcon(file(vcf_file, "rb")) else file(vcf_file, "r")
    on.exit(try(close(pipe_con), silent = TRUE), add = TRUE)
    all_lines  <- readLines(pipe_con, warn = FALSE)
    close(pipe_con)
    data_lines <- all_lines[!startsWith(all_lines, "#")]
    if (!length(data_lines))
      stop("No data lines found in VCF.", call. = FALSE)
    dt <- data.table::fread(
      text         = paste(data_lines, collapse = "\n"),
      header       = FALSE,
      sep          = "\t",
      colClasses   = "character",
      showProgress = FALSE,
      data.table   = FALSE
    )
  }

  nv <- nrow(dt)
  if (nv == 0L) stop("No variant rows found in phased VCF.", call. = FALSE)
  if (verbose) message("[read_phased_vcf] ", nv, " variants x ", ns, " samples")

  # Fixed column positions: 1=CHROM 2=POS 3=ID 4=REF 5=ALT 10..=sample GTs
  sc  <- .norm_chr_hap(dt[[1L]])
  sp  <- suppressWarnings(as.integer(dt[[2L]]))
  sid <- dt[[3L]]
  sr  <- dt[[4L]]
  sa  <- dt[[5L]]

  # -- Step 3: Vectorised GT parsing -----------------------------------------
  # Beagle 5.x always produces phased biallelic GT: "0|0" "0|1" "1|0" "1|1" ".|."
  # Alleles sit at fixed character positions 1 and 3 (pipe separator at 2).
  # substr() operates at C speed on the whole matrix - no strsplit loop needed.
  gt_mat  <- as.matrix(dt[, 10:(9L + ns), drop = FALSE])
  # Strip FORMAT subfields (e.g. "0|1:35:99" -> "0|1").
  # This makes the parser robust to VCFs with extra FORMAT fields.
  # Note: Beagle 5.x always writes FORMAT=GT only for the phased output,
  # so sub() is a no-op for Beagle-produced files.
  if (any(grepl(":", gt_mat, fixed = TRUE)))
    gt_mat <- matrix(sub(":.*$", "", gt_mat), nrow = nv, ncol = ns)
  # Parse alleles using sub() to support multi-character allele indices
  # (e.g. "10|2") as well as the standard single-character case.
  a1_char <- matrix(sub("\\|.*$", "", gt_mat), nrow = nv, ncol = ns)
  a2_char <- matrix(sub("^.*\\|",  "", gt_mat), nrow = nv, ncol = ns)

  missing_allele <- a1_char == "." | a1_char == "" |
    a2_char == "." | a2_char == ""
  h1 <- matrix(suppressWarnings(as.integer(a1_char)), nrow = nv, ncol = ns)
  h2 <- matrix(suppressWarnings(as.integer(a2_char)), nrow = nv, ncol = ns)
  h1[missing_allele] <- NA_real_
  h2[missing_allele] <- NA_real_

  rownames(h1) <- rownames(h2) <- sid
  colnames(h1) <- colnames(h2) <- sids
  dos <- h1 + h2

  # -- Step 4: Synthesise CHR_POS IDs for dot/empty ID fields -----------------
  dot_or_empty <- is.na(sid) | sid == "." | !nzchar(sid)
  if (any(dot_or_empty)) {
    sid[dot_or_empty] <- paste0(sc[dot_or_empty], "_", sp[dot_or_empty])
    rownames(h1) <- rownames(h2) <- rownames(dos) <- sid
    if (verbose)
      message("[read_phased_vcf] Synthesised ", sum(dot_or_empty),
              " SNP ID(s) as CHR_POS (original ID was '.' or empty)")
  }

  snp_info <- data.frame(SNP = sid, CHR = sc, POS = sp,
                         REF = sr,  ALT = sa,
                         stringsAsFactors = FALSE)

  # -- Step 5: Optional MAF filter --------------------------------------------
  if (min_maf > 0) {
    af  <- rowMeans(dos, na.rm = TRUE) / 2
    maf <- pmin(af, 1 - af)
    ok  <- !is.na(maf) & maf >= min_maf
    h1       <- h1[ok,  , drop = FALSE]
    h2       <- h2[ok,  , drop = FALSE]
    dos      <- dos[ok, , drop = FALSE]
    snp_info <- snp_info[ok, ]
    if (verbose)
      message("[read_phased_vcf] After MAF >= ", min_maf, ": ", sum(ok), " variants")
  }

  list(hap1 = h1, hap2 = h2, dosage = dos,
       snp_info     = snp_info,
       sample_ids   = sids,
       phased       = TRUE,
       phase_method = "vcf_phased")
}


#' Statistical Phasing via Beagle 5.x
#'
#' @description
#' Wrapper around Beagle 5.x for statistical haplotype phasing. Calls
#' \code{java -jar beagle.jar} and writes the phased VCF.gz to disk.
#' This is the recommended phasing method in LDxBlocks: Beagle performs
#' chromosome-level statistical phasing using population-LD across all
#' markers simultaneously, producing true inferred haplotypes suitable
#' for block-level frequency estimation and diversity analysis.
#'
#' @param input_vcf   Path to input VCF or VCF.gz (unphased).
#' @param out_prefix  Output path prefix. Beagle appends \code{.vcf.gz}.
#' @param beagle_jar  Path to \code{beagle.jar}. Default: searched in
#'   \code{dirname(out_prefix)}, then standard locations.
#' @param java_path   Java executable. Default \code{"java"}.
#' @param java_mem_gb Java heap size in GB (e.g. \code{8} sets \code{-Xmx8g}).
#'   Default \code{NULL} (uses JVM default). Increase for large VCFs.
#' @param nthreads    Beagle threads. Default \code{1L}.
#' @param ref_panel   Optional path to phased reference VCF. Default \code{NULL}.
#' @param map_file    Optional genetic map file for more accurate phasing
#'   in structured populations. Default \code{NULL}.
#' @param chrom       Optional chromosome string passed to Beagle
#'   (e.g. \code{"chr1"} or \code{"1"}). Default \code{NULL} (all chromosomes).
#' @param seed        Integer random seed for reproducibility. Default \code{NULL}.
#' @param burnin      Beagle burn-in iterations. Default \code{NULL} (Beagle default).
#' @param iterations  Beagle phasing iterations. Default \code{NULL} (Beagle default).
#' @param window      Beagle window size. Default \code{NULL} (Beagle default).
#' @param overlap     Beagle window overlap. Default \code{NULL} (Beagle default).
#' @param beagle_args Additional Beagle arguments string, space-separated.
#'   Default \code{""}.
#' @param verbose     Logical. Default \code{TRUE}.
#' @return Invisibly returns the path to the phased VCF.gz.
#'   Beagle stdout and stderr are written to \code{out_prefix.log}.
#' @note \code{Sys.which("beagle.jar")} typically fails to find \code{.jar}
#'   files on PATH. Supply \code{beagle_jar} explicitly or place
#'   \code{beagle.jar} next to \code{out_prefix}.
#'   Download: \url{https://faculty.washington.edu/browning/beagle/beagle.html}
#' @references Browning et al. (2018) Am J Hum Genet 103:338-348.
#' @export
phase_with_beagle <- function(
    input_vcf,
    out_prefix,
    beagle_jar  = NULL,
    java_path   = "java",
    java_mem_gb = NULL,
    nthreads    = 1L,
    ref_panel   = NULL,
    map_file    = NULL,
    chrom       = NULL,
    seed        = NULL,
    burnin      = NULL,
    iterations  = NULL,
    window      = NULL,
    overlap     = NULL,
    beagle_args = "",
    verbose     = TRUE
) {
  if (!is.character(input_vcf) || length(input_vcf) != 1L || !file.exists(input_vcf))
    stop("input_vcf not found: ", input_vcf, call. = FALSE)
  if (!grepl("\\.vcf(\\.gz)?$", input_vcf, ignore.case = TRUE))
    stop("input_vcf must be a .vcf or .vcf.gz file.", call. = FALSE)
  if (!is.null(ref_panel) && !file.exists(ref_panel))
    stop("ref_panel not found: ", ref_panel, call. = FALSE)
  if (!is.null(map_file) && !file.exists(map_file))
    stop("map_file not found: ", map_file, call. = FALSE)

  out_dir <- dirname(out_prefix)
  if (!dir.exists(out_dir))
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  if (is.null(beagle_jar)) {
    cands <- c(
      file.path(out_dir, "beagle.jar"),
      Sys.which("beagle.jar"),
      "/usr/local/bin/beagle.jar",
      file.path(Sys.getenv("HOME"), "beagle.jar")
    )
    found <- cands[nzchar(cands) & file.exists(cands)]
    if (!length(found))
      stop(
        "beagle.jar not found.\n",
        "Place beagle.jar in out_dir ('", out_dir, "') or pass beagle_jar explicitly.\n",
        "Download: https://faculty.washington.edu/browning/beagle/beagle.html",
        call. = FALSE
      )
    beagle_jar <- found[1]
  }

  log_file <- paste0(out_prefix, ".log")
  out_vcf  <- paste0(out_prefix, ".vcf.gz")

  args <- character(0)

  if (!is.null(java_mem_gb)) {
    java_mem_gb <- as.numeric(java_mem_gb)
    if (!is.finite(java_mem_gb) || java_mem_gb <= 0)
      stop("java_mem_gb must be a positive number.", call. = FALSE)
    args <- c(args, paste0("-Xmx", format(java_mem_gb, trim = TRUE), "g"))
  }

  args <- c(
    args,
    "-jar", normalizePath(beagle_jar, mustWork = TRUE),
    paste0("gt=",       normalizePath(input_vcf, mustWork = TRUE)),
    paste0("out=",      out_prefix),
    paste0("nthreads=", as.integer(nthreads))
  )

  if (!is.null(ref_panel))
    args <- c(args, paste0("ref=",        normalizePath(ref_panel, mustWork = TRUE)))
  if (!is.null(map_file))
    args <- c(args, paste0("map=",        normalizePath(map_file, mustWork = TRUE)))
  if (!is.null(chrom))
    args <- c(args, paste0("chrom=",      as.character(chrom)))
  if (!is.null(seed))
    args <- c(args, paste0("seed=",       as.integer(seed)))
  if (!is.null(burnin))
    args <- c(args, paste0("burnin=",     as.integer(burnin)))
  if (!is.null(iterations))
    args <- c(args, paste0("iterations=", as.integer(iterations)))
  if (!is.null(window))
    args <- c(args, paste0("window=",     window))
  if (!is.null(overlap))
    args <- c(args, paste0("overlap=",    overlap))

  if (nzchar(beagle_args)) {
    extra_args <- strsplit(beagle_args, "[[:space:]]+")[[1L]]
    args <- c(args, extra_args[nzchar(extra_args)])
  }

  if (verbose) {
    message("[phase_with_beagle] Running Beagle")
    message("[phase_with_beagle]   input_vcf : ", input_vcf)
    message("[phase_with_beagle]   out_prefix: ", out_prefix)
    message("[phase_with_beagle]   beagle_jar: ", beagle_jar)
    if (!is.null(map_file))  message("[phase_with_beagle]   map_file  : ", map_file)
    if (!is.null(ref_panel)) message("[phase_with_beagle]   ref_panel : ", ref_panel)
    if (!is.null(chrom))     message("[phase_with_beagle]   chrom     : ", chrom)
    message("[phase_with_beagle]   threads   : ", as.integer(nthreads))
    message("[phase_with_beagle]   log       : ", log_file)
  }

  status <- system2(command = java_path, args = args,
                    stdout = log_file, stderr = log_file)

  if (!identical(status, 0L)) {
    log_tail <- tryCatch({
      lns <- readLines(log_file, warn = FALSE)
      tail_lns <- tail(lns[nzchar(lns)], 10L)
      if (length(tail_lns)) paste0("\n  Log (last 10 lines):\n  ",
                                   paste(tail_lns, collapse = "\n  "))
      else ""
    }, error = function(e) "")
    stop("Beagle exited with status ", status, ". Log: ", log_file,
         log_tail, call. = FALSE)
  }
  if (!file.exists(out_vcf))
    stop("Beagle output not found: ", out_vcf, call. = FALSE)
  if (!file.exists(paste0(out_vcf, ".tbi")) && verbose)
    message("[phase_with_beagle] Note: phased VCF index (.tbi) not found.")
  if (verbose)
    message("[phase_with_beagle] Done: ", out_vcf)

  invisible(out_vcf)
}

#' Collapse Phased Gametes to 0/1/2 Dosage
#' @param phased_list List from \code{\link{read_phased_vcf}}.
#' @return Numeric matrix (individuals x SNPs, 0/1/2/NA).
#' @keywords internal
unphase_to_dosage <- function(phased_list) {
  if (!all(c("hap1", "hap2") %in% names(phased_list)))
    stop("Need hap1 and hap2 elements.", call. = FALSE)
  # hap1/hap2 are stored as SNPs x individuals; return individuals x SNPs
  t(phased_list$hap1 + phased_list$hap2)
}

# ==============================================================================
# Internal pipeline QC validator
# ==============================================================================

.validate_hap_output <- function(hap_matrix, hap_info, haplotypes = NULL) {
  issues <- character(0)

  # 1. NAs in feature matrix
  n_na <- sum(is.na(hap_matrix))
  if (n_na > 0L)
    issues <- c(issues, paste0("Feature matrix has ", n_na, " NA values"))

  if (!is.null(hap_info) && nrow(hap_info) > 0L) {
    # 2. Duplicated hap_ids
    n_dup <- sum(duplicated(hap_info$hap_id))
    if (n_dup > 0L)
      issues <- c(issues, paste0(n_dup, " duplicated hap_id values in hap_info"))

    # 3. Singleton coordinate / multi-SNP mismatch
    bad <- !is.na(hap_info$start_bp) & !is.na(hap_info$end_bp) &
      !is.na(hap_info$n_snps) &
      hap_info$start_bp == hap_info$end_bp & hap_info$n_snps > 1L
    if (any(bad))
      issues <- c(issues, paste0(sum(bad),
                                 " hap_info rows have start_bp==end_bp but n_snps>1 (retained_idx mismatch)"))

    # 4. hap_id / matrix column alignment
    if (!is.null(hap_matrix)) {
      missing_ids <- setdiff(hap_info$hap_id, colnames(hap_matrix))
      if (length(missing_ids) > 0L)
        issues <- c(issues, paste0(length(missing_ids),
                                   " hap_info hap_id(s) not found as matrix columns"))
    }

    # 5. n_snps range
    nsnps <- hap_info$n_snps[!is.na(hap_info$n_snps)]
    if (length(nsnps)) {
      message("[LDxBlocks] Hap QC: n_snps range [", min(nsnps), ", ", max(nsnps), "]",
              " | blocks=", nrow(hap_info),
              " | NA_matrix=", n_na)
    } else {
      message("[LDxBlocks] Hap QC: n_snps unavailable",
              " | blocks=", nrow(hap_info),
              " | NA_matrix=", n_na)
    }
  }

  if (length(issues) > 0L) {
    warning("[LDxBlocks] Pipeline QC warnings:\n",
            paste0("  - ", issues, collapse="\n"), call.=FALSE)
  } else {
    message("[LDxBlocks] Pipeline QC: all checks passed.")
  }
  invisible(issues)
}

.prefilter_blocks_by_span <- function(chr_blk_srt, snp_pos_sorted, min_snps) {
  # Keep only blocks whose bp interval [start.bp, end.bp] contains >= min_snps
  # positions in snp_pos_sorted (a sorted integer vector of SNP bp positions).
  # Returns a logical vector of length nrow(chr_blk_srt).
  lo <- findInterval(as.integer(chr_blk_srt$start.bp) - 1L, snp_pos_sorted) + 1L
  hi <- findInterval(as.integer(chr_blk_srt$end.bp),         snp_pos_sorted)
  pmax(0L, hi - lo + 1L) >= as.integer(min_snps)
}

#' Extract Haplotype Dosage Strings from LD Blocks
#'
#' @description
#' Builds per-block haplotype dosage strings for all individuals across the
#' LD blocks in \code{blocks}. Each block is processed by the C++ engine
#' \code{extract_chr_haplotypes_cpp()} (unphased) or
#' \code{extract_chr_haplotypes_phased_cpp()} (phased VCF input), which
#' assigns each individual a dosage string of 0/1/2 characters (one per SNP
#' in the block) and identifies the top haplotype alleles by frequency.
#'
#' @param geno One of:
#'   \itemize{
#'     \item An \code{LDxBlocks_backend} from \code{\link{read_geno}} or
#'       \code{\link{read_geno_bigmemory}} (streaming, one chromosome at a time).
#'     \item A named list with elements \code{hap1} and \code{hap2} (phased
#'       SNPs x individuals matrices from \code{\link{read_phased_vcf}}).
#'     \item A numeric matrix (individuals x SNPs, values 0/1/2/NA).
#'   }
#' @param snp_info Data frame with columns \code{SNP}, \code{CHR}, \code{POS}.
#' @param blocks Data frame of LD blocks from \code{\link{run_Big_LD_all_chr}},
#'   with columns \code{CHR}, \code{start.bp}, \code{end.bp}, \code{n_snps}.
#' @param chr Character vector of chromosomes to process. \code{NULL} (default)
#'   processes all chromosomes present in \code{blocks}.
#' @param min_snps Integer. Minimum number of SNPs a block must contain to be
#'   included. Default \code{3L}.
#' @param na_char Character. Symbol used to denote missing genotype in the
#'   dosage string. Default \code{"."}.
#'
#' @return A named list of per-block haplotype dosage matrices (individuals x
#'   haplotype alleles, values 0/1/2 for phased data or 0/1 for unphased).
#'   The list carries a \code{block_info} attribute (data frame with one row
#'   per block: \code{block_id}, \code{CHR}, \code{start_bp}, \code{end_bp},
#'   \code{n_snps}, \code{n_haplotypes}, \code{phased}).
#'
#' @seealso \code{\link{run_Big_LD_all_chr}}, \code{\link{build_haplotype_feature_matrix}},
#'   \code{\link{compute_haplotype_diversity}}, \code{\link{decode_haplotype_strings}}
#'
#' @examples
#' data(ldx_geno, ldx_snp_info, ldx_blocks)
#' haps <- extract_haplotypes(ldx_geno, ldx_snp_info, ldx_blocks, min_snps = 3L)
#' length(haps)                     # one element per block
#' names(haps)[1]                   # e.g. "block_1_1000_25000"
#' dim(haps[[1]])                   # individuals x haplotype alleles
#'
#' @export
extract_haplotypes <- function(geno, snp_info, blocks,
                               chr=NULL, min_snps=3L, na_char=".") {

  # -- Required-column validation -------------------------------------------
  req_snp  <- c("SNP", "CHR", "POS")
  miss_snp <- setdiff(req_snp, names(snp_info))
  if (length(miss_snp))
    stop("snp_info is missing required column(s): ",
         paste(miss_snp, collapse = ", "), call. = FALSE)

  req_blk  <- c("CHR", "start.bp", "end.bp")
  miss_blk <- setdiff(req_blk, names(blocks))
  if (length(miss_blk))
    stop("blocks is missing required column(s): ",
         paste(miss_blk, collapse = ", "), call. = FALSE)

  min_snps <- suppressWarnings(as.integer(min_snps))
  if (length(min_snps) != 1L || is.na(min_snps) || min_snps < 1L)
    stop("min_snps must be a positive integer.", call. = FALSE)

  # -- Input type detection --------------------------------------------------
  # Three paths:
  #  1. LDxBlocks_backend  - STREAMING: one chromosome at a time, never full
  #     genome in RAM. Each chromosome is extracted, all its blocks processed,
  #     then freed with gc(FALSE) before the next chromosome.
  #  2. Phased list (hap1/hap2) - already in RAM, process as before.
  #  3. Numeric dosage matrix - already in RAM, process as before.
  is_backend <- inherits(geno, "LDxBlocks_backend")
  isp <- !is_backend && is.list(geno) && all(c("hap1","hap2")%in%names(geno))

  if (is_backend) {
    # -- STREAMING PATH: chromosome by chromosome ---------------------------
    iids       <- geno$sample_ids
    si         <- geno$snp_info
    si$CHR     <- .norm_chr_hap(si$CHR)
    blocks$CHR <- .norm_chr_hap(blocks$CHR)
    if (!is.null(chr)) {
      chr    <- .norm_chr_hap(chr)
      blocks <- blocks[blocks$CHR %in% chr, ]
    }
    res      <- list()
    bi_rows  <- list()   # accumulate block_info rows; bind once at end
    bi_count <- 0L
    for (cb in unique(blocks$CHR)) {
      chr_blocks <- blocks[blocks$CHR == cb, ]
      if (!nrow(chr_blocks)) next
      chr_idx <- which(si$CHR == cb)
      if (!length(chr_idx)) next

      # Extract ONLY this chromosome - RAM proportional to one chromosome
      chr_geno  <- read_chunk(geno, chr_idx)   # individuals x chr_SNPs
      chr_si    <- si[chr_idx, ]

      # B7+B9+B10: single C++ call handles the entire chromosome:
      #   - interval lookup (no per-block R scanning)
      #   - string building with OpenMP parallel for
      #   - frequency tabulation and min_freq/top_n filtering
      # All returned in one List; no R loop over blocks needed.
      chr_pos_ord <- order(chr_si$POS)
      chr_si_srt  <- chr_si[chr_pos_ord, , drop=FALSE]
      blk_ord     <- order(chr_blocks$start.bp)
      chr_blk_srt <- chr_blocks[blk_ord, , drop=FALSE]

      # PRE-FILTER: remove blocks that cannot pass min_snps before calling C++.
      # This ensures the compact cpp_res arrays map 1:1 to chr_blk_srt rows,
      # eliminating any index mismatch between retained results and block metadata.
      snp_pos_sorted <- as.integer(chr_si_srt$POS)
      keep_blk <- .prefilter_blocks_by_span(chr_blk_srt, snp_pos_sorted, min_snps)
      chr_blk_srt <- chr_blk_srt[keep_blk, , drop=FALSE]
      if (!nrow(chr_blk_srt)) {
        rm(chr_geno, chr_si); gc(FALSE)
        next
      }

      geno_sorted <- chr_geno[, chr_pos_ord, drop=FALSE]
      cpp_res <- extract_chr_haplotypes_cpp(
        geno_chr  = matrix(as.integer(geno_sorted), nrow=nrow(geno_sorted)),
        snp_pos   = snp_pos_sorted,
        block_sb  = as.integer(chr_blk_srt$start.bp),
        block_eb  = as.integer(chr_blk_srt$end.bp),
        min_snps  = as.integer(min_snps),
        min_freq  = 0,
        top_n     = 0L,
        na_char   = na_char
      )
      n_ret   <- cpp_res$n_retained
      ret_idx <- as.integer(cpp_res$retained_idx)

      # Sanity check: retained_idx must exist and be the right length
      if (n_ret > 0L) {
        if (is.null(ret_idx) || length(ret_idx) != n_ret)
          stop("extract_chr_haplotypes_cpp(): retained_idx missing or wrong length. ",
               "Recompile the package with devtools::install().")
        if (any(ret_idx < 1L | ret_idx > nrow(chr_blk_srt)))
          stop("extract_chr_haplotypes_cpp(): retained_idx out of range.")

        # Use C++ bp coordinates directly when available (strongest contract);
        # fall back to R-side index lookup for backward compatibility.
        use_cpp_coords <- !is.null(cpp_res$retained_start_bp) &&
          length(cpp_res$retained_start_bp) == n_ret
        for (b in seq_len(n_ret)) {
          if (is.na(cpp_res$n_snps[b]) || cpp_res$n_snps[b] < min_snps) next

          if (use_cpp_coords) {
            sb <- as.numeric(cpp_res$retained_start_bp[b])
            eb <- as.numeric(cpp_res$retained_end_bp[b])
          } else {
            blk <- chr_blk_srt[ret_idx[b], , drop=FALSE]
            sb  <- as.numeric(blk$start.bp); eb <- as.numeric(blk$end.bp)
          }
          bid <- paste0("block_", cb, "_", sb, "_", eb)
          hs  <- cpp_res$hap_strings[[b]]
          names(hs) <- iids
          res[[bid]] <- hs
          bi_count <- bi_count + 1L
          bi_rows[[bi_count]] <- data.frame(
            block_id = bid, CHR = cb,
            start_bp = as.integer(sb), end_bp = as.integer(eb),
            n_snps   = as.integer(cpp_res$n_snps[b]),
            phased   = FALSE, stringsAsFactors = FALSE
          )
        }
        if (length(res) %% 1000L == 0L && length(res) > 0L)
          message("[extract_haplotypes] ", length(res), " blocks extracted...")
      }
      rm(chr_geno, chr_si, geno_sorted, cpp_res); gc(FALSE)
    }
    bi <- if (bi_count > 0L)
      data.table::rbindlist(bi_rows[seq_len(bi_count)], use.names=TRUE)
    else
      data.frame(block_id=character(),CHR=character(),start_bp=integer(),
                 end_bp=integer(),n_snps=integer(),phased=logical(),
                 stringsAsFactors=FALSE)
    attr(res,"block_info") <- as.data.frame(bi)
    return(res)
  }

  if (isp) {
    # Detect dosage orientation using sample_ids shape match.
    # bigmemory bm[] strips dimnames, so we cannot rely on rownames/colnames.
    # Priority order:
    #   1. Shape: if nrow(dosage) == length(sample_ids) -> [n_ind x n_snps]
    #      (pipeline build path, reattach path - bigmemory returns unnamed matrix)
    #   2. Shape: if ncol(dosage) == length(sample_ids) -> [n_snps x n_ind]
    #      (test/legacy path - t(ldx_geno) where individuals are columns)
    #   3. Dimname probe: rownames match snp_info$SNP -> [n_snps x n_ind]
    #      (legacy paths with named matrices)
    # After orientation, sg (SNP IDs) is filled from dimnames or snp_info$SNP.
    .dos  <- geno$dosage
    .sids <- geno$sample_ids
    .n_sids <- length(.sids)
    # Determine orientation
    .ind_are_rows <- if (.n_sids > 0L && nrow(.dos) == .n_sids) {
      TRUE   # [n_ind x n_snps]: shape match wins
    } else if (.n_sids > 0L && ncol(.dos) == .n_sids) {
      FALSE  # [n_snps x n_ind]: individuals are columns
    } else {
      # Final fallback: check if rownames match SNP IDs
      .snp_probe <- head(snp_info$SNP, 20L)
      !any(rownames(.dos) %in% .snp_probe)  # TRUE = ind are rows, FALSE = snps are rows
    }
    if (.ind_are_rows) {
      # dosage=[n_ind x n_snps]: individuals are rows, SNPs are columns
      iids <- .sids %||% rownames(.dos)
      sg   <- colnames(.dos) %||% snp_info$SNP  # bigmemory loses colnames; use snp_info
      gm   <- .dos
      h1m  <- geno$hap1
      h2m  <- geno$hap2
    } else {
      # dosage=[n_snps x n_ind]: SNPs are rows, individuals are columns
      iids <- .sids %||% colnames(.dos)
      sg   <- rownames(.dos) %||% snp_info$SNP
      gm   <- t(.dos)
      h1m  <- t(geno$hap1)
      h2m  <- t(geno$hap2)
    }
    # Safety checks: catch orientation/cache mismatches early with a clear error.
    if (length(sg) != ncol(gm))
      stop(
        "Internal error: SNP ID vector length does not match genotype columns: ",
        length(sg), " SNP IDs for ", ncol(gm), " genotype columns.",
        call. = FALSE
      )
    if (is.null(iids) || length(iids) != nrow(gm))
      iids <- paste0("ind", seq_len(nrow(gm)))
  } else {
    if (!is.matrix(geno)) geno <- as.matrix(geno)
    gm <- geno; h1m <- h2m <- NULL; iids <- rownames(geno); sg <- colnames(geno)
  }
  # Pre-build SNP name -> column index hash (O(n) once, O(1) per lookup).
  # Replaces match(snp_info$SNP[idx], sg) which rebuilds the hash every block.
  sg_idx <- stats::setNames(seq_along(sg), sg)
  snp_info$CHR <- .norm_chr_hap(snp_info$CHR)
  blocks$CHR   <- .norm_chr_hap(blocks$CHR)
  if (!is.null(chr)) { chr <- .norm_chr_hap(chr); blocks <- blocks[blocks$CHR%in%chr,] }
  res <- list()
  # Pre-filter: drop blocks smaller than min_snps using the index columns.
  # This is a vectorized O(n_blocks) operation that eliminates the expensive
  # per-block which() scan for singletons (which make up >90% of WGS blocks).
  if ("end" %in% names(blocks) && "start" %in% names(blocks)) {
    blocks <- blocks[(blocks$end - blocks$start + 1L) >= min_snps, , drop = FALSE]
  }
  if (!nrow(blocks)) return(structure(list(), block_info = data.frame(
    block_id=character(), CHR=character(), start_bp=integer(),
    end_bp=integer(), n_snps=integer(), phased=logical(),
    stringsAsFactors=FALSE)))

  # Pre-index snp_info by chromosome for O(log n) position lookup.
  # Without this, each block does which(snp_info$CHR==cb & ...) scanning all
  # 2.96M SNPs -- 291k blocks * 2.96M SNPs = 861 billion comparisons.
  snp_info$CHR <- .norm_chr_hap(snp_info$CHR)
  # B6: C++ single-pass interval lookup per chromosome.
  # block_snp_ranges_cpp() does O(p + n_blocks) work per chromosome
  # vs the previous O(p * n_blocks) repeated scanning.
  si_by_chr  <- split(seq_len(nrow(snp_info)), snp_info$CHR)
  blk_by_chr <- split(seq_len(nrow(blocks)),   blocks$CHR)

  bi_rows  <- vector("list", nrow(blocks))
  bi_count <- 0L
  for (cb in intersect(names(blk_by_chr), names(si_by_chr))) {
    chr_si_idx  <- si_by_chr[[cb]]
    chr_blk_idx <- blk_by_chr[[cb]]
    chr_si      <- snp_info[chr_si_idx, , drop=FALSE]
    chr_blks    <- blocks[chr_blk_idx, , drop=FALSE]
    # Index-span pre-filter (vectorized)
    if ("end" %in% names(chr_blks) && "start" %in% names(chr_blks))
      chr_blks <- chr_blks[(chr_blks$end - chr_blks$start + 1L) >= min_snps, , drop=FALSE]
    if (!nrow(chr_blks)) next
    # Sort SNP positions and block starts once per chromosome
    si_ord      <- order(chr_si$POS)
    chr_si_srt  <- chr_si[si_ord, , drop=FALSE]
    si_orig_idx <- chr_si_idx[si_ord]   # original snp_info row indices
    blk_ord     <- order(chr_blks$start.bp)
    chr_blk_srt <- chr_blks[blk_ord, , drop=FALSE]

    # PRE-FILTER: remove blocks that cannot reach min_snps before calling C++.
    # cpp_res arrays are compact (one entry per retained block); pre-filtering
    # ensures they map 1:1 to chr_blk_srt rows without any index mismatch.
    # ------------------------------------------------------------------------
    # Align SNP metadata to genotype columns by SNP ID (friend's clean fix).
    # Works uniformly for both phased and unphased paths.
    # gm is individuals x SNPs, so gm[, idx] is the direct column selection.
    # ------------------------------------------------------------------------
    chr_snps <- as.character(chr_si_srt$SNP)
    idx <- unname(sg_idx[chr_snps])

    valid <- !is.na(idx)

    if (!all(valid)) {
      warning(
        "[extract_haplotypes] Chromosome ", cb, ": dropping ",
        sum(!valid),
        " SNP(s) present in snp_info but absent from genotype matrix.",
        call. = FALSE
      )
    }

    chr_si_srt <- chr_si_srt[valid, , drop = FALSE]
    idx        <- as.integer(idx[valid])

    if (length(idx) < min_snps) next

    snp_pos_sorted_m <- as.integer(chr_si_srt$POS)

    keep_blk_m <- .prefilter_blocks_by_span(chr_blk_srt, snp_pos_sorted_m, min_snps)
    chr_blk_srt <- chr_blk_srt[keep_blk_m, , drop = FALSE]
    if (!nrow(chr_blk_srt)) next

    # Extract aligned chromosome genotype matrix: individuals x SNPs, pos-sorted.
    gm_chr_srt <- gm[, idx, drop = FALSE]

    if (isp) {
      # Phased path: use dedicated C++ phased extractor.
      h1_chr_srt <- h1m[, idx, drop = FALSE]
      h2_chr_srt <- h2m[, idx, drop = FALSE]
      cpp_res <- extract_chr_haplotypes_phased_cpp(
        hap1_chr = matrix(as.integer(h1_chr_srt), nrow=nrow(h1_chr_srt)),
        hap2_chr = matrix(as.integer(h2_chr_srt), nrow=nrow(h2_chr_srt)),
        snp_pos  = snp_pos_sorted_m,
        block_sb = as.integer(chr_blk_srt$start.bp),
        block_eb = as.integer(chr_blk_srt$end.bp),
        min_snps = as.integer(min_snps),
        min_freq = 0, top_n = 0L, na_char = na_char
      )
    } else {
      cpp_res <- extract_chr_haplotypes_cpp(
        geno_chr  = matrix(as.integer(gm_chr_srt), nrow=nrow(gm_chr_srt)),
        snp_pos   = snp_pos_sorted_m,
        block_sb  = as.integer(chr_blk_srt$start.bp),
        block_eb  = as.integer(chr_blk_srt$end.bp),
        min_snps  = as.integer(min_snps),
        min_freq  = 0, top_n = 0L, na_char = na_char
      )
    }
    n_ret   <- cpp_res$n_retained
    ret_idx <- as.integer(cpp_res$retained_idx)

    if (n_ret > 0L) {
      if (is.null(ret_idx) || length(ret_idx) != n_ret)
        stop("extract_chr_haplotypes_cpp(): retained_idx missing or wrong length. ",
             "Recompile the package with devtools::install().")
      if (any(ret_idx < 1L | ret_idx > nrow(chr_blk_srt)))
        stop("extract_chr_haplotypes_cpp(): retained_idx out of range.")
    }

    use_cpp_coords_m <- !is.null(cpp_res$retained_start_bp) &&
      length(cpp_res$retained_start_bp) == n_ret
    for (b in seq_len(n_ret)) {
      if (is.na(cpp_res$n_snps[b]) || cpp_res$n_snps[b] < min_snps) next

      if (use_cpp_coords_m) {
        sb <- as.numeric(cpp_res$retained_start_bp[b])
        eb <- as.numeric(cpp_res$retained_end_bp[b])
      } else {
        blk <- chr_blk_srt[ret_idx[b], , drop=FALSE]
        sb  <- as.numeric(blk$start.bp); eb <- as.numeric(blk$end.bp)
      }

      bid <- paste0("block_", cb, "_", sb, "_", eb)
      # Both phased and unphased: hap_strings come from C++ extractor
      hs <- cpp_res$hap_strings[[b]]
      names(hs) <- iids; res[[bid]] <- hs
      bi_count <- bi_count + 1L
      if (bi_count %% 1000L == 0L)
        message("[extract_haplotypes] ", bi_count, " blocks extracted...")
      bi_rows[[bi_count]] <- data.frame(
        block_id = bid, CHR = cb,
        start_bp = as.integer(sb), end_bp = as.integer(eb),
        n_snps   = as.integer(cpp_res$n_snps[b]),
        phased   = isp, stringsAsFactors = FALSE
      )
    }  # end for(b)
    rm(gm_chr_srt, cpp_res); gc(FALSE)
  }  # end for(cb)
  bi <- if (bi_count > 0L) data.table::rbindlist(bi_rows[seq_len(bi_count)],
                                                 use.names=TRUE) else data.frame(
                                                   block_id=character(),CHR=character(),start_bp=integer(),
                                                   end_bp=integer(),n_snps=integer(),phased=logical(),
                                                   stringsAsFactors=FALSE)
  attr(res,"block_info") <- as.data.frame(bi); res
}

#' Compute Haplotype Diversity Per Block
#'
#' @description
#' Calculates per-block haplotype diversity metrics: richness
#' (n_haplotypes), expected heterozygosity (He, Nei 1973 sample-size
#' corrected), Shannon entropy, effective number of alleles
#' (1/\eqn{\sum p_i^2}), dominant haplotype frequency, and a sweep
#' flag (TRUE when freq_dominant \eqn{\geq} 0.90). These metrics
#' directly correspond to those used to characterise block diversity
#' and identify selection signatures in Difabachew et al. (2023) and
#' Tong et al. (2024). Phased data contributes two gamete observations
#' per individual, doubling the effective sample size.
#'
#' @param haplotypes Named list from \code{\link{extract_haplotypes}}.
#' @param missing_string Missing data marker. Default \code{"."}.
#'
#' @return Data frame with one row per block: \code{block_id},
#'   \code{CHR}, \code{start_bp}, \code{end_bp}, \code{n_snps},
#'   \code{n_ind}, \code{n_haplotypes}, \code{He} (corrected),
#'   \code{Shannon}, \code{n_eff_alleles}, \code{freq_dominant},
#'   \code{sweep_flag}, \code{phased}.
#'
#' @references
#' Nei M (1973). Analysis of gene diversity in subdivided populations.
#' \emph{Proceedings of the National Academy of Sciences}
#' \strong{70}(12):3321-3323. \doi{10.1073/pnas.70.12.3321}
#'
#' Difabachew YF et al. (2023). Genomic prediction with haplotype
#' blocks in wheat. \emph{Frontiers in Plant Science} \strong{14}:1168547.
#' \doi{10.3389/fpls.2023.1168547}
#'
#' Tong J et al. (2024). Stacking beneficial haplotypes from the Vavilov
#' wheat collection to accelerate breeding for multiple disease resistance.
#' \emph{Theoretical and Applied Genetics} \strong{137}:274.
#' \doi{10.1007/s00122-024-04784-w}
#' @export
compute_haplotype_diversity <- function(haplotypes, missing_string=".") {
  bi <- attr(haplotypes,"block_info")
  rows <- lapply(seq_along(haplotypes), function(i) {
    bn <- names(haplotypes)[i]; hap <- haplotypes[[bn]]
    blk <- if(!is.null(bi)) bi[bi$block_id==bn,,drop=FALSE] else NULL
    phased <- isTRUE(blk$phased[1L])
    obs <- if(phased) unlist(strsplit(hap,"|",fixed=TRUE)) else hap
    obs <- obs[!grepl(missing_string,obs,fixed=TRUE)]
    ni  <- if(phased) length(obs)%/%2L else length(obs)
    if (!length(obs)) return(data.frame(block_id=bn,
                                        CHR=if(!is.null(blk))blk$CHR[1] else NA,
                                        start_bp=if(!is.null(blk))blk$start_bp[1] else NA,
                                        end_bp=if(!is.null(blk))blk$end_bp[1] else NA,
                                        n_snps=if(!is.null(blk))blk$n_snps[1] else NA,
                                        n_ind=0L, n_haplotypes=NA, He=NA, Shannon=NA,
                                        n_eff_alleles=NA, freq_dominant=NA, sweep_flag=NA,
                                        phased=phased, stringsAsFactors=FALSE))
    tbl  <- table(obs); freq <- as.numeric(tbl)/sum(tbl)
    He   <- 1 - sum(freq^2)
    # Nei (1973) sample-size correction.
    # The correction must use the number of HAPLOTYPE OBSERVATIONS (n_obs),
    # not the number of individuals (ni). For phased data n_obs = 2*ni
    # (two gametes per individual), which gives the correct correction.
    # Using ni for phased data would under-correct He.
    n_obs   <- length(obs)
    He_corr <- if (n_obs > 1L) (n_obs / (n_obs - 1L)) * He else He
    Sh  <- -sum(freq * log(pmax(freq, .Machine$double.eps)))
    # Effective number of alleles: reciprocal of homozygosity (Hill 1973)
    # Ranges from 1 (monomorphic) to n_haplotypes (equal frequencies)
    n_eff <- 1 / sum(freq^2)
    # Sweep flag: dominant haplotype > 90% suggests selective sweep or
    # strong founder effect in the region (Difabachew et al. 2023)
    sweep <- max(freq) >= 0.90
    data.frame(block_id=bn,
               CHR=if(!is.null(blk))blk$CHR[1] else NA,
               start_bp=if(!is.null(blk))blk$start_bp[1] else NA,
               end_bp=if(!is.null(blk))blk$end_bp[1] else NA,
               n_snps=if(!is.null(blk))blk$n_snps[1] else NA,
               n_ind=ni, n_haplotypes=length(tbl),
               He=He_corr, Shannon=Sh,
               n_eff_alleles=round(n_eff, 3),
               freq_dominant=max(freq),
               sweep_flag=sweep,
               phased=phased, stringsAsFactors=FALSE)
  })
  do.call(rbind,rows)
}

#' Map GWAS Hits to LD Blocks (Post-GWAS QTL Region Definition)
#' @description Maps significant GWAS markers onto LD blocks to define QTL
#'   regions. Blocks with significant markers from multiple traits are flagged
#'   pleiotropic. Implements the approach of Tong et al. (2024).
#' @param gwas_results Data frame with columns \code{SNP} (or \code{Marker}),
#'   \code{CHR}, \code{POS}. Optional columns: \code{P} (p-value),
#'   \code{BETA} (additive effect estimate), \code{trait}.
#' @param blocks LD block data frame from \code{\link{run_Big_LD_all_chr}}.
#' @param snp_info Full SNP metadata data frame.
#' @param p_threshold Significance threshold. Default \code{5e-8}.
#'   \code{NULL} uses all markers regardless of p-value.
#' @param trait_col Trait column name in \code{gwas_results}. Default
#'   \code{"trait"}.
#' @param min_snps Minimum SNPs per block. Default \code{3L}.
#' @param ld_decay Optional \code{LDxBlocks_decay} object from
#'   \code{\link{compute_ld_decay}}, or a data frame with columns
#'   \code{CHR} and \code{decay_dist_bp}. When supplied, candidate
#'   gene windows are extended by the chromosome-specific decay distance
#'   on both sides of each lead SNP position, adding columns
#'   \code{candidate_region_start}, \code{candidate_region_end}, and
#'   \code{candidate_region_size_kb} to the output. Default \code{NULL}.
#'
#' @section Block effect estimation:
#' When multiple SNPs within a block are GWAS-significant, their marginal
#' BETA values are correlated (due to LD) and cannot be summed directly.
#' LDxBlocks returns three complementary columns when \code{BETA} is
#' present in \code{gwas_results}:
#' \describe{
#'   \item{\code{lead_beta}}{BETA of the lead (lowest-p) SNP. The simplest
#'     proxy for the block effect -- assumes the lead SNP fully tags the
#'     causal variant. Underestimates the block effect when multiple
#'     independent signals exist.}
#'   \item{\code{sig_markers}}{Semicolon-separated IDs of all significant SNPs
#'     in the block. Pass these to a conditional/joint analysis tool (e.g.
#'     COJO, SuSiE, finemap) to obtain independent within-block effects.}
#'   \item{\code{sig_betas}}{Semicolon-separated marginal BETA values for
#'     all significant SNPs (same order as \code{sig_markers}). Their absolute
#'     values are an upper bound on the true block effect because they
#'     include LD-induced inflation; use \code{lead_beta} or joint analysis
#'     for calibrated estimates.}
#' }
#' For the biologically correct block-level effect, fit a haplotype model:
#' \code{build_haplotype_feature_matrix()} produces the 0/1/2 dosage columns
#' whose regression coefficients (alpha_h) are the true per-allele effects.
#'
#' @return Data frame with one row per block containing significant markers:
#'   \code{block_id}, \code{CHR}, \code{start_bp}, \code{end_bp},
#'   \code{n_snps_block}, \code{n_sig_markers}, \code{lead_marker},
#'   \code{lead_p}, \code{lead_beta} (if BETA supplied), \code{sig_markers},
#'   \code{sig_betas} (if BETA supplied), \code{traits}, \code{n_traits},
#'   \code{pleiotropic}.
#' @references
#' Tong J, Tarekegn ZT, Jambuthenne D, Alahmad S, Periyannan S,
#' Hickey L, Dinglasan E, Hayes B (2024). Stacking beneficial haplotypes
#' from the Vavilov wheat collection to accelerate breeding for multiple
#' disease resistance. \emph{Theoretical and Applied Genetics}
#' \strong{137}:274. \doi{10.1007/s00122-024-04784-w}
#'
#' Yang J et al. (2012). Conditional and joint multiple-SNP analysis of GWAS
#' summary statistics identifies additional variants influencing complex traits.
#' \emph{Nature Genetics} \strong{44}(4):369-375. \doi{10.1038/ng.2213}
#' @param verbose Logical. Print progress messages. Default \code{TRUE}.
#' @export
define_qtl_regions <- function(gwas_results, blocks, snp_info,
                               p_threshold = 5e-8, trait_col = "trait",
                               min_snps = 3L, ld_decay = NULL,
                               verbose = TRUE) {
  # Accept "Marker" as an alias for "SNP" (OptSLDP / GWAS convention)
  if (!"SNP" %in% names(gwas_results) && "Marker" %in% names(gwas_results))
    gwas_results$SNP <- gwas_results$Marker
  # Detect whether markers are LDxBlocks block IDs or external GWAS SNP IDs.
  # Block IDs always start with "block_"; rsIDs / SNP names never do.
  marker_source <- if (any(grepl("^block_", gwas_results$SNP, perl = TRUE)))
    "block_id" else "snp_id"
  miss <- setdiff(c("SNP","CHR","POS"),names(gwas_results))
  if (length(miss)) stop("gwas_results missing: ",paste(miss,collapse=","),call.=FALSE)
  gwas_results$CHR <- .norm_chr_hap(gwas_results$CHR)
  blocks$CHR <- .norm_chr_hap(blocks$CHR); snp_info$CHR <- .norm_chr_hap(snp_info$CHR)
  if (!is.null(p_threshold)&&"P"%in%names(gwas_results))
    gwas_results <- gwas_results[!is.na(gwas_results$P)&gwas_results$P<=p_threshold,]
  if (!nrow(gwas_results)){message("[define_qtl_regions] No significant markers.");return(data.frame())}
  if (!trait_col%in%names(gwas_results)) gwas_results[[trait_col]] <- "trait"
  rows <- list()
  # Prepare ld_decay lookup: named vector CHR -> decay_dist_bp
  # Accepts: LDxBlocks_decay object (from compute_ld_decay()) or
  #          data.frame with columns CHR and decay_dist_bp (or decay_bp).
  decay_map <- NULL
  if (!is.null(ld_decay)) {
    # Unwrap LDxBlocks_decay object
    if (inherits(ld_decay, "LDxBlocks_decay")) {
      if (is.null(ld_decay$decay_dist))
        stop("ld_decay$decay_dist is NULL -- rerun compute_ld_decay() ",
             "with a non-NULL r2_threshold.", call. = FALSE)
      ld_decay <- ld_decay$decay_dist
    }
    # Accept decay_dist_bp or legacy decay_bp column name
    bp_col <- if ("decay_dist_bp" %in% names(ld_decay)) "decay_dist_bp"
    else if ("decay_bp" %in% names(ld_decay)) "decay_bp"
    else stop("ld_decay must have a 'decay_dist_bp' column ",
              "(output of compute_ld_decay()).", call. = FALSE)
    if (!"CHR" %in% names(ld_decay))
      stop("ld_decay must have a 'CHR' column.", call. = FALSE)
    decay_map <- stats::setNames(ld_decay[[bp_col]], ld_decay$CHR)
    # Add genome-wide median as fallback for chromosomes not in decay_dist
    if (!"GENOME" %in% names(decay_map))
      decay_map["GENOME"] <- stats::median(decay_map, na.rm = TRUE)
  }

  for (b in seq_len(nrow(blocks))) {
    blk <- blocks[b,]; ch <- as.character(blk$CHR)
    sb <- as.numeric(blk$start.bp); eb <- as.numeric(blk$end.bp)

    # Extend search window by LD decay distance if supplied.
    # QTL window = block boundaries extended by the chromosome-specific
    # LD decay distance (or genome-wide if that chromosome is absent).
    # This captures candidate genes in LD with the GWAS hit that fall
    # outside the block boundary itself.
    if (!is.null(decay_map)) {
      ext <- if (ch %in% names(decay_map)) decay_map[[ch]]
      else if ("GENOME" %in% names(decay_map)) decay_map[["GENOME"]]
      else 0L
      sb_ext <- max(0, sb - ext)
      eb_ext <- eb + ext
    } else {
      sb_ext <- sb; eb_ext <- eb
    }

    hits <- gwas_results[
      gwas_results$CHR == ch &
        gwas_results$POS >= sb_ext &
        gwas_results$POS <= eb_ext, , drop = FALSE]
    if (!nrow(hits)) next
    nb <- sum(snp_info$CHR == ch & snp_info$POS >= sb & snp_info$POS <= eb)
    if (nb < min_snps) next
    traits <- unique(hits[[trait_col]])
    li <- if("P"%in%names(hits)) which.min(hits$P) else 1L
    has_beta <- "BETA" %in% names(hits)
    rows[[length(rows)+1L]] <- data.frame(
      block_id      = paste0("block_",ch,"_",sb,"_",eb),
      CHR           = ch,
      start_bp      = as.integer(sb),
      end_bp        = as.integer(eb),
      search_start  = as.integer(if (!is.null(decay_map)) sb_ext else sb),
      search_end    = as.integer(if (!is.null(decay_map)) eb_ext else eb),
      ld_decay_bp   = if (!is.null(decay_map)) as.integer(
        if (ch %in% names(decay_map)) decay_map[[ch]]
        else if ("GENOME" %in% names(decay_map)) decay_map[["GENOME"]]
        else 0L) else NA_integer_,
      n_snps_block  = nb,
      n_sig_markers = nrow(hits),
      lead_marker   = hits$SNP[li],
      lead_p        = if ("P" %in% names(hits)) hits$P[li] else NA_real_,
      candidate_region_start = if (!is.null(decay_map)) {
        as.integer(max(0, hits$POS[li] - (
          if (ch %in% names(decay_map)) decay_map[[ch]]
          else decay_map[["GENOME"]])))
      } else as.integer(sb),
      candidate_region_end   = if (!is.null(decay_map)) {
        as.integer(hits$POS[li] + (
          if (ch %in% names(decay_map)) decay_map[[ch]]
          else decay_map[["GENOME"]]))
      } else as.integer(eb),
      candidate_region_size_kb = if (!is.null(decay_map)) {
        round(2 * (if (ch %in% names(decay_map)) decay_map[[ch]]
                   else decay_map[["GENOME"]]) / 1000, 1)
      } else round((eb - sb) / 1000, 1),
      lead_beta     = if (has_beta) hits$BETA[li] else NA_real_,
      sig_markers   = paste(hits$SNP, collapse = ";"),
      sig_betas     = if (has_beta) paste(round(hits$BETA, 6), collapse = ";") else NA_character_,
      traits        = paste(sort(traits), collapse = ","),
      n_traits      = length(traits),
      pleiotropic   = length(traits) > 1L,
      marker_source = marker_source,
      stringsAsFactors = FALSE)
  }
  if (!length(rows)){message("[define_qtl_regions] No overlapping blocks.");return(data.frame())}
  out <- do.call(rbind, rows)

  # When marker_source = "block_id" (LDxBlocks haplotype input):
  # - lead_marker  = the significant haplotype block with lowest p-value
  #                  whose search window overlaps this block. When ld_decay=NULL
  #                  this equals block_id itself (the block IS the hit).
  # - sig_markers  = semicolon-separated block IDs of all significant blocks
  #                  whose search window overlaps this block.
  # - n_sig_markers = count of those significant blocks (not individual SNPs).
  # When marker_source = "snp_id" (external GWAS):
  # - lead_marker, sig_markers, n_sig_markers refer to individual SNPs.
  # The marker_source column always indicates which interpretation applies.
  if (all(out$marker_source == "block_id")) {
    # For clarity, expose n_sig_markers under a block-specific alias too.
    # Both columns are identical in value; n_sig_markers retained for API compat.
    out$n_sig_blocks_in_window <- out$n_sig_markers
    # Flag rows where block_id == lead_marker (block is its own lead hit).
    # True when ld_decay=NULL or when no other sig block falls in window.
    out$is_primary_hit <- out$block_id == out$lead_marker
  }

  out[order(out$CHR, out$start_bp), ]
}

#' Build Haplotype Dosage Matrix for Genomic Prediction
#' @description Converts haplotype strings to a numeric matrix for genomic
#'   prediction. Supports phased and unphased input with two encoding schemes.
#'
#'   \code{encoding="additive_012"} (default, recommended for GBLUP/rrBLUP/BGLR):
#'     Phased:   0=0 copies, 1=1 copy (het), 2=2 copies (hom)
#'     Unphased: 0=no match, 2=match (1 not identifiable without phase)
#'
#'   encoding="presence_02" (kernel methods, random forest):
#'     Phased:   2=either gamete matches, 0=neither, NA=missing
#'     Unphased: 2=match, 0=no match, NA=missing
#'
#' @param haplotypes List from extract_haplotypes().
#' @param top_n Integer or \code{NULL}. Maximum number of haplotype alleles
#'   to retain per block, ranked by frequency. \code{NULL} (default) retains
#'   all alleles that pass \code{min_freq} -- recommended for most analyses.
#'   Set an integer cap (e.g. \code{top_n = 5L}) only when you need to limit
#'   matrix width for memory reasons on panels with thousands of blocks and
#'   highly diverse haplotypes (many rare alleles above \code{min_freq}).
#' @param encoding     Dosage encoding for the feature matrix:
#'   \itemize{
#'     \item \code{"additive_012"} (default):
#'       \strong{Phased data}: values are 0 (neither gamete carries the allele),
#'       1 (one gamete carries it - heterozygous), or 2 (both gametes - homozygous).
#'       This gives true allele dosage on the standard 0/1/2 scale.
#'       \strong{Unphased data}: values are 0 or 1 only (presence/absence of the
#'       dosage-pattern haplotype). The value 2 is never produced for unphased
#'       data because it is impossible to confirm that both chromosomes carry
#'       the same haplotype string without phase information.
#'       Compatible with rrBLUP, BGLR, sommer, ASReml-R.
#'     \item \code{"presence_01"}: values are 0 or 1 for both phased and unphased
#'       data. For phased data: 1 if either gamete carries the allele (loses
#'       copy-number information vs \code{"additive_012"}). For unphased data:
#'       identical to \code{"additive_012"} since that already gives 0/1.
#'       May be preferable for Bayesian variable selection models (BayesB,
#'       BayesC) where the prior expects binary indicators.
#'   }
#'
#' @section Phased vs unphased haplotypes:
#' \strong{Phased data} (from \code{\link{read_phased_vcf}},
#' \code{\link{phase_with_beagle}}):
#' each individual's block string is \code{"g1|g2"} where \code{g1} and \code{g2}
#' are the two gametic sequences. Haplotype alleles are identified at the
#' \emph{gamete} level - true haplotypes in the biological sense. Frequencies
#' are gamete frequencies (each individual contributes two observations).
#' Dosage values 0, 1, 2 measure actual allele copy number.
#'
#' \strong{Unphased data} (from an unphased VCF or dosage matrix): each
#' individual's block string is a single multi-SNP dosage string (e.g.
#' \code{"021002"}). Haplotype \emph{alleles} are distinct dosage patterns,
#' not true gametic haplotypes. Frequencies measure the proportion of
#' \emph{individuals} carrying each dosage pattern. This is biologically
#' meaningful (distinct genotypic patterns at the block level) but should not
#' be equated with true haplotype allele frequency without phasing.
#' The recommended workflow for true haplotype analysis is:
#' phase first with Beagle, then call
#' \code{extract_haplotypes()} on the phased output.
#' @param missing_string Missing data marker. Default ".".
#' @param scale_features Center and scale columns. Default FALSE.
#' @param min_freq Minimum allele frequency to include. Default 0.01.
#' @return Numeric matrix (individuals x haplotype allele columns).
#' @references
#' Difabachew YF et al. (2023). Genomic prediction with haplotype
#' blocks in wheat. \emph{Frontiers in Plant Science} \strong{14}:1168547.
#' \doi{10.3389/fpls.2023.1168547}
#'
#' Weber SE, Frisch M, Snowdon RJ, Voss-Fels KP (2023). Haplotype
#' blocks for genomic prediction: a comparative evaluation in multiple
#' crop datasets. \emph{Frontiers in Plant Science} \strong{14}:1217589.
#' \doi{10.3389/fpls.2023.1217589}
#' @export
build_haplotype_feature_matrix <- function(haplotypes, top_n=NULL,
                                           encoding=c("additive_012","presence_01"),
                                           missing_string=".", scale_features=FALSE,
                                           min_freq=0.01) {
  encoding <- match.arg(encoding)
  if (identical(encoding, "presence_02")) encoding <- "presence_01"

  # Guard: empty haplotypes list (all blocks filtered by min_snps etc.)
  if (!length(haplotypes)) {
    info_empty <- data.frame(
      hap_id = character(), block_id = character(),
      CHR = character(), start_bp = integer(), end_bp = integer(),
      n_snps = integer(), hap_rank = integer(),
      hap_string = character(), frequency = numeric(),
      phased = logical(), alleles_export = character(),
      stringsAsFactors = FALSE
    )
    return(list(
      matrix = matrix(numeric(0), nrow = 0L, ncol = 0L),
      info   = info_empty
    ))
  }

  bi        <- attr(haplotypes, "block_info")
  inames    <- names(haplotypes[[1L]])
  mats      <- vector("list", length(haplotypes))
  info_rows <- list()
  info_i    <- 0L

  for (bk in seq_along(haplotypes)) {
    bn  <- names(haplotypes)[bk]
    hap <- haplotypes[[bk]]

    brow   <- if (!is.null(bi)) bi[bi$block_id == bn, , drop = FALSE] else data.frame()
    phased <- isTRUE(brow$phased[1L])

    if (phased) {
      parts <- strsplit(hap, "|", fixed = TRUE)
      g1 <- vapply(parts, `[`, character(1L), 1L)
      g2 <- vapply(parts, `[`, character(1L), 2L)
      miss <- grepl(missing_string, g1, fixed = TRUE) |
        grepl(missing_string, g2, fixed = TRUE)
      gam  <- c(g1[!grepl(missing_string, g1, fixed = TRUE)],
                g2[!grepl(missing_string, g2, fixed = TRUE)])
    } else {
      miss <- !nzchar(gsub(missing_string, "", hap, fixed = TRUE))
      gam  <- hap[!miss]
    }

    if (!length(gam)) { mats[[bk]] <- NULL; next }

    tbl  <- sort(table(gam), decreasing = TRUE)
    fv   <- as.numeric(tbl) / sum(tbl)
    tbl  <- tbl[fv >= min_freq]
    fv   <- fv[fv >= min_freq]
    tops <- if (is.null(top_n)) names(tbl) else
      names(tbl)[seq_len(min(as.integer(top_n), length(tbl)))]
    if (!length(tops)) { mats[[bk]] <- NULL; next }

    mat <- matrix(NA_real_, nrow = length(inames), ncol = length(tops),
                  dimnames = list(inames, paste0(bn, "_hap", seq_along(tops))))

    for (j in seq_along(tops)) {
      ref <- tops[j]
      if (phased) {
        parts <- strsplit(hap, "|", fixed = TRUE)
        g1 <- vapply(parts, `[`, character(1L), 1L)
        g2 <- vapply(parts, `[`, character(1L), 2L)
        miss <- grepl(missing_string, g1, fixed = TRUE) |
          grepl(missing_string, g2, fixed = TRUE)
        if (identical(encoding, "additive_012")) {
          mat[, j] <- ifelse(miss, NA_real_,
                             as.numeric(g1 == ref) + as.numeric(g2 == ref))
        } else {
          mat[, j] <- ifelse(miss, NA_real_,
                             as.numeric((g1 == ref) | (g2 == ref)))
        }
      } else {
        miss <- !nzchar(gsub(missing_string, "", hap, fixed = TRUE))
        mat[, j] <- ifelse(miss, NA_real_, as.numeric(hap == ref))
      }

      # Accumulate exact per-column metadata while we have the information.
      # This is the only place where the retained hap_string, rank, and
      # frequency are all known simultaneously. The writer uses this
      # directly -- no reconstruction needed later.
      info_i <- info_i + 1L
      info_rows[[info_i]] <- data.frame(
        hap_id        = colnames(mat)[j],
        block_id      = bn,
        CHR           = if (nrow(brow)) brow$CHR[1L]        else NA_character_,
        start_bp      = if (nrow(brow)) as.integer(brow$start_bp[1L]) else NA_integer_,
        end_bp        = if (nrow(brow)) as.integer(brow$end_bp[1L])   else NA_integer_,
        n_snps        = if (nrow(brow)) as.integer(brow$n_snps[1L])   else NA_integer_,
        hap_rank      = j,
        hap_string    = ref,
        frequency     = round(as.numeric(tbl[ref]) / sum(tbl), 4),
        phased        = phased,
        alleles_export = NA_character_,   # filled below after snp_info decode
        stringsAsFactors = FALSE
      )
    }
    mats[[bk]] <- mat
  }

  mats <- Filter(Negate(is.null), mats)

  # Build info data frame (empty case handled explicitly)
  info_empty <- data.frame(
    hap_id=character(), block_id=character(), CHR=character(),
    start_bp=integer(), end_bp=integer(), n_snps=integer(),
    hap_rank=integer(), hap_string=character(), frequency=numeric(),
    phased=logical(), stringsAsFactors=FALSE)

  if (!length(mats)) {
    out <- matrix(NA_real_, nrow=length(inames), ncol=0L,
                  dimnames=list(inames, character(0)))
    return(list(matrix=out, info=info_empty))
  }

  out <- do.call(cbind, mats)
  if (scale_features && ncol(out) > 0L) { out <- scale(out); out[is.nan(out)] <- 0 }

  n_pred <- ncol(out); n_blks <- length(mats)
  if (n_pred < 500L && n_blks >= 20L)
    warning("build_haplotype_feature_matrix: only ", n_pred,
            " predictor columns from ", n_blks,
            " blocks. Difabachew et al. (2023) show that < 500 haplotype ",
            "predictors reduce genomic prediction accuracy below the ",
            "single-SNP baseline. Consider lowering CLQcut or min_freq.",
            call. = FALSE)

  info <- if (info_i > 0L)
    as.data.frame(data.table::rbindlist(info_rows[seq_len(info_i)],
                                        use.names = TRUE, fill = TRUE))
  else info_empty

  list(matrix = out, info = info)
}


#' Compute Haplotype-Based Genomic Relationship Matrix
#'
#' @description
#' Computes the additive genomic relationship matrix (GRM) from a haplotype
#' feature matrix using the VanRaden (2008) method extended to multi-allelic
#' haplotype blocks. The resulting G matrix can be used directly in GBLUP,
#' rrBLUP, ASReml-R, or any software that accepts a realized relationship
#' matrix.
#'
#' @details
#' The GRM is computed as:
#' \deqn{G = \frac{ZZ^\top}{2\sum_j p_j(1-p_j)}}
#' where \eqn{Z} is the centred haplotype dosage matrix (columns centred by
#' \eqn{2p_j}) and \eqn{p_j} is the frequency of haplotype allele \eqn{j}.
#' This matches the standard G matrix of Weber et al. (2023) and Difabachew
#' et al. (2023), ensuring the relationship scale is compatible with
#' conventional SNP-based GRMs.
#'
#' Missing dosage values (\code{NA}) are mean-imputed per column before
#' centering.
#'
#' @param hap_matrix Numeric matrix (individuals x haplotype alleles) from
#'   \code{\link{build_haplotype_feature_matrix}}.
#' @param bend Logical. If \code{TRUE}, add a small constant to the diagonal
#'   to ensure positive-definiteness: \code{diag(G) + 0.001}. Useful when
#'   passing G to mixed model solvers. Default \code{FALSE}.
#' @param phased Logical or \code{NULL} (default). Whether the haplotype
#'   feature matrix uses \code{additive_012} encoding for phased data
#'   (dosage values 0/1/2, dose scale = 2) or unphased data (presence/absence
#'   0/1, dose scale = 1). \code{NULL} auto-detects from column maxima:
#'   any column with a value > 1 implies phased encoding. For guaranteed
#'   correctness when all haplotype alleles lack homozygous carriers (max=1
#'   even in phased data), supply \code{phased = TRUE} or \code{FALSE}
#'   explicitly.
#'
#' @return Symmetric n x n numeric matrix (individuals x individuals).
#'   Row and column names match \code{rownames(hap_matrix)}.
#'
#' @references
#' VanRaden PM (2008). Efficient methods to compute genomic predictions.
#' \emph{Journal of Dairy Science} \strong{91}(11):4414-4423.
#' \doi{10.3168/jds.2007-0980}
#'
#' Weber SE et al. (2023). Haplotype blocks for genomic prediction: a
#' comparative evaluation in multiple crop datasets.
#' \emph{Frontiers in Plant Science} \strong{14}:1217589.
#' \doi{10.3389/fpls.2023.1217589}
#'
#' @examples
#' data(ldx_geno, ldx_snp_info, ldx_blocks, package = "LDxBlocks")
#' haps <- extract_haplotypes(ldx_geno, ldx_snp_info, ldx_blocks, min_snps = 3)
#' feat <- build_haplotype_feature_matrix(haps, top_n = 5)$matrix
#' G    <- compute_haplotype_grm(feat)
#' dim(G)
#' round(range(diag(G)), 3)  # diagonal ~= 1 for typical populations
#'
#' @export
compute_haplotype_grm <- function(hap_matrix, bend = FALSE, phased = NULL) {
  if (!is.matrix(hap_matrix))
    hap_matrix <- as.matrix(hap_matrix)

  if (!ncol(hap_matrix))
    stop("hap_matrix has zero columns; cannot compute GRM.", call. = FALSE)

  # Drop all-NA columns (mean imputation would produce NaN)
  all_na <- colSums(!is.na(hap_matrix)) == 0L
  if (any(all_na)) hap_matrix <- hap_matrix[, !all_na, drop = FALSE]
  if (!ncol(hap_matrix))
    stop("All haplotype columns are NA; cannot compute GRM.", call. = FALSE)

  # Drop monomorphic columns (zero variance -> VanRaden scaling undefined)
  mono <- apply(hap_matrix, 2L, function(x) {
    ux <- unique(x[!is.na(x)]); length(ux) <= 1L
  })
  if (any(mono)) hap_matrix <- hap_matrix[, !mono, drop = FALSE]
  if (!ncol(hap_matrix))
    stop("All haplotype columns are monomorphic; cannot compute GRM.",
         call. = FALSE)

  # Mean-impute NA per column
  for (j in seq_len(ncol(hap_matrix))) {
    na_j <- is.na(hap_matrix[, j])
    if (any(na_j))
      hap_matrix[na_j, j] <- mean(hap_matrix[, j], na.rm = TRUE)
  }

  # Detect encoding scale:
  #   Phased data   -> additive_012 gives 0/1/2 dosage -> dose_scale = 2
  #   Unphased data -> additive_012 gives 0/1 presence -> dose_scale = 1
  #
  # VanRaden (2008) centres by 2p and scales by 2*sum(p*(1-p)), where p is
  # the allele frequency on the [0,1] scale. For 0/1/2 dosage p = mean/2;
  # for 0/1 presence/absence p = mean (no division needed).
  # Using the wrong dose_scale produces a GRM scaled by a constant factor,
  # which does not affect eigenvectors (PCs) but mis-scales relationship values.
  #
  # phased=NULL: auto-detect from column maxima.
  #   max > 1 in any column -> phased 0/1/2 encoding -> dose_scale = 2
  #   max <= 1 everywhere   -> unphased 0/1 encoding  -> dose_scale = 1
  #   Caveat: phased alleles with NO homozygous carriers also have max=1.
  #   Prefer supplying phased=TRUE/FALSE explicitly for guaranteed correctness.
  if (is.null(phased)) {
    col_max    <- apply(hap_matrix, 2, max, na.rm = TRUE)
    dose_scale <- if (any(col_max > 1 + .Machine$double.eps)) 2 else 1
  } else {
    dose_scale <- if (isTRUE(phased)) 2L else 1L
  }

  # Allele frequencies on [0,1] scale
  p <- colMeans(hap_matrix) / dose_scale
  p <- pmax(pmin(p, 1 - 1e-6), 1e-6)  # clamp for numerical safety

  # Centre: Z_ij = x_ij - dose_scale * p_j
  # (= x_ij - E[x_ij] for both encodings since E[x_j] = dose_scale * p_j)
  Z <- sweep(hap_matrix, 2, dose_scale * p, "-")

  # Scaling denominator: 2 * sum(p*(1-p))  [VanRaden 2008, eq. 3]
  # Correct for both encodings because p is on [0,1] scale in both cases.
  denom <- 2 * sum(p * (1 - p))
  if (denom < .Machine$double.eps)
    stop("All haplotype alleles are monomorphic. Cannot compute GRM.",
         call. = FALSE)

  G <- tcrossprod(Z) / denom

  if (bend) diag(G) <- diag(G) + 0.001

  dimnames(G) <- list(rownames(hap_matrix), rownames(hap_matrix))
  G
}

#' Write Haplotype Character (Nucleotide) Matrix
#'
#' Writes a matrix where each cell contains the nucleotide sequence of the
#' haplotype allele carried by each individual. Rows are haplotype alleles,
#' columns are individuals. This is the most interpretable format: you can
#' read directly which nucleotides define each haplotype allele and which
#' individuals carry it.
#'
#' @details
#' The cell value for individual i at haplotype allele h is:
#' \itemize{
#'   \item The nucleotide sequence (e.g. \code{"AGTTA"}) if the individual
#'     carries that allele (dosage = 2 for unphased, or present in either
#'     gamete for phased).
#'   \item \code{"-"} if the individual does not carry that allele.
#'   \item \code{"."} if the individual has missing data in that block.
#' }
#'
#' Heterozygous positions (dosage = 1, phased data only) are encoded using
#' IUPAC ambiguity codes: R=A/G, Y=C/T, S=G/C, W=A/T, K=G/T, M=A/C.
#' This keeps the nucleotide string the same length as \code{n_snps}
#' regardless of how many heterozygous positions are present.
#'
#' @param haplotypes List from \code{\link{extract_haplotypes}}.
#' @param snp_info   Data frame with \code{SNP}, \code{CHR}, \code{POS},
#'   \code{REF}, \code{ALT}.
#' @param out_file   Output file path (tab-delimited).
#' @param min_freq   Minimum haplotype frequency. Default \code{0.01}.
#' @param top_n      Integer or \code{NULL}. Cap alleles per block.
#'   \code{NULL} (default) keeps all above \code{min_freq}.
#' @param missing_string Missing genotype marker. Default \code{"."}.
#' @param verbose    Logical. Default \code{TRUE}.
#' @return Invisibly returns \code{out_file}.
#' @export
write_haplotype_character <- function(haplotypes, snp_info, out_file,
                                      min_freq       = 0.01,
                                      top_n          = NULL,
                                      missing_string = ".",
                                      verbose        = TRUE) {

  if (!all(c("CHR","POS","REF","ALT") %in% names(snp_info)))
    stop("snp_info must have columns CHR, POS, REF, ALT.", call. = FALSE)

  snp_info$CHR <- .norm_chr_hap(snp_info$CHR)
  bi           <- attr(haplotypes, "block_info")

  # Decode one dosage string to nucleotide sequence
  decode_str <- function(dstr, ref_v, alt_v) {
    chars <- strsplit(dstr, "", fixed = TRUE)[[1L]]
    n     <- min(length(chars), length(ref_v))
    vapply(seq_len(n), function(i) {
      switch(chars[i],
             "0" = ref_v[i],
             "2" = alt_v[i],
             "1" = {
               key <- paste0(ref_v[i], alt_v[i])
               iupac <- c(AG="R",GA="R",CT="Y",TC="Y",GC="S",CG="S",
                          AT="W",TA="W",GT="K",TG="K",AC="M",CA="M")
               if (!is.na(iupac[key])) unname(iupac[key]) else "N"
             },
             "N")
    }, character(1L)) |> paste(collapse = "")
  }

  rows_all <- vector("list", length(haplotypes))

  for (bk in seq_along(haplotypes)) {
    bn   <- names(haplotypes)[bk]
    hap  <- haplotypes[[bn]]
    brow <- if (!is.null(bi)) bi[bi$block_id == bn, , drop = FALSE] else NULL

    # Parse CHR, start_bp from block name
    parts_bn <- strsplit(bn, "_", fixed = TRUE)[[1L]]
    chr  <- if (!is.null(brow) && nrow(brow)) brow$CHR[1L] else parts_bn[2L]
    sb   <- if (!is.null(brow) && nrow(brow)) brow$start_bp[1L] else
      suppressWarnings(as.integer(parts_bn[3L]))
    eb   <- if (!is.null(brow) && nrow(brow)) brow$end_bp[1L] else
      suppressWarnings(as.integer(parts_bn[4L]))

    # SNPs in block
    blk_snps <- snp_info[snp_info$CHR == chr &
                           snp_info$POS >= sb &
                           snp_info$POS <= eb, , drop = FALSE]
    blk_snps <- blk_snps[order(blk_snps$POS), ]
    if (!nrow(blk_snps)) next

    ref_v <- as.character(blk_snps$REF)
    alt_v <- as.character(blk_snps$ALT)

    # Alleles string for metadata column: REF1/ALT1;REF2/ALT2;...
    alleles_str <- paste(paste0(ref_v, "/", alt_v), collapse = ";")

    # Frequency table - must split phased "g1|g2" strings into gametes before
    # tabulating. Without splitting, phased entries like "010|111" would be
    # treated as a single string and ranked as diplotypes, not haplotype alleles.
    brow_w    <- if (!is.null(bi)) bi[bi$block_id == bn, , drop = FALSE] else NULL
    is_phased_w <- isTRUE(brow_w$phased[1L])
    if (is_phased_w) {
      parts_w <- strsplit(hap, "|", fixed = TRUE)
      g1_w <- vapply(parts_w, `[`, character(1L), 1L)
      g2_w <- vapply(parts_w, `[`, character(1L), 2L)
      miss <- grepl(missing_string, g1_w, fixed = TRUE) |
        grepl(missing_string, g2_w, fixed = TRUE)
      gam  <- c(g1_w[!grepl(missing_string, g1_w, fixed = TRUE)],
                g2_w[!grepl(missing_string, g2_w, fixed = TRUE)])
    } else {
      miss <- grepl(missing_string, hap, fixed = TRUE)
      gam  <- hap[!miss]
    }
    if (!length(gam)) next
    tbl  <- sort(table(gam), decreasing = TRUE)
    fv   <- as.numeric(tbl) / sum(tbl)
    tbl  <- tbl[fv >= min_freq]
    tops <- if (is.null(top_n)) names(tbl)
    else names(tbl)[seq_len(min(as.integer(top_n), length(tbl)))]
    if (!length(tops)) next

    # For each top haplotype allele, build a row:
    # cols = hap_id | CHR | start_bp | end_bp | n_snps | Alleles | ind1 | ind2 | ...
    # cell = nucleotide_sequence if individual carries it, "-" if not, "." if missing
    for (r in seq_along(tops)) {
      dstr   <- tops[r]
      nuc_seq <- decode_str(dstr, ref_v, alt_v)
      hap_id  <- paste0(bn, "_hap", r)

      cell_vals <- ifelse(
        miss, missing_string,
        ifelse(hap == dstr, nuc_seq, "-")
      )

      # Build row: metadata cols + one col per individual
      ind_df <- as.data.frame(
        matrix(cell_vals, nrow = 1L,
               dimnames = list(NULL, names(hap))),
        check.names = FALSE, stringsAsFactors = FALSE
      )
      row_df <- data.frame(
        hap_id   = hap_id,
        CHR      = chr,
        start_bp = as.integer(sb),
        end_bp   = as.integer(eb),
        n_snps   = nrow(blk_snps),
        Alleles  = alleles_str,
        ind_df,
        check.names = FALSE, stringsAsFactors = FALSE,
        row.names = NULL
      )
      rows_all[[length(rows_all) + 1L]] <- row_df
    }
  }

  if (!length(rows_all)) {
    message("[write_haplotype_character] No haplotypes passed filters. File not written.")
    return(invisible(out_file))
  }

  out <- do.call(rbind, rows_all)
  data.table::fwrite(out, out_file, sep = "\t", quote = FALSE, na = ".")
  if (verbose) message("[write_haplotype_character] ", out_file,
                       " (", nrow(out), " haplotypes x ",
                       ncol(out) - 6L, " individuals)")
  invisible(out_file)
}

#' Write Haplotype Feature Matrix as Numeric Dosage Table
#'
#' @description
#' Writes the haplotype dosage matrix in a tab-delimited format with
#' haplotype alleles as rows and individuals as columns. Metadata columns
#' (\code{hap_id}, \code{CHR}, \code{start_bp}, \code{end_bp},
#' \code{n_snps}, \code{alleles}, \code{frequency}) precede the individual
#' columns. Individual cells contain 0/1/2/NA dosage values.
#'
#' @param hap_matrix Numeric matrix (individuals x haplotype alleles) from
#'   \code{\link{build_haplotype_feature_matrix}}.
#' @param out_file   Output file path.
#' @param haplotypes List from \code{\link{extract_haplotypes}}. When supplied
#'   together with \code{snp_info}, the \code{alleles} and \code{frequency}
#'   metadata columns are populated.
#' @param snp_info   Data frame with \code{CHR}, \code{POS}, \code{REF},
#'   \code{ALT}. Required for \code{alleles} column.
#' @param sep        Field separator. Default \code{","}.
#' @param na_str     NA string. Default \code{"NA"}.
#' @param hap_info   Data frame of exact per-column metadata from
#'   \code{\link{build_haplotype_feature_matrix}()\$info}. When supplied,
#'   the \code{alleles}, \code{frequency}, \code{n_snps}, \code{CHR},
#'   \code{start_bp}, and \code{end_bp} columns are written directly from
#'   this object without any reconstruction. Recommended - pass
#'   \code{hap_info = feat_out\$info} where \code{feat_out} is the return
#'   value of \code{build_haplotype_feature_matrix()}. Default \code{NULL}
#'   (falls back to legacy reconstruction from \code{haplotypes} and
#'   \code{snp_info}).
#' @param min_freq   Minimum frequency used when computing \code{alleles}.
#'   Default \code{0.01}.
#' @param missing_string Missing genotype marker. Default \code{"."}.
#' @param verbose    Logical. Default \code{TRUE}.
#' @return Invisibly returns \code{out_file}.
#' @export
write_haplotype_numeric <- function(hap_matrix, out_file,
                                    haplotypes     = NULL,
                                    snp_info       = NULL,
                                    hap_info       = NULL,
                                    sep            = "\t",
                                    na_str         = "NA",
                                    min_freq       = 0.01,
                                    missing_string = ".",
                                    verbose        = TRUE) {
  # Output orientation: haplotype alleles as ROWS, individuals as COLUMNS.
  # Metadata columns (before individual columns):
  #   hap_id, CHR, start_bp, end_bp, n_snps, alleles, frequency
  #
  # hap_info (from build_haplotype_feature_matrix()$info) is the authoritative
  # metadata source. When supplied, no reconstruction is performed. This fixes
  # the alleles/frequency inconsistency that occurred when the writer tried to
  # reverse-engineer metadata from the raw haplotypes list.

  hnm     <- colnames(hap_matrix)
  ind_ids <- rownames(hap_matrix)
  if (is.null(hnm))     hnm     <- paste0("hap", seq_len(ncol(hap_matrix)))
  if (is.null(ind_ids)) ind_ids <- paste0("Ind", seq_len(nrow(hap_matrix)))
  nhap <- length(hnm)

  # -- Metadata: use hap_info when available (authoritative path) ------------
  if (!is.null(hap_info)) {
    hap_info  <- as.data.frame(hap_info, stringsAsFactors = FALSE)
    idx       <- match(hnm, hap_info$hap_id)

    chr_col    <- hap_info$CHR[idx]
    spos_col   <- hap_info$start_bp[idx]
    epos_col   <- hap_info$end_bp[idx]
    n_snps_col <- hap_info$n_snps[idx]
    freq_col   <- round(hap_info$frequency[idx], 4)

    # Decode dosage strings to nucleotide sequences using snp_info REF/ALT.
    # Uses hap_info$hap_string (e.g. "021") + snp_info REF/ALT for the block's
    # SNP range. Falls back to the raw dosage string when snp_info is unavailable.
    # This is a direct per-row decode, no re-traversal of haplotypes list needed.
    iupac_map <- c(AG="R",GA="R",CT="Y",TC="Y",GC="S",CG="S",
                   AT="W",TA="W",GT="K",TG="K",AC="M",CA="M")
    if (!is.null(snp_info) &&
        all(c("CHR","POS","REF","ALT") %in% names(snp_info))) {
      snp_info_w <- snp_info
      snp_info_w$CHR <- sub("^chr","",as.character(snp_info_w$CHR),ignore.case=TRUE)
      alt_seq_col <- vapply(seq_along(hnm), function(k) {
        hi   <- hap_info[idx[k], ]
        if (is.na(hi$start_bp) || is.na(hi$end_bp) || is.na(hi$hap_string))
          return(hi$hap_string)
        blk_snps <- snp_info_w[snp_info_w$CHR == hi$CHR &
                                 snp_info_w$POS >= hi$start_bp &
                                 snp_info_w$POS <= hi$end_bp, , drop=FALSE]
        blk_snps <- blk_snps[order(blk_snps$POS), ]
        ds <- hi$hap_string
        if (!nzchar(ds) || nrow(blk_snps) == 0L) return(ds)
        chars <- strsplit(ds, "", fixed=TRUE)[[1L]]
        n <- min(length(chars), nrow(blk_snps))
        paste(vapply(seq_len(n), function(i) {
          switch(chars[i],
                 "0" = as.character(blk_snps$REF[i]),
                 "2" = as.character(blk_snps$ALT[i]),
                 "1" = { key <- paste0(blk_snps$REF[i], blk_snps$ALT[i])
                 if (!is.na(iupac_map[key])) unname(iupac_map[key]) else "N" },
                 "N")
        }, character(1L)), collapse="")
      }, character(1L))
    } else {
      alt_seq_col <- hap_info$hap_string[idx]
    }

  } else {
    # -- Fallback: reconstruct from column names (legacy path) ---------------
    # Used when hap_info is not supplied (e.g. direct call without pipeline).
    # NOTE: this path can produce wrong alleles/frequency for multi-SNP blocks
    # where multiple SNPs share the same bp coordinate. Prefer supplying hap_info.
    parts    <- strsplit(hnm, "_", fixed = TRUE)
    chr_col  <- vapply(parts, function(x) if (length(x)>=2L) x[2L] else "NA", character(1L))
    spos_col <- suppressWarnings(
      as.integer(vapply(parts, function(x) if (length(x)>=3L) x[3L] else "0", character(1L))))
    epos_col <- suppressWarnings(
      as.integer(vapply(parts, function(x) if (length(x)>=4L) x[4L] else "0", character(1L))))
    spos_col[is.na(spos_col)] <- 0L
    epos_col[is.na(epos_col)] <- 0L

    bi         <- if (!is.null(haplotypes)) attr(haplotypes, "block_info") else NULL
    alt_seq_col <- rep(NA_character_, nhap)
    freq_col    <- rep(NA_real_, nhap)
    n_snps_col  <- rep(NA_integer_, nhap)

    iupac <- c(AG="R",GA="R",CT="Y",TC="Y",GC="S",CG="S",
               AT="W",TA="W",GT="K",TG="K",AC="M",CA="M")
    decode_str <- function(dstr, ref_v, alt_v) {
      chars <- strsplit(dstr, "", fixed=TRUE)[[1L]]
      n     <- min(length(chars), length(ref_v))
      paste(vapply(seq_len(n), function(i) {
        switch(chars[i], "0"=ref_v[i], "2"=alt_v[i],
               "1"={ key <- paste0(ref_v[i],alt_v[i])
               if (!is.na(iupac[key])) unname(iupac[key]) else "N" }, "N")
      }, character(1L)), collapse="")
    }

    if (!is.null(haplotypes) && !is.null(snp_info) &&
        all(c("CHR","POS","REF","ALT") %in% names(snp_info))) {
      snp_info$CHR <- sub("^chr","",as.character(snp_info$CHR),ignore.case=TRUE)
      for (j in seq_len(nhap)) {
        bn       <- sub("_hap[0-9]+$","",hnm[j])
        hap_rank <- suppressWarnings(as.integer(sub(".*_hap","",hnm[j])))
        if (is.na(hap_rank)) next
        hap_strs <- haplotypes[[bn]]; if (is.null(hap_strs)) next
        blk_snps <- snp_info[snp_info$CHR==chr_col[j] &
                               snp_info$POS>=spos_col[j] &
                               snp_info$POS<=epos_col[j],,drop=FALSE]
        blk_snps <- blk_snps[order(blk_snps$POS),]
        # Use block_info n_snps as authoritative count
        if (!is.null(bi)) {
          brow <- bi[bi$block_id==bn,,drop=FALSE]
          if (nrow(brow)) n_snps_col[j] <- brow$n_snps[1L]
        }
        if (!nrow(blk_snps)) next
        if (is.na(n_snps_col[j])) n_snps_col[j] <- nrow(blk_snps)
        ref_v <- as.character(blk_snps$REF); alt_v <- as.character(blk_snps$ALT)
        miss  <- grepl(".",hap_strs,fixed=TRUE)
        gam   <- hap_strs[!miss]; if (!length(gam)) next
        tbl   <- sort(table(gam),decreasing=TRUE)
        fv    <- as.numeric(tbl)/sum(tbl); tbl <- tbl[fv>=min_freq]
        tops  <- names(tbl)
        if (hap_rank<=length(tops)) {
          alt_hapstr <- tops[hap_rank]
          # Guard against truncated decode when multiple SNPs share same bp
          if (length(ref_v) >= nchar(alt_hapstr))
            alt_seq_col[j] <- decode_str(alt_hapstr, ref_v, alt_v)
          freq_col[j] <- round(as.numeric(tbl[hap_rank])/sum(tbl), 4)
        }
      }
    }
    # Fill remaining n_snps from block_info
    if (!is.null(bi)) {
      for (j in seq_len(nhap)) {
        if (!is.na(n_snps_col[j])) next
        bn   <- sub("_hap[0-9]+$","",hnm[j])
        brow <- bi[bi$block_id==bn,,drop=FALSE]
        if (nrow(brow)) n_snps_col[j] <- brow$n_snps[1L]
      }
    }
  }

  # Transpose: haplotypes as rows, individuals as columns
  t_mat <- t(hap_matrix)
  out <- data.frame(
    hap_id    = hnm,
    CHR       = chr_col,
    start_bp  = spos_col,
    end_bp    = epos_col,
    n_snps    = n_snps_col,
    alleles   = alt_seq_col,
    frequency = freq_col,
    as.data.frame(t_mat, check.names=FALSE),
    check.names=FALSE, stringsAsFactors=FALSE, row.names=NULL
  )
  data.table::fwrite(out, out_file, sep=sep, quote=FALSE, na=na_str)
  if (verbose) message("[write_haplotype_numeric] ", out_file,
                       " (", nhap, " haplotypes x ", length(ind_ids), " individuals)")
  invisible(out_file)
}


#' Decode Haplotype Strings to Nucleotide Sequences
#'
#' Converts the dosage-encoded haplotype strings produced by
#' \code{\link{extract_haplotypes}} (e.g. \code{"02110"}) into
#' nucleotide sequences (e.g. \code{"AGTT?"}) using the REF and ALT
#' alleles of each SNP in the block.
#'
#' @details
#' Each character in a haplotype string is the dosage at one SNP in the block:
#' \itemize{
#'   \item \code{"0"} = homozygous REF  -> REF nucleotide (e.g. \code{A})
#'   \item \code{"1"} = heterozygous    -> IUPAC ambiguity code (e.g. \code{R} for A/G)
#'   \item \code{"2"} = homozygous ALT  -> ALT nucleotide (e.g. \code{G})
#'   \item \code{"."} = missing         -> \code{N}
#' }
#'
#' The result is a data frame with one row per unique haplotype allele per
#' block, showing its nucleotide sequence, frequency, and the REF/ALT at each
#' SNP position. This is the most interpretable representation of what each
#' haplotype allele actually encodes biologically.
#'
#' @param haplotypes List from \code{\link{extract_haplotypes}}.
#' @param snp_info   Data frame with columns \code{SNP}, \code{CHR},
#'   \code{POS}, \code{REF}, \code{ALT}. Must contain all SNPs in the blocks.
#' @param min_freq   Minimum haplotype frequency to include. Default \code{0.01}.
#' @param top_n      Integer or \code{NULL}. Maximum alleles per block.
#'   \code{NULL} (default) retains all above \code{min_freq}.
#' @param missing_string Missing genotype marker. Default \code{"."}.
#'
#' @return A data frame with columns:
#' \describe{
#'   \item{block_id}{Block identifier.}
#'   \item{CHR}{Chromosome.}
#'   \item{start_bp, end_bp}{Block boundaries.}
#'   \item{hap_rank}{Rank by frequency (1 = most common).}
#'   \item{hap_id}{Column name as it appears in the feature matrix.}
#'   \item{dosage_string}{Raw dosage string e.g. \code{"02110"}.}
#'   \item{nucleotide_sequence}{Decoded nucleotide string e.g. \code{"AGTT?"}.}
#'   \item{frequency}{Observed frequency across non-missing individuals.}
#'   \item{n_carriers}{Number of individuals carrying this haplotype (dosage > 0).}
#'   \item{snp_positions}{Semicolon-separated CHR:POS of each SNP in the block.}
#'   \item{snp_alleles}{Semicolon-separated REF/ALT for each SNP.}
#' }
#'
#' @examples
#' data(ldx_geno, ldx_snp_info, ldx_blocks, package = "LDxBlocks")
#' haps <- extract_haplotypes(ldx_geno, ldx_snp_info, ldx_blocks, min_snps = 3)
#' decoded <- decode_haplotype_strings(haps, ldx_snp_info)
#' head(decoded[, c("block_id","hap_rank","dosage_string",
#'                   "nucleotide_sequence","frequency")])
#'
#' @export
decode_haplotype_strings <- function(haplotypes, snp_info,
                                     min_freq       = 0.01,
                                     top_n          = NULL,
                                     missing_string = ".") {

  bi <- attr(haplotypes, "block_info")
  if (is.null(bi))
    stop("haplotypes must carry a block_info attribute from extract_haplotypes().",
         call. = FALSE)
  if (!all(c("CHR","POS","REF","ALT") %in% names(snp_info)))
    stop("snp_info must have columns CHR, POS, REF, ALT.", call. = FALSE)

  snp_info$CHR <- .norm_chr_hap(snp_info$CHR)
  rows_out     <- vector("list", length(haplotypes))

  for (bk in seq_along(haplotypes)) {

    bn  <- names(haplotypes)[bk]
    hap <- haplotypes[[bn]]
    brow <- bi[bi$block_id == bn, , drop = FALSE]
    if (!nrow(brow)) next

    chr <- brow$CHR[1L]
    sb  <- brow$start_bp[1L]
    eb  <- brow$end_bp[1L]

    # SNPs in this block from snp_info (ordered by position)
    blk_snps <- snp_info[snp_info$CHR == chr &
                           snp_info$POS >= sb &
                           snp_info$POS <= eb, , drop = FALSE]
    blk_snps <- blk_snps[order(blk_snps$POS), ]
    n_snps   <- nrow(blk_snps)

    if (!n_snps) next

    ref_v <- as.character(blk_snps$REF)
    alt_v <- as.character(blk_snps$ALT)
    pos_v <- blk_snps$POS

    # Frequency table of haplotype strings
    miss <- grepl(missing_string, hap, fixed = TRUE)
    gam  <- hap[!miss]
    if (!length(gam)) next

    tbl  <- sort(table(gam), decreasing = TRUE)
    fv   <- as.numeric(tbl) / sum(tbl)
    tbl  <- tbl[fv >= min_freq]
    if (!length(tbl)) next
    tops <- if (is.null(top_n)) names(tbl)
    else names(tbl)[seq_len(min(as.integer(top_n), length(tbl)))]

    # Decode each top haplotype string to nucleotides
    decode_one <- function(dstr) {
      chars <- strsplit(dstr, "", fixed = TRUE)[[1L]]
      n     <- min(length(chars), n_snps)
      iupac <- c(AG="R",GA="R",CT="Y",TC="Y",GC="S",CG="S",
                 AT="W",TA="W",GT="K",TG="K",AC="M",CA="M")
      nuc   <- character(n)
      for (i in seq_len(n)) {
        nuc[i] <- switch(chars[i],
                         "0" = ref_v[i],
                         "2" = alt_v[i],
                         "1" = { key <- paste0(ref_v[i], alt_v[i])
                         if (!is.na(iupac[key])) unname(iupac[key]) else "N" },
                         "." = "N",
                         "N"
        )
      }
      paste(nuc, collapse = "")
    }

    n_carriers_fn <- function(dstr) sum(hap == dstr, na.rm = TRUE)

    for (r in seq_along(tops)) {
      dstr <- tops[r]
      rows_out[[length(rows_out) + 1L]] <- data.frame(
        block_id           = bn,
        CHR                = chr,
        start_bp           = as.integer(sb),
        end_bp             = as.integer(eb),
        hap_rank           = r,
        hap_id             = paste0(bn, "_hap", r),
        dosage_string      = dstr,
        nucleotide_sequence = decode_one(dstr),
        frequency          = round(as.numeric(tbl[dstr]) / sum(tbl), 4),
        n_carriers         = n_carriers_fn(dstr),
        snp_positions      = paste(paste0(chr, ":", pos_v), collapse = ";"),
        snp_alleles        = paste(paste0(ref_v, "/", alt_v), collapse = ";"),
        stringsAsFactors   = FALSE
      )
    }
  }

  if (!length(rows_out))
    return(data.frame())

  do.call(rbind, rows_out)
}


#' Write Haplotype Diversity Table
#' @param diversity Data frame from compute_haplotype_diversity().
#' @param out_file Output CSV path.
#' @param append_summary Append genome-wide mean row. Default TRUE.
#' @param verbose Logical. Default TRUE.
#' @return Invisibly returns out_file.
#' @export
write_haplotype_diversity <- function(diversity, out_file,
                                      append_summary=TRUE, verbose=TRUE) {
  if (append_summary) {
    sr <- data.frame(block_id="GENOME",CHR="ALL",start_bp=NA_integer_,
                     end_bp=NA_integer_,n_snps=NA_integer_,n_ind=NA_integer_,
                     n_haplotypes=round(mean(diversity$n_haplotypes,na.rm=TRUE),1),
                     He=mean(diversity$He,na.rm=TRUE),
                     Shannon=mean(diversity$Shannon,na.rm=TRUE),
                     n_eff_alleles=round(mean(diversity$n_eff_alleles,na.rm=TRUE),3),
                     freq_dominant=mean(diversity$freq_dominant,na.rm=TRUE),
                     sweep_flag=NA,
                     phased=NA,stringsAsFactors=FALSE)
    diversity <- rbind(diversity,sr)
  }
  data.table::fwrite(diversity,out_file,sep=",",quote=FALSE,na="NA")
  if(verbose) message("[write_haplotype_diversity] ",out_file)
  invisible(out_file)
}


# ==============================================================================
# Tong et al. (2025) Haplotype Stacking Module
# Theor Appl Genet 138:267. https://doi.org/10.1007/s00122-025-05045-0
# ==============================================================================


#' Backsolve SNP Effects from GEBV (Tong et al. 2025)
#'
#' @description
#' Derives per-SNP additive effect estimates from genome-wide genomic estimated
#' breeding values (GEBV) without re-fitting a marker model. This implements
#' Step 2 of the Tong et al. (2025) haplotype stacking pipeline:
#'
#' \deqn{\hat{\alpha} = \frac{M^\top G^{-1} \hat{g}}{2 \sum_t p_t(1-p_t)}}
#'
#' where \eqn{M} is the centred genotype matrix, \eqn{G} is the VanRaden GRM,
#' \eqn{\hat{g}} are the GEBV, and \eqn{p_t} are allele frequencies.
#'
#' @details
#' This approach is preferred over direct marker-effect estimation when GEBV
#' are already available from a GBLUP run (e.g. from ASReml-R, sommer, or
#' rrBLUP). It avoids refitting the marker model and produces marker effects
#' on the same scale as the original GEBV.
#'
#' Missing genotype values are mean-imputed per column before computation.
#'
#' @param geno_matrix Numeric matrix (individuals x SNPs), values 0/1/2/NA.
#'   Row names must match names of \code{gebv}.
#' @param gebv Named numeric vector of GEBV, one per individual. Names must
#'   match \code{rownames(geno_matrix)}.
#' @param G Optional pre-computed VanRaden GRM (n x n). If \code{NULL}
#'   (default), computed internally via \code{compute_haplotype_grm()}-style
#'   logic. Supply your own if you used a bended or tuned G in the GBLUP.
#'
#' @return Named numeric vector of length p (SNPs), one effect per SNP.
#'   Names match \code{colnames(geno_matrix)}.
#'
#' @references
#' Tong J et al. (2025). Haplotype stacking to improve stability of stripe
#' rust resistance in wheat. \emph{Theoretical and Applied Genetics}
#' \strong{138}:267. \doi{10.1007/s00122-025-05045-0}
#'
#' VanRaden PM (2008). Efficient methods to compute genomic predictions.
#' \emph{Journal of Dairy Science} \strong{91}(11):4414-4423.
#' \doi{10.3168/jds.2007-0980}
#'
#' @examples
#' \dontrun{
#' # After fitting GBLUP with rrBLUP:
#' # fit  <- rrBLUP::kin.blup(data, geno="id", pheno="trait", K=G)
#' # gebv <- fit$g
#' snp_fx <- backsolve_snp_effects(geno_matrix = my_geno, gebv = gebv)
#' head(sort(abs(snp_fx), decreasing = TRUE))
#' }
#'
#' @export
backsolve_snp_effects <- function(geno_matrix, gebv, G = NULL) {
  if (!is.matrix(geno_matrix))
    geno_matrix <- as.matrix(geno_matrix)

  # Align individuals
  common <- intersect(rownames(geno_matrix), names(gebv))
  if (length(common) == 0L)
    stop("No matching individuals between geno_matrix rownames and gebv names.",
         call. = FALSE)
  geno_matrix <- geno_matrix[common, , drop = FALSE]
  g_hat       <- gebv[common]

  # Mean-impute NA per column
  for (j in seq_len(ncol(geno_matrix))) {
    na_j <- is.na(geno_matrix[, j])
    if (any(na_j))
      geno_matrix[na_j, j] <- mean(geno_matrix[, j], na.rm = TRUE)
  }

  # Allele frequencies and scaling denominator
  p     <- colMeans(geno_matrix) / 2
  p     <- pmax(pmin(p, 1 - 1e-8), 1e-8)
  denom <- 2 * sum(p * (1 - p))

  # Centre genotype matrix: M_ij = g_ij - 2*p_j
  M <- sweep(geno_matrix, 2, 2 * p, "-")

  # Compute or use supplied GRM
  if (is.null(G)) {
    G <- tcrossprod(M) / denom
  } else {
    G <- G[common, common, drop = FALSE]
  }

  # Invert G (add small ridge for numerical stability)
  G_inv <- tryCatch(
    solve(G + diag(1e-6, nrow(G))),
    error = function(e)
      stop("GRM inversion failed. Try supplying a pre-bended G via the G= argument.",
           call. = FALSE)
  )

  # Backsolve: alpha = M' G^{-1} g_hat / denom
  alpha <- as.numeric(crossprod(M, G_inv %*% g_hat)) / denom
  names(alpha) <- colnames(geno_matrix)
  alpha
}


#' Compute Local Haplotype GEBV per Block (Tong et al. 2025)
#'
#' @description
#' For each LD block, computes the local GEBV (haplotype effect) for every
#' individual by summing the per-SNP additive effects of the alleles they
#' carry within the block:
#'
#' \deqn{\text{local GEBV}_f = \sum_{t \in f} \left[ h_t \alpha_{1t} +
#'   (1-h_t)\alpha_{0t} \right]}
#'
#' where \eqn{h_t \in \{0, 0.5, 1\}} is the allele dosage at SNP \eqn{t}
#' scaled to [0,1], \eqn{\alpha_{1t}} is the ALT allele effect, and
#' \eqn{\alpha_{0t} = -\alpha_{1t}} is the REF allele effect (by symmetry
#' of the additive model).
#'
#' Blocks are ranked by \code{Var(local GEBV)} -- blocks with high variance
#' contribute strongly to trait differences among individuals and likely
#' harbour causal loci (Tong et al. 2025).
#'
#' @param geno_matrix Numeric matrix (individuals x SNPs), values 0/1/2/NA.
#' @param snp_info Data frame with columns \code{SNP}, \code{CHR}, \code{POS}.
#' @param blocks Block table from \code{\link{run_Big_LD_all_chr}}.
#' @param snp_effects Named numeric vector of per-SNP additive effects from
#'   \code{\link{backsolve_snp_effects}} or from a marker model directly.
#' @param scale Logical. If \code{TRUE} (default), scale \code{Var(local GEBV)}
#'   to [0,1] so blocks are comparable across traits and datasets.
#'
#' @return A list with two elements:
#' \describe{
#'   \item{\code{local_gebv}}{Numeric matrix (individuals x blocks) of per-block
#'     local GEBV values.}
#'   \item{\code{block_importance}}{Data frame with one row per block: block_id,
#'     CHR, start_bp, end_bp, n_snps, var_local_gebv, var_scaled, important
#'     (logical: scaled variance >= 0.9).}
#' }
#'
#' @references
#' Tong J et al. (2025). Haplotype stacking to improve stability of stripe
#' rust resistance in wheat. \emph{Theoretical and Applied Genetics}
#' \strong{138}:267. \doi{10.1007/s00122-025-05045-0}
#'
#' @examples
#' \dontrun{
#' snp_fx  <- backsolve_snp_effects(my_geno, gebv)
#' loc     <- compute_local_gebv(my_geno, snp_info, blocks, snp_fx)
#' # Blocks with scaled variance >= 0.9 are most important
#' head(loc$block_importance[loc$block_importance$important, ])
#' }
#'
#' @export
compute_local_gebv <- function(geno_matrix, snp_info, blocks,
                               snp_effects, scale = TRUE) {
  if (!is.matrix(geno_matrix))
    geno_matrix <- as.matrix(geno_matrix)

  snp_info$CHR <- .norm_chr_hap(as.character(snp_info$CHR))
  blocks$CHR   <- .norm_chr_hap(as.character(blocks$CHR))

  n_ind   <- nrow(geno_matrix)
  n_blk   <- nrow(blocks)
  ind_ids <- rownames(geno_matrix)

  local_mat  <- matrix(NA_real_, nrow = n_ind, ncol = n_blk,
                       dimnames = list(ind_ids, rep("", n_blk)))
  importance <- vector("list", n_blk)

  for (b in seq_len(n_blk)) {
    blk <- blocks[b, ]
    chr <- as.character(blk$CHR)
    sb  <- as.numeric(blk$start.bp)
    eb  <- as.numeric(blk$end.bp)

    # SNPs in this block
    blk_idx <- which(snp_info$CHR == chr &
                       snp_info$POS >= sb &
                       snp_info$POS <= eb)
    if (!length(blk_idx)) next

    blk_snps <- snp_info$SNP[blk_idx]
    # Only SNPs with known effects
    blk_snps <- intersect(blk_snps, names(snp_effects))
    if (!length(blk_snps)) next

    blk_snps_in_geno <- intersect(blk_snps, colnames(geno_matrix))
    if (!length(blk_snps_in_geno)) next

    G_sub  <- geno_matrix[, blk_snps_in_geno, drop = FALSE]
    alpha  <- snp_effects[blk_snps_in_geno]

    # Mean-impute NA
    for (j in seq_len(ncol(G_sub))) {
      na_j <- is.na(G_sub[, j])
      if (any(na_j)) G_sub[na_j, j] <- mean(G_sub[, j], na.rm = TRUE)
    }

    # local GEBV_f = sum_t [ h_t * alpha_1t + (1-h_t) * alpha_0t ]
    # With additive coding: alpha_0t = -alpha_1t (REF effect = -ALT effect)
    # dosage 0 -> contribution = -alpha, 1 -> 0, 2 -> +alpha
    # equivalently: local_gebv = (G_sub - 1) %*% alpha  (centred at het)
    # But Tong et al. use: h_t in {0, 0.5, 1}, so scale dosage by 0.5
    h    <- G_sub / 2          # scale to [0, 1]
    local_mat[, b] <- as.numeric(h %*% alpha + (1 - h) %*% (-alpha))

    bid <- paste0("block_", chr, "_", sb, "_", eb)
    importance[[b]] <- data.frame(
      block_id  = bid,
      CHR       = chr,
      start_bp  = as.integer(sb),
      end_bp    = as.integer(eb),
      n_snps    = length(blk_snps_in_geno),
      var_local_gebv = var(local_mat[, b], na.rm = TRUE),
      stringsAsFactors = FALSE
    )
    colnames(local_mat)[b] <- bid
  }

  importance <- do.call(rbind, Filter(Negate(is.null), importance))
  if (is.null(importance) || !nrow(importance))
    stop("No blocks could be evaluated. Check that SNP names match between ",
         "geno_matrix, snp_info, and snp_effects.", call. = FALSE)

  # Scale variance to [0,1]
  max_var <- max(importance$var_local_gebv, na.rm = TRUE)
  importance$var_scaled <- if (max_var > 0)
    importance$var_local_gebv / max_var else 0
  importance$important <- importance$var_scaled >= 0.9
  importance <- importance[order(importance$CHR, importance$start_bp), ]

  # Trim local_mat to blocks that were evaluated
  valid_blk_ids <- importance$block_id
  local_mat <- local_mat[, colnames(local_mat) %in% valid_blk_ids, drop = FALSE]

  list(local_gebv = local_mat, block_importance = importance)
}


#' Prepare Genomic Prediction Inputs for External GBLUP Software
#'
#' @description
#' Assembles the inputs required to fit a GBLUP model in external software
#' (rrBLUP, sommer, ASReml-R, BGLR) and subsequently run the Tong et al.
#' (2025) haplotype stacking pipeline. Returns the VanRaden GRM computed
#' from the haplotype feature matrix, aligned with a user-supplied phenotype
#' table.
#'
#' @section Workflow:
#' LDxBlocks handles genotype processing and block detection. Phenotype
#' handling and GBLUP fitting are intentionally left to dedicated R packages
#' because phenotype data requires preprocessing (multi-environment adjustment,
#' outlier removal, covariate inclusion) that is dataset-specific. The
#' handoff is:
#'
#' \enumerate{
#'   \item \strong{LDxBlocks} (this function): produce aligned G matrix and
#'     phenotype vector.
#'   \item \strong{External GBLUP} (rrBLUP / sommer / ASReml-R / BGLR):
#'     fit the model, obtain GEBV.
#'   \item \strong{LDxBlocks} (\code{\link{backsolve_snp_effects}} +
#'     \code{\link{compute_local_gebv}}): derive block-level haplotype effects
#'     from the GEBV.
#' }
#'
#' @section Example GBLUP calls after this function:
#' \preformatted{
#' # rrBLUP
#' library(rrBLUP)
#' fit  <- kin.blup(data = inp$pheno_df, geno = "id",
#'                  pheno = "trait", K = inp$G)
#' gebv <- fit$g
#'
#' # sommer
#' library(sommer)
#' fit  <- sommer::mmes(trait ~ 1, random = ~vsm(ism(id), Gu = inp$G),
#'              data = inp$pheno_df)
#' gebv <- fit$U$`u:id`$trait
#'
#' # BGLR
#' library(BGLR)
#' fit  <- BGLR(y = inp$y_vec, ETA = list(list(K = inp$G, model = "RKHS")))
#' gebv <- fit$yHat
#' }
#'
#' @param hap_matrix Numeric matrix (individuals x haplotype alleles) from
#'   \code{\link{build_haplotype_feature_matrix}}.
#' @param pheno_df Data frame of phenotypes. Must contain:
#'   \itemize{
#'     \item An ID column (set via \code{id_col}, default \code{"id"}).
#'       Values must match \code{rownames(hap_matrix)} exactly
#'       (case-sensitive, no leading/trailing spaces).
#'     \item One or more numeric trait columns (referenced via
#'       \code{trait_col}). Column names are arbitrary.
#'     \item \code{NA} values in trait columns are allowed.
#'   }
#'   Minimal single-trait format:
#'   \tabular{rr}{
#'     id \tab YLD \cr
#'     G001 \tab 4.21 \cr
#'     G002 \tab 3.87 \cr
#'     G003 \tab NA \cr
#'   }
#' @param id_col Name of the individual ID column in \code{pheno_df}.
#'   Default \code{"id"}.
#' @param trait_col Name of the trait column to extract as a numeric vector.
#'   Default \code{NULL} -- no \code{y_vec} is returned, only the aligned data
#'   frame.
#' @param bend Logical. Add 0.001 to diagonal of G for positive-definiteness.
#'   Default \code{TRUE} (recommended for mixed model solvers).
#'
#' @return A named list:
#' \describe{
#'   \item{\code{G}}{VanRaden GRM (n x n), aligned to individuals present
#'     in both \code{hap_matrix} and \code{pheno_df}.}
#'   \item{\code{pheno_df}}{Phenotype data frame subset and reordered to match
#'     rows of \code{G}.}
#'   \item{\code{y_vec}}{Named numeric vector of the requested trait (only if
#'     \code{trait_col} is supplied). \code{NA} values are preserved so the
#'     user can decide how to handle them (e.g. set to \code{NA} for
#'     prediction-only individuals in a training/validation split).}
#'   \item{\code{n_train}}{Number of individuals with non-missing trait values.}
#'   \item{\code{n_predict}}{Number of individuals with missing trait values
#'     (prediction candidates).}
#' }
#'
#' @references
#' Tong J et al. (2025). Haplotype stacking to improve stability of stripe
#' rust resistance in wheat. \emph{Theoretical and Applied Genetics}
#' \strong{138}:267. \doi{10.1007/s00122-025-05045-0}
#'
#' VanRaden PM (2008). Efficient methods to compute genomic predictions.
#' \emph{Journal of Dairy Science} \strong{91}(11):4414-4423.
#' \doi{10.3168/jds.2007-0980}
#'
#' @examples
#' \dontrun{
#' # After building haplotype feature matrix:
#' feat  <- build_haplotype_feature_matrix(haps, top_n = 5)$matrix
#' pheno <- read.csv("phenotypes.csv")   # columns: id, YLD, PHT, ...
#'
#' inp <- prepare_gblup_inputs(feat, pheno, id_col = "id",
#'                              trait_col = "YLD")
#'
#' # Fit GBLUP with rrBLUP
#' library(rrBLUP)
#' fit  <- kin.blup(data = inp$pheno_df, geno = "id",
#'                  pheno = "YLD", K = inp$G)
#' gebv <- fit$g
#'
#' # Then derive block-level haplotype effects
#' snp_fx <- backsolve_snp_effects(geno_matrix, gebv, G = inp$G)
#' loc    <- compute_local_gebv(geno_matrix, snp_info, blocks, snp_fx)
#' }
#'
#' @export
prepare_gblup_inputs <- function(hap_matrix, pheno_df,
                                 id_col    = "id",
                                 trait_col = NULL,
                                 bend      = TRUE) {
  if (!id_col %in% names(pheno_df))
    stop("id_col '", id_col, "' not found in pheno_df. ",
         "Set id_col to the column containing individual IDs.", call. = FALSE)

  if (!is.null(trait_col) && !trait_col %in% names(pheno_df))
    stop("trait_col '", trait_col, "' not found in pheno_df.", call. = FALSE)

  # Align individuals present in both hap_matrix and pheno_df
  geno_ids  <- rownames(hap_matrix)
  pheno_ids <- as.character(pheno_df[[id_col]])

  common <- intersect(geno_ids, pheno_ids)
  if (!length(common))
    stop("No matching individual IDs between hap_matrix rownames and ",
         "pheno_df[[id_col]]. Check that ID formats match exactly ",
         "(case-sensitive, no trailing spaces).", call. = FALSE)

  if (length(common) < length(geno_ids))
    message("[prepare_gblup_inputs] ", length(geno_ids) - length(common),
            " genotyped individuals not in pheno_df -- excluded from G.")
  if (length(common) < length(pheno_ids))
    message("[prepare_gblup_inputs] ", length(pheno_ids) - length(common),
            " phenotyped individuals not in hap_matrix -- excluded from output.")

  hap_sub   <- hap_matrix[common, , drop = FALSE]
  pheno_sub <- pheno_df[match(common, pheno_ids), , drop = FALSE]
  rownames(pheno_sub) <- NULL

  # Compute VanRaden GRM from aligned haplotype matrix
  G <- compute_haplotype_grm(hap_sub, bend = bend)

  # Build output
  out <- list(
    G         = G,
    pheno_df  = pheno_sub,
    y_vec     = NULL,
    n_train   = NA_integer_,
    n_predict = NA_integer_
  )

  if (!is.null(trait_col)) {
    y <- as.numeric(pheno_sub[[trait_col]])
    names(y)       <- common
    out$y_vec      <- y
    out$n_train    <- sum(!is.na(y))
    out$n_predict  <- sum(is.na(y))
    message("[prepare_gblup_inputs] ", out$n_train, " training individuals",
            " (non-missing ", trait_col, "), ",
            out$n_predict, " prediction candidates (NA).")
  }

  out
}


#' Haplotype Prediction and Block Importance from Pre-Adjusted Phenotypes
#'
#' @description
#' Runs the complete Tong et al. (2025) haplotype stacking pipeline using
#' pre-adjusted phenotype values (BLUEs, BLUPs, or adjusted entry means).
#' Accepts either a single trait or multiple traits simultaneously.
#'
#' When a single trait is supplied, GBLUP is fitted via
#' \code{rrBLUP::kin.blup()} and block importance is ranked by
#' \code{Var(local GEBV)} for that trait.
#'
#' When multiple traits are supplied, a single trait-agnostic GRM is computed
#' once and shared across all traits. \code{rrBLUP::kin.blup()} is fitted per
#' trait using this shared GRM. Block importance is summarised across traits
#' with per-trait columns plus cross-trait aggregates.
#'
#' @section Input format for \code{blues}:
#' Accepts any of three formats:
#' \itemize{
#'   \item A \strong{named numeric vector} -- single trait, names are genotype
#'     IDs: \code{c(G001 = 4.2, G002 = 3.8)}. \code{id_col}/\code{blue_col}
#'     are ignored.
#'   \item A \strong{data frame with one trait column} -- single trait,
#'     \code{id_col} and \code{blue_col} specify the columns.
#'   \item A \strong{data frame with multiple trait columns} -- multi-trait,
#'     \code{id_col} names the ID column, \code{blue_cols} names the trait
#'     columns (or \code{NULL} to use all numeric non-ID columns).
#'   \item A \strong{named list of named numeric vectors} -- multi-trait with
#'     potentially different individuals per trait:
#'     \code{list(YLD = c(G001=4.2, ...), DIS = c(G001=0.3, ...))}.
#' }
#'
#' @section Multi-trait GBLUP solver strategy:
#' All traits are fitted using \code{rrBLUP::kin.blup()} in a per-trait loop
#' sharing the same GRM. This produces numerically identical results to
#' \code{sommer::mmes()} for single-trait models (verified empirically to 4
#' decimal places), while avoiding the multi-trait \code{cbind()} formula
#' which fails in sommer <= 4.4.5 (Armadillo fixed-size matrix error in the
#' C++ coefficient matrix construction). Because all traits use the same GRM,
#' cross-trait block importance values are directly comparable.
#' \code{solver_used} is always \code{"rrBLUP"}.
#'
#' @section Block importance -- single trait:
#' \code{Var(local GEBV)} is computed per block after backsolving per-SNP
#' effects from GEBV. Blocks are scaled to [0,1]; those with scaled variance
#' >= 0.9 are flagged \code{important}. This follows Tong et al. (2024).
#'
#' @section Block importance -- multi-trait:
#' Per-trait columns \code{var_scaled_<trait>} and \code{important_<trait>}
#' are added for every trait. Cross-trait aggregates:
#' \describe{
#'   \item{\code{var_scaled_mean}}{Mean scaled variance across all traits.
#'     Primary ranking criterion -- blocks consistently important across traits
#'     are more robust stacking candidates.}
#'   \item{\code{n_traits_important}}{Count of traits for which this block
#'     is flagged important.}
#'   \item{\code{important_any}}{TRUE if important for at least one trait.}
#'   \item{\code{important_all}}{TRUE if important for all traits.}
#'   \item{\code{important}}{Controlled by \code{importance_rule}.}
#' }
#'
#' @param geno_matrix Numeric matrix (individuals x SNPs), values 0/1/2/NA.
#'   Row names must be genotype IDs.
#' @param snp_info    Data frame with columns \code{SNP}, \code{CHR},
#'   \code{POS}.
#' @param blocks      Block table from \code{\link{run_Big_LD_all_chr}}.
#' @param blues       Pre-adjusted phenotype values. See section above for
#'   accepted formats.
#' @param id_col      Name of the ID column when \code{blues} is a data
#'   frame. Default \code{"id"}.
#' @param blue_col    Name of the BLUE column for single-trait data frames.
#'   Default \code{"blue"}.
#' @param blue_cols   Character vector of trait column names for multi-trait
#'   data frames. Default \code{NULL} -- all numeric non-ID columns are used.
#' @param importance_rule How to set the combined \code{important} flag for
#'   multi-trait results:
#'   \itemize{
#'     \item \code{"any"} (default): TRUE if important for >= 1 trait.
#'     \item \code{"all"}: TRUE only if important for all traits.
#'     \item \code{"mean"}: TRUE if \code{var_scaled_mean} >= 0.9.
#'   }
#'   Ignored for single-trait runs (single-trait \code{important} is always
#'   scaled variance >= 0.9).
#' @param top_n       Integer or \code{NULL}. Max haplotype alleles per block.
#'   Default \code{NULL}.
#' @param min_freq    Minimum haplotype allele frequency. Default \code{0.01}.
#' @param min_snps    Minimum SNPs per block. Default \code{3L}.
#' @param bend        Logical. Add ridge to GRM diagonal. Default \code{TRUE}.
#' @param verbose     Logical. Print progress. Default \code{TRUE}.
#'
#' @return Named list. For single-trait runs, contains:
#' \describe{
#'   \item{\code{blocks}, \code{diversity}, \code{hap_matrix},
#'     \code{haplotypes}, \code{G}}{Core pipeline outputs.}
#'   \item{\code{n_blocks}, \code{n_hap_columns}}{Summary counts.}
#'   \item{\code{n_traits}}{Integer: 1 for single-trait.}
#'   \item{\code{traits}}{Character vector of trait names.}
#'   \item{\code{solver_used}}{Character: always \code{"rrBLUP"}.}
#'   \item{\code{gebv}}{Named numeric vector of GEBV.}
#'   \item{\code{snp_effects}}{Per-SNP additive effects.}
#'   \item{\code{local_gebv}}{Matrix (individuals x blocks).}
#'   \item{\code{block_importance}}{Data frame with \code{block_id},
#'     coordinates, \code{var_local_gebv}, \code{var_scaled},
#'     \code{important}.}
#'   \item{\code{n_train}, \code{n_predict}}{Training/prediction counts.}
#' }
#' For multi-trait runs, the same structure but \code{gebv},
#' \code{snp_effects}, \code{local_gebv} are named lists (one per trait),
#' and \code{block_importance} contains additional per-trait and cross-trait
#' columns as described in the \emph{Block importance} section.
#' \code{block_importance_list} is added as a named list of single-trait
#' block importance data frames.
#'
#' @references
#' Tong J et al. (2025). Haplotype stacking to improve stability of stripe
#' rust resistance in wheat. \emph{Theoretical and Applied Genetics}
#' \strong{138}:267. \doi{10.1007/s00122-025-05045-0}
#'
#' Tong J et al. (2024). Stacking beneficial haplotypes from the Vavilov
#' wheat collection. \emph{Theoretical and Applied Genetics} \strong{137}:274.
#' \doi{10.1007/s00122-024-04784-w}
#'
#' Endelman JB (2011). Ridge regression and other kernels for genomic
#' selection with R package rrBLUP. \emph{Plant Genome} \strong{4}:250-255.
#' \doi{10.3835/plantgenome2011.08.0024}
#'
#' Covarrubias-Pazaran G (2016). Genome-assisted prediction of quantitative
#' traits using the R package sommer. \emph{PLOS ONE} \strong{11}:e0156744.
#' \doi{10.1371/journal.pone.0156744}
#'
#' @examples
#' \dontrun{
#' library(LDxBlocks)
#' be     <- read_geno("mydata.vcf.gz")
#' blocks <- run_Big_LD_all_chr(be, CLQcut = 0.70)
#' geno   <- read_chunk(be, seq_len(be$n_snps))
#' rownames(geno) <- be$sample_ids
#' colnames(geno) <- be$snp_info$SNP
#'
#' # -- Single trait: named vector ---------------------------------------------
#' blues_vec <- c(G001 = 4.21, G002 = 3.87, G003 = 5.14)
#' res <- run_haplotype_prediction(geno, be$snp_info, blocks, blues = blues_vec)
#' res$block_importance[res$block_importance$important, ]
#' sort(res$gebv, decreasing = TRUE)
#'
#' # -- Single trait: data frame -----------------------------------------------
#' blues_df <- read.csv("blues.csv")   # columns: Genotype, YLD_BLUE
#' res <- run_haplotype_prediction(geno, be$snp_info, blocks,
#'                                  blues    = blues_df,
#'                                  id_col   = "Genotype",
#'                                  blue_col = "YLD_BLUE")
#'
#' # -- Multi-trait: data frame ------------------------------------------------
#' blues_mt <- read.csv("blues_mt.csv")  # columns: id, YLD, DIS, PHT
#' res_mt <- run_haplotype_prediction(geno, be$snp_info, blocks,
#'                                     blues           = blues_mt,
#'                                     id_col          = "id",
#'                                     blue_cols       = c("YLD","DIS","PHT"),
#'                                     importance_rule = "any")
#' res_mt$n_traits        # 3
#' res_mt$solver_used     # "rrBLUP"
#' res_mt$block_importance[
#'   res_mt$block_importance$important_any,
#'   c("block_id","var_scaled_YLD","var_scaled_DIS","var_scaled_mean",
#'     "n_traits_important")]
#'
#' # -- Multi-trait: named list (different individuals per trait) --------------
#' res_mt2 <- run_haplotype_prediction(geno, be$snp_info, blocks,
#'   blues = list(
#'     YLD = c(G001 = 4.2, G002 = 3.8),
#'     DIS = c(G001 = 0.3, G003 = 0.7)
#'   ))
#' }
#'
#' @export
run_haplotype_prediction <- function(geno_matrix,
                                     snp_info,
                                     blocks,
                                     blues,
                                     id_col          = "id",
                                     blue_col        = "blue",
                                     blue_cols       = NULL,
                                     importance_rule = c("any", "all", "mean"),
                                     top_n           = NULL,
                                     min_freq        = 0.01,
                                     min_snps        = 3L,
                                     bend            = TRUE,
                                     verbose         = TRUE) {

  importance_rule <- match.arg(importance_rule)
  .log <- function(...) if (verbose) message("[run_haplotype_prediction] ", ...)

  # -- Detect and normalise blues input into a named list of named vectors ----
  # Four accepted formats -- all normalised to list(trait = named_numeric_vec)
  blues_list <- .parse_blues(blues, id_col, blue_col, blue_cols)
  traits      <- names(blues_list)
  n_traits    <- length(traits)
  is_mt       <- n_traits > 1L

  if (is_mt) {
    .log(n_traits, " traits: ", paste(traits, collapse = ", "))
  } else {
    .log("Single-trait mode: ", traits[1L])
  }

  geno_matrix <- as.matrix(geno_matrix)
  geno_ids    <- rownames(geno_matrix)
  if (is.null(geno_ids))
    stop("geno_matrix must have row names (genotype IDs).", call. = FALSE)

  # -- Step 1: Haplotypes ----------------------------------------------------
  .log("Extracting haplotypes from ", nrow(blocks), " blocks ...")
  haps <- extract_haplotypes(geno_matrix, snp_info, blocks, min_snps = min_snps)
  .log(length(haps), " blocks extracted")

  # -- Step 2: Feature matrix ------------------------------------------------
  .log("Building haplotype feature matrix ...")
  feat_out <- build_haplotype_feature_matrix(haps, top_n = top_n, min_freq = min_freq)
  feat <- feat_out$matrix
  .log(nrow(feat), " ind x ", ncol(feat), " haplotype allele columns")

  # -- Step 3: Single shared GRM ---------------------------------------------
  G <- compute_haplotype_grm(feat, bend = bend)

  # -- Step 4: Build per-trait phenotype alignment table ---------------------
  n_train_vec   <- stats::setNames(integer(n_traits), traits)
  n_predict_vec <- stats::setNames(integer(n_traits), traits)
  pheno_wide    <- data.frame(gid = geno_ids, stringsAsFactors = FALSE)

  for (tr in traits) {
    bv_tr                  <- blues_list[[tr]]
    common_tr              <- intersect(geno_ids, names(bv_tr))
    y_tr                   <- rep(NA_real_, length(geno_ids))
    names(y_tr)            <- geno_ids
    y_tr[common_tr]        <- bv_tr[common_tr]
    pheno_wide[[tr]]       <- y_tr
    n_train_vec[tr]        <- sum(!is.na(y_tr))
    n_predict_vec[tr]      <- sum(is.na(y_tr))
    if (n_train_vec[tr] == 0L)
      stop("Trait '", tr, "': no matching IDs between geno_matrix and blues.",
           call. = FALSE)
  }

  # -- Step 5: GBLUP ---------------------------------------------------------
  gebv_list   <- stats::setNames(vector("list", n_traits), traits)
  solver_used <- "rrBLUP"

  # sommer::mmes() cbind() multi-trait formula crashes in sommer <= 4.4.5
  # (Armadillo fixed-size matrix error in C++ coefficient matrix construction).
  # Empirically, single-trait mmes(Gu=G) == rrBLUP::kin.blup() to 4 d.p.
  # so there is no accuracy benefit to using sommer here.
  # solver_used remains 'rrBLUP' for all inputs.

  # rrBLUP fallback (or single-trait path)
  if (solver_used == "rrBLUP") {
    if (is_mt) .log("Fitting rrBLUP per-trait loop (", n_traits, " traits) ...")
    for (tr in traits) {
      if (is_mt) .log("  Fitting: ", tr)
      pheno_tr <- data.frame(gid = pheno_wide$gid, y = pheno_wide[[tr]],
                             stringsAsFactors = FALSE)
      fit_tr <- tryCatch(
        rrBLUP::kin.blup(data = pheno_tr, geno = "gid", pheno = "y", K = G),
        error = function(e)
          stop("rrBLUP::kin.blup() failed for trait '", tr, "': ",
               conditionMessage(e),
               "\nTry bend = TRUE if GRM is not positive definite.", call. = FALSE)
      )
      gebv_list[[tr]] <- stats::setNames(fit_tr$g, names(fit_tr$g))
    }
  }

  # -- Step 6 & 7: Per-trait SNP effects and local GEBV ----------------------
  snp_fx_list     <- stats::setNames(vector("list", n_traits), traits)
  local_gebv_list <- stats::setNames(vector("list", n_traits), traits)
  bi_list         <- stats::setNames(vector("list", n_traits), traits)

  for (tr in traits) {
    if (is_mt) .log("  Local GEBV: ", tr, " ...")
    common_tr         <- intersect(geno_ids, names(blues_list[[tr]]))
    snp_fx_list[[tr]] <- backsolve_snp_effects(
      geno_matrix[common_tr, , drop = FALSE],
      gebv_list[[tr]][common_tr],
      G = G[common_tr, common_tr, drop = FALSE]
    )
    loc_tr                <- compute_local_gebv(geno_matrix, snp_info, blocks,
                                                snp_fx_list[[tr]])
    local_gebv_list[[tr]] <- loc_tr$local_gebv
    bi_list[[tr]]         <- loc_tr$block_importance
  }

  # -- Step 8: Aggregate block importance ------------------------------------
  .log("Computing haplotype diversity ...")
  diversity <- compute_haplotype_diversity(haps)

  if (!is_mt) {
    # -- Single-trait return ------------------------------------------------
    bi <- bi_list[[traits[1L]]]
    .log("Done. ", sum(bi$important, na.rm = TRUE),
         " important blocks (scaled variance >= 0.9).")
    return(list(
      blocks            = blocks,
      diversity         = diversity,
      hap_matrix        = feat,
      haplotypes        = haps,
      snp_info_filtered = snp_info,
      n_blocks          = nrow(blocks),
      n_hap_columns     = ncol(feat),
      n_traits          = 1L,
      traits            = traits,
      solver_used       = solver_used,
      gebv              = gebv_list[[traits[1L]]],
      snp_effects       = snp_fx_list[[traits[1L]]],
      local_gebv        = local_gebv_list[[traits[1L]]],
      block_importance  = bi,
      G                 = G,
      n_train           = n_train_vec[[traits[1L]]],
      n_predict         = n_predict_vec[[traits[1L]]]
    ))
  }

  # -- Multi-trait aggregation ------------------------------------------------
  .log("Aggregating block importance across ", n_traits, " traits ...")
  ref_bi <- bi_list[[traits[1L]]]
  agg <- data.frame(
    block_id = ref_bi$block_id, CHR = ref_bi$CHR,
    start_bp = ref_bi$start_bp, end_bp = ref_bi$end_bp,
    n_snps   = ref_bi$n_snps, stringsAsFactors = FALSE
  )
  for (tr in traits) {
    bi_tr        <- bi_list[[tr]]
    m_idx        <- match(agg$block_id, bi_tr$block_id)
    agg[[paste0("var_scaled_",  tr)]] <- bi_tr$var_scaled[m_idx]
    imp                               <- bi_tr$important[m_idx]
    imp[is.na(imp)]                   <- FALSE
    agg[[paste0("important_", tr)]]   <- imp
  }
  vs_cols  <- paste0("var_scaled_", traits)
  imp_cols <- paste0("important_", traits)
  agg$var_scaled_mean    <- rowMeans(agg[, vs_cols,  drop = FALSE], na.rm = TRUE)
  agg$n_traits_important <- rowSums( agg[, imp_cols, drop = FALSE], na.rm = TRUE)
  agg$important_any      <- agg$n_traits_important >= 1L
  agg$important_all      <- agg$n_traits_important == n_traits
  agg$important <- switch(importance_rule,
                          "any"  = agg$important_any,
                          "all"  = agg$important_all,
                          "mean" = !is.na(agg$var_scaled_mean) & agg$var_scaled_mean >= 0.9
  )
  agg <- agg[order(agg$CHR, agg$start_bp), ]
  rownames(agg) <- NULL

  .log("Done. ", sum(agg$important, na.rm = TRUE),
       " blocks flagged important (", importance_rule, " rule, ",
       n_traits, " traits).")

  list(
    blocks                = blocks,
    diversity             = diversity,
    hap_matrix            = feat,
    haplotypes            = haps,
    snp_info_filtered     = snp_info,
    n_blocks              = nrow(blocks),
    n_hap_columns         = ncol(feat),
    n_traits              = n_traits,
    traits                = traits,
    solver_used           = solver_used,
    gebv                  = gebv_list,
    snp_effects           = snp_fx_list,
    local_gebv            = local_gebv_list,
    block_importance      = agg,
    block_importance_list = bi_list,
    G                     = G,
    n_train               = n_train_vec,
    n_predict             = n_predict_vec
  )
}

# Internal helper: normalise all blues input formats into a named list
# of named numeric vectors, one entry per trait.
.parse_blues <- function(blues, id_col, blue_col, blue_cols) {
  if (is.numeric(blues) && !is.null(names(blues))) {
    # Named numeric vector -- single trait
    return(list(trait = blues))
  }
  if (is.data.frame(blues)) {
    if (!id_col %in% names(blues))
      stop("id_col '", id_col, "' not found in blues.", call. = FALSE)
    ids <- as.character(blues[[id_col]])
    # Determine trait columns
    if (!is.null(blue_cols)) {
      # Explicit multi-trait columns
      miss <- setdiff(blue_cols, names(blues))
      if (length(miss))
        stop("Columns not found in blues: ", paste(miss, collapse=", "), call. = FALSE)
      trts <- blue_cols
    } else if (blue_col %in% names(blues)) {
      # Single trait via blue_col
      trts <- blue_col
    } else if (blue_col != "blue") {
      # blue_col was explicitly supplied by the user but not found -> error
      stop("blue_col '", blue_col, "' not found in blues. ",
           "Available columns: ", paste(names(blues), collapse = ", "),
           call. = FALSE)
    } else {
      # blue_col is the default 'blue' (user did not specify it) and it is
      # absent -> auto-detect all numeric non-ID columns
      trts <- names(blues)[vapply(blues, is.numeric, logical(1L))]
      trts <- setdiff(trts, id_col)
      if (!length(trts))
        stop("No numeric trait columns found in blues.", call. = FALSE)
    }
    return(stats::setNames(
      lapply(trts, function(tr) {
        v <- as.numeric(blues[[tr]]); names(v) <- ids; v
      }), trts
    ))
  }
  if (is.list(blues) && all(vapply(blues, is.numeric, logical(1L)))) {
    if (is.null(names(blues)))
      stop("blues list must be named (one name per trait).", call. = FALSE)
    return(blues)
  }
  stop("blues must be a named numeric vector, a data frame, or a named list of named numeric vectors.",
       call. = FALSE)
}

#' Integrate GWAS QTL Regions with Haplotype Prediction Results
#'
#' @description
#' Links the output of \code{\link{define_qtl_regions}} (biological evidence
#' from GWAS) with the output of \code{\link{run_haplotype_prediction}}
#' (statistical evidence from haplotype variance) to identify blocks that
#' are supported by both lines of evidence. These are the priority candidates
#' for haplotype stacking in breeding.
#'
#' @details
#' Three complementary evidence layers are combined per block:
#' \describe{
#'   \item{Biological (GWAS)}{Does the block contain a genome-wide significant
#'     marker? Sourced from \code{define_qtl_regions()}.}
#'   \item{Statistical (variance)}{Does the block explain substantial variance
#'     in the trait? Blocks with scaled \code{Var(local GEBV)} >= 0.9 are
#'     flagged \code{important} by \code{run_haplotype_prediction()}.}
#'   \item{Diversity (He)}{Does the block have enough haplotype diversity to
#'     stack favourable alleles? Sourced from
#'     \code{compute_haplotype_diversity()}.}
#' }
#'
#' The output \code{priority_score} is the sum of binary flags for each layer
#' (0-3). Blocks scoring 3 are supported by all three lines of evidence and
#' are the strongest candidates for haplotype stacking. Blocks scoring 2 are
#' worth investigating. Blocks scoring 1 or 0 require caution.
#'
#' @section Interpretation guide:
#' \tabular{rll}{
#'   \strong{Score} \tab \strong{Meaning} \tab \strong{Action} \cr
#'   3 \tab GWAS hit + high variance + diverse \tab Top priority for stacking \cr
#'   2 (GWAS + var) \tab Real effect, low diversity \tab Select across populations \cr
#'   2 (GWAS + div) \tab Real locus, small effect \tab Include if trait is oligogenic \cr
#'   2 (var + div)  \tab Variance explained, no GWAS \tab May be pop. structure -- verify \cr
#'   1 \tab Single evidence only \tab Use with caution \cr
#'   0 \tab No evidence \tab Exclude from stacking \cr
#' }
#'
#' @param qtl_regions   Data frame from \code{\link{define_qtl_regions}}.
#' @param pred_result   List from \code{\link{run_haplotype_prediction}}.
#' @param diversity     Data frame from \code{\link{compute_haplotype_diversity}}.
#'   Optional. If supplied, adds He and sweep_flag to the output.
#' @param He_threshold  Minimum expected heterozygosity to flag a block as
#'   sufficiently diverse for haplotype stacking. Default \code{0.3}.
#'
#' @return Data frame with one row per block that appears in at least one of
#'   the input sources, with columns:
#'   \code{block_id}, \code{CHR}, \code{start_bp}, \code{end_bp},
#'   \code{n_snps}, \code{has_gwas_hit} (logical),
#'   \code{lead_marker}, \code{lead_p}, \code{lead_beta},
#'   \code{n_sig_markers}, \code{is_important} (logical, scaled var >= 0.9),
#'   \code{var_scaled}, \code{He}, \code{sweep_flag},
#'   \code{is_diverse} (logical, He >= He_threshold),
#'   \code{priority_score} (0-3),
#'   \code{recommendation} (character label).
#'   Sorted by \code{priority_score} descending, then \code{var_scaled}
#'   descending.
#'
#' @references
#' Difabachew YF et al. (2023). Genomic prediction with haplotype
#' blocks in wheat. \emph{Frontiers in Plant Science} \strong{14}:1168547.
#' \doi{10.3389/fpls.2023.1168547}
#'
#' Weber SE, Frisch M, Snowdon RJ, Voss-Fels KP (2023). Haplotype
#' blocks for genomic prediction: a comparative evaluation in multiple
#' crop datasets. \emph{Frontiers in Plant Science} \strong{14}:1217589.
#' \doi{10.3389/fpls.2023.1217589}
#'
#' Tong J et al. (2024). Stacking beneficial haplotypes from the Vavilov
#' wheat collection. \emph{Theoretical and Applied Genetics} \strong{137}:274.
#' \doi{10.1007/s00122-024-04784-w}
#'
#' Tong J et al. (2025). Haplotype stacking to improve stability of stripe
#' rust resistance in wheat. \emph{Theoretical and Applied Genetics}
#' \strong{138}:267. \doi{10.1007/s00122-025-05045-0}
#'
#' @examples
#' \dontrun{
#' # 1. GWAS integration
#' gwas   <- read.csv("gwas_results.csv")   # SNP, CHR, POS, P, BETA
#' qtl    <- define_qtl_regions(gwas, blocks, snp_info, p_threshold = 5e-8)
#'
#' # 2. Haplotype prediction
#' blues  <- read.csv("blues.csv")
#' pred   <- run_haplotype_prediction(geno, snp_info, blocks,
#'                                     blues    = blues,
#'                                     id_col   = "id",
#'                                     blue_col = "YLD")
#'
#' # 3. Haplotype diversity
#' haps   <- extract_haplotypes(geno, snp_info, blocks, min_snps = 3)
#' div    <- compute_haplotype_diversity(haps)
#'
#' # 4. Integrate all three
#' priority <- integrate_gwas_haplotypes(qtl, pred, diversity = div)
#'
#' # Top priority blocks for haplotype stacking
#' priority[priority$priority_score == 3, ]
#' }
#'
#' @export
integrate_gwas_haplotypes <- function(qtl_regions,
                                      pred_result,
                                      diversity    = NULL,
                                      He_threshold = 0.3) {

  # -- Validate inputs --------------------------------------------------------
  if (!is.data.frame(qtl_regions))
    stop("qtl_regions must be a data frame from define_qtl_regions().",
         call. = FALSE)
  if (!is.list(pred_result) || !"block_importance" %in% names(pred_result))
    stop("pred_result must be the list returned by run_haplotype_prediction().",
         call. = FALSE)

  bi <- pred_result$block_importance

  # -- Build unified block table from prediction results ----------------------
  out <- data.frame(
    block_id      = bi$block_id,
    CHR           = bi$CHR,
    start_bp      = bi$start_bp,
    end_bp        = bi$end_bp,
    n_snps        = bi$n_snps,
    has_gwas_hit  = FALSE,
    lead_marker      = NA_character_,
    lead_p        = NA_real_,
    lead_beta     = NA_real_,
    n_sig_markers = NA_integer_,
    is_important  = bi$important,
    var_scaled    = bi$var_scaled,
    He            = NA_real_,
    sweep_flag    = NA,
    stringsAsFactors = FALSE
  )

  # -- Merge GWAS evidence ----------------------------------------------------
  if (nrow(qtl_regions) > 0L) {
    for (i in seq_len(nrow(qtl_regions))) {
      qr  <- qtl_regions[i, ]
      idx <- which(out$block_id == qr$block_id)
      if (!length(idx)) {
        # Block has GWAS evidence but was not in pred_result (e.g. < min_snps)
        # Add it with NA prediction columns
        new_row <- data.frame(
          block_id      = qr$block_id,
          CHR           = qr$CHR,
          start_bp      = qr$start_bp,
          end_bp        = qr$end_bp,
          n_snps        = qr$n_snps_block,
          has_gwas_hit  = TRUE,
          lead_marker      = qr$lead_marker,
          lead_p        = if ("lead_p"    %in% names(qr)) qr$lead_p    else NA_real_,
          lead_beta     = if ("lead_beta" %in% names(qr)) qr$lead_beta else NA_real_,
          n_sig_markers = qr$n_sig_markers,
          is_important  = FALSE,
          var_scaled    = NA_real_,
          He            = NA_real_,
          sweep_flag    = NA,
          stringsAsFactors = FALSE
        )
        out <- rbind(out, new_row)
      } else {
        out$has_gwas_hit [idx] <- TRUE
        out$lead_marker     [idx] <- qr$lead_marker
        out$lead_p       [idx] <- if ("lead_p"    %in% names(qr)) qr$lead_p    else NA_real_
        out$lead_beta    [idx] <- if ("lead_beta" %in% names(qr)) qr$lead_beta else NA_real_
        out$n_sig_markers[idx] <- qr$n_sig_markers
      }
    }
  }

  # -- Merge diversity evidence -----------------------------------------------
  if (!is.null(diversity) && nrow(diversity) > 0L) {
    div_match <- match(out$block_id, diversity$block_id)
    if ("He"         %in% names(diversity)) out$He         <- diversity$He        [div_match]
    if ("sweep_flag" %in% names(diversity)) out$sweep_flag <- diversity$sweep_flag[div_match]
  }

  # -- Compute priority score and recommendation ------------------------------
  out$is_diverse <- !is.na(out$He) & out$He >= He_threshold

  out$priority_score <- as.integer(out$has_gwas_hit) +
    as.integer(!is.na(out$is_important) & out$is_important) +
    as.integer(out$is_diverse)

  out$recommendation <- dplyr::case_when(
    out$priority_score == 3L                                       ~ "Top priority: stack",
    out$priority_score == 2L &  out$has_gwas_hit & out$is_important  ~ "Strong: real effect, low diversity",
    out$priority_score == 2L &  out$has_gwas_hit & out$is_diverse    ~ "Moderate: real locus, small effect",
    out$priority_score == 2L & !out$has_gwas_hit                  ~ "Verify: variance explained, no GWAS hit",
    out$priority_score == 1L &  out$has_gwas_hit                  ~ "GWAS only: effect may be small",
    out$priority_score == 1L &  out$is_important                  ~ "Variance only: check for structure",
    out$priority_score == 1L &  out$is_diverse                    ~ "Diverse only: no trait evidence",
    TRUE                                                           ~ "Low priority"
  )

  # Sort: priority descending, then variance descending
  out <- out[order(-out$priority_score,
                   -ifelse(is.na(out$var_scaled), -1, out$var_scaled)), ]
  rownames(out) <- NULL
  out
}


#' Rank Haplotype Blocks by Evidence Strength
#'
#' @description
#' Provides a unified block ranking that works across three use cases,
#' depending on what data the user has available. In all cases,
#' \code{min_freq} filtering inside
#' \code{\link{build_haplotype_feature_matrix}} is applied first as a hard
#' population-level filter -- equivalent to MAF filtering for single SNPs.
#' Blocks are ranked only among those that survive this filter.
#'
#' @section The three use cases:
#' \enumerate{
#'   \item \strong{Genotype only} (no GWAS, no phenotype): blocks ranked by
#'     haplotype diversity (He, effective number of alleles). Blocks with
#'     high diversity have the most potential for haplotype stacking.
#'   \item \strong{Genotype + GWAS} (no phenotype): blocks are binary-flagged
#'     by whether they contain a GWAS-significant marker. Within the
#'     GWAS-hit group and within the non-hit group, blocks are further
#'     ordered by He. The p-value is not used for ranking -- a marker either
#'     crosses the significance threshold or it does not.
#'   \item \strong{Genotype + phenotype} (? GWAS): blocks ranked by scaled
#'     \code{Var(local GEBV)} from \code{\link{run_haplotype_prediction}}.
#'     When GWAS results are also available, the binary GWAS flag is added
#'     as a secondary layer via \code{\link{integrate_gwas_haplotypes}}.
#' }
#'
#' @section Link to \code{min_freq}:
#' \code{min_freq} is a hard population-level pre-filter identical in
#' purpose to MAF filtering for single SNPs: haplotype alleles observed
#' at frequency below this threshold cannot have their effects reliably
#' estimated regardless of trait association, and are dropped before the
#' dosage matrix is built. \code{rank_haplotype_blocks} operates entirely
#' downstream of this filter -- it ranks blocks that have already passed
#' \code{min_freq}, not individual alleles. A block survives as long as
#' at least one of its alleles passes \code{min_freq}.
#'
#' @param diversity     Data frame from \code{\link{compute_haplotype_diversity}}.
#'   Required for all three use cases.
#' @param qtl_regions   Optional. Data frame from
#'   \code{\link{define_qtl_regions}}. When supplied, blocks are
#'   binary-flagged as containing a GWAS hit or not. The p-value is not
#'   used for ranking. Default \code{NULL}.
#' @param pred_result   Optional. List from
#'   \code{\link{run_haplotype_prediction}}. When supplied, blocks are
#'   ranked by \code{Var(local GEBV)}. Default \code{NULL}.
#' @param He_threshold  Minimum He to consider a block diverse enough for
#'   haplotype stacking. Default \code{0.3}.
#' @param top_n_blocks  Return only the top n blocks. Default \code{NULL}
#'   (return all).
#'
#' @return Data frame with one row per block, sorted by evidence strength,
#'   with columns: \code{block_id}, \code{CHR}, \code{start_bp},
#'   \code{end_bp}, \code{n_snps}, \code{He}, \code{n_eff_alleles},
#'   \code{freq_dominant}, \code{sweep_flag}, \code{is_diverse},
#'   \code{has_gwas_hit} (if \code{qtl_regions} supplied),
#'   \code{lead_marker}, \code{lead_beta}, \code{n_sig_markers}
#'   (if \code{qtl_regions} supplied),
#'   \code{var_scaled}, \code{is_important}
#'   (if \code{pred_result} supplied),
#'   \code{use_case}, \code{rank_score}, \code{recommendation}.
#'
#' @references
#' Difabachew YF et al. (2023). Genomic prediction with haplotype
#' blocks in wheat. \emph{Frontiers in Plant Science} \strong{14}:1168547.
#' \doi{10.3389/fpls.2023.1168547}
#'
#' Weber SE, Frisch M, Snowdon RJ, Voss-Fels KP (2023). Haplotype
#' blocks for genomic prediction: a comparative evaluation in multiple
#' crop datasets. \emph{Frontiers in Plant Science} \strong{14}:1217589.
#' \doi{10.3389/fpls.2023.1217589}
#'
#' Tong J et al. (2024). Stacking beneficial haplotypes from the Vavilov
#' wheat collection. \emph{Theoretical and Applied Genetics} \strong{137}:274.
#' \doi{10.1007/s00122-024-04784-w}
#'
#' Tong J et al. (2025). Haplotype stacking to improve stability of
#' stripe rust resistance in wheat. \emph{Theoretical and Applied
#' Genetics} \strong{138}:267. \doi{10.1007/s00122-025-05045-0}
#'
#' @examples
#' \dontrun{
#' haps <- extract_haplotypes(geno, snp_info, blocks)
#' div  <- compute_haplotype_diversity(haps)
#'
#' # Use case 1: genotype only -- rank by diversity
#' ranked <- rank_haplotype_blocks(div)
#'
#' # Use case 2: genotype + GWAS -- binary flag, then diversity within groups
#' qtl    <- define_qtl_regions(gwas_df, blocks, snp_info)
#' ranked <- rank_haplotype_blocks(div, qtl_regions = qtl)
#'
#' # Use case 3: genotype + phenotype -- rank by Var(local GEBV)
#' pred   <- run_haplotype_prediction(geno, snp_info, blocks, blues = blues)
#' ranked <- rank_haplotype_blocks(div, qtl_regions = qtl, pred_result = pred)
#'
#' # Top 20 blocks
#' ranked <- rank_haplotype_blocks(div, pred_result = pred, top_n_blocks = 20)
#' }
#'
#' @export
rank_haplotype_blocks <- function(diversity,
                                  qtl_regions  = NULL,
                                  pred_result  = NULL,
                                  He_threshold = 0.3,
                                  top_n_blocks = NULL) {

  if (!is.data.frame(diversity))
    stop("diversity must be a data frame from compute_haplotype_diversity().",
         call. = FALSE)

  # -- Base table from diversity (all blocks that passed min_freq) ------------
  out <- data.frame(
    block_id      = diversity$block_id,
    CHR           = diversity$CHR,
    start_bp      = diversity$start_bp,
    end_bp        = diversity$end_bp,
    n_snps        = diversity$n_snps,
    He            = diversity$He,
    n_eff_alleles = if ("n_eff_alleles" %in% names(diversity))
      diversity$n_eff_alleles else NA_real_,
    freq_dominant = diversity$freq_dominant,
    sweep_flag    = if ("sweep_flag" %in% names(diversity))
      diversity$sweep_flag else NA,
    stringsAsFactors = FALSE
  )
  out$is_diverse <- !is.na(out$He) & out$He >= He_threshold

  has_gwas <- !is.null(qtl_regions) && is.data.frame(qtl_regions) &&
    nrow(qtl_regions) > 0L
  has_pred <- !is.null(pred_result) && is.list(pred_result) &&
    "block_importance" %in% names(pred_result)

  # -- GWAS flag (binary: hit or no hit, p-value not used for ranking) --------
  out$has_gwas_hit  <- FALSE
  out$lead_marker      <- NA_character_
  out$lead_beta     <- NA_real_
  out$n_sig_markers <- NA_integer_

  if (has_gwas) {
    qtl_match             <- match(out$block_id, qtl_regions$block_id)
    out$has_gwas_hit      <- !is.na(qtl_match)
    out$lead_marker          <- qtl_regions$lead_marker     [qtl_match]
    out$n_sig_markers     <- qtl_regions$n_sig_markers[qtl_match]
    if ("lead_beta" %in% names(qtl_regions))
      out$lead_beta       <- qtl_regions$lead_beta[qtl_match]
  }

  # -- Use case 3: phenotype available ---------------------------------------
  if (has_pred) {
    bi           <- pred_result$block_importance
    bi_match     <- match(out$block_id, bi$block_id)
    out$var_scaled   <- bi$var_scaled[bi_match]
    out$is_important <- bi$important [bi_match]
    out$is_important[is.na(out$is_important)] <- FALSE
    out$use_case     <- "phenotype"

    # Primary sort: var_scaled (trait evidence)
    # Secondary sort within tied var_scaled: gwas_hit flag, then He
    out$rank_score <- ifelse(is.na(out$var_scaled), 0, out$var_scaled) +
      ifelse(out$has_gwas_hit, 0.001, 0) +   # tiny tiebreaker only
      ifelse(out$is_diverse,   0.001, 0)

    out$recommendation <- dplyr::case_when(
      out$is_important &  out$has_gwas_hit &  out$is_diverse ~
        "Top priority: high variance + GWAS hit + diverse",
      out$is_important &  out$has_gwas_hit & !out$is_diverse ~
        "Strong: high variance + GWAS hit, low diversity",
      out$is_important & !out$has_gwas_hit &  out$is_diverse ~
        "Good: high variance + diverse, no GWAS hit",
      out$is_important & !out$has_gwas_hit & !out$is_diverse ~
        "Moderate: high variance, no GWAS hit, low diversity",
      !out$is_important &  out$has_gwas_hit &  out$is_diverse ~
        "Moderate: GWAS hit + diverse, low variance",
      !out$is_important &  out$has_gwas_hit & !out$is_diverse ~
        "Low: GWAS hit only",
      TRUE ~
        "Low priority"
    )

    # -- Use case 2: GWAS only, no phenotype -----------------------------------
  } else if (has_gwas) {
    out$use_case <- "gwas"

    # Primary sort: GWAS hit (binary)
    # Secondary sort within each group: He (diversity)
    out$rank_score <- as.integer(out$has_gwas_hit) +
      ifelse(is.na(out$He), 0, out$He / 10)  # He as tiebreaker only

    out$recommendation <- dplyr::case_when(
      out$has_gwas_hit &  out$is_diverse  ~ "GWAS hit + diverse: good stacking candidate",
      out$has_gwas_hit & !out$is_diverse  ~ "GWAS hit, low diversity: limited stacking potential",
      !out$has_gwas_hit &  out$is_diverse  ~ "No GWAS hit, diverse: background block",
      TRUE                                ~ "No evidence"
    )

    # -- Use case 1: diversity only ---------------------------------------------
  } else {
    out$use_case <- "diversity_only"

    # Rank by He, with n_eff_alleles as tiebreaker
    max_neff <- max(out$n_eff_alleles, na.rm = TRUE)
    if (is.infinite(max_neff) || max_neff == 0) max_neff <- 1
    out$rank_score <- ifelse(is.na(out$He), 0, out$He) +
      ifelse(is.na(out$n_eff_alleles), 0,
             out$n_eff_alleles / max_neff) / 10

    out$recommendation <- dplyr::case_when(
      out$is_diverse & !is.na(out$sweep_flag) & !out$sweep_flag ~
        "Diverse, no sweep: good candidate",
      out$is_diverse  ~
        "Diverse: candidate (verify sweep status)",
      TRUE            ~
        "Low diversity: limited stacking potential"
    )
  }

  # -- Sort and trim ----------------------------------------------------------
  out <- out[order(-out$rank_score), ]
  rownames(out) <- NULL

  if (!is.null(top_n_blocks))
    out <- utils::head(out, as.integer(top_n_blocks))

  # Return ranked table plus all standard outputs for downstream use.
  # hap_matrix is the dimensionality-reduced genotype matrix for GP.
  # blocks, diversity, haplotypes are passed through unchanged so the
  # user does not need to keep multiple result objects.
  result <- list(ranked_blocks = out)

  # Pass through outputs from run_haplotype_prediction if supplied
  if (has_pred) {
    pred_passthrough <- c("blocks", "diversity", "hap_matrix", "haplotypes",
                          "snp_info_filtered", "n_blocks", "n_hap_columns",
                          "gebv", "snp_effects", "local_gebv",
                          "block_importance", "G", "n_train", "n_predict")
    for (nm in pred_passthrough)
      if (nm %in% names(pred_result)) result[[nm]] <- pred_result[[nm]]
  }

  # Always attach diversity for use cases 1 and 2
  if (!"diversity" %in% names(result)) result$diversity <- diversity

  result
}
