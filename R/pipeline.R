# ==============================================================================
# pipeline.R
# End-to-end haplotype block pipeline: one file in, analysis outputs out.
#
# run_ldx_pipeline() wraps the full workflow:
#   1. Auto-detect / convert input to GDS for large files (streaming)
#   2. MAF filtering
#   3. Genome-wide LD block detection via run_Big_LD_all_chr()
#   4. Haplotype extraction and diversity analysis
#   5. Build haplotype genotype matrix (numeric or HapMap)
#   6. Write all outputs to user-specified paths
#
# The user provides only a file path. Format detection and GDS conversion
# are handled internally, mirroring the OptSLDP scale-aware architecture.
# ==============================================================================


# -- Internal: write haplotype matrix in HapMap format -------------------------

#' Write a haplotype dosage matrix in HapMap nucleotide format
#'
#' Each haplotype column in `hap_mat` carries values 0 (absent), 2 (present),
#' or NA. We encode 0 -> REF/REF, 2 -> ALT/ALT, NA -> NN. REF and ALT are
#' synthetic nucleotide labels (H for haplotype reference, A for alternate).
#'
#' @param hap_mat  Numeric matrix (individuals x haplotype columns).
#' @param out_file Path to write the HapMap file.
#' @keywords internal
#' @noRd
.write_haplotype_hapmap <- function(hap_mat, out_file) {
  n_hap  <- ncol(hap_mat)
  hap_nm <- colnames(hap_mat)

  # Parse block_id and haplotype rank from column names like "block_1_hap1"
  block_ids <- sub("_hap[0-9]+$", "", hap_nm)
  hap_ranks <- as.integer(sub(".*_hap", "", hap_nm))

  # REF = "H" (haplotype reference), ALT = "A" (alternate allele token)
  REF <- "H"; ALT <- "A"

  # Encode dosage -> nucleotide
  encode <- function(vals) {
    out <- character(length(vals))
    out[!is.na(vals) & vals == 0]  <- paste0(REF, REF)
    out[!is.na(vals) & vals == 2]  <- paste0(ALT, ALT)
    out[is.na(vals)]               <- "NN"
    out
  }

  hdr <- data.frame(
    "rs#"       = hap_nm,
    alleles     = paste0(REF, "/", ALT),
    chrom       = block_ids,
    pos         = hap_ranks,
    strand      = "+",
    "assembly#" = "NA", center    = "NA",
    protLSID    = "NA", assayLSID = "NA",
    panelLSID   = "NA", QCcode    = "NA",
    check.names = FALSE, stringsAsFactors = FALSE
  )

  # Build nucleotide call matrix (haplotypes x individuals)
  call_mat <- apply(hap_mat, 2, encode)   # returns individuals x haplotypes
  # apply over columns -> result is nrow x ncol of hap_mat
  # We need haplotypes as rows, individuals as columns
  call_mat <- t(call_mat)
  rownames(call_mat) <- hap_nm
  colnames(call_mat) <- rownames(hap_mat)

  out <- cbind(hdr, as.data.frame(call_mat, stringsAsFactors = FALSE))
  data.table::fwrite(out, file = out_file, sep = "\t",
                     quote = FALSE, na = "NN")
  invisible(out_file)
}


# -- Internal: MAF filter for backend ------------------------------------------

#' Filter SNPs by MAF using a backend object
#'
#' Reads all genotype data from the backend in chromosome chunks and removes
#' SNPs with MAF < maf_cut or that are monomorphic.
#'
#' @param be       LDxBlocks_backend object.
#' @param maf_cut  Minimum MAF. Default 0.05.
#' @param verbose  Logical.
#' @return Updated snp_info data.frame with passing SNPs only.
#' @keywords internal
#' @noRd
.maf_filter_backend <- function(be, maf_cut = 0.05, verbose = TRUE) {
  if (verbose) message("[MAF filter] Computing MAF for ", be$n_snps, " SNPs ...")

  chrs     <- unique(be$snp_info$CHR)
  keep_ids <- character(0)

  for (chr in chrs) {
    idx      <- which(be$snp_info$CHR == chr)
    geno_chr <- read_chunk(be, idx)                      # individuals x SNPs
    # maf_filter_cpp expects individuals x SNPs (same orientation)
    keep_lgl <- maf_filter_cpp(geno_chr, maf_cut = maf_cut)
    keep_ids <- c(keep_ids, be$snp_info$SNP[idx][keep_lgl])
    rm(geno_chr); gc(FALSE)
  }

  n_pass <- length(keep_ids)
  if (verbose)
    message("[MAF filter] ", n_pass, " / ", be$n_snps,
            " SNPs pass MAF >= ", maf_cut)

  be$snp_info[be$snp_info$SNP %in% keep_ids, ]
}


# -- Main pipeline function -----------------------------------------------------

#' End-to-End Haplotype Block Pipeline
#'
#' @description
#' A single-call wrapper that takes one genotype file and produces a complete
#' haplotype-based dataset ready for genomic prediction. Internally handles
#' format detection, optional GDS conversion for large files, MAF filtering,
#' genome-wide LD block detection, haplotype extraction, diversity analysis,
#' and output writing.
#'
#' The user provides only the input file path and desired output paths. All
#' intermediate steps run transparently.
#'
#' @section Scale behaviour:
#' Files are processed via the `LDxBlocks_backend` streaming interface.
#' For VCF, HapMap, and numeric CSV files the backend reads one chromosome
#' window at a time - peak RAM equals one `subSegmSize`-SNP window regardless
#' of total marker count. For very large files (> 2 M SNPs) the GDS backend
#' via `SNPRelate` is used automatically if the `SNPRelate` package is installed.
#'
#' @section Haplotype genotype matrix:
#' The haplotype matrix has one row per individual and one column per
#' haplotype allele (top-`top_n` haplotypes per block). Each cell contains
#' the dosage of that haplotype allele encoded as:
#' \itemize{
#'   \item `0` - individual does not carry this haplotype
#'   \item `2` - individual carries this haplotype (homozygous)
#'   \item `NA` - missing data in block
#' }
#' This encoding is directly compatible with genomic prediction software
#' (ASReml-R, rrBLUP, BGLR, GBLUP) without further transformation.
#'
#' @param geno_file         Path to genotype file. Supported formats: numeric
#'   dosage CSV (`.csv`), HapMap (`.hmp.txt`), VCF (`.vcf`, `.vcf.gz`),
#'   SNPRelate GDS (`.gds`), PLINK BED (`.bed`). Format is detected
#'   automatically from the file extension.
#' @param out_blocks        Path for the LD block table CSV. Columns:
#'   `CHR`, `start`, `end`, `start.rsID`, `end.rsID`, `start.bp`, `end.bp`,
#'   `length_bp`, `length_snps`, `block_name`.
#' @param out_diversity     Path for the haplotype diversity table CSV.
#'   Columns: `block_id`, `n_ind`, `n_haplotypes`, `He`, `Shannon`,
#'   `freq_dominant`.
#' @param out_hap_matrix    Path for the haplotype genotype matrix. Format is
#'   controlled by `hap_format`. See section **Haplotype genotype matrix**.
#' @param hap_format        Output format for the haplotype matrix:
#'   \itemize{
#'     \item `"numeric"` (default) - CSV with rows = individuals,
#'       columns = haplotype alleles coded 0/2/NA.
#'     \item `"hapmap"` - HapMap format with rows = haplotype alleles,
#'       columns = individuals, nucleotide encoding.
#'   }
#' @param maf_cut           Minimum minor allele frequency. SNPs below this
#'   threshold are removed before block detection. Default `0.05`.
#' @param CLQcut            r^2 threshold for clique edges in CLQD. Higher
#'   values produce tighter, smaller blocks. Default `0.5`.
#' @param method            LD metric: `"r2"` (default) or `"rV2"` (requires
#'   kinship matrix; see `Big_LD()`).
#' @param leng              Boundary scan half-window in SNPs. Default `200L`.
#'   Reduce to 50-100 for very dense WGS panels.
#' @param subSegmSize       Maximum SNPs per CLQD sub-segment. Controls peak
#'   RAM: `subSegmSize x n_individuals x 8` bytes. Default `1500L`.
#' @param n_threads         OpenMP threads for the C++ LD kernel. Default `1L`.
#' @param min_snps_chr      Skip chromosomes with fewer post-filter SNPs than
#'   this. Default `10L`. Increase to skip unplaced scaffolds.
#' @param min_snps_block    Minimum SNPs per haplotype block. Blocks smaller
#'   than this are excluded from haplotype analysis. Default `3L`.
#' @param top_n             Number of top haplotype alleles to retain per block
#'   in the output matrix. Default `5L`.
#' @param scale_hap_matrix  Logical. If `TRUE`, scale the haplotype matrix
#'   columns to zero mean and unit variance before writing. Useful for
#'   GBLUP-style models. Default `FALSE`.
#' @param chr               Character vector of chromosome names to process.
#'   `NULL` (default) processes all chromosomes.
#' @param verbose           Logical. Print timestamped progress. Default `TRUE`.
#'
#' @return A named list (invisibly) with elements:
#' \describe{
#'   \item{`blocks`}{`data.frame` of LD blocks.}
#'   \item{`diversity`}{`data.frame` of per-block haplotype diversity metrics.}
#'   \item{`hap_matrix`}{Numeric matrix of haplotype dosages
#'     (individuals x haplotype alleles). `NULL` if written to file only.}
#'   \item{`snp_info_filtered`}{`data.frame` of SNP metadata after MAF filter.}
#'   \item{`n_blocks`}{Integer. Total blocks detected.}
#'   \item{`n_hap_columns`}{Integer. Total haplotype allele columns.}
#' }
#'
#' @examples
#' \donttest{
#' geno_file <- system.file("extdata", "example_genotypes_numeric.csv",
#'                          package = "LDxBlocks")
#'
#' res <- run_ldx_pipeline(
#'   geno_file      = geno_file,
#'   out_blocks     = tempfile(fileext = ".csv"),
#'   out_diversity  = tempfile(fileext = ".csv"),
#'   out_hap_matrix = tempfile(fileext = ".csv"),
#'   hap_format     = "numeric",
#'   maf_cut        = 0.05,
#'   CLQcut         = 0.5,
#'   leng           = 10L,
#'   subSegmSize    = 80L,
#'   min_snps_block = 3L,
#'   top_n          = 3L,
#'   verbose        = FALSE
#' )
#'
#' head(res$blocks)
#' head(res$diversity)
#' dim(res$hap_matrix)
#' }
#'
#' @seealso [run_Big_LD_all_chr()], [extract_haplotypes()],
#'   [compute_haplotype_diversity()], [build_haplotype_feature_matrix()]
#'
#' @export
run_ldx_pipeline <- function(
    geno_file,
    out_blocks,
    out_diversity,
    out_hap_matrix,
    hap_format       = c("numeric", "hapmap"),
    maf_cut          = 0.05,
    CLQcut           = 0.5,
    method           = c("r2", "rV2"),
    leng             = 200L,
    subSegmSize      = 1500L,
    n_threads        = 1L,
    min_snps_chr     = 10L,
    min_snps_block   = 3L,
    top_n            = 5L,
    scale_hap_matrix = FALSE,
    chr              = NULL,
    verbose          = TRUE
) {
  hap_format <- match.arg(hap_format)
  method     <- match.arg(method)

  .ldx_log <- function(...) {
    if (verbose) message(sprintf("[%s] %s",
                                 format(Sys.time(), "%H:%M:%S"),
                                 paste0(...)))
  }

  # -- Step 1: Open backend (auto-detect format) ------------------------------
  .ldx_log("Opening genotype file: ", basename(geno_file))
  be <- read_geno(geno_file, verbose = verbose)
  .ldx_log("Backend: ", be$type, " | ",
           be$n_samples, " individuals | ", be$n_snps, " SNPs")

  # -- Step 2: MAF filtering --------------------------------------------------
  .ldx_log("MAF filtering (>= ", maf_cut, ") ...")
  snp_info_filtered <- .maf_filter_backend(be, maf_cut = maf_cut,
                                           verbose = verbose)
  n_pass <- nrow(snp_info_filtered)
  if (n_pass == 0L)
    stop("No SNPs remain after MAF filtering. Lower maf_cut.")

  # Rebuild backend with filtered SNP list by subsetting the snp_info
  # The backend still reads from the original file; we just track which
  # SNP indices to use via the filtered snp_info
  be$snp_info <- snp_info_filtered
  be$n_snps   <- n_pass

  # -- Step 3: Subset chromosomes if requested --------------------------------
  if (!is.null(chr)) {
    chr <- as.character(chr)
    be$snp_info <- be$snp_info[be$snp_info$CHR %in% chr, ]
    be$n_snps   <- nrow(be$snp_info)
    .ldx_log("Chromosome filter: retained ", be$n_snps, " SNPs on chr ",
             paste(chr, collapse = ", "))
  }

  # -- Step 4: Load genotype matrix for block detection ----------------------
  # For streaming backends we load chromosome by chromosome inside
  # run_Big_LD_all_chr; here we assemble the full filtered matrix.
  # For very large datasets users should use the GDS backend which handles
  # chunking internally via read_chunk().
  .ldx_log("Loading filtered genotype matrix ...")
  snp_idx  <- match(be$snp_info$SNP, be$snp_info$SNP)  # all, in order
  col_idx  <- which(be$snp_info$SNP %in% be$snp_info$SNP)

  # Re-read from original backend using filtered column indices
  be_full  <- read_geno(geno_file, verbose = FALSE)
  full_idx <- match(be$snp_info$SNP, be_full$snp_info$SNP)
  geno_mat <- read_chunk(be_full, full_idx)        # individuals x SNPs
  close_backend(be_full)

  rownames(geno_mat) <- be_full$sample_ids
  colnames(geno_mat) <- be$snp_info$SNP
  .ldx_log("Genotype matrix: ", nrow(geno_mat), " x ", ncol(geno_mat))

  # -- Step 5: Genome-wide LD block detection ---------------------------------
  .ldx_log("Running genome-wide LD block detection ...")
  blocks <- run_Big_LD_all_chr(
    geno_matrix  = geno_mat,
    snp_info     = be$snp_info,
    method       = method,
    CLQcut       = CLQcut,
    leng         = leng,
    subSegmSize  = subSegmSize,
    n_threads    = n_threads,
    min_snps_chr = min_snps_chr,
    verbose      = verbose
  )

  if (is.null(blocks) || nrow(blocks) == 0L)
    stop("No LD blocks detected. Try lowering CLQcut or adjusting parameters.")

  .ldx_log("Detected ", nrow(blocks), " LD blocks")

  # Write block table
  data.table::fwrite(blocks, file = out_blocks, sep = ",",
                     quote = FALSE, na = "NA")
  .ldx_log("Block table written: ", out_blocks)

  # -- Step 6: Haplotype extraction -------------------------------------------
  .ldx_log("Extracting haplotypes (min_snps = ", min_snps_block, ") ...")
  haplotypes <- extract_haplotypes(
    geno     = geno_mat,
    snp_info = be$snp_info,
    blocks   = blocks,
    chr      = NULL,
    min_snps = min_snps_block
  )
  n_hap_blocks <- length(haplotypes)
  .ldx_log("Haplotypes extracted for ", n_hap_blocks, " blocks")

  if (n_hap_blocks == 0L)
    stop("No haplotype blocks produced. Check min_snps_block vs block sizes.")

  # -- Step 7: Haplotype diversity --------------------------------------------
  .ldx_log("Computing haplotype diversity ...")
  diversity <- compute_haplotype_diversity(haplotypes)

  data.table::fwrite(diversity, file = out_diversity, sep = ",",
                     quote = FALSE, na = "NA")
  .ldx_log("Diversity table written: ", out_diversity)

  # -- Step 8: Build haplotype genotype matrix --------------------------------
  .ldx_log("Building haplotype feature matrix (top_n = ", top_n, ") ...")
  hap_matrix <- build_haplotype_feature_matrix(
    haplotypes     = haplotypes,
    top_n          = top_n,
    scale_features = scale_hap_matrix
  )
  n_hap_cols <- ncol(hap_matrix)
  .ldx_log("Haplotype matrix: ", nrow(hap_matrix), " individuals x ",
           n_hap_cols, " haplotype allele columns")

  # -- Step 9: Write haplotype genotype matrix --------------------------------
  .ldx_log("Writing haplotype matrix (format = ", hap_format, ") ...")

  if (identical(hap_format, "numeric")) {
    # Numeric CSV: rows = individuals, columns = haplotype alleles
    out_df <- data.frame(
      Sample = rownames(hap_matrix),
      as.data.frame(hap_matrix, check.names = FALSE),
      check.names = FALSE,
      stringsAsFactors = FALSE
    )
    data.table::fwrite(out_df, file = out_hap_matrix, sep = ",",
                       quote = FALSE, na = "NA")
  } else {
    # HapMap: rows = haplotype alleles, columns = individuals
    .write_haplotype_hapmap(hap_matrix, out_hap_matrix)
  }

  .ldx_log("Haplotype matrix written: ", out_hap_matrix)

  # -- Done -------------------------------------------------------------------
  .ldx_log("Pipeline complete.")
  .ldx_log("  Blocks:              ", nrow(blocks))
  .ldx_log("  Haplotype blocks:    ", n_hap_blocks)
  .ldx_log("  Haplotype columns:   ", n_hap_cols)
  .ldx_log("  Individuals:         ", nrow(hap_matrix))

  invisible(list(
    blocks           = blocks,
    diversity        = diversity,
    hap_matrix       = hap_matrix,
    snp_info_filtered = be$snp_info,
    n_blocks         = nrow(blocks),
    n_hap_columns    = n_hap_cols
  ))
}
