# ==============================================================================
# pipeline.R
# End-to-end haplotype block pipeline: one file in, analysis outputs out.
#
# run_ldx_pipeline() wraps the full workflow:
#   1. Auto-detect / convert input to GDS for large files (streaming)
#   2. MAF filtering
#   3. Genome-wide LD block detection via run_Big_LD_all_chr()
#   4. Haplotype extraction and diversity analysis
#   5. Build haplotype genotype matrix (numeric or character/nucleotide)
#   6. Write all outputs to user-specified paths
#
# The user provides only a file path. Format detection and GDS conversion
# are handled internally, mirroring the OptSLDP scale-aware architecture.
# ==============================================================================



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
#' Both output formats have haplotype alleles as rows and individuals as
#' columns, with four metadata columns before the individual columns:
#' \code{hap_id}, \code{CHR}, \code{start_bp}, \code{end_bp},
#' \code{n_snps}, \code{alleles} (nucleotide sequence), \code{frequency}.
#' \itemize{
#'   \item \code{"numeric"}: individual cells are 0/1/2/NA dosage.
#'     0 = does not carry this allele, 2 = carries it (unphased), 1 = one
#'     gamete carries it (phased only). Compatible with rrBLUP, BGLR, ASReml-R.
#'   \item \code{"character"}: individual cells are the nucleotide sequence
#'     (e.g. \code{"AGTTA"}) if the individual carries this allele, \code{"-"}
#'     if absent, \code{"."} if missing.
#' }
#'
#' @param clean_malformed   Logical. Stream-clean the input file before reading
#'   by removing any lines whose column count does not match the header. Needed
#'   for files from NGSEP and some variant callers. Default \code{FALSE}.
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
#'     \item \code{"numeric"} (default) - Tab-delimited. Rows = haplotype
#'       alleles, columns = individuals. Metadata columns: \code{hap_id},
#'       \code{CHR}, \code{start_bp}, \code{end_bp}, \code{n_snps},
#'       \code{alleles} (nucleotide sequence of this allele), \code{frequency}.
#'       Individual columns contain 0/1/2/NA dosage.
#'     \item \code{"character"} - Tab-delimited. Same row/column orientation.
#'       Individual cells contain the nucleotide sequence of the haplotype
#'       allele if the individual carries it, \code{"-"} if absent, \code{"."}
#'       if missing.
#'   }
#' @param maf_cut           Minimum minor allele frequency. SNPs below this
#'   threshold are removed before block detection. Default `0.05`.
#' @param CLQcut            r^2 threshold for clique edges in CLQD. Higher
#'   values produce tighter, smaller blocks. Default `0.5`.
#' @param method            LD metric: \code{"r2"} (default) or \code{"rV2"}
#'   (requires kinship matrix; see \code{\link{run_Big_LD_all_chr}}).
#' @param leng              Boundary scan half-window in SNPs. Default `200L`.
#'   Reduce to 50-100 for very dense WGS panels.
#' @param subSegmSize       Maximum SNPs per CLQD sub-segment. Controls peak
#'   RAM: `subSegmSize x n_individuals x 8` bytes. Default `1500L`.
#' @param n_threads         OpenMP threads for the C++ LD kernel. Default `1L`.
#' @param min_snps_chr      Skip chromosomes with fewer post-filter SNPs than
#'   this. Default `10L`. Increase to skip unplaced scaffolds.
#' @param min_snps_block    Minimum SNPs per haplotype block. Blocks smaller
#'   than this are excluded from haplotype analysis. Default `3L`.
#' @param top_n             Integer or \code{NULL}. Maximum haplotype alleles per
#'   block in the output matrix. \code{NULL} (default) retains all alleles
#'   above \code{min_freq} -- no cap. Set an integer (e.g. \code{5L}) only to
#'   limit column count for memory reasons.
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
#'   \item{hap_matrix}{Numeric matrix of haplotype dosages
#'     (individuals x haplotype alleles, 0/1/2/NA). Always returned
#'     regardless of \code{hap_format}. Use this for downstream modelling.}
#'   \item{haplotypes}{Named list of haplotype strings from
#'     \code{extract_haplotypes()}. Useful for \code{decode_haplotype_strings()}.}
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
#' #'   # top_n       = 5L,          # optional integer cap; NULL = all above min_freq
#'   # hap_format  = "character", # or "numeric" (default)
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
    hap_format       = c("numeric", "character"),
    maf_cut          = 0.05,
    CLQcut           = 0.5,
    method           = c("r2", "rV2"),
    leng             = 200L,
    subSegmSize      = 1500L,
    n_threads        = 1L,
    min_snps_chr     = 10L,
    min_snps_block   = 3L,
    top_n            = NULL,
    scale_hap_matrix = FALSE,
    chr              = NULL,
    verbose          = TRUE,
    clean_malformed  = FALSE
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
  be <- read_geno(geno_file, clean_malformed = clean_malformed, verbose = verbose)
  .ldx_log("Backend: ", be$type, " | ",
           be$n_samples, " individuals | ", be$n_snps, " SNPs")

  # -- Step 2: MAF filtering --------------------------------------------------
  # Store the original full SNP info before filtering — needed in Step 4
  # to map filtered SNP IDs back to correct column indices in the backend.
  orig_snp_info <- be$snp_info

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
  # Read the MAF-filtered SNPs from the already-open backend `be`.
  # We do NOT open a second backend — on Windows, opening a second connection
  # to the same GDS file while `be` is open causes a file-lock error.
  # The backend's snp_info was updated in Step 2 to the filtered set;
  # read_chunk() with the original full-genome SNP indices loads only those.
  .ldx_log("Loading filtered genotype matrix ...")

  # Map filtered SNP IDs back to column positions in the original full backend
  # We need the original (pre-filter) snp_info to get correct column indices.
  # Since be$snp_info is now filtered, we stored the original before filtering.
  full_idx <- match(be$snp_info$SNP, orig_snp_info$SNP)
  geno_mat <- read_chunk(be, full_idx)             # individuals x filtered SNPs
  rownames(geno_mat) <- be$sample_ids
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
  # Both formats: rows = haplotype alleles, cols = individuals.
  # Metadata columns: hap_id, CHR, start_bp, end_bp, n_snps, alleles, frequency.
  .ldx_log("Writing haplotype matrix (format = ", hap_format, ") ...")

  if (identical(hap_format, "numeric")) {
    # Numeric: individual cells = 0/1/2/NA dosage
    # write_haplotype_numeric() handles transposition and metadata internally
    write_haplotype_numeric(
      hap_matrix = hap_matrix,
      out_file   = out_hap_matrix,
      haplotypes = haplotypes,
      snp_info   = be$snp_info,
      sep        = "\t",
      verbose    = verbose
    )
  } else {
    # Character: individual cells = nucleotide sequence or "-" or "."
    write_haplotype_character(
      haplotypes = haplotypes,
      snp_info   = be$snp_info,
      out_file   = out_hap_matrix,
      verbose    = verbose
    )
  }

  .ldx_log("Haplotype matrix written: ", out_hap_matrix)

  # -- Done -------------------------------------------------------------------
  .ldx_log("Pipeline complete.")
  .ldx_log("  Blocks:              ", nrow(blocks))
  .ldx_log("  Haplotype blocks:    ", n_hap_blocks)
  .ldx_log("  Haplotype columns:   ", n_hap_cols)
  .ldx_log("  Individuals:         ", nrow(hap_matrix))

  # res$hap_matrix: always return the numeric dosage matrix (individuals x
  # haplotype alleles) regardless of hap_format, so downstream R code can
  # use it directly for modelling. The file on disk matches hap_format.
  invisible(list(
    blocks            = blocks,
    diversity         = diversity,
    hap_matrix        = hap_matrix,          # raw numeric, always
    haplotypes        = haplotypes,           # list of dosage strings
    snp_info_filtered = be$snp_info,
    n_blocks          = nrow(blocks),
    n_hap_columns     = n_hap_cols
  ))
}
