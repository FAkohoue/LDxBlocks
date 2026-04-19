# ==============================================================================
# pipeline.R
# End-to-end haplotype block pipeline: one file in, analysis outputs out.
#
# run_ldx_pipeline() wraps the full workflow:
#   1. Open genotype backend (file path, bigmemory, or pre-built backend)
#   2. MAF filtering in C++ chromosome by chromosome
#   3. Optional chromosome subset
#   4. Load MAF-filtered genotype matrix
#   5. Genome-wide LD block detection via run_Big_LD_all_chr()
#   6. Haplotype extraction via extract_haplotypes() (C++ string builder)
#   7. Haplotype diversity via compute_haplotype_diversity()
#   8. Build haplotype feature matrix via build_haplotype_feature_matrix()
#   9. Write all outputs to user-specified paths
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
#' @section File-backed memory-mapped genotype store (\code{use_bigmemory}):
#' When \code{use_bigmemory = TRUE} the pipeline converts the source file into
#' a \code{bigmemory::big.matrix} backed by two binary files on disk before
#' running any analysis. Only the OS pages needed for each sub-segment window
#' are loaded into RAM; the rest of the genome stays on disk. Peak RAM is
#' proportional to \code{n_samples x subSegmSize x bytes_per_cell} rather than
#' the full genome matrix.
#'
#' \strong{Backing files created in \code{bigmemory_path}:}
#' \describe{
#'   \item{\code{ldxblocks_bm.bin}}{Raw binary genotype data.
#'     Size = n_samples x n_snps x bytes_per_cell.
#'     With \code{type = "char"} and 204 samples x 3M SNPs: ~0.6 GB.}
#'   \item{\code{ldxblocks_bm.desc}}{Tiny text descriptor that bigmemory
#'     uses to memory-map the \code{.bin} file. Never edit this manually.}
#'   \item{\code{ldxblocks_bm_snpinfo.rds}}{Cached SNP metadata (SNP, CHR,
#'     POS, REF, ALT). bigmemory does not store metadata, so this file is
#'     saved separately to enable restart without re-reading the source file.}
#' }
#'
#' \strong{Restart behaviour:} if all three files already exist in
#' \code{bigmemory_path} when the pipeline is called, the \code{.bin} is
#' reattached instantly via \code{bigmemory::attach.big.matrix()} -- the
#' source VCF is not touched. To force a rebuild, delete the \code{.bin}
#' and \code{.desc} files.
#'
#' \strong{Choosing \code{bigmemory_type}:}
#' \describe{
#'   \item{\code{"char"} (default, recommended)}{1 signed byte per cell
#'     (range --128..127). Genotype dosage values 0, 1, 2 fit without loss.
#'     8x smaller than \code{"double"}.}
#'   \item{\code{"short"}}{2 bytes per cell. Not needed for 0/1/2 dosage;
#'     use only if your pipeline stores values outside --128..127 in the
#'     same matrix.}
#'   \item{\code{"double"}}{8 bytes per cell. Same as a standard R
#'     \code{matrix()}. Use only if downstream code requires \code{double}
#'     precision and the RAM saving is not needed.}
#' }
#'
#' \strong{When to use \code{use_bigmemory = TRUE}:}
#' \itemize{
#'   \item Peak RAM from GDS streaming exceeds available node memory.
#'   \item You need fast restart: the \code{.bin} persists across R sessions
#'     (supply a persistent \code{bigmemory_path}, not \code{tempdir()}).
#'   \item Multiple pipeline runs share the same panel (one \code{.bin},
#'     many readers -- bigmemory memory-maps are safe for concurrent access).
#' }
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
#' columns, preceded by metadata columns: \code{hap_id}, \code{CHR},
#' \code{start_bp}, \code{end_bp}, \code{n_snps}, \code{alleles}
#' (nucleotide sequence of this allele), \code{frequency}.
#' \itemize{
#'   \item \code{"numeric"}: individual cells are haplotype dosage values:
#'     \itemize{
#'       \item \strong{Phased data}: 0/1/2/NA -- 0 = neither gamete carries
#'         this allele, 1 = one gamete carries it (heterozygous),
#'         2 = both gametes carry it (homozygous).
#'       \item \strong{Unphased data}: 0/1/NA -- 0 = absent, 1 = present
#'         (individual's block-level string matches this allele exactly).
#'         The value 2 is not used for unphased data because the two
#'         chromosomes cannot be distinguished -- an individual homozygous
#'         for this allele and one heterozygous for it produce different
#'         observable strings and are treated as different alleles.
#'     }
#'     Compatible with rrBLUP, BGLR, sommer, ASReml-R.
#'   \item \code{"character"}: individual cells are the full nucleotide
#'     sequence of the allele (e.g. \code{"AGTTA"}) if the individual
#'     carries it, \code{"-"} if absent, \code{"."} if missing.
#'     Heterozygous positions use IUPAC ambiguity codes
#'     (R=A/G, Y=C/T, S=G/C, W=A/T, K=G/T, M=A/C).
#' }
#'
#' @param clean_malformed   Logical. Stream-clean the input file before reading
#'   by removing any lines whose column count does not match the header. Needed
#'   for files from NGSEP and some variant callers. Default \code{FALSE}.
#' @param use_bigmemory     Logical. If \code{TRUE} and \code{geno_source} is a
#'   file path, the source file is first loaded into a file-backed
#'   \code{bigmemory::big.matrix} before block detection. Only the OS pages
#'   needed for each sub-segment window are loaded into RAM, keeping peak memory
#'   proportional to \code{n_samples x subSegmSize} rather than the full genome.
#'   Requires the \pkg{bigmemory} package. Default \code{FALSE}.
#'   Ignored when \code{geno_source} is already an
#'   \code{LDxBlocks_backend} (the supplied backend is used as-is).
#' @param bigmemory_path    Directory where the bigmemory backing files
#'   (\code{.bin} and \code{.desc}) and SNP info cache
#'   (\code{_snpinfo.rds}) are written or read from. Defaults to
#'   \code{tempdir()} so files are cleaned up at the end of the R session.
#'   Supply a persistent directory (e.g. next to the output CSVs) when you
#'   want to reattach the backing file on a restart without re-reading the VCF.
#' @param bigmemory_type    Storage type for the \code{big.matrix}: \code{"char"}
#'   (1 byte per cell, 8x smaller than double, sufficient for 0/1/2 dosage),
#'   \code{"short"} (2 bytes), or \code{"double"} (8 bytes). Default
#'   \code{"char"}.
#' @param geno_source  Path to a genotype file, OR an
#'   \code{LDxBlocks_backend} object already created by \code{\link{read_geno}}
#'   or \code{\link{read_geno_bigmemory}}. Passing a backend skips the
#'   internal \code{read_geno()} call, allowing you to use any pre-built
#'   backend including file-backed \code{bigmemory} stores.
#'   Supported file formats when a path is supplied:
#'   numeric dosage CSV (`.csv`), HapMap (`.hmp.txt`),
#'   VCF (`.vcf`, `.vcf.gz`), SNPRelate GDS (`.gds`),
#'   PLINK BED (`.bed`). Format is detected automatically
#'   from the file extension.
#' @param out_blocks        Path for the LD block table CSV. Columns:
#'   `CHR`, `start`, `end`, `start.rsID`, `end.rsID`, `start.bp`, `end.bp`,
#'   `length_bp`, `length_snps`, `block_name`.
#' @param out_diversity     Path for the haplotype diversity table CSV.
#'   Columns: \code{block_id}, \code{CHR}, \code{start_bp},
#'   \code{end_bp}, \code{n_snps}, \code{n_ind},
#'   \code{n_haplotypes}, \code{He} (Nei 1973 corrected),
#'   \code{Shannon}, \code{n_eff_alleles} (effective number of
#'   alleles = 1/\eqn{\sum p_i^2}), \code{freq_dominant},
#'   \code{sweep_flag} (\code{TRUE} when freq_dominant >= 0.90,
#'   indicating a possible selective sweep or strong founder effect),
#'   \code{phased}.
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
#' @param impute            Character. Missing genotype imputation strategy
#'   applied to the genotype matrix after MAF filtering and before LD block
#'   detection and haplotype extraction. One of:
#'   \describe{
#'     \item{\code{"mean_rounded"} (default)}{Per-SNP mean imputation rounded
#'       to the nearest valid dosage (0, 1, or 2). For each SNP, compute the
#'       mean dosage over non-missing samples and round: mean 0.08 -> 0,
#'       mean 0.94 -> 1, mean 1.82 -> 2. Recommended for WGS data with
#'       partial missingness (e.g. \code{miss20} VCF filter). Ensures
#'       haplotype strings never contain \code{"."}  so all individuals
#'       contribute 0/1 values to the feature matrix.}
#'     \item{\code{"mode"}}{Per-SNP mode imputation. Imputes each missing
#'       value with the most common observed dosage (0, 1, or 2) at that
#'       SNP.}
#'     \item{\code{"none"}}{No imputation. If any missing values remain
#'       after MAF filtering the function stops with an informative error.
#'       Use only when the input data are already fully imputed.}
#'   }
#' @param min_callrate      Numeric in \code{[0, 1]}. Minimum per-SNP call
#'   rate (proportion of non-missing samples) required to retain a SNP.
#'   Applied as the first filtering step (before MAF re-filter and imputation).
#'   SNPs with call rate below this threshold are removed. Default \code{0.0}
#'   (disabled). The three-step order is: (1) call-rate filter, (2) MAF
#'   re-filter on the remaining SNPs, (3) imputation. Example: set
#'   \code{min_callrate = 0.8} to require calls in at least 80\% of samples.
#' @param CLQcut            r^2 threshold for clique edges in CLQD. Higher
#'   values produce tighter, smaller blocks. Default `0.5`.
#' @param method            LD metric: \code{"r2"} (default) or \code{"rV2"}
#'   (requires kinship matrix; see \code{\link{run_Big_LD_all_chr}}).
#' @param kin_method        Whitening method for \code{rV2}: \code{"chol"}
#'   (Cholesky, default, faster) or \code{"eigen"} (eigendecomposition,
#'   more stable for near-singular kinship matrices).
#' @param CLQmode           Clique scoring mode: \code{"Density"} (default,
#'   prefers compact high-density cliques -- recommended for most analyses)
#'   or \code{"Maximal"} (prefers the largest cliques regardless of span).
#' @param clstgap           Maximum base-pair gap allowed within a clique
#'   before it is split. Default \code{40000L} (40 kb). Increase for
#'   populations with long-range LD (e.g. inbred lines); decrease for
#'   high-recombination panels.
#' @param split             Logical. If \code{TRUE}, split cliques whose
#'   SNP span exceeds \code{clstgap} bp at the largest internal gap.
#'   Default \code{FALSE}.
#' @param appendrare        Logical. If \code{TRUE}, SNPs that fail MAF
#'   filtering are appended to the nearest block after detection. Default
#'   \code{FALSE}.
#' @param singleton_as_block Logical. If \code{TRUE}, SNPs that receive
#'   no clique assignment (singletons at recombination hotspots) are
#'   returned as single-SNP blocks in the block table. Default \code{FALSE}
#'   (original Big-LD behaviour -- singletons are silently dropped).
#' @param checkLargest      Logical. If \code{TRUE}, apply a dense-core
#'   pre-pass before clique enumeration on sub-segments with >=500 SNPs to
#'   prevent exponential blowup. Default \code{FALSE}.
#' @param digits            Integer. Round r^2 values to this many decimal
#'   places before clique detection. \code{-1L} (default) disables rounding.
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
#' @param min_freq          Minimum haplotype allele frequency (0--1). Alleles
#'   observed at lower frequency than this threshold are dropped before building
#'   the feature matrix and output files. Default \code{0.01} (1\%). Rare
#'   alleles below this threshold cannot be estimated reliably in typical
#'   training sets and add noise. Lower values retain more rare alleles;
#'   higher values (e.g. \code{0.05}) match the MAF filter applied to SNPs.
#' @param scale_hap_matrix  Logical. If `TRUE`, scale the haplotype matrix
#'   columns to zero mean and unit variance before writing. Useful for
#'   GBLUP-style models. Default `FALSE`.
#' @param max_bp_distance  Integer. Maximum bp distance between a SNP pair
#'   for its r\eqn{^2} to be computed in the LD graph. \code{0L} (default)
#'   computes all pairs. Recommended for WGS panels: \code{500000L} (500 kb).
#'   Reduces O(p\eqn{^2}) LD computation to near-O(p) when set.
#' @param chr               Character vector of chromosome names to process.
#'   `NULL` (default) processes all chromosomes.
#' @param verbose           Logical. Print timestamped progress. Default `TRUE`.
#'
#' @return A named list (invisibly) with elements:
#' \describe{
#'   \item{\code{blocks}}{Data frame of LD blocks from \code{run_Big_LD_all_chr}.}
#'   \item{\code{diversity}}{Data frame of per-block haplotype diversity
#'     metrics: \code{block_id}, \code{CHR}, \code{start_bp},
#'     \code{end_bp}, \code{n_snps}, \code{n_ind}, \code{n_haplotypes},
#'     \code{He} (sample-size corrected expected heterozygosity),
#'     \code{Shannon}, \code{n_eff_alleles}, \code{freq_dominant},
#'     \code{sweep_flag}, \code{phased}.}
#'   \item{\code{hap_matrix}}{Numeric matrix (individuals x haplotype
#'     allele columns) -- the dimensionality-reduced genotype matrix for
#'     genomic prediction. Always returned as a numeric R matrix
#'     regardless of \code{hap_format} (which only controls the
#'     \emph{file} written to \code{out_hap_matrix}). Dosage values:
#'     0/1/2/NA for phased data; 0/1/NA for unphased data.
#'     Compatible with rrBLUP, BGLR, sommer, ASReml-R.}
#'   \item{\code{haplotypes}}{Named list of per-block haplotype dosage strings
#'     from \code{extract_haplotypes()}. Pass to
#'     \code{decode_haplotype_strings()} for nucleotide sequences, or to
#'     \code{rank_haplotype_blocks()} for evidence-based ranking.}
#'   \item{\code{snp_info_filtered}}{Data frame of SNP metadata after
#'     MAF filtering: \code{SNP}, \code{CHR}, \code{POS}, and any
#'     additional columns from the input file.}
#'   \item{\code{geno_matrix}}{Numeric matrix (individuals x SNPs) of
#'     MAF-filtered genotypes (0/1/2/NA). Needed directly by
#'     \code{\link{tune_LD_params}} and
#'     \code{\link{run_haplotype_prediction}} -- avoids reloading
#'     the genotype file after the pipeline completes.}
#'   \item{\code{n_blocks}}{Integer. Total LD blocks detected genome-wide.}
#'   \item{\code{n_hap_columns}}{Integer. Total haplotype allele columns
#'     after \code{min_freq} filtering -- the effective number of predictors
#'     for genomic prediction.}
#' }
#'
#' @examples
#' \donttest{
#' # Path A: supply a file path (original interface, unchanged)
#' geno_file <- system.file("extdata", "example_genotypes_numeric.csv",
#'                          package = "LDxBlocks")
#'
#' res <- run_ldx_pipeline(
#'   geno_source    = geno_file,
#'   out_blocks     = tempfile(fileext = ".csv"),
#'   out_diversity  = tempfile(fileext = ".csv"),
#'   out_hap_matrix = tempfile(fileext = ".csv"),
#'   hap_format     = "numeric",
#'   maf_cut        = 0.05,
#'   CLQcut         = 0.5,
#'   leng           = 10L,
#'   subSegmSize    = 80L,
#'   min_snps_block = 3L,
#'   # top_n       = 5L,   # optional integer cap; NULL = all above min_freq
#'   # hap_format  = "character",  # alternative to "numeric"
#'   verbose        = FALSE
#' )
#'
#' head(res$blocks)
#' head(res$diversity)
#' dim(res$hap_matrix)
#'
#' # Path B: supply a pre-built bigmemory backend
#' if (requireNamespace("bigmemory", quietly = TRUE)) {
#'   # read_geno_bigmemory() accepts a file path directly:
#'   # it calls read_geno() internally and wraps the result.
#'   be_bm <- read_geno_bigmemory(
#'     source      = geno_file,
#'     backingfile = tempfile("ldxbm"),
#'     type        = "char"
#'   )
#'   res2 <- run_ldx_pipeline(
#'     geno_source    = be_bm,
#'     out_blocks     = tempfile(fileext = ".csv"),
#'     out_diversity  = tempfile(fileext = ".csv"),
#'     out_hap_matrix = tempfile(fileext = ".csv"),
#'     CLQcut         = 0.5, leng = 10L, subSegmSize = 70L
#'   )
#'   close_backend(be_bm)
#' }
#' }
#'
#' @seealso \code{\link{run_Big_LD_all_chr}},
#'   \code{\link{extract_haplotypes}},
#'   \code{\link{compute_haplotype_diversity}},
#'   \code{\link{build_haplotype_feature_matrix}},
#'   \code{\link{rank_haplotype_blocks}},
#'   \code{\link{run_haplotype_prediction}},
#'   \code{\link{tune_LD_params}}
#'
#' @export
run_ldx_pipeline <- function(
    geno_source,
    out_blocks,
    out_diversity,
    out_hap_matrix,
    hap_format       = c("numeric", "character"),
    maf_cut          = 0.05,
    CLQcut           = 0.5,
    method           = c("r2", "rV2"),
    kin_method       = "chol",
    CLQmode          = c("Density", "Maximal", "Louvain", "Leiden"),
    leng             = 200L,
    subSegmSize      = 1500L,
    clstgap          = 40000L,
    split            = FALSE,
    appendrare       = FALSE,
    singleton_as_block = FALSE,
    checkLargest     = FALSE,
    digits           = -1L,
    n_threads        = 1L,
    min_snps_chr     = 10L,
    min_snps_block   = 3L,
    top_n            = NULL,
    min_freq         = 0.01,
    scale_hap_matrix = FALSE,
    chr              = NULL,
    verbose          = TRUE,
    max_bp_distance  = 0L,
    clean_malformed  = FALSE,
    use_bigmemory    = FALSE,
    bigmemory_path   = tempdir(),
    bigmemory_type   = "char",
    impute           = c("mean_rounded", "mode", "none"),
    min_callrate     = 0.0
) {
  hap_format    <- match.arg(hap_format)
  method        <- match.arg(method)
  CLQmode       <- match.arg(CLQmode)
  impute        <- match.arg(impute)
  bigmemory_type <- match.arg(bigmemory_type,
                              choices = c("char", "short", "double"))

  .ldx_log <- function(...) {
    if (verbose) message(sprintf("[%s] %s",
                                 format(Sys.time(), "%H:%M:%S"),
                                 paste0(...)))
  }

  # -- Step 1: Open backend ------------------------------------------------
  # Three input paths, in priority order:
  #   (a) geno_source is already an LDxBlocks_backend -> use as-is
  #   (b) use_bigmemory = TRUE -> open file, convert to bigmemory backend
  #   (c) default -> open file via read_geno() (GDS streaming)
  if (inherits(geno_source, "LDxBlocks_backend")) {
    # Path (a): pre-built backend supplied by caller
    be <- geno_source
    .ldx_log("Using pre-built backend: ", be$type, " | ",
             be$n_samples, " ind | ", be$n_snps, " SNPs")
    # Caller owns this backend -- do NOT register on.exit(close_backend)

  } else {
    if (!is.character(geno_source) || length(geno_source) != 1L)
      stop("geno_source must be a file path (character) or an ",
           "LDxBlocks_backend object.", call. = FALSE)

    if (isTRUE(use_bigmemory)) {
      # Path (b): file-backed bigmemory backend
      if (!requireNamespace("bigmemory", quietly = TRUE))
        stop("use_bigmemory = TRUE requires the 'bigmemory' package. ",
             "Install with: install.packages('bigmemory')", call. = FALSE)

      # Normalize and convert to forward slashes. bigmemory's C code
      # uses sprintf("%s/%s", backingpath, file) so backslashes from
      # Windows normalizePath produce mixed separators that fail.
      bigmemory_path <- normalizePath(bigmemory_path, mustWork = FALSE)
      bigmemory_path <- gsub("\\\\", "/", bigmemory_path, fixed = TRUE)
      bm_stem     <- file.path(bigmemory_path, "ldxblocks_bm")
      bm_bin_file <- paste0(bm_stem, ".bin")
      bm_desc_file <- paste0(bm_stem, ".desc")
      bm_si_file  <- paste0(bm_stem, "_snpinfo.rds")

      # -- Determine bigmemory state and act accordingly ------------------
      # Three files must ALL exist for a valid reattach:
      #   .bin  -- raw genotype bytes
      #   .desc -- bigmemory memory-map descriptor
      #   _snpinfo.rds -- SNP metadata (not stored in .desc)
      # Any other combination (partial write, interrupted job, manual
      # deletion) is treated as a stale/corrupt set: all three files
      # are removed automatically and the matrix is rebuilt from scratch.
      # bigmemory < 1.4.7 creates an extensionless backing file (no .bin).
      # Accept either form when checking for existing state.
      bm_bin_noext <- sub("\\.bin$", "", bm_bin_file)  # e.g. .../ldxblocks_bm
      bm_bin_ok  <- file.exists(bm_bin_file) || file.exists(bm_bin_noext)
      bm_desc_ok <- file.exists(bm_desc_file)
      bm_si_ok   <- file.exists(bm_si_file)
      bm_all_ok  <- bm_bin_ok && bm_desc_ok && bm_si_ok
      bm_partial <- (bm_bin_ok || bm_desc_ok || bm_si_ok) && !bm_all_ok

      if (bm_partial) {
        # Partial write detected (e.g. job was killed mid-build).
        # Remove all three files and rebuild cleanly.
        partial <- c(bm_bin_file, bm_desc_file, bm_si_file)[c(bm_bin_ok, bm_desc_ok, bm_si_ok)]
        .ldx_log("[bigmemory] Partial backing files detected (incomplete ",
                 "previous run). Removing and rebuilding:")
        for (f in partial) {
          .ldx_log("[bigmemory]   removing ", basename(f))
          file.remove(f)
        }
        bm_all_ok <- FALSE
      }

      if (bm_all_ok) {
        # State 1: all three files present -> reattach instantly
        .ldx_log("[bigmemory] Reattaching existing backing file: ",
                 basename(bm_bin_file))
        si_bm <- readRDS(bm_si_file)
        be <- read_geno_bigmemory(
          source      = bm_desc_file,
          snp_info    = si_bm,
          backingfile = basename(bm_stem),
          backingpath = bigmemory_path
        )
      } else {
        # State 2: no files (or just cleaned) -> build fresh
        .ldx_log("[bigmemory] Building file-backed matrix (type = '",
                 bigmemory_type, "') ...")
        .ldx_log("[bigmemory]   ", bm_bin_file)
        be_tmp <- read_geno(geno_source,
                            clean_malformed = clean_malformed,
                            verbose = verbose)
        be <- read_geno_bigmemory(
          source      = be_tmp,
          snp_info    = be_tmp$snp_info,
          backingfile = basename(bm_stem),
          backingpath = bigmemory_path,
          type        = bigmemory_type,
          verbose     = verbose
        )
        close_backend(be_tmp)  # release temporary GDS/BED handle
        saveRDS(be$snp_info, bm_si_file)
        .ldx_log("[bigmemory] SNP info cached: ", basename(bm_si_file))
        # Save sample IDs separately (bigmemory rownames unreliable on file-backed)
        bm_sid_file <- paste0(bm_stem, "_sampleids.rds")
        saveRDS(be$sample_ids, bm_sid_file)
        .ldx_log("[bigmemory] Sample IDs cached: ", basename(bm_sid_file))
        .ldx_log("[bigmemory] Rerun with same bigmemory_path to reattach ",
                 "without re-reading the source file.")
      }
      on.exit(close_backend(be), add = TRUE)

    } else {
      # Path (c): default GDS streaming backend
      .ldx_log("Opening genotype file: ", basename(geno_source))
      be <- read_geno(geno_source, clean_malformed = clean_malformed,
                      verbose = verbose)
      on.exit(close_backend(be), add = TRUE)
    }

    .ldx_log("Backend: ", be$type, " | ",
             be$n_samples, " ind | ", be$n_snps, " SNPs")
  }

  # -- Step 2: MAF filtering --------------------------------------------------
  # Store the original full SNP info before filtering -- needed in Step 4
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
  # We do NOT open a second backend -- on Windows, opening a second connection
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

  # -- Step 2.5: Call-rate filter -> MAF re-filter -> Impute (C++) -----------
  # Correct biological order:
  #   1. Call-rate filter : drop SNPs with too much missing data.
  #      Computed on the actual loaded matrix (not the streaming approximation).
  #   2. MAF re-filter    : recompute MAF after call-rate filter.
  #      A SNP with 40% missing that passes the backend MAF threshold may fail
  #      once low-callrate samples are excluded from the denominator.
  #      Uses maf_filter_cpp() on the subset matrix.
  #   3. Impute           : fill remaining NAs in passing SNPs.
  #      Uses impute_and_filter_cpp() with min_callrate = 0 (filter already done).
  # All three steps use C++; applied to geno_mat only (disk data unchanged).

  # -- 2.5a: Call-rate filter ------------------------------------------------
  n_ind  <- nrow(geno_mat)
  n_snps_before_cr <- ncol(geno_mat)
  if (min_callrate > 0) {
    cr_res <- impute_and_filter_cpp(
      geno         = matrix(as.integer(geno_mat), nrow = n_ind),
      min_callrate = min_callrate,
      method       = 0L   # method irrelevant here; we only use $keep
    )
    n_cr_removed <- cr_res$n_filtered
    # Call-rate distribution summary for the log
    cr_vals   <- as.numeric(cr_res$call_rates)
    cr_pct    <- round(cr_vals * 100, 1)
    cr_below  <- sum(cr_vals < min_callrate, na.rm = TRUE)
    cr_min    <- round(min(cr_vals, na.rm = TRUE) * 100, 1)
    cr_median <- round(median(cr_vals, na.rm = TRUE) * 100, 1)
    cr_mean   <- round(mean(cr_vals, na.rm = TRUE) * 100, 1)
    .ldx_log(sprintf(
      "Call-rate filter: threshold = %.0f%% | SNP call rates: min=%.1f%% mean=%.1f%% median=%.1f%%",
      min_callrate * 100, cr_min, cr_mean, cr_median))
    if (n_cr_removed > 0L) {
      keep_cr           <- as.logical(cr_res$keep)
      geno_mat          <- geno_mat[, keep_cr, drop = FALSE]
      snp_info_filtered <- snp_info_filtered[keep_cr, , drop = FALSE]
      be$snp_info       <- snp_info_filtered
      be$n_snps         <- nrow(snp_info_filtered)
      .ldx_log(sprintf(
        "Call-rate filter: removed %d / %d SNPs (%.1f%%) with call rate < %.0f%% | Remaining: %d SNPs",
        n_cr_removed, n_snps_before_cr,
        100 * n_cr_removed / n_snps_before_cr,
        min_callrate * 100, ncol(geno_mat)))
    } else {
      .ldx_log(sprintf(
        "Call-rate filter: all %d SNPs passed (>= %.0f%% call rate).",
        n_snps_before_cr, min_callrate * 100))
    }
  } else {
    .ldx_log("Call-rate filter: disabled (min_callrate = 0).")
  }

  # -- 2.5b: MAF re-filter after call-rate filter ----------------------------
  if (min_callrate > 0 && ncol(geno_mat) > 0L) {
    n_before_maf2 <- ncol(geno_mat)
    keep_maf <- maf_filter_cpp(geno_mat, maf_cut = maf_cut)
    n_maf2_removed <- sum(!keep_maf)
    if (n_maf2_removed > 0L) {
      geno_mat          <- geno_mat[, keep_maf, drop = FALSE]
      snp_info_filtered <- snp_info_filtered[keep_maf, , drop = FALSE]
      be$snp_info       <- snp_info_filtered
      be$n_snps         <- nrow(snp_info_filtered)
      .ldx_log(sprintf(
        "MAF re-filter (>= %.2f after call-rate filter): removed %d / %d SNPs (%.1f%%) | Remaining: %d SNPs",
        maf_cut, n_maf2_removed, n_before_maf2,
        100 * n_maf2_removed / n_before_maf2, ncol(geno_mat)))
    } else {
      .ldx_log(sprintf(
        "MAF re-filter (>= %.2f): all %d SNPs passed.", maf_cut, n_before_maf2))
    }
    if (ncol(geno_mat) == 0L)
      stop("No SNPs remain after call-rate and MAF filtering.")
  }

  # -- 2.5c: Imputation ------------------------------------------------------
  n_na_before   <- sum(is.na(geno_mat))
  n_cells_total <- as.numeric(nrow(geno_mat)) * ncol(geno_mat)
  pct_missing   <- round(100 * n_na_before / n_cells_total, 3)
  if (impute == "none") {
    if (n_na_before > 0L)
      stop("impute = 'none' but ", n_na_before, " missing genotypes remain (",
           pct_missing, "% of ", nrow(geno_mat), " ind x ", ncol(geno_mat),
           " SNPs). Supply fully imputed data or set impute = 'mean_rounded' or 'mode'.")
    .ldx_log(sprintf(
      "Imputation: none | matrix %d ind x %d SNPs | 0 missing values (data complete).",
      nrow(geno_mat), ncol(geno_mat)))
  } else if (n_na_before == 0L) {
    .ldx_log(sprintf(
      "Imputation: skipped | matrix %d ind x %d SNPs | 0 missing values.",
      nrow(geno_mat), ncol(geno_mat)))
  } else {
    .ldx_log(sprintf(
      "Imputation (%s): %d ind x %d SNPs | %d missing values (%.3f%% of matrix) -> filling ...",
      impute, nrow(geno_mat), ncol(geno_mat), n_na_before, pct_missing))
    method_int <- if (impute == "mode") 1L else 0L
    imp_res <- impute_and_filter_cpp(
      geno         = matrix(as.integer(geno_mat), nrow = nrow(geno_mat)),
      min_callrate = 0.0,   # call-rate already handled in step 2.5a
      method       = method_int
    )
    geno_mat <- imp_res$geno_imputed
    storage.mode(geno_mat) <- "numeric"
    rownames(geno_mat) <- be$sample_ids
    colnames(geno_mat) <- snp_info_filtered$SNP
    n_na_after <- sum(is.na(geno_mat))
    .ldx_log(sprintf(
      "Imputation (%s): filled %d / %d missing values (%.3f%%) | Remaining NA: %d",
      impute, imp_res$n_imputed, n_na_before, pct_missing, n_na_after))
  }

  # -- Step 5: Genome-wide LD block detection ---------------------------------
  .ldx_log("Running genome-wide LD block detection ...")
  blocks <- run_Big_LD_all_chr(
    geno_matrix        = geno_mat,
    snp_info           = be$snp_info,
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
    stop("No LD blocks detected. Try lowering CLQcut or adjusting parameters.")

  .ldx_log("Detected ", nrow(blocks), " LD blocks")

  # Write block table
  data.table::fwrite(blocks, file = out_blocks, sep = ",",
                     quote = FALSE, na = "NA")
  .ldx_log("Block table written: ", out_blocks)

  # -- Step 6: Haplotype extraction -------------------------------------------
  # extract_haplotypes() uses build_hap_strings_cpp() (C++) internally,
  # replacing an R vapply() loop. ~20-50x faster for large panels.
  .ldx_log("Extracting haplotypes (min_snps = ", min_snps_block, ") ...")
  # Pass the backend (be) not geno_mat so extract_haplotypes() uses
  # the chromosome-wise streaming path. This avoids repeated O(n_snps)
  # scans of the full 2.96M SNP metadata for each of the 291k blocks.
  haplotypes <- extract_haplotypes(
    geno     = be,
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

  write_haplotype_diversity(
    diversity,
    out_diversity,
    append_summary = TRUE,   # adds genome-wide mean row at the bottom
    verbose        = FALSE
  )
  .ldx_log("Diversity table written: ", out_diversity)

  # -- Step 8: Build haplotype genotype matrix --------------------------------
  .ldx_log("Building haplotype feature matrix (top_n = ", top_n,
           ", min_freq = ", min_freq, ") ...")
  hap_matrix <- build_haplotype_feature_matrix(
    haplotypes     = haplotypes,
    top_n          = top_n,
    min_freq       = min_freq,
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
      min_freq   = min_freq,
      top_n      = top_n,
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
    hap_matrix        = hap_matrix,          # individuals x haplotype allele columns
    haplotypes        = haplotypes,           # raw dosage strings per block
    geno_matrix       = geno_mat,            # individuals x SNPs (MAF-filtered)
    snp_info_filtered = be$snp_info,
    n_blocks          = nrow(blocks),
    n_hap_columns     = n_hap_cols
  ))
}
