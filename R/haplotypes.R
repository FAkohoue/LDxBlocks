# ─────────────────────────────────────────────────────────────────────────────
# haplotypes.R  –  Phase-free haplotype reconstruction and diversity analysis
# ─────────────────────────────────────────────────────────────────────────────
#
# These functions extend LDxBlocks beyond block detection into biologically
# actionable haplotype-level analysis:
#
#   extract_haplotypes()         -- phase-free haplotype strings per block
#   compute_haplotype_diversity() -- richness, He, Shannon entropy per block
#   build_haplotype_feature_matrix() -- numeric dosage matrix for prediction
#
# Rationale
# ---------
# LD blocks define genomic segments of correlated variation. Within each block
# the set of SNP alleles co-inherited on the same chromosome represents a
# *haplotype*. Because most livestock / plant breeding datasets are unphased,
# we use a pragmatic *allele-string* representation (concatenated 0/1/2 codes)
# rather than true gametic phases. This is equivalent to a multi-locus genotype
# code and captures most of the haplotype diversity signal without requiring
# statistical phasing.
#
# For prediction, each unique haplotype string is encoded as an integer dosage
# (0/1/2), producing a feature matrix that can be fed directly into GBLUP,
# RR-BLUP, or machine learning models — often outperforming SNP-based models
# because haplotypes capture multi-locus epistatic effects implicitly.
# ─────────────────────────────────────────────────────────────────────────────


# ── Internal: collapse genotype row to allele string ─────────────────────────
.geno_to_string <- function(row) paste(row, collapse = "")


#' Extract Phase-Free Haplotypes Within LD Blocks
#'
#' @description
#' For each LD block in \code{blocks}, extracts the SNP columns from
#' \code{geno} that fall within the block's bp interval and constructs a
#' *phase-free haplotype* for each individual by concatenating their allele
#' codes (0, 1, or 2) into a character string. Missing genotypes are
#' represented as \code{"."}.
#'
#' These haplotype strings serve as multi-locus genotype identifiers and can be
#' used directly with \code{\link{compute_haplotype_diversity}} or
#' \code{\link{build_haplotype_feature_matrix}}.
#'
#' @section Note on phasing:
#' True gametic haplotypes require statistical or read-based phasing (e.g.
#' SHAPEIT, Beagle). The strings produced here are *diploid allele strings*,
#' not gametic phases. They are nonetheless highly informative for diversity
#' metrics and genomic prediction because within a high-LD block, unphased
#' multi-locus codes are nearly 1:1 with true haplotype classes.
#'
#' @param geno Numeric matrix (individuals × SNPs; 0/1/2). Column names must
#'   be SNP IDs matching \code{snp_info$SNP}. Row names are used as individual
#'   IDs.
#' @param snp_info Data frame with columns \code{SNP}, \code{CHR}, \code{POS}.
#' @param blocks Data frame of LD blocks as returned by
#'   \code{\link{run_Big_LD_all_chr}} or \code{\link{Big_LD}}. Must include
#'   columns \code{start.bp}, \code{end.bp}, and (if \code{chr} is given)
#'   \code{CHR}.
#' @param chr Character. If supplied, only the blocks on this chromosome and
#'   the SNPs on this chromosome (from \code{snp_info}) are processed.
#'   \code{NULL} processes all blocks jointly.
#' @param min_snps Integer. Minimum number of SNPs required for a block to be
#'   included. Blocks with fewer SNPs are skipped. Default \code{2}.
#' @param na_char Character. Replacement for \code{NA} alleles in strings.
#'   Default \code{"."}.
#'
#' @return A named list with one element per block (named
#'   \code{"block_<start.bp>_<end.bp>"}). Each element is a character vector
#'   of length \code{nrow(geno)}, giving the haplotype string for each
#'   individual. The list has an attribute \code{"block_info"} — a data frame
#'   summarising block ID, chromosome, positions, and SNP count.
#'
#' @seealso \code{\link{compute_haplotype_diversity}},
#'   \code{\link{build_haplotype_feature_matrix}}
#'
#' @examples
#' \donttest{
#' set.seed(1)
#' geno <- matrix(sample(0:2, 50 * 40, replace = TRUE), 50, 40)
#' rownames(geno) <- paste0("ind", 1:50)
#' colnames(geno) <- paste0("rs", 1:40)
#' snp_info <- data.frame(
#'   SNP = colnames(geno),
#'   CHR = "chr1",
#'   POS = seq(1000, by = 5000, length.out = 40)
#' )
#' # Toy blocks
#' blocks <- data.frame(
#'   start.bp = c(1000, 80000),
#'   end.bp   = c(70000, 200000),
#'   CHR      = "chr1"
#' )
#' haps <- extract_haplotypes(geno, snp_info, blocks)
#' head(haps[[1]])
#' }
#'
#' @export
extract_haplotypes <- function(
    geno,
    snp_info,
    blocks,
    chr      = NULL,
    min_snps = 2L,
    na_char  = "."
) {
  if (!is.matrix(geno)) geno <- as.matrix(geno)

  # ── Subset to chromosome ───────────────────────────────────────────────────
  if (!is.null(chr)) {
    chr      <- as.character(chr)
    if ("CHR" %in% names(snp_info)) snp_info <- snp_info[snp_info$CHR == chr, ]
    if ("CHR" %in% names(blocks))   blocks   <- blocks[blocks$CHR   == chr, ]
  }

  snp_pos <- as.numeric(snp_info$POS)
  snp_ids <- as.character(snp_info$SNP)

  # match SNP columns in geno
  col_idx <- match(snp_ids, colnames(geno))
  if (all(is.na(col_idx))) {
    # fallback: positional match
    col_idx <- seq_len(min(nrow(snp_info), ncol(geno)))
    warning("SNP IDs in snp_info not found in colnames(geno). Using positional matching.")
  }

  result     <- list()
  block_info <- vector("list", nrow(blocks))

  for (b in seq_len(nrow(blocks))) {
    s_bp <- as.numeric(blocks$start.bp[b])
    e_bp <- as.numeric(blocks$end.bp[b])
    in_block <- which(snp_pos >= s_bp & snp_pos <= e_bp)

    if (length(in_block) < min_snps) next

    # columns of geno inside this block
    block_cols <- col_idx[in_block]
    block_cols <- block_cols[!is.na(block_cols)]
    if (length(block_cols) < min_snps) next

    sub_geno   <- geno[, block_cols, drop = FALSE]
    sub_geno[is.na(sub_geno)] <- na_char

    hap_strings <- apply(sub_geno, 1L, function(row) {
      paste(ifelse(is.na(row), na_char, as.character(row)), collapse = "")
    })

    block_name       <- paste0("block_", s_bp, "_", e_bp)
    result[[block_name]] <- hap_strings

    block_info[[b]] <- data.frame(
      block_id  = block_name,
      CHR       = if ("CHR" %in% names(blocks)) as.character(blocks$CHR[b]) else NA_character_,
      start_bp  = s_bp,
      end_bp    = e_bp,
      n_snps    = length(block_cols),
      stringsAsFactors = FALSE
    )
  }

  block_info <- do.call(rbind, Filter(Negate(is.null), block_info))
  attr(result, "block_info") <- block_info
  result
}


#' Compute Haplotype Diversity Metrics Per LD Block
#'
#' @description
#' For each block produced by \code{\link{extract_haplotypes}}, computes:
#' \describe{
#'   \item{Richness (\eqn{k})}{Number of unique haplotype strings.}
#'   \item{Expected heterozygosity (\eqn{H_e})}{
#'     Nei's (1973) gene diversity:
#'     \eqn{H_e = \frac{n}{n-1}\left(1 - \sum_i p_i^2\right)}
#'     where \eqn{p_i} is the frequency of haplotype \eqn{i} and \eqn{n} is
#'     the number of individuals (with non-missing haplotypes).}
#'   \item{Shannon entropy (\eqn{H'})}{
#'     \eqn{H' = -\sum_i p_i \log_2(p_i)}.
#'     Sensitive to both richness and evenness.}
#'   \item{Dominant haplotype frequency (\eqn{f_{\max}})}{
#'     Frequency of the most common haplotype; high values indicate a selective
#'     sweep or strong founder effect.}
#' }
#'
#' @param haplotypes List as returned by \code{\link{extract_haplotypes}}.
#' @param missing_string Character. Haplotype strings containing this pattern
#'   are treated as missing and excluded from frequency calculations.
#'   Default \code{"."}.
#'
#' @return A \code{data.frame} with one row per block and columns:
#'   \code{block_id}, \code{n_ind} (non-missing individuals),
#'   \code{n_haplotypes} (richness), \code{He} (expected heterozygosity),
#'   \code{Shannon} (entropy in bits), \code{freq_dominant} (frequency of most
#'   common haplotype).
#'
#' @references
#' Nei M (1973) Analysis of gene diversity in subdivided populations.
#' \emph{PNAS} \strong{70}(12):3321–3323.
#'
#' @seealso \code{\link{extract_haplotypes}}, \code{\link{build_haplotype_feature_matrix}}
#'
#' @examples
#' \donttest{
#' # (Continuing extract_haplotypes example)
#' set.seed(1)
#' geno <- matrix(sample(0:2, 50 * 40, replace = TRUE), 50, 40)
#' rownames(geno) <- paste0("ind", 1:50)
#' colnames(geno) <- paste0("rs", 1:40)
#' snp_info <- data.frame(SNP = colnames(geno), CHR = "chr1",
#'   POS = seq(1000, by = 5000, length.out = 40))
#' blocks <- data.frame(start.bp = c(1000, 80000), end.bp = c(70000, 200000), CHR = "chr1")
#' haps <- extract_haplotypes(geno, snp_info, blocks)
#' diversity <- compute_haplotype_diversity(haps)
#' print(diversity)
#' }
#'
#' @export
compute_haplotype_diversity <- function(haplotypes, missing_string = ".") {
  block_names <- names(haplotypes)
  res <- lapply(block_names, function(bn) {
    hap     <- haplotypes[[bn]]
    # exclude missing
    is_miss <- grepl(missing_string, hap, fixed = TRUE)
    hap     <- hap[!is_miss]
    n       <- length(hap)
    if (n == 0L) return(data.frame(block_id = bn, n_ind = 0L, n_haplotypes = NA,
                                   He = NA, Shannon = NA, freq_dominant = NA))
    tbl   <- table(hap)
    freqs <- as.numeric(tbl) / n
    k     <- length(tbl)
    He    <- if (n > 1L) (n / (n - 1L)) * (1 - sum(freqs^2)) else NA_real_
    Sh    <- -sum(freqs * log2(freqs + .Machine$double.eps))
    fmax  <- max(freqs)
    data.frame(block_id       = bn,
               n_ind          = n,
               n_haplotypes   = k,
               He             = He,
               Shannon        = Sh,
               freq_dominant  = fmax,
               stringsAsFactors = FALSE)
  })
  do.call(rbind, res)
}


#' Build a Haplotype Dosage Feature Matrix for Genomic Prediction
#'
#' @description
#' Converts the phase-free haplotype strings from
#' \code{\link{extract_haplotypes}} into a numeric dosage matrix suitable for
#' genomic prediction (GBLUP, BayesB, random forests, etc.).
#'
#' For each block, the \code{top_n} most frequent haplotypes are selected as
#' *reference haplotypes*. Each individual then receives an integer dosage
#' count (0, 1, or 2 for diploids) for each reference haplotype, analogous to
#' allele dosage in single-SNP analysis. Rare haplotypes (not in the top
#' \code{top_n}) are collapsed into an \code{"OTHER"} category.
#'
#' @section Why haplotype features?:
#' Single-SNP models assume additive effects and miss multi-locus interactions
#' implicit in haplotype structure. Haplotype dosages capture these
#' interactions without requiring explicit interaction terms, often improving
#' prediction accuracy in structured populations (Calus et al. 2008; de Roos
#' et al. 2009).
#'
#' @param haplotypes List as returned by \code{\link{extract_haplotypes}}.
#' @param top_n Integer. Number of most frequent haplotypes per block to
#'   include as features. Remaining haplotypes are pooled as
#'   \code{"hap_OTHER"}. Default \code{5}.
#' @param missing_string Character. Haplotype strings matching this pattern are
#'   coded as \code{NA} in the output matrix. Default \code{"."}.
#' @param scale_features Logical. If \code{TRUE}, each haplotype dosage column
#'   is centred and scaled (mean 0, sd 1). Recommended when combining blocks
#'   in a single model. Default \code{FALSE}.
#'
#' @return A numeric matrix of dimension
#'   \code{n_individuals × n_features}, where features are haplotype dosages
#'   named \code{"<block_id>__hap<rank>"} (rank 1 = most frequent). Row names
#'   are individual IDs (from the input haplotype list). Columns for the
#'   \code{"OTHER"} category are omitted (it is the implicit reference level).
#'
#' @references
#' Calus MPL et al. (2008) Accuracy of genomic selection using different
#' methods to define haplotypes. \emph{Genetics} \strong{178}(1):553–561.\cr
#' de Roos APW et al. (2009) Linkage disequilibrium and persistence of phase in
#' Holstein–Friesian, Jersey and Angus cattle. \emph{Genetics}
#' \strong{179}(3):1503–1512.
#'
#' @seealso \code{\link{extract_haplotypes}}, \code{\link{compute_haplotype_diversity}}
#'
#' @examples
#' \donttest{
#' set.seed(1)
#' geno <- matrix(sample(0:2, 50 * 40, replace = TRUE), 50, 40)
#' rownames(geno) <- paste0("ind", 1:50)
#' colnames(geno) <- paste0("rs", 1:40)
#' snp_info <- data.frame(SNP = colnames(geno), CHR = "chr1",
#'   POS = seq(1000, by = 5000, length.out = 40))
#' blocks <- data.frame(start.bp = c(1000, 80000), end.bp = c(70000, 200000), CHR = "chr1")
#' haps <- extract_haplotypes(geno, snp_info, blocks)
#' feat <- build_haplotype_feature_matrix(haps, top_n = 3)
#' dim(feat)
#' feat[1:5, 1:min(6, ncol(feat))]
#' }
#'
#' @export
build_haplotype_feature_matrix <- function(
    haplotypes,
    top_n           = 5L,
    missing_string  = ".",
    scale_features  = FALSE
) {
  ind_names   <- names(haplotypes[[1L]])
  block_names <- names(haplotypes)
  all_mats    <- vector("list", length(block_names))

  for (bi in seq_along(block_names)) {
    bn   <- block_names[bi]
    hap  <- haplotypes[[bn]]

    # treat missing
    is_miss              <- grepl(missing_string, hap, fixed = TRUE)
    hap_clean            <- hap
    hap_clean[is_miss]   <- NA_character_

    # frequency table on non-missing
    tbl      <- sort(table(hap_clean[!is_miss]), decreasing = TRUE)
    top_haps <- names(tbl)[seq_len(min(top_n, length(tbl)))]

    # dosage matrix: for each reference haplotype, count occurrences in diploid string
    # Since strings are concatenated diploid codes (length = 2 * n_snps for phased
    # or 1 * n_snps for unphased), we treat exact string match as the haplotype
    # identity and count how many times it occurs across the two "copies" inferred
    # by simple majority assignment.
    # For unphased 0/1/2 codes, a single string per individual is the diploid
    # representation, so dosage is: 2 if individual == top_hap, else 0.
    # The "1" (heterozygous equivalent) arises when blocks are defined
    # across phased data. We default to diploid string identity.

    mat <- matrix(NA_real_, nrow = length(hap_clean), ncol = length(top_haps),
                  dimnames = list(ind_names,
                                  paste0(bn, "__hap", seq_along(top_haps))))
    for (j in seq_along(top_haps)) {
      ref    <- top_haps[j]
      # dosage: 2 = homozygous match, 1 = het (only possible with phased strings),
      # 0 = no match. For unphased strings we only produce 0 or 2.
      dosage          <- ifelse(is.na(hap_clean), NA_real_,
                                ifelse(hap_clean == ref, 2, 0))
      mat[, j]        <- dosage
    }
    all_mats[[bi]] <- mat
  }

  feat_matrix <- do.call(cbind, all_mats)

  if (isTRUE(scale_features)) {
    feat_matrix <- scale(feat_matrix)
    feat_matrix[is.nan(feat_matrix)] <- 0   # zero-variance columns -> 0
  }
  feat_matrix
}
