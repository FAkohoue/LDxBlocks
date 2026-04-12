#' LDxBlocks: Genome-Wide LD Block Detection and Haplotype Analysis
#'
#' @description
#' Fast, scalable linkage disequilibrium (LD) block detection for genome-wide
#' SNP data using a C++/Armadillo computational core and two LD metrics:
#' standard r^2 and kinship-adjusted rV^2. Includes haplotype reconstruction,
#' diversity metrics, prediction feature matrices, and multi-format genotype I/O.
#'
#' @section Choosing between r^2 and rV^2:
#' \describe{
#'   \item{r^2 (\code{method = "r2"})}{Standard squared Pearson correlation.
#'     No kinship matrix required. C++ kernel via RcppArmadillo + OpenMP.
#'     Suitable for unstructured populations and large datasets (10 M+ markers).}
#'   \item{rV^2 (\code{method = "rV2"})}{Kinship-whitened squared correlation
#'     (Mangin et al. 2012, Heredity 108:285-291). Corrects LD inflation from relatedness. For livestock,
#'     inbred lines, or family-based cohorts with < 200 k markers per chromosome.}
#' }
#'
#' @section Pipeline overview:
#' \enumerate{
#'   \item Read genotype data: \code{\link{read_geno}} - auto-detects CSV,
#'     HapMap, VCF, GDS, BED, or plain matrix.
#'   \item Detect LD blocks chromosome-wise: \code{\link{run_Big_LD_all_chr}}.
#'   \item Optionally auto-tune parameters: \code{\link{tune_LD_params}}.
#'   \item Reconstruct haplotypes: \code{\link{extract_haplotypes}}.
#'   \item Compute diversity metrics: \code{\link{compute_haplotype_diversity}}.
#'   \item Build prediction feature matrix:
#'     \code{\link{build_haplotype_feature_matrix}}.
#' }
#'
#' @section Example data:
#' \describe{
#'   \item{\code{\link{ldx_geno}}}{120 x 200 genotype matrix, 3 chromosomes,
#'     9 simulated LD blocks.}
#'   \item{\code{\link{ldx_snp_info}}}{SNP metadata (SNP, CHR, POS, REF, ALT).}
#'   \item{\code{\link{ldx_blocks}}}{Reference block table for \code{ldx_geno}.}
#'   \item{\code{\link{ldx_gwas}}}{20 toy GWAS markers for tuning demos.}
#' }
#'
#' @references
#' Kim S-A et al. (2018) Bioinformatics 34(4):588-596.
#' \doi{10.1093/bioinformatics/btx609}
#'
#' Mangin B et al. (2012) Heredity 108(3):285-291.
#' \doi{10.1038/hdy.2011.73}
#'
#' VanRaden PM (2008) J. Dairy Sci. 91(11):4414-4423.
#' \doi{10.3168/jds.2007-0980}
#'
#' @docType package
#' @name LDxBlocks-package
#' @aliases LDxBlocks
#' @keywords package
#'
#' @import Rcpp
#' @importFrom data.table fread fwrite rbindlist setnames setorder :=
#' @importFrom stats cor median quantile na.omit setNames
#' @importFrom igraph graph_from_adjacency_matrix coreness max_cliques cliques components
#'   cliques components
#' @importFrom utils globalVariables read.table
"_PACKAGE"

utils::globalVariables(c(
  "CHR", "start.bp", "end.bp", "start", "end",
  "block_idx", ".N", "block_name", "length_snps", "length_bp",
  "n_unassigned", "n_forced", "n_blocks", "penalty_bp",
  "REF", "ALT", "SNP", "POS", "..samp_cols", ".",
  "start_x", "end_x", "chr_int"
))
