#' LDxBlocks: Genome-wide LD Block Detection using Kinship-Adjusted LD
#'
#' The \code{LDxBlocks} package provides functions to detect linkage disequilibrium (LD) blocks
#' from genome-wide SNP data, incorporating kinship structure using adjusted \eqn{rV^2} correlations.
#'
#' The main function for users is \code{\link{run_Big_LD_all_chr}}, which performs chromosome-wise
#' LD block detection and outputs genomic intervals representing LD blocks.
#'
#' @section Key Functions:
#' \describe{
#'   \item{\code{\link{run_Big_LD_all_chr}}}{Main wrapper to run \code{Big_LD} by chromosome.}
#'   \item{\code{\link{tune_LD_params}}}{Auto-tunes \code{Big_LD} parameters via grid search to minimize unassigned/forced GWAS marker assignments and produce final LD blocks.}
#'   \item{\code{\link{Big_LD}}}{Core LD block segmentation logic with kinship correction.}
#'   \item{\code{\link{CLQD}}}{Clique detection from \eqn{rV^2} matrix.}
#'   \item{\code{\link{compute_rV2}}}{Computes kinship-adjusted \eqn{rV^2} LD matrix.}
#'   \item{\code{\link{get_V_inv_sqrt}}}{Computes inverse square root of a kinship matrix.}
#' }
#'
#' @docType package
#' @name LDxBlocks-package
#' @keywords package
#'
#' @import data.table
#' @importFrom stats cov median quantile na.omit
#' @importFrom igraph graph_from_adjacency_matrix coreness max_cliques cliques components
#' @importFrom utils globalVariables
"_PACKAGE"

## Silence data.table NSE / scoring column NOTES
utils::globalVariables(c(
  "CHR","start.bp","end.bp","start","end",
  "block_idx",".N","block_name","length_snps","length_bp",
  "n_unassigned","n_forced","n_blocks","penalty_bp"
))
