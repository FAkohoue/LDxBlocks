#' LDxBlocks: Genome-wide LD Block Detection using Kinship-Adjusted LD
#'
#' @description
#' The `LDxBlocks` package provides functions to detect linkage disequilibrium (LD) blocks 
#' from genome-wide SNP data, incorporating kinship structure using adjusted rV² correlations.
#'
#' The main function for users is \code{\link{run_Big_LD_all_chr}}, which performs chromosome-wise
#' LD block detection and outputs genomic intervals representing LD blocks.
#'
#' @section Key Functions:
#' \describe{
#'   \item{\code{\link{run_Big_LD_all_chr}}}{Main wrapper function to run Big_LD by chromosome}
#'   \item{\code{\link{Big_LD}}}{Core LD block segmentation logic with kinship correction}
#'   \item{\code{\link{CLQD}}}{Clique detection from rV² matrix}
#'   \item{\code{\link{compute_rV2}}}{Computes kinship-adjusted rV² LD matrix}
#'   \item{\code{\link{get_V_inv_sqrt}}}{Computes inverse square root of kinship matrix}
#' }
#'
#' @author
#' Félicien Akohoue
#'
#' @keywords package
"_PACKAGE"

