#' @title Compute Kinship-Adjusted Squared Correlations (rV²)
#'
#' @description
#' Computes the squared correlation matrix \eqn{rV^2} between SNPs from a kinship-adjusted genotype matrix.
#' This matrix quantifies the linkage disequilibrium (LD) between SNPs after removing individual-level covariance
#' using a transformation such as \eqn{V^{-1/2} \cdot G}, where \eqn{V} is the kinship matrix.
#'
#' @param X A numeric matrix of dimension \eqn{n \times p}, where rows are individuals and columns are SNPs.
#'          This should be a kinship-adjusted and mean-centered genotype matrix.
#'
#' @details
#' The function computes the squared Pearson correlation for each SNP pair using:
#' \deqn{rV^2_{ij} = \frac{\text{Cov}(X_i, X_j)^2}{\text{Var}(X_i) \cdot \text{Var}(X_j)}}
#'
#' The resulting matrix is symmetric, with zero on the diagonal. Any \code{NA} values due to missing data
#' are replaced with zero. Floating-point asymmetries are corrected by enforcing perfect symmetry.
#'
#' This measure is used in LD block detection pipelines where individual relatedness (kinship)
#' can inflate traditional \eqn{r^2} LD estimates. Using \eqn{rV^2} offers a more accurate LD signal
#' for structured populations.
#'
#' @return A symmetric numeric matrix of dimension \eqn{p \times p}, where \eqn{p} is the number of SNPs.
#'         The matrix contains squared correlations between columns of \code{X}, with zeros on the diagonal.
#'
#' @examples
#' # Simulate adjusted genotype matrix
#' set.seed(1)
#' G <- matrix(rnorm(100 * 50), nrow = 100, ncol = 50)
#' G_centered <- scale(G, center = TRUE, scale = FALSE)
#'
#' # Simulate kinship and adjust genotypes
#' kin <- crossprod(G_centered) / ncol(G)
#' L_inv <- backsolve(chol(kin), diag(nrow(kin)))  # fast V^{-1/2}
#' G_adj <- L_inv %*% G_centered
#'
#' # Compute rV²
#' rV2_mat <- compute_rV2(G_adj)
#'
#' # Check properties
#' dim(rV2_mat)
#' range(rV2_mat)
#'
#' @author
#' Félicien Akohoue
#'
#' @seealso \code{\link{cov}}, \code{\link{scale}}, \code{\link{get_V_inv_sqrt}}
#'
#' @export
#'
compute_rV2 <- function(X) {
  # Pairwise covariance
  cov_mat <- cov(X, use = "pairwise.complete.obs")

  # Variances of each SNP (on diagonal)
  diag_cov <- diag(cov_mat)

  # rV²: squared correlation
  rV2 <- (cov_mat^2) / outer(diag_cov, diag_cov, `*`)

  # Clean and format
  rV2[is.na(rV2)] <- 0
  diag(rV2) <- 0
  rV2 <- (rV2 + t(rV2)) / 2  # ensure perfect symmetry

  return(rV2)
}
