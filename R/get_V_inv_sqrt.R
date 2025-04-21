#' @title Compute the Inverse Square Root of a Kinship Matrix
#'
#' @name get_V_inv_sqrt
#' @alias get_V_inv_sqrt
#'
#' @description
#' Computes the inverse square root of a symmetric positive-definite matrix, commonly used to adjust genotype data for kinship.
#' This transformation removes individual-level covariance from genotype matrices prior to downstream analyses such as kinship-adjusted LD estimation (\eqn{rV^2}) or linear mixed model fitting.
#'
#' @param V A symmetric, numeric matrix of dimension \eqn{n \times n}. Typically a kinship matrix derived from additive genetic relatedness (e.g., using the VanRaden method).
#' @param rV2method Character string specifying the matrix decomposition method to use. Options are:
#'   \itemize{
#'     \item \code{"chol"}: Cholesky decomposition (default), fast and numerically stable.
#'     \item \code{"eigen"}: Eigenvalue decomposition, more general but potentially slower.
#'   }
#'
#' @details
#' The function returns \eqn{V^{-1/2}}, satisfying:
#' \deqn{V^{-1/2} \cdot V \cdot V^{-1/2} = I}
#'
#' This matrix is used to transform centered genotype matrices prior to LD calculations involving kinship-adjusted squared correlation (rV²) or mixed model inference.
#'
#' \strong{Cholesky Decomposition ("chol"):}
#' \itemize{
#'   \item Computes the Cholesky factor \eqn{L} such that \eqn{V = L L^\top}.
#'   \item Then returns \eqn{V^{-1/2} = L^{-1}}.
#' }
#'
#' \strong{Eigen Decomposition ("eigen"):}
#' \itemize{
#'   \item Computes \eqn{V = Q \Lambda Q^\top}, where \eqn{\Lambda} is diagonal.
#'   \item Then returns \eqn{V^{-1/2} = Q \Lambda^{-1/2} Q^\top}.
#'   \item Small eigenvalues (less than 1e-6) are stabilized to avoid numerical artifacts.
#' }
#'
#' @return A numeric matrix of the same dimension as \code{V}, representing the inverse square root \eqn{V^{-1/2}}.
#'
#' @examples
#' # Simulate genotype and kinship matrix
#' set.seed(1)
#' G <- matrix(rnorm(100 * 20), nrow = 100, ncol = 20)
#' G_centered <- scale(G, center = TRUE, scale = FALSE)
#' kinship <- crossprod(G_centered) / ncol(G)
#'
#' # Compute V^(-1/2) using Cholesky (default)
#' V_inv_sqrt_chol <- get_V_inv_sqrt(kinship, rV2method = "chol")
#'
#' # Compute V^(-1/2) using Eigen decomposition
#' V_inv_sqrt_eig <- get_V_inv_sqrt(kinship, rV2method = "eigen")
#'
#' # Use result to adjust genotype matrix
#' G_adjusted <- V_inv_sqrt_chol %*% G_centered
#'
#' @seealso \code{\link{compute_rV2}}, \code{\link{CLQD}}, \code{\link{Big_LD}}, \code{\link{eigen}}, \code{\link{chol}}
#'
#' @author
#' Félicien Akohoue
#'
#' @keywords matrix decomposition kinship normalization LD
#'
#' @export
#' 
get_V_inv_sqrt <- function(V, rV2method = c("chol", "eigen")) {
  rV2method <- match.arg(rV2method)
  if (rV2method == "chol") {
    L <- chol(V)
    return(backsolve(L, diag(nrow(V))))  # Fast V^{-1/2}
  } else {
    eig <- eigen(V, symmetric = TRUE)
    vals <- pmax(eig$values, 1e-6)  # stabilize small eigenvalues
    return(eig$vectors %*% diag(1 / sqrt(vals)) %*% t(eig$vectors))  # General V^{-1/2}
  }
}
