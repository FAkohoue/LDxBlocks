#' Compute the Inverse Square Root (Whitening) of a Kinship Matrix
#'
#' @description
#' Returns a matrix \eqn{A} that whitens a positive-definite matrix \eqn{V}:
#' \deqn{A\, V\, A^\top = I.}
#' This is used to decorrelate centered genotype data prior to kinship-adjusted LD
#' estimation (\eqn{rV^2}) or mixed-model transforms.
#'
#' @param V A symmetric, numeric matrix of dimension \eqn{n \times n} (typically a kinship matrix).
#' @param rV2method Character string specifying the decomposition:
#'   \itemize{
#'     \item \code{"chol"} (default): Cholesky factorization. If \eqn{V = R^\top R} (R from \code{chol}),
#'           this returns \eqn{A = R^{-1}} (via \code{backsolve}). Fast and stable.
#'     \item \code{"eigen"}: Eigen decomposition. If \eqn{V = Q \Lambda Q^\top}, returns
#'           \eqn{A = Q\,\Lambda^{-1/2}\,Q^\top}. Eigenvalues are thresholded at \code{1e-6}.
#'   }
#'
#' @details
#' The Cholesky-based \eqn{A} is not symmetric but satisfies \eqn{A\,V\,A^\top=I} and is standard for
#' whitening vectors by left-multiplication (\eqn{X^* = A X}). The eigen-based \eqn{A} is symmetric
#' and also satisfies \eqn{A\,V\,A^\top=I}.
#'
#' @return A numeric matrix \eqn{A} of the same dimension as \code{V} such that \eqn{A\,V\,A^\top=I}.
#'
#' @examples
#' \donttest{
#' set.seed(1)
#' G <- matrix(rnorm(100 * 20), nrow = 100, ncol = 20)
#' Gc <- scale(G, center = TRUE, scale = FALSE)
#' kin <- crossprod(Gc) / ncol(G)  # simple kinship
#' A_chol <- get_V_inv_sqrt(kin, rV2method = "chol")
#' A_eig  <- get_V_inv_sqrt(kin, rV2method = "eigen")
#' }
#'
#' @seealso \code{\link{compute_rV2}}, \code{\link{CLQD}}, \code{\link{Big_LD}}, \code{\link{eigen}}, \code{\link{chol}}
#' @export
get_V_inv_sqrt <- function(V, rV2method = c("chol", "eigen")) {
  rV2method <- match.arg(rV2method)
  if (rV2method == "chol") {
    R <- chol(V)                                # V = R^T R (R upper-triangular)
    return(backsolve(R, diag(nrow(V))))         # A = R^{-1}; A V A^T = I
  } else {
    eig  <- eigen(V, symmetric = TRUE)
    vals <- pmax(eig$values, 1e-6)              # stabilize tiny eigenvalues
    return(eig$vectors %*% diag(1 / sqrt(vals)) %*% t(eig$vectors))  # A = Q Λ^{-1/2} Q^T
  }
}
