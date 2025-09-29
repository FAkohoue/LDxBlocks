#' Compute Kinship-Adjusted Squared Correlations (\eqn{rV^2})
#'
#' @description
#' Computes the squared correlation matrix \eqn{rV^2} between SNPs from a kinship-adjusted genotype matrix.
#' Optionally rounds the result to \code{digits} decimal places to reduce floating-point jitter.
#'
#' @param X A numeric matrix of dimension \eqn{n \times p}, where rows are
#'   individuals and columns are SNPs. This should be a kinship-adjusted,
#'   mean-centered genotype matrix (i.e., \eqn{V^{-1/2}} times the centered
#'   genotype matrix).
#'
#' @param digits Integer or \code{NULL}. If non-\code{NULL}, round the \eqn{rV^2} matrix to this many
#'   decimal places for reproducibility. Default \code{NULL} (no rounding).
#'
#' @return A symmetric \eqn{p \times p} numeric matrix with zeros on the diagonal.
#' @seealso \code{\link{cov}}, \code{\link{scale}}, \code{\link{get_V_inv_sqrt}}
#' @export
compute_rV2 <- function(X, digits = NULL) {
  cov_mat <- stats::cov(X, use = "pairwise.complete.obs")
  diag_cov <- diag(cov_mat)

  # rV²: squared correlation
  rV2 <- (cov_mat^2) / outer(diag_cov, diag_cov, `*`)

  # Clean / format
  rV2[is.na(rV2)] <- 0
  diag(rV2) <- 0
  rV2 <- (rV2 + t(rV2)) / 2  # enforce symmetry

  if (!is.null(digits)) {
    rV2 <- round(rV2, digits = as.integer(digits))
  }
  rV2
}

