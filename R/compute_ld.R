# -----------------------------------------------------------------------------
# compute_ld.R  -  Unified LD computation: standard r^2 or kinship-adjusted rV^2
# -----------------------------------------------------------------------------
# Design:
#   compute_ld()      -- dispatcher: routes to r^2 or rV^2 path
#   compute_r2()      -- standard r^2 via C++ (fast, no kinship)
#   compute_rV2()     -- kinship-adjusted rV^2 via C++ (structured populations)
#   get_V_inv_sqrt()  -- whitening factor (unchanged, needed for rV^2)
#   prepare_geno()    -- centre, optionally whiten; returns list(X, V_inv_sqrt)
# -----------------------------------------------------------------------------


#' Compute LD Matrix: Standard r^2 or Kinship-Adjusted rV^2
#'
#' @description
#' Unified dispatcher. With \code{method = "r2"} (default and recommended
#' for large datasets) it calls the C++ Armadillo back-end for fast standard
#' squared Pearson correlations. With \code{method = "rV2"} it applies kinship
#' whitening first and then the same C++ correlation kernel.
#'
#' @param X NumericMatrix (individuals x SNPs).
#'   For \code{method = "r2"}: raw or mean-centred genotypes (0/1/2).
#'   For \code{method = "rV2"}: the pre-whitened matrix
#'   \eqn{V^{-1/2} \tilde{G}} produced by \code{prepare_geno()}.
#' @param method Character. \code{"r2"} (default) or \code{"rV2"}.
#' @param digits Integer. Rounding precision. \code{-1} skips rounding.
#' @param n_threads Integer. OpenMP threads for the C++ kernel. Default 1.
#'   Use \code{parallel::detectCores()} to choose automatically.
#'
#' @return Symmetric p x p numeric matrix, diagonal 0, values in [0, 1].
#'
#' @seealso \code{\link{compute_r2}}, \code{\link{compute_rV2}},
#'   \code{\link{prepare_geno}}, \code{Big_LD()}
#'
#' @examples
#' \dontrun{
#' # Internal function -- use compute_r2() or compute_rV2() instead
#' set.seed(1)
#' G <- matrix(sample(0:2, 60 * 20, replace = TRUE), 60, 20)
#' ld_r2  <- LDxBlocks:::compute_ld(G, method = "r2")
#' range(ld_r2)
#' }
#'
#' @keywords internal
compute_ld <- function(X, method = c("r2", "rV2"), digits = -1L, n_threads = 1L) {
  method <- match.arg(method)
  if (!is.matrix(X)) X <- as.matrix(X)
  if (method == "r2") {
    compute_r2_cpp(X, digits = as.integer(digits), n_threads = as.integer(n_threads))
  } else {
    compute_rV2_cpp(X, digits = as.integer(digits), n_threads = as.integer(n_threads))
  }
}


#' Compute Standard r^2 LD Matrix
#'
#' @description
#' Fast C++/Armadillo implementation of the standard pairwise squared Pearson
#' correlation (r^2) for a window of SNP columns. Missing genotypes (NA) are
#' mean-imputed per column before computation. This is the default LD metric
#' in \pkg{LDxBlocks} and is 10-50x faster than \code{stats::cor()}.
#'
#' @section When to use r^2 vs rV^2:
#' \describe{
#'   \item{r^2}{Use for large unstructured datasets (> 500 k markers), random
#'     mating populations, or whenever speed matters. The standard estimator
#'     is inflated in related populations (i.e. will over-estimate LD) but
#'     this usually leads to slightly more conservative (larger) blocks rather
#'     than catastrophically wrong ones.}
#'   \item{rV^2}{Use for highly structured / related populations (livestock,
#'     inbred lines, family-based human cohorts) where kinship inflation
#'     would meaningfully distort block boundaries. Requires computing and
#'     inverting the GRM - prohibitive beyond ~5 k individuals.}
#' }
#'
#' @param X NumericMatrix (individuals x SNPs). Values 0/1/2. NA allowed.
#' @param digits Integer. Rounding decimal places. \code{-1} (default) = none.
#' @param n_threads Integer. OpenMP threads. Default 1.
#'
#' @return Symmetric p x p NumericMatrix, diagonal 0, values in [0, 1].
#'
#' @seealso \code{\link{compute_rV2}}, \code{\link{compute_ld}}
#'
#' @examples
#' set.seed(1)
#' G <- matrix(sample(0:2, 80 * 30, replace = TRUE), 80, 30)
#' r2 <- compute_r2(G, digits = 6L)
#' range(r2)   # [0, 1]
#'
#' @export
compute_r2 <- function(X, digits = -1L, n_threads = 1L) {
  if (!is.matrix(X)) X <- as.matrix(X)
  compute_r2_cpp(X, digits = as.integer(digits), n_threads = as.integer(n_threads))
}


#' Compute Kinship-Adjusted rV^2 LD Matrix
#'
#' @description
#' Computes the kinship-adjusted squared correlation (rV^2) for a
#' pre-whitened genotype matrix \eqn{X = V^{-1/2} \tilde{G}}.
#' The whitening step must be performed first via \code{\link{prepare_geno}}.
#'
#' Mathematically identical to \code{\link{compute_r2}} applied to \eqn{X}
#' because the whitening already removes kinship structure. The C++ kernel
#' is the same - the distinction is purely in the preparation step.
#'
#' @param X NumericMatrix (n x p). Pre-whitened, mean-centred genotype matrix.
#' @param digits Integer. Rounding. \code{-1} = none.
#' @param n_threads Integer. OpenMP threads.
#'
#' @return Symmetric p x p NumericMatrix, diagonal 0.
#'
#' @references Kim S-A et al. (2018) GENETICS 209(3):855-868.
#'
#' @seealso \code{\link{compute_r2}}, \code{\link{prepare_geno}}
#'
#' @export
compute_rV2 <- function(X, digits = -1L, n_threads = 1L) {
  if (!is.matrix(X)) X <- as.matrix(X)
  compute_rV2_cpp(X, digits = as.integer(digits), n_threads = as.integer(n_threads))
}


#' Compute the Inverse Square Root (Whitening Factor) of a Kinship Matrix
#'
#' @description
#' Returns \eqn{A} such that \eqn{A V A^\top = I}. Used by
#' \code{\link{prepare_geno}} when \code{method = "rV2"}.
#'
#' @param V Symmetric positive-definite matrix (n x n). Typically a VanRaden
#'   GRM after bending/tuning.
#' @param method Character. \code{"chol"} (default, faster) or
#'   \code{"eigen"} (more robust for near-singular matrices).
#'
#' @return Numeric matrix A (n x n).
#'
#' @export
get_V_inv_sqrt <- function(V, method = c("chol", "eigen")) {
  method <- match.arg(method)
  if (method == "chol") {
    R <- chol(V)
    backsolve(R, diag(nrow(V)))
  } else {
    eig  <- eigen(V, symmetric = TRUE)
    vals <- pmax(eig$values, 1e-6)
    eig$vectors %*% diag(1 / sqrt(vals)) %*% t(eig$vectors)
  }
}


#' Prepare Genotype Matrix for LD Computation
#'
#' @description
#' Central preparation function called inside \code{Big_LD()} and
#' \code{\link{run_Big_LD_all_chr}}. Depending on \code{method}:
#'
#' \describe{
#'   \item{\code{"r2"}}{Mean-centres each SNP column. No kinship needed.
#'     Returns the centred matrix directly. Fast - O(np).}
#'   \item{\code{"rV2"}}{Computes VanRaden GRM via AGHmatrix, tunes it
#'     via ASRgenomics, inverts via Cholesky/eigen, and left-multiplies
#'     the centred genotype matrix. Returns the whitened matrix and the
#'     whitening factor for use in \code{appendSGTs}.}
#' }
#'
#' @param geno Numeric matrix (individuals x SNPs), 0/1/2.
#' @param method Character. \code{"r2"} or \code{"rV2"}.
#' @param kin_method Character. Whitening decomposition for rV^2:
#'   \code{"chol"} or \code{"eigen"}.
#' @param verbose Logical. Print progress.
#'
#' @return Named list:
#'   \item{adj_geno}{n x p numeric matrix ready for \code{compute_ld()}.}
#'   \item{V_inv_sqrt}{n x n whitening matrix, or \code{NULL} for r^2.}
#'
#' @export
prepare_geno <- function(
    geno,
    method     = c("r2", "rV2"),
    kin_method = c("chol", "eigen"),
    verbose    = FALSE
) {
  method     <- match.arg(method)
  kin_method <- match.arg(kin_method)

  # Ensure individual IDs (required by AGHmatrix)
  ids <- rownames(geno)
  if (is.null(ids)) {
    ids <- sprintf("ind%04d", seq_len(nrow(geno)))
    rownames(geno) <- ids
  }

  geno_centered <- scale(geno, center = TRUE, scale = FALSE)

  if (method == "r2") {
    return(list(adj_geno = geno_centered, V_inv_sqrt = NULL))
  }

  # rV^2 path: compute GRM, whiten
  if (!requireNamespace("AGHmatrix", quietly = TRUE))
    stop("AGHmatrix is required for method = 'rV2'. ",
         "Install with: install.packages('AGHmatrix')", call. = FALSE)
  if (!requireNamespace("ASRgenomics", quietly = TRUE))
    stop("ASRgenomics is required for method = 'rV2'. ",
         "Install with: install.packages('ASRgenomics')", call. = FALSE)
  if (isTRUE(verbose)) cat("[prepare_geno] Computing VanRaden GRM...\n")
  kin <- AGHmatrix::Gmatrix(geno, method = "VanRaden")
  dimnames(kin) <- list(ids, ids)

  if (isTRUE(verbose)) cat("[prepare_geno] Tuning GRM (bend + rcn)...\n")
  kin_tuned  <- ASRgenomics::G.tuneup(kin, bend = TRUE, rcn = TRUE)$Gb
  V_inv_sqrt <- get_V_inv_sqrt(kin_tuned, method = kin_method)

  adj_geno <- V_inv_sqrt %*% geno_centered
  list(adj_geno = adj_geno, V_inv_sqrt = V_inv_sqrt)
}
