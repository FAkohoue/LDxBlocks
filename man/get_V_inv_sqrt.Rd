\name{get_V_inv_sqrt}
\alias{get_V_inv_sqrt}
\title{
Compute the Inverse Square Root of a Kinship Matrix
}
\description{
Computes the matrix square root of the inverse of a symmetric, positive-definite matrix, such as a kinship matrix. This transformation is used to decorrelate genotype data before downstream analyses including LD estimation using kinship-adjusted metrics such as rV².
}
\usage{
get_V_inv_sqrt(V, rV2method = c("chol", "eigen"))
}
\arguments{
  \item{V}{A symmetric, numeric matrix of dimension \eqn{n \times n}. Typically a kinship matrix estimated from genotype data (e.g., using VanRaden's method).}

  \item{rV2method}{Decomposition method to compute \eqn{V^{-1/2}}. One of \code{"chol"} for Cholesky decomposition (fast, default), or \code{"eigen"} for eigenvalue decomposition (more general).}
}
\details{
This function returns the matrix square root of the inverse of \code{V}, satisfying:

\deqn{V^{-1/2} \cdot V \cdot V^{-1/2} = I}

\strong{Method "chol"} uses the Cholesky decomposition:
\itemize{
  \item Computes \eqn{L} such that \eqn{V = L L^\top}.
  \item Returns \eqn{V^{-1/2} = L^{-1}}.
}

\strong{Method "eigen"} uses the eigen decomposition:
\itemize{
  \item Computes \eqn{V = Q \Lambda Q^\top}, then \eqn{V^{-1/2} = Q \Lambda^{-1/2} Q^\top}.
  \item Stabilizes eigenvalues below 1e-6 to avoid inversion errors.
}

This inverse root matrix is used for kinship-adjusted transformation of genotype data before computing \eqn{rV^2} or fitting LMMs.
}
\value{
A numeric matrix of the same dimensions as \code{V}, representing \eqn{V^{-1/2}}.
}

\keyword{matrix}
\keyword{kinship}
\keyword{LD}

\references{Félicien Akohoue
%% ~put references to the literature/web site here ~
}
\author{Félicien Akohoue
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{\code{\link{compute_rV2}}, \code{\link{CLQD}}, \code{\link{Big_LD}}}

\examples{
set.seed(1)
G <- matrix(rnorm(100 * 20), nrow = 100, ncol = 20)
G_centered <- scale(G, center = TRUE, scale = FALSE)
kinship <- crossprod(G_centered) / ncol(G)

# Using Cholesky
V_inv_sqrt_chol <- get_V_inv_sqrt(kinship, rV2method = "chol")

# Using Eigenvalue
V_inv_sqrt_eig <- get_V_inv_sqrt(kinship, rV2method = "eigen")

# Transform genotype
G_adjusted <- V_inv_sqrt_chol \%*\% G_centered

##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or standard data sets, see data().

## The function is currently defined as

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
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
