\name{compute_rV2}
\alias{compute_rV2}
\title{Compute Kinship-Adjusted Squared Correlations (rV²)}
\description{
Computes the squared correlation matrix \eqn{rV^2} between SNPs from a kinship-adjusted genotype matrix.
This measure quantifies linkage disequilibrium (LD) after removing covariance induced by population structure or relatedness.
}
\usage{
compute_rV2(X)
}
\arguments{
  \item{X}{A numeric matrix of dimension \eqn{n \times p}, where rows are individuals and columns are SNPs. 
  This should be a kinship-adjusted and mean-centered genotype matrix, typically obtained by transforming the genotype matrix 
  \eqn{G} using \eqn{V^{-1/2} \cdot G}, where \eqn{V} is the kinship matrix.}
}
\details{
This function estimates the pairwise squared correlation (similar to \eqn{r^2}) between SNPs after correcting for kinship structure.
The formula used is:

\deqn{rV^2_{ij} = \frac{\mathrm{Cov}(X_i, X_j)^2}{\mathrm{Var}(X_i) \cdot \mathrm{Var}(X_j)}}

The result is a symmetric matrix of \eqn{rV^2} values, with diagonal elements set to zero. 
This method is used in LD block detection pipelines where individual-level covariance (kinship) inflates traditional LD estimates.

\strong{Data handling:}
\itemize{
  \item Missing values are handled using \code{use = "pairwise.complete.obs"} in the covariance computation.
  \item Any \code{NA} values in the result are replaced with 0.
  \item Floating point asymmetry is corrected by enforcing symmetry with \code{(rV2 + t(rV2)) / 2}.
}
}
\value{
A symmetric numeric matrix of size \eqn{p \times p}, where \eqn{p} is the number of SNPs.
Each entry contains the squared correlation \eqn{rV^2} between two SNPs after kinship correction.
}

\keyword{LD}
\keyword{genetics}
\keyword{correlation}


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

\seealso{
\code{\link{cov}}, \code{\link{scale}}, \code{\link{get_V_inv_sqrt}}
}

\examples{
# Simulate adjusted genotype matrix
set.seed(1)
G <- matrix(rnorm(100 * 50), nrow = 100, ncol = 50)
G_centered <- scale(G, center = TRUE, scale = FALSE)

# Simulate kinship and adjust genotypes
kin <- crossprod(G_centered) / ncol(G)
L_inv <- backsolve(chol(kin), diag(nrow(kin)))  # fast V^{-1/2}
G_adj <- L_inv \%*\% G_centered

# Compute rV²
rV2_mat <- compute_rV2(G_adj)

# Check properties
dim(rV2_mat)
range(rV2_mat)

##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or standard data sets, see data().

## The function is currently defined as

compute_rV2 <- function(X) {
  # Pairwise covariance
  cov_mat <- cov(X, use = "pairwise.complete.obs")
  
  # Variances of each SNP (on diagonal)
  diag_cov <- diag(cov_mat)
  
  # rV²: squared correlation
  rV2 <- (cov_mat^2) / (diag_cov %*% t(diag_cov))
  
  # Clean and format
  rV2[is.na(rV2)] <- 0
  diag(rV2) <- 0
  rV2 <- (rV2 + t(rV2)) / 2  # ensure perfect symmetry
  
  return(rV2)
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
