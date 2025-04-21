\name{CLQD}
\alias{CLQD}
\title{Clique-Based LD Block Detection Using Kinship-Adjusted Squared Correlations (rV²)}
\description{
Partitions SNPs into linkage disequilibrium (LD) blocks (cliques) based on the kinship-adjusted squared correlation matrix \eqn{rV^2}.
This function builds cliques from SNPs that exceed a correlation threshold, allowing for physical distance-based splitting and clique optimization heuristics.
}
\usage{
CLQD(
  subgeno, subSNPinfo, adj_subgeno,
  CLQcut = 0.5,
  clstgap = 40000,
  CLQmode = c("Density", "Maximal"),
  codechange = FALSE,
  checkLargest = FALSE,
  split = FALSE,
  verbose = TRUE
)
}
\arguments{
  \item{subgeno}{A numeric matrix of raw genotype values (individuals × SNPs).}
  \item{subSNPinfo}{A data frame with SNP metadata. The first column must contain SNP IDs (rsID), and the second column must contain base-pair positions.}
  \item{adj_subgeno}{A numeric matrix containing kinship-adjusted, mean-centered genotype values. Typically computed as \eqn{V^{-1/2} \cdot G}, where \eqn{G} is the centered genotype matrix and \eqn{V} is a kinship matrix.}
  \item{CLQcut}{Threshold for the squared correlation (\eqn{rV^2}) to form edges in the LD graph. Default is 0.5.}
  \item{clstgap}{Maximum allowed physical distance (in base pairs) between consecutive SNPs within a clique. If exceeded, the clique is split. Default is 40000.}
  \item{CLQmode}{Clique scoring method used to prioritize cliques. Options are \code{"Density"} (default), which prefers tightly clustered SNPs, or \code{"Maximal"}, which favors cliques with more members.}
  \item{codechange}{Logical. If \code{TRUE}, attempts to flip correlation signs within cliques to maximize clique inclusion (less common). Default is \code{FALSE}.}
  \item{checkLargest}{Logical. If \code{TRUE}, enables a fast pre-processing step using dense-core decomposition for very dense LD matrices (SNPs > 500). Default is \code{FALSE}.}
  \item{split}{Logical. If \code{TRUE}, split cliques that span long genomic distances (> \code{clstgap}). Default is \code{FALSE}.}
  \item{verbose}{Logical. Print diagnostic and progress messages. Default is \code{TRUE}.}
}
\details{
This function applies a graph-based approach to identify LD blocks (cliques) using the kinship-adjusted squared correlation matrix (\eqn{rV^2}) as a similarity measure.

The \code{igraph} package is used to construct a graph from the adjacency matrix defined by SNPs with \eqn{rV^2} values above \code{CLQcut}.
Maximal cliques are identified from this graph and optionally split based on physical distances between SNPs (\code{clstgap}).
If \code{checkLargest} is enabled and the SNP count exceeds 500, a fast heuristic reduces the search space.

The final result is a vector assigning SNPs to block IDs, where \code{NA} indicates singletons (SNPs not clustered into any clique).
}
\value{
A numeric vector of length equal to the number of SNPs in \code{subgeno}. Each element contains the ID of the assigned LD block, or \code{NA} if the SNP is a singleton.
}

\keyword{LD}
\keyword{genetics}
\keyword{clustering}
\keyword{linkage}


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
\code{\link{compute_rV2}}, \code{\link{get_V_inv_sqrt}}, \code{\link[igraph]{max_cliques}}, \code{\link{Big_LD}}
}

\examples{
# Example assumes kinship-adjusted genotype matrix is already computed
set.seed(123)
G <- matrix(rnorm(100 * 50), nrow = 100)
G_centered <- scale(G, center = TRUE, scale = FALSE)
kin <- crossprod(G_centered) / ncol(G)
L_inv <- backsolve(chol(kin), diag(nrow(kin)))
G_adj <- L_inv \%*\% G_centered

snp_info <- data.frame(SNP = paste0("SNP", 1:50), POS = sort(sample(1e6, 50)))
block_assignments <- CLQD(G, snp_info, G_adj, CLQcut = 0.4, clstgap = 100000, CLQmode = "Density")
table(block_assignments)



##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or standard data sets, see data().

## The function is currently defined as

CLQD <- function(
    subgeno, subSNPinfo, adj_subgeno,
    CLQcut = 0.5,
    clstgap = 40000,
    CLQmode = c("Density", "Maximal"),
    codechange = FALSE,
    checkLargest = FALSE,
    split = FALSE,
    verbose = TRUE
) {
  # Ensure mode is singular
  CLQmode <- match.arg(CLQmode)
  
  # Get SNP base-pair positions
  SNPbps <- as.numeric(subSNPinfo[[2]])
  
  # Clique scoring
  clique_score <- function(bp_positions) {
    if (CLQmode == "Maximal") {
      return(length(bp_positions))
    } else { # "Density"
      return(length(bp_positions) / (diff(range(bp_positions)) / 1000 + 1))  # avoid /0
    }
  }
  
  # Sign correction for cliques
  code_sign_opt <- function(vt, corr_mat) {
    rin <- corr_mat[vt, vt]
    nr <- length(vt)
    code <- rep(1, nr)
    change <- TRUE
    iter <- 0
    max_iter <- 2^nr
    
    while (change && iter < max_iter) {
      change <- FALSE
      neg_counts <- rowSums(rin < 0)
      idx <- which.max(neg_counts)
      if (neg_counts[idx] > (nr - 1) / 2) {
        rin[idx, ] <- -rin[idx, ]
        rin[, idx] <- -rin[, idx]
        code[idx] <- 1 - code[idx]
        change <- TRUE
      }
      iter <- iter + 1
    }
    return(list(matrix = rin, code = code))
  }
  
  # Physically split cliques by SNP gaps
  split_cliques <- function(cliques_bp, gapdist) {
    cliques_bp <- lapply(cliques_bp, function(x) sort(as.numeric(na.omit(x))))
    final <- list()
    
    while (length(cliques_bp) > 0) {
      newlist <- list()
      for (c in cliques_bp) {
        gaps <- diff(c)
        if (any(gaps > gapdist)) {
          split_at <- which(gaps > gapdist)[1]
          front <- c[1:split_at]
          rest <- c[(split_at + 1):length(c)]
          final <- c(final, list(front))
          if (length(rest) > 1) newlist <- c(newlist, list(rest))
        } else {
          final <- c(final, list(c))
        }
      }
      cliques_bp <- newlist
    }
    return(final)
  }
  
  # Begin computation
  rV2_mat <- compute_rV2(adj_subgeno)
  adj_bin <- ifelse(rV2_mat >= CLQcut, 1, 0)
  diag(adj_bin) <- 0
  
  binvector <- rep(NA, ncol(rV2_mat))
  re_SNPbps <- SNPbps
  bin_id <- 1
  
  # Dense decomposition (optional)
  if (checkLargest && ncol(rV2_mat) >= 500) {
    repeat {
      g <- igraph::graph_from_adjacency_matrix(adj_bin, mode = "undirected", diag = FALSE)
      comp <- igraph::components(g)
      largest <- which.max(comp$csize)
      comp_idx <- which(comp$membership == largest)
      sub_mat <- rV2_mat[comp_idx, comp_idx]
      
      cg <- igraph::graph_from_adjacency_matrix(sub_mat, mode = "undirected", diag = FALSE)
      core_vals <- igraph::coreness(cg)
      if (median(core_vals) > 80 && max(core_vals) > 100) {
        high_deg <- which(rowSums(adj_bin) >= quantile(rowSums(adj_bin), 0.7))
        densest <- high_deg[which.max(sapply(high_deg, function(v) {
          neighbors <- which(adj_bin[v, ] > 0)
          nb_mat <- adj_bin[neighbors, neighbors]
          sum(nb_mat) / (length(neighbors) * (length(neighbors) - 1))
        }))]
        idx_to_remove <- which(adj_bin[densest, ] > 0)
        binvector[match(re_SNPbps[idx_to_remove], SNPbps)] <- bin_id
        bin_id <- bin_id + 1
        
        keep_idx <- setdiff(seq_along(re_SNPbps), idx_to_remove)
        re_SNPbps <- re_SNPbps[keep_idx]
        rV2_mat <- rV2_mat[keep_idx, keep_idx]
        adj_bin <- adj_bin[keep_idx, keep_idx]
        
        if (length(re_SNPbps) < 500) break
      } else {
        break
      }
    }
  }
  
  if (verbose) message("End of pre-processing")
  
  g <- igraph::graph_from_adjacency_matrix(adj_bin, mode = "undirected", diag = FALSE)
  cliques <- igraph::max_cliques(g, min = 2)
  
  if (verbose) {
    cat("Detected cliques:", length(cliques), "\n")
    print(table(sapply(cliques, length)))
  }
  
  bp_cliques <- lapply(cliques, function(x) re_SNPbps[x])
  if (split) bp_cliques <- split_cliques(bp_cliques, clstgap)
  
  while (any(is.na(binvector)) && length(bp_cliques) > 0) {
    scores <- sapply(bp_cliques, clique_score)
    best <- bp_cliques[[which.max(scores)]]
    if (codechange) {
      idx <- match(best, re_SNPbps)
      corr_adj <- code_sign_opt(idx, rV2_mat)$matrix
      # optionally reselect clique here from corrected matrix
    }
    binvector[match(best, SNPbps)] <- bin_id
    bin_id <- bin_id + 1
    bp_cliques <- lapply(bp_cliques, function(x) setdiff(x, best))
    bp_cliques <- Filter(function(x) length(x) > 1, bp_cliques)
  }
  
  if (all(is.na(binvector))) binvector <- seq_along(binvector)
  return(binvector)
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
