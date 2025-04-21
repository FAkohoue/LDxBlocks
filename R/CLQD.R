#' @title Clique-based LD block detection using rV²
#' @description Partition SNPs into LD cliques based on kinship-adjusted squared correlations (rV²).
#'
#' @param subgeno A numeric genotype matrix (individuals × SNPs).
#' @param subSNPinfo A data.frame with SNP metadata (first column = rsID, second = base-pair positions).
#' @param adj_subgeno Kinship-adjusted centered genotype matrix (V^(-1/2) * centered G).
#' @param CLQcut Threshold for rV² correlation (0–1).
#' @param clstgap Maximum bp gap allowed within a clique.
#' @param CLQmode Clique priority mode: "Density" (default) or "Maximal".
#' @param codechange Logical. Flip negative LD signs to maximize cliques.
#' @param checkLargest Logical. Use dense-core decomposition when SNPs > 500.
#' @param split Logical. Apply physical-distance-based clique splitting.
#' @param verbose Logical. Print process updates.
#'
#' @return A vector assigning each SNP to an LD block (NA = singleton).
#' @export
#' 
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
