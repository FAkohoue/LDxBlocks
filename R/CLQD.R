#' Clique-based LD block detection using \eqn{rV^2}
#'
#' @description Partition SNPs into LD cliques based on kinship-adjusted squared correlations (\eqn{rV^2}).
#'
#' @param subgeno Numeric genotype matrix (individuals x SNPs).
#' @param subSNPinfo Data frame with SNP metadata (col1 = rsID, col2 = base-pair positions).
#' @param adj_subgeno Kinship-adjusted centered genotype matrix (\eqn{V^{-1/2}} %*% centered G).
#' @param CLQcut Threshold for \eqn{rV^2} (0–1).
#' @param clstgap Maximum bp gap allowed within a clique.
#' @param CLQmode Clique priority mode: \code{"Density"} (default) or \code{"Maximal"}.
#' @param codechange Logical. Flip negative LD signs to maximize cliques.
#' @param checkLargest Logical. Use dense-core decomposition when SNPs > 500.
#' @param split Logical. Apply distance-based clique splitting.
#' @param digits Integer; rounding precision passed to \code{compute_rV2()}.
#' @param verbose Logical; if \code{TRUE}, print progress messages.
#'
#' @return Integer vector assigning each SNP to an LD block (NA = singleton).
#' @export
CLQD <- function(
    subgeno, subSNPinfo, adj_subgeno,
    CLQcut = 0.5,
    clstgap = 40000,
    CLQmode = c("Density", "Maximal"),
    codechange = FALSE,
    checkLargest = FALSE,
    split = FALSE,
    digits = 6,
    verbose = TRUE
) {
  CLQmode <- match.arg(CLQmode)
  SNPbps  <- as.numeric(subSNPinfo[[2]])

  # rV^2 -> adjacency
  rV2_mat <- compute_rV2(adj_subgeno, digits = digits)
  adj_bin <- ifelse(rV2_mat >= CLQcut, 1L, 0L)
  diag(adj_bin) <- 0L

  binvector  <- rep(NA_integer_, ncol(rV2_mat))
  re_SNPbps  <- SNPbps
  bin_id     <- 1L

  # optional dense decomposition
  if (isTRUE(checkLargest) && ncol(rV2_mat) >= 500L) {
    repeat {
      g <- igraph::graph_from_adjacency_matrix(adj_bin, mode = "undirected", diag = FALSE)
      comp <- igraph::components(g)
      largest <- which.max(comp$csize)
      comp_idx <- which(comp$membership == largest)

      sub_mat <- rV2_mat[comp_idx, comp_idx, drop = FALSE]
      cg <- igraph::graph_from_adjacency_matrix(sub_mat >= CLQcut, mode = "undirected", diag = FALSE)
      core_vals <- igraph::coreness(cg)

      if (stats::median(core_vals) > 80 && max(core_vals) > 100) {
        deg <- rowSums(adj_bin)
        high_deg <- which(deg >= stats::quantile(deg, 0.7, names = FALSE))
        densest <- high_deg[which.max(sapply(high_deg, function(v) {
          nb <- which(adj_bin[v, ] > 0L)
          if (length(nb) < 2) return(0)
          nb_mat <- adj_bin[nb, nb, drop = FALSE]
          sum(nb_mat) / (length(nb) * (length(nb) - 1))
        }))]

        idx_to_remove <- which(adj_bin[densest, ] > 0L)
        binvector[match(re_SNPbps[idx_to_remove], SNPbps)] <- bin_id
        bin_id <- bin_id + 1L

        keep_idx  <- setdiff(seq_along(re_SNPbps), idx_to_remove)
        re_SNPbps <- re_SNPbps[keep_idx]
        rV2_mat   <- rV2_mat[keep_idx, keep_idx, drop = FALSE]
        adj_bin   <- adj_bin[keep_idx, keep_idx, drop = FALSE]

        if (length(re_SNPbps) < 500L) break
      } else {
        break
      }
    }
  }
  if (verbose) message("[CLQD] End of pre-processing")

  g <- igraph::graph_from_adjacency_matrix(adj_bin, mode = "undirected", diag = FALSE)
  cliques <- igraph::max_cliques(g, min = 2)

  if (verbose) {
    message(sprintf("[CLQD] Detected cliques: %d", length(cliques)))
  }

  # helper: scoring
  clique_score <- function(bp_positions) {
    if (CLQmode == "Maximal") {
      length(bp_positions)
    } else {
      length(bp_positions) / (diff(range(bp_positions)) / 1000 + 1)  # avoid div by zero
    }
  }

  # optional sign correction heuristic inside clique (kept but not reselecting)
  code_sign_opt <- function(vt, corr_mat) {
    rin <- corr_mat[vt, vt, drop = FALSE]
    nr <- length(vt)
    code <- rep(1L, nr)
    change <- TRUE
    iter <- 0L
    max_iter <- 2^nr
    while (change && iter < max_iter) {
      change <- FALSE
      neg_counts <- rowSums(rin < 0)
      idx <- which.max(neg_counts)
      if (neg_counts[idx] > (nr - 1) / 2) {
        rin[idx, ] <- -rin[idx, ]
        rin[, idx] <- -rin[, idx]
        code[idx] <- 1L - code[idx]
        change <- TRUE
      }
      iter <- iter + 1L
    }
    list(matrix = rin, code = code)
  }

  # split long cliques by genomic gaps
  split_cliques <- function(cliques_bp, gapdist) {
    cliques_bp <- lapply(cliques_bp, function(x) sort(as.numeric(stats::na.omit(x))))
    final <- list()
    while (length(cliques_bp) > 0) {
      newlist <- list()
      for (c in cliques_bp) {
        gaps <- diff(c)
        if (any(gaps > gapdist)) {
          split_at <- which(gaps > gapdist)[1]
          front <- c[1:split_at]
          rest  <- c[(split_at + 1):length(c)]
          final <- c(final, list(front))
          if (length(rest) > 1) newlist <- c(newlist, list(rest))
        } else {
          final <- c(final, list(c))
        }
      }
      cliques_bp <- newlist
    }
    final
  }

  bp_cliques <- lapply(cliques, function(x) re_SNPbps[x])
  if (isTRUE(split)) bp_cliques <- split_cliques(bp_cliques, clstgap)

  while (any(is.na(binvector)) && length(bp_cliques) > 0) {
    scores <- sapply(bp_cliques, clique_score)
    best   <- bp_cliques[[which.max(scores)]]

    if (isTRUE(codechange)) {
      idx <- match(best, re_SNPbps)
      tmp <- code_sign_opt(idx, rV2_mat)
      rV2_mat[idx, idx] <- tmp$matrix
    }

    binvector[match(best, SNPbps)] <- bin_id
    bin_id <- bin_id + 1L
    bp_cliques <- lapply(bp_cliques, function(x) setdiff(x, best))
    bp_cliques <- Filter(function(x) length(x) > 1, bp_cliques)
  }

  if (all(is.na(binvector))) binvector <- seq_along(binvector)
  binvector
}
