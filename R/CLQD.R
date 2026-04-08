# ─────────────────────────────────────────────────────────────────────────────
# CLQD.R  –  Clique-based LD partitioning  (C++ accelerated)
# ─────────────────────────────────────────────────────────────────────────────

#' Clique-Based LD Block Detection
#'
#' @description
#' Partitions SNPs in a genomic window into LD cliques. The LD matrix is
#' computed by the C++ Armadillo kernel; adjacency is built by
#' \code{build_adj_matrix_cpp}; cliques are found by igraph.
#'
#' @param subgeno Numeric matrix (individuals x SNPs), 0/1/2 coded.
#' @param subSNPinfo Data frame: col 1 = rsID, col 2 = bp position.
#' @param adj_subgeno n x p matrix from \code{prepare_geno()}.
#' @param CLQcut Numeric in (0,1]. LD threshold for edges. Default 0.5.
#' @param clstgap Integer. Max bp gap within clique (split=TRUE).
#' @param CLQmode "Density" (default) or "Maximal".
#' @param codechange Logical. Sign-flip heuristic. Default FALSE.
#' @param checkLargest Logical. Dense-core pre-pass for >= 500 SNPs.
#' @param split Logical. Split cliques over gaps. Default FALSE.
#' @param digits Integer. Rounding for LD. -1 = none (default).
#' @param n_threads Integer. OpenMP threads. Default 1.
#' @param verbose Logical.
#'
#' @return Integer vector length p: clique ID per SNP, NA = singleton.
#' @export
CLQD <- function(
    subgeno,
    subSNPinfo,
    adj_subgeno,
    CLQcut       = 0.5,
    clstgap      = 40000,
    CLQmode      = c("Density", "Maximal"),
    codechange   = FALSE,
    checkLargest = FALSE,
    split        = FALSE,
    digits       = -1L,
    n_threads    = 1L,
    verbose      = FALSE
) {
  CLQmode <- match.arg(CLQmode)
  SNPbps  <- as.numeric(subSNPinfo[[2L]])
  p       <- ncol(adj_subgeno)

  ld_mat  <- compute_r2_cpp(adj_subgeno,
                             digits    = as.integer(digits),
                             n_threads = as.integer(n_threads))
  adj_bin <- build_adj_matrix_cpp(ld_mat, CLQcut)

  binvector <- rep(NA_integer_, p)
  re_SNPbps <- SNPbps
  re_idx    <- seq_len(p)
  bin_id    <- 1L

  if (isTRUE(checkLargest) && p >= 500L) {
    repeat {
      g        <- igraph::graph_from_adjacency_matrix(adj_bin, mode = "undirected", diag = FALSE)
      comp     <- igraph::components(g)
      largest  <- which.max(comp$csize)
      comp_idx <- which(comp$membership == largest)
      sub_mat  <- ld_mat[comp_idx, comp_idx, drop = FALSE]
      cg       <- igraph::graph_from_adjacency_matrix(
        build_adj_matrix_cpp(sub_mat, CLQcut), mode = "undirected", diag = FALSE)
      core_vals <- igraph::coreness(cg)
      if (stats::median(core_vals) > 80 && max(core_vals) > 100) {
        deg      <- rowSums(adj_bin)
        high_deg <- which(deg >= stats::quantile(deg, 0.7, names = FALSE))
        densest  <- high_deg[which.max(vapply(high_deg, function(v) {
          nb <- which(adj_bin[v, ] > 0L)
          if (length(nb) < 2L) return(0)
          sum(adj_bin[nb, nb]) / (length(nb) * (length(nb) - 1L))
        }, numeric(1)))]
        rm_idx            <- which(adj_bin[densest, ] > 0L)
        binvector[re_idx[rm_idx]] <- bin_id
        bin_id    <- bin_id + 1L
        keep      <- setdiff(seq_along(re_SNPbps), rm_idx)
        re_SNPbps <- re_SNPbps[keep]
        re_idx    <- re_idx[keep]
        ld_mat    <- ld_mat[keep, keep, drop = FALSE]
        adj_bin   <- adj_bin[keep, keep, drop = FALSE]
        if (length(re_SNPbps) < 500L) break
      } else { break }
    }
  }

  if (isTRUE(verbose))
    message(sprintf("[CLQD] Building graph on %d SNPs...", length(re_SNPbps)))

  g       <- igraph::graph_from_adjacency_matrix(adj_bin, mode = "undirected", diag = FALSE)
  cliques <- igraph::max_cliques(g, min = 2L)

  if (isTRUE(verbose))
    message(sprintf("[CLQD] Found %d maximal cliques.", length(cliques)))

  clique_score <- function(bp_pos) {
    if (CLQmode == "Maximal") length(bp_pos)
    else                      length(bp_pos) / (diff(range(bp_pos)) / 1000 + 1)
  }

  cliques_bp <- lapply(cliques, function(x) re_SNPbps[x])
  if (isTRUE(split))
    cliques_bp <- .split_cliques_by_gap(cliques_bp, clstgap)

  while (any(is.na(binvector)) && length(cliques_bp) > 0L) {
    scores <- vapply(cliques_bp, clique_score, numeric(1L))
    best   <- cliques_bp[[which.max(scores)]]
    binvector[match(best, SNPbps)] <- bin_id
    bin_id     <- bin_id + 1L
    cliques_bp <- Filter(function(x) length(x) > 1L,
                         lapply(cliques_bp, function(x) setdiff(x, best)))
  }

  if (all(is.na(binvector))) binvector <- seq_along(binvector)
  binvector
}

.split_cliques_by_gap <- function(cliques_bp, gapdist) {
  cliques_bp <- lapply(cliques_bp, function(x) sort(as.numeric(stats::na.omit(x))))
  final <- list()
  while (length(cliques_bp) > 0L) {
    newlist <- list()
    for (cl in cliques_bp) {
      gaps <- diff(cl)
      if (any(gaps > gapdist)) {
        sp    <- which(gaps > gapdist)[1L]
        final <- c(final, list(cl[seq_len(sp)]))
        rest  <- cl[(sp + 1L):length(cl)]
        if (length(rest) > 1L) newlist <- c(newlist, list(rest))
      } else { final <- c(final, list(cl)) }
    }
    cliques_bp <- newlist
  }
  final
}
