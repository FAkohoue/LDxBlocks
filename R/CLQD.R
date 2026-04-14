# -----------------------------------------------------------------------------
# CLQD.R  -  Clique-based LD partitioning  (C++ accelerated)
# -----------------------------------------------------------------------------

#' Clique-Based LD Block Detection
#'
#' @description
#' Partitions SNPs in a genomic window into LD cliques. Three algorithms are
#' available via \code{CLQmode}:
#'
#' \describe{
#'   \item{\code{"Density"} (default)}{Standard Big-LD clique scoring -- edges
#'     built from the full r? matrix, cliques found by Bron-Kerbosch
#'     (\code{igraph::max_cliques}), scored by density = size / span.
#'     Exact, but exponential worst-case. Use \code{checkLargest = TRUE} to
#'     guard against blowup on dense WGS panels.}
#'   \item{\code{"Maximal"}}{Same algorithm as Density but scored by clique
#'     size only (original Big-LD maximality criterion).}
#'   \item{\code{"Louvain"}}{Louvain community detection
#'     (\code{igraph::cluster_louvain}). Polynomial time O(n log n) -- never
#'     exponential. Recommended for WGS panels (>500 SNPs per window) or
#'     any time \code{CLQmode = "Density"} is slow. Communities are returned
#'     as bin IDs matching the Density output format so all downstream code
#'     is unaffected.}
#'   \item{\code{"Leiden"}}{Leiden community detection
#'     (\code{igraph::cluster_leiden}). Stricter than Louvain -- produces
#'     well-connected communities without disconnected nodes. Requires
#'     igraph >= 1.3.0. Slightly slower than Louvain but higher quality
#'     communities. Recommended over Louvain when resolution matters.}
#' }
#'
#' When \code{max_bp_distance} is supplied (> 0), only SNP pairs within that
#' physical distance have their r? computed (\code{compute_r2_sparse_cpp}).
#' This reduces the O(p?) LD computation to near-O(p) for WGS panels where
#' long-range LD is negligible (typically > 500 kb). Applied before any
#' clique algorithm.
#'
#' @param subgeno Numeric matrix (individuals x SNPs), 0/1/2 coded.
#' @param subSNPinfo Data frame: col 1 = rsID, col 2 = bp position.
#' @param adj_subgeno n x p matrix from \code{prepare_geno()}.
#' @param CLQcut Numeric in (0,1]. LD threshold for edges. Default 0.5.
#' @param clstgap Integer. Max bp gap within clique (split=TRUE).
#' @param CLQmode Algorithm:
#'   \code{"Density"} (default, exact, exponential worst-case),
#'   \code{"Maximal"} (exact, maximality criterion),
#'   \code{"Louvain"} (polynomial, recommended for WGS),
#'   \code{"Leiden"} (polynomial, stricter communities, igraph >= 1.3.0).
#' @param codechange Logical. Sign-flip heuristic. Default FALSE.
#' @param checkLargest Logical. Dense-core pre-pass for >= 500 SNPs when
#'   using Density or Maximal mode. Default FALSE. Has no effect when
#'   CLQmode is Louvain or Leiden (already polynomial).
#' @param split Logical. Split cliques over gaps. Default FALSE.
#' @param digits Integer. Rounding for LD. -1 = none (default).
#' @param n_threads Integer. OpenMP threads. Default 1.
#' @param max_bp_distance Integer. Maximum base-pair distance between a
#'   SNP pair for its r? to be computed. Pairs beyond this distance are
#'   assumed to be in negligible LD and set to 0. Default \code{0L} =
#'   disabled (compute all pairs, original behaviour). Recommended value
#'   for WGS panels: \code{500000L} (500 kb). Has no effect when the
#'   window spans less than \code{max_bp_distance}.
#' @param verbose Logical.
#'
#' @return Integer vector length p: clique/community ID per SNP, NA = singleton.
#' @keywords internal
CLQD <- function(
    subgeno,
    subSNPinfo,
    adj_subgeno,
    CLQcut          = 0.5,
    clstgap         = 40000,
    CLQmode         = c("Density", "Maximal", "Louvain", "Leiden"),
    codechange      = FALSE,
    checkLargest    = FALSE,
    split           = FALSE,
    digits          = -1L,
    n_threads       = 1L,
    max_bp_distance = 0L,
    verbose         = FALSE
) {
  CLQmode <- match.arg(CLQmode)
  SNPbps  <- as.numeric(subSNPinfo[[2L]])
  p       <- ncol(adj_subgeno)

  # -- LD matrix computation --------------------------------------------------
  # When max_bp_distance > 0 and at least one pair exceeds that distance,
  # use sparse computation: only pairs within max_bp_distance bp are computed,
  # all others are set to 0 in the adjacency matrix. This replaces O(p?) with
  # near-O(p) for WGS panels where long-range LD is negligible.
  use_sparse <- (max_bp_distance > 0L) &&
    (diff(range(SNPbps)) > max_bp_distance)

  if (use_sparse) {
    # compute_r2_sparse_cpp returns triplet (row, col, r2) for pairs within
    # max_bp_distance that also pass `threshold = CLQcut` (avoids storing zeros).
    # bp must be an integer vector sorted in the same order as subgeno columns.
    bp_int   <- as.integer(SNPbps)
    sp_out   <- compute_r2_sparse_cpp(adj_subgeno,
                                      bp        = bp_int,
                                      max_bp_dist = as.integer(max_bp_distance),
                                      threshold   = CLQcut)
    # Rebuild binary adjacency from sparse triplet -- symmetric
    adj_bin  <- matrix(0L, nrow = p, ncol = p)
    if (length(sp_out$row) > 0L) {
      for (k in seq_along(sp_out$row)) {
        i <- sp_out$row[k]; j <- sp_out$col[k]
        adj_bin[i, j] <- 1L; adj_bin[j, i] <- 1L
      }
    }
    # ld_mat is not needed for Louvain/Leiden -- only build it for Density/Maximal
    ld_mat <- if (CLQmode %in% c("Density", "Maximal")) {
      m <- matrix(0, nrow = p, ncol = p)
      if (length(sp_out$row) > 0L) {
        for (k in seq_along(sp_out$row)) {
          m[sp_out$row[k], sp_out$col[k]] <- sp_out$r2[k]
          m[sp_out$col[k], sp_out$row[k]] <- sp_out$r2[k]
        }
      }
      diag(m) <- 1
      m
    } else NULL
  } else {
    ld_mat  <- compute_r2_cpp(adj_subgeno,
                              digits    = as.integer(digits),
                              n_threads = as.integer(n_threads))
    adj_bin <- build_adj_matrix_cpp(ld_mat, CLQcut)
  }

  # -- Route to algorithm -----------------------------------------------------
  if (CLQmode %in% c("Louvain", "Leiden")) {
    return(.clqd_community(adj_bin, SNPbps, CLQmode, CLQcut, clstgap, split, verbose))
  }

  # -- Density / Maximal path (original Bron-Kerbosch) -----------------------
  binvector <- rep(NA_integer_, p)
  re_SNPbps <- SNPbps
  re_idx    <- seq_len(p)
  bin_id    <- 1L

  # Dense-core pre-pass (checkLargest guard against exponential blowup)
  if (isTRUE(checkLargest) && p >= 500L) {
    repeat {
      g        <- igraph::graph_from_adjacency_matrix(adj_bin, mode = "undirected",
                                                      diag = FALSE)
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

  g       <- igraph::graph_from_adjacency_matrix(adj_bin, mode = "undirected",
                                                 diag = FALSE)
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

  # Only apply the all-singleton fallback when the graph had actual edges.
  # If the graph was empty (no pairs above CLQcut, or max_bp_distance pruned all
  # pairs), all-NA is the correct result -- every SNP is a true singleton.
  # The original fallback seq_along() was designed to prevent downstream crashes
  # when no cliques form despite edges existing (rare degenerate case).
  n_edges <- sum(adj_bin) / 2L
  if (all(is.na(binvector)) && n_edges > 0L)
    binvector <- seq_along(binvector)
  binvector
}


# -- Community-detection path (Louvain / Leiden) ----------------------------
# Polynomial time O(n log n) -- safe for WGS panels.
# Returns the same format as the Bron-Kerbosch path: integer vector of bin IDs.
.clqd_community <- function(adj_bin, SNPbps, CLQmode, CLQcut, clstgap, split,
                            verbose) {
  p <- nrow(adj_bin)

  g <- igraph::graph_from_adjacency_matrix(adj_bin, mode = "undirected",
                                           diag = FALSE, weighted = NULL)
  # Weight edges by physical proximity (inverse bp distance) so the community
  # detector respects local LD structure -- SNPs closer in bp space are
  # preferentially grouped.
  el    <- igraph::as_edgelist(g, names = FALSE)
  if (nrow(el) > 0L) {
    wt <- 1 / (abs(SNPbps[el[, 1L]] - SNPbps[el[, 2L]]) / 1000 + 1)
    igraph::E(g)$weight <- wt
  }

  clust <- if (CLQmode == "Leiden") {
    if (!igraph::is_igraph(g))
      stop("cluster_leiden requires igraph >= 1.3.0", call. = FALSE)
    igraph::cluster_leiden(g, objective_function = "modularity",
                           weights = igraph::E(g)$weight %||% NULL)
  } else {
    igraph::cluster_louvain(g, weights = igraph::E(g)$weight %||% NULL)
  }

  mem <- igraph::membership(clust)

  if (isTRUE(verbose))
    message(sprintf("[CLQD] %s found %d communities from %d SNPs.",
                    CLQmode, max(mem), p))

  # -- Louvain connectivity fix -------------------------------------------
  # Traag et al. (2019) showed that Louvain can produce internally
  # disconnected communities: two subgraphs grouped together despite having
  # no path between them within the community. For LD block detection this
  # is a correctness problem -- a block must be a contiguous genomic interval.
  # Leiden guarantees connected communities by design; Louvain does not.
  # Post-processing: for each Louvain community, extract the subgraph and
  # check connectivity. If disconnected, split into connected components and
  # assign each component a new unique community ID.
  if (CLQmode == "Louvain") {
    next_id <- max(mem) + 1L
    for (cid in unique(mem)) {
      nodes_c <- which(mem == cid)
      if (length(nodes_c) < 2L) next
      sg <- igraph::induced_subgraph(g, nodes_c)
      if (!igraph::is_connected(sg)) {
        comps <- igraph::components(sg)$membership
        # First component keeps the original community ID
        for (comp_id in unique(comps)) {
          comp_nodes <- nodes_c[comps == comp_id]
          if (comp_id == 1L) next  # keep original ID
          mem[comp_nodes] <- next_id
          next_id <- next_id + 1L
        }
      }
    }
    if (isTRUE(verbose)) {
      n_split <- max(mem) - length(unique(igraph::membership(clust)))
      if (n_split > 0L)
        message(sprintf("[CLQD] Louvain connectivity fix: split %d disconnected community/ies.",
                        n_split))
    }
  }

  # Communities with only 1 member -> NA (singleton, same as Density path)
  comm_sizes <- tabulate(mem)
  binvector  <- ifelse(comm_sizes[mem] > 1L, as.integer(mem), NA_integer_)

  # Optional: split communities at physical gaps > clstgap
  if (isTRUE(split) && clstgap > 0) {
    binvector <- .split_community_bins(binvector, SNPbps, clstgap)
  }

  binvector
}

# Re-split community bins at physical gaps (mirrors .split_cliques_by_gap logic)
.split_community_bins <- function(binvector, SNPbps, gapdist) {
  result    <- binvector
  next_id   <- max(binvector, na.rm = TRUE) + 1L
  for (bid in unique(stats::na.omit(binvector))) {
    idx    <- which(binvector == bid)
    bps    <- SNPbps[idx]
    ord    <- order(bps)
    idx    <- idx[ord]; bps <- bps[ord]
    gaps   <- diff(bps)
    splits <- which(gaps > gapdist)
    if (length(splits) == 0L) next
    # Reassign sub-segments after each gap
    segs   <- c(0L, splits, length(idx))
    for (s in seq_len(length(segs) - 1L)) {
      seg_idx <- idx[(segs[s] + 1L):segs[s + 1L]]
      if (s == 1L) {
        result[seg_idx] <- bid        # first segment keeps original ID
      } else {
        result[seg_idx] <- next_id
        next_id <- next_id + 1L
      }
    }
  }
  result
}

`%||%` <- function(a, b) if (!is.null(a)) a else b

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
