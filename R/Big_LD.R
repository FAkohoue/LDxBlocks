# -----------------------------------------------------------------------------
# Big_LD.R  -  Core per-chromosome LD block segmentation
#              Supports standard r^2 (default) and kinship-adjusted rV^2
#              C++ kernels used throughout for speed
# -----------------------------------------------------------------------------

#' LD Block Segmentation (r^2 or rV^2, C++ accelerated)
#'
#' @description
#' Core per-chromosome LD block detection. Two LD metrics are supported:
#'
#' \describe{
#'   \item{\code{method = "r2"} (default)}{Standard squared Pearson
#'     correlation, computed by the C++ Armadillo kernel. No kinship matrix
#'     is required. Suitable for large datasets (hundreds of thousands to
#'     millions of markers) and unstructured or mildly structured populations.}
#'   \item{\code{method = "rV2"}}{Kinship-adjusted squared correlation
#'     (Kim et al. 2018). Requires computing and inverting a GRM via
#'     AGHmatrix + ASRgenomics. Recommended for highly related populations
#'     (livestock, inbred lines) with moderate marker counts (< 200 k per chr).}
#' }
#'
#' For genome-wide analyses use \code{\link{run_Big_LD_all_chr}}.
#' For automatic parameter tuning use \code{\link{tune_LD_params}}.
#'
#' @param geno Numeric matrix (individuals x SNPs), 0/1/2. Row names used
#'   as individual IDs; auto-generated if absent.
#' @param SNPinfo Data frame: col 1 = rsID, col 2 = bp position.
#' @param method Character. \code{"r2"} (default) or \code{"rV2"}.
#' @param CLQcut Numeric in (0,1]. LD threshold for clique edges. Default 0.5.
#' @param clstgap Integer. Max bp gap within clique (split=TRUE). Default 40000.
#' @param leng Integer. Boundary-scan half-window (SNPs). Default 200.
#' @param subSegmSize Integer. Max SNPs per CLQD call. Default 1500.
#' @param MAFcut Numeric. Minor allele frequency minimum. Default 0.05.
#' @param appendrare Logical. Append rare SNPs after block detection. Default FALSE.
#' @param singleton_as_block Logical. If \code{TRUE}, every SNP that passes MAF
#'   filtering but has pairwise r^2 below \code{CLQcut} with all neighbours
#'   (i.e. is not assigned to any clique) is returned as a single-SNP block with
#'   \code{start == end} and \code{length_bp == 1}. Default \code{FALSE}.
#'   These blocks are excluded from haplotype analysis by the default
#'   \code{min_snps = 3} threshold in \code{extract_haplotypes()}, but are
#'   useful for auditing coverage and for single-SNP feature engineering.
#' @param checkLargest Logical. Dense-core pre-pass for large windows. Default FALSE.
#' @param CLQmode \code{"Density"} (default) or \code{"Maximal"}.
#' @param kin_method Character. GRM whitening: \code{"chol"} (default) or
#'   \code{"eigen"}. Ignored when \code{method = "r2"}.
#' @param split Logical. Split cliques at gaps > clstgap. Default FALSE.
#' @param digits Integer. LD rounding: \code{-1} (default, no rounding)
#'   or a positive integer.
#' @param n_threads Integer. OpenMP threads for C++ LD kernel. Default 1.
#' @param seed Integer or NULL. Set for reproducibility.
#' @param verbose Logical. Progress messages.
#'
#' @return data.frame with columns: start, end, start.rsID, end.rsID,
#'   start.bp, end.bp.
#'
#' @seealso \code{\link{run_Big_LD_all_chr}}, \code{\link{CLQD}},
#'   \code{\link{compute_r2}}, \code{\link{compute_rV2}},
#'   \code{\link{tune_LD_params}}
#'
#' @references
#' Kim S-A et al. (2018) GENETICS 209(3):855-868.\cr
#' VanRaden PM (2008) J. Dairy Sci. 91(11):4414-4423.
#'
#' @examples
#' \dontrun{
#' set.seed(1)
#' m <- 80; b1 <- 60; b2 <- 60
#' make_block <- function(m, size, p, flip = 0.03) {
#'   seed_snp <- rbinom(m, 2, p)
#'   M <- matrix(seed_snp, nrow = m, ncol = size)
#'   for (j in seq_len(size)) {
#'     idx <- sample.int(m, max(1L, floor(flip * m)))
#'     M[idx, j] <- pmin(2L, pmax(0L, M[idx, j] + sample(c(-1L, 1L), length(idx), TRUE)))
#'   }
#'   M
#' }
#' G <- cbind(make_block(m, b1, 0.30), make_block(m, b2, 0.60))
#' colnames(G) <- paste0("rs", seq_len(ncol(G)))
#' rownames(G) <- paste0("ind", seq_len(nrow(G)))
#' pos <- c(seq(1, by = 1000, length.out = b1), seq(5e6, by = 1000, length.out = b2))
#' SNPinfo <- data.frame(SNP = colnames(G), POS = pos)
#' blocks <- LDxBlocks:::Big_LD(G, SNPinfo, method = "r2", CLQcut = 0.6,
#'                  leng = 30, subSegmSize = 120, verbose = FALSE)
#' head(blocks)
#' }
#' @keywords internal
# ─────────────────────────────────────────────────────────────────────────────
# LD-informed overlap resolution (package-scope so tests can access via :::)
# ─────────────────────────────────────────────────────────────────────────────
#' @noRd
.resolve_overlap <- function(blocks, adj_mat, k_rep = 10L) {
  i <- 1L
  while (i < nrow(blocks)) {
    sA <- blocks[i,   1L]; eA <- blocks[i,   2L]
    sB <- blocks[i+1L,1L]; eB <- blocks[i+1L,2L]

    if (eA < sB) { i <- i + 1L; next }  # no overlap, advance

    # Three zones (use seq.int with explicit length to avoid R's
    # descending-sequence behaviour when start > end).
    left_core  <- if (sB > sA) seq.int(sA,   sB - 1L) else integer(0L)
    overlap    <- if (eA >= sB) seq.int(sB,  eA)       else integer(0L)
    right_core <- if (eB > eA) seq.int(eA + 1L, eB)   else integer(0L)

    # Fallback: right core empty means block B is fully inside block A
    # (eB <= eA), OR left core empty means blocks start at same position.
    # Either case: no valid LD anchor on one side -> union merge.
    has_left  <- length(left_core)  > 0L
    has_right <- length(right_core) > 0L
    if (!has_left || !has_right) {
      blocks[i,] <- c(min(sA, sB), max(eA, eB))
      blocks     <- blocks[-(i + 1L), , drop = FALSE]
      next  # recheck position i after merge (may now overlap i+1)
    }

    # Representatives: boundary-adjacent SNPs from each core.
    # Left:  last  k_rep of left_core  (closest to the overlap zone)
    # Right: first k_rep of right_core (closest to the overlap zone)
    k_L       <- min(k_rep, length(left_core))
    k_R       <- min(k_rep, length(right_core))
    left_reps  <- tail(left_core,  k_L)   # last  k_L indices of left_core
    right_reps <- head(right_core, k_R)   # first k_R indices of right_core

    # Compute signed LD score for each disputed SNP:
    #   score[k] = r2_left[k] - r2_right[k]
    # Positive = SNP prefers block A; negative = SNP prefers block B.
    # NA: scored as 0 (neutral -- constant/monomorphic column).
    # Ties are handled by the >= 0 boundary test below: a zero
    # contribution keeps cum_score non-negative if preceding SNPs
    # were positive, giving block A the benefit on a tie. A pure
    # all-zero sequence (no LD signal at all) gives cum_score = 0
    # everywhere, last_left = length(overlap), and the overlap goes
    # entirely to block A -- consistent with block A being detected
    # first in the sub-segment processing order.
    scores <- vapply(overlap, function(s) {
      r2v  <- as.numeric(col_r2_cpp(adj_mat, s))  # r2 of s with all cols
      r2_L <- mean(r2v[left_reps],  na.rm = TRUE)
      r2_R <- mean(r2v[right_reps], na.rm = TRUE)
      if (is.na(r2_L) || is.na(r2_R)) 0 else r2_L - r2_R
    }, numeric(1L))

    # Cumulative-score boundary rule (superior to last-left):
    #   cum_score[k] = sum(scores[1..k])
    #   Split at the LAST position where cum_score > 0.
    #   This is the point after which the remaining overlap SNPs are
    #   collectively in stronger LD with block B than with block A.
    #
    # Example: scores = [+0.3, -0.1, +0.5, -0.3, -0.5, -0.7]
    #   cum_score      = [+0.3, +0.2, +0.7, +0.4, -0.1, -0.8]
    #   last non-negative (>= 0) = position 4
    # All-negative: scores=[-0.3,-0.5,-0.2], cumsum=[-0.3,-0.8,-1.0]
    #   which(cumsum >= 0) = empty -> last_left=0 -> all-right path
    #   Block A ends at sB-1, block B keeps full overlap [sB..eB]
    # All-zero (ties): scores=[0,0,0], cumsum=[0,0,0]
    #   which(cumsum >= 0) = [1,2,3] -> last_left=3=length(overlap)
    #   -> all-left path: block A keeps full overlap, block B starts eA+1
    #   -> Block A gets overlap SNPs 1-4, block B gets 5-6.
    #
    # Compare last-left rule on same data:
    #   assignments    = [L, R, L, R, R, R]  (based on individual scores)
    #   last_left      = 3  (SNP 3 is last L)
    #   -> Block A gets 1-3, block B gets 4-6.
    #   SNP 4 (cum_score still +0.4, collectively A leads) is misassigned.
    #
    # The cumulative rule correctly uses all preceding evidence rather
    # than reacting to each SNP individually.
    cum_score <- cumsum(scores)
    last_left <- max(c(0L, which(cum_score >= 0)))

    if (last_left == 0L) {
      # Every disputed SNP belongs to block B.
      # Block A shrinks to [sA .. sB-1], block B stays [sB .. eB].
      blocks[i,   ] <- c(sA, sB - 1L)
      blocks[i+1L, ] <- c(sB, eB)
    } else if (last_left == length(overlap)) {
      # Every disputed SNP belongs to block A.
      # Block A stays [sA .. eA], block B shrinks to [eA+1 .. eB].
      # Both blocks survive -- this is NOT a union merge.
      blocks[i,   ] <- c(sA, eA)
      blocks[i+1L, ] <- c(eA + 1L, eB)
    } else {
      # Mixed assignment: split at last-left boundary.
      split_snp      <- overlap[last_left]
      blocks[i,   ]  <- c(sA, split_snp)
      blocks[i+1L, ] <- c(split_snp + 1L, eB)
    }

    i <- i + 1L  # both blocks are now non-overlapping; advance
  }
  blocks
}


Big_LD <- function(
    geno,
    SNPinfo,
    method       = c("r2", "rV2"),
    CLQcut       = 0.5,
    clstgap      = 40000,
    leng         = 200,
    subSegmSize  = 1500,
    MAFcut       = 0.05,
    appendrare   = FALSE,
    singleton_as_block = FALSE,
    checkLargest = FALSE,
    CLQmode         = c("Density", "Maximal", "Louvain", "Leiden"),
    kin_method      = "chol",
    split           = FALSE,
    max_bp_distance = 0L,
    digits       = -1L,
    n_threads    = 1L,
    seed         = NULL,
    verbose      = FALSE
) {
  if (!is.null(seed)) set.seed(seed)
  method <- match.arg(method)

  if (!is.matrix(geno)) geno <- as.matrix(geno)
  if (ncol(geno) != nrow(SNPinfo))
    stop("ncol(geno) [", ncol(geno), "] != nrow(SNPinfo) [", nrow(SNPinfo), "]")
  if (ncol(SNPinfo) != 2L)
    stop("SNPinfo must have exactly 2 columns: [rsID, bp_position]")

  ids <- rownames(geno)
  if (is.null(ids)) { ids <- sprintf("ind%04d", seq_len(nrow(geno))); rownames(geno) <- ids }

  # -- Prepare (centre or whiten) --------------------------------------------
  prep       <- prepare_geno(geno, method = method, kin_method = kin_method, verbose = verbose)
  adj_geno   <- prep$adj_geno
  V_inv_sqrt <- prep$V_inv_sqrt   # NULL for r^2

  Ogeno    <- geno
  OSNPinfo <- SNPinfo

  # -- MAF + monomorphic filter (C++) ---------------------------------------
  keep_poly <- maf_filter_cpp(geno, maf_cut = MAFcut)
  monoSNPs  <- OSNPinfo[!keep_poly, , drop = FALSE]
  geno      <- geno[, keep_poly, drop = FALSE]
  adjN      <- adj_geno[, keep_poly, drop = FALSE]
  SNPinfo   <- SNPinfo[keep_poly, , drop = FALSE]

  if (ncol(geno) < 2L) {
    warning("[Big_LD] Fewer than 2 polymorphic SNPs after MAF filter. Returning empty.")
    return(data.frame(start = integer(), end = integer(),
                      start.rsID = character(), end.rsID = character(),
                      start.bp = numeric(), end.bp = numeric()))
  }

  # ==========================================================================
  # Internal helpers
  # ==========================================================================

  # -- cutsequence.modi: find weak-LD boundaries using C++ boundary_scan ----
  cutsequence.modi <- function(geno, adj_geno, leng, subSegmSize, CLQcut, digits, n_threads) {
    p <- ncol(geno)
    if (p <= subSegmSize) return(list(p, NULL))
    if (isTRUE(verbose)) message("[Big_LD] Subsegmenting via C++ boundary scan...")

    modeNum <- 1L; lastnum <- 0L; cutpoints <- NULL; i <- leng

    while (i <= (p - leng)) {
      if ((i - lastnum) > 5L * subSegmSize) { modeNum <- 2L; break }

      found_cut <- FALSE
      for (j in unique(c(1L, min(10L, leng), leng))) {
        l_start <- i - j + 1L; l_end <- i
        r_start <- i + 1L;     r_end <- i + j
        if (r_end > p) next

        # C++ cross-boundary r^2 check
        sub_cols  <- l_start:r_end
        sub_adj   <- adj_geno[, sub_cols, drop = FALSE]
        cross_r2  <- compute_r2_cpp(sub_adj,
                                    digits    = as.integer(digits),
                                    n_threads = as.integer(n_threads))
        nl <- j; nr <- j
        cross_block <- cross_r2[seq_len(nl), (nl + 1L):(nl + nr), drop = FALSE]
        cross_block[cross_block < CLQcut] <- 0

        if (sum(cross_block, na.rm = TRUE) > 0) { i <- i + 1L; break }
        if (j == leng) found_cut <- TRUE
      }

      if (found_cut) {
        cutpoints <- c(cutpoints, i)
        lastnum   <- i
        if (isTRUE(verbose)) message(sprintf("[Big_LD]   Cut at SNP %d", i))
        i <- i + floor(leng / 2L)
      }
    }

    if (modeNum == 1L) {
      cutpoints <- c(0L, cutpoints, p)
      atfcut    <- NULL
      while (max(diff(cutpoints)) > subSegmSize) {
        diffseq    <- diff(cutpoints)
        recutpoint <- which(diffseq > subSegmSize)
        tt         <- cbind(cutpoints[recutpoint] + 1L, cutpoints[recutpoint + 1L])
        numvec     <- NULL
        tick       <- max(1L, as.integer(leng / 5L))

        for (k in seq_len(nrow(tt))) {
          st <- tt[k, 1L]; ed <- tt[k, 2L]
          if (ed > (p - leng)) ed <- p - leng
          cands <- (st + leng):(ed - leng)
          if (length(cands) == 0L) next

          # C++ boundary scan over candidate positions
          scan_cols <- st:(ed + leng)
          scan_cols <- scan_cols[scan_cols >= 1L & scan_cols <= p]
          scan_adj  <- adj_geno[, scan_cols, drop = FALSE]
          weak_flags <- boundary_scan_cpp(scan_adj,
                                          start     = tick + 1L,
                                          end       = ncol(scan_adj) - tick,
                                          half_w    = tick,
                                          threshold = CLQcut)
          # find positions with fewest cross-LD (weak_flags == 1)
          weakcount    <- as.integer(!weak_flags)   # 0 = weak, good for cut
          weakcount.s  <- sort(weakcount)
          weaks        <- weakcount.s[min(10L, length(weakcount.s))]
          rel_idx      <- which(weakcount <= weaks)
          abs_idx      <- rel_idx + st + leng - 1L
          abs_idx      <- abs_idx[abs_idx >= st & abs_idx <= ed]
          if (length(abs_idx) == 0L) abs_idx <- as.integer(round((st + ed) / 2L))
          nearcenter   <- vapply(abs_idx, function(x) abs((ed - x) - (x - st)), numeric(1L))
          addcut       <- abs_idx[which.min(nearcenter)]
          numvec       <- c(numvec, addcut)
          atfcut       <- c(atfcut, addcut)
        }
        cutpoints <- sort(c(cutpoints, numvec))
        if (!any(diff(cutpoints) > subSegmSize)) break
      }
    } else {
      cutpoints <- seq(subSegmSize, p, subSegmSize / 2L)
      atfcut    <- cutpoints
      if (max(cutpoints) == p) atfcut <- atfcut[-length(atfcut)]
      else                     cutpoints <- c(cutpoints, p)
    }
    list(cutpoints, atfcut)
  }

  # -- intervalCliqueList ----------------------------------------------------
  intervalCliqueList <- function(clstlist, allsnps, onlybp) {
    bp.clstlist <- lapply(clstlist, function(x) onlybp[x])
    bp.allsnps  <- lapply(allsnps,  function(x) onlybp[x])
    IMsize      <- length(bp.clstlist)
    adjacencyM  <- matrix(0L, IMsize, IMsize)
    for (i in seq_len(IMsize))
      for (j in seq_len(IMsize))
        adjacencyM[i, j] <- length(intersect(bp.allsnps[[i]], bp.allsnps[[j]]))
    diag(adjacencyM) <- 0L
    ig <- igraph::graph_from_adjacency_matrix(adjacencyM, mode = "undirected",
                                              weighted = TRUE, diag = FALSE)
    icliques <- if (max(igraph::coreness(ig)) > 10L) igraph::max_cliques(ig, min = 1L)
    else                                  igraph::cliques(ig, min = 1L)
    icliques  <- icliques[order(vapply(icliques, min, numeric(1L)))]
    intervals  <- lapply(icliques, function(x) sort(unique(unlist(bp.clstlist[x]))))
    weight.itv <- vapply(intervals, length, integer(1L))
    ir         <- t(vapply(intervals, range, numeric(2L)))
    ur         <- unique(ir)
    ri         <- cbind(ir, weight.itv)
    info_idx   <- apply(ur, 1L, function(x) {
      same <- intersect(which(ri[,1L] == x[1L]), which(ri[,2L] == x[2L]))
      same[which.max(ri[same, 3L])]
    })
    list(intervals[info_idx], ri[info_idx, 3L])
  }

  # -- find.maximum.indept ---------------------------------------------------
  find.maximum.indept <- function(sample.itv, sample.weight) {
    n              <- length(sample.itv)
    interval.range <- t(vapply(sample.itv, range, numeric(2L)))
    pre.range      <- vector("list", n)
    for (i in seq_len(n)) {
      idx <- which(interval.range[, 2L] < interval.range[i, 1L])
      pre.range[[i]] <- if (length(idx) > 0L) idx else NA_integer_
    }
    sources <- which(vapply(pre.range, function(x) all(is.na(x)), logical(1L)))
    if (length(sources) < n) {
      not.s <- setdiff(seq_len(n), sources)
      for (i in not.s) {
        pp <- sort(unique(unlist(pre.range[pre.range[[i]]])))
        pre.range[[i]] <- setdiff(pre.range[[i]], pp)
      }
      rw  <- rep(0, n); rw[sources] <- sample.weight[sources]
      ptr <- rep(0L, n); ptr[sources] <- NA_integer_
      info <- cbind(seq_len(n), rw, ptr, expl = replace(rep(0L,n), sources, 1L))
      for (i in not.s) {
        mp   <- pre.range[[i]]
        ni   <- info[mp, , drop = FALSE]
        mi   <- ni[ni[,2L] == max(ni[,2L]), , drop = FALSE][1L,]
        info[i,2L] <- sample.weight[i] + mi[2L]
        info[i,3L] <- mi[1L]; info[i,4L] <- 1L
      }
      si  <- which(info[,2L] == max(info[,2L]))[1L]
      pre <- info[si,3L]; is <- c(pre, si)
      while (!is.na(pre)) { pre <- info[pre,3L]; is <- c(pre, is) }
      list(indept.set = is[!is.na(is)], indept.set.weight = max(info[,2L]))
    } else {
      list(indept.set = which(sample.weight == max(sample.weight)),
           indept.set.weight = max(sample.weight))
    }
  }

  # -- constructLDblock ------------------------------------------------------
  constructLDblock <- function(clstlist, subSNPinfo) {
    safe_range <- function(v) {
      y <- as.numeric(stats::na.omit(v))
      if (length(y) < 2L) return(NULL)
      st <- min(y); ed <- max(y)
      if (!is.finite(st) || !is.finite(ed) || st > ed) return(NULL)
      c(start = st, end = ed)
    }
    SP   <- as.numeric(as.character(subSNPinfo[[2L]]))
    clst <- Filter(function(x) length(x) >= 2L,
                   lapply(clstlist, function(x) as.numeric(stats::na.omit(unlist(x)))))
    if (length(clst) == 0L) return(NULL)
    Total <- NULL
    while (length(clst) > 0L) {
      cr   <- Filter(Negate(is.null), lapply(clst, safe_range))
      asnp <- lapply(cr, function(rg) seq.int(rg["start"], rg["end"]))
      if (length(asnp) == 0L) break
      cand <- intervalCliqueList(clst, asnp, SP)
      ivs  <- cand[[1L]]; wts <- cand[[2L]]
      if (length(ivs) == 0L || length(wts) == 0L) break
      MWIS <- find.maximum.indept(ivs, wts)
      iset <- ivs[MWIS[[1L]]]
      sblk <- Filter(Negate(is.null), lapply(iset, function(x) {
        safe_range(as.numeric(stats::na.omit(match(x, SP))))
      }))
      if (length(sblk) == 0L) break
      sb   <- do.call(rbind, sblk)
      Total <- rbind(Total, sb)
      used  <- unique(unlist(lapply(sblk, function(rg) seq.int(rg["start"], rg["end"]))))
      clst  <- Filter(function(x) length(x) >= 2L, lapply(clst, function(x) setdiff(x, used)))
    }
    Total
  }

  # -- subBigLD --------------------------------------------------------------
  subBigLD <- function(sg, si, ag, CLQcut, clstgap, CLQmode, checkLargest, split, digits, n_threads, max_bp_distance = 0L) {
    bv   <- CLQD(sg, si, ag, CLQcut = CLQcut, clstgap = clstgap, CLQmode = CLQmode,
                 codechange = FALSE, checkLargest = checkLargest, split = split,
                 digits = digits, n_threads = n_threads, max_bp_distance = max_bp_distance, verbose = FALSE)
    if (all(is.na(bv))) return(matrix(integer(0), nrow=0L, ncol=2L))
    bins <- seq_len(max(bv[!is.na(bv)]))
    cl   <- lapply(bins, function(x) sort(which(bv == x)))
    cl   <- Filter(length, cl)  # drop empty bins (Louvain IDs may not be 1..n)
    if (length(cl)) cl <- cl[order(vapply(cl, min, numeric(1L)))]
    nowLD <- constructLDblock(cl, si)
    if (is.null(nowLD)) return(matrix(integer(0), nrow=0L, ncol=2L))
    if (!is.matrix(nowLD)) nowLD <- matrix(nowLD, nrow=1L)
    nowLD[order(nowLD[,1L]), , drop = FALSE]
  }

  # -- appendSGTs ------------------------------------------------------------
  appendSGTs <- function(LDblocks, Ogeno, OSNPinfo, CLQcut, clstgap,
                         checkLargest, CLQmode, prep_full, split, digits, n_threads,
                         max_bp_distance = 0L) {
    if (isTRUE(verbose)) message("[Big_LD] appendSGTs: assigning rare SNPs.")
    expandB <- NULL

    # Recompute adjusted geno for full original set
    Ogeno_c <- scale(Ogeno, center = TRUE, scale = FALSE)
    ag_full <- if (!is.null(prep_full$V_inv_sqrt)) {
      prep_full$V_inv_sqrt %*% Ogeno_c
    } else {
      Ogeno_c
    }

    snp1 <- which(OSNPinfo[,2L] < LDblocks[1L, 5L])
    if (length(snp1) > 2L) {
      OSNPs      <- seq_len(max(snp1))
      firstB     <- LDblocks[1L, ]
      secondSNPs <- which(OSNPinfo[,2L] >= firstB$start.bp & OSNPinfo[,2L] <= firstB$end.bp)
      blk_adj    <- ag_full[, c(secondSNPs, OSNPs), drop = FALSE]
      rv2        <- compute_r2_cpp(blk_adj, digits = as.integer(digits), n_threads = as.integer(n_threads))
      rv2        <- rv2[seq_along(secondSNPs), (length(secondSNPs)+1L):ncol(rv2), drop = FALSE]
      cor2ratio  <- colSums(rv2 > CLQcut) / length(secondSNPs)
      cor2numT   <- c(cor2ratio > 0.6, 1L)
      points2    <- min(which(cor2numT > 0))
      NsecondSNPs <- points2:max(secondSNPs)
      reOSNPs    <- setdiff(seq_len(max(NsecondSNPs)), NsecondSNPs)
      if (length(reOSNPs) > 1L) {
        sb <- subBigLD(Ogeno[,reOSNPs], OSNPinfo[reOSNPs,], ag_full[,reOSNPs,drop=FALSE],
                       CLQcut, clstgap, CLQmode, checkLargest, split, digits, n_threads)
        expandB <- rbind(expandB, sb + min(reOSNPs) - 1L)
      }
      firstSNPs <- NsecondSNPs
    } else {
      firstB    <- LDblocks[1L, ]
      firstSNPs <- which(OSNPinfo[,2L] >= firstB$start.bp & OSNPinfo[,2L] <= firstB$end.bp)
    }

    if (nrow(LDblocks) > 1L) {
      for (i in seq_len(nrow(LDblocks) - 1L)) {
        if (!length(firstSNPs)) next  # skip if first block has no SNPs in OSNPinfo
        secondB    <- LDblocks[i+1L, ]
        secondSNPs <- which(OSNPinfo[,2L] >= secondB$start.bp & OSNPinfo[,2L] <= secondB$end.bp)
        if (!length(secondSNPs)) {
          expandB <- rbind(expandB, range(firstSNPs))
          firstSNPs <- firstSNPs  # no next block found; keep and skip
          next
        }
        OSNPs      <- setdiff(max(firstSNPs):min(secondSNPs), c(max(firstSNPs), min(secondSNPs)))
        if (length(OSNPs) == 0L) {
          expandB <- rbind(expandB, range(firstSNPs)); firstSNPs <- secondSNPs
        } else {
          rv1 <- compute_r2_cpp(ag_full[,c(firstSNPs,OSNPs),drop=FALSE],
                                digits=as.integer(digits), n_threads=as.integer(n_threads))
          rv1 <- rv1[seq_along(firstSNPs),(length(firstSNPs)+1L):ncol(rv1),drop=FALSE]
          rv2 <- compute_r2_cpp(ag_full[,c(secondSNPs,OSNPs),drop=FALSE],
                                digits=as.integer(digits), n_threads=as.integer(n_threads))
          rv2 <- rv2[seq_along(secondSNPs),(length(secondSNPs)+1L):ncol(rv2),drop=FALSE]
          cor1ratio <- colSums(rv1 > CLQcut) / length(firstSNPs)
          cor2ratio <- colSums(rv2 > CLQcut) / length(secondSNPs)
          cor1numT  <- c(1L, cor1ratio > 0.6, 0L)
          cor2numT  <- c(0L, cor2ratio > 0.6, 1L)
          p1  <- min(max(firstSNPs) + max(which(cor1numT > 0)) - 1L, ncol(Ogeno))
          p2  <- min(max(firstSNPs) + max(which(cor2numT > 0)) - 1L, ncol(Ogeno))
          p2  <- max(p2, min(secondSNPs))  # ensure p2 <= max(secondSNPs) direction
          NF  <- min(firstSNPs):p1
          NS  <- min(p2, max(secondSNPs)):max(secondSNPs)
          if (max(NF) < min(NS)) {
            expandB <- rbind(expandB, range(NF))
            reO <- setdiff(min(NF):max(NS), c(NF, NS))
            if (length(reO) > 1L) {
              sb <- subBigLD(Ogeno[,reO], OSNPinfo[reO,], ag_full[,reO,drop=FALSE],
                             CLQcut, clstgap, CLQmode, checkLargest, split, digits, n_threads)
              if (nrow(sb) > 0L) expandB <- rbind(expandB, sb + min(reO) - 1L)
            }
            firstSNPs <- NS
          } else {
            rng <- min(firstSNPs):max(secondSNPs)
            sb  <- subBigLD(Ogeno[,rng], OSNPinfo[rng,], ag_full[,rng,drop=FALSE],
                            CLQcut, clstgap, CLQmode, checkLargest, split, digits, n_threads)
            if (nrow(sb) == 0L) { firstSNPs <- secondSNPs; next }
            sb  <- sb + min(rng) - 1L
            if (nrow(sb) == 1L) firstSNPs <- sb[1L,1L]:sb[1L,2L]
            else { expandB <- rbind(expandB, sb[-nrow(sb),]); firstSNPs <- sb[nrow(sb),1L]:sb[nrow(sb),2L] }
          }
        }
      }
    }

    if (length(firstSNPs) > 0L && max(firstSNPs) < (ncol(Ogeno) - 1L)) {
      OSNPs <- (max(firstSNPs)+1L):ncol(Ogeno)
      rv2   <- compute_r2_cpp(ag_full[,c(firstSNPs,OSNPs),drop=FALSE],
                              digits=as.integer(digits), n_threads=as.integer(n_threads))
      rv2   <- rv2[seq_along(firstSNPs),(length(firstSNPs)+1L):ncol(rv2),drop=FALSE]
      cor1ratio <- colSums(rv2 > CLQcut) / length(firstSNPs)
      cor1numT  <- c(1L, cor1ratio > 0.6, 0L)
      p1 <- max(firstSNPs) + max(which(cor1numT > 0)) - 1L
      NF <- min(firstSNPs):p1
      expandB <- rbind(expandB, range(NF))
      reO <- setdiff(min(NF):ncol(Ogeno), NF)
      if (length(reO) > 1L) {
        sb <- subBigLD(Ogeno[,reO], OSNPinfo[reO,], ag_full[,reO,drop=FALSE],
                       CLQcut, clstgap, CLQmode, checkLargest, split, digits, n_threads)
        expandB <- rbind(expandB, sb + min(reO) - 1L)
      }
    } else if (length(firstSNPs) > 0L) {
      expandB <- rbind(expandB, range(firstSNPs))
    }

    if (is.null(expandB) || nrow(expandB) == 0L) return(LDblocks)
    expandB <- expandB[expandB[,1L] != expandB[,2L], , drop = FALSE]
    data.frame(start      = expandB[,1L],
               end        = expandB[,2L],
               start.rsID = as.character(OSNPinfo[[1L]][expandB[,1L]]),
               end.rsID   = as.character(OSNPinfo[[1L]][expandB[,2L]]),
               start.bp   = as.numeric(OSNPinfo[[2L]][expandB[,1L]]),
               end.bp     = as.numeric(OSNPinfo[[2L]][expandB[,2L]]))
  }

  # ==========================================================================
  # Main loop
  # ==========================================================================

  cuts      <- cutsequence.modi(geno, adjN, leng, subSegmSize, CLQcut, digits, n_threads)
  cutpoints <- setdiff(cuts[[1L]], 0L)
  atfcut    <- cuts[[2L]]
  if (!is.null(atfcut)) atfcut <- sort(atfcut)

  cutblock <- cbind(c(1L, cutpoints + 1L), c(cutpoints, ncol(geno)))
  if (nrow(cutblock) > 1L) cutblock <- cutblock[-nrow(cutblock), , drop = FALSE]

  LDblocks <- matrix(NA_integer_, nrow(SNPinfo), 2L)

  for (i in seq_len(nrow(cutblock))) {
    nowst <- cutblock[i,1L]; nowed <- cutblock[i,2L]
    sg    <- geno[, nowst:nowed, drop = FALSE]
    ag    <- adjN[, nowst:nowed, drop = FALSE]
    si    <- SNPinfo[nowst:nowed, , drop = FALSE]

    bv   <- CLQD(sg, si, ag, CLQcut = CLQcut, clstgap = clstgap, CLQmode = CLQmode,
                 codechange = FALSE, checkLargest = checkLargest, split = split,
                 digits = digits, n_threads = n_threads, max_bp_distance = max_bp_distance, verbose = verbose)
    # Collect singleton indices (NA in bin vector = no clique partner found)
    if (isTRUE(singleton_as_block)) {
      sing_local <- which(is.na(bv))           # local indices within segment
      sing_global <- sing_local + (nowst - 1L) # map back to full-chromosome index
      if (!exists("singleton_idx")) singleton_idx <- integer(0L)
      singleton_idx <- c(singleton_idx, sing_global)
    }
    non_empty <- !all(is.na(bv))
    if (!non_empty) { nowLD <- matrix(integer(0), nrow=0L, ncol=2L) } else {
      bins <- seq_len(max(bv[!is.na(bv)]))
      cl   <- lapply(bins, function(x) sort(which(bv == x)))
      cl   <- Filter(length, cl)  # drop empty bins (Louvain IDs may not be 1..n)
      if (length(cl)) cl <- cl[order(vapply(cl, min, numeric(1L)))]
      nowLD <- constructLDblock(cl, si)
      if (!is.null(nowLD)) {
        if (!is.matrix(nowLD)) nowLD <- matrix(nowLD, nrow = 1L)
        nowLD <- nowLD + (nowst - 1L)
        nowLD <- nowLD[order(nowLD[,1L]), , drop = FALSE]
      } else { nowLD <- matrix(integer(0), nrow = 0L, ncol = 2L) }
    }  # end else(!non_empty)

    pre <- sum(!is.na(LDblocks[,1L]))
    if (nrow(nowLD) > 0L)
      LDblocks[(pre+1L):(pre+nrow(nowLD)), ] <- nowLD
    if (isTRUE(verbose)) message(sprintf("[Big_LD] Segment %d/%d done.", i, nrow(cutblock)))
  }

  done <- LDblocks[!is.na(LDblocks[,1L]), , drop = FALSE]

  # -- Optional re-merge across forced cut-points ----------------------------
  if (length(atfcut)) {
    newLDblocks <- matrix(NA_integer_, nrow(SNPinfo), 2L)
    consec <- 0L
    for (i in seq_len(nrow(done) - 1L)) {
      eb <- done[i,]; nb <- done[i+1L,]
      if (length(intersect(eb[2L]:nb[1L], atfcut)) > 0L) {
        consec <- consec + 1L
        if (consec > 1L) { addi <- max(which(!is.na(newLDblocks[,1L]))); eb <- newLDblocks[addi,] }
        rng <- range(c(eb, nb))
        bv  <- CLQD(geno[,rng[1L]:rng[2L]], SNPinfo[rng[1L]:rng[2L],],
                    adjN[,rng[1L]:rng[2L]],
                    CLQcut=CLQcut, clstgap=clstgap, CLQmode=CLQmode,
                    codechange=FALSE, checkLargest=checkLargest, split=split,
                    digits=digits, n_threads=n_threads, verbose=FALSE)
        bins <- seq_len(max(bv[!is.na(bv)]))
        cl   <- lapply(bins, function(x) sort(which(bv == x)))
        cl   <- Filter(length, cl)  # drop empty bins (Louvain IDs non-contiguous)
        if (length(cl)) cl <- cl[order(vapply(cl, min, numeric(1L)))]
        tmp  <- constructLDblock(cl, SNPinfo[rng[1L]:rng[2L],])
        tmp  <- tmp + (rng[1L]-1L); tmp <- tmp[order(tmp[,1L]),,drop=FALSE]
        ii   <- sum(!is.na(newLDblocks[,1L])) + 1L
        newLDblocks[ii:(ii+nrow(tmp)-1L), ] <- tmp
      } else {
        consec <- 0L
        ii <- min(which(is.na(newLDblocks[,1L]))); newLDblocks[ii,] <- eb
        if (i == nrow(done)-1L) { ii <- min(which(is.na(newLDblocks[,1L]))); newLDblocks[ii,] <- nb }
      }
    }
    LDblocks <- newLDblocks
  } else { LDblocks <- done }

  LDblocks <- LDblocks[!is.na(LDblocks[,1L]),,drop=FALSE]
  LDblocks <- LDblocks[order(LDblocks[,1L]),,drop=FALSE]

  # -- Resolve overlapping blocks (LD-informed split) ------------------------
  # Overlaps arise at sub-segment seams: a boundary SNP can be assigned to
  # a block in segment A and also to a block in segment B, producing two
  # blocks whose SNP index ranges overlap.
  #
  # OLD behaviour: blind union merge -- take [min(all), max(all)].
  #   Problem: merges blocks with low inter-block LD, creating an
  #   artificially large block with a deflated He and a high count of
  #   low-frequency haplotype alleles.
  #
  # NEW behaviour: LD-informed split --
  #   For each disputed SNP in the overlap zone, compute its mean r2
  #   with a sample of SNPs from the LEFT core (A-exclusive region) and
  #   from the RIGHT core (B-exclusive region). Assign it to whichever
  #   core it is in stronger LD with. The block boundary is set at the
  #   last left-assigned SNP, preserving contiguity of both blocks.
  #
  #   Fallback to union merge when either core is empty (one block fully
  #   contained in the other), where there is no LD anchor to compare.
  #
  # Parameters:
  #   k_rep = max representatives sampled from each core (bounds cost)
  #   Uses col_r2_cpp() -- O(n_ind) per call, negligible overhead.

  # -- LD-informed overlap resolution ----------------------------------------
  # Replaces the blind union merge [min(starts), max(ends)] with a per-SNP
  # LD-based assignment of disputed SNPs to one of the two overlapping blocks.
  #
  # GEOMETRY (after sorting LDblocks by start index):
  #   Block A: [sA ......... sB-1 | sB ......... eA]
  #   Block B:              [sB ......... eA | eA+1 ......... eB]
  #                          ^--- overlap zone ---^
  #   left_core:  [sA .. sB-1]   A-exclusive (always non-empty after sort)
  #   overlap:    [sB .. eA]     disputed SNPs

  LDblocks <- .resolve_overlap(LDblocks, adjN, k_rep = 10L)

  LDblocks <- .resolve_overlap(LDblocks, adjN, k_rep = 10L)

  # -- Build output data.frame -----------------------------------------------
  out <- data.frame(
    start      = LDblocks[,1L],
    end        = LDblocks[,2L],
    start.rsID = as.character(SNPinfo[[1L]][LDblocks[,1L]]),
    end.rsID   = as.character(SNPinfo[[1L]][LDblocks[,2L]]),
    start.bp   = as.numeric(SNPinfo[[2L]][LDblocks[,1L]]),
    end.bp     = as.numeric(SNPinfo[[2L]][LDblocks[,2L]])
  )

  # Re-index over full set including monomorphic SNPs
  all_bp  <- sort(c(as.numeric(monoSNPs[[2L]]), as.numeric(OSNPinfo[[2L]])))
  out$start <- match(out$start.bp, all_bp)
  out$end   <- match(out$end.bp,   all_bp)
  out <- out[order(out$start),]
  if (nrow(out) > 0L) {
    out[,c("start","end")]       <- t(apply(out[,c("start","end")],       1L, sort))
    out[,c("start.bp","end.bp")] <- t(apply(out[,c("start.bp","end.bp")], 1L, sort))
  }

  if (isTRUE(appendrare)) {
    prep_full <- list(V_inv_sqrt = V_inv_sqrt)
    out <- appendSGTs(out, Ogeno, OSNPinfo, CLQcut, clstgap, checkLargest,
                      CLQmode, prep_full, split, digits, n_threads)
  }

  # -- Optional: include singletons as single-SNP blocks ---------------------
  # Each SNP that passed MAF filtering but was not assigned to any clique
  # (no r^2 >= CLQcut partner within its sub-segment) becomes its own block.
  # start == end, length_bp == 1, start.rsID == end.rsID.
  # These are biologically meaningful: they mark low-LD regions, recombination
  # hotspots, or rapidly-evolving loci. They are excluded from haplotype
  # analysis by the default min_snps = 3 threshold in extract_haplotypes().
  if (isTRUE(singleton_as_block) && exists("singleton_idx") &&
      length(singleton_idx) > 0L) {
    # Map singleton indices to OSNPinfo positions (already over full SNP set)
    si_bp  <- as.numeric(SNPinfo[[2L]][singleton_idx])
    si_id  <- as.character(SNPinfo[[1L]][singleton_idx])
    # Re-index to full SNP set (including monomorphics) using bp position
    all_bp  <- sort(c(as.numeric(monoSNPs[[2L]]), as.numeric(OSNPinfo[[2L]])))
    si_full <- match(si_bp, all_bp)
    sing_df <- data.frame(
      start      = si_full,
      end        = si_full,
      start.rsID = si_id,
      end.rsID   = si_id,
      start.bp   = si_bp,
      end.bp     = si_bp,
      CHR        = if ("CHR" %in% names(out)) out$CHR[1L] else NA_character_,
      length_bp  = 1L,
      stringsAsFactors = FALSE
    )
    # Merge with multi-SNP blocks and re-sort by position
    out <- rbind(out, sing_df[, names(out), drop = FALSE])
    out <- out[order(out$start), ]
    if (isTRUE(verbose))
      message("[Big_LD] Added ", nrow(sing_df), " singleton SNP block(s).")
  }
  out
}
