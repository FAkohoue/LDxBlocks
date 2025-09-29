#' Estimation of LD Block Regions Using Kinship-Adjusted Correlation (\eqn{rV^2})
#'
#' @name Big_LD
#'
#' @description
#' Detects linkage disequilibrium (LD) blocks in SNP genotype data using a Big-LD style approach
#' with kinship-adjusted squared correlation (\eqn{rV^2}). Suitable for related or structured populations.
#'
#' @param geno Numeric matrix/data frame of additive genotypes (individuals x SNPs).
#' @param SNPinfo Data frame with 2 columns: col1 = SNP ID (rsID), col2 = base-pair position.
#' @param CLQcut Numeric threshold for \eqn{rV^2} (0–1) used to define LD cliques. Default: 0.5.
#' @param clstgap Integer. Maximum physical gap (bp) allowed within a clique.
#' @param leng Integer. Window size for weak-LD boundary checks during segmentation. Default: 200.
#' @param MAFcut Numeric. Minor allele frequency threshold (exclude below unless \code{appendrare=TRUE}). Default: 0.05.
#' @param subSegmSize Integer. Max SNPs in a sub-segment (segments above are split). Default: 1500.
#' @param appendrare Logical. If TRUE, appends rare SNPs after block construction. Default: FALSE.
#' @param checkLargest Logical. Fast clique detection in large regions via dense-core decomposition. Default: FALSE.
#' @param CLQmode Character. Clique prioritization: \code{"Density"} (default) or \code{"Maximal"}.
#' @param rV2method \code{"chol"} (default) or \code{"eigen"} for computing \eqn{V^{-1/2}}.
#' @param split Logical. If TRUE, allow splitting cliques over long distances. Default: FALSE.
#' @param digits Integer; decimals for rounding \eqn{rV^2} inside \code{compute_rV2()}. Default: \code{6}.
#' @param seed Integer or \code{NULL}. If non-\code{NULL}, \code{set.seed(seed)} is called. Default: \code{NULL}.
#' @param verbose Logical; if \code{TRUE}, print progress messages. Default: \code{FALSE}.
#'
#' @return \code{data.frame} of detected LD blocks with columns:
#' \describe{
#'   \item{start}{SNP index of block start (including monomorphic SNPs).}
#'   \item{end}{SNP index of block end.}
#'   \item{start.rsID, end.rsID}{SNP IDs at block boundaries.}
#'   \item{start.bp, end.bp}{Base-pair positions of block boundaries.}
#' }
#'
#' @seealso \code{\link{CLQD}}, \code{\link{compute_rV2}}, \code{\link{get_V_inv_sqrt}}, \code{\link{run_Big_LD_all_chr}}
#'
#' @examples
#' \donttest{
#' ## No test:
#' set.seed(1)
#'
#' # Small synthetic data with two strong LD blocks
#' m  <- 80          # individuals
#' b1 <- 60          # SNPs in block 1
#' b2 <- 60          # SNPs in block 2
#'
#' make_block <- function(m, size, p, flip = 0.03) {
#'   # Start from a common "founder" SNP and create near-duplicates
#'   seed <- rbinom(m, 2, p)
#'   M <- matrix(seed, nrow = m, ncol = size)
#'   for (j in seq_len(size)) {
#'     idx <- sample.int(m, max(1, floor(flip * m)))
#'     delta <- sample(c(-1, 1), length(idx), replace = TRUE)
#'     M[idx, j] <- pmin(2, pmax(0, M[idx, j] + delta))  # keep in 0/1/2
#'   }
#'   M
#' }
#'
#' G <- cbind(
#'   make_block(m, b1, p = 0.30, flip = 0.03),
#'   make_block(m, b2, p = 0.60, flip = 0.03)
#' )
#'
#' colnames(G) <- paste0("rs", seq_len(ncol(G)))
#' rownames(G) <- paste0("ind", seq_len(nrow(G)))  # required so G.tuneup has names
#'
#' # Positions: two blocks separated by a large genomic gap
#' pos <- c(
#'   seq(1,    by = 1000, length.out = b1),
#'   seq(5e6,  by = 1000, length.out = b2)
#' )
#' SNPinfo <- data.frame(SNP = colnames(G), POS = pos)
#'
#' # Run with small windows so the example is quick & stable
#' blocks <- Big_LD(G, SNPinfo, CLQcut = 0.6, leng = 30, subSegmSize = 120, verbose = FALSE)
#' head(blocks)
#' ## End(No test)


#' }
#'
#' @export
Big_LD <- function(
    geno, SNPinfo, CLQcut = 0.5, clstgap = 40000, leng = 200, subSegmSize = 1500, MAFcut = 0.05,
    appendrare = FALSE, checkLargest = FALSE, CLQmode = "Density", rV2method = "chol", split = FALSE,
    digits = 6, seed = NULL, verbose = FALSE
) {
  if (!is.null(seed)) set.seed(seed)

  # ───────────────────────────────
  # Computing kinship and whitening
  # ───────────────────────────────

  # Ensure individual IDs exist (ASRgenomics::G.tuneup requires dimnames)
  ids <- rownames(geno)
  if (is.null(ids)) ids <- sprintf("ind%04d", seq_len(nrow(geno)))
  rownames(geno) <- ids

  if (isTRUE(verbose)) cat("Computing kinship matrix using VanRaden...\n")
  kinship <- AGHmatrix::Gmatrix(geno, method = "VanRaden")

  # Ensure names are present on G
  dimnames(kinship) <- list(ids, ids)

  if (isTRUE(verbose)) cat("Computing V^(-1/2)...\n")
  kinship_tuned <- ASRgenomics::G.tuneup(kinship, bend = TRUE, rcn = TRUE)$Gb
  V_inv_sqrt    <- get_V_inv_sqrt(kinship_tuned, rV2method)

  # adjust genotypes
  geno_centered <- scale(geno, center = TRUE, scale = FALSE)
  adj_geno      <- V_inv_sqrt %*% geno_centered

  #--- helpers ------------------------------------------------------------
  cutsequence.modi <- function(geno, adj_geno, leng, subSegmSize, CLQcut, digits) {
    if (isTRUE(verbose)) message("[Big_LD] Splitting whole sequence into subsegments using rV2.")
    modeNum <- 1L; lastnum <- 0L
    if (ncol(geno) <= subSegmSize) return(list(ncol(geno), NULL))

    calterms <- c(1L, 10L, leng)
    cutpoints <- NULL
    i <- leng

    while (i <= (ncol(geno) - leng)) {
      if ((i - lastnum) > 5 * subSegmSize) { modeNum <- 2L; break }
      cutnow <- FALSE
      for (j in calterms) {
        left_cols <- (i - j + 1):i
        right_cols <- (i + 1):(i + j)
        combined_cols <- c(left_cols, right_cols)
        if (max(combined_cols) > ncol(adj_geno)) next

        block_adj <- adj_geno[, combined_cols, drop = FALSE]
        nowr2 <- compute_rV2(block_adj, digits = digits)
        nowr2 <- nowr2[1:length(left_cols), (length(left_cols) + 1):(length(combined_cols)), drop = FALSE]
        nowr2[nowr2 < CLQcut] <- 0

        if (sum(nowr2, na.rm = TRUE) > 0) { i <- i + 1L; cutnow <- FALSE; break }
        if (j == leng) cutnow <- TRUE
      }
      if (cutnow) {
        cutpoints <- c(cutpoints, i)
        lastnum <- i
        if (isTRUE(verbose)) message(sprintf("[Big_LD] add cutpoint %d", i))
        i <- i + floor(leng / 2)
      }
    }

    if (modeNum == 1L) {
      cutpoints <- c(0, cutpoints, ncol(geno))
      atfcut <- NULL
      while (max(diff(cutpoints)) > subSegmSize) {
        diffseq <- diff(cutpoints)
        recutpoint <- which(diffseq > subSegmSize)
        tt <- cbind((cutpoints[recutpoint] + 1), cutpoints[recutpoint + 1])
        numvec <- NULL

        for (k in seq_len(nrow(tt))) {
          st <- tt[k, 1]; ed <- tt[k, 2]
          if (ed > (ncol(geno) - leng)) ed <- ncol(geno) - leng
          weakcount <- sapply((st + leng):(ed - leng), function(x) {
            tick <- as.integer(leng / 5)
            cols <- c((x - tick + 1):(x), (x + 1):(x + tick))
            if (min(cols) < 1 || max(cols) > ncol(adj_geno)) return(Inf)
            rv2 <- compute_rV2(adj_geno[, cols, drop = FALSE], digits = digits)
            rv2 <- rv2[1:tick, (tick + 1):(2 * tick), drop = FALSE]
            diag(rv2) <- 0
            length(which(rv2 >= CLQcut))
          })
          weakcount.s <- sort(weakcount)
          weaks <- weakcount.s[min(10L, length(weakcount.s))]
          weakpoint <- which(weakcount <= weaks) + st + leng - 1L
          nearcenter <- sapply(weakpoint, function(x) abs((ed - x) - (x - st)))
          addcut <- weakpoint[which.min(nearcenter)]
          numvec <- c(numvec, addcut)
          atfcut <- c(atfcut, addcut)
        }
        cutpoints <- sort(c(cutpoints, numvec))
        if (!any(diff(cutpoints) > subSegmSize)) break
      }
    } else {
      cutpoints <- seq(subSegmSize, ncol(geno), subSegmSize / 2)
      atfcut <- cutpoints
      if (max(cutpoints) == ncol(geno)) atfcut <- atfcut[-length(atfcut)] else cutpoints <- c(cutpoints, ncol(geno))
    }
    list(cutpoints, atfcut)
  }

  intervalCliqueList <- function(clstlist, allsnps, onlybp) {
    bp.clstlist <- lapply(clstlist, function(x) onlybp[x])
    bp.allsnps  <- lapply(allsnps,  function(x) onlybp[x])
    IMsize <- length(bp.clstlist)
    adjacencyM <- matrix(0, IMsize, IMsize)
    for (i in seq_len(IMsize)) for (j in seq_len(IMsize))
      adjacencyM[i, j] <- length(intersect(bp.allsnps[[i]], bp.allsnps[[j]]))
    diag(adjacencyM) <- 0
    interval.graph <- igraph::graph_from_adjacency_matrix(adjacencyM, mode = "undirected", weighted = TRUE, diag = FALSE)
    if (max(igraph::coreness(interval.graph)) > 10) {
      interval.cliques <- igraph::max_cliques(interval.graph, min = 1)
    } else {
      interval.cliques <- igraph::cliques(interval.graph, min = 1)
    }
    interval.cliques <- interval.cliques[order(sapply(interval.cliques, min))]
    intervals <- lapply(interval.cliques, function(x) sort(unique(unlist(bp.clstlist[x]))))
    weight.itv <- sapply(intervals, length)
    intervals.range <- t(sapply(intervals, range))
    unique.range <- unique(intervals.range)
    rangeinfo <- cbind(intervals.range, weight.itv)
    interval.info <- apply(unique.range, 1, function(x) {
      same <- intersect(which(rangeinfo[, 1] == x[1]), which(rangeinfo[, 2] == x[2]))
      same[which.max(rangeinfo[same, 3])]
    })
    list(intervals[interval.info], rangeinfo[interval.info, 3])
  }

  find.maximum.indept <- function(sample.itv, sample.weight) {
    n <- length(sample.itv)
    interval.range <- t(sapply(sample.itv, range))
    pre.range <- vector("list", n)
    for (i in seq_len(n)) {
      nowstart <- interval.range[i, 1]
      idx <- which(interval.range[, 2] < nowstart)
      if (length(idx) > 0) pre.range[[i]] <- idx
    }
    sources <- which(sapply(pre.range, function(x) all(is.na(x))))
    if (length(sources) < n) {
      not.s <- setdiff(seq_len(n), sources)
      for (i in not.s) {
        pre.pre <- sort(unique(unlist(pre.range[pre.range[[i]]])))
        pre.range[[i]] <- setdiff(pre.range[[i]], pre.pre)
      }
      route.weights <- rep(0, n)
      route.weights[sources] <- sample.weight[sources]
      pointers <- rep(0, n); pointers[sources] <- NA
      info <- cbind(seq_len(n), route.weights, pointers, explored = c(replace(rep(0, n), sources, 1)))
      for (i in not.s) {
        maybe.pred <- pre.range[[i]]
        now.info <- info[maybe.pred, , drop = FALSE]
        max.info <- now.info[which(now.info[, 2] == max(now.info[, 2])), , drop = FALSE]
        if (nrow(max.info) > 1) max.info <- max.info[1, , drop = FALSE]
        info[i, 2] <- sample.weight[i] + max.info[2]
        info[i, 3] <- max.info[1]
        info[i, 4] <- 1
      }
      start.itv <- which(info[, 2] == max(info[, 2]))[1]
      predecessor <- info[start.itv, 3]
      indept.set <- c(predecessor, start.itv)
      while (!is.na(predecessor)) {
        predecessor <- info[predecessor, 3]
        indept.set <- c(predecessor, indept.set)
      }
      indept.set <- indept.set[!is.na(indept.set)]
      indept.weight <- max(info[, 2])
    } else {
      indept.set <- which(sample.weight == max(sample.weight))
      indept.weight <- max(sample.weight)
    }
    list(indept.set = indept.set, indept.set.weight = indept.weight)
  }

  constructLDblock <- function(clstlist, subSNPinfo) {
    safe_range <- function(v) {
      y <- as.numeric(stats::na.omit(v))
      if (length(y) < 2) return(NULL)
      st <- min(y); ed <- max(y)
      if (!is.finite(st) || !is.finite(ed) || st > ed) return(NULL)
      c(start = st, end = ed)
    }
    Totalblocks <- NULL
    SNP_positions <- as.numeric(as.character(subSNPinfo[[2]]))
    clstlist <- lapply(clstlist, function(x) as.numeric(stats::na.omit(unlist(x))))
    clstlist <- Filter(function(x) length(x) >= 2, clstlist)
    if (length(clstlist) == 0) return(NULL)

    while (length(clstlist) > 0) {
      clst_ranges <- lapply(clstlist, safe_range)
      clst_ranges <- Filter(Negate(is.null), clst_ranges)
      allsnps <- lapply(clst_ranges, function(rg) seq.int(rg["start"], rg["end"]))
      if (length(allsnps) == 0) break

      candi <- intervalCliqueList(clstlist, allsnps, SNP_positions)
      intervals  <- candi[[1]]
      weight.itv <- candi[[2]]
      if (length(intervals) == 0 || length(weight.itv) == 0) break

      MWIS <- find.maximum.indept(intervals, weight.itv)
      indept.set <- intervals[MWIS[[1]]]

      subLDblocks_list <- lapply(indept.set, function(x) {
        idx <- as.numeric(stats::na.omit(match(x, SNP_positions)))
        safe_range(idx)
      })
      subLDblocks_list <- Filter(Negate(is.null), subLDblocks_list)
      if (length(subLDblocks_list) == 0) break

      subLDblocks <- do.call(rbind, subLDblocks_list)
      Totalblocks <- rbind(Totalblocks, subLDblocks)

      used <- unique(unlist(lapply(subLDblocks_list, function(rg) seq.int(rg["start"], rg["end"]))))
      clstlist <- lapply(clstlist, function(x) setdiff(x, used))
      clstlist <- Filter(function(x) length(x) >= 2, clstlist)
    }
    Totalblocks
  }

  subBigLD <- function(subgeno, subSNPinfo, adj_subgeno, CLQcut, clstgap, CLQmode, checkLargest, split, digits, verbose) {
    subbinvec <- CLQD(subgeno, subSNPinfo, adj_subgeno,
                      CLQcut = CLQcut, clstgap = clstgap, CLQmode = CLQmode,
                      codechange = FALSE, checkLargest = checkLargest,
                      split = split, digits = digits, verbose = verbose)
    bins <- seq_len(max(subbinvec[!is.na(subbinvec)]))
    clstlist <- lapply(bins, function(x) sort(which(subbinvec == x)))
    clstlist <- clstlist[order(sapply(clstlist, min))]
    nowLDblocks <- constructLDblock(clstlist, subSNPinfo)
    nowLDblocks[order(nowLDblocks[, 1]), , drop = FALSE]
  }

  appendSGTs <- function(LDblocks, Ogeno, OSNPinfo, CLQcut, clstgap, checkLargest, CLQmode, V_inv_sqrt, split, digits, verbose) {
    expandB <- NULL
    if (isTRUE(verbose)) message("[Big_LD] appendSGTs using rV2")
    Ogeno_centered <- scale(Ogeno, center = TRUE, scale = FALSE)
    adj_geno <- V_inv_sqrt %*% Ogeno_centered

    snp1 <- which(OSNPinfo[, 2] < LDblocks[1, 5])
    if (length(snp1) > 2) {
      OSNPs <- 1:max(snp1)
      firstB <- LDblocks[1, ]
      secondSNPs <- which(OSNPinfo[, 2] >= firstB$start.bp & OSNPinfo[, 2] <= firstB$end.bp)

      block_adj <- adj_geno[, c(secondSNPs, OSNPs), drop = FALSE]
      rv2 <- compute_rV2(block_adj, digits = digits)
      rv2 <- rv2[1:length(secondSNPs), (length(secondSNPs) + 1):ncol(rv2), drop = FALSE]

      cor2num   <- colSums(rv2 > CLQcut)
      cor2ratio <- cor2num / length(secondSNPs)
      cor2numT  <- c(cor2ratio > 0.6, 1)
      points2   <- min(which(cor2numT > 0))
      NsecondSNPs <- points2:max(secondSNPs)
      reOSNPs     <- setdiff(1:max(NsecondSNPs), NsecondSNPs)

      if (length(reOSNPs) > 1) {
        subBlocks <- subBigLD(Ogeno[, reOSNPs], OSNPinfo[reOSNPs, ], adj_geno[, reOSNPs, drop = FALSE],
                              CLQcut, clstgap, CLQmode, checkLargest, split, digits, verbose)
        subBlocks <- subBlocks + min(reOSNPs) - 1
        expandB <- rbind(expandB, subBlocks)
      }
      firstSNPs <- NsecondSNPs
    } else {
      firstB <- LDblocks[1, ]
      firstSNPs <- which(OSNPinfo[, 2] >= firstB$start.bp & OSNPinfo[, 2] <= firstB$end.bp)
    }

    if (nrow(LDblocks) > 1) {
      for (i in 1:(nrow(LDblocks) - 1)) {
        secondB <- LDblocks[i + 1, ]
        secondSNPs <- which(OSNPinfo[, 2] >= secondB$start.bp & OSNPinfo[, 2] <= secondB$end.bp)
        OSNPs <- setdiff(max(firstSNPs):min(secondSNPs), c(max(firstSNPs), min(secondSNPs)))

        if (length(OSNPs) == 0) {
          expandB <- rbind(expandB, range(firstSNPs))
          firstSNPs <- secondSNPs
        } else {
          rv1 <- compute_rV2(adj_geno[, c(firstSNPs, OSNPs), drop = FALSE], digits = digits)
          rv1 <- rv1[1:length(firstSNPs), (length(firstSNPs) + 1):ncol(rv1), drop = FALSE]
          cor1num   <- colSums(rv1 > CLQcut)
          cor1ratio <- cor1num / length(firstSNPs)

          rv2 <- compute_rV2(adj_geno[, c(secondSNPs, OSNPs), drop = FALSE], digits = digits)
          rv2 <- rv2[1:length(secondSNPs), (length(secondSNPs) + 1):ncol(rv2), drop = FALSE]
          cor2num   <- colSums(rv2 > CLQcut)
          cor2ratio <- cor2num / length(secondSNPs)

          cor1numT <- c(1, cor1ratio > 0.6, 0)
          cor2numT <- c(0, cor2ratio > 0.6, 1)

          points1 <- max(firstSNPs) + max(which(cor1numT > 0)) - 1
          NfirstSNPs <- min(firstSNPs):points1

          points2 <- max(firstSNPs) + max(which(cor2numT > 0)) - 1
          NsecondSNPs <- points2:max(secondSNPs)

          if (max(NfirstSNPs) < min(NsecondSNPs)) {
            expandB <- rbind(expandB, range(NfirstSNPs))
            reOSNPs <- setdiff(min(NfirstSNPs):max(NsecondSNPs), c(NfirstSNPs, NsecondSNPs))
            if (length(reOSNPs) > 1) {
              subBlocks <- subBigLD(Ogeno[, reOSNPs], OSNPinfo[reOSNPs, ], adj_geno[, reOSNPs, drop = FALSE],
                                    CLQcut, clstgap, CLQmode, checkLargest, split, digits, verbose)
              subBlocks <- subBlocks + min(reOSNPs) - 1
              expandB <- rbind(expandB, subBlocks)
            }
            firstSNPs <- NsecondSNPs
          } else {
            subBlocks <- subBigLD(Ogeno[, min(firstSNPs):max(secondSNPs)],
                                  OSNPinfo[min(firstSNPs):max(secondSNPs), ],
                                  adj_geno[, min(firstSNPs):max(secondSNPs), drop = FALSE],
                                  CLQcut, clstgap, CLQmode, checkLargest, split, digits, verbose)
            subBlocks <- subBlocks + min(firstSNPs) - 1
            if (nrow(subBlocks) == 1) {
              firstSNPs <- subBlocks[1, 1]:subBlocks[1, 2]
            } else {
              expandB <- rbind(expandB, subBlocks[-nrow(subBlocks), ])
              firstSNPs <- subBlocks[nrow(subBlocks), 1]:subBlocks[nrow(subBlocks), 2]
            }
          }
        }
      }
    }

    if (max(firstSNPs) < (ncol(Ogeno) - 1)) {
      OSNPs <- (max(firstSNPs) + 1):ncol(Ogeno)
      rv2 <- compute_rV2(adj_geno[, c(firstSNPs, OSNPs), drop = FALSE], digits = digits)
      rv2 <- rv2[1:length(firstSNPs), (length(firstSNPs) + 1):ncol(rv2), drop = FALSE]
      cor1num   <- colSums(rv2 > CLQcut)
      cor1ratio <- cor1num / length(firstSNPs)
      cor1numT  <- c(1, cor1ratio > 0.6, 0)
      points1   <- max(firstSNPs) + max(which(cor1numT > 0)) - 1
      NfirstSNPs <- min(firstSNPs):points1
      expandB <- rbind(expandB, range(NfirstSNPs))

      reOSNPs <- setdiff(min(NfirstSNPs):ncol(Ogeno), NfirstSNPs)
      if (length(reOSNPs) > 1) {
        subBlocks <- subBigLD(Ogeno[, reOSNPs], OSNPinfo[reOSNPs, ], adj_geno[, reOSNPs, drop = FALSE],
                              CLQcut, clstgap, CLQmode, checkLargest, split, digits, verbose)
        subBlocks <- subBlocks + min(reOSNPs) - 1
        expandB <- rbind(expandB, subBlocks)
      }
    } else {
      expandB <- rbind(expandB, range(firstSNPs))
    }

    expandB <- expandB[expandB[, 1] != expandB[, 2], , drop = FALSE]
    start.bp <- OSNPinfo[, 2][expandB[, 1]]
    end.bp   <- OSNPinfo[, 2][expandB[, 2]]
    start.rsID <- as.character(OSNPinfo[, 1][expandB[, 1]])
    end.rsID   <- as.character(OSNPinfo[, 1][expandB[, 2]])
    TexpandB <- data.frame(expandB, start.rsID, end.rsID, start.bp, end.bp)
    colnames(TexpandB) <- c("start", "end", "start.rsID", "end.rsID", "start.bp", "end.bp")
    TexpandB
  }
  #--- /helpers -----------------------------------------------------------

  # input checks
  Ogeno <- geno; OSNPinfo <- SNPinfo
  if (ncol(Ogeno) != nrow(OSNPinfo)) stop("Mismatch: ncol(geno) != nrow(SNPinfo)")
  if (ncol(OSNPinfo) != 2) stop("SNPinfo must have exactly 2 columns.")

  # drop monomorphic, MAF prune
  keep_poly <- apply(Ogeno, 2, function(x) { y <- x[!is.na(x)]; length(unique(y)) != 1 })
  monoSNPs  <- OSNPinfo[!keep_poly, , drop = FALSE]
  Ogeno     <- Ogeno[, keep_poly, drop = FALSE]
  OSNPinfo  <- OSNPinfo[keep_poly, , drop = FALSE]

  maf <- apply(Ogeno, 2, function(x) mean(x, na.rm = TRUE) / 2)
  maf <- ifelse(maf >= 0.5, 1 - maf, maf)
  mafprun <- which(maf >= MAFcut)

  geno  <- Ogeno[, mafprun, drop = FALSE]
  adjN  <- adj_geno[, mafprun, drop = FALSE]
  SNPinfo <- OSNPinfo[mafprun, , drop = FALSE]

  # subsegmenting
  cuts <- cutsequence.modi(geno, adjN, leng, subSegmSize, CLQcut = CLQcut, digits = digits)
  cutpoints <- setdiff(cuts[[1]], 0)
  atfcut <- cuts[[2]]
  if (!is.null(atfcut)) atfcut <- sort(atfcut)

  cutblock <- cbind(c(1, cutpoints + 1), c(cutpoints, ncol(geno)))
  if (nrow(cutblock) > 1) cutblock <- cutblock[-nrow(cutblock), , drop = FALSE]

  LDblocks <- matrix(NA_integer_, nrow(SNPinfo), 2)
  for (i in seq_len(nrow(cutblock))) {
    nowst <- cutblock[i, 1]; nowed <- cutblock[i, 2]
    subgeno     <- geno[, nowst:nowed, drop = FALSE]
    adj_subgeno <- adjN[,  nowst:nowed, drop = FALSE]
    subinfo     <- SNPinfo[nowst:nowed, , drop = FALSE]

    subbinvec <- CLQD(subgeno, subinfo, adj_subgeno,
                      CLQcut = CLQcut, clstgap = clstgap, CLQmode = CLQmode,
                      codechange = FALSE, checkLargest = checkLargest,
                      split = split, digits = digits, verbose = verbose)

    bins <- seq_len(max(subbinvec[!is.na(subbinvec)]))
    clstlist <- lapply(bins, function(x) sort(which(subbinvec == x)))
    clstlist <- clstlist[order(sapply(clstlist, min))]
    nowLDblocks <- constructLDblock(clstlist, subinfo)
    nowLDblocks <- nowLDblocks + (nowst - 1L)
    nowLDblocks <- nowLDblocks[order(nowLDblocks[, 1]), , drop = FALSE]

    pre <- sum(!is.na(LDblocks[, 1]))
    LDblocks[(pre + 1):(pre + nrow(nowLDblocks)), ] <- nowLDblocks
    if (isTRUE(verbose)) message(sprintf("[Big_LD] segment %d/%d", i, nrow(cutblock)))
  }

  done <- LDblocks[!is.na(LDblocks[, 1]), , drop = FALSE]

  # optional re-merge across atf cut points
  if (length(atfcut)) {
    newLDblocks <- matrix(NA_integer_, nrow(SNPinfo), 2)
    consecutive <- 0L
    for (i in 1:(nrow(done) - 1)) {
      endblock <- done[i, ]
      nextblock <- done[i + 1, ]
      gap <- c(endblock[2]:nextblock[1])
      if (length(intersect(gap, atfcut)) > 0) {
        consecutive <- consecutive + 1L
        if (consecutive > 1L) {
          addi <- max(which(!is.na(newLDblocks[, 1])))
          endblock <- newLDblocks[addi, ]
        }
        rng <- range(c(endblock, nextblock))
        subgeno     <- geno[, rng[1]:rng[2], drop = FALSE]
        adj_subgeno <- adjN[,  rng[1]:rng[2], drop = FALSE]
        subinfo     <- SNPinfo[rng[1]:rng[2], , drop = FALSE]
        subbinvec <- CLQD(subgeno, subinfo, adj_subgeno,
                          CLQcut = CLQcut, clstgap = clstgap, CLQmode = CLQmode,
                          codechange = FALSE, checkLargest = checkLargest,
                          split = split, digits = digits, verbose = verbose)
        bins <- seq_len(max(subbinvec[!is.na(subbinvec)]))
        clstlist <- lapply(bins, function(x) sort(which(subbinvec == x)))
        clstlist <- clstlist[order(sapply(clstlist, min))]
        tmp <- constructLDblock(clstlist, subinfo)
        tmp <- tmp + (rng[1] - 1L)
        tmp <- tmp[order(tmp[, 1]), , drop = FALSE]
        ii <- sum(!is.na(newLDblocks[, 1])) + 1L
        newLDblocks[ii:(ii + nrow(tmp) - 1), ] <- tmp
      } else {
        consecutive <- 0L
        ii <- min(which(is.na(newLDblocks[, 1])))
        newLDblocks[ii, ] <- endblock
        if (i == (nrow(done) - 1)) {
          ii <- min(which(is.na(newLDblocks[, 1])))
          newLDblocks[ii, ] <- nextblock
        }
      }
    }
    LDblocks <- newLDblocks
  } else {
    LDblocks <- done
  }

  LDblocks <- LDblocks[!is.na(LDblocks[, 1]), , drop = FALSE]
  LDblocks <- LDblocks[order(LDblocks[, 1]), , drop = FALSE]

  # merge overlaps
  i <- 1L
  while (i < nrow(LDblocks)) {
    nowb <- LDblocks[i, ]; nextb <- LDblocks[i + 1, ]
    if (nowb[2] < nextb[1]) { i <- i + 1L; next } else {
      newb <- c(min(c(nowb, nextb)), max(c(nowb, nextb)))
      LDblocks[i, ] <- newb
      LDblocks[i + 1, ] <- c(NA_integer_, NA_integer_)
      LDblocks <- LDblocks[!is.na(LDblocks[, 1]), , drop = FALSE]
    }
  }

  # map to bp and rsID
  start.bp  <- SNPinfo[, 2][LDblocks[, 1]]
  end.bp    <- SNPinfo[, 2][LDblocks[, 2]]
  start.rs  <- as.character(SNPinfo[, 1][LDblocks[, 1]])
  end.rs    <- as.character(SNPinfo[, 1][LDblocks[, 2]])
  out <- data.frame(LDblocks, start.rsID = start.rs, end.rsID = end.rs, start.bp = start.bp, end.bp = end.bp)
  colnames(out) <- c("start", "end", "start.rsID", "end.rsID", "start.bp", "end.bp")

  # place blocks on full set (incl. monomorphic) by bp
  mono_bp <- as.numeric(monoSNPs[[2]])
  poly_bp <- as.numeric(OSNPinfo[[2]])
  all_bp  <- sort(c(mono_bp, poly_bp))

  out$start.bp <- as.numeric(out$start.bp)
  out$end.bp   <- as.numeric(out$end.bp)
  out$start <- match(out$start.bp, all_bp)
  out$end   <- match(out$end.bp,   all_bp)
  out <- out[order(out$start), ]
  out[, c("start","end")] <- t(apply(out[, c("start","end")], 1, sort))
  out[, c("start.bp","end.bp")] <- t(apply(out[, c("start.bp","end.bp")], 1, sort))

  if (isTRUE(appendrare)) {
    out <- appendSGTs(out, Ogeno, OSNPinfo, CLQcut, clstgap, checkLargest, CLQmode, V_inv_sqrt, split, digits, verbose)
  }
  return(out)
}
