\name{Big_LD}
\alias{Big_LD}
\title{
Estimation of LD Block Regions Using Kinship-Adjusted Correlation (rV²)
}
\description{
Detects LD blocks in genotype data using kinship-adjusted squared correlation (\eqn{rV^2}) and a graph-theoretic clique method derived from Big-LD (Kim et al., 2017).
}
\usage{
Big_LD(geno, SNPinfo, CLQcut = 0.5, clstgap = 40000, leng = 200,
       subSegmSize = 1500, MAFcut = 0.05, appendrare = FALSE,
       checkLargest = FALSE, CLQmode = "Density", rV2method = "chol", split = FALSE)
}
\arguments{
  \item{geno}{Additive genotype matrix (individuals × SNPs).}
  \item{SNPinfo}{A data frame with 2 columns: rsID and base-pair position.}
  \item{CLQcut}{rV² threshold for defining LD cliques.}
  \item{clstgap}{Max base-pair distance allowed between SNPs within same clique.}
  \item{leng}{Window size (number of SNPs) for weak LD detection.}
  \item{subSegmSize}{Max SNPs per segment before forced subdivision.}
  \item{MAFcut}{Minor allele frequency cutoff. SNPs below this are removed unless \code{appendrare = TRUE}.}
  \item{appendrare}{Logical. Append rare SNPs post-block detection?}
  \item{checkLargest}{Enable dense-core decomposition in high-density LD regions.}
  \item{CLQmode}{Clique scoring method: "Density" or "Maximal".}
  \item{rV2method}{Method to compute V^{-1/2}: "chol" or "eigen".}
  \item{split}{Allow block splitting based on long physical gaps.}
}
\details{
This function integrates LD structure, physical positions, allele frequency, and kinship matrix into a unified pipeline.
It segments chromosomes into weak-LD intervals, identifies cliques, constructs interval graphs, and outputs maximal independent sets.
}
\value{
A \code{data.frame} with columns:
\itemize{
  \item \code{start}, \code{end}: Index of SNPs defining the block (relative to all SNPs including monomorphic ones).
  \item \code{start.rsID}, \code{end.rsID}: SNP IDs at the start and end.
  \item \code{start.bp}, \code{end.bp}: Base-pair positions of block boundaries.
}
}
\references{
Kim S, Scheet P, Kang EY. (2017). \emph{Bioinformatics}, 33(16):2526–2534. \doi{10.1093/bioinformatics/btx609}
}

\keyword{LD}
\keyword{genomics}

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
\code{\link{CLQD}}, \code{\link{compute_rV2}}, \code{\link{get_V_inv_sqrt}}, \code{\link{run_Big_LD_all_chr}}
}

\examples{
data(geno)
data(SNPinfo)
blocks <- Big_LD(geno, SNPinfo)

##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or standard data sets, see data().

## The function is currently defined as

Big_LD <- function(geno, SNPinfo, CLQcut = 0.5, clstgap = 40000, leng = 200, subSegmSize = 1500, MAFcut = 0.05, 
                   appendrare = FALSE, checkLargest = FALSE, CLQmode="Density", rV2method = "chol", split = FALSE) {
  # packages
  # library(igraph)
  #######################################################################################################
  # sub-Functions 1. cutsequence.modi, 2.intervalCliqueList, 3. find.maximum.indept, 4. constructLDblock, 5. CLQ
  
  # ───────────────────────────────
  # Compute kinship
  # ───────────────────────────────
  cat("Computing kinship matrix using VanRaden...\n")
  kinship <- AGHmatrix::Gmatrix(geno, method = "VanRaden")
  kinship_tuned <- ASRgenomics::G.tuneup(kinship, bend = TRUE, rcn = TRUE)$Gb
  
  cat("Computing V^(-1/2)...\n")
  V_inv_sqrt <- get_V_inv_sqrt(kinship_tuned, rV2method)
  
  
  # Genotype adjustment
  geno_centered <- scale(geno, center = TRUE, scale = FALSE)
  adj_geno <- V_inv_sqrt %*% geno_centered
  
  #-----------------------------------------
  # Execute cut sequence.modi
  #-----------------------------------------
  cutsequence.modi <- function(geno, adj_geno, leng, subSegmSize, CLQcut = 0.5) {
    print("Split whole sequence into subsegments using rV²")
    modeNum <- 1
    lastnum <- 0 
    
    if (ncol(geno) <= subSegmSize) {
      print("There is only one sub-region!")
      return(list(ncol(geno), NULL))
    } else {
      calterms <- c(1, 10, leng)
      cutpoints <- NULL
      i <- leng
      
      while (i <= (ncol(geno) - leng)) {
        if ((i - lastnum) > 5 * subSegmSize) {
          modeNum <- 2
          break
        }
        
        for (j in calterms) {
          left_cols <- (i - j + 1):i
          right_cols <- (i + 1):(i + j)
          combined_cols <- c(left_cols, right_cols)
          
          if (max(combined_cols) > ncol(adj_geno)) next
          
          block_adj <- adj_geno[, combined_cols, drop = FALSE]
          nowr2 <- compute_rV2(block_adj)
          nowr2 <- nowr2[1:length(left_cols), (length(left_cols) + 1):(length(combined_cols))]
          nowr2[nowr2 < CLQcut] <- 0
          
          if (sum(nowr2, na.rm = TRUE) > 0) {
            i <- i + 1
            cutnow <- FALSE
            break
          }
          
          if (j == leng) cutnow <- TRUE
        }
        
        if (cutnow) {
          cutpoints <- c(cutpoints, i)
          lastnum <- i
          print(i)
          i <- i + (leng / 2)
          cutnow <- FALSE
        }
      }
      
      if (modeNum == 1) {
        cutpoints <- c(0, cutpoints, ncol(geno))
        atfcut <- NULL
        
        while (max(diff(cutpoints)) > subSegmSize) {
          diffseq <- diff(cutpoints)
          recutpoint <- which(diffseq > subSegmSize)
          tt <- cbind((cutpoints[recutpoint] + 1), cutpoints[recutpoint + 1])
          numvec <- NULL
          
          for (i in 1:nrow(tt)) {
            st <- tt[i, 1]
            ed <- tt[i, 2]
            if (ed > (ncol(geno) - leng)) ed <- ncol(geno) - leng
            
            weakcount <- sapply((st + leng):(ed - leng), function(x) {
              tick <- as.integer(leng / 5)
              cols <- c((x - tick + 1):(x), (x + 1):(x + tick))
              if (min(cols) < 1 || max(cols) > ncol(adj_geno)) return(Inf)
              block_adj <- adj_geno[, cols, drop = FALSE]
              rv2 <- compute_rV2(block_adj)
              rv2 <- rv2[1:tick, (tick + 1):(2 * tick)]
              nowr2 <- rv2
              diag(nowr2) <- 0
              length(which(nowr2>= CLQcut))
            })
            
            weakcount.s <- sort(weakcount)
            weaks <- weakcount.s[10]
            weakpoint <- which(weakcount <= weaks)
            weakpoint <- weakpoint + st + leng - 1
            nearcenter <- sapply(weakpoint, function(x) abs((ed - x) - (x - st)), simplify = TRUE)
            addcut <- weakpoint[which.min(nearcenter)]
            print(paste("Add cutpoint", addcut))
            numvec <- c(numvec, addcut)
            atfcut <- c(atfcut, addcut)
          }
          
          cutpoints <- sort(c(cutpoints, numvec))
          newcandi <- which(diff(cutpoints) > subSegmSize)
          if (length(newcandi) == 0) break
        }
      }
      
      if (modeNum == 2) {
        cutpoints <- seq(subSegmSize, ncol(geno), subSegmSize / 2)
        atfcut <- cutpoints
        if (max(cutpoints) == ncol(geno)) {
          atfcut <- atfcut[-length(atfcut)]
        } else {
          cutpoints <- c(cutpoints, ncol(geno))
        }
      }
    }
    
    print("Cutting sequence, done")
    return(list(cutpoints, atfcut))
  }
  
  
  intervalCliqueList = function(clstlist, allsnps, onlybp) {
    bp.clstlist <- lapply(clstlist, function(x) onlybp[x])  ###
    bp.allsnps <- lapply(allsnps, function(x) onlybp[x])
    
    IMsize <- length(bp.clstlist)  ## adjacency matrix of intervals in interval graph
    adjacencyM <- matrix(0, IMsize, IMsize)
    for (i in 1:IMsize) {
      for (j in 1:IMsize) {
        adjacencyM[i, j] <- length(intersect(bp.allsnps[[i]], bp.allsnps[[j]]))
      }
    }
    diag(adjacencyM) <- 0
    interval.graph <- graph_from_adjacency_matrix(adjacencyM, mode = "undirected", weighted = TRUE, diag = FALSE, add.colnames = NULL)
    # print(paste("max coreness", max(coreness(interval.graph))))
    # print(paste("ecount", ecount(interval.graph), "vertex*5 ", 5*IMsize))
    if(max(coreness(interval.graph))>10){ #ecount(interval.graph)> 5*IMsize|
      interval.cliques <- max_cliques(interval.graph, min = 1)  
    }else{
      interval.cliques <- cliques(interval.graph, min = 1)  
    }
    interval.cliques <- interval.cliques[order(sapply(interval.cliques, min))]
    
    intervals <- lapply(interval.cliques, function(x) unlist(bp.clstlist[x]))
    intervals <- lapply(intervals, sort)
    intervals <- lapply(intervals, unique)
    weight.itv <- sapply(intervals, length)
    
    intervals.range <- t(sapply(intervals, range))
    unique.intervals.range <- unique(intervals.range)
    
    rangeinfo <- cbind(intervals.range, weight.itv)
    
    interval.info <- apply(unique.intervals.range, 1, function(x) {
      info1 = which(rangeinfo[, 1] == x[1])
      info2 = which(rangeinfo[, 2] == x[2])
      sameitv <- intersect(info1, info2)
      maxweight <- max(rangeinfo[sameitv, 3])
      sameitv[which(maxweight == rangeinfo[sameitv, 3])][1]
    })
    
    final.intervals <- intervals[interval.info]
    final.intervals.w <- rangeinfo[interval.info, 3]
    return(list(final.intervals, final.intervals.w))
  }
  # find maximum weight independent set input: clique interval list, clique weights
  find.maximum.indept = function(sample.itv, sample.weight) {
    n.of.sample <- length(sample.itv)
    interval.range <- t(sapply(sample.itv, range))
    pre.range <- as.list(rep(NA, n.of.sample))  #pre.range : range of predecessor
    ## pre.range : n by 2, i row (x,y) : possible predecessors of i interval are from x interval to y interval
    for (i in 1:n.of.sample) {
      nowstart <- interval.range[i, 1]
      if (length(which(interval.range[, 2] < nowstart)) > 0) {
        pre.range[[i]] <- which(interval.range[, 2] < nowstart)
      }
    }
    sources <- c(1:n.of.sample)[(sapply(pre.range, function(x) all(is.na(x)) == TRUE))]
    
    ## source of comparability graph of complement of Interval graph
    if (length(sources) < n.of.sample) {
      not.s <- setdiff(c(1:n.of.sample), sources)
      for (i in not.s) {
        pre.pre <- sort(unique(unlist(pre.range[pre.range[[i]]])))
        pre.range[[i]] <- setdiff(pre.range[[i]], pre.pre)
      }
      names(pre.range) <- sample.weight
      n.interval <- c(1:n.of.sample)
      route.weights <- rep(0, n.of.sample)  ##cumulative weights
      route.weights[sources] <- sample.weight[sources]
      pointers <- rep(0, n.of.sample)  ## predecessor of current interval
      pointers[sources] <- NA
      explored <- rep(0, n.of.sample)
      explored[sources] <- 1
      info <- cbind(n.interval, route.weights, pointers, explored)
      
      for (i in not.s) {
        maybe.pred <- pre.range[[i]]
        now.info <- info[maybe.pred, , drop = FALSE]
        max.info <- now.info[which(now.info[, 2] == max(now.info[, 2])), , drop = FALSE]
        if (dim(max.info)[1] > 1)
          max.info <- max.info[1, , drop = FALSE]
        info[i, 2] <- sample.weight[i] + max.info[2]
        info[i, 3] <- max.info[1]
        info[i, 4] <- 1
      }
      
      #### trace maximum independent set
      start.itv <- which(info[, 2] == max(info[, 2]))[1]
      predecessor <- info[start.itv, 3]
      indept.set <- c(predecessor, start.itv)
      while (!is.na(predecessor)) {
        predecessor <- info[predecessor, 3]
        indept.set <- c(predecessor, indept.set)
      }
      
      indept.set <- as.vector(indept.set)
      indept.set <- indept.set[-which(is.na(indept.set))]
      indept.set.weight <- max(info[, 2])
    } else {
      indept.set = which(sample.weight == max(sample.weight))
      indept.set.weight = max(sample.weight)
    }
    
    final.result <- list(indept.set = indept.set, indept.set.weight = indept.set.weight)
    return(final.result)
  }

  constructLDblock <- function(clstlist, subSNPinfo) {
    # ----------------------------------------------------------
    # Helper: only return a valid c(start,end) pair or NULL
    # ----------------------------------------------------------
    safe_range <- function(v) {
      y <- as.numeric(na.omit(v))
      if (length(y) < 2) return(NULL)
      st <- min(y, na.rm = TRUE)
      ed <- max(y, na.rm = TRUE)
      if (!is.finite(st) || !is.finite(ed) || st > ed) return(NULL)
      c(start = st, end = ed)
    }
    
    Totalblocks <- NULL
    SNP_positions <- as.numeric(as.character(subSNPinfo[[2]]))
    
    # Normalize and filter cliques to numeric vectors ≥ 1
    clstlist <- lapply(clstlist, function(x) {
      x <- as.numeric(na.omit(unlist(x)))
      x
    })
    clstlist <- Filter(function(x) length(x) >= 2, clstlist)
    if (length(clstlist) == 0) return(NULL)
    
    while (length(clstlist) > 0) {
      # 1) build only valid start–end pairs
      clst_ranges <- lapply(clstlist, safe_range)
      clst_ranges <- Filter(Negate(is.null), clst_ranges)
      
      # 2) expand each range into a full sequence of positions
      allsnps <- lapply(clst_ranges, function(rg) {
        seq.int(rg["start"], rg["end"])
      })
      if (length(allsnps) == 0) break
      
      # 3) find cliques of overlapping intervals
      candi.interval <- intervalCliqueList(clstlist, allsnps, SNP_positions)
      intervals  <- candi.interval[[1]]
      weight.itv <- candi.interval[[2]]
      if (length(intervals) == 0 || length(weight.itv) == 0) break
      
      # 4) pick maximum‐weight independent set
      MWIS      <- find.maximum.indept(intervals, weight.itv)
      indept.set <- intervals[MWIS[[1]]]
      
      # 5) turn each independent set into start–end blocks
      subLDblocks_list <- lapply(indept.set, function(x) {
        idx <- as.numeric(na.omit(match(x, SNP_positions)))
        safe_range(idx)
      })
      subLDblocks_list <- Filter(Negate(is.null), subLDblocks_list)
      if (length(subLDblocks_list) == 0) break
      
      subLDblocks <- do.call(rbind, subLDblocks_list)
      Totalblocks <- rbind(Totalblocks, subLDblocks)
      
      # 6) remove used SNPs from all cliques
      used <- unique(unlist(lapply(subLDblocks_list, function(rg) {
        seq.int(rg["start"], rg["end"])
      })))
      clstlist <- lapply(clstlist, function(x) setdiff(x, used))
      clstlist <- Filter(function(x) length(x) >= 2, clstlist)
    }
    
    return(Totalblocks)
  }
  
  subBigLD = function(subgeno, subSNPinfo,  CLQcut, clstgap, CLQmode, checkLargest){
    subbinvec <- CLQD(subgeno, subSNPinfo, adj_subgeno, CLQcut, clstgap, CLQmode, codechange = FALSE, checkLargest)
    # print('CLQ done!')
    bins <- c(1:max(subbinvec[which(!is.na(subbinvec))]))
    clstlist <- sapply(bins, function(x) which(subbinvec == x), simplify = FALSE)
    clstlist <- lapply(clstlist, sort)  ###
    clstlist <- clstlist[order(sapply(clstlist, min))]  ###
    nowLDblocks <- constructLDblock(clstlist, subSNPinfo)
    # print('constructLDblock done!')
    nowLDblocks <- nowLDblocks[order(nowLDblocks[, 1]), , drop = FALSE]
    return(nowLDblocks)
    
  }
  appendSGTs <- function(LDblocks, Ogeno, OSNPinfo, CLQcut, clstgap, checkLargest, CLQmode, V_inv_sqrt) {
    expandB <- NULL
    print("[INFO] Starting appendSGTs with rV² integration")
    
    # Step 1: Center and kinship-adjust genotype matrix
    print("[INFO] Centering and adjusting genotype matrix using V⁻¹ᐟ²")
    Ogeno_centered <- scale(Ogeno, center = TRUE, scale = FALSE)
    adj_geno <- V_inv_sqrt %*% Ogeno_centered
    
    snp1 <- which(OSNPinfo[, 2] < LDblocks[1, 5])
    
    if (length(snp1) > 2) {
      print("[INFO] Checking SNPs before the first LD block")
      OSNPs <- 1:max(snp1)
      firstB <- LDblocks[1, ]
      secondSNPs <- which(OSNPinfo[, 2] >= firstB$start.bp & OSNPinfo[, 2] <= firstB$end.bp)
      
      print("here")
      
      # Compute rV² between secondSNPs and OSNPs
      combined_idx <- c(secondSNPs, OSNPs)
      block_adj <- adj_geno[, combined_idx, drop = FALSE]
      rv2 <- compute_rV2(block_adj)
      rv2 <- rv2[1:length(secondSNPs), (length(secondSNPs) + 1):ncol(rv2)]
      
      # Threshold and count LD connections
      cor2num <- apply(rv2, 2, function(x) sum(x > CLQcut))
      cor2ratio <- cor2num / length(secondSNPs)
      cor2numT <- c(cor2ratio > 0.6, 1)
      
      # Identify SNPs in weak LD
      points2 <- min(which(cor2numT > 0))
      NsecondSNPs <- points2:max(secondSNPs)
      reOSNPs <- setdiff(1:max(NsecondSNPs), NsecondSNPs)
      
      if (length(reOSNPs) > 1) {
        print("[INFO] Creating subblocks from rare SNPs before first block")
        subgeno <- Ogeno[, reOSNPs]
        adj_subgeno <- adj_geno[, reOSNPs]
        subSNPinfo <- OSNPinfo[reOSNPs, ]
        subBlocks <- subBigLD(subgeno, subSNPinfo, CLQcut, clstgap, CLQmode, checkLargest)
        subBlocks <- subBlocks + min(reOSNPs) - 1
        expandB <- rbind(expandB, subBlocks)
      }
      
      firstSNPs <- NsecondSNPs
    } else {
      print("[INFO] No rare SNPs before first LD block")
      firstB <- LDblocks[1, ]
      firstSNPs <- which(OSNPinfo[, 2] >= firstB$start.bp & OSNPinfo[, 2] <= firstB$end.bp)
    }
    
    if (nrow(LDblocks) > 1) {
      for (i in 1:(nrow(LDblocks) - 1)) {
        print(paste("[INFO] Processing LD block", i, "of", nrow(LDblocks)))
        secondB <- LDblocks[i + 1, ]
        secondSNPs <- which(OSNPinfo[, 2] >= secondB$start.bp & OSNPinfo[, 2] <= secondB$end.bp)
        OSNPs <- setdiff(max(firstSNPs):min(secondSNPs), c(max(firstSNPs), min(secondSNPs)))
        
        if (length(OSNPs) == 0) {
          expandB <- rbind(expandB, range(firstSNPs))
          firstSNPs <- secondSNPs
        } else {
          # rV² between firstSNPs and OSNPs
          block1_adj <- adj_geno[, c(firstSNPs, OSNPs), drop = FALSE]
          rv1 <- compute_rV2(block1_adj)
          rv1 <- rv1[1:length(firstSNPs), (length(firstSNPs) + 1):ncol(rv1)]
          cor1num <- apply(rv1, 2, function(x) sum(x > CLQcut))
          cor1ratio <- cor1num / length(firstSNPs)
          
          # rV² between secondSNPs and OSNPs
          block2_adj <- adj_geno[, c(secondSNPs, OSNPs), drop = FALSE]
          rv2 <- compute_rV2(block2_adj)
          rv2 <- rv2[1:length(secondSNPs), (length(secondSNPs) + 1):ncol(rv2)]
          cor2num <- apply(rv2, 2, function(x) sum(x > CLQcut))
          cor2ratio <- cor2num / length(secondSNPs)
          
          # Logical thresholds
          cor1numT <- c(1, cor1ratio > 0.6, 0)
          cor2numT <- c(0, cor2ratio > 0.6, 1)
          
          # Define LD extension points
          points1 <- max(firstSNPs) + max(which(cor1numT > 0)) - 1
          NfirstSNPs <- min(firstSNPs):points1
          
          points2 <- max(firstSNPs) + max(which(cor2numT > 0)) - 1
          NsecondSNPs <- points2:max(secondSNPs)
          
          if (max(NfirstSNPs) < min(NsecondSNPs)) {
            expandB <- rbind(expandB, range(NfirstSNPs))
            reOSNPs <- setdiff(min(NfirstSNPs):max(NsecondSNPs), c(NfirstSNPs, NsecondSNPs))
            
            if (length(reOSNPs) > 1) {
              subgeno <- Ogeno[, reOSNPs]
              adj_subgeno <- adj_geno[, reOSNPs]
              subSNPinfo <- OSNPinfo[reOSNPs, ]
              subBlocks <- subBigLD(subgeno, subSNPinfo, CLQcut, clstgap, CLQmode, checkLargest)
              subBlocks <- subBlocks + min(reOSNPs) - 1
              expandB <- rbind(expandB, subBlocks)
            }
            
            firstSNPs <- NsecondSNPs
          } else {
            print("[INFO] Merging two blocks due to overlap")
            subgeno <- Ogeno[, min(firstSNPs):max(secondSNPs)]
            adj_subgeno <- adj_geno[, min(firstSNPs):max(secondSNPs)]
            subSNPinfo <- OSNPinfo[min(firstSNPs):max(secondSNPs), ]
            subBlocks <- subBigLD(subgeno, subSNPinfo, CLQcut, clstgap, CLQmode, checkLargest)
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
    
    # Handle last segment
    if (max(firstSNPs) < (ncol(Ogeno) - 1)) {
      OSNPs <- (max(firstSNPs) + 1):ncol(Ogeno)
      block_adj <- adj_geno[, c(firstSNPs, OSNPs), drop = FALSE]
      rv2 <- compute_rV2(block_adj)
      rv2 <- rv2[1:length(firstSNPs), (length(firstSNPs) + 1):ncol(rv2)]
      
      cor1num <- apply(rv2, 2, function(x) sum(x > CLQcut))
      cor1ratio <- cor1num / length(firstSNPs)
      cor1numT <- c(1, cor1ratio > 0.6, 0)
      points1 <- max(firstSNPs) + max(which(cor1numT > 0)) - 1
      NfirstSNPs <- min(firstSNPs):points1
      
      expandB <- rbind(expandB, range(NfirstSNPs))
      
      reOSNPs <- setdiff(min(NfirstSNPs):ncol(Ogeno), NfirstSNPs)
      if (length(reOSNPs) > 1) {
        subgeno <- Ogeno[, reOSNPs]
        adj_subgeno <- adj_geno[, reOSNPs]
        subSNPinfo <- OSNPinfo[reOSNPs, ]
        subBlocks <- subBigLD(subgeno, subSNPinfo, CLQcut, clstgap, CLQmode, checkLargest)
        subBlocks <- subBlocks + min(reOSNPs) - 1
        expandB <- rbind(expandB, subBlocks)
      }
    } else {
      expandB <- rbind(expandB, range(firstSNPs))
    }
    
    # Final clean-up and return
    expandB <- expandB[expandB[, 1] != expandB[, 2], , drop = FALSE]
    start.bp <- OSNPinfo[, 2][expandB[, 1]]
    end.bp <- OSNPinfo[, 2][expandB[, 2]]
    start.rsID <- as.character(OSNPinfo[, 1][expandB[, 1]])
    end.rsID <- as.character(OSNPinfo[, 1][expandB[, 2]])
    TexpandB <- data.frame(expandB, start.rsID, end.rsID, start.bp, end.bp)
    colnames(TexpandB) <- c("start", "end", "start.rsID", "end.rsID", "start.bp", "end.bp")
    
    print("[INFO] appendSGTs complete")
    return(TexpandB)
  }
  
  #######################################################################################################
  # Main part input data check!!!!!!!!!!!!!!!!!
  Ogeno = geno
  OSNPinfo = SNPinfo
  if (dim(Ogeno)[2] != dim(OSNPinfo)[1]) {
    stop("N of SNPs in geno data and N of SNPs in SNPinfo data Do Not Match!!")
    
  } else if (dim(OSNPinfo)[2] != 2) {
    stop("SNPinfo data Must Contain 2 columns!!")
  }
  Omono = apply(Ogeno, 2, function(x) {
    y<- x[!is.na(x)]
    length(unique(y))!=1
  })
  Ogeno <- Ogeno[,Omono]
  monoSNPs = OSNPinfo[!Omono,]
  OSNPinfo <- OSNPinfo[Omono,]
  maf = apply(Ogeno, 2, function(x) mean(x,na.rm=TRUE)/2)
  maf_ok=ifelse(maf>=0.5,1-maf,maf)
  maf=maf_ok
  mafprun <- which(maf >= MAFcut)
  geno <- Ogeno[,mafprun]
  adj_genoN <- adj_geno[,mafprun]
  SNPinfo <- OSNPinfo[mafprun,]
  
  # print("split whole sequence into subsegments")
  cutpoints.all <- cutsequence.modi(geno, adj_geno, leng, subSegmSize, CLQcut = 0.5)
  cutpoints <- cutpoints.all[[1]]
  atfcut <- (cutpoints.all[[2]])
  if (!is.null(atfcut)){
    atfcut <- sort(atfcut)
  }
  cutpoints = setdiff(cutpoints, 0)
  cutblock <- cbind(c(1, cutpoints + 1), c(cutpoints, dim(geno)[2]))
  cutblock <- cutblock[-(dim(cutblock)[1]), , drop = FALSE]
  LDblocks <- matrix(NA, dim(SNPinfo)[1], 2)
  # partition each segment into LD blocks
  for (i in 1:dim(cutblock)[1]) {
    nowst <- cutblock[i, 1]
    nowed <- cutblock[i, 2]
    subgeno <- geno[, nowst:nowed]
    adj_subgeno <- adj_genoN[, nowst:nowed]
    subSNPinfo <- SNPinfo[nowst:nowed, ]
    # subbinvec <- CLQD(subgeno, subSNPinfo, CLQcut, clstgap, CLQmode = "Density", codechange = FALSE)
    subbinvec <- CLQD(subgeno, subSNPinfo, adj_subgeno, CLQcut, clstgap, CLQmode, codechange = FALSE, checkLargest)
    print('CLQ done!')
    bins <- c(1:max(subbinvec[which(!is.na(subbinvec))]))
    clstlist <- sapply(bins, function(x) which(subbinvec == x), simplify = FALSE)
    clstlist <- lapply(clstlist, sort)  ###
    clstlist <- clstlist[order(sapply(clstlist, min))]  ###
    nowLDblocks <- constructLDblock(clstlist, subSNPinfo)
    # print('constructLDblock done!')
    nowLDblocks <- nowLDblocks + (cutblock[i, 1] - 1)
    nowLDblocks <- nowLDblocks[order(nowLDblocks[, 1]), , drop = FALSE]
    preleng1 <- length(which(!is.na(LDblocks[, 1])))
    LDblocks[(preleng1 + 1):(preleng1 + dim(nowLDblocks)[1]), ] <- nowLDblocks
    print(c(i, dim(cutblock)[1]))
    print(Sys.time())
  }
  doneLDblocks <- LDblocks[which(!is.na(LDblocks[, 1])), , drop = FALSE]
  if (length(atfcut) != 0) {
    newLDblocks <- matrix(NA, dim(SNPinfo)[1], 2)
    consecutive.atf = 0
    for(i in 1:(dim(doneLDblocks)[1]-1)){
      # if(i==1080) break;
      print(paste(i, dim(doneLDblocks)[1]))
      #if(i == 1){
      endblock = doneLDblocks[i,]
      #}
      # if(nowblock[1]<end) next;
      nextblock = doneLDblocks[(i+1),]
      gap = c(endblock[2]:nextblock[1])
      if(length(intersect(gap, atfcut))>0){
        consecutive.atf = consecutive.atf+1
        ## merge 
        if(consecutive.atf>1){
          addlinei = max(which(!is.na(newLDblocks[,1])==TRUE))
          endblock = newLDblocks[addlinei,] 
        }
        nowatfcut = intersect(gap, atfcut)
        newbigblock = range(c(endblock, nextblock))
        newbigblocksize = diff(newbigblock)+1
        nowst = newbigblock[1]
        nowed = newbigblock[2]
        subgeno <- geno[, nowst:nowed]
        adj_subgeno <- adj_genoN[, nowst:nowed]
        subSNPinfo <- SNPinfo[nowst:nowed, ]
        subbinvec <- CLQD(subgeno, subSNPinfo, adj_subgeno, CLQcut, clstgap, CLQmode, codechange = FALSE, checkLargest)
        # print('CLQ done!')
        bins <- c(1:max(subbinvec[which(!is.na(subbinvec))]))
        clstlist <- sapply(bins, function(x) which(subbinvec == x), simplify = FALSE)
        clstlist <- lapply(clstlist, sort)  ###
        clstlist <- clstlist[order(sapply(clstlist, min))]  ###
        nowLDblocks <- constructLDblock(clstlist, SNPinfo[nowst:nowed, ])
        # print('constructLDblock done!')
        nowLDblocks <- nowLDblocks + (newbigblock[1] - 1)
        nowLDblocks <- nowLDblocks[order(nowLDblocks[, 1]), , drop = FALSE]
        preleng1 <- length(which(!is.na(newLDblocks[, 1])))
        nowLDbleng = dim(nowLDblocks)[1]
        #if(nowLDbleng != 1){
        newLDblocks[(preleng1 + 1):(preleng1 + dim(nowLDblocks)[1]), ] <- nowLDblocks
        #}
        #endblock <- nowLDblocks[nowLDbleng, ,drop = FALSE]
        # end <- max(nowLDblocks)
        # if(diff(newbigblock)+1 < subSegmSize)
        print(Sys.time())
      }else{
        consecutive.atf = 0
        addlinei = min(which(is.na(newLDblocks[,1])==TRUE))
        newLDblocks[addlinei,] <-endblock
        #endblock <- nextblock
        if(i == (dim(doneLDblocks)[1]-1)){
          addlinei = min(which(is.na(newLDblocks[,1])==TRUE))
          newLDblocks[addlinei,] <-nextblock
        }
      }
    }
    LDblocks = newLDblocks
  }else{
    LDblocks = doneLDblocks
  }
  LDblocks <- LDblocks[which(!is.na(LDblocks[, 1])), , drop = FALSE]
  LDblocks <- LDblocks[order(LDblocks[, 1]), , drop = FALSE]
  #overlapping LD block merging
  i = 1
  while(i <dim(LDblocks)[1]){
    nowb = LDblocks[i,]
    nextb = LDblocks[(i+1),]
    if(nowb[2]<nextb[1]){
      i = i+1
      next
    }else{
      newb = c(min(c(nowb, nextb)), max(c(nowb, nextb)))
      LDblocks[i,]<-newb
      LDblocks[(i+1),]<-c(NA, NA)
      LDblocks<-LDblocks[!is.na(LDblocks[,1]),]
    }
  }
  start.bp <- SNPinfo[, 2][LDblocks[, 1]]
  end.bp <- SNPinfo[, 2][LDblocks[, 2]]
  start.rsID <- as.character(SNPinfo[, 1][LDblocks[, 1]])
  end.rsID <- as.character(SNPinfo[, 1][LDblocks[, 2]])
  LDblocks <- data.frame(LDblocks, start.rsID, end.rsID, start.bp, end.bp)
  colnames(LDblocks) <- c("start", "end", "start.rsID", "end.rsID", "start.bp", "end.bp")
  if(appendrare==TRUE){
    LDblocks<-appendSGTs(LDblocks, Ogeno, OSNPinfo, CLQcut=CLQcut, clstgap = clstgap, checkLargest)
  }
  #allSNPbp = sort(c(monoSNPs[,2], OSNPinfo[,2]))
  
  # pull out column 2 as a simple numeric vector:
  mono_bp <- as.numeric(monoSNPs[[2]])
  poly_bp <- as.numeric(OSNPinfo[[2]])
  allSNPbp <- sort(c(mono_bp, poly_bp))
  
  LDblocks$start.bp <- as.numeric(LDblocks$start.bp)
  LDblocks$end.bp   <- as.numeric(LDblocks$end.bp)
  
  LDblocks$start <- match(LDblocks$start.bp, allSNPbp)
  LDblocks$end   <- match(LDblocks$end.bp,   allSNPbp)
  
  LDblocks <- LDblocks[order(LDblocks[,1]),]
  
  # ensure start/end indices and bp‐coordinates are in ascending order
  LDblocks[, c("start","end")] <- t(apply(LDblocks[, c("start","end")], 1, sort))
  LDblocks[, c("start.bp","end.bp")] <- t(apply(LDblocks[, c("start.bp","end.bp")], 1, sort))
  
  return(LDblocks)
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
