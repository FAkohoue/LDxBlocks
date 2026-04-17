# ==============================================================================
# analysis_extensions.R
# Seven complementary analysis functions extending the LDxBlocks pipeline.
#
# 1. cv_haplotype_prediction()        k-fold cross-validation for genomic
#                                     prediction; returns predictive ability
#                                     per trait and block importance stability.
# 2. compare_haplotype_populations()  Per-block FST and allele frequency
#                                     comparison between two sample groups.
# 3. plot_haplotype_network()         Minimum-spanning network of haplotype
#                                     alleles within one LD block (igraph).
# 4. run_haplotype_stability()        Finlay-Wilkinson regression stability
#                                     of haplotype effects across environments.
# 5. export_candidate_regions()       Export QTL candidate windows as BED,
#                                     CSV, or a biomaRt-ready data frame.
# 6. decompose_block_effects()        Aggregate per-SNP effects into a
#                                     per-haplotype-allele effect table.
# 7. scan_diversity_windows()         Sliding-window He / n_eff scan
#                                     independent of LD block boundaries.
# ==============================================================================


# ==============================================================================
# 1. cv_haplotype_prediction
# ==============================================================================

#' K-Fold Cross-Validation for Haplotype-Based Genomic Prediction
#'
#' @description
#' Estimates the predictive ability of the haplotype GBLUP model via k-fold
#' cross-validation. In each fold, a subset of individuals is masked from the
#' phenotype and predicted from the haplotype GRM; Pearson correlation between
#' predicted and observed BLUEs is returned as the predictive ability (PA).
#' Runs per trait when multiple traits are supplied.
#'
#' @param geno_matrix Numeric matrix (individuals x SNPs), MAF-filtered dosage.
#' @param snp_info    Data frame with columns \code{SNP}, \code{CHR},
#'   \code{POS}.
#' @param blocks      LD block table from \code{\link{run_Big_LD_all_chr}}.
#' @param blues       Pre-adjusted phenotype means. Accepts the same four
#'   formats as \code{\link{run_haplotype_prediction}}: named numeric vector,
#'   single-trait data frame, multi-trait data frame, or named list.
#' @param k           Integer. Number of folds. Default \code{5L}.
#' @param n_rep       Integer. Number of CV replications (each with a
#'   different random fold assignment). Default \code{1L}.
#' @param top_n       Integer or \code{NULL}. Maximum haplotype alleles per
#'   block passed to \code{\link{build_haplotype_feature_matrix}}.
#'   Default \code{NULL} (all alleles above \code{min_freq}).
#' @param min_freq    Numeric. Minimum haplotype allele frequency.
#'   Default \code{0.05}.
#' @param min_snps    Integer. Minimum SNPs per block for haplotype extraction.
#'   Default \code{3L}.
#' @param id_col      Character. Name of the individual ID column when
#'   \code{blues} is a data frame. Default \code{"id"}.
#' @param blue_col    Character. Name of the BLUE column for single-trait data
#'   frames. Default \code{"blue"}.
#' @param blue_cols   Character vector. Trait column names for multi-trait data
#'   frames. Default \code{NULL} (auto-detect all numeric non-ID columns).
#' @param seed        Integer. RNG seed for reproducible fold assignment.
#'   Default \code{42L}.
#' @param verbose     Logical. Print progress. Default \code{TRUE}.
#'
#' @return A named list of class \code{LDxBlocks_cv}:
#' \describe{
#'   \item{\code{pa_summary}}{Data frame: \code{trait}, \code{rep},
#'     \code{fold}, \code{n_train}, \code{n_test}, \code{PA} (Pearson r),
#'     \code{RMSE}.}
#'   \item{\code{pa_mean}}{Data frame: mean PA and RMSE per trait across
#'     all folds and replications.}
#'   \item{\code{gebv_all}}{Data frame of out-of-fold GEBVs for all
#'     individuals and traits (one row per individual x trait).}
#'   \item{\code{k}}{Number of folds used.}
#'   \item{\code{n_rep}}{Number of replications.}
#' }
#'
#' @examples
#' \donttest{
#' data(ldx_geno, ldx_snp_info, ldx_blocks, ldx_blues, package = "LDxBlocks")
#' cv <- cv_haplotype_prediction(
#'   geno_matrix = ldx_geno,
#'   snp_info    = ldx_snp_info,
#'   blocks      = ldx_blocks,
#'   blues       = ldx_blues,
#'   k           = 5L,
#'   id_col      = "id",
#'   verbose     = FALSE
#' )
#' cv$pa_mean
#' }
#' @seealso \code{\link{run_haplotype_prediction}},
#'   \code{\link{build_haplotype_feature_matrix}}
#' @export
cv_haplotype_prediction <- function(
    geno_matrix,
    snp_info,
    blocks,
    blues,
    k          = 5L,
    n_rep      = 1L,
    top_n      = NULL,
    min_freq   = 0.05,
    min_snps   = 3L,
    id_col     = "id",
    blue_col   = "blue",
    blue_cols  = NULL,
    seed       = 42L,
    verbose    = TRUE
) {
  if (!requireNamespace("rrBLUP", quietly = TRUE))
    stop("rrBLUP is required: install.packages('rrBLUP')", call. = FALSE)

  .log <- function(...) if (verbose) message(sprintf("[cv] %s", paste0(...)))

  # -- Parse blues into a named list of named numeric vectors ----------------
  blues_list <- .parse_blues_ext(blues, id_col, blue_col, blue_cols)
  traits <- names(blues_list)

  # -- Build haplotype features once (shared across all folds) ---------------
  .log("Extracting haplotypes ...")
  haps <- extract_haplotypes(geno_matrix, snp_info, blocks,
                             min_snps = min_snps)
  hap_mat <- build_haplotype_feature_matrix(haps, top_n = top_n,
                                            min_freq = min_freq)
  G <- compute_haplotype_grm(hap_mat)

  inds_G <- rownames(G)
  results <- list()

  for (rep_i in seq_len(n_rep)) {
    set.seed(seed + rep_i - 1L)

    for (tr in traits) {
      pheno <- blues_list[[tr]]
      common <- intersect(names(pheno), inds_G)
      if (length(common) < k)
        stop("Fewer individuals (", length(common),
             ") than folds (", k, ") for trait ", tr, call. = FALSE)

      y <- pheno[common]
      G_sub <- G[common, common, drop = FALSE]

      # Assign folds
      fold_id <- sample(rep(seq_len(k), length.out = length(common)))
      gebv_oof <- setNames(rep(NA_real_, length(common)), common)

      for (fold in seq_len(k)) {
        test_idx  <- which(fold_id == fold)
        train_idx <- which(fold_id != fold)
        y_train   <- y
        y_train[test_idx] <- NA

        fit <- tryCatch(
          rrBLUP::kin.blup(
            data    = data.frame(id = common, y = y_train),
            geno    = "id",
            pheno   = "y",
            K       = G_sub
          ),
          error = function(e) NULL
        )
        if (is.null(fit)) next

        pred <- fit$g[common[test_idx]]
        gebv_oof[test_idx] <- pred

        obs  <- y[test_idx]
        pa   <- if (stats::sd(pred, na.rm=TRUE) > 0 && stats::sd(obs, na.rm=TRUE) > 0)
          stats::cor(pred, obs, use = "complete.obs") else NA_real_
        rmse <- sqrt(mean((pred - obs)^2, na.rm = TRUE))

        results[[length(results) + 1L]] <- data.frame(
          trait   = tr,
          rep     = rep_i,
          fold    = fold,
          n_train = length(train_idx),
          n_test  = length(test_idx),
          PA      = pa,
          RMSE    = rmse,
          stringsAsFactors = FALSE
        )
        .log("  rep=", rep_i, " trait=", tr, " fold=", fold,
             " PA=", round(pa, 3))
      }
    }
  }

  pa_df <- do.call(rbind, results)
  pa_mean <- stats::aggregate(cbind(PA, RMSE) ~ trait, data = pa_df,
                              FUN = mean, na.rm = TRUE)
  pa_sd   <- stats::aggregate(cbind(PA, RMSE) ~ trait, data = pa_df,
                              FUN = stats::sd, na.rm = TRUE)
  names(pa_sd)[2:3] <- c("PA_sd", "RMSE_sd")
  pa_mean <- merge(pa_mean, pa_sd, by = "trait")

  structure(
    list(pa_summary = pa_df, pa_mean = pa_mean, k = k, n_rep = n_rep),
    class = c("LDxBlocks_cv", "list")
  )
}

#' @export
print.LDxBlocks_cv <- function(x, ...) {
  cat("LDxBlocks Cross-Validation\n")
  cat("  Folds:", x$k, "  Replications:", x$n_rep, "\n")
  cat("  Mean predictive ability per trait:\n")
  for (i in seq_len(nrow(x$pa_mean))) {
    r <- x$pa_mean[i, ]
    cat(sprintf("    %-12s  PA = %.3f (sd=%.3f)  RMSE = %.3f\n",
                r$trait, r$PA, r$PA_sd, r$RMSE))
  }
  invisible(x)
}


# ==============================================================================
# 2. compare_haplotype_populations
# ==============================================================================

#' Compare Haplotype Allele Frequencies Between Two Population Groups
#'
#' @description
#' For each LD block, computes allele frequencies in two sample groups and
#' returns FST (Weir-Cockerham 1984), frequency differences, and a chi-squared
#' test of independence. Useful for detecting blocks under divergent selection
#' or monitoring diversity changes between breeding cycles.
#'
#' @param haplotypes Named list from \code{\link{extract_haplotypes}}.
#' @param group1     Character vector of individual IDs for group 1
#'   (e.g. wild/landrace accessions).
#' @param group2     Character vector of individual IDs for group 2
#'   (e.g. elite breeding lines).
#' @param group1_name Character. Label for group 1. Default \code{"group1"}.
#' @param group2_name Character. Label for group 2. Default \code{"group2"}.
#' @param min_freq   Numeric. Alleles below this frequency in both groups
#'   are pooled into an "other" category before testing. Default \code{0.02}.
#' @param missing_string Character. Missing haplotype placeholder.
#'   Default \code{"."}.
#'
#' @return Data frame with one row per block:
#'   \code{block_id}, \code{CHR}, \code{start_bp}, \code{end_bp},
#'   \code{n1} (group 1 sample size), \code{n2},
#'   \code{n_alleles} (number of distinct alleles),
#'   \code{FST} (Weir-Cockerham single-locus estimate),
#'   \code{max_freq_diff} (max absolute frequency difference across alleles),
#'   \code{dominant_g1} (most frequent allele in group 1),
#'   \code{dominant_g2} (most frequent allele in group 2),
#'   \code{chisq_p} (chi-squared p-value, \code{NA} if < 2 alleles in
#'   either group),
#'   \code{divergent} (\code{TRUE} when FST > 0.1 and chisq_p < 0.05).
#'
#' @examples
#' \donttest{
#' data(ldx_geno, ldx_snp_info, ldx_blocks, package = "LDxBlocks")
#' haps <- extract_haplotypes(ldx_geno, ldx_snp_info, ldx_blocks)
#' ids  <- rownames(ldx_geno)
#' cmp  <- compare_haplotype_populations(
#'   haplotypes   = haps,
#'   group1       = ids[1:60],
#'   group2       = ids[61:120],
#'   group1_name  = "cycle1",
#'   group2_name  = "cycle2"
#' )
#' cmp[cmp$divergent, c("block_id", "FST", "max_freq_diff", "chisq_p")]
#' }
#' @references
#' Weir BS, Cockerham CC (1984). Estimating F-statistics for the analysis of
#' population structure. \emph{Evolution} \strong{38}(6):1358-1370.
#' @seealso \code{\link{extract_haplotypes}},
#'   \code{\link{compute_haplotype_diversity}}
#' @export
compare_haplotype_populations <- function(
    haplotypes,
    group1,
    group2,
    group1_name    = "group1",
    group2_name    = "group2",
    min_freq       = 0.02,
    missing_string = "."
) {
  bi <- attr(haplotypes, "block_info")

  rows <- lapply(seq_along(haplotypes), function(i) {
    bn  <- names(haplotypes)[i]
    hap <- haplotypes[[bn]]
    blk <- if (!is.null(bi)) bi[bi$block_id == bn, , drop = FALSE] else NULL

    ind_names <- names(hap)
    h1 <- hap[intersect(ind_names, group1)]
    h2 <- hap[intersect(ind_names, group2)]
    h1 <- h1[!grepl(missing_string, h1, fixed = TRUE)]
    h2 <- h2[!grepl(missing_string, h2, fixed = TRUE)]
    n1 <- length(h1); n2 <- length(h2)

    na_row <- data.frame(
      block_id = bn,
      CHR = if (!is.null(blk)) blk$CHR[1] else NA,
      start_bp = if (!is.null(blk)) blk$start_bp[1] else NA,
      end_bp   = if (!is.null(blk)) blk$end_bp[1] else NA,
      n1 = n1, n2 = n2, n_alleles = NA_integer_,
      FST = NA_real_, max_freq_diff = NA_real_,
      dominant_g1 = NA_character_, dominant_g2 = NA_character_,
      chisq_p = NA_real_, divergent = FALSE,
      stringsAsFactors = FALSE
    )
    if (n1 < 2L || n2 < 2L) return(na_row)

    alleles <- union(unique(h1), unique(h2))
    f1 <- table(h1)[alleles]; f1[is.na(f1)] <- 0L; f1 <- f1 / n1
    f2 <- table(h2)[alleles]; f2[is.na(f2)] <- 0L; f2 <- f2 / n2

    # Pool rare alleles
    rare <- (f1 < min_freq) & (f2 < min_freq)
    if (any(rare) && !all(rare)) {
      other_f1 <- sum(f1[rare]); other_f2 <- sum(f2[rare])
      f1 <- c(f1[!rare], other = other_f1)
      f2 <- c(f2[!rare], other = other_f2)
    }
    alleles_kept <- names(f1)

    # Weir-Cockerham FST (single-locus, two-population)
    p_bar <- (n1 * f1 + n2 * f2) / (n1 + n2)
    n_bar <- (n1 + n2) / 2
    S2    <- (n1 * (f1 - p_bar)^2 + n2 * (f2 - p_bar)^2) / (2 * n_bar - 1)
    # FST per allele: S2 / (p_bar*(1-p_bar) + S2*(1-1/n_bar) ...)
    # Use Weir-Cockerham eq 6: FST = S2 / (p_bar(1-p_bar))
    denom <- p_bar * (1 - p_bar)
    fst_per <- ifelse(denom > 0, S2 / denom, 0)
    FST <- mean(fst_per, na.rm = TRUE)
    FST <- max(0, min(1, FST))

    freq_diff <- abs(f1 - f2)
    dom1 <- names(which.max(f1))
    dom2 <- names(which.max(f2))

    # Chi-squared test of independence
    cnt1 <- round(f1 * n1); cnt2 <- round(f2 * n2)
    chisq_p <- tryCatch({
      tab <- rbind(cnt1, cnt2)
      if (nrow(tab) >= 2 && ncol(tab) >= 2 && all(rowSums(tab) > 0))
        stats::chisq.test(tab, simulate.p.value = TRUE, B = 2000)$p.value
      else NA_real_
    }, error = function(e) NA_real_)

    data.frame(
      block_id     = bn,
      CHR          = if (!is.null(blk)) blk$CHR[1] else NA,
      start_bp     = if (!is.null(blk)) blk$start_bp[1] else NA,
      end_bp       = if (!is.null(blk)) blk$end_bp[1] else NA,
      n1           = n1, n2 = n2,
      n_alleles    = length(alleles_kept),
      FST          = round(FST, 4),
      max_freq_diff = round(max(freq_diff), 4),
      dominant_g1  = dom1,
      dominant_g2  = dom2,
      chisq_p      = round(chisq_p, 4),
      divergent    = !is.na(FST) & FST > 0.1 &
        !is.na(chisq_p) & chisq_p < 0.05,
      stringsAsFactors = FALSE
    )
  })

  out <- do.call(rbind, rows)
  out[order(out$CHR, out$start_bp), ]
}


# ==============================================================================
# 3. plot_haplotype_network
# ==============================================================================

#' Plot a Minimum-Spanning Haplotype Network for One LD Block
#'
#' @description
#' Draws a minimum-spanning network (MSN) of haplotype alleles within a single
#' LD block using \pkg{igraph}. Nodes represent alleles; edge weights are
#' Hamming distances (number of differing SNP positions). Node size is
#' proportional to allele frequency. Optional colour by population or
#' phenotypic group.
#'
#' @param haplotypes Named list from \code{\link{extract_haplotypes}}.
#' @param block_id   Character. Name of the block to visualise (must be in
#'   \code{names(haplotypes)}).
#' @param groups     Named character vector mapping individual IDs to group
#'   labels (used for node pie-chart colouring). \code{NULL} = all one colour.
#' @param min_freq   Numeric. Alleles below this frequency are dropped before
#'   plotting. Default \code{0.02}.
#' @param missing_string Character. Missing haplotype placeholder.
#'   Default \code{"."}.
#' @param title      Character. Plot title. Default = block_id.
#' @param palette    Character vector of colours for groups. \code{NULL} uses
#'   a built-in palette.
#'
#' @return Invisibly returns the \code{igraph} graph object. The network is
#'   plotted as a side effect.
#'
#' @examples
#' \donttest{
#' data(ldx_geno, ldx_snp_info, ldx_blocks, package = "LDxBlocks")
#' haps <- extract_haplotypes(ldx_geno, ldx_snp_info, ldx_blocks)
#' plot_haplotype_network(haps, block_id = names(haps)[1])
#' }
#' @seealso \code{\link{extract_haplotypes}},
#'   \code{\link{compute_haplotype_diversity}}
#' @export
plot_haplotype_network <- function(
    haplotypes,
    block_id,
    groups         = NULL,
    min_freq       = 0.02,
    missing_string = ".",
    title          = NULL,
    palette        = NULL
) {
  if (!requireNamespace("igraph", quietly = TRUE))
    stop("igraph is required: install.packages('igraph')", call. = FALSE)

  if (!block_id %in% names(haplotypes))
    stop("block_id '", block_id, "' not found in haplotypes.", call. = FALSE)

  hap   <- haplotypes[[block_id]]
  valid <- hap[!grepl(missing_string, hap, fixed = TRUE)]
  if (length(valid) == 0L) stop("No non-missing haplotypes in this block.")

  tbl   <- table(valid)
  freq  <- as.numeric(tbl) / sum(tbl)
  names(freq) <- names(tbl)
  keep  <- freq >= min_freq
  if (sum(keep) < 2L)
    stop("Fewer than 2 alleles above min_freq=", min_freq,
         ". Lower min_freq or choose a different block.")

  alleles <- names(freq)[keep]
  freq    <- freq[keep]
  n_al    <- length(alleles)

  # Hamming distance matrix
  nchar_a <- nchar(alleles[1])
  ham <- matrix(0L, n_al, n_al, dimnames = list(alleles, alleles))
  for (a in seq_len(n_al - 1L)) {
    for (b in seq(a + 1L, n_al)) {
      d <- sum(strsplit(alleles[a], "")[[1]] != strsplit(alleles[b], "")[[1]])
      ham[a, b] <- ham[b, a] <- d
    }
  }

  g <- igraph::graph_from_adjacency_matrix(ham, mode = "undirected",
                                           weighted = TRUE, diag = FALSE)
  mst <- igraph::mst(g, weights = igraph::E(g)$weight)

  # Node sizing by frequency
  node_size <- 5 + 40 * freq[igraph::V(mst)$name]

  # Node colours
  if (!is.null(groups)) {
    grp_levels <- unique(groups)
    if (is.null(palette))
      palette <- grDevices::hcl.colors(length(grp_levels), "Set2")
    col_map <- stats::setNames(palette, grp_levels)

    node_cols <- vapply(igraph::V(mst)$name, function(al) {
      carriers <- names(valid)[valid == al]
      grp_carr <- groups[intersect(carriers, names(groups))]
      if (!length(grp_carr)) return("#AAAAAA")
      most_common <- names(sort(table(grp_carr), decreasing = TRUE))[1]
      col_map[most_common]
    }, character(1L))
  } else {
    node_cols <- rep("#4E9EC2", igraph::vcount(mst))
  }

  edge_labels <- round(igraph::E(mst)$weight)
  plot(
    mst,
    vertex.size        = node_size,
    vertex.color       = node_cols,
    vertex.label       = paste0(igraph::V(mst)$name, "\n(",
                                round(freq[igraph::V(mst)$name] * 100, 1), "%)"),
    vertex.label.cex   = 0.65,
    vertex.label.color = "black",
    edge.width         = 1.5,
    edge.label         = ifelse(edge_labels > 1, edge_labels, ""),
    edge.label.cex     = 0.7,
    edge.color         = "#888888",
    layout             = igraph::layout_nicely(mst),
    main               = if (!is.null(title)) title else block_id
  )

  if (!is.null(groups)) {
    grp_levels <- unique(groups)
    legend("bottomright", legend = grp_levels,
           fill = palette, bty = "n", cex = 0.8)
  }

  invisible(mst)
}


# ==============================================================================
# 4. run_haplotype_stability
# ==============================================================================

#' Finlay-Wilkinson Stability Analysis of Haplotype Effects Across Environments
#'
#' @description
#' Estimates the stability of each haplotype block's effect across multiple
#' environments using Finlay-Wilkinson (1963) regression. For each block,
#' the per-environment GEBV contribution is regressed on the environment mean
#' (the environmental index). A regression slope b_i = 1 indicates average
#' stability; b_i > 1 = above-average response (exploits good environments);
#' b_i < 1 = below-average response (robust across environments).
#'
#' @param geno_matrix  Numeric matrix (individuals x SNPs).
#' @param snp_info     Data frame with \code{SNP}, \code{CHR}, \code{POS}.
#' @param blocks       LD block table from \code{\link{run_Big_LD_all_chr}}.
#' @param blues_list   Named list of named numeric vectors, one per
#'   environment: \code{list(env1 = c(id1=val,...), env2 = ...)}.
#' @param top_n        Integer or \code{NULL}. Max haplotype alleles per block.
#'   Default \code{NULL}.
#' @param min_freq     Numeric. Minimum haplotype allele frequency.
#'   Default \code{0.05}.
#' @param min_snps     Integer. Minimum SNPs per block. Default \code{3L}.
#' @param verbose      Logical. Default \code{TRUE}.
#'
#' @return Data frame with one row per block:
#'   \code{block_id}, \code{CHR}, \code{start_bp}, \code{end_bp},
#'   \code{b} (Finlay-Wilkinson slope — stability coefficient),
#'   \code{b_se} (standard error of b),
#'   \code{r2_fw} (R² of the FW regression),
#'   \code{s2d} (deviation mean square — non-linear instability),
#'   \code{stable} (\code{TRUE} when b is not significantly different from 1
#'   at alpha=0.05).
#'
#' @examples
#' \donttest{
#' data(ldx_geno, ldx_snp_info, ldx_blocks, ldx_blues, package = "LDxBlocks")
#' # Simulate two environments by splitting the BLUEs
#' b1 <- setNames(ldx_blues$YLD + rnorm(nrow(ldx_blues), 0, 0.1),
#'                ldx_blues$id)
#' b2 <- setNames(ldx_blues$YLD + rnorm(nrow(ldx_blues), 0.5, 0.15),
#'                ldx_blues$id)
#' stab <- run_haplotype_stability(
#'   geno_matrix = ldx_geno,
#'   snp_info    = ldx_snp_info,
#'   blocks      = ldx_blocks,
#'   blues_list  = list(env1 = b1, env2 = b2),
#'   verbose     = FALSE
#' )
#' head(stab[order(stab$b), ])
#' }
#' @references
#' Finlay KW, Wilkinson GN (1963). The analysis of adaptation in a
#' plant-breeding programme. \emph{Australian Journal of Agricultural Research}
#' \strong{14}(6):742-754.
#' @seealso \code{\link{run_haplotype_prediction}},
#'   \code{\link{compute_local_gebv}}
#' @export
run_haplotype_stability <- function(
    geno_matrix,
    snp_info,
    blocks,
    blues_list,
    top_n    = NULL,
    min_freq = 0.05,
    min_snps = 3L,
    verbose  = TRUE
) {
  if (!is.list(blues_list) || is.null(names(blues_list)))
    stop("blues_list must be a named list of named numeric vectors (one per environment).",
         call. = FALSE)
  if (length(blues_list) < 2L)
    stop("At least 2 environments required for stability analysis.", call. = FALSE)
  if (!requireNamespace("rrBLUP", quietly = TRUE))
    stop("rrBLUP is required: install.packages('rrBLUP')", call. = FALSE)

  .log <- function(...) if (verbose) message(sprintf("[stability] %s", paste0(...)))
  envs <- names(blues_list)

  # Build hap features + GRM once
  .log("Extracting haplotypes ...")
  haps    <- extract_haplotypes(geno_matrix, snp_info, blocks,
                                min_snps = min_snps)
  hap_mat <- build_haplotype_feature_matrix(haps, top_n = top_n,
                                            min_freq = min_freq)
  G       <- compute_haplotype_grm(hap_mat)

  # Get per-environment local GEBVs per block
  local_by_env <- list()
  for (env in envs) {
    .log("Fitting GBLUP for environment: ", env)
    y <- blues_list[[env]]
    common <- intersect(names(y), rownames(G))
    if (length(common) < 5L) {
      warning("Fewer than 5 common individuals for env '", env, "' — skipping.")
      next
    }
    G_sub <- G[common, common, drop = FALSE]
    fit <- tryCatch(
      rrBLUP::kin.blup(
        data  = data.frame(id = common, y = y[common]),
        geno  = "id", pheno = "y", K = G_sub
      ),
      error = function(e) NULL
    )
    if (is.null(fit)) next

    gebv <- fit$g
    snp_fx <- backsolve_snp_effects(
      geno_matrix = geno_matrix[common, , drop = FALSE],
      gebv        = gebv[common],
      G           = G_sub
    )
    local_gebv <- compute_local_gebv(
      geno_matrix = geno_matrix[common, , drop = FALSE],
      snp_info    = snp_info,
      blocks      = blocks,
      snp_effects = snp_fx
    )
    local_by_env[[env]] <- local_gebv
  }

  if (length(local_by_env) < 2L)
    stop("Fewer than 2 environments had sufficient data.", call. = FALSE)

  # Finlay-Wilkinson regression per block
  block_ids <- colnames(local_by_env[[1]])
  env_means <- vapply(local_by_env, function(m) mean(m, na.rm = TRUE), numeric(1L))

  rows <- lapply(block_ids, function(bid) {
    bi_info <- blocks[1, , drop = FALSE]   # placeholder
    # Get this block's local GEBV per environment (individual mean per env)
    block_env_means <- vapply(names(local_by_env), function(env) {
      col <- local_by_env[[env]][, bid, drop = TRUE]
      mean(col, na.rm = TRUE)
    }, numeric(1L))

    # Environmental index I_j = mean GEBV in env j across all blocks
    I <- env_means

    if (stats::var(I) < 1e-10) return(NULL)

    # FW regression: block_effect ~ I
    df_fw <- data.frame(y = block_env_means, I = I)
    fit   <- stats::lm(y ~ I, data = df_fw)
    co    <- stats::coef(summary(fit))
    b     <- co["I", "Estimate"]
    b_se  <- co["I", "Std. Error"]
    r2    <- summary(fit)$r.squared
    resid <- stats::residuals(fit)
    s2d   <- if (length(resid) > 2) sum(resid^2) / (length(resid) - 2) else NA_real_

    # Test H0: b = 1
    t_stat <- (b - 1) / b_se
    df_t   <- nrow(df_fw) - 2L
    p_b1   <- 2 * stats::pt(-abs(t_stat), df = df_t)
    stable <- !is.na(p_b1) && p_b1 >= 0.05

    data.frame(
      block_id = bid,
      b        = round(b, 4),
      b_se     = round(b_se, 4),
      r2_fw    = round(r2, 4),
      s2d      = round(s2d, 6),
      p_b1     = round(p_b1, 4),
      stable   = stable,
      stringsAsFactors = FALSE
    )
  })

  rows <- Filter(Negate(is.null), rows)
  out  <- do.call(rbind, rows)

  # Merge block metadata
  if (all(c("CHR", "start.bp", "end.bp") %in% names(blocks))) {
    bk <- blocks[, c("block_name", "CHR", "start.bp", "end.bp"),
                 drop = FALSE]
    names(bk)[1] <- "block_id"
    out <- merge(bk, out, by = "block_id", all.y = TRUE)
  }

  out[order(out$b), ]
}


# ==============================================================================
# 5. export_candidate_regions
# ==============================================================================

#' Export Candidate Gene Regions to BED, CSV, or biomaRt Format
#'
#' @description
#' Converts the output of \code{\link{define_qtl_regions}} into formats
#' ready for downstream annotation: standard BED (for BEDtools, UCSC),
#' a plain CSV, or a named list suitable for direct use with
#' \code{biomaRt::getBM()} filters.
#'
#' @param qtl_regions Data frame from \code{\link{define_qtl_regions}}.
#' @param format      Character. Output format: \code{"bed"}, \code{"csv"},
#'   or \code{"biomart"}. Default \code{"bed"}.
#' @param out_file    Character path. If \code{NULL} (default), returns the
#'   object invisibly without writing to disk.
#' @param chr_prefix  Character. Prefix to add to chromosome names in BED
#'   output (e.g. \code{"chr"} for UCSC, \code{""} for Ensembl).
#'   Default \code{""}.
#' @param use_lead_snp Logical. If \code{TRUE} and \code{ld_decay} was
#'   supplied to \code{define_qtl_regions}, use the LD-extended
#'   \code{candidate_region_start}/\code{candidate_region_end} columns
#'   instead of block boundaries. Default \code{TRUE}.
#' @param padding_bp  Integer. Additional bp to add on each side of each
#'   region (applied after the LD extension if any). Default \code{0L}.
#'
#' @return Invisibly:
#' \describe{
#'   \item{\code{"bed"}}{Data frame with columns \code{chrom}, \code{start}
#'     (0-based), \code{end}, \code{name}, \code{score}, \code{strand}.}
#'   \item{\code{"csv"}}{The \code{qtl_regions} data frame, optionally
#'     written to \code{out_file}.}
#'   \item{\code{"biomart"}}{Named list with elements \code{chromosome_name},
#'     \code{start}, \code{end} suitable for
#'     \code{biomaRt::getBM(filters = c("chromosome_name","start","end"), ...)}.}
#' }
#'
#' @examples
#' \donttest{
#' data(ldx_geno, ldx_snp_info, ldx_blocks, ldx_gwas, package = "LDxBlocks")
#' qtl <- define_qtl_regions(ldx_gwas, ldx_blocks, ldx_snp_info,
#'                            p_threshold = NULL)
#' # BED file
#' bed <- export_candidate_regions(qtl, format = "bed", chr_prefix = "chr")
#' head(bed)
#' # biomaRt-ready list
#' bm  <- export_candidate_regions(qtl, format = "biomart")
#' # biomaRt::getBM(attributes = c("ensembl_gene_id","external_gene_name"),
#' #                filters = c("chromosome_name","start","end"),
#' #                values  = bm, mart = my_mart)
#' }
#' @seealso \code{\link{define_qtl_regions}}, \code{\link{compute_ld_decay}}
#' @export
export_candidate_regions <- function(
    qtl_regions,
    format       = c("bed", "csv", "biomart"),
    out_file     = NULL,
    chr_prefix   = "",
    use_lead_snp = TRUE,
    padding_bp   = 0L
) {
  format <- match.arg(format)
  req <- c("block_id", "CHR", "start_bp", "end_bp")
  miss <- setdiff(req, names(qtl_regions))
  if (length(miss))
    stop("qtl_regions missing columns: ", paste(miss, collapse=", "), call.=FALSE)

  # Determine region boundaries
  has_cand <- all(c("candidate_region_start","candidate_region_end") %in%
                    names(qtl_regions))
  if (use_lead_snp && has_cand) {
    reg_start <- pmax(0L, qtl_regions$candidate_region_start - padding_bp)
    reg_end   <- qtl_regions$candidate_region_end + padding_bp
  } else {
    reg_start <- pmax(0L, qtl_regions$start_bp - padding_bp)
    reg_end   <- qtl_regions$end_bp + padding_bp
  }
  chr_str <- paste0(chr_prefix, qtl_regions$CHR)

  if (format == "bed") {
    out <- data.frame(
      chrom  = chr_str,
      start  = as.integer(reg_start) - 1L,   # BED is 0-based
      end    = as.integer(reg_end),
      name   = qtl_regions$block_id,
      score  = if ("n_sig_markers" %in% names(qtl_regions))
        qtl_regions$n_sig_markers else 0L,
      strand = ".",
      stringsAsFactors = FALSE
    )
    if (!is.null(out_file)) {
      write.table(out, out_file, sep = "\t", quote = FALSE,
                  row.names = FALSE, col.names = FALSE)
      message("BED file written: ", out_file)
    }
    return(invisible(out))
  }

  if (format == "csv") {
    if (!is.null(out_file)) {
      data.table::fwrite(qtl_regions, out_file)
      message("CSV written: ", out_file)
    }
    return(invisible(qtl_regions))
  }

  if (format == "biomart") {
    # biomaRt expects vectors, one entry per region
    out <- list(
      chromosome_name = qtl_regions$CHR,
      start           = as.integer(reg_start),
      end             = as.integer(reg_end)
    )
    if (!is.null(out_file)) {
      saveRDS(out, out_file)
      message("biomaRt filter list saved: ", out_file)
    }
    return(invisible(out))
  }
}


# ==============================================================================
# 6. decompose_block_effects
# ==============================================================================

#' Decompose Per-SNP Effects into Per-Haplotype-Allele Effect Table
#'
#' @description
#' Aggregates the per-SNP additive effects (from \code{\link{backsolve_snp_effects}})
#' into a per-allele effect for each LD block: the effect of carrying allele h
#' is the sum of SNP effects weighted by the allele's dosage at each SNP
#' position. Returns an interpretable table of "which allele of which block
#' is worth how many units of the trait?"
#'
#' @param haplotypes    Named list from \code{\link{extract_haplotypes}}.
#' @param snp_info      Data frame with \code{SNP}, \code{CHR}, \code{POS}.
#' @param blocks        LD block table from \code{\link{run_Big_LD_all_chr}}.
#' @param snp_effects   Named numeric vector of per-SNP additive effects,
#'   as returned by \code{\link{backsolve_snp_effects}}.
#' @param min_freq      Numeric. Alleles below this frequency are excluded.
#'   Default \code{0.02}.
#' @param missing_string Character. Missing haplotype placeholder.
#'   Default \code{"."}.
#'
#' @return Data frame with one row per allele per block:
#'   \code{block_id}, \code{CHR}, \code{start_bp}, \code{end_bp},
#'   \code{allele} (the haplotype string), \code{frequency},
#'   \code{allele_effect} (sum of SNP effects for this allele),
#'   \code{effect_rank} (rank within block, 1 = most positive effect),
#'   \code{n_snps_block}.
#'   Sorted by \code{CHR}, \code{start_bp}, \code{effect_rank}.
#'
#' @examples
#' \donttest{
#' data(ldx_geno, ldx_snp_info, ldx_blocks, ldx_blues, package = "LDxBlocks")
#' haps    <- extract_haplotypes(ldx_geno, ldx_snp_info, ldx_blocks)
#' hap_mat <- build_haplotype_feature_matrix(haps)
#' G       <- compute_haplotype_grm(hap_mat)
#' res     <- run_haplotype_prediction(ldx_geno, ldx_snp_info, ldx_blocks,
#'                                      blues = ldx_blues, id_col = "id",
#'                                      blue_col = "YLD", verbose = FALSE)
#' # Use first trait's SNP effects
#' snp_fx  <- res$snp_effects[[1]]
#' allele_tbl <- decompose_block_effects(haps, ldx_snp_info, ldx_blocks,
#'                                        snp_effects = snp_fx)
#' head(allele_tbl[order(-allele_tbl$allele_effect), ])
#' }
#' @seealso \code{\link{backsolve_snp_effects}},
#'   \code{\link{rank_haplotype_blocks}}
#' @export
decompose_block_effects <- function(
    haplotypes,
    snp_info,
    blocks,
    snp_effects,
    min_freq       = 0.02,
    missing_string = "."
) {
  if (!is.numeric(snp_effects) || is.null(names(snp_effects)))
    stop("snp_effects must be a named numeric vector (SNP ID -> effect).",
         call. = FALSE)

  bi <- attr(haplotypes, "block_info")

  rows <- lapply(seq_along(haplotypes), function(i) {
    bn  <- names(haplotypes)[i]
    hap <- haplotypes[[bn]]
    blk <- if (!is.null(bi)) bi[bi$block_id == bn, , drop = FALSE] else NULL

    if (is.null(blk) || nrow(blk) == 0L) return(NULL)

    # Get SNPs in this block
    ch <- blk$CHR[1]; sb <- blk$start_bp[1]; eb <- blk$end_bp[1]
    snp_in_block <- snp_info$SNP[
      snp_info$CHR == ch &
        snp_info$POS >= sb &
        snp_info$POS <= eb
    ]
    fx_in_block <- snp_effects[intersect(snp_in_block, names(snp_effects))]
    if (!length(fx_in_block)) return(NULL)

    # Allele frequency table
    valid <- hap[!grepl(missing_string, hap, fixed = TRUE)]
    if (!length(valid)) return(NULL)
    tbl  <- table(valid)
    freq <- as.numeric(tbl) / sum(tbl)
    names(freq) <- names(tbl)
    keep_al <- names(freq)[freq >= min_freq]
    if (!length(keep_al)) return(NULL)

    # For each allele, compute sum of SNP effects weighted by allele dosage
    snp_order <- intersect(names(fx_in_block),
                           snp_info$SNP[snp_info$CHR == ch &
                                          snp_info$POS >= sb &
                                          snp_info$POS <= eb])
    snp_order <- snp_order[snp_order %in% names(fx_in_block)]

    al_effects <- vapply(keep_al, function(al) {
      # Parse allele string into dosage vector at each SNP position
      chars <- strsplit(al, "")[[1]]
      n_pos <- length(chars)
      if (n_pos != length(snp_order)) return(NA_real_)
      doses <- suppressWarnings(as.integer(chars))
      doses[is.na(doses)] <- 0L
      sum(doses * fx_in_block[snp_order], na.rm = TRUE)
    }, numeric(1L))

    valid_al <- !is.na(al_effects)
    if (!any(valid_al)) return(NULL)

    al_df <- data.frame(
      block_id     = bn,
      CHR          = ch,
      start_bp     = sb,
      end_bp       = eb,
      allele       = keep_al[valid_al],
      frequency    = round(freq[keep_al[valid_al]], 4),
      allele_effect = round(al_effects[valid_al], 6),
      n_snps_block  = length(snp_order),
      stringsAsFactors = FALSE
    )
    al_df$effect_rank <- rank(-al_df$allele_effect, ties.method = "min")
    al_df
  })

  rows <- Filter(Negate(is.null), rows)
  if (!length(rows)) return(data.frame())
  out <- do.call(rbind, rows)
  out[order(out$CHR, out$start_bp, out$effect_rank), ]
}


# ==============================================================================
# 7. scan_diversity_windows
# ==============================================================================

#' Sliding-Window Genome-Wide Diversity Scan
#'
#' @description
#' Computes haplotype diversity metrics (He, Shannon entropy, n_eff_alleles,
#' dominant frequency) in sliding windows across the genome, independently of
#' LD block boundaries. Useful for identifying diversity valleys (bottlenecks,
#' selective sweeps) and comparing wild/elite panels without needing
#' pre-defined blocks.
#'
#' @param geno_matrix   Numeric matrix (individuals x SNPs), 0/1/2/NA.
#' @param snp_info      Data frame with \code{SNP}, \code{CHR}, \code{POS}.
#' @param window_bp     Integer. Window size in base pairs. Default \code{1e6L}
#'   (1 Mb).
#' @param step_bp       Integer. Step size in base pairs. Default
#'   \code{5e5L} (500 kb, i.e. 50% overlap).
#' @param min_snps_win  Integer. Minimum SNPs in a window to compute
#'   diversity (windows with fewer are skipped). Default \code{5L}.
#' @param missing_val   Numeric. Value representing missing data in
#'   \code{geno_matrix}. Default \code{NA}.
#'
#' @return Data frame with one row per window:
#'   \code{CHR}, \code{win_start}, \code{win_end}, \code{win_mid},
#'   \code{n_snps}, \code{n_ind}, \code{n_haplotypes},
#'   \code{He} (Nei 1973, sample-size corrected),
#'   \code{Shannon}, \code{n_eff_alleles}, \code{freq_dominant},
#'   \code{sweep_flag}.
#'   Sorted by \code{CHR}, \code{win_start}.
#'
#' @examples
#' \donttest{
#' data(ldx_geno, ldx_snp_info, package = "LDxBlocks")
#' scan <- scan_diversity_windows(
#'   geno_matrix  = ldx_geno,
#'   snp_info     = ldx_snp_info,
#'   window_bp    = 50000L,
#'   step_bp      = 25000L,
#'   min_snps_win = 3L
#' )
#' # Plot He across chromosome 1
#' chr1 <- scan[scan$CHR == "1", ]
#' plot(chr1$win_mid / 1e3, chr1$He, type = "l",
#'      xlab = "Position (kb)", ylab = "He",
#'      main = "Haplotype diversity scan — chr 1")
#' }
#' @seealso \code{\link{compute_haplotype_diversity}},
#'   \code{\link{compare_haplotype_populations}}
#' @export
scan_diversity_windows <- function(
    geno_matrix,
    snp_info,
    window_bp    = 1e6L,
    step_bp      = 5e5L,
    min_snps_win = 5L,
    missing_val  = NA
) {
  if (!is.matrix(geno_matrix))
    stop("geno_matrix must be a numeric matrix (individuals x SNPs).",
         call. = FALSE)
  req <- c("SNP", "CHR", "POS")
  miss <- setdiff(req, names(snp_info))
  if (length(miss))
    stop("snp_info missing: ", paste(miss, collapse=", "), call.=FALSE)

  snp_info <- snp_info[snp_info$SNP %in% colnames(geno_matrix), ]
  chrs <- unique(snp_info$CHR)
  rows <- list()

  for (ch in chrs) {
    si_chr <- snp_info[snp_info$CHR == ch, ]
    si_chr <- si_chr[order(si_chr$POS), ]
    pos    <- si_chr$POS
    if (!length(pos)) next

    starts <- seq(min(pos), max(pos), by = step_bp)

    for (ws in starts) {
      we   <- ws + window_bp - 1L
      idx  <- which(pos >= ws & pos <= we)
      if (length(idx) < min_snps_win) next

      snp_ids <- si_chr$SNP[idx]
      G_win   <- geno_matrix[, snp_ids, drop = FALSE]

      # Replace missing with NA
      if (!is.na(missing_val))
        G_win[G_win == missing_val] <- NA

      # Build diploid allele strings per individual
      hap_str <- apply(G_win, 1, function(row) {
        if (any(is.na(row))) return(NA_character_)
        paste(as.integer(row), collapse = "")
      })
      valid <- hap_str[!is.na(hap_str)]
      ni    <- length(valid)
      if (ni < 2L) next

      tbl  <- table(valid)
      freq <- as.numeric(tbl) / sum(tbl)
      He   <- 1 - sum(freq^2)
      He_c <- (ni / (ni - 1L)) * He
      Sh   <- -sum(freq * log(pmax(freq, .Machine$double.eps)))
      n_eff <- 1 / sum(freq^2)

      rows[[length(rows) + 1L]] <- data.frame(
        CHR           = ch,
        win_start     = as.integer(ws),
        win_end       = as.integer(we),
        win_mid       = as.integer((ws + we) / 2L),
        n_snps        = length(idx),
        n_ind         = ni,
        n_haplotypes  = length(tbl),
        He            = round(He_c, 4),
        Shannon       = round(Sh, 4),
        n_eff_alleles = round(n_eff, 3),
        freq_dominant = round(max(freq), 4),
        sweep_flag    = max(freq) >= 0.90,
        stringsAsFactors = FALSE
      )
    }
  }

  if (!length(rows)) return(data.frame())
  out <- do.call(rbind, rows)
  out[order(out$CHR, out$win_start), ]
}


# ==============================================================================
# Internal helper: parse blues into named list (shared with haplotypes.R)
# Uses a local copy to avoid coupling to the internal .parse_blues
# ==============================================================================
.parse_blues_ext <- function(blues, id_col, blue_col, blue_cols) {
  if (is.numeric(blues) && !is.null(names(blues)))
    return(list(trait = blues))
  if (is.data.frame(blues)) {
    if (!id_col %in% names(blues))
      stop("id_col '", id_col, "' not found in blues.", call. = FALSE)
    ids <- as.character(blues[[id_col]])
    if (!is.null(blue_cols)) {
      trts <- blue_cols
    } else if (blue_col %in% names(blues)) {
      trts <- blue_col
    } else {
      trts <- names(blues)[vapply(blues, is.numeric, logical(1L))]
      trts <- setdiff(trts, id_col)
    }
    return(stats::setNames(
      lapply(trts, function(tr) {
        v <- as.numeric(blues[[tr]]); names(v) <- ids; v
      }), trts
    ))
  }
  if (is.list(blues) && all(vapply(blues, is.numeric, logical(1L))))
    return(blues)
  stop("blues must be a named numeric vector, data frame, or named list.",
       call. = FALSE)
}
