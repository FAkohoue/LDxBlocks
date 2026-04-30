# ==============================================================================
# epistasis.R
# Within-block and between-block haplotype epistasis detection.
#
# Three exported functions:
#
#   scan_block_epistasis()
#     Pairwise SNP x SNP interaction scan within significant haplotype blocks.
#     Tests H0: aa_ij = 0 in y = mu + a_i*x_i + a_j*x_j + aa_ij*(x_i*x_j) + e
#     on GRM-corrected REML residuals. Restricted to blocks flagged as
#     significant by test_block_haplotypes() to avoid genome-wide explosion.
#
#   scan_block_by_block_epistasis()
#     Haplotype allele x haplotype allele interaction scan between blocks.
#     Tests interactions between dosage columns from different blocks using
#     the same REML residuals. Restricted to: significant alleles x all blocks,
#     making the scan O(n_sig_alleles x n_total_alleles) -- tractable even at
#     WGS scale.
#
#   fine_map_epistasis_block()
#     Fine-maps interacting SNP pairs within a single block. Dispatches to:
#       "pairwise" -- exhaustive C(p,2) scan (recommended for p <= 200)
#       "lasso"    -- LASSO with pairwise interaction terms (p > 200)
#     Returns ranked interacting pairs with effect, SE, p-value.
#
# Statistical model (all functions)
# ----------------------------------
# All tests operate on REML residuals from the same null model as
# test_block_haplotypes():
#   y_resid = y - X_null*beta_hat - u_hat
# where u ~ MVN(0, sigma_g^2 * G_hap) is the polygenic random effect.
# This ensures all epistasis tests are population-structure-corrected and
# consistent with the marginal association results.
#
# Multiple testing correction
# ---------------------------
# scan_block_epistasis():     Bonferroni AND simpleM Sidak within each block.
#                             Meff estimated from eigenspectrum of the pairwise
#                             interaction column correlation matrix. Both
#                             p_bonf and p_simplem_sidak are always returned;
#                             significant flag uses sig_metric parameter.
# scan_block_by_block():      Bonferroni over n_sig x n_total_alleles tests.
# fine_map_epistasis_block(): Bonferroni over C(p,2) pairs (pairwise method).
# ==============================================================================


# ==============================================================================
# Internal helper: fit REML null model and return de-regressed residuals.
# Shared by all three epistasis functions.
# ==============================================================================
.fit_null_reml <- function(y, G, n_pcs = 0L, verbose = FALSE) {
  if (!requireNamespace("rrBLUP", quietly = TRUE))
    stop("rrBLUP is required.", call. = FALSE)

  n_pcs <- as.integer(n_pcs)
  X_null <- if (n_pcs > 0L) {
    eig <- tryCatch(eigen(G, symmetric = TRUE), error = function(e) NULL)
    if (is.null(eig)) {
      matrix(1, nrow = length(y), ncol = 1L)
    } else {
      cbind(1, eig$vectors[, seq_len(min(n_pcs, ncol(eig$vectors))),
                           drop = FALSE])
    }
  } else {
    matrix(1, nrow = length(y), ncol = 1L)
  }

  fit <- tryCatch(
    rrBLUP::mixed.solve(y = y, X = X_null, K = G, method = "REML"),
    error = function(e) NULL
  )
  if (is.null(fit)) {
    # Fallback to OLS residuals
    lm_fit <- stats::lm.fit(x = X_null, y = y)
    return(stats::residuals(lm_fit))
  }
  y_fixed <- as.numeric(X_null %*% fit$beta)
  u       <- fit$u; u[is.na(u)] <- 0
  y - y_fixed - u
}


# ==============================================================================
# Internal helper: vectorised pairwise interaction Wald test.
# Tests H0: aa_ij = 0 for all pairs (i,j) in columns of X.
# Uses REML residuals y_resid.
# Returns data.frame: col_i, col_j, aa_effect, SE, t_stat, p_wald
# ==============================================================================
.pairwise_interaction_scan <- function(X, y_resid, col_names = NULL) {
  n    <- length(y_resid)
  p    <- ncol(X)
  if (p < 2L) return(data.frame())

  yb <- mean(y_resid)
  ss_tot <- sum((y_resid - yb)^2)

  rows <- vector("list", p * (p - 1L) / 2L)
  k    <- 0L

  for (i in seq_len(p - 1L)) {
    xi <- X[, i]
    xi_c <- xi - mean(xi)
    for (j in seq(i + 1L, p)) {
      xj  <- X[, j]
      # Interaction column: (x_i - mean(x_i)) * (x_j - mean(x_j))
      # Centering reduces collinearity with main effects
      xij <- xi_c * (xj - mean(xj))
      xij_var <- sum(xij^2) - sum(xij)^2 / n
      if (!is.finite(xij_var) || xij_var < 1e-10) next

      # Partial regression of xij on y_resid after projecting out xi, xj
      # For speed: use simple OLS on [xi, xj, xi*xj] matrix
      Xm   <- cbind(1, xi, xj, xij)
      fit_m <- tryCatch(
        stats::lm.fit(x = Xm, y = y_resid),
        error = function(e) NULL
      )
      if (is.null(fit_m)) next
      coefs <- stats::coef(fit_m)
      if (length(coefs) < 4L || is.na(coefs[4L])) next

      rss   <- sum(stats::residuals(fit_m)^2)
      df_r  <- n - 4L
      if (df_r < 1L) next
      sigma2 <- rss / df_r

      # SE of interaction coefficient via XtX inverse diagonal
      xtx_inv <- tryCatch(solve(crossprod(Xm)), error = function(e) NULL)
      if (is.null(xtx_inv)) next
      se_int <- sqrt(sigma2 * xtx_inv[4L, 4L])
      if (!is.finite(se_int) || se_int <= 0) next

      aa   <- coefs[4L]
      t_aa <- aa / se_int
      p_aa <- 2 * stats::pt(-abs(t_aa), df = df_r)

      k <- k + 1L
      rows[[k]] <- data.frame(
        col_i     = if (!is.null(col_names)) col_names[i] else as.character(i),
        col_j     = if (!is.null(col_names)) col_names[j] else as.character(j),
        aa_effect = round(aa,    6),
        SE        = round(se_int, 6),
        t_stat    = round(t_aa,  4),
        p_wald    = p_aa,
        stringsAsFactors = FALSE
      )
    }
  }
  if (k == 0L) return(data.frame())
  do.call(rbind, rows[seq_len(k)])
}


# ==============================================================================
# 1. scan_block_epistasis()
# ==============================================================================

#' Within-Block Pairwise SNP Epistasis Scan
#'
#' @description
#' Tests pairwise SNP x SNP interactions within haplotype blocks that were
#' identified as significant by \code{\link{test_block_haplotypes}}. For each
#' significant block, all C(p, 2) SNP pairs are tested for interaction using
#' the model:
#'
#' \deqn{y = \mu + a_i x_i + a_j x_j + aa_{ij}(x_i x_j) + \varepsilon}
#'
#' on GRM-corrected REML residuals (same null model as
#' \code{\link{test_block_haplotypes}}), ensuring all tests are
#' population-structure-corrected.
#'
#' Restricting to significant blocks avoids the genome-wide explosion of
#' pairwise tests: for 15 significant blocks with ~200 SNPs each, the total
#' number of tests is ~300,000, compared to ~4.4 billion for an unrestricted
#' genome-wide scan.
#'
#' @param assoc Output of \code{\link{test_block_haplotypes}}. Used to
#'   identify significant blocks and obtain pre-computed GRM residuals.
#' @param geno_matrix Numeric matrix (individuals x SNPs) or
#'   \code{LDxBlocks_backend}. The imputed, MAF-filtered genotype matrix
#'   from Job 1 (\code{res$geno_matrix}).
#' @param snp_info Data frame with columns \code{SNP}, \code{CHR}, \code{POS}.
#' @param blocks LD block table from \code{run_Big_LD_all_chr}.
#' @param blues Pre-adjusted phenotype means (same format as
#'   \code{test_block_haplotypes}). Named numeric vector or named list.
#' @param haplotypes Named list from \code{extract_haplotypes}. Used to
#'   re-build the GRM for REML residual computation.
#' @param trait Character. Which trait to use for residual computation when
#'   \code{blues} is a named list. Default \code{NULL} uses the first trait.
#' @param sig_blocks Character vector. Block IDs to scan. \code{NULL} (default)
#'   uses all blocks with \code{significant_omnibus = TRUE} in
#'   \code{assoc$block_tests}.
#' @param min_freq Numeric. Minimum MAF for SNPs within the block.
#'   Default \code{0.05}.
#' @param max_snps_per_block Integer. Maximum SNPs per block before switching
#'   to random subsampling of pairs. Default \code{300L} (C(300,2)=44,850
#'   pairs per block). Set \code{NULL} to always use all SNPs.
#' @param sig_threshold Numeric. Significance threshold. Default \code{0.05}.
#' @param sig_metric Character. Which correction drives the primary
#'   \code{significant} flag. One of:
#'   \itemize{
#'     \item \code{"p_simplem_sidak"} (default, recommended) -- simpleM Sidak.
#'     \item \code{"p_simplem"} -- simpleM Bonferroni.
#'     \item \code{"p_bonf"} -- plain Bonferroni (p x n_pairs).
#'   }
#'   All three p-value columns are always present regardless of this choice.
#' @param meff_percent_cut Numeric. Variance cutoff for simpleM Meff.
#'   Default \code{0.995}.
#' @param id_col Character. ID column when blues is a data frame.
#' @param blue_col Character. Phenotype column when blues is a data frame.
#' @param verbose Logical. Default \code{TRUE}.
#'
#' @return A named list of class \code{LDxBlocks_epistasis}:
#' \describe{
#'   \item{\code{results}}{Data frame. One row per tested SNP pair per block
#'     per trait. Always present columns: \code{block_id}, \code{CHR},
#'     \code{start_bp}, \code{end_bp}, \code{trait}, \code{SNP_i},
#'     \code{SNP_j}, \code{POS_i}, \code{POS_j}, \code{dist_bp},
#'     \code{aa_effect}, \code{SE}, \code{t_stat}, \code{p_wald},
#'     \code{Meff} (simpleM effective test count from interaction eigenspectrum),
#'     \code{p_bonf} (Bonferroni: p x n_pairs),
#'     \code{p_simplem} (simpleM Bonferroni: p x Meff),
#'     \code{p_simplem_sidak} (simpleM Sidak: 1-(1-p)^Meff),
#'     \code{significant} (driven by \code{sig_metric}),
#'     \code{significant_bonf}, \code{significant_simplem},
#'     \code{significant_simplem_sidak}.}
#'   \item{\code{scan_summary}}{Data frame. One row per block: number of pairs
#'     tested, number significant, minimum p-value.}
#'   \item{\code{n_blocks_scanned}}{Integer.}
#'   \item{\code{n_pairs_total}}{Integer. Total pairwise tests performed.}
#' }
#'
#' @references
#' Cordell HJ (2009). Detecting gene-gene interactions that underlie human
#' diseases. \emph{Nature Reviews Genetics} \strong{10}:392-404.
#'
#' @seealso \code{\link{test_block_haplotypes}},
#'   \code{\link{scan_block_by_block_epistasis}},
#'   \code{\link{fine_map_epistasis_block}}
#' @importFrom utils combn
#' @export
scan_block_epistasis <- function(
    assoc,
    geno_matrix,
    snp_info,
    blocks,
    blues,
    haplotypes,
    trait            = NULL,
    sig_blocks       = NULL,
    min_freq         = 0.05,
    max_snps_per_block = 300L,
    sig_threshold    = 0.05,
    sig_metric       = c("p_simplem_sidak", "p_simplem", "p_bonf", "p_fdr"),
    meff_percent_cut = 0.995,
    id_col           = "id",
    blue_col         = "blue",
    verbose          = TRUE
) {
  if (!inherits(assoc, "LDxBlocks_haplotype_assoc"))
    stop("assoc must be output of test_block_haplotypes().", call. = FALSE)

  sig_metric <- match.arg(sig_metric)
  .log <- function(...) if (verbose) message(sprintf("[epi_block] %s", paste0(...)))

  # -- Resolve trait and blues -------------------------------------------------
  blues_list <- .parse_blues_assoc(blues, id_col, blue_col, NULL)
  tr <- if (is.null(trait)) names(blues_list)[1L] else trait
  if (!tr %in% names(blues_list))
    stop("trait '", tr, "' not found in blues.", call. = FALSE)
  pheno <- blues_list[[tr]]

  # -- Significant blocks to scan ----------------------------------------------
  if (is.null(sig_blocks)) {
    bt <- assoc$block_tests
    if (!"significant_omnibus" %in% names(bt))
      stop("assoc$block_tests missing 'significant_omnibus'. ",
           "Run test_block_haplotypes() with sig_metric set.", call. = FALSE)
    sig_blocks <- unique(bt$block_id[bt$significant_omnibus &
                                       !is.na(bt$significant_omnibus)])
  }
  if (!length(sig_blocks)) {
    .log("No significant blocks to scan.")
    return(structure(
      list(results = data.frame(), scan_summary = data.frame(),
           n_blocks_scanned = 0L, n_pairs_total = 0L),
      class = "LDxBlocks_epistasis"))
  }
  .log(length(sig_blocks), " significant blocks to scan for within-block epistasis")

  # -- Build GRM and REML residuals (shared across all blocks) ----------------
  .log("Building haplotype GRM for REML null model ...")
  hap_mat <- build_haplotype_feature_matrix(
    haplotypes, min_freq = min_freq, encoding = "additive_012")$matrix
  G <- tryCatch(
    compute_haplotype_grm(hap_mat),
    error = function(e) {
      diag_G <- diag(nrow(hap_mat))
      rownames(diag_G) <- colnames(diag_G) <- rownames(hap_mat)
      diag_G
    }
  )

  common <- intersect(names(pheno), rownames(G))
  common <- common[!is.na(pheno[common]) & is.finite(pheno[common])]
  if (length(common) < 10L)
    stop("Fewer than 10 common individuals.", call. = FALSE)

  y      <- pheno[common]
  G_sub  <- G[common, common, drop = FALSE]
  y_resid <- .fit_null_reml(y, G_sub, n_pcs = assoc$n_pcs_used,
                            verbose = verbose)
  names(y_resid) <- common
  .log("REML residuals computed (n=", length(y_resid), ")")

  # -- Normalise inputs --------------------------------------------------------
  snp_info$CHR  <- .norm_chr_hap(snp_info$CHR)
  blocks$CHR    <- .norm_chr_hap(blocks$CHR)
  if (!"start_bp" %in% names(blocks)) {
    blocks$start_bp <- blocks$start.bp
    blocks$end_bp   <- blocks$end.bp
  }
  if (!"block_id" %in% names(blocks))
    blocks$block_id <- paste0("block_", blocks$CHR, "_",
                              blocks$start_bp, "_", blocks$end_bp)

  # -- Per-block pairwise scan -------------------------------------------------
  result_rows  <- list()
  summary_rows <- list()
  n_pairs_total <- 0L

  # Handle both matrix and backend inputs
  is_backend <- inherits(geno_matrix, "LDxBlocks_backend")

  for (bn in sig_blocks) {
    brow <- blocks[blocks$block_id == bn, , drop = FALSE]
    if (!nrow(brow)) next
    chr_b <- brow$CHR[1]; start_b <- brow$start_bp[1]; end_b <- brow$end_bp[1]

    # SNPs within this block
    si_b <- snp_info[snp_info$CHR == chr_b &
                       snp_info$POS >= start_b &
                       snp_info$POS <= end_b, , drop = FALSE]
    if (nrow(si_b) < 2L) next

    # Extract genotype submatrix for block SNPs
    if (is_backend) {
      snp_idx <- which(snp_info$CHR == chr_b &
                         snp_info$POS >= start_b &
                         snp_info$POS <= end_b)
      G_blk_raw <- read_chunk(geno_matrix, snp_idx)
    } else {
      snp_idx <- which(colnames(geno_matrix) %in% si_b$SNP)
      if (!length(snp_idx))
        snp_idx <- which(snp_info$CHR == chr_b &
                           snp_info$POS >= start_b &
                           snp_info$POS <= end_b)
      G_blk_raw <- geno_matrix[, snp_idx, drop = FALSE]
    }

    # Align to individuals with residuals
    G_blk <- G_blk_raw[rownames(G_blk_raw) %in% common, , drop = FALSE]
    G_blk <- G_blk[common[common %in% rownames(G_blk)], , drop = FALSE]
    if (nrow(G_blk) < 10L) next

    y_r_bn <- y_resid[rownames(G_blk)]

    # MAF filter -- keep si_b aligned with G_blk columns
    maf_b    <- apply(G_blk, 2L, function(x) {
      x <- x[!is.na(x)]; af <- mean(x, na.rm = TRUE) / 2
      min(af, 1 - af)
    })
    maf_keep <- maf_b >= min_freq
    G_blk    <- G_blk[, maf_keep, drop = FALSE]
    if (ncol(G_blk) < 2L) next
    # Re-align si_b to surviving columns by column name
    if (!is.null(colnames(G_blk)) && all(colnames(G_blk) %in% si_b$SNP))
      si_b <- si_b[match(colnames(G_blk), si_b$SNP), , drop = FALSE]

    # Mean-impute NAs
    for (j in seq_len(ncol(G_blk))) {
      na_j <- is.na(G_blk[, j])
      if (any(na_j)) G_blk[na_j, j] <- mean(G_blk[, j], na.rm = TRUE)
    }

    # Remove monomorphic columns -- keep si_b aligned
    col_var  <- apply(G_blk, 2L, stats::var)
    var_keep <- col_var > 1e-10
    G_blk    <- G_blk[, var_keep, drop = FALSE]
    if (ncol(G_blk) < 2L) next
    if (!is.null(colnames(G_blk)) && all(colnames(G_blk) %in% si_b$SNP))
      si_b <- si_b[match(colnames(G_blk), si_b$SNP), , drop = FALSE]

    p_blk     <- ncol(G_blk)
    snp_names <- colnames(G_blk)
    if (is.null(snp_names)) snp_names <- si_b$SNP[seq_len(p_blk)]
    # pos_map is derived from the now-aligned si_b
    pos_map_blk <- if (!is.null(snp_names) && all(snp_names %in% si_b$SNP))
      stats::setNames(si_b$POS, si_b$SNP) else
        stats::setNames(si_b$POS[seq_len(p_blk)], snp_names)

    # Subsample SNPs if block is very large; keep pos_map_blk aligned
    if (!is.null(max_snps_per_block) && p_blk > max_snps_per_block) {
      .log("  ", bn, ": ", p_blk, " SNPs > max (", max_snps_per_block,
           ") -- subsampling")
      sel       <- sort(sample(p_blk, max_snps_per_block))
      G_blk     <- G_blk[, sel, drop = FALSE]
      snp_names <- snp_names[sel]
      p_blk     <- max_snps_per_block
      # Re-derive pos_map_blk for subsampled SNPs
      pos_map_blk <- pos_map_blk[snp_names]
    }

    n_pairs <- p_blk * (p_blk - 1L) / 2L
    n_pairs_total <- n_pairs_total + n_pairs
    .log("  Scanning ", bn, ": ", p_blk, " SNPs, ", n_pairs, " pairs")

    # Run pairwise interaction scan
    scan_df <- .pairwise_interaction_scan(G_blk, y_r_bn, col_names = snp_names)
    if (!nrow(scan_df)) {
      summary_rows[[bn]] <- data.frame(
        block_id = bn, CHR = chr_b, trait = tr,
        n_snps = p_blk, n_pairs = n_pairs,
        n_significant = 0L, min_p = NA_real_,
        stringsAsFactors = FALSE)
      next
    }

    # -- Multiple testing correction within block -----------------------------
    # Bonferroni (always computed)
    scan_df$p_bonf <- pmin(scan_df$p_wald * n_pairs, 1)
    scan_df$p_fdr  <- stats::p.adjust(scan_df$p_wald, method = "BH")

    # simpleM on the pairwise interaction columns (Meff from eigenspectrum)
    # Build the interaction column matrix for Meff estimation
    pair_cols <- utils::combn(seq_len(ncol(G_blk)), 2L)
    n_pairs_c <- ncol(pair_cols)
    X_pair    <- matrix(NA_real_, nrow = nrow(G_blk), ncol = n_pairs_c)
    for (kk in seq_len(n_pairs_c)) {
      ii <- pair_cols[1L, kk]; jj <- pair_cols[2L, kk]
      xi_c <- G_blk[, ii] - mean(G_blk[, ii])
      xj_c <- G_blk[, jj] - mean(G_blk[, jj])
      X_pair[, kk] <- xi_c * xj_c
    }
    meff_blk <- tryCatch(
      .compute_simplem_meff_matrix(
        X_pair, percent_cut = meff_percent_cut, max_cols = 1000L)$meff,
      error = function(e) n_pairs
    )
    if (!is.finite(meff_blk) || meff_blk < 1L) meff_blk <- n_pairs

    scan_df$Meff                    <- meff_blk
    scan_df$p_simplem               <- .adjust_p_simplem_bonf(scan_df$p_wald, meff_blk)
    scan_df$p_simplem_sidak         <- .adjust_p_simplem_sidak(scan_df$p_wald, meff_blk)
    scan_df$significant_bonf        <- scan_df$p_bonf        < sig_threshold
    scan_df$significant_fdr         <- scan_df$p_fdr         < sig_threshold
    scan_df$significant_simplem     <- scan_df$p_simplem      < sig_threshold
    scan_df$significant_simplem_sidak <- scan_df$p_simplem_sidak < sig_threshold
    # Primary significant flag: driven by sig_metric
    scan_df$significant <- switch(sig_metric,
                                  p_simplem_sidak = !is.na(scan_df$significant_simplem_sidak) & scan_df$significant_simplem_sidak,
                                  p_simplem       = !is.na(scan_df$significant_simplem)       & scan_df$significant_simplem,
                                  p_bonf          = !is.na(scan_df$significant_bonf)          & scan_df$significant_bonf,
                                  p_fdr           = !is.na(scan_df$significant_fdr)           & scan_df$significant_fdr
    )

    # Add metadata
    pos_map <- stats::setNames(si_b$POS, si_b$SNP)
    scan_df$block_id <- bn
    scan_df$CHR      <- chr_b
    scan_df$start_bp <- start_b
    scan_df$end_bp   <- end_b
    scan_df$trait    <- tr
    scan_df$POS_i    <- pos_map_blk[scan_df$col_i]
    scan_df$POS_j    <- pos_map_blk[scan_df$col_j]
    scan_df$dist_bp  <- abs(scan_df$POS_j - scan_df$POS_i)
    names(scan_df)[names(scan_df) == "col_i"] <- "SNP_i"
    names(scan_df)[names(scan_df) == "col_j"] <- "SNP_j"

    # Reorder columns
    col_order <- c("block_id","CHR","start_bp","end_bp","trait",
                   "SNP_i","SNP_j","POS_i","POS_j","dist_bp",
                   "aa_effect","SE","t_stat","p_wald",
                   "Meff","p_bonf","p_fdr","p_simplem","p_simplem_sidak",
                   "significant","significant_bonf","significant_fdr",
                   "significant_simplem","significant_simplem_sidak")
    scan_df <- scan_df[, intersect(col_order, names(scan_df)), drop = FALSE]
    scan_df <- scan_df[order(scan_df$p_wald), ]

    result_rows[[bn]]  <- scan_df
    summary_rows[[bn]] <- data.frame(
      block_id = bn, CHR = chr_b, trait = tr,
      n_snps = p_blk, n_pairs = n_pairs,
      n_significant = sum(scan_df$significant, na.rm = TRUE),
      min_p = min(scan_df$p_wald, na.rm = TRUE),
      stringsAsFactors = FALSE)
  }

  results_df <- if (length(result_rows))
    do.call(rbind, result_rows) else data.frame()
  summary_df <- if (length(summary_rows))
    do.call(rbind, summary_rows) else data.frame()
  rownames(results_df) <- NULL
  rownames(summary_df) <- NULL

  .log("Done. Blocks scanned: ", length(sig_blocks),
       " | Total pairs: ", format(n_pairs_total, big.mark = ","),
       " | Significant (", sig_metric, "): ",
       sum(results_df$significant, na.rm = TRUE))

  structure(
    list(results          = results_df,
         scan_summary     = summary_df,
         n_blocks_scanned = length(sig_blocks),
         n_pairs_total    = n_pairs_total,
         sig_metric        = sig_metric,
         sig_threshold     = sig_threshold),
    class = "LDxBlocks_epistasis"
  )
}


# ==============================================================================
# 2. scan_block_by_block_epistasis()
# ==============================================================================

#' Between-Block Haplotype Allele Epistasis Scan
#'
#' @description
#' Tests interactions between haplotype allele dosage columns from different
#' LD blocks. For each significant allele (from
#' \code{\link{test_block_haplotypes}}), tests its interaction with every
#' haplotype allele at all other blocks:
#'
#' \deqn{y = \mu + \alpha_i x_i + \alpha_j x_j + \gamma_{ij}(x_i x_j) + \varepsilon}
#'
#' Restricting to significant alleles x all blocks gives
#' O(n_sig x n_total_alleles) tests -- tractable at WGS scale.
#' For your panel: ~25 significant alleles x 17,943 alleles = ~450,000 tests.
#'
#' This is a \strong{trans-haplotype epistasis scan}: it tests whether the
#' effect of a resistance haplotype at one block depends on the haplotype
#' background at another block -- a form of genetic background dependence
#' that single-SNP and even single-block analyses cannot detect.
#'
#' @param assoc Output of \code{\link{test_block_haplotypes}}.
#' @param haplotypes Named list from \code{extract_haplotypes}.
#' @param blues Pre-adjusted phenotype means.
#' @param blocks LD block table.
#' @param trait Character. Trait to use. Default \code{NULL} uses first trait.
#' @param sig_alleles Data frame with columns \code{block_id} and \code{allele}
#'   identifying the query alleles. \code{NULL} (default) uses all
#'   \code{significant = TRUE} rows from \code{assoc$allele_tests}.
#' @param min_freq Numeric. Minimum allele frequency. Default \code{0.05}.
#' @param top_n Integer or \code{NULL}. Alleles per block in feature matrix.
#' @param sig_threshold Numeric. Significance threshold applied to the
#'   p-value chosen by \code{sig_metric}. Default \code{0.05}.
#' @param sig_metric Character. Which p-value drives the \code{significant}
#'   flag. One of \code{"p_simplem_sidak"} (default, recommended),
#'   \code{"p_simplem"}, \code{"p_bonf"}, or \code{"p_fdr"}.
#'   All four p-value columns are always present in the output regardless
#'   of this choice.
#' @param meff_percent_cut Numeric in (0, 1). Variance threshold for simpleM
#'   eigendecomposition when estimating \eqn{M_{\mathrm{eff}}}.
#'   Default \code{0.995} (99.5\%), following Gao et al. (2008).
#' @param id_col Character. ID column when blues is a data frame.
#' @param blue_col Character. Phenotype column when blues is a data frame.
#' @param verbose Logical. Default \code{TRUE}.
#'
#' @return A named list of class \code{LDxBlocks_block_epistasis}:
#' \describe{
#'   \item{\code{results}}{Data frame sorted by p_wald. Columns:
#'     \code{block_i}, \code{allele_i}, \code{block_j}, \code{allele_j},
#'     \code{CHR_i}, \code{CHR_j}, \code{same_chr},
#'     \code{aa_effect}, \code{SE}, \code{t_stat},
#'     \code{p_wald}, \code{p_bonf}, \code{significant}.}
#'   \item{\code{n_tests}}{Integer. Total interaction tests performed.}
#'   \item{\code{n_sig_alleles}}{Integer. Number of query alleles.}
#' }
#'
#' @references
#' Mackay TFC (2014). Epistasis and quantitative traits: using model
#' organisms to study gene-gene interactions. \emph{Nature Reviews Genetics}
#' \strong{15}:22-33.
#'
#' @seealso \code{\link{scan_block_epistasis}},
#'   \code{\link{test_block_haplotypes}},
#'   \code{\link{fine_map_epistasis_block}}
#' @export
scan_block_by_block_epistasis <- function(
    assoc,
    haplotypes,
    blues,
    blocks,
    trait         = NULL,
    sig_alleles   = NULL,
    min_freq      = 0.05,
    top_n         = NULL,
    sig_threshold = 0.05,
    sig_metric    = c("p_simplem_sidak", "p_simplem", "p_bonf", "p_fdr"),
    meff_percent_cut = 0.995,
    id_col        = "id",
    blue_col      = "blue",
    verbose       = TRUE
) {
  if (!inherits(assoc, "LDxBlocks_haplotype_assoc"))
    stop("assoc must be output of test_block_haplotypes().", call. = FALSE)

  sig_metric <- match.arg(sig_metric)
  .log <- function(...) if (verbose) message(sprintf("[epi_block_x_block] %s", paste0(...)))

  # -- Resolve trait and blues -------------------------------------------------
  blues_list <- .parse_blues_assoc(blues, id_col, blue_col, NULL)
  tr <- if (is.null(trait)) names(blues_list)[1L] else trait
  pheno <- blues_list[[tr]]

  # -- Significant query alleles -----------------------------------------------
  if (is.null(sig_alleles)) {
    at <- assoc$allele_tests
    at_sig <- at[!is.na(at$significant) & at$significant &
                   at$trait == tr, , drop = FALSE]
    if (!nrow(at_sig)) {
      .log("No significant alleles for trait '", tr, "'.")
      return(structure(
        list(results = data.frame(), n_tests = 0L, n_sig_alleles = 0L),
        class = "LDxBlocks_block_epistasis"))
    }
    sig_alleles <- at_sig[, c("block_id", "allele"), drop = FALSE]
  }
  n_sig <- nrow(sig_alleles)
  .log(n_sig, " query alleles x all blocks (trans-haplotype epistasis scan)")

  # -- Build full haplotype feature matrix ------------------------------------
  .log("Building haplotype feature matrix ...")
  feat_obj <- build_haplotype_feature_matrix(
    haplotypes, top_n = top_n, min_freq = min_freq, encoding = "additive_012")
  hap_mat  <- feat_obj$matrix
  col_blocks <- sub("_hap[0-9]+$", "", colnames(hap_mat))

  # -- Build GRM and REML residuals -------------------------------------------
  .log("Computing GRM and REML residuals ...")
  G <- tryCatch(
    compute_haplotype_grm(hap_mat),
    error = function(e) {
      diag_G <- diag(nrow(hap_mat))
      rownames(diag_G) <- colnames(diag_G) <- rownames(hap_mat)
      diag_G
    }
  )
  common <- intersect(names(pheno), rownames(G))
  common <- common[!is.na(pheno[common]) & is.finite(pheno[common])]
  y       <- pheno[common]
  G_sub   <- G[common, common, drop = FALSE]
  y_resid <- .fit_null_reml(y, G_sub, n_pcs = assoc$n_pcs_used)
  names(y_resid) <- common

  hap_sub <- hap_mat[common, , drop = FALSE]
  # Mean-impute NAs
  for (j in seq_len(ncol(hap_sub))) {
    na_j <- is.na(hap_sub[, j])
    if (any(na_j)) hap_sub[na_j, j] <- mean(hap_sub[, j], na.rm = TRUE)
  }
  # Remove constant columns
  col_var    <- apply(hap_sub, 2L, stats::var)
  valid_cols <- col_var > 1e-10
  hap_sub    <- hap_sub[, valid_cols, drop = FALSE]

  if (ncol(hap_sub) < 2L) {
    .log("Fewer than 2 valid haplotype allele columns after filtering -- cannot scan.")
    return(structure(
      list(results = data.frame(), n_tests = 0L, n_sig_alleles = n_sig),
      class = "LDxBlocks_block_epistasis"))
  }

  col_names  <- colnames(hap_sub)
  col_blks   <- col_blocks[valid_cols]

  # CHR lookup for block metadata
  if (!"start_bp" %in% names(blocks)) {
    blocks$start_bp <- blocks$start.bp
    blocks$end_bp   <- blocks$end.bp
  }
  if (!"block_id" %in% names(blocks))
    blocks$block_id <- paste0("block_", .norm_chr_hap(blocks$CHR), "_",
                              blocks$start_bp, "_", blocks$end_bp)
  blocks$CHR <- .norm_chr_hap(blocks$CHR)
  block_chr_map <- stats::setNames(blocks$CHR, blocks$block_id)

  n_total  <- ncol(hap_sub)
  n_tests  <- n_sig * (n_total - 1L)
  .log("Total interaction tests: ", format(n_tests, big.mark = ","),
       " (", n_sig, " x ", n_total - 1L, ")")

  yb    <- mean(y_resid)
  ss_tot <- sum((y_resid - yb)^2)

  result_rows <- vector("list", n_sig * 50L)   # pre-allocate estimate
  k <- 0L

  for (qi in seq_len(n_sig)) {
    q_block  <- sig_alleles$block_id[qi]
    q_allele <- sig_alleles$allele[qi]

    # Find the column index of this query allele.
    # Strategy: try exact match "block_id_allele" first, then escape the
    # allele string for regex and search for it as a suffix.
    # Never fall back to a different allele from the same block.
    q_col_name <- paste0(q_block, "_", q_allele)
    q_idx <- which(col_names == q_col_name)
    if (!length(q_idx)) {
      # Allele label may have been stored differently; try regex-escaped suffix
      q_allele_esc <- gsub("([.|()*+?\\^${}\\[\\]])", "\\\1", q_allele)
      q_idx <- grep(paste0("^", q_block, "_", q_allele_esc, "$"), col_names)
    }
    if (!length(q_idx)) next   # allele not in feature matrix -- skip, do NOT substitute

    xi   <- hap_sub[, q_idx]
    xi_c <- xi - mean(xi)

    # Test against all other blocks (skip same block)
    other_idx <- which(col_blks != q_block)
    if (!length(other_idx)) next

    X_other <- hap_sub[, other_idx, drop = FALSE]
    n_ind   <- length(y_resid)

    # Vectorised interaction scan: for each column j in X_other
    # Interaction column: xi_c * (xj - mean(xj))
    # Partial regression coefficient via OLS on [1, xi, xj, xi*xj]
    # Use the Frisch-Waugh-Lovell theorem shortcut:
    # Project xi and xj out of y_resid and the interaction column,
    # then regress residuals on the interaction residuals.
    # Full shortcut for speed:
    # beta_int = Cov(y_resid*, xi_c * xj_c) / Var(xi_c * xj_c)
    # where y_resid* = resid(y_resid ~ xi + xj),
    #       xi_c*xj_c = centred interaction column

    # For each j: build [1, xi, xj, xi*xj], OLS, extract interaction coef
    # Do this as a vectorised loop (no pure R per-pair loop overhead for
    # the inner arithmetic; only the OLS solve is per-pair)
    x_means <- colMeans(X_other)
    X_c     <- sweep(X_other, 2L, x_means, "-")   # centred partner columns
    X_int   <- xi_c * X_c                         # n x n_other interaction matrix

    # Pre-compute projection of y_resid onto [1, xi]
    X_main <- cbind(1, xi)
    XtX_main <- crossprod(X_main)
    XtX_inv_main <- tryCatch(solve(XtX_main), error = function(e) NULL)
    if (is.null(XtX_inv_main)) next
    beta_main <- XtX_inv_main %*% crossprod(X_main, y_resid)
    y_star    <- y_resid - as.numeric(X_main %*% beta_main)

    for (jj in seq_len(ncol(X_other))) {
      xj    <- X_other[, jj]
      xj_c  <- X_c[, jj]
      x_int <- X_int[, jj]

      # Full [1, xi, xj, xi*xj] OLS
      Xm <- cbind(1, xi, xj, x_int)
      fit_m <- tryCatch(stats::lm.fit(x = Xm, y = y_resid),
                        error = function(e) NULL)
      if (is.null(fit_m)) next
      coefs <- stats::coef(fit_m)
      if (length(coefs) < 4L || is.na(coefs[4L])) next
      rss   <- sum(stats::residuals(fit_m)^2)
      df_r  <- n_ind - 4L
      if (df_r < 1L) next
      sigma2 <- rss / df_r
      xtx_inv <- tryCatch(solve(crossprod(Xm)), error = function(e) NULL)
      if (is.null(xtx_inv)) next
      se_int <- sqrt(sigma2 * xtx_inv[4L, 4L])
      if (!is.finite(se_int) || se_int <= 0) next
      aa   <- coefs[4L]
      t_aa <- aa / se_int
      p_aa <- 2 * stats::pt(-abs(t_aa), df = df_r)

      k <- k + 1L
      if (k > length(result_rows))
        result_rows <- c(result_rows, vector("list", 1000L))

      j_block <- col_blks[other_idx[jj]]
      result_rows[[k]] <- data.frame(
        block_i    = q_block,
        allele_i   = q_allele,
        block_j    = j_block,
        allele_j   = sub(paste0("^", j_block, "_"), "", col_names[other_idx[jj]]),
        CHR_i      = unname(block_chr_map[q_block]),
        CHR_j      = unname(block_chr_map[j_block]),
        same_chr   = identical(unname(block_chr_map[q_block]),
                               unname(block_chr_map[j_block])),
        aa_effect  = round(aa,    6),
        SE         = round(se_int, 6),
        t_stat     = round(t_aa,  4),
        p_wald     = p_aa,
        stringsAsFactors = FALSE
      )
    }
  }

  results_df <- if (k > 0L) {
    df <- do.call(rbind, result_rows[seq_len(k)])
    df$p_bonf  <- pmin(df$p_wald * n_tests, 1)
    df$p_fdr   <- stats::p.adjust(df$p_wald, method = "BH")
    # simpleM Meff from eigenspectrum of all tested interaction columns
    # Use the haplotype dosage matrix already in scope (hap_sub)
    meff_bb <- tryCatch(
      .compute_simplem_meff_matrix(
        hap_sub, percent_cut = meff_percent_cut, max_cols = 1000L)$meff,
      error = function(e) n_tests
    )
    if (!is.finite(meff_bb) || meff_bb < 1L) meff_bb <- n_tests
    df$Meff            <- meff_bb
    df$p_simplem       <- .adjust_p_simplem_bonf(df$p_wald, meff_bb)
    df$p_simplem_sidak <- .adjust_p_simplem_sidak(df$p_wald, meff_bb)
    df$significant <- switch(sig_metric,
                             p_simplem_sidak = !is.na(df$p_simplem_sidak) & df$p_simplem_sidak < sig_threshold,
                             p_simplem       = !is.na(df$p_simplem)       & df$p_simplem       < sig_threshold,
                             p_bonf          = !is.na(df$p_bonf)          & df$p_bonf          < sig_threshold,
                             p_fdr           = !is.na(df$p_fdr)           & df$p_fdr           < sig_threshold
    )
    df[order(df$p_wald), ]
  } else data.frame()
  rownames(results_df) <- NULL

  n_sig_out <- sum(results_df$significant, na.rm = TRUE)
  .log("Done. Tests: ", format(n_tests, big.mark=","),
       " | Significant (", sig_metric, "): ", n_sig_out)

  structure(
    list(results       = results_df,
         n_tests       = n_tests,
         n_sig_alleles = n_sig,
         sig_metric    = sig_metric,
         sig_threshold = sig_threshold),
    class = "LDxBlocks_block_epistasis"
  )
}


# ==============================================================================
# 3. fine_map_epistasis_block()
# ==============================================================================

#' Fine-Map Epistatic SNP Pairs Within a Single Block
#'
#' @description
#' Exhaustively or adaptively identifies the specific SNP pairs within a
#' single LD block that drive the block-level epistatic signal. Three methods
#' are available:
#'
#' \describe{
#'   \item{\code{"pairwise"}}{Exhaustive C(p,2) scan. Tests all pairs using
#'     the same model as \code{\link{scan_block_epistasis}}. Recommended for
#'     blocks with p <= 200 SNPs.}
#'   \item{\code{"lasso"}}{Fits a LASSO model with main effects and all
#'     pairwise interaction terms using \code{glmnet}. Identifies the
#'     interaction terms with non-zero coefficients at the cross-validated
#'     lambda. Recommended for blocks with p > 200 SNPs.}
#'   \item{\code{"auto"}}{Dispatches to \code{"pairwise"} when p <= 200,
#'     \code{"lasso"} otherwise.}
#' }
#'
#' @param block_id Character. The block to fine-map.
#' @param geno_matrix Numeric matrix (individuals x SNPs) or
#'   \code{LDxBlocks_backend}.
#' @param snp_info Data frame with columns \code{SNP}, \code{CHR}, \code{POS}.
#' @param blocks LD block table.
#' @param y_resid Named numeric vector of GRM-corrected REML residuals (from
#'   \code{test_block_haplotypes} null model). Names must match individual IDs.
#' @param method Character. One of \code{"pairwise"}, \code{"lasso"},
#'   \code{"auto"}. Default \code{"auto"}.
#' @param min_freq Numeric. Minimum MAF. Default \code{0.05}.
#' @param sig_threshold Numeric. Significance threshold applied to the
#'   p-value chosen by \code{sig_metric}. Default \code{0.05}.
#' @param sig_metric Character. Which p-value drives the \code{significant}
#'   flag in pairwise output. One of \code{"p_simplem_sidak"} (default,
#'   recommended), \code{"p_simplem"}, \code{"p_bonf"}, or \code{"p_fdr"}.
#'   All four p-value columns are always present regardless of this choice.
#'   Ignored for the \code{"lasso"} method (which returns \code{lasso_coef}
#'   and \code{selected} instead of p-values).
#' @param meff_percent_cut Numeric in (0, 1). Variance threshold for simpleM
#'   eigendecomposition. Default \code{0.995}.
#' @param lasso_nfolds Integer. CV folds for glmnet lambda selection.
#'   Default \code{5L}.
#' @param verbose Logical. Default \code{TRUE}.
#'
#' @return A data frame sorted by p_wald (pairwise) or |coefficient|
#'   (lasso). Columns: \code{SNP_i}, \code{SNP_j}, \code{POS_i},
#'   \code{POS_j}, \code{dist_bp}, \code{aa_effect} or \code{lasso_coef},
#'   \code{SE}, \code{t_stat}, \code{p_wald}, \code{p_bonf},
#'   \code{significant} (pairwise) or \code{selected} (lasso).
#'
#' @seealso \code{\link{scan_block_epistasis}},
#'   \code{\link{scan_block_by_block_epistasis}}
#' @export
fine_map_epistasis_block <- function(
    block_id,
    geno_matrix,
    snp_info,
    blocks,
    y_resid,
    method        = c("auto", "pairwise", "lasso"),
    min_freq      = 0.05,
    sig_threshold = 0.05,
    sig_metric    = c("p_simplem_sidak", "p_simplem", "p_bonf", "p_fdr"),
    meff_percent_cut = 0.995,
    lasso_nfolds  = 5L,
    verbose       = TRUE
) {
  method     <- match.arg(method)
  sig_metric <- match.arg(sig_metric)
  .log  <- function(...) if (verbose) message(sprintf("[fine_map_epi] %s", paste0(...)))

  # Normalise block table
  snp_info$CHR <- .norm_chr_hap(snp_info$CHR)
  blocks$CHR   <- .norm_chr_hap(blocks$CHR)
  if (!"start_bp" %in% names(blocks)) {
    blocks$start_bp <- blocks$start.bp
    blocks$end_bp   <- blocks$end.bp
  }
  if (!"block_id" %in% names(blocks))
    blocks$block_id <- paste0("block_", blocks$CHR, "_",
                              blocks$start_bp, "_", blocks$end_bp)

  brow <- blocks[blocks$block_id == block_id, , drop = FALSE]
  if (!nrow(brow)) stop("block_id '", block_id, "' not found.", call. = FALSE)
  chr_b <- brow$CHR[1]; start_b <- brow$start_bp[1]; end_b <- brow$end_bp[1]

  # SNPs in block
  si_b <- snp_info[snp_info$CHR == chr_b &
                     snp_info$POS >= start_b &
                     snp_info$POS <= end_b, , drop = FALSE]
  if (nrow(si_b) < 2L) stop("Fewer than 2 SNPs in block.", call. = FALSE)

  # Extract genotype matrix
  is_backend <- inherits(geno_matrix, "LDxBlocks_backend")
  if (is_backend) {
    snp_idx <- which(snp_info$CHR == chr_b &
                       snp_info$POS >= start_b &
                       snp_info$POS <= end_b)
    G_blk <- read_chunk(geno_matrix, snp_idx)
  } else {
    snp_idx <- which(colnames(geno_matrix) %in% si_b$SNP)
    if (!length(snp_idx))
      snp_idx <- which(snp_info$CHR == chr_b &
                         snp_info$POS >= start_b &
                         snp_info$POS <= end_b)
    G_blk <- geno_matrix[, snp_idx, drop = FALSE]
  }

  # Align individuals
  common <- intersect(rownames(G_blk), names(y_resid))
  if (length(common) < 10L) stop("Fewer than 10 common individuals.", call. = FALSE)
  G_blk  <- G_blk[common, , drop = FALSE]
  y_r    <- y_resid[common]

  # MAF filter -- keep si_b aligned with G_blk
  maf_b    <- apply(G_blk, 2L, function(x) {
    x <- x[!is.na(x)]; af <- mean(x) / 2; min(af, 1 - af)
  })
  maf_keep <- maf_b >= min_freq
  G_blk    <- G_blk[, maf_keep, drop = FALSE]
  if (ncol(G_blk) < 2L) stop("Fewer than 2 SNPs pass MAF filter.", call. = FALSE)
  if (!is.null(colnames(G_blk)) && all(colnames(G_blk) %in% si_b$SNP))
    si_b <- si_b[match(colnames(G_blk), si_b$SNP), , drop = FALSE]

  for (j in seq_len(ncol(G_blk))) {
    na_j <- is.na(G_blk[, j])
    if (any(na_j)) G_blk[na_j, j] <- mean(G_blk[, j], na.rm = TRUE)
  }
  col_var  <- apply(G_blk, 2L, stats::var)
  var_keep <- col_var > 1e-10
  G_blk    <- G_blk[, var_keep, drop = FALSE]
  if (!is.null(colnames(G_blk)) && all(colnames(G_blk) %in% si_b$SNP))
    si_b <- si_b[match(colnames(G_blk), si_b$SNP), , drop = FALSE]

  p_blk     <- ncol(G_blk)
  snp_names <- colnames(G_blk)
  if (is.null(snp_names)) snp_names <- si_b$SNP[seq_len(p_blk)]

  # Dispatch method
  if (method == "auto") method <- if (p_blk <= 200L) "pairwise" else "lasso"
  .log("Block: ", block_id, " | ", p_blk, " SNPs | method: ", method)

  # pos_map derived from aligned si_b (same row count as p_blk after filtering)
  pos_map <- if (nrow(si_b) == p_blk)
    stats::setNames(si_b$POS, snp_names) else
      stats::setNames(si_b$POS[seq_len(p_blk)], snp_names)

  if (method == "pairwise") {
    n_pairs <- p_blk * (p_blk - 1L) / 2L
    .log("Exhaustive pairwise scan: ", n_pairs, " pairs")
    scan_df <- .pairwise_interaction_scan(G_blk, y_r, col_names = snp_names)
    if (!nrow(scan_df)) {
      .log("No pairs passed variance filter.")
      return(data.frame())
    }
    scan_df$p_bonf  <- pmin(scan_df$p_wald * n_pairs, 1)
    scan_df$p_fdr   <- stats::p.adjust(scan_df$p_wald, method = "BH")
    # simpleM Meff from eigenspectrum of interaction columns
    meff_fm <- tryCatch(
      .compute_simplem_meff_matrix(
        G_blk, percent_cut = meff_percent_cut, max_cols = 1000L)$meff,
      error = function(e) n_pairs
    )
    if (!is.finite(meff_fm) || meff_fm < 1L) meff_fm <- n_pairs
    scan_df$Meff            <- meff_fm
    scan_df$p_simplem       <- .adjust_p_simplem_bonf(scan_df$p_wald, meff_fm)
    scan_df$p_simplem_sidak <- .adjust_p_simplem_sidak(scan_df$p_wald, meff_fm)
    scan_df$significant <- switch(sig_metric,
                                  p_simplem_sidak = !is.na(scan_df$p_simplem_sidak) & scan_df$p_simplem_sidak < sig_threshold,
                                  p_simplem       = !is.na(scan_df$p_simplem)       & scan_df$p_simplem       < sig_threshold,
                                  p_bonf          = !is.na(scan_df$p_bonf)          & scan_df$p_bonf          < sig_threshold,
                                  p_fdr           = !is.na(scan_df$p_fdr)           & scan_df$p_fdr           < sig_threshold
    )
    scan_df$POS_i   <- pos_map[scan_df$col_i]
    scan_df$POS_j   <- pos_map[scan_df$col_j]
    scan_df$dist_bp <- abs(scan_df$POS_j - scan_df$POS_i)
    names(scan_df)[names(scan_df) == "col_i"] <- "SNP_i"
    names(scan_df)[names(scan_df) == "col_j"] <- "SNP_j"
    col_order <- c("SNP_i","SNP_j","POS_i","POS_j","dist_bp",
                   "aa_effect","SE","t_stat","p_wald",
                   "Meff","p_bonf","p_fdr","p_simplem","p_simplem_sidak",
                   "significant")
    scan_df <- scan_df[, intersect(col_order, names(scan_df)), drop = FALSE]
    .log("Significant pairs (", sig_metric, "): ",
         sum(scan_df$significant, na.rm = TRUE))
    return(scan_df[order(scan_df$p_wald), ])

  } else {
    # LASSO with pairwise interactions
    if (!requireNamespace("glmnet", quietly = TRUE))
      stop("glmnet required for method='lasso': install.packages('glmnet')",
           call. = FALSE)

    .log("Building interaction design matrix (", p_blk,
         " SNPs, ", p_blk*(p_blk-1L)/2L, " interaction terms) ...")

    # Build [main effects | interaction terms] matrix
    pair_names <- character(p_blk * (p_blk - 1L) / 2L)
    pair_i     <- integer(length(pair_names))
    pair_j     <- integer(length(pair_names))
    k <- 0L
    for (i in seq_len(p_blk - 1L))
      for (j in seq(i + 1L, p_blk)) {
        k <- k + 1L
        pair_names[k] <- paste0(snp_names[i], ":", snp_names[j])
        pair_i[k] <- i; pair_j[k] <- j
      }

    X_int_mat <- matrix(0, nrow = nrow(G_blk), ncol = k)
    for (kk in seq_len(k))
      X_int_mat[, kk] <- G_blk[, pair_i[kk]] * G_blk[, pair_j[kk]]
    colnames(X_int_mat) <- pair_names[seq_len(k)]

    X_full <- cbind(G_blk, X_int_mat)
    X_full <- scale(X_full, center = TRUE, scale = TRUE)
    X_full[!is.finite(X_full)] <- 0

    .log("Fitting LASSO (", ncol(X_full), " terms, ", lasso_nfolds, "-fold CV) ...")
    cv_fit <- tryCatch(
      glmnet::cv.glmnet(x = X_full, y = y_r, alpha = 1,
                        nfolds = lasso_nfolds, intercept = TRUE),
      error = function(e) {
        stop("glmnet::cv.glmnet failed: ", conditionMessage(e), call. = FALSE)
      }
    )

    coef_cv <- as.numeric(
      glmnet::coef.glmnet(cv_fit, s = "lambda.1se")[-1L])   # drop intercept
    names(coef_cv) <- colnames(X_full)

    # Extract non-zero interaction terms
    int_coefs <- coef_cv[(p_blk + 1L):length(coef_cv)]
    selected  <- names(int_coefs)[abs(int_coefs) > 1e-10]

    if (!length(selected)) {
      .log("No interaction terms selected by LASSO at lambda.1se.")
      return(data.frame())
    }

    .log(length(selected), " interaction terms selected by LASSO")

    # Build output table
    out_rows <- lapply(selected, function(nm) {
      parts <- strsplit(nm, ":", fixed = TRUE)[[1L]]
      if (length(parts) != 2L) return(NULL)
      si <- parts[1L]; sj <- parts[2L]
      data.frame(
        SNP_i      = si, SNP_j = sj,
        POS_i      = pos_map[si], POS_j = pos_map[sj],
        dist_bp    = abs(pos_map[sj] - pos_map[si]),
        lasso_coef = round(int_coefs[nm], 6),
        selected   = TRUE,
        stringsAsFactors = FALSE
      )
    })
    out_df <- do.call(rbind, out_rows[!vapply(out_rows, is.null, logical(1L))])
    rownames(out_df) <- NULL
    return(out_df[order(-abs(out_df$lasso_coef)), ])
  }
}


# ==============================================================================
# print methods
# ==============================================================================

#' @export
print.LDxBlocks_epistasis <- function(x, ...) {
  cat("LDxBlocks Within-Block Epistasis Scan\n")
  cat("  Blocks scanned  :", x$n_blocks_scanned, "\n")
  cat("  Total pairs     :", format(x$n_pairs_total, big.mark = ","), "\n")
  if (nrow(x$results) > 0) {
    sm_label <- if (!is.null(x$sig_metric)) x$sig_metric else "sig_metric"
    cat("  Significant pairs (", sm_label, "):",
        sum(x$results$significant, na.rm = TRUE), "\n", sep = "")
    if (nrow(x$scan_summary) > 0) {
      cat("  Per-block summary:\n")
      for (i in seq_len(nrow(x$scan_summary))) {
        cat(sprintf("    %-40s  %d pairs  %d significant\n",
                    x$scan_summary$block_id[i],
                    x$scan_summary$n_pairs[i],
                    x$scan_summary$n_significant[i]))
      }
    }
  } else {
    cat("  No significant interactions found.\n")
  }
  invisible(x)
}

#' @export
print.LDxBlocks_block_epistasis <- function(x, ...) {
  cat("LDxBlocks Between-Block Haplotype Epistasis Scan\n")
  cat("  Query alleles   :", x$n_sig_alleles, "\n")
  cat("  Total tests     :", format(x$n_tests, big.mark = ","), "\n")
  if (nrow(x$results) > 0) {
    sm_label_bb <- if (!is.null(x$sig_metric)) x$sig_metric else "sig_metric"
    cat("  Significant (", sm_label_bb, "):",
        sum(x$results$significant, na.rm = TRUE), "\n", sep = "")
    cat("  Top 5 interactions by p_wald:\n")
    top5 <- head(x$results[, c("block_i","allele_i","block_j","allele_j",
                               "aa_effect","p_wald","p_bonf")], 5)
    print(top5, row.names = FALSE)
  } else {
    cat("  No significant interactions found.\n")
  }
  invisible(x)
}
