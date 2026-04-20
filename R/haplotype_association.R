# ==============================================================================
# haplotype_association.R
# Block-level haplotype association testing and diplotype effect estimation.
#
# Statistical model (test_block_haplotypes)
# -----------------------------------------
# Unified mixed linear model (Q+K, EMMAX/GAPIT3 formulation):
#
#   y = mu + alpha * x_hap + sum(beta_k * PC_k) + g + e
#
#   PC_k = k-th eigenvector of the haplotype GRM (GRM-derived PCA)
#   g    ~ MVN(0, sigma_g^2 * G_hap)
#   n_pcs = 0: pure GRM correction (EMMAX)
#   n_pcs > 0: Q+K model (Yu et al. 2006)
#
# Scaling design
# --------------
# The GRM inversion via rrBLUP::mixed.solve() is O(n^3) but done ONCE per
# trait. The per-allele scan on de-regressed residuals is then fully
# vectorised: all haplotype allele tests across all blocks are performed in
# a single crossprod() call (analogous to the BLAS DGEMV trick in
# screening.R). This makes the scan step O(n_ind * n_hap_cols) regardless
# of the number of blocks, with no R-level function call overhead per allele.
#
# For 5000 individuals and 51,000 haplotype allele columns (17,000 blocks
# x 3 alleles), the scan completes in seconds. The dominant cost is the
# O(n^3) GRM inversion - done once, not per block.
# ==============================================================================


# ==============================================================================
# Internal helper: parse blues argument into a named list of numeric vectors.
# Defined here so haplotype_association.R is self-contained and does not depend
# on the copy in haplotype_analysis.R being loaded first.
# ==============================================================================
.parse_blues_assoc <- function(blues, id_col, blue_col, blue_cols) {
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


# ==============================================================================
# 1. test_block_haplotypes
# ==============================================================================

#' Block-Level Haplotype Association Testing (Q+K Mixed Linear Model)
#'
#' @description
#' Performs genome-wide haplotype block association tests for one or more
#' quantitative traits. Each LD block is tested as a unit: per-allele Wald
#' tests identify which specific haplotype alleles drive association, and an
#' omnibus F-test evaluates the block as a whole. Population structure and
#' kinship are corrected jointly through a unified mixed linear model.
#'
#' \strong{Statistical model (Q+K / EMMAX formulation):}
#'
#' \deqn{y = \mu + \alpha \cdot x_{\mathrm{hap}} +
#'        \sum_{k=1}^{K} \beta_k \, PC_k + g + \varepsilon}
#'
#' \describe{
#'   \item{\eqn{x_{\mathrm{hap}}}}{Haplotype allele dosage (0, 1, or 2 copies
#'     for phased data; 0 or 1 for unphased) - the quantity being tested for
#'     association.}
#'   \item{\eqn{PC_k}}{The k-th eigenvector of the haplotype GRM, included as
#'     a fixed-effect covariate to capture discrete population structure
#'     explicitly. Derived from \code{eigen(G_hap)}, so both fixed-effect PCs
#'     and the random-effect GRM use the same kinship model.}
#'   \item{\eqn{g \sim MVN(0,\,\sigma_g^2 G)}}{Polygenic background - captures
#'     residual continuous kinship (within-family, cryptic relatedness) as a
#'     random effect after PC removal.}
#'   \item{\eqn{\varepsilon \sim MVN(0,\,\sigma_e^2 I)}}{Residual error.}
#' }
#'
#' \strong{Implementation and scaling:}
#' The GRM is inverted once per trait via \code{rrBLUP::mixed.solve()} (dominant
#' cost, O(n^3)). Per-allele tests on the de-regressed residuals are then fully
#' vectorised across all blocks simultaneously using a single
#' \code{crossprod()} call (O(n x p), analogous to the BLAS DGEMV trick in
#' marginal SNP screening). For 5,000 individuals and 51,000 haplotype allele
#' columns across 17,000 blocks, the scan step takes seconds after the
#' one-time GRM inversion (~30 s for n = 5,000).
#'
#' @param haplotypes Named list produced by \code{\link{extract_haplotypes}}.
#'   Each element is a named character vector (one haplotype dosage string per
#'   individual), and the list must carry a \code{block_info} attribute (added
#'   automatically by \code{extract_haplotypes}). The number of elements equals
#'   the number of qualifying LD blocks.
#'
#' @param blues Pre-adjusted phenotype means (BLUEs or BLUPs from a field
#'   trial mixed model). Accepted in four formats:
#'   \itemize{
#'     \item Named numeric vector: \code{c(ind1 = 2.3, ind2 = 1.8, ...)}.
#'       Names must match individual IDs (row names of the genotype matrix).
#'     \item Single-trait data frame: columns \code{id_col} (individual IDs)
#'       and \code{blue_col} (numeric phenotype values).
#'     \item Multi-trait data frame: columns \code{id_col} plus one column per
#'       trait specified in \code{blue_cols}. All named traits are tested in a
#'       single call sharing the same GRM.
#'     \item Named list of named numeric vectors: one element per trait,
#'       e.g. \code{list(YLD = c(ind1=2.3,...), RES = c(ind1=0.8,...))}.
#'   }
#'   Individuals in \code{blues} not present in \code{haplotypes} are silently
#'   dropped. At least 10 common individuals are required per trait.
#'
#' @param blocks LD block table returned by \code{\link{run_Big_LD_all_chr}}
#'   or \code{\link{run_ldx_pipeline}}. Required columns: \code{CHR},
#'   \code{start.bp}, \code{end.bp}. Used only for block metadata annotation
#'   in the output; the actual haplotype allele columns come from
#'   \code{haplotypes}.
#'
#' @param n_pcs Integer (\code{>= 0}) or \code{NULL}. Number of haplotype-GRM
#'   eigenvectors to include as fixed-effect population structure covariates
#'   in the null model. Controls the trade-off between correction strategies:
#'   \itemize{
#'     \item \code{0L} (default): Pure GRM correction - EMMAX / P3D
#'       approximation. The GRM random effect absorbs all structure and kinship.
#'       Appropriate for populations with diffuse continuous kinship
#'       (e.g. livestock half-sib families, advanced inbred lines) where
#'       there are no sharp discrete subpopulation boundaries.
#'     \item \code{1} to \code{10}: Q+K model (Yu et al. 2006) - top-k GRM
#'       eigenvectors as fixed effects plus the GRM random effect. The fixed-
#'       effect PCs capture discrete subpopulation membership (e.g. breeds,
#'       ecotypes, geographic clusters); the GRM random effect captures
#'       within-subpopulation continuous kinship. Use when strong cluster
#'       structure inflates the Q-Q plot under \code{n_pcs = 0}.
#'       Typically 3-5 PCs are sufficient; using more than 10 risks
#'       over-correction.
#'     \item \code{NULL}: Auto-select via the elbow of the GRM eigenvalue scree
#'       plot (first position where the marginal gain in variance explained
#'       drops below 1\%), capped at 10.
#'   }
#'   PCs are derived from \code{eigen(G_hap)} - the same GRM that enters as
#'   the random effect - ensuring mathematical consistency.
#'
#' @param top_n Integer or \code{NULL}. Maximum number of haplotype alleles per
#'   block to include in the feature matrix before testing. \code{NULL}
#'   (default) retains all alleles above \code{min_freq}. Supply an integer
#'   (e.g. \code{5L}) only to cap column count for very large panels where
#'   memory is constrained. Alleles are ranked by frequency; the most common
#'   \code{top_n} are retained.
#'
#' @param min_freq Numeric in (0, 1). Minimum haplotype allele frequency in
#'   the panel. Alleles with frequency below this threshold are excluded from
#'   both the feature matrix and all tests. Default \code{0.05}. Lower values
#'   (e.g. \code{0.02}) include rarer alleles at the cost of reduced power and
#'   increased multiple testing burden. Values below \code{0.02} are not
#'   recommended for panels smaller than 200 individuals.
#'
#' @param id_col Character. Name of the individual-ID column when \code{blues}
#'   is a data frame. Must exactly match a column name in the data frame.
#'   Default \code{"id"}.
#'
#' @param blue_col Character. Name of the phenotype column when \code{blues}
#'   is a single-trait data frame. Must be numeric. Default \code{"blue"}.
#'
#' @param blue_cols Character vector. Names of phenotype columns when
#'   \code{blues} is a multi-trait wide data frame. Each named column is
#'   treated as a separate trait and tested independently (but using the same
#'   GRM). Default \code{NULL} (ignored unless \code{blues} is a data frame
#'   with more than one numeric column beyond \code{id_col}).
#'
#' @param alpha Numeric or \code{NULL}. Significance threshold for the
#'   \code{significant} flag in \code{allele_tests} and
#'   \code{significant_omnibus} in \code{block_tests}. \code{NULL} (default)
#'   applies genome-wide Bonferroni correction:
#'   \eqn{\alpha = 0.05 / n_{\mathrm{tests}}} where \eqn{n_{\mathrm{tests}}}
#'   is the total number of allele-level tests across all blocks and traits.
#'   Supply a fixed value (e.g. \code{0.05}) to use a less conservative
#'   threshold, or \code{0.05 / nrow(blocks)} to apply Bonferroni at the block
#'   level rather than the allele level.
#'
#' @param verbose Logical. If \code{TRUE} (default), prints timestamped
#'   progress messages: model type, trait name, number of alleles scanned, and
#'   a final summary. Set \code{FALSE} for batch use or inside loops.
#'
#' @return A named list of class \code{c("LDxBlocks_haplotype_assoc", "list")}
#'   with the following elements:
#'
#' \describe{
#'   \item{\code{allele_tests}}{Data frame. One row per haplotype allele per
#'     LD block per trait. Contains 13 columns:
#'     \itemize{
#'       \item \code{block_id} (character) - Block identifier string matching
#'         \code{names(haplotypes)}, e.g. \code{"block_1_1000_103000"}.
#'       \item \code{CHR} (character) - Chromosome label.
#'       \item \code{start_bp} (integer) - Block start coordinate (base pairs).
#'       \item \code{end_bp} (integer) - Block end coordinate (base pairs).
#'       \item \code{trait} (character) - Trait name.
#'       \item \code{allele} (character) - Haplotype allele identifier string.
#'       \item \code{frequency} (numeric, (0,1]) - Allele frequency in the
#'         panel (proportion of individuals carrying the allele).
#'       \item \code{effect} (numeric) - Estimated additive effect: mean
#'         phenotype difference per unit increase in allele dosage on the
#'         de-regressed residual scale. Positive = favourable if higher trait
#'         values are desirable.
#'       \item \code{SE} (numeric) - Standard error of the effect estimate.
#'       \item \code{t_stat} (numeric) - t-statistic (\code{effect / SE}).
#'       \item \code{p_wald} (numeric, (0,1]) - Two-sided Wald p-value (raw,
#'         uncorrected).
#'       \item \code{p_wald_adj} (numeric, (0,1]) - Bonferroni-adjusted p-value:
#'         \code{min(p_wald * n_tests, 1)}.
#'       \item \code{significant} (logical) - \code{TRUE} when
#'         \code{p_wald <= alpha}.
#'     }
#'     Rows are sorted ascending by \code{CHR}, \code{start_bp}, then
#'     \code{p_wald} within each trait.}
#'
#'   \item{\code{block_tests}}{Data frame. One row per LD block per trait.
#'     Contains 12 columns:
#'     \itemize{
#'       \item \code{block_id}, \code{CHR}, \code{start_bp}, \code{end_bp},
#'         \code{trait} - as above.
#'       \item \code{n_alleles_tested} (integer) - Number of alleles above
#'         \code{min_freq} in this block that were included in the omnibus test.
#'       \item \code{F_stat} (numeric) - Omnibus F-statistic testing all
#'         \code{n_alleles_tested} allele columns jointly against the de-
#'         regressed residual.
#'       \item \code{df_LRT} (integer) - Numerator degrees of freedom
#'         (= \code{n_alleles_tested}).
#'       \item \code{p_omnibus} (numeric, (0,1]) - Raw omnibus p-value from
#'         the F-distribution.
#'       \item \code{p_omnibus_adj} (numeric, (0,1]) - Bonferroni-adjusted
#'         across all blocks per trait: \code{min(p_omnibus * n_blocks, 1)}.
#'       \item \code{var_explained} (numeric, [0,1]) - Proportion of
#'         de-regressed phenotypic variance explained by all haplotype alleles
#'         of this block jointly: \code{1 - RSS_full / RSS_null}.
#'       \item \code{significant_omnibus} (logical) - \code{TRUE} when
#'         \code{p_omnibus_adj < 0.05}.
#'     }
#'     Rows sorted ascending by \code{CHR}, \code{start_bp}, \code{p_omnibus}.}
#'
#'   \item{\code{traits}}{Character vector of trait names that were tested.}
#'
#'   \item{\code{n_pcs_used}}{Integer. Number of GRM eigenvectors included as
#'     fixed-effect covariates in the null model. \code{0L} when
#'     \code{n_pcs = 0L} (pure EMMAX).}
#'
#'   \item{\code{alpha}}{Numeric. The significance threshold actually used
#'     (either supplied or computed by Bonferroni).}
#'
#'   \item{\code{n_tests}}{Integer. Total number of allele-level Wald tests
#'     performed across all blocks and traits (the Bonferroni denominator).}
#' }
#'
#' @examples
#' \donttest{
#' data(ldx_geno, ldx_snp_info, ldx_blocks, ldx_blues, package = "LDxBlocks")
#' haps <- extract_haplotypes(ldx_geno, ldx_snp_info, ldx_blocks)
#'
#' # Pure GRM correction (EMMAX, default - n_pcs = 0)
#' assoc <- test_block_haplotypes(
#'   haplotypes = haps,
#'   blues      = setNames(ldx_blues$YLD, ldx_blues$id),
#'   blocks     = ldx_blocks,
#'   verbose    = FALSE
#' )
#' head(assoc$block_tests[order(assoc$block_tests$p_omnibus), ])
#'
#' # Q+K model (3 GRM-derived PCs + GRM, for structured populations)
#' assoc_qk <- test_block_haplotypes(
#'   haplotypes = haps,
#'   blues      = setNames(ldx_blues$YLD, ldx_blues$id),
#'   blocks     = ldx_blocks,
#'   n_pcs      = 3L,
#'   verbose    = FALSE
#' )
#'
#' # Multi-trait
#' assoc_mt <- test_block_haplotypes(
#'   haplotypes = haps,
#'   blues      = ldx_blues,
#'   blocks     = ldx_blocks,
#'   id_col     = "id",
#'   blue_cols  = c("YLD", "RES"),
#'   n_pcs      = 3L,
#'   verbose    = FALSE
#' )
#' }
#' @references
#' Endelman JB (2011). Ridge regression and other kernels for genomic
#' selection with R package rrBLUP. \emph{Plant Genome} \strong{4}:250-255.
#' \doi{10.3835/plantgenome2011.08.0024}
#'
#' Yu J, Pressoir G, Briggs WH, et al. (2006). A unified mixed-model method
#' for association mapping that accounts for multiple levels of relatedness.
#' \emph{Nature Genetics} \strong{38}(2):203-208. \doi{10.1038/ng1702}
#'
#' Kang HM, Sul JH, Service SK, et al. (2010). Variance component model to
#' account for sample structure in genome-wide association studies.
#' \emph{Nature Genetics} \strong{42}(4):348-354. \doi{10.1038/ng.548}
#' @seealso \code{\link{extract_haplotypes}},
#'   \code{\link{compute_haplotype_grm}},
#'   \code{\link{run_haplotype_prediction}},
#'   \code{\link{estimate_diplotype_effects}}
#' @export
test_block_haplotypes <- function(
    haplotypes,
    blues,
    blocks,
    n_pcs    = 0L,
    top_n    = NULL,
    min_freq = 0.05,
    id_col   = "id",
    blue_col = "blue",
    blue_cols = NULL,
    alpha    = NULL,
    verbose  = TRUE
) {
  if (!requireNamespace("rrBLUP", quietly = TRUE))
    stop("rrBLUP is required: install.packages('rrBLUP')", call. = FALSE)

  .log <- function(...) if (verbose) message(sprintf("[assoc] %s", paste0(...)))

  blues_list <- .parse_blues_assoc(blues, id_col, blue_col, blue_cols)
  traits     <- names(blues_list)
  bi         <- attr(haplotypes, "block_info")

  # -- Build haplotype feature matrix and GRM (once, shared across all traits)
  .log("Building haplotype feature matrix ...")
  hap_mat <- build_haplotype_feature_matrix(
    haplotypes, top_n = top_n, min_freq = min_freq, encoding = "additive_012"
  )$matrix
  # Column names are "blockID__hapN" - extract block ID prefix for grouping
  col_blocks <- sub("__hap\\d+$", "", colnames(hap_mat))

  .log("Computing haplotype GRM ...")
  G <- compute_haplotype_grm(hap_mat)

  # -- Derive GRM PCs once (shared across all traits) -------------------------
  # G = Q Lambda Q^T - columns of Q are PCs of individuals in haplotype space.
  # Mathematically equivalent to PCA on genotypes (SVD duality) but guaranteed
  # consistent with the GRM random effect - same kinship model for both.
  n_pcs_used <- 0L
  pc_scores  <- NULL
  auto_pcs   <- is.null(n_pcs)
  n_pcs_req  <- if (auto_pcs) NA_integer_ else as.integer(n_pcs)

  if (auto_pcs || n_pcs_req > 0L) {
    .log("Eigendecomposing GRM for PC fixed-effect covariates ...")
    eig      <- eigen(G, symmetric = TRUE)
    eig_vals <- pmax(eig$values, 0)
    var_prop <- eig_vals / sum(eig_vals)
    if (auto_pcs) {
      marginal   <- diff(-var_prop)
      elbow      <- which(marginal > 0.01)[1]
      n_pcs_used <- if (is.na(elbow)) 3L else min(elbow, 10L)
    } else {
      n_pcs_used <- min(n_pcs_req, ncol(eig$vectors))
    }
    if (n_pcs_used > 0L) {
      pc_scores <- eig$vectors[, seq_len(n_pcs_used), drop = FALSE]
      rownames(pc_scores) <- rownames(G)
      colnames(pc_scores) <- paste0("PC", seq_len(n_pcs_used))
      cum_var <- sum(var_prop[seq_len(n_pcs_used)]) * 100
      .log("  ", n_pcs_used, " GRM-derived PCs (", round(cum_var, 1),
           "% of GRM variance)")
    }
  }

  allele_rows <- list()
  block_rows  <- list()

  for (tr in traits) {
    model_str <- if (n_pcs_used > 0L)
      paste0("GRM + ", n_pcs_used, " PC(s) [Q+K]") else "GRM [EMMAX]"
    .log("Testing trait: ", tr, " [", model_str, "]")

    pheno  <- blues_list[[tr]]
    common <- intersect(names(pheno), rownames(G))
    if (length(common) < 10L) {
      warning("Fewer than 10 common individuals for trait '", tr,
              "' - skipping.", call. = FALSE)
      next
    }

    y     <- pheno[common]
    G_sub <- G[common, common, drop = FALSE]

    # -- Null model: y = mu + PC_fixed_effects + g + e ----------------------
    # X_null includes intercept + optional PC covariates.
    # mixed.solve() with X argument fits: y = X*beta + g + e
    if (n_pcs_used > 0L) {
      X_null <- cbind(1, pc_scores[common, , drop = FALSE])
    } else {
      X_null <- matrix(1, nrow = length(common), ncol = 1L)
    }

    null_fit <- tryCatch(
      rrBLUP::mixed.solve(y = y, X = X_null, K = G_sub, method = "REML"),
      error = function(e) NULL
    )
    if (is.null(null_fit)) {
      warning("Null model failed for trait '", tr, "' - skipping.",
              call. = FALSE)
      next
    }

    # De-regressed phenotype: remove fixed effects (mu + PCs) and polygenic g
    y_fixed <- as.numeric(X_null %*% null_fit$beta)
    y_resid <- y - y_fixed - null_fit$u[common]
    n_ind   <- length(y_resid)

    # -- VECTORISED per-allele scan across ALL blocks simultaneously ---------
    # This is the key scaling improvement over the original lm() loop.
    # Identical in spirit to screening.R: one BLAS crossprod() for all tests.
    #
    # Model per allele j:  y_resid ~ x_j
    # Normal equations:    beta_j = (x_j' y_resid) / (x_j' x_j)
    # RSS_j = y_resid' y_resid - beta_j^2 * (x_j' x_j)
    # sigma2_j = RSS_j / (n - 2)
    # SE_j = sqrt(sigma2_j / (x_j' x_j))

    X_all <- hap_mat[common, , drop = FALSE]

    # Mean-impute NA across all columns at once
    for (j in seq_len(ncol(X_all))) {
      na_j <- is.na(X_all[, j])
      if (any(na_j)) X_all[na_j, j] <- mean(X_all[, j], na.rm = TRUE)
    }

    # Remove constant columns
    XtX_vec <- colSums(X_all^2) - (colSums(X_all))^2 / n_ind  # variance * (n-1)
    valid_cols <- XtX_vec > 1e-10

    X_valid     <- X_all[, valid_cols, drop = FALSE]
    XtX_valid   <- XtX_vec[valid_cols]           # x_j'x_j (centred)
    col_bl_valid <- col_blocks[valid_cols]

    # Single BLAS call for all allele effects
    Xty    <- as.numeric(crossprod(X_valid, y_resid))   # x_j' y_resid
    # For centred regression: need crossprod with centred X
    # Approximate: beta = (X'y - n*xbar*ybar) / (X'X - n*xbar^2)
    # Equivalent to the screening.R vectorised formula
    x_bars  <- colMeans(X_valid)
    y_bar_r <- mean(y_resid)
    Xty_c   <- Xty - n_ind * x_bars * y_bar_r   # centred cross-product
    beta_all <- Xty_c / XtX_valid                # all betas in one step

    ss_tot   <- sum((y_resid - y_bar_r)^2)
    ss_reg   <- beta_all^2 * XtX_valid           # regression SS per allele
    rss_all  <- ss_tot - ss_reg
    rss_all  <- pmax(rss_all, 0)                 # numerical floor
    sigma2   <- rss_all / (n_ind - 2L)
    SE_all   <- sqrt(sigma2 / XtX_valid)
    t_all    <- beta_all / SE_all
    p_all    <- 2 * stats::pt(-abs(t_all), df = n_ind - 2L)
    freq_all <- (colMeans(X_valid, na.rm = TRUE)) / 2

    .log("  Allele scan: ", sum(valid_cols), " alleles across ",
         length(unique(col_bl_valid)), " blocks")

    # -- Collect per-allele results ------------------------------------------
    valid_idx <- which(valid_cols)
    for (ci in seq_along(valid_idx)) {
      orig_ci <- valid_idx[ci]
      bn      <- col_blocks[orig_ci]
      blk_meta <- if (!is.null(bi)) bi[bi$block_id == bn, , drop=FALSE]
      else NULL
      gm <- function(col, def = NA) {
        if (!is.null(blk_meta) && nrow(blk_meta)>0) blk_meta[[col]][1] else def
      }
      allele_label <- sub(paste0("^", bn, "__hap\\d+__?"), "",
                          colnames(hap_mat)[orig_ci])
      if (allele_label == "" || allele_label == colnames(hap_mat)[orig_ci])
        allele_label <- colnames(hap_mat)[orig_ci]

      allele_rows[[length(allele_rows)+1L]] <- data.frame(
        block_id  = bn,
        CHR       = gm("CHR",      NA_character_),
        start_bp  = gm("start_bp", NA_integer_),
        end_bp    = gm("end_bp",   NA_integer_),
        trait     = tr,
        allele    = allele_label,
        frequency = round(freq_all[ci], 4),
        effect    = round(beta_all[ci],  6),
        SE        = round(SE_all[ci],    6),
        t_stat    = round(t_all[ci],     4),
        p_wald    = p_all[ci],
        stringsAsFactors = FALSE
      )
    }

    # -- Per-block omnibus F-test (vectorised per block) ---------------------
    # After the vectorised scan, omnibus F-test is a per-block loop but cheap:
    # only iterates over ~17,000 blocks with 3-5 columns each.
    unique_blocks <- unique(col_bl_valid)
    for (bn in unique_blocks) {
      blk_idx <- which(col_bl_valid == bn)
      X_blk   <- X_valid[, blk_idx, drop = FALSE]
      if (!ncol(X_blk)) next

      blk_meta <- if (!is.null(bi)) bi[bi$block_id == bn, , drop=FALSE]
      else NULL
      gm2 <- function(col, def=NA) {
        if (!is.null(blk_meta) && nrow(blk_meta)>0) blk_meta[[col]][1] else def
      }

      df_full <- as.data.frame(X_blk)
      names(df_full) <- paste0("h", seq_len(ncol(df_full)))
      df_full$y_r <- y_resid
      fit_full <- tryCatch(stats::lm(y_r ~ ., data=df_full),
                           error=function(e) NULL)
      if (is.null(fit_full)) next

      RSS_null <- ss_tot
      RSS_full <- sum(stats::residuals(fit_full)^2)
      df_lrt   <- ncol(X_blk)
      df_res   <- n_ind - df_lrt - 1L
      if (RSS_null < 1e-15 || df_lrt < 1L || df_res < 1L) next

      F_stat  <- ((RSS_null - RSS_full)/df_lrt) / (RSS_full/df_res)
      p_omni  <- stats::pf(F_stat, df1=df_lrt, df2=df_res, lower.tail=FALSE)
      var_exp <- max(0, 1 - RSS_full/RSS_null)

      block_rows[[length(block_rows)+1L]] <- data.frame(
        block_id         = bn,
        CHR              = gm2("CHR",      NA_character_),
        start_bp         = gm2("start_bp", NA_integer_),
        end_bp           = gm2("end_bp",   NA_integer_),
        trait            = tr,
        n_alleles_tested = ncol(X_blk),
        F_stat           = round(F_stat,  4),
        df_LRT           = df_lrt,
        p_omnibus        = p_omni,
        var_explained    = round(var_exp, 4),
        stringsAsFactors = FALSE
      )
    }
  }

  # -- Assemble output and apply multiple testing correction ------------------
  if (!length(allele_rows) && !length(block_rows)) {
    message("[test_block_haplotypes] No results produced.")
    return(structure(
      list(allele_tests=data.frame(), block_tests=data.frame(),
           traits=traits, n_pcs_used=n_pcs_used, alpha=alpha, n_tests=0L),
      class=c("LDxBlocks_haplotype_assoc","list")
    ))
  }

  allele_df <- do.call(rbind, allele_rows)
  block_df  <- do.call(rbind, block_rows)
  n_tests   <- nrow(allele_df)
  if (is.null(alpha)) alpha <- 0.05 / max(n_tests, 1L)

  if (nrow(allele_df) > 0) {
    allele_df$p_wald_adj  <- pmin(allele_df$p_wald * n_tests, 1)
    allele_df$significant <- allele_df$p_wald <= alpha
    allele_df <- allele_df[order(allele_df$CHR, allele_df$start_bp,
                                 allele_df$p_wald), ]
  }
  if (nrow(block_df) > 0) {
    n_bt <- nrow(block_df)
    block_df$p_omnibus_adj       <- pmin(block_df$p_omnibus * n_bt, 1)
    block_df$significant_omnibus <- block_df$p_omnibus_adj < 0.05
    block_df <- block_df[order(block_df$CHR, block_df$start_bp,
                               block_df$p_omnibus), ]
  }

  .log("Done. Allele tests: ", nrow(allele_df),
       " | Significant: ", sum(allele_df$significant, na.rm=TRUE),
       " | Alpha: ", formatC(alpha, format="e", digits=2))

  structure(
    list(allele_tests=allele_df, block_tests=block_df,
         traits=traits, n_pcs_used=n_pcs_used, alpha=alpha, n_tests=n_tests),
    class=c("LDxBlocks_haplotype_assoc","list")
  )
}

#' @export
print.LDxBlocks_haplotype_assoc <- function(x, ...) {
  cat("LDxBlocks Haplotype Association Results\n")
  cat("  Model:       GRM")
  if (x$n_pcs_used > 0L)
    cat(" + ", x$n_pcs_used, " GRM-derived PC(s) [Q+K]", sep="")
  else
    cat(" [EMMAX/P3D]")
  cat("\n")
  cat("  Traits:      ", paste(x$traits, collapse=", "), "\n")
  cat("  Allele tests:", nrow(x$allele_tests), "\n")
  cat("  Block tests: ", nrow(x$block_tests), "\n")
  cat("  Alpha:       ", formatC(x$alpha, format="e", digits=3), "\n")
  if (nrow(x$allele_tests) > 0)
    cat("  Significant alleles:", sum(x$allele_tests$significant, na.rm=TRUE), "\n")
  if (nrow(x$block_tests) > 0)
    cat("  Significant blocks (omnibus):",
        sum(x$block_tests$significant_omnibus, na.rm=TRUE), "\n")
  invisible(x)
}


# ==============================================================================
# 2. estimate_diplotype_effects
# ==============================================================================

#' Estimate Diplotype Effects and Dominance Deviations Per LD Block
#'
#' @description
#' Fits a per-block diplotype model to decompose phenotypic variation into
#' additive and dominance components at each LD block, after correcting for
#' population structure and kinship via the haplotype GRM. This function
#' provides direct evidence for non-additive gene action: dominance,
#' overdominance, and heterosis.
#'
#' \strong{Method:} For each LD block, unique diplotypes (canonical
#' haplotype-pair strings inferred by \code{\link{infer_block_haplotypes}},
#' e.g. \code{"010/110"}) are treated as levels of a fixed factor. A mixed
#' model with the haplotype GRM as a kinship random effect is fitted via
#' \code{rrBLUP::mixed.solve()} to obtain de-regressed phenotype residuals.
#' Diplotype class means are computed on these residuals. For each pair of
#' common alleles (A, B), classical quantitative genetics parameters are
#' derived from the three diplotype class means:
#' \itemize{
#'   \item Additive effect: \eqn{a = (\bar{y}_{BB} - \bar{y}_{AA}) / 2}
#'   \item Dominance deviation: \eqn{d = \bar{y}_{AB} -
#'     (\bar{y}_{AA} + \bar{y}_{BB}) / 2}
#'   \item Dominance ratio: \eqn{d/a} - the key quantity for interpreting
#'     gene action (see \code{d_over_a} column description).
#' }
#'
#' @param haplotypes Named list produced by \code{\link{extract_haplotypes}}.
#'   Phase-ambiguous diplotypes (unphased heterozygotes) are excluded from
#'   the dominance analysis but included in the omnibus F-test if their
#'   diplotype string is unambiguous. For accurate dominance decomposition,
#'   phased input (from \code{read_phased_vcf}, \code{phase_with_beagle}, or
#'   \code{phase_with_pedigree}) is strongly recommended.
#'
#' @param blues Pre-adjusted phenotype means. Accepts the same four formats
#'   as \code{\link{test_block_haplotypes}}: named numeric vector,
#'   single-trait data frame, multi-trait data frame, or named list. All
#'   formats require individual IDs matching the names of haplotype strings.
#'
#' @param blocks LD block table from \code{\link{run_Big_LD_all_chr}}.
#'   Used only to supply block coordinate metadata (\code{CHR},
#'   \code{start_bp}, \code{end_bp}) to the output tables.
#'
#' @param min_freq Numeric in (0, 1). Minimum allele frequency. Alleles
#'   below this threshold in the full panel are excluded from the dominance
#'   decomposition (\code{dominance_table}) but their diplotypes still
#'   contribute to \code{diplotype_means} and \code{omnibus_tests} if
#'   sufficient individuals carry them. Default \code{0.05}.
#'
#' @param min_n_diplotype Integer. Minimum number of individuals that must
#'   carry a given diplotype class for it to be included in \code{diplotype_means}
#'   and the dominance decomposition. Diplotype classes with fewer individuals
#'   are silently excluded. The omnibus F-test uses only the retained classes.
#'   Default \code{3L}. Increasing to 5-10 improves mean estimate stability
#'   for small panels; decreasing below 3 is not recommended.
#'
#' @param id_col Character. Name of the individual-ID column when \code{blues}
#'   is a data frame. Default \code{"id"}.
#'
#' @param blue_col Character. Name of the phenotype column for single-trait
#'   data frames. Default \code{"blue"}.
#'
#' @param blue_cols Character vector. Phenotype column names for multi-trait
#'   data frames. Default \code{NULL}.
#'
#' @param verbose Logical. \code{TRUE} (default) prints progress per trait.
#'
#' @return A named list of class \code{c("LDxBlocks_diplotype", "list")}
#'   with three elements:
#'
#' \describe{
#'   \item{\code{diplotype_means}}{Data frame. One row per diplotype class per
#'     LD block per trait, for diplotype classes with \code{>= min_n_diplotype}
#'     individuals. Contains 9 columns:
#'     \itemize{
#'       \item \code{block_id} (character) - Block identifier.
#'       \item \code{CHR} (character) - Chromosome.
#'       \item \code{start_bp}, \code{end_bp} (integer) - Block coordinates.
#'       \item \code{trait} (character) - Trait name.
#'       \item \code{diplotype} (character) - Canonical diplotype string: two
#'         haplotype allele strings sorted alphabetically and joined by
#'         \code{"/"}, e.g. \code{"010/110"}. Homozygotes have identical
#'         halves: \code{"010/010"}. This format ensures that \code{"010/110"}
#'         and \code{"110/010"} are always represented identically.
#'       \item \code{n} (integer) - Number of individuals with this diplotype.
#'       \item \code{mean_blue} (numeric) - Mean de-regressed phenotype value
#'         for this diplotype class (on the residual scale after GRM correction).
#'       \item \code{se_mean} (numeric) - Standard error of the mean
#'         (\code{sd / sqrt(n)}).
#'     }}
#'
#'   \item{\code{dominance_table}}{Data frame. One row per ordered allele pair
#'     (A, B) per block per trait, for pairs where all three diplotype classes
#'     (AA, AB, BB) each have \code{>= min_n_diplotype} individuals.
#'     Contains 14 columns:
#'     \itemize{
#'       \item \code{block_id}, \code{CHR}, \code{start_bp}, \code{end_bp},
#'         \code{trait} - as above.
#'       \item \code{allele_A}, \code{allele_B} (character) - The two alleles
#'         being compared (alphabetically ordered: A comes before B).
#'       \item \code{mean_AA}, \code{mean_AB}, \code{mean_BB} (numeric) -
#'         Diplotype class means on the de-regressed phenotype scale.
#'       \item \code{a} (numeric) - Additive effect:
#'         \eqn{(\bar{y}_{BB} - \bar{y}_{AA}) / 2}. Positive means allele B
#'         increases trait value relative to A.
#'       \item \code{d} (numeric) - Dominance deviation:
#'         \eqn{\bar{y}_{AB} - (\bar{y}_{AA} + \bar{y}_{BB}) / 2}.
#'         Positive: heterozygote advantage (overdominance when \eqn{|d| > |a|}).
#'         Negative: heterozygote disadvantage (underdominance).
#'       \item \code{d_over_a} (numeric or \code{NA}) - Dominance ratio
#'         \eqn{d/a}. Interpretation: 0 = purely additive gene action;
#'         \eqn{\pm 0.5} = partial dominance; \eqn{\pm 1} = complete dominance
#'         (heterozygote equal to one homozygote); \eqn{|d/a| > 1} = overdominance
#'         (heterozygote exceeds both homozygotes - evidence for heterosis at
#'         this locus). \code{NA} when \eqn{|a| < 10^{-10}} (no additive
#'         variation between homozygotes).
#'       \item \code{overdominance} (logical) - \code{TRUE} when
#'         \eqn{|d/a| > 1} and \code{d_over_a} is not \code{NA}.
#'     }}
#'
#'   \item{\code{omnibus_tests}}{Data frame. One row per LD block per trait,
#'     for blocks where at least \code{min_n_diplotype} individuals carry at
#'     least 2 distinct diplotype classes. Contains 9 columns:
#'     \itemize{
#'       \item \code{block_id}, \code{trait} - as above.
#'       \item \code{n_diplotypes} (integer) - Number of diplotype classes
#'         with \code{>= min_n_diplotype} individuals, used as factor levels.
#'       \item \code{F_stat} (numeric) - F-statistic from one-way ANOVA
#'         (diplotype factor) on de-regressed phenotype residuals.
#'       \item \code{df1} (integer) - Numerator df = \code{n_diplotypes - 1}.
#'       \item \code{df2} (integer) - Denominator df = \code{n_obs - n_diplotypes}.
#'       \item \code{p_omnibus} (numeric, (0,1]) - Raw p-value.
#'       \item \code{p_omnibus_adj} (numeric, (0,1]) - Bonferroni-adjusted
#'         across all tested blocks per trait.
#'       \item \code{significant} (logical) - \code{TRUE} when
#'         \code{p_omnibus_adj < 0.05}.
#'     }
#'     Sorted ascending by \code{p_omnibus}.}
#' }
#'
#' @examples
#' \donttest{
#' data(ldx_geno, ldx_snp_info, ldx_blocks, ldx_blues, package = "LDxBlocks")
#' haps <- extract_haplotypes(ldx_geno, ldx_snp_info, ldx_blocks, min_snps = 5)
#' dip_res <- estimate_diplotype_effects(
#'   haplotypes = haps,
#'   blues      = setNames(ldx_blues$YLD, ldx_blues$id),
#'   blocks     = ldx_blocks,
#'   verbose    = FALSE
#' )
#' # Overdominant blocks
#' dip_res$dominance_table[dip_res$dominance_table$overdominance, ]
#' # Significant diplotype effects
#' dip_res$omnibus_tests[dip_res$omnibus_tests$significant, ]
#' }
#' @seealso \code{\link{infer_block_haplotypes}},
#'   \code{\link{test_block_haplotypes}}
#' @export
estimate_diplotype_effects <- function(
    haplotypes,
    blues,
    blocks,
    min_freq         = 0.05,
    min_n_diplotype  = 3L,
    id_col           = "id",
    blue_col         = "blue",
    blue_cols        = NULL,
    verbose          = TRUE
) {
  if (!requireNamespace("rrBLUP", quietly = TRUE))
    stop("rrBLUP is required: install.packages('rrBLUP')", call. = FALSE)

  .log <- function(...) if (verbose) message(sprintf("[diplotype] %s", paste0(...)))

  blues_list <- .parse_blues_assoc(blues, id_col, blue_col, blue_cols)
  traits     <- names(blues_list)
  bi         <- attr(haplotypes, "block_info")

  .log("Building haplotype GRM ...")
  hap_mat <- build_haplotype_feature_matrix(haplotypes, min_freq = min_freq,
                                            encoding = "additive_012")$matrix
  G <- compute_haplotype_grm(hap_mat)

  .log("Inferring diplotypes ...")
  # resolve_unphased=TRUE: heterozygous individuals get an arbitrary phase
  # assignment, which is acceptable for dominance decomposition since
  # mean(y_AA), mean(y_AB), mean(y_BB) depend only on diplotype identity,
  # not on which chromosome carries which haplotype.
  dip_tbl <- infer_block_haplotypes(haplotypes, resolve_unphased = TRUE)
  dip_tbl <- dip_tbl[!dip_tbl$missing, ]

  means_rows <- list(); dom_rows <- list(); omnibus_rows <- list()

  for (tr in traits) {
    .log("Trait: ", tr)
    pheno  <- blues_list[[tr]]
    common <- intersect(names(pheno), rownames(G))
    if (length(common) < 10L) next

    G_sub    <- G[common, common, drop = FALSE]
    y        <- pheno[common]
    null_fit <- tryCatch(rrBLUP::mixed.solve(y=y, K=G_sub, method="REML"),
                         error=function(e) NULL)
    if (is.null(null_fit)) next
    y_resid  <- y - as.numeric(null_fit$beta) - null_fit$u[common]

    for (bn in names(haplotypes)) {
      dip_bn  <- dip_tbl[dip_tbl$block_id == bn, ]
      if (!nrow(dip_bn)) next
      dip_ind <- dip_bn[dip_bn$id %in% common, ]
      if (nrow(dip_ind) < min_n_diplotype) next

      dip_tbl2    <- table(dip_ind$diplotype)
      keep_dip    <- names(dip_tbl2)[dip_tbl2 >= min_n_diplotype]
      dip_ind_sub <- dip_ind[dip_ind$diplotype %in% keep_dip, ]
      if (nrow(dip_ind_sub) < min_n_diplotype) next

      blk_meta <- if (!is.null(bi)) bi[bi$block_id == bn, , drop=FALSE] else NULL
      gm <- function(col, def=NA) {
        if (!is.null(blk_meta) && nrow(blk_meta)>0) blk_meta[[col]][1] else def
      }

      for (dp in keep_dip) {
        idx <- intersect(dip_ind_sub$id[dip_ind_sub$diplotype==dp], names(y_resid))
        if (length(idx) < 2L) next
        vals <- y_resid[idx]
        means_rows[[length(means_rows)+1L]] <- data.frame(
          block_id=bn, CHR=gm("CHR",NA_character_),
          start_bp=gm("start_bp",NA_integer_), end_bp=gm("end_bp",NA_integer_),
          trait=tr, diplotype=dp, n=length(vals),
          mean_blue=round(mean(vals),6),
          se_mean=round(stats::sd(vals)/sqrt(length(vals)),6),
          stringsAsFactors=FALSE)
      }

      # Build allele frequency table from resolved hap1/hap2 gamete strings.
      # Using dip_tbl (not raw haplotypes[[bn]]) ensures phased input works:
      # haplotypes[[bn]] contains 'A|B' strings for phased data, while
      # dip_tbl$hap1/hap2 already contain the individual gamete strings.
      gametes_bn <- c(dip_ind$hap1, dip_ind$hap2)
      gametes_bn <- gametes_bn[!is.na(gametes_bn) & !grepl("\\.", gametes_bn)]
      if (!length(gametes_bn)) next
      freq_tbl <- table(gametes_bn) / length(gametes_bn)
      com_al   <- names(freq_tbl)[freq_tbl >= min_freq]
      all_al   <- intersect(
        unique(unlist(lapply(keep_dip, function(dp) strsplit(dp,"/")[[1]]))),
        com_al)
      if (length(all_al) < 2L) next

      for (ai in seq_len(length(all_al)-1L)) {
        for (bi2 in seq(ai+1L, length(all_al))) {
          A <- all_al[ai]; B <- all_al[bi2]
          gm2 <- function(dp) {
            idx <- intersect(dip_ind_sub$id[dip_ind_sub$diplotype==dp],
                             names(y_resid))
            if (length(idx)<2L) NA_real_ else mean(y_resid[idx])
          }
          mAA <- gm2(paste(sort(c(A,A)),collapse="/"))
          mAB <- gm2(paste(sort(c(A,B)),collapse="/"))
          mBB <- gm2(paste(sort(c(B,B)),collapse="/"))
          if (any(is.na(c(mAA,mAB,mBB)))) next
          a_e  <- (mBB-mAA)/2
          d_e  <- mAB-(mAA+mBB)/2
          doa  <- if (abs(a_e)>1e-10) d_e/a_e else NA_real_
          dom_rows[[length(dom_rows)+1L]] <- data.frame(
            block_id=bn, CHR=gm("CHR",NA_character_),
            start_bp=gm("start_bp",NA_integer_), end_bp=gm("end_bp",NA_integer_),
            trait=tr, allele_A=A, allele_B=B,
            mean_AA=round(mAA,6), mean_AB=round(mAB,6), mean_BB=round(mBB,6),
            a=round(a_e,6), d=round(d_e,6),
            d_over_a=if(!is.na(doa)) round(doa,4) else NA_real_,
            overdominance=!is.na(doa) & abs(doa)>1,
            stringsAsFactors=FALSE)
        }
      }

      df_dip  <- data.frame(
        y_r=y_resid[dip_ind_sub$id], dip=factor(dip_ind_sub$diplotype))
      df_dip  <- df_dip[!is.na(df_dip$y_r),]
      if (nrow(df_dip) < min_n_diplotype+1L) next
      fit_dip <- tryCatch(stats::lm(y_r~dip,data=df_dip), error=function(e) NULL)
      if (is.null(fit_dip)) next
      an <- stats::anova(fit_dip)
      omnibus_rows[[length(omnibus_rows)+1L]] <- data.frame(
        block_id=bn, trait=tr, n_diplotypes=length(keep_dip),
        F_stat=round(an["dip","F value"],4),
        df1=an["dip","Df"], df2=an["Residuals","Df"],
        p_omnibus=an["dip","Pr(>F)"], stringsAsFactors=FALSE)
    }
  }

  means_df   <- if (length(means_rows))   do.call(rbind,means_rows)   else data.frame()
  dom_df     <- if (length(dom_rows))     do.call(rbind,dom_rows)     else data.frame()
  omnibus_df <- if (length(omnibus_rows)) {
    df <- do.call(rbind,omnibus_rows)
    df$p_omnibus_adj <- pmin(df$p_omnibus*nrow(df),1)
    df$significant   <- df$p_omnibus_adj < 0.05
    df[order(df$p_omnibus),]
  } else data.frame()

  .log("Done. Diplotype means:",nrow(means_df),"| Allele pairs:",nrow(dom_df),
       "| Blocks:",nrow(omnibus_df))
  structure(list(diplotype_means=means_df,dominance_table=dom_df,
                 omnibus_tests=omnibus_df),
            class=c("LDxBlocks_diplotype","list"))
}

#' @export
print.LDxBlocks_diplotype <- function(x, ...) {
  cat("LDxBlocks Diplotype Effect Results\n")
  cat("  Diplotype means rows:", nrow(x$diplotype_means), "\n")
  cat("  Allele pair comparisons:", nrow(x$dominance_table), "\n")
  if (nrow(x$dominance_table)>0)
    cat("  Overdominant pairs (|d/a|>1):",
        sum(x$dominance_table$overdominance,na.rm=TRUE),"\n")
  if (nrow(x$omnibus_tests)>0)
    cat("  Significant blocks:",sum(x$omnibus_tests$significant,na.rm=TRUE),"\n")
  invisible(x)
}
