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
# simpleM helpers (Gao et al. 2008, 2010, 2011)
# .MeFF_PCA          : effective-test count from eigenvalue spectrum
# .compute_simplem_*  : Meff from a dosage matrix, optionally by group
# .make_block_summary_matrix : one PC1 per block for block-level Meff
# .adjust_p_simplem_* : Bonferroni-style and Sidak-style p-value adjustment
# ==============================================================================

.MeFF_PCA <- function(eigenValues, percentCut = 0.995) {
  eigenValues <- as.numeric(eigenValues)
  eigenValues <- eigenValues[is.finite(eigenValues) & eigenValues > 0]
  if (!length(eigenValues)) return(0L)
  totalEigenValues <- sum(eigenValues)
  myCut            <- percentCut * totalEigenValues
  myEigenSum       <- 0
  index_Eigen      <- 0L
  for (i in seq_along(eigenValues)) {
    if (myEigenSum <= myCut) {
      myEigenSum  <- myEigenSum + eigenValues[i]
      index_Eigen <- i
    } else break
  }
  as.integer(index_Eigen)
}

.compute_simplem_meff_matrix <- function(X,
                                         percent_cut = 0.995,
                                         max_cols    = 1000L,
                                         cor_use     = "pairwise.complete.obs") {
  X    <- as.matrix(X)
  if (!nrow(X) || !ncol(X)) return(list(meff = 0L, n_tests = 0L))
  keep <- apply(X, 2L, function(z) {
    z <- z[is.finite(z)]; length(z) > 1L && stats::sd(z) > 0
  })
  X <- X[, keep, drop = FALSE]
  if (!ncol(X)) return(list(meff = 0L, n_tests = 0L))

  if (ncol(X) <= max_cols) {
    R <- suppressWarnings(stats::cor(X, use = cor_use))
    R[!is.finite(R)] <- 0; diag(R) <- 1
    eig  <- eigen(R, symmetric = TRUE, only.values = TRUE)$values
    meff <- .MeFF_PCA(abs(eig), percentCut = percent_cut)
    return(list(meff = as.integer(meff), n_tests = ncol(X)))
  }

  # chunk-and-sum for large matrices
  idx      <- split(seq_len(ncol(X)), ceiling(seq_len(ncol(X)) / max_cols))
  meff_sum <- 0L
  for (ii in idx) {
    R <- suppressWarnings(stats::cor(X[, ii, drop = FALSE], use = cor_use))
    R[!is.finite(R)] <- 0; diag(R) <- 1
    eig      <- eigen(R, symmetric = TRUE, only.values = TRUE)$values
    meff_sum <- meff_sum + .MeFF_PCA(abs(eig), percentCut = percent_cut)
  }
  list(meff = as.integer(meff_sum), n_tests = ncol(X))
}

.compute_simplem_meff_by_group <- function(X,
                                           group,
                                           percent_cut = 0.995,
                                           max_cols    = 1000L) {
  X <- as.matrix(X)
  if (!ncol(X)) return(integer(0))
  if (length(group) != ncol(X))
    stop("group length must match ncol(X).", call. = FALSE)
  split_idx <- split(seq_len(ncol(X)), as.character(group))
  vapply(split_idx, function(ii) {
    .compute_simplem_meff_matrix(
      X[, ii, drop = FALSE],
      percent_cut = percent_cut,
      max_cols    = max_cols
    )$meff
  }, integer(1L))
}

.make_block_summary_matrix <- function(X, block_ids) {
  X <- as.matrix(X)
  if (!ncol(X)) return(matrix(numeric(0), nrow = nrow(X), ncol = 0L))
  ublocks <- unique(block_ids)
  out <- matrix(NA_real_, nrow = nrow(X), ncol = length(ublocks))
  colnames(out) <- ublocks; rownames(out) <- rownames(X)
  for (i in seq_along(ublocks)) {
    ii <- which(block_ids == ublocks[i])
    Xi <- X[, ii, drop = FALSE]
    if (!ncol(Xi)) next
    if (ncol(Xi) == 1L) {
      out[, i] <- Xi[, 1L]
    } else {
      Xi <- scale(Xi, center = TRUE, scale = FALSE)
      Xi[!is.finite(Xi)] <- 0
      pc1 <- tryCatch(
        stats::prcomp(Xi, center = FALSE, scale. = FALSE)$x[, 1L],
        error = function(e) rowMeans(Xi, na.rm = TRUE)
      )
      out[, i] <- pc1
    }
  }
  out
}

.adjust_p_simplem_bonf <- function(p, meff) {
  out <- rep(NA_real_, length(p))
  ok  <- is.finite(p) & is.finite(meff) & meff > 0
  out[ok] <- pmin(p[ok] * meff[ok], 1)
  out
}

.adjust_p_simplem_sidak <- function(p, meff) {
  out <- rep(NA_real_, length(p))
  ok  <- is.finite(p) & is.finite(meff) & meff > 0
  out[ok] <- pmin(1 - (1 - p[ok])^meff[ok], 1)
  out
}

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

#' Block-Level Haplotype Association Testing (Q+K Mixed Linear Model
#' with simpleM Multiple-Testing Correction)
#'
#' @description
#' Performs genome-wide haplotype block association tests for one or more
#' quantitative traits. Each LD block is tested as a unit: per-allele Wald
#' tests identify which specific haplotype alleles drive association, and an
#' omnibus F-test evaluates the block as a whole. Population structure and
#' kinship are corrected jointly through a unified mixed linear model.
#'
#' In addition to raw and FDR-adjusted p-values, this function applies
#' \strong{simpleM}-style multiple-testing correction to account for
#' correlation among tested haplotype alleles. simpleM estimates the effective
#' number of independent tests (\eqn{M_{\mathrm{eff}}}) from the eigenvalues
#' of the haplotype allele correlation matrix and derives either a
#' Bonferroni-style or Sidak-style adjusted threshold.
#'
#' \strong{Statistical model (Q+K / EMMAX formulation):}
#'
#' \deqn{y = \mu + \alpha \cdot x_{\mathrm{hap}} +
#'        \sum_{k=1}^{K} \beta_k \, PC_k + g + \varepsilon}
#'
#' \describe{
#'   \item{\eqn{x_{\mathrm{hap}}}}{Haplotype allele dosage (0, 1, or 2 copies
#'     for phased data; 0 or 1 for unphased).}
#'   \item{\eqn{PC_k}}{k-th eigenvector of the haplotype GRM, as a
#'     fixed-effect covariate for population structure.}
#'   \item{\eqn{g \sim MVN(0,\,\sigma_g^2 G)}}{Polygenic background as a
#'     random effect.}
#'   \item{\eqn{\varepsilon \sim MVN(0,\,\sigma_e^2 I)}}{Residual error.}
#' }
#'
#' \strong{simpleM extension (Gao et al. 2008, 2010, 2011):}
#' \eqn{M_{\mathrm{eff}}} is estimated from the eigenspectrum of the
#' haplotype allele dosage correlation matrix - the number of principal
#' components needed to explain \code{meff_percent_cut} of variance (default
#' 99.5\%). Two simpleM adjusted p-value paths are always computed:
#' \itemize{
#'   \item \code{p_simplem} - Bonferroni-style:
#'     \eqn{\min(p \times M_{\mathrm{eff}},\,1)}.
#'   \item \code{p_simplem_sidak} - Sidak-style:
#'     \eqn{1 - (1 - p)^{M_{\mathrm{eff}}}}.
#' }
#' Because haplotype alleles are correlated within and across nearby blocks,
#' simpleM is less conservative than raw Bonferroni while still providing
#' family-wise error control.
#'
#' @param haplotypes Named list produced by \code{\link{extract_haplotypes}}.
#'   Each element is a named character vector (one haplotype dosage string per
#'   individual), and the list must carry a \code{block_info} attribute.
#'
#' @param blues Pre-adjusted phenotype means. Accepted in four formats:
#'   named numeric vector, single-trait data frame, multi-trait data frame,
#'   or named list of named numeric vectors. At least 10 common individuals
#'   are required per trait.
#'
#' @param blocks LD block table from \code{\link{run_Big_LD_all_chr}}.
#'   Used for block metadata annotation in the output only.
#'
#' @param n_pcs Integer (\code{>= 0}) or \code{NULL}. Number of haplotype-GRM
#'   eigenvectors to include as fixed-effect population structure covariates.
#'   \code{0L} (default): pure GRM / EMMAX. \code{1-10}: Q+K model.
#'   \code{NULL}: auto-selected from the GRM scree-plot elbow, capped at 10.
#'
#' @param top_n Integer or \code{NULL}. Maximum haplotype alleles per block
#'   in the feature matrix. \code{NULL} retains all alleles above
#'   \code{min_freq}.
#'
#' @param min_freq Numeric in (0, 1). Minimum haplotype allele frequency.
#'   Default \code{0.05}.
#'
#' @param id_col Character. Individual-ID column name when \code{blues} is a
#'   data frame. Default \code{"id"}.
#'
#' @param blue_col Character. Phenotype column name for single-trait data
#'   frames. Default \code{"blue"}.
#'
#' @param blue_cols Character vector. Phenotype column names for multi-trait
#'   data frames. Default \code{NULL}.
#'
#' @param sig_threshold Numeric in (0, 1]. Significance cutoff applied to the
#'   p-value chosen by \code{sig_metric}. Default \code{0.05}. Common choices:
#'   \itemize{
#'     \item \code{0.05} - recommended with \code{sig_metric = "p_fdr"},
#'       \code{"p_simplem"}, or \code{"p_simplem_sidak"}.
#'     \item \code{0.05 / n_tests} - raw Bonferroni with
#'       \code{sig_metric = "p_wald"}.
#'   }
#'
#' @param sig_metric Character. Which p-value to use for the
#'   \code{significant} and \code{significant_omnibus} flags. One of:
#'   \itemize{
#'     \item \code{"p_wald"} - raw Wald p-value. Use with a pre-corrected
#'       \code{sig_threshold} (e.g. \code{0.05 / n_tests}).
#'     \item \code{"p_fdr"} - Benjamini-Hochberg FDR. Recommended for
#'       discovery-oriented analyses.
#'     \item \code{"p_simplem"} - simpleM Bonferroni-style adjusted p-value.
#'     \item \code{"p_simplem_sidak"} - simpleM Sidak-style adjusted p-value.
#'       Recommended default when family-wise error control is desired with
#'       correlated haplotype predictors.
#'   }
#'   All four p-value columns are always present in the output regardless of
#'   this choice.
#'
#' @param meff_scope Character. Scope at which \eqn{M_{\mathrm{eff}}} is
#'   estimated. One of:
#'   \itemize{
#'     \item \code{"chromosome"} (recommended) - separate \eqn{M_{\mathrm{eff}}}
#'       per chromosome; best balance of scalability and LD awareness.
#'     \item \code{"global"} - one genome-wide \eqn{M_{\mathrm{eff}}}; suitable
#'       for moderate-sized scans.
#'     \item \code{"block"} - one \eqn{M_{\mathrm{eff}}} per LD block at the
#'       allele level; for block omnibus tests \eqn{M_{\mathrm{eff}} = 1}.
#'   }
#'
#' @param meff_percent_cut Numeric in (0, 1). Proportion of variance explained
#'   by the retained PCs in the simpleM eigendecomposition. Default \code{0.995}
#'   (99.5\%), following the original simpleM recommendation.
#'
#' @param meff_max_cols Integer. Maximum columns per eigendecomposition block
#'   when computing \eqn{M_{\mathrm{eff}}}. Larger groups are chunked and
#'   summed. Default \code{1000L}.
#'
#' @param plot Logical. If \code{TRUE}, save Manhattan and Q-Q plots.
#'   Default \code{FALSE}. Significance markers reflect \code{sig_metric}.
#'
#' @param out_dir Character. Directory for optional plots. Default \code{"."}.
#'
#' @param verbose Logical. Print progress messages. Default \code{TRUE}.
#'
#' @return A named list of class \code{c("LDxBlocks_haplotype_assoc", "list")}:
#'
#' \describe{
#'   \item{\code{allele_tests}}{Data frame with one row per allele per block
#'     per trait. Always includes: \code{p_wald}, \code{p_fdr},
#'     \code{Meff}, \code{alpha_simplem}, \code{alpha_simplem_sidak},
#'     \code{p_simplem}, \code{p_simplem_sidak}, \code{significant}.}
#'   \item{\code{block_tests}}{Data frame with one row per block per trait.
#'     Always includes: \code{p_omnibus}, \code{p_omnibus_fdr},
#'     \code{p_omnibus_adj} (backward-compat Bonferroni),
#'     \code{Meff}, \code{alpha_simplem}, \code{alpha_simplem_sidak},
#'     \code{p_omnibus_simplem}, \code{p_omnibus_simplem_sidak},
#'     \code{significant_omnibus}.}
#'   \item{\code{traits}}{Character vector of trait names tested.}
#'   \item{\code{n_pcs_used}}{Integer. GRM PCs used in the null model.}
#'   \item{\code{sig_threshold}}{Numeric. Significance cutoff used.}
#'   \item{\code{sig_metric}}{Character. P-value used for significance flags.}
#'   \item{\code{n_tests}}{Integer. Total allele-level tests performed.}
#'   \item{\code{meff_scope}}{Character. simpleM estimation scope.}
#'   \item{\code{meff_percent_cut}}{Numeric. Variance cutoff for simpleM.}
#'   \item{\code{meff}}{Named list of \eqn{M_{\mathrm{eff}}} summaries per
#'     trait, each with \code{$allele} (global/chromosome/block) and
#'     \code{$block} (global/chromosome) components.}
#' }
#'
#' @examples
#' \donttest{
#' data(ldx_geno, ldx_snp_info, ldx_blocks, ldx_blues, package = "LDxBlocks")
#' haps <- extract_haplotypes(ldx_geno, ldx_snp_info, ldx_blocks)
#'
#' # FDR-based discovery (default EMMAX model)
#' assoc_fdr <- test_block_haplotypes(
#'   haplotypes = haps, blues = setNames(ldx_blues$YLD, ldx_blues$id),
#'   blocks = ldx_blocks, sig_metric = "p_fdr", verbose = FALSE
#' )
#'
#' # simpleM Sidak with 3 GRM PCs (Q+K), chromosome-wise Meff
#' assoc_sm <- test_block_haplotypes(
#'   haplotypes = haps, blues = setNames(ldx_blues$YLD, ldx_blues$id),
#'   blocks = ldx_blocks, n_pcs = 3L,
#'   sig_metric = "p_simplem_sidak", meff_scope = "chromosome",
#'   verbose = FALSE
#' )
#'
#' # Check per-trait Meff
#' assoc_sm$meff$trait$allele$chromosome
#'
#' # Multi-trait with simpleM
#' assoc_mt <- test_block_haplotypes(
#'   haplotypes = haps, blues = ldx_blues,
#'   blocks = ldx_blocks, id_col = "id", blue_cols = c("YLD","RES"),
#'   sig_metric = "p_simplem", meff_scope = "chromosome", verbose = FALSE
#' )
#' }
#'
#' @references
#' Endelman JB (2011). Ridge regression and other kernels for genomic
#' selection with R package rrBLUP. \emph{Plant Genome} \strong{4}:250-255.
#' \doi{10.3835/plantgenome2011.08.0024}
#'
#' Yu J et al. (2006). A unified mixed-model method for association mapping.
#' \emph{Nature Genetics} \strong{38}(2):203-208. \doi{10.1038/ng1702}
#'
#' Kang HM et al. (2010). Variance component model to account for sample
#' structure in genome-wide association studies.
#' \emph{Nature Genetics} \strong{42}(4):348-354. \doi{10.1038/ng.548}
#'
#' Gao X et al. (2008). A multiple testing correction method for genetic
#' association studies using correlated single nucleotide polymorphisms.
#' \emph{Genetic Epidemiology} \strong{32}:361-369.
#' \doi{10.1002/gepi.20310}
#'
#' Gao X et al. (2010). Avoiding the high Bonferroni penalty in genome-wide
#' association studies. \emph{Genetic Epidemiology} \strong{34}:100-105.
#' \doi{10.1002/gepi.20430}
#'
#' Gao X (2011). Multiple testing corrections for imputed SNPs.
#' \emph{Genetic Epidemiology} \strong{35}:154-158.
#' \doi{10.1002/gepi.20563}
#'
#' @seealso \code{\link{extract_haplotypes}},
#'   \code{\link{compute_haplotype_grm}},
#'   \code{\link{run_haplotype_prediction}},
#'   \code{\link{estimate_diplotype_effects}}
#' @export

test_block_haplotypes <- function(
    haplotypes,
    blues,
    blocks,
    n_pcs            = 0L,
    top_n            = NULL,
    min_freq         = 0.05,
    id_col           = "id",
    blue_col         = "blue",
    blue_cols        = NULL,
    sig_threshold    = 0.05,
    sig_metric       = c("p_wald", "p_fdr", "p_simplem", "p_simplem_sidak"),
    meff_scope       = c("chromosome", "global", "block"),
    meff_percent_cut = 0.995,
    meff_max_cols    = 1000L,
    plot             = FALSE,
    out_dir          = ".",
    verbose          = TRUE
) {
  if (!requireNamespace("rrBLUP", quietly = TRUE))
    stop("rrBLUP is required: install.packages('rrBLUP')", call. = FALSE)

  sig_metric <- match.arg(sig_metric)
  meff_scope <- match.arg(meff_scope)
  .log <- function(...) if (verbose) message(sprintf("[assoc] %s", paste0(...)))

  blues_list <- .parse_blues_assoc(blues, id_col, blue_col, blue_cols)
  traits     <- names(blues_list)
  bi         <- attr(haplotypes, "block_info")

  # -- Build haplotype feature matrix and GRM (once, shared across all traits)
  .log("Building haplotype feature matrix ...")
  hap_mat <- build_haplotype_feature_matrix(
    haplotypes, top_n = top_n, min_freq = min_freq, encoding = "additive_012"
  )$matrix
  # Column names are "blockID_hapN" - extract block ID prefix for grouping
  col_blocks <- sub("_hap\\d+$", "", colnames(hap_mat))

  # Detect phasing from block_info attribute set by extract_haplotypes().
  # This is unambiguous: the 'phased' column is TRUE for phased VCF input,
  # FALSE for unphased dosage matrix input.
  # Do NOT detect from data values (col_max trick) - a phased allele with
  # no homozygous carriers has max dosage = 1 even though the scale is 0/1/2.
  .bi_ph  <- attr(haplotypes, "block_info")
  is_phased_data <- isTRUE(.bi_ph$phased[1L]) ||
    (is.data.frame(.bi_ph) && nrow(.bi_ph) > 0L &&
       any(isTRUE(.bi_ph$phased)))
  # Scale factor: 2 for phased (0/1/2 dosage), 1 for unphased (0/1 presence)
  dose_scale <- if (is_phased_data) 2 else 1

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
  meff_cache  <- vector("list", length(traits))
  names(meff_cache) <- traits

  # block_id -> CHR lookup for chromosome-wise Meff grouping
  block_chr_map <- if (!is.null(bi) && nrow(bi) > 0L) {
    stats::setNames(as.character(bi$CHR), bi$block_id)
  } else {
    NULL
  }

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

    X_valid      <- X_all[, valid_cols, drop = FALSE]
    XtX_valid    <- XtX_vec[valid_cols]           # x_j'x_j (centred)
    col_bl_valid <- col_blocks[valid_cols]
    # Per-column CHR label for chromosome-wise Meff grouping
    col_chr_valid <- if (!is.null(block_chr_map)) {
      unname(block_chr_map[col_bl_valid])
    } else {
      rep(NA_character_, length(col_bl_valid))
    }

    # Cache the exact tested matrix + grouping vectors for post-scan Meff
    meff_cache[[tr]] <- list(
      X_valid      = X_valid,
      col_blocks   = col_bl_valid,
      col_chr      = col_chr_valid
    )

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
    # Allele frequency among phenotyped individuals actually used in the test.
    # dose_scale is set once above from block_info$phased:
    #   dose_scale = 2: phased data, dosage values are 0/1/2
    #                   freq = colMeans(X) / 2
    #   dose_scale = 1: unphased data, dosage values are 0/1 (presence/absence)
    #                   freq = colMeans(X)
    # Using block_info$phased is reliable even for rare alleles with no
    # homozygous carriers (where col_max would be 1 even in phased data).
    freq_all <- colMeans(X_valid, na.rm = TRUE) / dose_scale

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
      allele_label <- sub(paste0("^", bn, "_hap\\d+_?"), "",
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
        allele_freq_tested = round(freq_all[ci], 4),
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
      list(allele_tests     = data.frame(), block_tests = data.frame(),
           traits           = traits,       n_pcs_used  = n_pcs_used,
           sig_threshold    = sig_threshold, sig_metric = sig_metric,
           n_tests          = 0L,
           meff_scope       = meff_scope, meff_percent_cut = meff_percent_cut,
           meff             = list()),
      class = c("LDxBlocks_haplotype_assoc", "list")
    ))
  }

  allele_df    <- do.call(rbind, allele_rows)
  block_df     <- do.call(rbind, block_rows)
  n_tests      <- nrow(allele_df)
  meff_summary <- list()

  # -- FDR adjustment (always computed, all traits at once) -------------------
  if (nrow(allele_df) > 0)
    allele_df$p_fdr <- stats::p.adjust(allele_df$p_wald, method = "BH")
  if (nrow(block_df) > 0) {
    block_df$p_omnibus_fdr <- stats::p.adjust(block_df$p_omnibus, method = "BH")
    # Bonferroni-style block column kept for backward compatibility
    block_df$p_omnibus_adj <- NA_real_
    for (tr in unique(block_df$trait)) {
      idx_tr <- which(block_df$trait == tr)
      block_df$p_omnibus_adj[idx_tr] <- pmin(block_df$p_omnibus[idx_tr] * length(idx_tr), 1)
    }
  }

  # -- simpleM Meff estimation and adjustment (per trait) --------------------
  # Initialise simpleM columns (NA until filled below)
  if (nrow(allele_df) > 0) {
    allele_df$Meff               <- NA_real_
    allele_df$alpha_simplem       <- NA_real_
    allele_df$alpha_simplem_sidak <- NA_real_
    allele_df$p_simplem           <- NA_real_
    allele_df$p_simplem_sidak     <- NA_real_
  }
  if (nrow(block_df) > 0) {
    block_df$Meff                    <- NA_real_
    block_df$alpha_simplem            <- NA_real_
    block_df$alpha_simplem_sidak      <- NA_real_
    block_df$p_omnibus_simplem        <- NA_real_
    block_df$p_omnibus_simplem_sidak  <- NA_real_
  }

  for (tr in traits) {
    cache <- meff_cache[[tr]]
    if (is.null(cache)) next

    X_valid    <- cache$X_valid
    col_bl     <- cache$col_blocks
    col_chr    <- cache$col_chr

    # -- Allele-level Meff estimates ------------------------------------------
    meff_global_a <- .compute_simplem_meff_matrix(
      X_valid, percent_cut = meff_percent_cut, max_cols = meff_max_cols)$meff

    meff_chr_a <- if (length(col_chr) && !all(is.na(col_chr))) {
      .compute_simplem_meff_by_group(
        X_valid, group = col_chr,
        percent_cut = meff_percent_cut, max_cols = meff_max_cols)
    } else integer(0)

    meff_block_a <- .compute_simplem_meff_by_group(
      X_valid, group = col_bl,
      percent_cut = meff_percent_cut, max_cols = meff_max_cols)

    # -- Block-level Meff: one PC1 summary per block --------------------------
    X_bsum <- .make_block_summary_matrix(X_valid, col_bl)
    bids_sum <- colnames(X_bsum)
    bchr_sum <- if (!is.null(block_chr_map)) unname(block_chr_map[bids_sum]) else character(0)

    meff_global_b <- if (ncol(X_bsum)) {
      .compute_simplem_meff_matrix(
        X_bsum, percent_cut = meff_percent_cut, max_cols = meff_max_cols)$meff
    } else 0L

    meff_chr_b <- if (length(bchr_sum) && !all(is.na(bchr_sum))) {
      .compute_simplem_meff_by_group(
        X_bsum, group = bchr_sum,
        percent_cut = meff_percent_cut, max_cols = meff_max_cols)
    } else integer(0)

    meff_summary[[tr]] <- list(
      allele = list(global = meff_global_a, chromosome = meff_chr_a, block = meff_block_a),
      block  = list(global = meff_global_b, chromosome = meff_chr_b)
    )

    # -- Assign Meff to allele_df rows for this trait -------------------------
    if (nrow(allele_df) > 0) {
      idx_tr <- which(allele_df$trait == tr)
      meff_a <- switch(meff_scope,
                       global     = rep(meff_global_a, length(idx_tr)),
                       chromosome = if (length(meff_chr_a))
                         unname(meff_chr_a[as.character(allele_df$CHR[idx_tr])]) else meff_global_a,
                       block      = unname(meff_block_a[as.character(allele_df$block_id[idx_tr])])
      )
      allele_df$Meff[idx_tr]               <- meff_a
      allele_df$alpha_simplem[idx_tr]       <- sig_threshold / meff_a
      allele_df$alpha_simplem_sidak[idx_tr] <- 1 - (1 - sig_threshold)^(1 / meff_a)
      allele_df$p_simplem[idx_tr] <-
        .adjust_p_simplem_bonf(allele_df$p_wald[idx_tr], meff_a)
      allele_df$p_simplem_sidak[idx_tr] <-
        .adjust_p_simplem_sidak(allele_df$p_wald[idx_tr], meff_a)
    }

    # -- Assign Meff to block_df rows for this trait --------------------------
    if (nrow(block_df) > 0) {
      idx_tb <- which(block_df$trait == tr)
      meff_b <- switch(meff_scope,
                       global     = rep(meff_global_b, length(idx_tb)),
                       chromosome = if (length(meff_chr_b))
                         unname(meff_chr_b[as.character(block_df$CHR[idx_tb])]) else meff_global_b,
                       block      = rep(1L, length(idx_tb))  # one omnibus test per block -> Meff=1
      )
      meff_b[is.na(meff_b) | meff_b < 1L] <- 1L
      block_df$Meff[idx_tb]                    <- meff_b
      block_df$alpha_simplem[idx_tb]            <- sig_threshold / meff_b
      block_df$alpha_simplem_sidak[idx_tb]      <- 1 - (1 - sig_threshold)^(1 / meff_b)
      block_df$p_omnibus_simplem[idx_tb] <-
        .adjust_p_simplem_bonf(block_df$p_omnibus[idx_tb], meff_b)
      block_df$p_omnibus_simplem_sidak[idx_tb] <-
        .adjust_p_simplem_sidak(block_df$p_omnibus[idx_tb], meff_b)
    }
  }  # end per-trait Meff loop

  # -- Significance flags using chosen sig_metric ----------------------------
  if (nrow(allele_df) > 0) {
    allele_df$significant <- switch(sig_metric,
                                    p_wald          = allele_df$p_wald          <= sig_threshold,
                                    p_fdr           = allele_df$p_fdr            <= sig_threshold,
                                    p_simplem       = allele_df$p_simplem        <= sig_threshold,
                                    p_simplem_sidak = allele_df$p_simplem_sidak  <= sig_threshold
    )
    allele_df <- allele_df[order(allele_df$CHR, allele_df$start_bp, allele_df$p_wald), ]
    rownames(allele_df) <- NULL
  }
  if (nrow(block_df) > 0) {
    block_df$significant_omnibus <- switch(sig_metric,
                                           p_wald          = block_df$p_omnibus                 <= sig_threshold,
                                           p_fdr           = block_df$p_omnibus_fdr             <= sig_threshold,
                                           p_simplem       = block_df$p_omnibus_simplem         <= sig_threshold,
                                           p_simplem_sidak = block_df$p_omnibus_simplem_sidak   <= sig_threshold
    )
    block_df <- block_df[order(block_df$CHR, block_df$start_bp, block_df$p_omnibus), ]
    rownames(block_df) <- NULL
  }

  n_sig <- sum(allele_df$significant, na.rm = TRUE)
  .log("Done. Allele tests: ", n_tests,
       " | Significant (", sig_metric, " <= ", sig_threshold, "): ", n_sig)

  # -- Optional Manhattan and QQ plots ----------------------------------------
  if (isTRUE(plot) && nrow(allele_df) > 0) {
    .plot_assoc_results(allele_df, sig_threshold = sig_threshold,
                        out_dir = out_dir, verbose = verbose)
  }

  structure(
    list(
      allele_tests     = allele_df,
      block_tests      = block_df,
      traits           = traits,
      n_pcs_used       = n_pcs_used,
      sig_threshold    = sig_threshold,
      sig_metric       = sig_metric,
      n_tests          = n_tests,
      meff_scope       = meff_scope,
      meff_percent_cut = meff_percent_cut,
      meff             = meff_summary
    ),
    class = c("LDxBlocks_haplotype_assoc", "list")
  )
}


# ==============================================================================
# Internal: Manhattan and QQ plots for test_block_haplotypes results
# Style: discrete chromosome x-axis, jitter, per-chromosome colour,
#        shape legend for significant multi-trait hits.
# ==============================================================================

.plot_assoc_results <- function(allele_df, sig_threshold = 0.05,
                                out_dir = ".", verbose = TRUE) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    if (verbose) message("[plot] ggplot2 not installed. Install with: install.packages('ggplot2')")
    return(invisible(NULL))
  }
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  traits <- unique(allele_df$trait)

  # -- Shared chromosome ordering ----------------------------------------------
  all_chrs <- unique(allele_df$CHR)
  chr_num  <- suppressWarnings(as.integer(all_chrs))
  chr_order <- if (!any(is.na(chr_num))) {
    all_chrs[order(chr_num)]
  } else {
    sort(all_chrs)
  }
  chrom_colors <- stats::setNames(
    rep(c("#6495ED","#800000","#f0b27a","#008080","#8e44ad"),
        length.out = length(chr_order)),
    chr_order
  )

  # -- Per-trait Manhattan and QQ -----------------------------------------------
  for (tr in traits) {
    df <- allele_df[allele_df$trait == tr & !is.na(allele_df$p_wald) &
                      allele_df$p_wald > 0, ]
    if (!nrow(df)) next

    df$CHR <- factor(as.character(df$CHR), levels = chr_order)
    df$log10p      <- -log10(df$p_wald)
    df$Significant <- ifelse(df$significant, "Above threshold", "Below threshold")
    sig_line  <- -log10(sig_threshold)
    sugg_line <- max(0, sig_line - 1)  # one order of magnitude less stringent

    # -- Manhattan plot (your style) --------------------------------------------
    p_manh <- ggplot2::ggplot(
      df, ggplot2::aes(x = CHR, y = log10p, colour = CHR)) +
      ggplot2::geom_jitter(
        data    = df[df$Significant == "Below threshold", ],
        alpha   = 0.5, width = 0.42, size = 0.8) +
      ggplot2::geom_jitter(
        data    = df[df$Significant == "Above threshold", ],
        ggplot2::aes(shape = allele),
        alpha   = 0.85, width = 0.42, size = 2.2) +
      ggplot2::scale_colour_manual(values = chrom_colors, guide = "none") +
      ggplot2::scale_shape_manual(
        name   = "Haplotype allele",
        values = seq(0, max(0, length(unique(df$allele[df$significant])) - 1))) +
      ggplot2::geom_hline(yintercept = sig_line,  colour = "red",
                          linetype = "solid",  linewidth = 0.4) +
      ggplot2::geom_hline(yintercept = sugg_line, colour = "blue",
                          linetype = "dotted", linewidth = 0.3) +
      ggplot2::scale_x_discrete(expand = ggplot2::expansion(add = 0.5)) +
      ggplot2::scale_y_continuous(
        expand = ggplot2::expansion(0),
        limits = c(0, max(df$log10p, sig_line + 1) * 1.05),
        breaks = scales::pretty_breaks(n = 6)
      ) +
      ggplot2::labs(
        title    = paste0("Manhattan plot \u2014 ", tr),
        subtitle = paste0(nrow(df), " haplotype allele tests | ",
                          "threshold p = ", sig_threshold, " | ",
                          sum(df$significant, na.rm = TRUE), " significant"),
        x = "Chromosome",
        y = expression(-log[10](italic(p)))
      ) +
      ggplot2::theme_classic(base_size = 11) +
      ggplot2::guides(shape = ggplot2::guide_legend(nrow = 2, override.aes = list(size = 3))) +
      ggplot2::theme(
        plot.title        = ggplot2::element_text(face = "bold"),
        plot.subtitle     = ggplot2::element_text(size = 9, colour = "grey40"),
        legend.title      = ggplot2::element_text(size = 9, face = "bold"),
        legend.position   = "bottom",
        axis.title.x      = ggplot2::element_blank(),
        axis.text.x       = ggplot2::element_text(size = 8, face = "bold"),
        axis.ticks.x      = ggplot2::element_blank(),
        axis.line.x       = ggplot2::element_blank(),
        axis.title.y      = ggplot2::element_text(size = 12, face = "bold"),
        plot.margin       = ggplot2::unit(c(0.2, 0.2, 0.1, 0.2), "cm"),
        panel.grid.major  = ggplot2::element_blank(),
        panel.grid.minor  = ggplot2::element_blank()
      )

    # -- QQ plot ----------------------------------------------------------------
    n_p  <- nrow(df)
    obs  <- sort(df$log10p, decreasing = TRUE)
    exp  <- -log10((seq_len(n_p)) / (n_p + 1))
    lambda <- stats::median(
      stats::qchisq(df$p_wald, df = 1, lower.tail = FALSE), na.rm = TRUE
    ) / stats::qchisq(0.5, df = 1, lower.tail = FALSE)

    qq_df <- data.frame(expected = rev(exp), observed = rev(obs),
                        significant = rev(sort(df$significant)))

    p_qq <- ggplot2::ggplot(qq_df, ggplot2::aes(x = expected, y = observed)) +
      ggplot2::geom_abline(intercept = 0, slope = 1,
                           colour = "red", linetype = "solid", linewidth = 0.5) +
      ggplot2::geom_point(
        data  = qq_df[!qq_df$significant, ],
        colour = "#6495ED", size = 0.7, alpha = 0.6) +
      ggplot2::geom_point(
        data  = qq_df[qq_df$significant, ],
        colour = "#800000", size = 1.5, alpha = 0.9) +
      ggplot2::labs(
        title    = paste0("QQ plot \u2014 ", tr),
        subtitle = sprintf("\u03bb (GC) = %.3f | %d tests", lambda, n_p),
        x = expression("Expected " * -log[10](italic(p))),
        y = expression("Observed " * -log[10](italic(p)))
      ) +
      ggplot2::theme_classic(base_size = 11) +
      ggplot2::theme(
        plot.title     = ggplot2::element_text(face = "bold"),
        plot.subtitle  = ggplot2::element_text(size = 9, colour = "grey40"),
        axis.title     = ggplot2::element_text(size = 11, face = "bold"),
        panel.grid.major = ggplot2::element_blank(),
        panel.grid.minor = ggplot2::element_blank()
      )

    # -- Save -------------------------------------------------------------------
    safe_tr   <- gsub("[^A-Za-z0-9_]", "_", tr)
    manh_file <- file.path(out_dir, paste0("manhattan_", safe_tr, ".png"))
    qq_file   <- file.path(out_dir, paste0("qq_",        safe_tr, ".png"))

    ggplot2::ggsave(manh_file, p_manh, width = 14, height = 5.5, dpi = 300)
    ggplot2::ggsave(qq_file,   p_qq,   width  =  5, height = 5,   dpi = 300)

    if (verbose) {
      message("[plot] Manhattan saved: ", manh_file)
      message("[plot] QQ saved:        ", qq_file)
      message(sprintf("[plot] lambda (GC) for %s: %.3f", tr, lambda))
    }
  }
  invisible(NULL)
}


#' @export
print.LDxBlocks_haplotype_assoc <- function(x, ...) {
  cat("LDxBlocks Haplotype Association Results\n")
  cat("  Model:        GRM")
  if (x$n_pcs_used > 0L)
    cat(" + ", x$n_pcs_used, " GRM-derived PC(s) [Q+K]", sep = "")
  else
    cat(" [EMMAX/P3D]")
  cat("\n")
  cat("  Traits:       ", paste(x$traits, collapse = ", "), "\n")
  cat("  Allele tests: ", nrow(x$allele_tests), "\n")
  cat("  Block tests:  ", nrow(x$block_tests), "\n")
  cat("  sig_metric:   ", x$sig_metric,
      " (significant when ", x$sig_metric, " <= ", x$sig_threshold, ")\n", sep = "")
  if (!is.null(x$meff_scope))
    cat("  meff_scope:   ", x$meff_scope,
        " (cut = ", x$meff_percent_cut, ")\n", sep = "")
  cat("  Columns always present:\n")
  cat("    allele_tests: p_wald, p_fdr, p_simplem, p_simplem_sidak, Meff,\n")
  cat("                  alpha_simplem, alpha_simplem_sidak\n")
  cat("    block_tests:  p_omnibus, p_omnibus_fdr, p_omnibus_adj,\n")
  cat("                  p_omnibus_simplem, p_omnibus_simplem_sidak, Meff,\n")
  cat("                  alpha_simplem, alpha_simplem_sidak\n")
  if (nrow(x$allele_tests) > 0)
    cat("  Significant alleles (", x$sig_metric, " <= ", x$sig_threshold, "): ",
        sum(x$allele_tests$significant, na.rm = TRUE), "\n", sep = "")
  if (nrow(x$block_tests) > 0)
    cat("  Significant blocks  (", x$sig_metric, " <= ", x$sig_threshold, "): ",
        sum(x$block_tests$significant_omnibus, na.rm = TRUE), "\n", sep = "")
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
#'   via \code{\link{phase_with_beagle}}) is strongly recommended.
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
#'       \item \code{n_class} (integer) - Number of individuals carrying
#'         this specific diplotype combination.
#'       \item \code{n_total} (integer) - Total individuals with non-missing
#'         data in this block (sum of n_class across all diplotype classes).
#'       \item \code{mean_blue} (numeric) - Mean de-regressed phenotype value
#'         for this diplotype class (on the residual scale after GRM correction).
#'       \item \code{se_mean} (numeric) - Standard error of the mean
#'         (\code{sd / sqrt(n_class)}).
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
          trait=tr, diplotype=dp,
          n=length(vals),               # individuals carrying THIS diplotype (n per class)
          n_total=nrow(dip_ind),         # all individuals with a valid diplotype call
          # at this block AND with phenotype data;
          # n_total >= sum(n_class) because rare
          # diplotype classes (< min_n_diplotype)
          # contribute to n_total but not to any row
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


# ==============================================================================
# compare_block_effects
# ==============================================================================


# ==============================================================================
# .match_blocks_by_position
# Internal: match blocks from two populations by genomic interval overlap.
#
# When LD block boundaries differ between populations (different ancestry,
# different CLQcut, different segmentation), the same causal QTL region will
# have different block_id strings. .match_blocks_by_position() finds the
# best positional match for every Pop1 block in Pop2 using
# Intersection-over-Union (IoU) in base pairs.
#
# Returns a data.frame:
#   block_id_pop1, chr_pop1, start_pop1, end_pop1
#   block_id_pop2, chr_pop2, start_pop2, end_pop2  (NA when no match)
#   overlap_ratio  (IoU; NA when no match)
#   match_type:  "exact"     -- same block_id string
#                "position"  -- matched by positional overlap >= overlap_min
#                "pop1_only" -- no Pop2 block overlaps at >= overlap_min
# ==============================================================================

.match_blocks_by_position <- function(blocks1, blocks2,
                                      overlap_min = 0.50) {
  .nb <- function(b) {
    if (!"CHR" %in% names(b))
      stop("Block table must contain a CHR column.", call. = FALSE)
    if (!"start.bp" %in% names(b)) {
      if ("start_bp" %in% names(b)) {
        b$start.bp <- as.integer(b$start_bp)
        b$end.bp   <- as.integer(b$end_bp)
      } else stop("Block table must contain start.bp/end.bp or start_bp/end_bp.",
                  call. = FALSE)
    }
    if (!"block_id" %in% names(b))
      b$block_id <- paste0("block_", .norm_chr_hap(b$CHR), "_",
                           b$start.bp, "_", b$end.bp)
    b$CHR <- .norm_chr_hap(as.character(b$CHR))
    b[!duplicated(b$block_id),
      c("block_id", "CHR", "start.bp", "end.bp"), drop = FALSE]
  }

  b1 <- .nb(blocks1); b2 <- .nb(blocks2)
  rows <- vector("list", nrow(b1))

  for (i in seq_len(nrow(b1))) {
    id1  <- b1$block_id[i];  chr1 <- b1$CHR[i]
    s1   <- b1$start.bp[i];  e1   <- b1$end.bp[i]

    cands <- b2[b2$CHR == chr1 & b2$end.bp >= s1 & b2$start.bp <= e1,
                , drop = FALSE]

    if (!nrow(cands)) {
      rows[[i]] <- data.frame(
        block_id_pop1 = id1,  chr_pop1 = chr1, start_pop1 = s1, end_pop1 = e1,
        block_id_pop2 = NA_character_, chr_pop2 = NA_character_,
        start_pop2 = NA_integer_, end_pop2 = NA_integer_,
        overlap_ratio = NA_real_, match_type = "pop1_only",
        stringsAsFactors = FALSE)
      next
    }

    inter  <- pmax(0L, pmin(cands$end.bp, e1) - pmax(cands$start.bp, s1))
    union_ <- pmax(1L, pmax(cands$end.bp, e1) - pmin(cands$start.bp, s1))
    iou    <- inter / union_
    best_i <- which.max(iou); best_r <- iou[best_i]

    if (best_r < overlap_min) {
      rows[[i]] <- data.frame(
        block_id_pop1 = id1, chr_pop1 = chr1, start_pop1 = s1, end_pop1 = e1,
        block_id_pop2 = NA_character_, chr_pop2 = NA_character_,
        start_pop2 = NA_integer_, end_pop2 = NA_integer_,
        overlap_ratio = round(best_r, 4), match_type = "pop1_only",
        stringsAsFactors = FALSE)
    } else {
      id2 <- cands$block_id[best_i]
      rows[[i]] <- data.frame(
        block_id_pop1 = id1,  chr_pop1 = chr1,  start_pop1 = s1, end_pop1 = e1,
        block_id_pop2 = id2,  chr_pop2 = cands$CHR[best_i],
        start_pop2 = cands$start.bp[best_i], end_pop2 = cands$end.bp[best_i],
        overlap_ratio = round(best_r, 4),
        match_type = if (id1 == id2) "exact" else "position",
        stringsAsFactors = FALSE)
    }
  }
  do.call(rbind, rows)
}

#' Cross-Population Haplotype Effect Concordance
#'
#' @description
#' Given two \code{allele_tests} data frames produced by
#' \code{\link{test_block_haplotypes}} on two independent populations,
#' computes per-block (and optionally per-trait) statistics that quantify
#' how consistently haplotype allele effects replicate across populations.
#'
#' \strong{Statistics computed per block:}
#' \describe{
#'   \item{Effect correlation (\eqn{r})}{Pearson correlation of per-allele
#'     effect sizes between populations across all shared alleles. Requires
#'     at least 3 shared alleles; \code{NA} otherwise.}
#'   \item{Direction agreement}{Proportion of shared alleles where both
#'     populations assign the same sign to the effect. A value >= 0.75
#'     (i.e. at least 3 out of 4 alleles agree in direction) is considered
#'     strong directional replication.}
#'   \item{Inverse-variance weighted (IVW) meta-analytic effect}{The weighted
#'     mean effect \eqn{\hat\beta_{\mathrm{IVW}} = \sum w_i \beta_i /
#'     \sum w_i} where \eqn{w_i = 1/\mathrm{SE}_i^2}.  Computed separately
#'     for each shared allele and summarised as the mean of per-allele IVW
#'     estimates. The IVW meta-analysis is the same framework used in
#'     two-sample Mendelian randomisation and cross-population GWAS
#'     meta-analysis (see Borenstein et al. 2009).}
#'   \item{Cochran's Q heterogeneity statistic}{Tests whether effect sizes
#'     differ significantly between the two populations:
#'     \eqn{Q = \sum w_i (\beta_i - \hat\beta_{\mathrm{IVW}})^2}.
#'     Under the null of no heterogeneity, \eqn{Q \sim \chi^2_{n-1}}.
#'     Significant Q (large \eqn{Q_p}) indicates that effect sizes differ
#'     between populations - a sign of GxE interaction, LD structure
#'     differences, or population-specific allelic action.}
#'   \item{\eqn{I^2} inconsistency}{
#'     \eqn{I^2 = 100 \times \max(0,\, (Q - df)/Q)}.
#'     Values > 50\% indicate substantial between-population heterogeneity.}
#'   \item{Block boundary concordance}{When \code{blocks_pop1} and
#'     \code{blocks_pop2} are both supplied, the bp overlap ratio of the two
#'     populations' block definitions is computed for each block. A ratio
#'     < 0.8 flags blocks where LD structure likely differs between
#'     populations and effect comparisons should be interpreted cautiously.}
#' }
#'
#' \strong{Addressing the two principal limitations of cross-population validation:}
#'
#' \enumerate{
#'   \item \emph{Population structure confounding}: The Q+K mixed model in
#'     \code{test_block_haplotypes()} already corrects for within-population
#'     structure via the haplotype GRM.  Cross-population validation
#'     inherently controls between-population confounding by design: the
#'     same haplotype allele must associate with the phenotype in a genetically
#'     distinct background, making false-positive carry-over from Pop A's
#'     stratification implausible.  The meta-analytic \eqn{I^2} statistic
#'     additionally flags cases where the effect sizes are so heterogeneous
#'     that a shared biological mechanism is unlikely.
#'
#'   \item \emph{LD block boundary differences}: Block boundaries can differ
#'     between populations due to different historical recombination rates
#'     and ancestral haplotype structure.  The \code{boundary_overlap_ratio}
#'     column quantifies this directly for every block.  Low overlap (< 0.8)
#'     means the two populations carve the same genomic region into blocks of
#'     different sizes; in that case a haplotype allele tested over a 50 kb
#'     window in Pop A is compared to an allele over a 30 kb window in Pop B
#'     and the strings will not match well.  The \code{n_shared_alleles}
#'     column is the empirical consequence: very low shared-allele counts
#'     despite adequate allele frequencies in both populations are a direct
#'     symptom of block boundary mismatch.  The recommended remedy - passing
#'     Pop A's block table as the \code{blocks} argument to
#'     \code{test_block_haplotypes()} in Pop B - eliminates this issue
#'     entirely by forcing both runs to use identical coordinates.
#' }
#'
#' @param assoc_pop1 Output of \code{\link{test_block_haplotypes}} for
#'   population 1 (discovery). Must contain an \code{allele_tests} data frame
#'   with columns \code{block_id}, \code{CHR}, \code{start_bp}, \code{end_bp},
#'   \code{trait}, \code{allele}, \code{effect}, \code{SE}, \code{p_wald}.
#' @param assoc_pop2 Output of \code{\link{test_block_haplotypes}} for
#'   population 2 (validation). Same structure as \code{assoc_pop1}.
#' @param pop1_name Character. Label for population 1 in the output.
#'   Default \code{"pop1"}.
#' @param pop2_name Character. Label for population 2 in the output.
#'   Default \code{"pop2"}.
#' @param traits Character vector or \code{NULL}. Traits to compare. If
#'   \code{NULL} (default), all traits present in both result objects are
#'   included.
#' @param min_shared_alleles Integer. Minimum number of alleles shared
#'   between the two populations for a block to be included in the output.
#'   Blocks with fewer shared alleles are retained but marked with
#'   \code{enough_shared = FALSE}. Default \code{2L}.
#' @param blocks_pop1 Optional data frame of LD blocks from population 1
#'   (output of \code{\link{run_Big_LD_all_chr}}). When both
#'   \code{blocks_pop1} and \code{blocks_pop2} are supplied, a
#'   \code{boundary_overlap_ratio} column is computed for every block.
#'   Required columns: \code{block_id} (or constructible from
#'   \code{CHR}/\code{start.bp}/\code{end.bp}), \code{start.bp},
#'   \code{end.bp}.
#' @param blocks_pop2 Optional data frame of LD blocks from population 2.
#'   Same structure as \code{blocks_pop1}.
#' @param block_match Character. How to match blocks between populations.
#'   \itemize{
#'     \item \code{"id"} (default) - match by exact \code{block_id} string.
#'       Fast and backward-compatible. Use when both populations were analysed
#'       with the same block table (recommended workflow: run
#'       \code{test_block_haplotypes()} on both populations using Pop A's block
#'       boundaries for Pop B as well).
#'     \item \code{"position"} - match by genomic interval overlap
#'       (Intersection-over-Union, IoU, in base pairs). For each Pop1 block,
#'       the best-matching Pop2 block on the same chromosome is found. Blocks
#'       with IoU \eqn{\geq} \code{overlap_min} are matched; those below are
#'       labelled \code{"pop1_only"}. Use when block boundaries genuinely
#'       differ between populations (different ancestral LD structures,
#'       different \code{CLQcut} used, or independent block detection runs).
#'   }
#' @param overlap_min Numeric in (0, 1]. Minimum Intersection-over-Union (IoU)
#'   in base pairs for two blocks to be considered the same region when
#'   \code{block_match = "position"}. Default \code{0.50}. Blocks below this
#'   threshold are assigned \code{match_type = "pop1_only"} and not compared.
#'   Lower values (e.g. \code{0.30}) tolerate more boundary discordance;
#'   higher values (e.g. \code{0.80}) require tighter boundary agreement.
#'   Ignored when \code{block_match = "id"}.
#' @param direction_threshold Numeric in (0.5, 1]. Minimum direction-agreement
#'   proportion to consider a block directionally concordant. Default
#'   \code{0.75}.
#' @param boundary_overlap_warn Numeric in (0, 1). Threshold used to raise a
#'   warning flag in the output. When the automatically computed
#'   \code{boundary_overlap_ratio} (an output column, not a user-set value) is
#'   below this threshold, the corresponding \code{boundary_warning} output
#'   column is set to \code{TRUE}. Default \code{0.80}. Raise this value
#'   (e.g. \code{0.90}) to flag more conservatively, or lower it (e.g.
#'   \code{0.60}) to flag only severely mismatched blocks.
#' @param verbose Logical. Print progress. Default \code{TRUE}.
#'
#' @return A named list of class \code{c("LDxBlocks_effect_concordance", "list")}:
#' \describe{
#'   \item{\code{concordance}}{Data frame with one row per block per trait.
#'     Columns:
#'     \itemize{
#'       \item \code{block_id}, \code{CHR}, \code{start_bp}, \code{end_bp},
#'         \code{trait} - block coordinates and trait name.
#'       \item \code{n_alleles_pop1}, \code{n_alleles_pop2} - alleles tested
#'         in each population (before intersection).
#'       \item \code{n_shared_alleles} - alleles present in both populations.
#'       \item \code{enough_shared} - logical; \code{TRUE} when
#'         \code{n_shared_alleles >= min_shared_alleles}.
#'       \item \code{effect_correlation} - Pearson r of per-allele effects
#'         across populations (NA when n_shared < 3).
#'       \item \code{direction_agreement} - proportion of shared alleles with
#'         concordant effect signs.
#'       \item \code{directionally_concordant} - logical; direction_agreement
#'         >= \code{direction_threshold}.
#'       \item \code{meta_effect} - IVW meta-analytic effect (mean over shared
#'         alleles, each weighted by combined inverse-variance).
#'       \item \code{meta_SE} - SE of the IVW estimate.
#'       \item \code{meta_z} - meta-analytic z-score.
#'       \item \code{meta_p} - two-sided p-value of meta-analytic effect.
#'       \item \code{Q_stat} - Cochran Q heterogeneity statistic.
#'       \item \code{Q_df} - degrees of freedom of Q (n_shared - 1).
#'       \item \code{Q_p} - p-value of Q under chi-squared distribution.
#'       \item \code{I2} - \eqn{I^2} inconsistency (0-100\%).
#'       \item \code{replicated} - logical; \code{TRUE} when
#'         \code{directionally_concordant} AND \code{Q_p > 0.05} AND
#'         \code{enough_shared}.
#'       \item \code{boundary_overlap_ratio} - bp intersection / bp union of
#'         the two populations' block boundaries. \code{NA} when block tables
#'         not supplied.
#'       \item \code{boundary_warning} - logical; \code{TRUE} when
#'         \code{boundary_overlap_ratio < boundary_overlap_warn}. These blocks
#'         should be interpreted cautiously because different LD structures
#'         likely produce non-comparable haplotype strings.
#'       \item \code{match_type} - character; how the block was matched between
#'         populations: \code{"exact"} (same \code{block_id} string),
#'         \code{"position"} (matched by genomic IoU \eqn{\geq}
#'         \code{overlap_min} when \code{block_match = "position"}),
#'         \code{"pop1_only"} (no Pop2 block overlapped at the threshold),
#'         or \code{NA} when no block tables were supplied.
#'     }
#'     Sorted by \code{CHR}, \code{start_bp}, \code{meta_p} (ascending).}
#'   \item{\code{shared_alleles}}{Data frame. One row per shared allele per
#'     block per trait. Contains \code{effect_pop1}, \code{SE_pop1},
#'     \code{p_wald_pop1}, \code{effect_pop2}, \code{SE_pop2},
#'     \code{p_wald_pop2}, \code{direction_agree}, \code{ivw_effect},
#'     \code{ivw_SE}, for detailed per-allele inspection.}
#'   \item{\code{pop1_name}, \code{pop2_name}}{Character. Population labels.}
#'   \item{\code{block_match}}{Character. Matching strategy used.}
#'   \item{\code{overlap_min}}{Numeric. IoU threshold used (relevant when
#'     \code{block_match = "position"}).}
#'   \item{\code{traits}}{Character vector of traits compared.}
#'   \item{\code{direction_threshold}}{Numeric. Threshold used.}
#'   \item{\code{boundary_overlap_warn}}{Numeric. Boundary warning threshold.}
#' }
#'
#' @examples
#' \donttest{
#' data(ldx_geno, ldx_snp_info, ldx_blocks, ldx_blues, package = "LDxBlocks")
#'
#' # Simulate two populations by splitting the example dataset
#' set.seed(1L)
#' n     <- nrow(ldx_geno)
#' idx_1 <- sample(n, round(n * 0.6))
#' idx_2 <- setdiff(seq_len(n), idx_1)
#'
#' haps_1 <- extract_haplotypes(ldx_geno[idx_1, ], ldx_snp_info, ldx_blocks)
#' haps_2 <- extract_haplotypes(ldx_geno[idx_2, ], ldx_snp_info, ldx_blocks)
#'
#' blues_1 <- setNames(ldx_blues$YLD[idx_1], ldx_blues$id[idx_1])
#' blues_2 <- setNames(ldx_blues$YLD[idx_2], ldx_blues$id[idx_2])
#'
#' res_1 <- test_block_haplotypes(haps_1, blues = blues_1,
#'                                blocks = ldx_blocks, verbose = FALSE)
#' res_2 <- test_block_haplotypes(haps_2, blues = blues_2,
#'                                blocks = ldx_blocks, verbose = FALSE)
#'
#' conc <- compare_block_effects(res_1, res_2,
#'                               pop1_name = "pop1", pop2_name = "pop2",
#'                               blocks_pop1 = ldx_blocks,
#'                               blocks_pop2 = ldx_blocks)
#'
#' # Replicated blocks
#' conc$concordance[conc$concordance$replicated, ]
#'
#' # Full per-allele details
#' head(conc$shared_alleles)
#'
#' # Summary
#' print(conc)
#' }
#'
#' @references
#' Borenstein M, Hedges LV, Higgins JPT, Rothstein HR (2009).
#' \emph{Introduction to Meta-Analysis}. Wiley.
#'
#' Higgins JPT, Thompson SG (2002). Quantifying heterogeneity in a
#' meta-analysis. \emph{Statistics in Medicine} \strong{21}(11):1539-1558.
#' \doi{10.1002/sim.1186}
#'
#' @seealso \code{\link{test_block_haplotypes}},
#'   \code{\link{harmonize_haplotypes}},
#'   \code{\link{score_favorable_haplotypes}}
#' @export
compare_block_effects <- function(
    assoc_pop1,
    assoc_pop2,
    pop1_name             = "pop1",
    pop2_name             = "pop2",
    traits                = NULL,
    min_shared_alleles    = 2L,
    blocks_pop1           = NULL,
    blocks_pop2           = NULL,
    block_match           = c("id", "position"),
    overlap_min           = 0.50,
    direction_threshold   = 0.75,
    boundary_overlap_warn = 0.80,
    verbose               = TRUE
) {
  block_match <- match.arg(block_match)
  .log <- function(...) if (verbose) message("[compare_block_effects] ", ...)

  # -- Input validation --------------------------------------------------------
  if (!inherits(assoc_pop1, "LDxBlocks_haplotype_assoc"))
    stop("assoc_pop1 must be of class LDxBlocks_haplotype_assoc (output of test_block_haplotypes()).", call. = FALSE)
  if (!inherits(assoc_pop2, "LDxBlocks_haplotype_assoc"))
    stop("assoc_pop2 must be of class LDxBlocks_haplotype_assoc (output of test_block_haplotypes()).", call. = FALSE)

  at1 <- assoc_pop1$allele_tests
  at2 <- assoc_pop2$allele_tests
  req_cols <- c("block_id","CHR","start_bp","end_bp","trait",
                "allele","effect","SE","p_wald")
  for (nm in c("assoc_pop1","assoc_pop2")) {
    df <- get(paste0("at", if (nm == "assoc_pop1") 1 else 2))
    miss <- setdiff(req_cols, names(df))
    if (length(miss))
      stop(nm, " allele_tests missing columns: ", paste(miss, collapse=", "),
           call. = FALSE)
  }

  # -- Trait selection ---------------------------------------------------------
  if (is.null(traits)) {
    traits <- intersect(unique(at1$trait), unique(at2$trait))
    if (!length(traits))
      stop("No common traits between the two result objects.", call. = FALSE)
    .log("Comparing traits: ", paste(traits, collapse=", "))
  } else {
    miss1 <- setdiff(traits, unique(at1$trait))
    miss2 <- setdiff(traits, unique(at2$trait))
    if (length(miss1)) warning("Traits not in pop1: ", paste(miss1, collapse=", "),
                               call. = FALSE)
    if (length(miss2)) warning("Traits not in pop2: ", paste(miss2, collapse=", "),
                               call. = FALSE)
    traits <- intersect(traits, intersect(unique(at1$trait), unique(at2$trait)))
    if (!length(traits))
      stop("No usable common traits after filtering.", call. = FALSE)
  }

  # -- Pre-compute block boundary overlap if tables supplied -------------------
  # -- Block matching: by ID (fast) or by genomic position (handles different boundaries)
  # block_map: data.frame mapping Pop1 block_id -> Pop2 block_id + overlap_ratio
  # overlap_map: named vector block_id_pop1 -> overlap_ratio (for output column)
  block_map   <- NULL
  overlap_map <- NULL

  if (!is.null(blocks_pop1) && !is.null(blocks_pop2)) {
    if (block_match == "position") {
      block_map <- .match_blocks_by_position(blocks_pop1, blocks_pop2,
                                             overlap_min = overlap_min)
      overlap_map <- stats::setNames(block_map$overlap_ratio,
                                     block_map$block_id_pop1)
      n_pos <- sum(block_map$match_type == "position", na.rm = TRUE)
      n_ex  <- sum(block_map$match_type == "exact",    na.rm = TRUE)
      n_mis <- sum(block_map$match_type == "pop1_only", na.rm = TRUE)
      .log("Position matching: ", n_ex, " exact, ", n_pos,
           " position-matched, ", n_mis, " Pop1-only (IoU < ", overlap_min, ")")
    } else {
      # block_match == "id": compute overlap only for blocks that share the same ID
      block_map <- .match_blocks_by_position(blocks_pop1, blocks_pop2,
                                             overlap_min = 0)
      block_map$match_type <- ifelse(
        block_map$block_id_pop1 == block_map$block_id_pop2 &
          !is.na(block_map$block_id_pop2), "exact", "pop1_only")
      block_map$block_id_pop2[block_map$match_type == "pop1_only"] <- NA_character_
      overlap_map <- stats::setNames(
        ifelse(block_map$match_type == "exact", block_map$overlap_ratio, NA_real_),
        block_map$block_id_pop1)
    }
  }

  # -- Main loop: per-trait, per-block comparison ------------------------------
  conc_rows  <- list()
  shared_rows <- list()

  for (tr in traits) {
    a1_tr <- at1[at1$trait == tr, ]
    a2_tr <- at2[at2$trait == tr, ]

    # When block_match = "position", remap Pop2 block IDs to their Pop1 equivalents
    # so that blocks with different boundaries but overlapping regions are compared.
    if (!is.null(block_map) && block_match == "position") {
      # For each Pop2 allele, look up its block's Pop1 equivalent ID
      remap <- stats::setNames(block_map$block_id_pop1, block_map$block_id_pop2)
      remap <- remap[!is.na(names(remap))]
      a2_tr$block_id_orig <- a2_tr$block_id   # preserve original for output
      a2_tr$block_id <- ifelse(a2_tr$block_id %in% names(remap),
                               remap[a2_tr$block_id],
                               a2_tr$block_id)
    }

    all_blocks <- union(a1_tr$block_id, a2_tr$block_id)

    for (bn in all_blocks) {
      a1_b <- a1_tr[a1_tr$block_id == bn, ]
      a2_b <- a2_tr[a2_tr$block_id == bn, ]

      n1 <- nrow(a1_b); n2 <- nrow(a2_b)
      # Block metadata from whichever population has data
      meta_src <- if (nrow(a1_b) > 0) a1_b else a2_b
      blk_chr   <- meta_src$CHR[1]
      blk_start <- meta_src$start_bp[1]
      blk_end   <- meta_src$end_bp[1]

      # Shared alleles: same allele string in both populations
      shared_alleles <- intersect(a1_b$allele, a2_b$allele)
      n_shared <- length(shared_alleles)
      enough   <- n_shared >= min_shared_alleles

      # Block boundary info
      overlap_ratio    <- if (!is.null(overlap_map) && bn %in% names(overlap_map))
        overlap_map[[bn]] else NA_real_
      boundary_warning <- !is.na(overlap_ratio) && overlap_ratio < boundary_overlap_warn

      if (!n_shared) {
        # No shared alleles - still report the block
        conc_rows[[length(conc_rows) + 1L]] <- data.frame(
          block_id             = bn,
          CHR                  = blk_chr,
          start_bp             = blk_start,
          end_bp               = blk_end,
          trait                = tr,
          n_alleles_pop1       = n1,
          n_alleles_pop2       = n2,
          n_shared_alleles     = 0L,
          enough_shared        = FALSE,
          effect_correlation   = NA_real_,
          direction_agreement  = NA_real_,
          directionally_concordant = FALSE,
          meta_effect          = NA_real_,
          meta_SE              = NA_real_,
          meta_z               = NA_real_,
          meta_p               = NA_real_,
          Q_stat               = NA_real_,
          Q_df                 = NA_integer_,
          Q_p                  = NA_real_,
          I2                   = NA_real_,
          replicated           = FALSE,
          boundary_overlap_ratio = overlap_ratio,
          boundary_warning     = boundary_warning,
          match_type           = if (!is.null(block_map)) {
            mt <- block_map$match_type[block_map$block_id_pop1 == bn]
            if (length(mt)) mt[1] else NA_character_
          } else NA_character_,
          stringsAsFactors     = FALSE
        )
        next
      }

      # Extract matched effects for shared alleles
      e1  <- a1_b$effect[match(shared_alleles, a1_b$allele)]
      se1 <- a1_b$SE[match(shared_alleles, a1_b$allele)]
      p1  <- a1_b$p_wald[match(shared_alleles, a1_b$allele)]
      e2  <- a2_b$effect[match(shared_alleles, a2_b$allele)]
      se2 <- a2_b$SE[match(shared_alleles, a2_b$allele)]
      p2  <- a2_b$p_wald[match(shared_alleles, a2_b$allele)]

      # Guard: drop rows where either SE is zero or NA (degenerate)
      valid <- is.finite(e1) & is.finite(e2) &
        is.finite(se1) & is.finite(se2) &
        se1 > 0 & se2 > 0
      e1 <- e1[valid]; se1 <- se1[valid]; p1 <- p1[valid]
      e2 <- e2[valid]; se2 <- se2[valid]; p2 <- p2[valid]
      shared_alleles_v <- shared_alleles[valid]
      n_v <- length(shared_alleles_v)

      # -- Effect correlation -------------------------------------------------
      eff_cor <- if (n_v >= 3L) {
        tryCatch(stats::cor(e1, e2, use = "complete.obs"), error = function(e) NA_real_)
      } else NA_real_

      # -- Direction agreement -----------------------------------------------
      dir_agree <- if (n_v > 0L) mean(sign(e1) == sign(e2)) else NA_real_
      dir_concordant <- !is.na(dir_agree) && dir_agree >= direction_threshold

      # -- IVW meta-analysis per allele, then average -----------------------
      # For each shared allele: w_i = 1/(se1_i^2 + se2_i^2)
      # IVW combined effect for allele i: b_ivw_i = (e1*w1 + e2*w2) / (w1 + w2)
      # where w1 = 1/se1^2, w2 = 1/se2^2
      w1  <- 1 / se1^2
      w2  <- 1 / se2^2
      b_ivw  <- (e1 * w1 + e2 * w2) / (w1 + w2)
      se_ivw <- sqrt(1 / (w1 + w2))

      # Combined weight across both populations for block-level summary
      w_comb    <- w1 + w2
      meta_eff  <- if (n_v > 0L) sum(b_ivw * w_comb) / sum(w_comb) else NA_real_
      meta_SE_v <- if (n_v > 0L) sqrt(1 / sum(w_comb)) else NA_real_
      meta_z_v  <- if (!is.na(meta_eff) && !is.na(meta_SE_v) && meta_SE_v > 0)
        meta_eff / meta_SE_v else NA_real_
      meta_p_v  <- if (!is.na(meta_z_v))
        2 * stats::pnorm(-abs(meta_z_v)) else NA_real_

      # -- Cochran Q heterogeneity -------------------------------------------
      # Q = sum over alleles of [ sum_pop w_{ij} (e_{ij} - b_ivw_i)^2 ]
      # where i = allele, j = population
      Q_v  <- sum(w1 * (e1 - b_ivw)^2 + w2 * (e2 - b_ivw)^2)
      Q_df_v <- as.integer(max(0L, n_v - 1L))  # df = n_alleles - 1
      Q_p_v  <- if (Q_df_v > 0L) stats::pchisq(Q_v, df = Q_df_v,
                                               lower.tail = FALSE) else NA_real_
      I2_v   <- if (!is.na(Q_v) && !is.na(Q_df_v) && Q_v > 0)
        100 * max(0, (Q_v - Q_df_v) / Q_v) else 0

      replicated <- enough & dir_concordant & !is.na(Q_p_v) & Q_p_v > 0.05

      conc_rows[[length(conc_rows) + 1L]] <- data.frame(
        block_id             = bn,
        CHR                  = blk_chr,
        start_bp             = blk_start,
        end_bp               = blk_end,
        trait                = tr,
        n_alleles_pop1       = n1,
        n_alleles_pop2       = n2,
        n_shared_alleles     = n_v,
        enough_shared        = enough,
        effect_correlation   = round(eff_cor,    4),
        direction_agreement  = round(dir_agree,  4),
        directionally_concordant = dir_concordant,
        meta_effect          = round(meta_eff,   6),
        meta_SE              = round(meta_SE_v,  6),
        meta_z               = round(meta_z_v,   4),
        meta_p               = meta_p_v,
        Q_stat               = round(Q_v,        4),
        Q_df                 = Q_df_v,
        Q_p                  = Q_p_v,
        I2                   = round(I2_v,       2),
        replicated           = replicated,
        boundary_overlap_ratio = round(overlap_ratio, 4),
        boundary_warning     = boundary_warning,
        match_type           = if (!is.null(block_map)) {
          mt <- block_map$match_type[block_map$block_id_pop1 == bn]
          if (length(mt)) mt[1] else NA_character_
        } else NA_character_,
        stringsAsFactors     = FALSE
      )

      # -- Per-allele detail rows -------------------------------------------
      if (n_v > 0L) {
        for (ai in seq_len(n_v)) {
          shared_rows[[length(shared_rows) + 1L]] <- data.frame(
            block_id       = bn,
            CHR            = blk_chr,
            start_bp       = blk_start,
            end_bp         = blk_end,
            trait          = tr,
            allele         = shared_alleles_v[ai],
            effect_pop1    = round(e1[ai],  6),
            SE_pop1        = round(se1[ai], 6),
            p_wald_pop1    = p1[ai],
            effect_pop2    = round(e2[ai],  6),
            SE_pop2        = round(se2[ai], 6),
            p_wald_pop2    = p2[ai],
            direction_agree = sign(e1[ai]) == sign(e2[ai]),
            ivw_effect     = round(b_ivw[ai],  6),
            ivw_SE         = round(se_ivw[ai], 6),
            stringsAsFactors = FALSE
          )
        }
      }
    }  # end block loop
  }  # end trait loop

  # -- Assemble and sort -------------------------------------------------------
  conc_df   <- if (length(conc_rows))
    do.call(rbind, conc_rows) else data.frame()
  shared_df <- if (length(shared_rows))
    do.call(rbind, shared_rows) else data.frame()

  if (nrow(conc_df) > 0) {
    conc_df <- conc_df[order(conc_df$CHR, conc_df$start_bp,
                             ifelse(is.na(conc_df$meta_p), 1, conc_df$meta_p)), ]
    rownames(conc_df) <- NULL
  }
  if (nrow(shared_df) > 0) rownames(shared_df) <- NULL

  n_rep <- sum(conc_df$replicated, na.rm = TRUE)
  .log("Done. Blocks compared: ", nrow(conc_df),
       " | Replicated (concordant, Q_p > 0.05): ", n_rep)

  structure(
    list(
      concordance           = conc_df,
      shared_alleles        = shared_df,
      pop1_name             = pop1_name,
      pop2_name             = pop2_name,
      traits                = traits,
      direction_threshold   = direction_threshold,
      boundary_overlap_warn = boundary_overlap_warn,
      block_match           = block_match,
      overlap_min           = overlap_min
    ),
    class = c("LDxBlocks_effect_concordance", "list")
  )
}


#' @export
print.LDxBlocks_effect_concordance <- function(x, ...) {
  cat("LDxBlocks Cross-Population Effect Concordance\n")
  cat("  Populations: ", x$pop1_name, " vs ", x$pop2_name, "\n", sep = "")
  cat("  Traits:      ", paste(x$traits, collapse = ", "), "\n")
  if (nrow(x$concordance) > 0) {
    cat("  Blocks compared:          ", nrow(x$concordance), "\n")
    cat("  With enough shared alleles:", sum(x$concordance$enough_shared, na.rm = TRUE), "\n")
    cat("  Directionally concordant:  ",
        sum(x$concordance$directionally_concordant, na.rm = TRUE), "\n")
    cat("  Replicated (dir + Q_p>0.05):",
        sum(x$concordance$replicated, na.rm = TRUE), "\n")
    cat("  Boundary warnings:         ",
        sum(x$concordance$boundary_warning, na.rm = TRUE),
        "(overlap ratio <", x$boundary_overlap_warn, ")\n")
    if (any(!is.na(x$concordance$I2)))
      cat("  Median I2 (heterogeneity): ",
          round(median(x$concordance$I2, na.rm = TRUE), 1), "%\n")
  }
  if (nrow(x$shared_alleles) > 0)
    cat("  Shared allele comparisons: ", nrow(x$shared_alleles), "\n")
  invisible(x)
}


# ==============================================================================
# compare_gwas_effects
# Cross-population GWAS effect concordance from external association results
# ==============================================================================

#' Cross-Population GWAS Effect Concordance from External Results
#'
#' @description
#' Compares haplotype block-level association signals between two populations
#' when GWAS was run externally (e.g. in GAPIT, TASSEL, FarmCPU, PLINK, or
#' any other tool) rather than through \code{\link{test_block_haplotypes}}.
#'
#' The function accepts external GWAS results in two ways:
#' \enumerate{
#'   \item \strong{Pre-mapped} (recommended): supply the output of
#'     \code{\link{define_qtl_regions}} for each population. The QTL table
#'     already maps each marker to its LD block, providing the lead SNP,
#'     its effect size (BETA), and its p-value per block. This is the most
#'     reliable input because the block assignment is explicit.
#'   \item \strong{Raw GWAS + blocks}: supply raw GWAS data frames plus block
#'     tables; the function calls \code{\link{define_qtl_regions}} internally
#'     to map markers to blocks before comparing.
#' }
#'
#' \strong{What is compared and how:}
#'
#' External GWAS produces marker-level effects, not haplotype-allele-level
#' effects. At the block level each population contributes one observation:
#' the lead SNP (smallest p-value in the block). The comparison is therefore
#' analogous to two-sample Mendelian randomisation - one effect estimate
#' per exposure (block) per population, combined via IVW meta-analysis.
#'
#' Because there is only one "allele" per block (the lead SNP), several
#' statistics differ from \code{\link{compare_block_effects}}:
#' \itemize{
#'   \item \code{n_shared_alleles} is always 1 (the lead SNP concept).
#'   \item \code{effect_correlation} is always \code{NA} (needs >= 3 alleles).
#'   \item \code{direction_agreement} is 0 or 1 (sign of lead SNP effect agrees
#'     or not).
#'   \item \code{Q_df} is always 0 (one observation per population; Cochran Q
#'     is undefined with two populations and one allele per block).
#'   \item \code{meta_p} tests whether the IVW combined effect differs from
#'     zero - the primary replication test.
#' }
#'
#' \strong{Standard error derivation:}
#'
#' Many GWAS tools report BETA and P but not SE. When SE is absent,
#' \code{compare_gwas_effects()} derives it from the z-score:
#' \deqn{SE = \frac{|\beta|}{|z|},\quad z = \Phi^{-1}(P/2)}
#' where \eqn{\Phi^{-1}} is the standard normal quantile function. This
#' is exact for large-sample Wald tests (the dominant GWAS framework) and
#' provides a valid SE for IVW weighting.
#'
#' \strong{Interpreting pleiotropic blocks:}
#'
#' When \code{trait_col} is supplied and GWAS results contain multiple traits,
#' the function compares trait by trait. Blocks significant for multiple traits
#' in both populations (\code{pleiotropic = TRUE} in both QTL tables) are
#' flagged with \code{both_pleiotropic = TRUE} in the output, identifying
#' the most robust cross-population replication targets.
#'
#' @param qtl_pop1 Data frame. Output of \code{\link{define_qtl_regions}} for
#'   population 1 (discovery). Required columns: \code{block_id}, \code{CHR},
#'   \code{start_bp}, \code{end_bp}, \code{lead_snp}, \code{lead_p},
#'   \code{lead_beta}. Optional: \code{traits}, \code{pleiotropic}.
#'   Ignored when \code{gwas_pop1} is supplied instead.
#' @param qtl_pop2 Data frame. Output of \code{\link{define_qtl_regions}} for
#'   population 2 (validation). Same structure as \code{qtl_pop1}.
#'   Ignored when \code{gwas_pop2} is supplied instead.
#' @param gwas_pop1 Data frame. Raw GWAS results for population 1. Required
#'   when \code{qtl_pop1} is \code{NULL}. Must contain \code{SNP} (or
#'   \code{Marker}), \code{CHR}, \code{POS}, and the columns named by
#'   \code{beta_col}, \code{se_col} (optional), \code{p_col}.
#' @param gwas_pop2 Data frame. Raw GWAS results for population 2. Same
#'   structure as \code{gwas_pop1}.
#' @param blocks_pop1 Data frame. LD block table for population 1 (output of
#'   \code{\link{run_Big_LD_all_chr}}). Required when \code{gwas_pop1} is
#'   supplied. Also used for \code{boundary_overlap_ratio} computation when
#'   \code{qtl_pop1} is the input.
#' @param blocks_pop2 Data frame. LD block table for population 2. Same
#'   structure as \code{blocks_pop1}.
#' @param snp_info_pop1 Data frame. SNP metadata for population 1 (columns:
#'   \code{SNP}, \code{CHR}, \code{POS}). Required when \code{gwas_pop1} is
#'   supplied.
#' @param snp_info_pop2 Data frame. SNP metadata for population 2. Same
#'   structure as \code{snp_info_pop1}. If \code{NULL} (default), \code{snp_info_pop1}
#'   is reused (appropriate when both populations share the same marker panel).
#' @param pop1_name Character. Label for population 1. Default \code{"pop1"}.
#' @param pop2_name Character. Label for population 2. Default \code{"pop2"}.
#' @param beta_col Character. Name of the effect-size column in raw GWAS data
#'   frames. Default \code{"BETA"}.
#' @param se_col Character or \code{NULL}. Name of the standard-error column.
#'   When \code{NULL} or absent, SE is derived from \code{beta_col} and
#'   \code{p_col} via the z-score formula. Default \code{"SE"}.
#' @param p_col Character. Name of the p-value column. Default \code{"P"}.
#' @param trait_col Character. Name of the trait column when GWAS results
#'   contain multiple traits. Default \code{"trait"}.
#' @param p_threshold Numeric. Significance threshold applied when mapping raw
#'   GWAS results to blocks via \code{\link{define_qtl_regions}}. Ignored when
#'   \code{qtl_pop1} / \code{qtl_pop2} are supplied directly. Default
#'   \code{5e-8}.
#' @param min_snps Integer. Minimum number of SNPs in a block for it to be
#'   included in the QTL mapping step. Default \code{3L}.
#' @param block_match Character. How to match blocks between populations.
#'   \code{"id"} (default) matches by exact \code{block_id} string - fast and
#'   backward-compatible, appropriate when both populations were mapped against
#'   the same block table. \code{"position"} matches by genomic interval
#'   overlap (Intersection-over-Union \eqn{\geq} \code{overlap_min}) -
#'   recommended when block boundaries differ between populations due to
#'   different ancestral LD structures or independent block-detection runs.
#' @param overlap_min Numeric in (0, 1]. Minimum Intersection-over-Union (IoU)
#'   in base pairs for two blocks to be considered the same region when
#'   \code{block_match = "position"}. Default \code{0.50}. Ignored when
#'   \code{block_match = "id"}.
#' @param direction_threshold Numeric in (0.5, 1]. Minimum direction agreement
#'   to consider a block directionally concordant. Default \code{0.75}.
#'   Note: with one lead SNP per block this is effectively a 0/1 flag at
#'   this threshold, but the parameter is kept for consistency with
#'   \code{\link{compare_block_effects}}.
#' @param boundary_overlap_warn Numeric in (0, 1). \code{boundary_overlap_ratio}
#'   below this value triggers \code{boundary_warning = TRUE}. The
#'   \code{boundary_overlap_ratio} is an **output column** automatically
#'   computed from the block tables - it cannot be set by the user.
#'   Default \code{0.80}.
#' @param verbose Logical. Print progress. Default \code{TRUE}.
#'
#' @return A named list of class \code{c("LDxBlocks_effect_concordance", "list")}
#'   with the same structure as \code{\link{compare_block_effects}}. Additional
#'   columns in \code{$concordance} specific to GWAS input:
#'   \itemize{
#'     \item \code{lead_snp_pop1}, \code{lead_snp_pop2} - lead SNP ID from
#'       each population (same SNP = same tag; different SNP = different LD
#'       tag, same region).
#'     \item \code{lead_p_pop1}, \code{lead_p_pop2} - lead SNP p-values.
#'     \item \code{se_derived_pop1}, \code{se_derived_pop2} - logical; \code{TRUE}
#'       when SE was derived from BETA and P rather than read directly.
#'     \item \code{both_pleiotropic} - logical; \code{TRUE} when the block is
#'       pleiotropic in both populations (requires \code{pleiotropic} column
#'       in QTL tables).
#'   }
#'   The \code{$shared_alleles} data frame contains one row per block per trait
#'   with \code{lead_snp} instead of \code{allele}, and \code{effect_pop1},
#'   \code{SE_pop1}, \code{effect_pop2}, \code{SE_pop2}, \code{direction_agree},
#'   \code{ivw_effect}, \code{ivw_SE}.
#'
#' @examples
#' \donttest{
#' data(ldx_geno, ldx_snp_info, ldx_blocks, ldx_gwas, package = "LDxBlocks")
#'
#' # Simulate two populations
#' set.seed(1L)
#' n    <- nrow(ldx_geno)
#' idx1 <- sample(n, 70L)
#' idx2 <- setdiff(seq_len(n), idx1)
#'
#' # Suppose we have external GWAS results for both populations
#' # (here we use ldx_gwas as a stand-in)
#' gwas_A <- ldx_gwas; gwas_A$BETA <- rnorm(nrow(gwas_A), 0.5, 0.1)
#' gwas_B <- ldx_gwas; gwas_B$BETA <- rnorm(nrow(gwas_B), 0.4, 0.15)
#'
#' # Path 1: raw GWAS + blocks (convenience)
#' conc <- compare_gwas_effects(
#'   gwas_pop1     = gwas_A,
#'   gwas_pop2     = gwas_B,
#'   blocks_pop1   = ldx_blocks,
#'   blocks_pop2   = ldx_blocks,
#'   snp_info_pop1 = ldx_snp_info,
#'   pop1_name     = "PopA",
#'   pop2_name     = "PopB",
#'   p_threshold   = NULL   # keep all markers for this demo
#' )
#' print(conc)
#'
#' # Path 2: pre-mapped QTL tables (recommended)
#' qtl_A <- define_qtl_regions(gwas_A, ldx_blocks, ldx_snp_info,
#'                              p_threshold = NULL)
#' qtl_B <- define_qtl_regions(gwas_B, ldx_blocks, ldx_snp_info,
#'                              p_threshold = NULL)
#' conc2 <- compare_gwas_effects(
#'   qtl_pop1      = qtl_A,
#'   qtl_pop2      = qtl_B,
#'   blocks_pop1   = ldx_blocks,
#'   blocks_pop2   = ldx_blocks,
#'   pop1_name     = "PopA",
#'   pop2_name     = "PopB"
#' )
#' conc2$concordance[conc2$concordance$replicated, ]
#' }
#'
#' @references
#' Borenstein M, Hedges LV, Higgins JPT, Rothstein HR (2009).
#' \emph{Introduction to Meta-Analysis}. Wiley.
#'
#' Higgins JPT, Thompson SG (2002). Quantifying heterogeneity in a
#' meta-analysis. \emph{Statistics in Medicine} \strong{21}(11):1539-1558.
#' \doi{10.1002/sim.1186}
#'
#' @seealso \code{\link{define_qtl_regions}},
#'   \code{\link{compare_block_effects}},
#'   \code{\link{test_block_haplotypes}}
#' @export
compare_gwas_effects <- function(
    qtl_pop1              = NULL,
    qtl_pop2              = NULL,
    gwas_pop1             = NULL,
    gwas_pop2             = NULL,
    blocks_pop1           = NULL,
    blocks_pop2           = NULL,
    snp_info_pop1         = NULL,
    snp_info_pop2         = NULL,
    pop1_name             = "pop1",
    pop2_name             = "pop2",
    beta_col              = "BETA",
    se_col                = "SE",
    p_col                 = "P",
    trait_col             = "trait",
    p_threshold           = 5e-8,
    min_snps              = 3L,
    block_match           = c("id", "position"),
    overlap_min           = 0.50,
    direction_threshold   = 0.75,
    boundary_overlap_warn = 0.80,
    verbose               = TRUE
) {
  block_match <- match.arg(block_match)
  .log <- function(...) if (verbose) message("[compare_gwas_effects] ", ...)

  # -- Input validation and routing --------------------------------------------
  using_raw <- !is.null(gwas_pop1) || !is.null(gwas_pop2)
  using_qtl <- !is.null(qtl_pop1) || !is.null(qtl_pop2)

  if (!using_raw && !using_qtl)
    stop("Supply either (qtl_pop1 + qtl_pop2) or (gwas_pop1 + gwas_pop2).",
         call. = FALSE)
  if (using_raw && using_qtl)
    stop("Supply either qtl_pop1/qtl_pop2 OR gwas_pop1/gwas_pop2, not both.",
         call. = FALSE)

  # -- Internal: derive SE from BETA and P via z-score ---------------------
  .derive_se <- function(beta, p) {
    # SE = |beta| / |z|  where z = qnorm(p/2) (two-sided)
    # Guard: p = 0 -> z = Inf -> SE = 0 (degenerate); use p = .Machine$double.xmin
    p   <- pmax(p, .Machine$double.xmin)
    z   <- abs(stats::qnorm(p / 2))
    se  <- abs(beta) / z
    se[!is.finite(se) | se <= 0] <- NA_real_
    se
  }

  # -- Internal: normalise a raw GWAS data frame ---------------------------
  .normalise_gwas <- function(gwas, pop_label) {
    if (!is.data.frame(gwas))
      stop(pop_label, ": gwas must be a data frame.", call. = FALSE)
    # SNP ID column
    if (!"SNP" %in% names(gwas) && "Marker" %in% names(gwas))
      gwas$SNP <- gwas$Marker
    req <- c("SNP", "CHR", "POS")
    miss <- setdiff(req, names(gwas))
    if (length(miss))
      stop(pop_label, ": missing columns: ", paste(miss, collapse = ", "),
           call. = FALSE)
    # Try specified column names; fall back to standard "BETA"/"SE"/"P" if not found
    # This allows pop1 and pop2 to have different column naming conventions
    eff_col_use <- if (beta_col %in% names(gwas)) beta_col else
      if ("BETA" %in% names(gwas))   "BETA"   else
        stop(pop_label, ": beta_col '", beta_col, "' not found and 'BETA' also absent.", call. = FALSE)
    p_col_use   <- if (p_col %in% names(gwas))    p_col    else
      if ("P" %in% names(gwas))       "P"      else
        stop(pop_label, ": p_col '", p_col, "' not found and 'P' also absent.", call. = FALSE)
    # Standardise column names
    gwas$BETA <- as.numeric(gwas[[eff_col_use]])
    gwas$P    <- as.numeric(gwas[[p_col_use]])
    # SE: use supplied column or derive from z-score
    se_col_use  <- if (!is.null(se_col) && se_col %in% names(gwas)) se_col else
      if ("SE" %in% names(gwas)) "SE" else NULL
    se_present <- !is.null(se_col_use) && !all(is.na(gwas[[se_col_use]]))
    if (se_present) {
      gwas$SE          <- as.numeric(gwas[[se_col_use]])
      gwas$.se_derived <- FALSE
    } else {
      .log("  ", pop_label, ": SE not found - deriving from BETA and ", p_col,
           " via z-score")
      gwas$SE          <- .derive_se(gwas$BETA, gwas$P)
      gwas$.se_derived <- TRUE
    }
    gwas
  }

  # -- Internal: normalise a QTL table (output of define_qtl_regions) ------
  .normalise_qtl <- function(qtl, pop_label) {
    req <- c("block_id", "CHR", "start_bp", "end_bp", "lead_snp",
             "lead_p", "lead_beta")
    miss <- setdiff(req, names(qtl))
    if (length(miss))
      stop(pop_label, " qtl table missing columns: ", paste(miss, collapse = ", "),
           "\n  Run define_qtl_regions() with a GWAS table that includes BETA.",
           call. = FALSE)
    # Derive SE if not present (lead_beta + lead_p)
    if (!"lead_se" %in% names(qtl) || all(is.na(qtl$lead_se))) {
      .log("  ", pop_label, ": lead_se absent - deriving from lead_beta + lead_p")
      qtl$lead_se       <- .derive_se(qtl$lead_beta, qtl$lead_p)
      qtl$.se_derived   <- TRUE
    } else {
      qtl$.se_derived <- FALSE
    }
    if (!trait_col %in% names(qtl)) qtl[[trait_col]] <- "trait"
    if (!"pleiotropic" %in% names(qtl)) qtl$pleiotropic <- FALSE
    qtl
  }

  # -- Build QTL tables if raw GWAS supplied -------------------------------
  if (using_raw) {
    if (is.null(gwas_pop1) || is.null(gwas_pop2))
      stop("Both gwas_pop1 and gwas_pop2 must be supplied when using raw GWAS.",
           call. = FALSE)
    if (is.null(blocks_pop1))
      stop("blocks_pop1 is required when gwas_pop1 is supplied.", call. = FALSE)
    if (is.null(snp_info_pop1))
      stop("snp_info_pop1 is required when gwas_pop1 is supplied.", call. = FALSE)

    # Use snp_info_pop1 for pop2 if not separately supplied (shared marker panel)
    if (is.null(snp_info_pop2)) snp_info_pop2 <- snp_info_pop1
    if (is.null(blocks_pop2))   blocks_pop2   <- blocks_pop1

    gwas1 <- .normalise_gwas(gwas_pop1, pop1_name)
    gwas2 <- .normalise_gwas(gwas_pop2, pop2_name)

    # Capture SE-derivation flag as plain scalars immediately while gwas1/gwas2
    # are in scope. SE derivation is a population-level property (either ALL
    # markers in a population had SE or none did), so a scalar is correct.
    # Storing as scalar avoids all per-row column extraction complexity.
    se_derived_flag_pop1 <- isTRUE(gwas1$.se_derived[1L])
    se_derived_flag_pop2 <- isTRUE(gwas2$.se_derived[1L])

    .log("Mapping ", pop1_name, " GWAS results to blocks ...")
    qtl1 <- define_qtl_regions(
      gwas_results = gwas1,
      blocks       = blocks_pop1,
      snp_info     = snp_info_pop1,
      p_threshold  = p_threshold,
      trait_col    = if (trait_col %in% names(gwas1)) trait_col else "trait",
      min_snps     = min_snps,
      verbose      = FALSE
    )
    if (nrow(qtl1) == 0L)
      stop("No blocks with significant markers in ", pop1_name, ". ",
           "Try a less stringent p_threshold.", call. = FALSE)
    # Carry over SE values (needed for IVW weighting)
    qtl1$lead_se <- vapply(seq_len(nrow(qtl1)), function(i) {
      idx <- which(gwas1$SNP == qtl1$lead_snp[i])
      if (length(idx)) gwas1$SE[idx[1]] else NA_real_
    }, numeric(1L))

    .log("Mapping ", pop2_name, " GWAS results to blocks ...")
    qtl2 <- define_qtl_regions(
      gwas_results = gwas2,
      blocks       = blocks_pop2,
      snp_info     = snp_info_pop2,
      p_threshold  = p_threshold,
      trait_col    = if (trait_col %in% names(gwas2)) trait_col else "trait",
      min_snps     = min_snps,
      verbose      = FALSE
    )
    if (nrow(qtl2) == 0L)
      stop("No blocks with significant markers in ", pop2_name, ". ",
           "Try a less stringent p_threshold.", call. = FALSE)
    qtl2$lead_se <- vapply(seq_len(nrow(qtl2)), function(i) {
      idx <- which(gwas2$SNP == qtl2$lead_snp[i])
      if (length(idx)) gwas2$SE[idx[1]] else NA_real_
    }, numeric(1L))

    qtl_pop1 <- qtl1; qtl_pop2 <- qtl2
  } else {
    # Pre-mapped QTL tables supplied directly
    qtl_pop1 <- .normalise_qtl(qtl_pop1, pop1_name)
    qtl_pop2 <- .normalise_qtl(qtl_pop2, pop2_name)
    # For pre-mapped path, se_derived flag comes from .normalise_qtl
    # which sets $.se_derived when it derives SE from lead_beta + lead_p
    se_derived_flag_pop1 <- if (".se_derived" %in% names(qtl_pop1) && nrow(qtl_pop1))
      isTRUE(as.logical(qtl_pop1$.se_derived)[1L]) else FALSE
    se_derived_flag_pop2 <- if (".se_derived" %in% names(qtl_pop2) && nrow(qtl_pop2))
      isTRUE(as.logical(qtl_pop2$.se_derived)[1L]) else FALSE
  }

  .log("Comparing effects for blocks with hits in both populations ...")

  # -- Block matching and boundary overlap ------------------------------------
  # block_map_gwas maps QTL block IDs from Pop2 -> Pop1 when block_match="position"
  block_map_gwas <- NULL
  overlap_map    <- NULL

  if (!is.null(blocks_pop1) && !is.null(blocks_pop2)) {
    if (block_match == "position") {
      block_map_gwas <- .match_blocks_by_position(blocks_pop1, blocks_pop2,
                                                  overlap_min = overlap_min)
      overlap_map <- stats::setNames(block_map_gwas$overlap_ratio,
                                     block_map_gwas$block_id_pop1)
      n_pos <- sum(block_map_gwas$match_type == "position", na.rm = TRUE)
      n_ex  <- sum(block_map_gwas$match_type == "exact",    na.rm = TRUE)
      .log("Position matching: ", n_ex, " exact, ", n_pos, " position-matched")
    } else {
      block_map_gwas <- .match_blocks_by_position(blocks_pop1, blocks_pop2,
                                                  overlap_min = 0)
      block_map_gwas$match_type <- ifelse(
        block_map_gwas$block_id_pop1 == block_map_gwas$block_id_pop2 &
          !is.na(block_map_gwas$block_id_pop2), "exact", "pop1_only")
      block_map_gwas$block_id_pop2[block_map_gwas$match_type == "pop1_only"] <- NA_character_
      overlap_map <- stats::setNames(
        ifelse(block_map_gwas$match_type == "exact", block_map_gwas$overlap_ratio, NA_real_),
        block_map_gwas$block_id_pop1)
    }
  }

  # -- Per-trait, per-block comparison -------------------------------------
  traits <- intersect(unique(qtl_pop1[[trait_col]]), unique(qtl_pop2[[trait_col]]))
  if (!length(traits)) {
    # Fall back: treat as single-trait
    qtl_pop1[[trait_col]] <- "trait"
    qtl_pop2[[trait_col]] <- "trait"
    traits <- "trait"
  }
  .log("Traits: ", paste(traits, collapse = ", "))

  conc_rows  <- list()
  shared_rows <- list()

  for (tr in traits) {
    q1_tr <- qtl_pop1[qtl_pop1[[trait_col]] == tr, , drop = FALSE]
    q2_tr <- qtl_pop2[qtl_pop2[[trait_col]] == tr, , drop = FALSE]
    all_blocks <- union(q1_tr$block_id, q2_tr$block_id)
    # Note: q2_tr$block_id may have been remapped above; recompute union

    # When block_match = "position", remap Pop2 QTL block IDs to Pop1 equivalents
    if (!is.null(block_map_gwas) && block_match == "position") {
      remap_gwas <- stats::setNames(block_map_gwas$block_id_pop1,
                                    block_map_gwas$block_id_pop2)
      remap_gwas <- remap_gwas[!is.na(names(remap_gwas))]
      q2_tr$block_id_orig <- q2_tr$block_id
      q2_tr$block_id <- ifelse(q2_tr$block_id %in% names(remap_gwas),
                               remap_gwas[q2_tr$block_id],
                               q2_tr$block_id)
    }

    for (bn in all_blocks) {
      r1 <- q1_tr[q1_tr$block_id == bn, , drop = FALSE]
      r2 <- q2_tr[q2_tr$block_id == bn, , drop = FALSE]

      # Block metadata
      meta_src  <- if (nrow(r1) > 0) r1 else r2
      blk_chr   <- meta_src$CHR[1]
      blk_start <- meta_src$start_bp[1]
      blk_end   <- meta_src$end_bp[1]

      overlap_ratio    <- if (!is.null(overlap_map) && bn %in% names(overlap_map))
        overlap_map[[bn]] else NA_real_
      boundary_warning <- !is.na(overlap_ratio) && overlap_ratio < boundary_overlap_warn

      # Base row with NA results for blocks present in only one population
      base_row <- data.frame(
        block_id             = bn,
        CHR                  = blk_chr,
        start_bp             = blk_start,
        end_bp               = blk_end,
        trait                = tr,
        lead_snp_pop1        = if (nrow(r1)) r1$lead_snp[1] else NA_character_,
        lead_snp_pop2        = if (nrow(r2)) r2$lead_snp[1] else NA_character_,
        lead_p_pop1          = if (nrow(r1)) r1$lead_p[1] else NA_real_,
        lead_p_pop2          = if (nrow(r2)) r2$lead_p[1] else NA_real_,
        # se_derived_pop1/pop2 is a POPULATION-LEVEL flag (did pop1/pop2 lack SE?)
        # Use the scalar captured at the start of the raw/pre-mapped path
        # to avoid platform-specific data.frame column extraction issues.
        se_derived_pop1      = se_derived_flag_pop1,
        se_derived_pop2      = se_derived_flag_pop2,
        n_alleles_pop1       = as.integer(nrow(r1) > 0),
        n_alleles_pop2       = as.integer(nrow(r2) > 0),
        n_shared_alleles     = 0L,
        enough_shared        = FALSE,
        effect_correlation   = NA_real_,
        direction_agreement  = NA_real_,
        directionally_concordant = FALSE,
        meta_effect          = NA_real_,
        meta_SE              = NA_real_,
        meta_z               = NA_real_,
        meta_p               = NA_real_,
        Q_stat               = NA_real_,
        Q_df                 = NA_integer_,
        Q_p                  = NA_real_,
        I2                   = NA_real_,
        replicated           = FALSE,
        both_pleiotropic     = if ("pleiotropic" %in% names(r1) &&
                                   "pleiotropic" %in% names(r2) &&
                                   nrow(r1) > 0 && nrow(r2) > 0)
          isTRUE(r1$pleiotropic[1]) && isTRUE(r2$pleiotropic[1]) else NA,
        boundary_overlap_ratio = round(overlap_ratio, 4),
        boundary_warning       = boundary_warning,
        match_type             = if (!is.null(block_map_gwas)) {
          mt <- block_map_gwas$match_type[block_map_gwas$block_id_pop1 == bn]
          if (length(mt)) mt[1] else NA_character_
        } else NA_character_,
        stringsAsFactors = FALSE
      )

      if (!nrow(r1) || !nrow(r2)) {
        conc_rows[[length(conc_rows) + 1L]] <- base_row
        next
      }

      # Both populations have this block
      e1  <- r1$lead_beta[1];  e2  <- r2$lead_beta[1]
      se1 <- r1$lead_se[1];   se2 <- r2$lead_se[1]

      if (!is.finite(e1) || !is.finite(e2) ||
          !is.finite(se1) || !is.finite(se2) || se1 <= 0 || se2 <= 0) {
        conc_rows[[length(conc_rows) + 1L]] <- base_row
        next
      }

      # Direction agreement (0 or 1 with one observation per population)
      dir_agree     <- as.numeric(sign(e1) == sign(e2))
      dir_concordant <- dir_agree >= direction_threshold

      # IVW meta-analysis (two populations, one lead SNP each)
      w1 <- 1 / se1^2;  w2 <- 1 / se2^2
      b_ivw  <- (e1 * w1 + e2 * w2) / (w1 + w2)
      se_ivw <- sqrt(1 / (w1 + w2))
      meta_z_v <- b_ivw / se_ivw
      meta_p_v <- 2 * stats::pnorm(-abs(meta_z_v))

      # Cochran Q: df = n_alleles - 1 = 0 with one lead SNP
      # Q is undefined; report NA (consistent documentation)
      Q_v    <- NA_real_
      Q_df_v <- NA_integer_
      Q_p_v  <- NA_real_
      I2_v   <- NA_real_

      # Replicated: directionally concordant AND meta_p significant
      # (Q not available for single-allele comparison; use meta_p instead)
      replicated <- dir_concordant && !is.na(meta_p_v) && meta_p_v <= 0.05

      base_row$n_shared_alleles     <- 1L
      base_row$enough_shared        <- TRUE
      base_row$effect_correlation   <- NA_real_  # needs >= 3 alleles
      base_row$direction_agreement  <- dir_agree
      base_row$directionally_concordant <- dir_concordant
      base_row$meta_effect          <- round(b_ivw,   6)
      base_row$meta_SE              <- round(se_ivw,  6)
      base_row$meta_z               <- round(meta_z_v, 4)
      base_row$meta_p               <- meta_p_v
      base_row$Q_stat               <- Q_v
      base_row$Q_df                 <- Q_df_v
      base_row$Q_p                  <- Q_p_v
      base_row$I2                   <- I2_v
      base_row$replicated           <- replicated

      conc_rows[[length(conc_rows) + 1L]] <- base_row

      shared_rows[[length(shared_rows) + 1L]] <- data.frame(
        block_id       = bn,
        CHR            = blk_chr,
        start_bp       = blk_start,
        end_bp         = blk_end,
        trait          = tr,
        lead_snp       = paste(r1$lead_snp[1], r2$lead_snp[1], sep = " / "),
        lead_snp_pop1  = r1$lead_snp[1],
        lead_snp_pop2  = r2$lead_snp[1],
        effect_pop1    = round(e1,  6),
        SE_pop1        = round(se1, 6),
        p_wald_pop1    = r1$lead_p[1],
        effect_pop2    = round(e2,  6),
        SE_pop2        = round(se2, 6),
        p_wald_pop2    = r2$lead_p[1],
        direction_agree = sign(e1) == sign(e2),
        ivw_effect     = round(b_ivw,  6),
        ivw_SE         = round(se_ivw, 6),
        stringsAsFactors = FALSE
      )
    }
  }

  conc_df   <- if (length(conc_rows))  do.call(rbind, conc_rows)  else data.frame()
  shared_df <- if (length(shared_rows)) do.call(rbind, shared_rows) else data.frame()

  if (nrow(conc_df) > 0) {
    conc_df <- conc_df[order(conc_df$CHR, conc_df$start_bp,
                             ifelse(is.na(conc_df$meta_p), 1, conc_df$meta_p)), ]
    rownames(conc_df) <- NULL
  }
  if (nrow(shared_df) > 0) rownames(shared_df) <- NULL

  n_both <- sum(!is.na(conc_df$meta_p))
  n_rep  <- sum(conc_df$replicated, na.rm = TRUE)
  .log("Done. Blocks in both pops: ", n_both,
       " | Replicated (dir concordant + meta_p <= 0.05): ", n_rep)

  structure(
    list(
      concordance           = conc_df,
      shared_alleles        = shared_df,
      pop1_name             = pop1_name,
      pop2_name             = pop2_name,
      traits                = traits,
      direction_threshold   = direction_threshold,
      boundary_overlap_warn = boundary_overlap_warn
    ),
    class = c("LDxBlocks_effect_concordance", "list")
  )
}
