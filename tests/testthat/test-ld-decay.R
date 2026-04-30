## tests/testthat/test-ld-decay.R
## Tests for compute_ld_decay() and plot_ld_decay().
## Covers: input validation, all sampling modes, backend types,
## parametric threshold, censored flag, model fitting, and
## integration with define_qtl_regions().

library(testthat)
library(LDxBlocks)

data(ldx_geno,     package = "LDxBlocks")
data(ldx_snp_info, package = "LDxBlocks")
data(ldx_blocks,   package = "LDxBlocks")
data(ldx_gwas,     package = "LDxBlocks")

# ── Input validation ──────────────────────────────────────────────────────────

test_that("compute_ld_decay: rejects invalid r2_threshold values", {
  expect_error(
    compute_ld_decay(ldx_geno, ldx_snp_info, r2_threshold = -0.1,
                     verbose = FALSE),
    "\\[0, 1\\]"
  )
  expect_error(
    compute_ld_decay(ldx_geno, ldx_snp_info, r2_threshold = 1.5,
                     verbose = FALSE),
    "\\[0, 1\\]"
  )
  expect_error(
    compute_ld_decay(ldx_geno, ldx_snp_info, r2_threshold = NA_real_,
                     verbose = FALSE),
    "\\[0, 1\\]"
  )
  expect_error(
    compute_ld_decay(ldx_geno, ldx_snp_info, r2_threshold = "invalid",
                     verbose = FALSE),
    "parametric"
  )
})

test_that("compute_ld_decay: accepts NULL r2_threshold (no threshold)", {
  decay <- compute_ld_decay(ldx_geno, ldx_snp_info,
                            r2_threshold = NULL, fit_model = "none",
                            n_pairs = 500L, verbose = FALSE)
  expect_null(decay$critical_r2)
  expect_null(decay$decay_dist)
  expect_s3_class(decay, "LDxBlocks_decay")
})

test_that("compute_ld_decay: rejects geno without snp_info when matrix", {
  expect_error(
    compute_ld_decay(ldx_geno, snp_info = NULL, verbose = FALSE),
    "snp_info"
  )
})

# ── Return structure ──────────────────────────────────────────────────────────

test_that("compute_ld_decay: returns LDxBlocks_decay object with required fields", {
  decay <- compute_ld_decay(ldx_geno, ldx_snp_info,
                            sampling     = "random",
                            r2_threshold = 0.1,
                            fit_model    = "loess",
                            n_pairs      = 1000L,
                            verbose      = FALSE)
  expect_s3_class(decay, "LDxBlocks_decay")
  req <- c("pairs","decay_curve","decay_dist","decay_dist_genome",
           "critical_r2","critical_r2_fixed","critical_r2_param",
           "unlinked_r2","model_params","n_pairs_used","method",
           "sampling","pctile","call")
  expect_true(all(req %in% names(decay)))
})

test_that("compute_ld_decay: pairs data.frame has CHR, dist_bp, r2 columns", {
  decay <- compute_ld_decay(ldx_geno, ldx_snp_info,
                            n_pairs = 500L, r2_threshold = NULL,
                            fit_model = "none", verbose = FALSE)
  expect_true(all(c("CHR","dist_bp","r2") %in% names(decay$pairs)))
  expect_true(all(decay$pairs$r2 >= 0 & decay$pairs$r2 <= 1 + 1e-8))
  expect_true(all(decay$pairs$dist_bp >= 0))
})

test_that("compute_ld_decay: decay_curve has expected columns", {
  decay <- compute_ld_decay(ldx_geno, ldx_snp_info,
                            n_pairs = 500L, r2_threshold = NULL,
                            fit_model = "loess", verbose = FALSE)
  expect_true(all(c("CHR","dist_bp","r2_mean","r2_median","n_pairs")
                  %in% names(decay$decay_curve)))
  expect_true("r2_loess" %in% names(decay$decay_curve))
})

test_that("compute_ld_decay: decay_dist has censored column", {
  decay <- compute_ld_decay(ldx_geno, ldx_snp_info,
                            r2_threshold = 0.1, n_pairs = 500L,
                            fit_model = "loess", verbose = FALSE)
  expect_true("censored" %in% names(decay$decay_dist))
  expect_type(decay$decay_dist$censored, "logical")
})

test_that("compute_ld_decay: n_pairs_used named vector matches chromosomes", {
  decay <- compute_ld_decay(ldx_geno, ldx_snp_info,
                            n_pairs = 500L, r2_threshold = NULL,
                            fit_model = "none", verbose = FALSE)
  expect_equal(sort(names(decay$n_pairs_used)),
               sort(unique(ldx_snp_info$CHR)))
})

# ── Sampling modes ────────────────────────────────────────────────────────────

test_that("compute_ld_decay: random sampling produces non-empty pairs", {
  decay <- compute_ld_decay(ldx_geno, ldx_snp_info,
                            sampling = "random", n_pairs = 1000L,
                            r2_threshold = NULL, fit_model = "none",
                            verbose = FALSE)
  expect_gt(nrow(decay$pairs), 0L)
  expect_equal(decay$sampling, "random")
})

test_that("compute_ld_decay: sliding_window sampling produces non-empty pairs", {
  decay <- compute_ld_decay(ldx_geno, ldx_snp_info,
                            sampling = "sliding_window",
                            window_snps = 10L, n_pairs = 500L,
                            r2_threshold = NULL, fit_model = "none",
                            verbose = FALSE)
  expect_gt(nrow(decay$pairs), 0L)
  expect_equal(decay$sampling, "sliding_window")
})

test_that("compute_ld_decay: both sampling mode uses sliding for pairs", {
  decay <- compute_ld_decay(ldx_geno, ldx_snp_info,
                            sampling = "both", window_snps = 10L,
                            n_pairs = 500L, r2_threshold = NULL,
                            fit_model = "none", verbose = FALSE)
  expect_equal(decay$sampling, "both")
  expect_gt(nrow(decay$pairs), 0L)
})

test_that("compute_ld_decay: chr parameter restricts to named chromosomes", {
  decay <- compute_ld_decay(ldx_geno, ldx_snp_info,
                            chr = "1", n_pairs = 500L,
                            r2_threshold = NULL, fit_model = "none",
                            verbose = FALSE)
  expect_true(all(decay$pairs$CHR == "1"))
  if (!is.null(decay$decay_curve))
    expect_true(all(decay$decay_curve$CHR == "1"))
})

# ── SNP position sorting ──────────────────────────────────────────────────────

test_that("compute_ld_decay: produces valid pairs when snp_info is unsorted", {
  # Shuffle snp_info rows -- function must sort internally
  set.seed(7)
  shuffled <- ldx_snp_info[sample(nrow(ldx_snp_info)), ]
  decay <- compute_ld_decay(ldx_geno[, match(shuffled$SNP, colnames(ldx_geno))],
                            shuffled, n_pairs = 500L,
                            r2_threshold = NULL, fit_model = "none",
                            verbose = FALSE)
  # All dist_bp must be non-negative (would fail if positions are mis-ordered)
  expect_true(all(decay$pairs$dist_bp >= 0L))
})

# ── Threshold types ───────────────────────────────────────────────────────────

test_that("compute_ld_decay: fixed numeric threshold populates decay_dist", {
  decay <- compute_ld_decay(ldx_geno, ldx_snp_info,
                            r2_threshold = 0.1, n_pairs = 1000L,
                            fit_model = "loess", verbose = FALSE)
  expect_equal(decay$critical_r2_fixed, 0.1)
  expect_null(decay$critical_r2_param)
  if (!is.null(decay$decay_dist))
    expect_true(all(decay$decay_dist$threshold_used == 0.1))
})

test_that("compute_ld_decay: parametric threshold is in [0,1]", {
  decay <- compute_ld_decay(ldx_geno, ldx_snp_info,
                            r2_threshold = "parametric",
                            n_unlinked = 1000L, n_pairs = 500L,
                            fit_model = "none", verbose = FALSE)
  expect_false(is.null(decay$critical_r2_param))
  expect_true(decay$critical_r2_param >= 0 & decay$critical_r2_param <= 1)
  expect_false(is.null(decay$unlinked_r2))
  expect_true(length(decay$unlinked_r2) > 0L)
})

test_that("compute_ld_decay: 'both' returns both thresholds", {
  decay <- compute_ld_decay(ldx_geno, ldx_snp_info,
                            r2_threshold = "both", n_pairs = 500L,
                            n_unlinked = 500L, fit_model = "none",
                            verbose = FALSE)
  expect_equal(decay$critical_r2_fixed, 0.1)
  expect_false(is.null(decay$critical_r2_param))
  # Active threshold should be the parametric one
  expect_equal(decay$critical_r2, decay$critical_r2_param)
})

# ── Model fitting ─────────────────────────────────────────────────────────────

test_that("compute_ld_decay: fit_model='none' skips smoothing", {
  decay <- compute_ld_decay(ldx_geno, ldx_snp_info,
                            n_pairs = 500L, fit_model = "none",
                            r2_threshold = NULL, verbose = FALSE)
  if (!is.null(decay$decay_curve))
    expect_false("r2_loess" %in% names(decay$decay_curve))
})

test_that("compute_ld_decay: fit_model='loess' produces r2_loess column", {
  decay <- compute_ld_decay(ldx_geno, ldx_snp_info,
                            n_pairs = 1000L, fit_model = "loess",
                            r2_threshold = NULL, verbose = FALSE)
  if (!is.null(decay$decay_curve))
    expect_true("r2_loess" %in% names(decay$decay_curve))
})

test_that("compute_ld_decay: fit_model='nonlinear' attempts Hill-Weir fit", {
  decay <- suppressWarnings(
    compute_ld_decay(ldx_geno, ldx_snp_info,
                     n_pairs = 2000L, fit_model = "nonlinear",
                     r2_threshold = NULL, verbose = FALSE)
  )
  # Either model_params is populated or a warning was emitted (small data)
  # -- either outcome is acceptable; just check no error
  expect_s3_class(decay, "LDxBlocks_decay")
})

# ── Decay distance ────────────────────────────────────────────────────────────

test_that("compute_ld_decay: decay_dist_bp is positive for each chromosome", {
  decay <- compute_ld_decay(ldx_geno, ldx_snp_info,
                            r2_threshold = 0.1, n_pairs = 1000L,
                            fit_model = "loess", verbose = FALSE)
  if (!is.null(decay$decay_dist))
    expect_true(all(decay$decay_dist$decay_dist_bp > 0))
})

test_that("compute_ld_decay: decay_dist_genome is median of per-chr distances", {
  decay <- compute_ld_decay(ldx_geno, ldx_snp_info,
                            r2_threshold = 0.1, n_pairs = 1000L,
                            fit_model = "loess", verbose = FALSE)
  if (!is.null(decay$decay_dist) && nrow(decay$decay_dist) > 0L) {
    expected_median <- stats::median(decay$decay_dist$decay_dist_bp)
    expect_equal(decay$decay_dist_genome, expected_median)
  }
})

test_that("compute_ld_decay: censored=TRUE when threshold never crossed", {
  # Use a very high threshold (1.0) that will never be crossed
  decay <- suppressWarnings(
    compute_ld_decay(ldx_geno, ldx_snp_info,
                     r2_threshold = 0.999, n_pairs = 1000L,
                     fit_model = "loess", verbose = FALSE)
  )
  if (!is.null(decay$decay_dist))
    expect_true(any(decay$decay_dist$censored))
})

# ── Backend types ─────────────────────────────────────────────────────────────

test_that("compute_ld_decay: accepts LDxBlocks_backend (matrix type)", {
  be <- read_geno(ldx_geno, format = "matrix", snp_info = ldx_snp_info)
  on.exit(close_backend(be))
  decay <- compute_ld_decay(be, n_pairs = 500L,
                            r2_threshold = NULL, fit_model = "none",
                            verbose = FALSE)
  expect_s3_class(decay, "LDxBlocks_decay")
  expect_gt(nrow(decay$pairs), 0L)
})

test_that("compute_ld_decay: matrix and backend produce consistent pairs", {
  be <- read_geno(ldx_geno, format = "matrix", snp_info = ldx_snp_info)
  on.exit(close_backend(be))
  set.seed(1)
  d_mat <- compute_ld_decay(ldx_geno, ldx_snp_info,
                            n_pairs = 300L, r2_threshold = NULL,
                            fit_model = "none", seed = 1L, verbose = FALSE)
  set.seed(1)
  d_be  <- compute_ld_decay(be, n_pairs = 300L, r2_threshold = NULL,
                            fit_model = "none", seed = 1L, verbose = FALSE)
  # Same number of chromosomes processed
  expect_equal(sum(d_mat$n_pairs_used > 0),
               sum(d_be$n_pairs_used  > 0))
})

test_that("compute_ld_decay: accepts bigmemory backend", {
  skip_if_not_installed("bigmemory")
  be_mat <- read_geno(ldx_geno, format = "matrix", snp_info = ldx_snp_info)
  bm_path <- file.path(tempdir(), paste0("ldxbm_decay_", as.integer(Sys.time())))
  dir.create(bm_path, showWarnings = FALSE)
  on.exit({ unlink(bm_path, recursive = TRUE); close_backend(be_mat) },
          add = TRUE)
  be_bm <- read_geno_bigmemory(be_mat,
                               backingfile = tempfile("bm_", tmpdir = bm_path),
                               backingpath = bm_path, type = "char",
                               verbose = FALSE)
  on.exit(close_backend(be_bm), add = TRUE)
  decay <- compute_ld_decay(be_bm, n_pairs = 300L, r2_threshold = NULL,
                            fit_model = "none", verbose = FALSE)
  expect_s3_class(decay, "LDxBlocks_decay")
  expect_gt(nrow(decay$pairs), 0L)
})

# ── print method ─────────────────────────────────────────────────────────────

test_that("print.LDxBlocks_decay: prints without error", {
  decay <- compute_ld_decay(ldx_geno, ldx_snp_info,
                            r2_threshold = 0.1, n_pairs = 500L,
                            fit_model = "loess", verbose = FALSE)
  expect_output(print(decay), "LDxBlocks LD Decay Analysis")
})

# ── plot_ld_decay ─────────────────────────────────────────────────────────────

test_that("plot_ld_decay: returns ggplot when fit_model produces curve", {
  skip_if_not_installed("ggplot2")
  decay <- compute_ld_decay(ldx_geno, ldx_snp_info,
                            r2_threshold = 0.1, n_pairs = 1000L,
                            fit_model = "loess", verbose = FALSE)
  p <- plot_ld_decay(decay)
  expect_s3_class(p, "ggplot")
})

test_that("plot_ld_decay: errors on non-LDxBlocks_decay input", {
  skip_if_not_installed("ggplot2")
  expect_error(plot_ld_decay(list(a = 1)), "LDxBlocks_decay")
})

test_that("plot_ld_decay: errors when fit_model='none' (no curve)", {
  skip_if_not_installed("ggplot2")
  decay <- compute_ld_decay(ldx_geno, ldx_snp_info,
                            r2_threshold = NULL, n_pairs = 500L,
                            fit_model = "none", verbose = FALSE)
  expect_error(plot_ld_decay(decay), "fit_model")
})

test_that("plot_ld_decay: facet=TRUE produces facet_wrap", {
  skip_if_not_installed("ggplot2")
  decay <- compute_ld_decay(ldx_geno, ldx_snp_info,
                            r2_threshold = 0.1, n_pairs = 1000L,
                            fit_model = "loess", verbose = FALSE)
  p <- plot_ld_decay(decay, facet = TRUE)
  # Check FacetWrap layer is present
  expect_true(inherits(p$facet, "FacetWrap"))
})

# ── define_qtl_regions integration ───────────────────────────────────────────

test_that("define_qtl_regions: accepts LDxBlocks_decay object in ld_decay", {
  decay <- compute_ld_decay(ldx_geno, ldx_snp_info,
                            r2_threshold = 0.1, n_pairs = 1000L,
                            fit_model = "loess", verbose = FALSE)
  qtl <- suppressWarnings(
    define_qtl_regions(ldx_gwas, ldx_blocks, ldx_snp_info,
                       ld_decay = decay, p_threshold = NULL,
                       trait_col = "trait")
  )
  expect_s3_class(qtl, "data.frame")
  # With ld_decay, output should have candidate region columns
  expect_true("candidate_region_start" %in% names(qtl))
  expect_true("candidate_region_end"   %in% names(qtl))
  expect_true("candidate_region_size_kb" %in% names(qtl))
})

test_that("define_qtl_regions with ld_decay: candidate regions are wider than blocks", {
  decay <- compute_ld_decay(ldx_geno, ldx_snp_info,
                            r2_threshold = 0.1, n_pairs = 1000L,
                            fit_model = "loess", verbose = FALSE)
  qtl_decay <- suppressWarnings(
    define_qtl_regions(ldx_gwas, ldx_blocks, ldx_snp_info,
                       ld_decay = decay, p_threshold = NULL)
  )
  qtl_plain <- define_qtl_regions(ldx_gwas, ldx_blocks, ldx_snp_info,
                                  p_threshold = NULL)
  # Candidate regions should be at least as wide as block boundaries
  if (nrow(qtl_decay) > 0L && !all(is.na(qtl_decay$candidate_region_start))) {
    sizes <- qtl_decay$candidate_region_end - qtl_decay$candidate_region_start
    block_sizes <- qtl_plain$end_bp - qtl_plain$start_bp
    # Check against matching blocks by block_id
    common_ids <- intersect(qtl_decay$block_id, qtl_plain$block_id)
    if (length(common_ids) > 0L) {
      for (bid in common_ids) {
        cr_size <- qtl_decay$candidate_region_end[qtl_decay$block_id == bid] -
          qtl_decay$candidate_region_start[qtl_decay$block_id == bid]
        blk_size <- qtl_plain$end_bp[qtl_plain$block_id == bid] -
          qtl_plain$start_bp[qtl_plain$block_id == bid]
        expect_true(cr_size >= blk_size - 1L,
                    label = paste("candidate region >= block for", bid))
      }
    }
  }
})

test_that("define_qtl_regions: without ld_decay candidate_region_start uses block boundary", {
  qtl <- define_qtl_regions(ldx_gwas, ldx_blocks, ldx_snp_info,
                            p_threshold = NULL)
  # When ld_decay = NULL, fallback is block start/end boundary (not NA)
  expect_true("candidate_region_start" %in% names(qtl))
  expect_true(all(!is.na(qtl$candidate_region_start)))
  expect_true(all(!is.na(qtl$candidate_region_end)))
})
