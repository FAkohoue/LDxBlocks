## tests/testthat/test-prediction.R
## -----------------------------------------------------------------------------
## Thorough tests for run_haplotype_prediction():
##   - All four blues input formats (named vector, ST df, MT df, named list)
##   - Single-trait: full return structure verification
##   - Multi-trait: all three importance_rule values, auto blue_cols, list format
##   - Edge cases: id mismatch, NA blues, partial overlap
##   - run_haplotype_prediction via ldx_blues canonical example
## -----------------------------------------------------------------------------

library(testthat)
library(LDxBlocks)

data(ldx_geno,     package = "LDxBlocks")
data(ldx_snp_info, package = "LDxBlocks")
data(ldx_blocks,   package = "LDxBlocks")
data(ldx_blues,    package = "LDxBlocks")

# Shared synthetic blues (reproducible)
set.seed(7L)
blues_vec <- structure(
  round(ldx_geno[, 1L] * 0.4 + rnorm(120L), 4L),
  names = rownames(ldx_geno)
)
blues_st_df <- data.frame(
  Genotype = rownames(ldx_geno),
  YLD_BLUE = blues_vec,
  stringsAsFactors = FALSE
)
blues_mt_df <- data.frame(
  id  = rownames(ldx_geno),
  YLD = blues_vec,
  RES = round(ldx_geno[, 2L] * 0.3 + rnorm(120L), 4L),
  stringsAsFactors = FALSE
)

# -- Format 1: Named numeric vector --------------------------------------------

test_that("ST named vector: runs without error", {
  expect_no_error(
    run_haplotype_prediction(ldx_geno, ldx_snp_info, ldx_blocks,
                             blues = blues_vec, verbose = FALSE)
  )
})

test_that("ST named vector: return structure has all required keys", {
  res <- run_haplotype_prediction(ldx_geno, ldx_snp_info, ldx_blocks,
                                  blues = blues_vec, verbose = FALSE)
  required <- c("blocks", "diversity", "hap_matrix", "haplotypes",
                "n_blocks", "n_hap_columns", "n_traits", "traits",
                "solver_used", "gebv", "snp_effects", "local_gebv",
                "block_importance", "G", "n_train", "n_predict")
  expect_true(all(required %in% names(res)),
              info = paste("Missing:", paste(setdiff(required, names(res)), collapse=", ")))
})

test_that("ST named vector: n_traits == 1, traits is length-1 character", {
  res <- run_haplotype_prediction(ldx_geno, ldx_snp_info, ldx_blocks,
                                  blues = blues_vec, verbose = FALSE)
  expect_equal(res$n_traits, 1L)
  expect_type(res$traits, "character")
  expect_equal(length(res$traits), 1L)
})

test_that("ST named vector: solver_used is 'rrBLUP' (single-trait always uses rrBLUP)", {
  res <- run_haplotype_prediction(ldx_geno, ldx_snp_info, ldx_blocks,
                                  blues = blues_vec, verbose = FALSE)
  expect_equal(res$solver_used, "rrBLUP")
})

test_that("ST named vector: gebv is a named numeric vector covering training individuals", {
  res <- run_haplotype_prediction(ldx_geno, ldx_snp_info, ldx_blocks,
                                  blues = blues_vec, verbose = FALSE)
  expect_type(res$gebv, "double")
  expect_false(is.null(names(res$gebv)))
  expect_equal(res$n_train, sum(!is.na(blues_vec)))
})

test_that("ST named vector: local_gebv is a matrix (individuals x blocks)", {
  res <- run_haplotype_prediction(ldx_geno, ldx_snp_info, ldx_blocks,
                                  blues = blues_vec, verbose = FALSE)
  expect_true(is.matrix(res$local_gebv))
  expect_equal(ncol(res$local_gebv), res$n_blocks)
})

test_that("ST named vector: block_importance has required columns and correct row count", {
  res <- run_haplotype_prediction(ldx_geno, ldx_snp_info, ldx_blocks,
                                  blues = blues_vec, verbose = FALSE)
  bi <- res$block_importance
  req_cols <- c("block_id", "CHR", "start_bp", "end_bp", "n_snps",
                "var_local_gebv", "var_scaled", "important")
  expect_true(all(req_cols %in% names(bi)),
              info = paste("Missing:", paste(setdiff(req_cols, names(bi)), collapse=", ")))
  expect_equal(nrow(bi), res$n_blocks)
})

test_that("ST named vector: var_scaled is in [0,1]", {
  res <- run_haplotype_prediction(ldx_geno, ldx_snp_info, ldx_blocks,
                                  blues = blues_vec, verbose = FALSE)
  vs <- res$block_importance$var_scaled
  expect_true(all(vs >= 0 & vs <= 1 + 1e-8, na.rm = TRUE))
})

test_that("ST named vector: important column is logical", {
  res <- run_haplotype_prediction(ldx_geno, ldx_snp_info, ldx_blocks,
                                  blues = blues_vec, verbose = FALSE)
  expect_type(res$block_importance$important, "logical")
})

test_that("ST named vector: snp_effects named by SNP IDs", {
  res <- run_haplotype_prediction(ldx_geno, ldx_snp_info, ldx_blocks,
                                  blues = blues_vec, verbose = FALSE)
  expect_type(res$snp_effects, "double")
  expect_false(is.null(names(res$snp_effects)))
  # All SNP names should be in snp_info
  expect_true(all(names(res$snp_effects) %in% ldx_snp_info$SNP))
})

test_that("ST named vector: G is a square symmetric matrix", {
  res <- run_haplotype_prediction(ldx_geno, ldx_snp_info, ldx_blocks,
                                  blues = blues_vec, verbose = FALSE)
  expect_true(is.matrix(res$G))
  expect_equal(nrow(res$G), ncol(res$G))
  expect_true(isSymmetric(res$G, tol = 1e-8))
})

# -- Format 2: Data frame, single trait, explicit id_col + blue_col ------------

test_that("ST data.frame (id_col + blue_col): runs without error", {
  expect_no_error(
    run_haplotype_prediction(ldx_geno, ldx_snp_info, ldx_blocks,
                             blues    = blues_st_df,
                             id_col   = "Genotype",
                             blue_col = "YLD_BLUE",
                             verbose  = FALSE)
  )
})

test_that("ST data.frame: n_traits == 1 and solver is rrBLUP", {
  res <- run_haplotype_prediction(ldx_geno, ldx_snp_info, ldx_blocks,
                                  blues    = blues_st_df,
                                  id_col   = "Genotype",
                                  blue_col = "YLD_BLUE",
                                  verbose  = FALSE)
  expect_equal(res$n_traits, 1L)
  expect_equal(res$solver_used, "rrBLUP")
})

test_that("ST data.frame: gebv matches output from named vector (same data)", {
  res_vec <- run_haplotype_prediction(ldx_geno, ldx_snp_info, ldx_blocks,
                                      blues = blues_vec, verbose = FALSE)
  res_df  <- run_haplotype_prediction(ldx_geno, ldx_snp_info, ldx_blocks,
                                      blues    = blues_st_df,
                                      id_col   = "Genotype",
                                      blue_col = "YLD_BLUE",
                                      verbose  = FALSE)
  # GEBVs on same data should be identical regardless of input format
  common <- intersect(names(res_vec$gebv), names(res_df$gebv))
  expect_gt(length(common), 0L)
  expect_equal(res_vec$gebv[common], res_df$gebv[common], tolerance = 1e-6)
})

# -- Format 2b: Data frame, auto-detect single numeric column ------------------

test_that("ST data.frame auto-detect: single numeric column detected without blue_col", {
  blues_auto <- data.frame(
    id   = rownames(ldx_geno),
    YLD  = blues_vec,
    stringsAsFactors = FALSE
  )
  # blue_col default is "blue" - not present, so auto-detect should kick in
  expect_no_error(
    res <- run_haplotype_prediction(ldx_geno, ldx_snp_info, ldx_blocks,
                                    blues   = blues_auto,
                                    id_col  = "id",
                                    verbose = FALSE)
  )
  expect_equal(res$n_traits, 1L)
})

# -- Format 3: Data frame, multiple traits, explicit blue_cols -----------------

test_that("MT data.frame (blue_cols): runs without error", {
  expect_no_error(
    run_haplotype_prediction(ldx_geno, ldx_snp_info, ldx_blocks,
                             blues     = blues_mt_df,
                             id_col    = "id",
                             blue_cols = c("YLD", "RES"),
                             verbose   = FALSE)
  )
})

test_that("MT data.frame: n_traits == 2, traits has correct names", {
  res <- run_haplotype_prediction(ldx_geno, ldx_snp_info, ldx_blocks,
                                  blues     = blues_mt_df,
                                  id_col    = "id",
                                  blue_cols = c("YLD", "RES"),
                                  verbose   = FALSE)
  expect_equal(res$n_traits, 2L)
  expect_equal(sort(res$traits), sort(c("YLD", "RES")))
})

test_that("MT data.frame: solver_used is 'rrBLUP'", {
  res <- run_haplotype_prediction(ldx_geno, ldx_snp_info, ldx_blocks,
                                  blues     = blues_mt_df,
                                  id_col    = "id",
                                  blue_cols = c("YLD", "RES"),
                                  verbose   = FALSE)
  expect_equal(res$solver_used, "rrBLUP")
})

test_that("MT data.frame: block_importance has per-trait columns", {
  res <- run_haplotype_prediction(ldx_geno, ldx_snp_info, ldx_blocks,
                                  blues     = blues_mt_df,
                                  id_col    = "id",
                                  blue_cols = c("YLD", "RES"),
                                  verbose   = FALSE)
  bi <- res$block_importance
  # Per-trait scaled variance columns
  expect_true("var_scaled_YLD" %in% names(bi))
  expect_true("var_scaled_RES" %in% names(bi))
  # Per-trait importance flags
  expect_true("important_YLD" %in% names(bi))
  expect_true("important_RES" %in% names(bi))
  # Cross-trait aggregates
  expect_true("var_scaled_mean"    %in% names(bi))
  expect_true("n_traits_important" %in% names(bi))
  expect_true("important_any"      %in% names(bi))
  expect_true("important_all"      %in% names(bi))
})

test_that("MT data.frame: n_traits_important in 0..n_traits", {
  res <- run_haplotype_prediction(ldx_geno, ldx_snp_info, ldx_blocks,
                                  blues     = blues_mt_df,
                                  id_col    = "id",
                                  blue_cols = c("YLD", "RES"),
                                  verbose   = FALSE)
  n <- res$block_importance$n_traits_important
  expect_true(all(n >= 0L & n <= 2L))
})

test_that("MT data.frame: important_any TRUE iff n_traits_important >= 1", {
  res <- run_haplotype_prediction(ldx_geno, ldx_snp_info, ldx_blocks,
                                  blues           = blues_mt_df,
                                  id_col          = "id",
                                  blue_cols       = c("YLD", "RES"),
                                  importance_rule = "any",
                                  verbose         = FALSE)
  bi <- res$block_importance
  expect_equal(bi$important, bi$n_traits_important >= 1L)
  expect_equal(bi$important, bi$important_any)
})

test_that("MT data.frame: importance_rule='all' sets important = important_all", {
  res <- run_haplotype_prediction(ldx_geno, ldx_snp_info, ldx_blocks,
                                  blues           = blues_mt_df,
                                  id_col          = "id",
                                  blue_cols       = c("YLD", "RES"),
                                  importance_rule = "all",
                                  verbose         = FALSE)
  bi <- res$block_importance
  expect_equal(bi$important, bi$important_all)
})

test_that("MT data.frame: importance_rule='mean' sets important via var_scaled_mean >= 0.9", {
  res <- run_haplotype_prediction(ldx_geno, ldx_snp_info, ldx_blocks,
                                  blues           = blues_mt_df,
                                  id_col          = "id",
                                  blue_cols       = c("YLD", "RES"),
                                  importance_rule = "mean",
                                  verbose         = FALSE)
  bi <- res$block_importance
  expect_equal(bi$important, bi$var_scaled_mean >= 0.9)
})

test_that("MT data.frame: block_importance_list is a named list (one per trait)", {
  res <- run_haplotype_prediction(ldx_geno, ldx_snp_info, ldx_blocks,
                                  blues     = blues_mt_df,
                                  id_col    = "id",
                                  blue_cols = c("YLD", "RES"),
                                  verbose   = FALSE)
  expect_type(res$block_importance_list, "list")
  expect_equal(sort(names(res$block_importance_list)), sort(c("YLD", "RES")))
  for (tr in c("YLD", "RES")) {
    bi_tr <- res$block_importance_list[[tr]]
    expect_true("var_scaled" %in% names(bi_tr))
    expect_true("important"  %in% names(bi_tr))
  }
})

test_that("MT data.frame: gebv and snp_effects are named lists (one per trait)", {
  res <- run_haplotype_prediction(ldx_geno, ldx_snp_info, ldx_blocks,
                                  blues     = blues_mt_df,
                                  id_col    = "id",
                                  blue_cols = c("YLD", "RES"),
                                  verbose   = FALSE)
  expect_type(res$gebv,        "list")
  expect_type(res$snp_effects, "list")
  expect_equal(sort(names(res$gebv)),        sort(c("YLD", "RES")))
  expect_equal(sort(names(res$snp_effects)), sort(c("YLD", "RES")))
})

# -- Format 3b: Data frame, auto-detect blue_cols (blue_cols = NULL) -----------

test_that("MT data.frame auto-detect blue_cols=NULL: finds YLD and RES", {
  res <- run_haplotype_prediction(ldx_geno, ldx_snp_info, ldx_blocks,
                                  blues      = blues_mt_df,
                                  id_col     = "id",
                                  blue_cols  = NULL,   # auto-detect
                                  verbose    = FALSE)
  expect_equal(res$n_traits, 2L)
  expect_equal(sort(res$traits), sort(c("YLD", "RES")))
})

# -- Format 4: Named list of named numeric vectors -----------------------------

test_that("MT named list: runs without error", {
  blues_list <- list(
    YLD = blues_vec,
    RES = setNames(blues_mt_df$RES, rownames(ldx_geno))
  )
  expect_no_error(
    run_haplotype_prediction(ldx_geno, ldx_snp_info, ldx_blocks,
                             blues   = blues_list,
                             verbose = FALSE)
  )
})

test_that("MT named list: n_traits == 2, traits correct", {
  blues_list <- list(
    YLD = blues_vec,
    RES = setNames(blues_mt_df$RES, rownames(ldx_geno))
  )
  res <- run_haplotype_prediction(ldx_geno, ldx_snp_info, ldx_blocks,
                                  blues   = blues_list,
                                  verbose = FALSE)
  expect_equal(res$n_traits, 2L)
  expect_equal(sort(res$traits), sort(c("YLD", "RES")))
})

test_that("MT named list: partial overlap per trait handled gracefully", {
  # RES measured on only first 80 individuals
  blues_partial <- list(
    YLD = blues_vec,                             # 120 individuals
    RES = blues_vec[1:80]                        # 80 individuals (subset)
  )
  expect_no_error(
    res <- run_haplotype_prediction(ldx_geno, ldx_snp_info, ldx_blocks,
                                    blues   = blues_partial,
                                    verbose = FALSE)
  )
  expect_equal(res$n_traits, 2L)
  # YLD training count should be 120, RES should be 80
  expect_equal(res$n_train[["YLD"]], 120L)
  expect_equal(res$n_train[["RES"]], 80L)
})

test_that("MT named list: unnamed list throws an error", {
  blues_unnamed <- list(blues_vec, blues_vec)  # no names
  expect_error(
    run_haplotype_prediction(ldx_geno, ldx_snp_info, ldx_blocks,
                             blues = blues_unnamed, verbose = FALSE),
    "named"
  )
})

# -- Canonical example: ldx_blues dataset -------------------------------------

test_that("ldx_blues single-trait (YLD): complete pipeline runs end-to-end", {
  res <- run_haplotype_prediction(
    ldx_geno, ldx_snp_info, ldx_blocks,
    blues    = ldx_blues,
    id_col   = "id",
    blue_col = "YLD",
    verbose  = FALSE
  )
  expect_equal(res$n_traits, 1L)
  expect_equal(res$n_train, 120L)
  expect_true(all(c("var_scaled", "important") %in% names(res$block_importance)))
  expect_equal(nrow(res$block_importance), res$n_blocks)
})

test_that("ldx_blues multi-trait (YLD + RES): complete pipeline runs end-to-end", {
  res <- run_haplotype_prediction(
    ldx_geno, ldx_snp_info, ldx_blocks,
    blues           = ldx_blues,
    id_col          = "id",
    blue_cols       = c("YLD", "RES"),
    importance_rule = "any",
    verbose         = FALSE
  )
  expect_equal(res$n_traits, 2L)
  expect_equal(sort(res$traits), sort(c("YLD", "RES")))
  expect_true("var_scaled_YLD"    %in% names(res$block_importance))
  expect_true("var_scaled_RES"    %in% names(res$block_importance))
  expect_true("var_scaled_mean"   %in% names(res$block_importance))
  expect_true("important_any"     %in% names(res$block_importance))
  expect_true("important_all"     %in% names(res$block_importance))
  expect_equal(nrow(res$block_importance), res$n_blocks)
})

# -- Edge cases ----------------------------------------------------------------

test_that("ID mismatch (zero overlap): throws informative error", {
  blues_bad <- c(NOBODY1 = 1.0, NOBODY2 = 2.0, NOBODY3 = 3.0)
  expect_error(
    run_haplotype_prediction(ldx_geno, ldx_snp_info, ldx_blocks,
                             blues = blues_bad, verbose = FALSE)
  )
})

test_that("Partial ID overlap: proceeds with common individuals only", {
  # Use only first 60 individuals - valid partial overlap
  blues_partial <- blues_vec[1:60]
  expect_no_error(
    res <- run_haplotype_prediction(ldx_geno, ldx_snp_info, ldx_blocks,
                                    blues = blues_partial, verbose = FALSE)
  )
  expect_equal(res$n_train, 60L)
})

test_that("NA blues values: treated as prediction candidates (n_predict > 0)", {
  blues_na      <- blues_vec
  blues_na[1:5] <- NA  # 5 individuals become prediction candidates
  res <- run_haplotype_prediction(ldx_geno, ldx_snp_info, ldx_blocks,
                                  blues = blues_na, verbose = FALSE)
  expect_equal(res$n_train,   115L)
  expect_equal(res$n_predict, 5L)
  # GEBVs should exist for all 120 individuals (training + prediction)
  expect_equal(length(res$gebv), 120L)
})

test_that("wrong id_col name throws informative error", {
  expect_error(
    run_haplotype_prediction(ldx_geno, ldx_snp_info, ldx_blocks,
                             blues  = blues_st_df,
                             id_col = "WRONG_COL",
                             verbose = FALSE),
    "id_col"
  )
})

test_that("wrong blue_col name throws informative error", {
  # When the specified blue_col doesn't exist, it should not silently
  # fall through - it should error or warn clearly.
  expect_error(
    run_haplotype_prediction(ldx_geno, ldx_snp_info, ldx_blocks,
                             blues    = blues_st_df,
                             id_col   = "Genotype",
                             blue_col = "WRONG_TRAIT",
                             verbose  = FALSE)
  )
})
