## tests/testthat/test-association.R
## ─────────────────────────────────────────────────────────────────────────────
## Tests for haplotype_association.R and breeding_decisions.R:
##   test_block_haplotypes()
##   estimate_diplotype_effects()
##   score_favorable_haplotypes()
##   summarize_parent_haplotypes()
##
## All tests use the ldx_* example datasets plus small synthetic fixtures
## from helper.R. rrBLUP is required for the first two functions; tests
## are skipped automatically when the package is unavailable.
## ─────────────────────────────────────────────────────────────────────────────

library(testthat)
library(LDxBlocks)

data(ldx_geno,     package = "LDxBlocks")
data(ldx_snp_info, package = "LDxBlocks")
data(ldx_blocks,   package = "LDxBlocks")
data(ldx_blues,    package = "LDxBlocks")

# ── Shared fixtures ───────────────────────────────────────────────────────────

.haps <- extract_haplotypes(ldx_geno, ldx_snp_info, ldx_blocks, min_snps = 5L)
.blues_vec  <- setNames(ldx_blues$YLD, ldx_blues$id)
.blues_df1  <- ldx_blues[, c("id", "YLD")]           # single-trait data.frame
.blues_dfmt <- ldx_blues[, c("id", "YLD", "RES")]    # multi-trait data.frame

# Small synthetic fixture for fast tests
.G_s   <- make_geno(n = 50, p = 30, seed = 9L)
.si_s  <- make_snpinfo(p = 30)
.blk_s <- make_blocks(.si_s, n_blocks = 3L)
.haps_s <- extract_haplotypes(.G_s, .si_s, .blk_s, min_snps = 3L)
.blues_s <- make_blues(.G_s, seed = 9L)

# ── Deterministic diplotype fixture (exercises estimate_diplotype_effects) ───
# 18 individuals x 1 block x 5 SNPs.
# Exactly 6 AA, 6 AB, 6 BB diplotypes, ensuring dominance_table is non-empty.
# Blues are set so: a = 0.5, d = 1.0, d/a = 2.0 -> overdominance.
local({
  n_per   <- 6L
  n_snps  <- 5L
  n_ind   <- 3L * n_per

  # Dosage matrix: AA=all-0, AB=all-1, BB=all-2
  G_dip <- matrix(
    c(rep(0L, n_per * n_snps),   # AA individuals
      rep(1L, n_per * n_snps),   # AB individuals
      rep(2L, n_per * n_snps)),  # BB individuals
    nrow = n_ind, ncol = n_snps, byrow = FALSE
  )
  G_dip <- t(matrix(
    c(rep(0L, n_snps * n_per),
      rep(1L, n_snps * n_per),
      rep(2L, n_snps * n_per)),
    nrow = n_snps, ncol = n_ind
  ))
  rownames(G_dip) <- paste0("dip_ind", seq_len(n_ind))
  colnames(G_dip) <- paste0("dsnp", seq_len(n_snps))

  si_dip <- data.frame(
    SNP = paste0("dsnp", seq_len(n_snps)),
    CHR = "1",
    POS = seq(1000L, by = 2000L, length.out = n_snps),
    REF = "A", ALT = "T",
    stringsAsFactors = FALSE
  )
  blk_dip <- data.frame(
    start      = 1L, end = n_snps,
    start.rsID = "dsnp1", end.rsID = paste0("dsnp", n_snps),
    start.bp   = 1000L, end.bp = 1000L + (n_snps - 1L) * 2000L,
    CHR = "1",
    length_bp  = (n_snps - 1L) * 2000L + 1L,
    stringsAsFactors = FALSE
  )

  # Extract haplotypes (haplotype strings: "00000", "11111", "22222")
  .haps_dip <<- extract_haplotypes(G_dip, si_dip, blk_dip, min_snps = 3L)

  # Blues: AA mean=1, AB mean=3 (overdominance d/a=2>1), BB mean=2
  # a=(2-1)/2=0.5, d=3-1.5=1.5, d/a=3.0 (overdominance)
  set.seed(123L)
  blues_dip <- c(
    rnorm(n_per, mean = 1.0, sd = 0.2),   # AA
    rnorm(n_per, mean = 3.0, sd = 0.2),   # AB
    rnorm(n_per, mean = 2.0, sd = 0.2)    # BB
  )
  .blues_dip <<- setNames(blues_dip, rownames(G_dip))
})

# ══════════════════════════════════════════════════════════════════════════════
# 1. test_block_haplotypes
# ══════════════════════════════════════════════════════════════════════════════

test_that("test_block_haplotypes: returns LDxBlocks_haplotype_assoc", {
  skip_if_not_installed("rrBLUP")
  res <- test_block_haplotypes(.haps, blues = .blues_vec,
                               blocks = ldx_blocks, n_pcs = 0L, verbose = FALSE)
  expect_s3_class(res, "LDxBlocks_haplotype_assoc")
  expect_true(is.list(res))
})

test_that("test_block_haplotypes: allele_tests has required columns", {
  skip_if_not_installed("rrBLUP")
  res <- test_block_haplotypes(.haps, blues = .blues_vec,
                               blocks = ldx_blocks, n_pcs = 0L, verbose = FALSE)
  req <- c("block_id","CHR","start_bp","end_bp","trait",
           "allele","frequency","effect","SE","t_stat",
           "p_wald","p_wald_adj","significant")
  expect_true(all(req %in% names(res$allele_tests)),
              info = paste("Missing:", paste(setdiff(req, names(res$allele_tests)),
                                             collapse=", ")))
})

test_that("test_block_haplotypes: block_tests has required columns", {
  skip_if_not_installed("rrBLUP")
  res <- test_block_haplotypes(.haps, blues = .blues_vec,
                               blocks = ldx_blocks, n_pcs = 0L, verbose = FALSE)
  req <- c("block_id","CHR","start_bp","end_bp","trait",
           "n_alleles_tested","F_stat","df_LRT",
           "p_omnibus","p_omnibus_adj","var_explained","significant_omnibus")
  expect_true(all(req %in% names(res$block_tests)),
              info = paste("Missing:", paste(setdiff(req, names(res$block_tests)),
                                             collapse=", ")))
})

test_that("test_block_haplotypes: p_wald values in (0, 1]", {
  skip_if_not_installed("rrBLUP")
  res <- test_block_haplotypes(.haps, blues = .blues_vec,
                               blocks = ldx_blocks, n_pcs = 0L, verbose = FALSE)
  p <- res$allele_tests$p_wald
  expect_true(all(!is.na(p) & p > 0 & p <= 1))
})

test_that("test_block_haplotypes: p_omnibus values in (0, 1]", {
  skip_if_not_installed("rrBLUP")
  res <- test_block_haplotypes(.haps, blues = .blues_vec,
                               blocks = ldx_blocks, n_pcs = 0L, verbose = FALSE)
  p <- res$block_tests$p_omnibus
  expect_true(all(!is.na(p) & p > 0 & p <= 1))
})

test_that("test_block_haplotypes: var_explained in [0, 1]", {
  skip_if_not_installed("rrBLUP")
  res <- test_block_haplotypes(.haps, blues = .blues_vec,
                               blocks = ldx_blocks, n_pcs = 0L, verbose = FALSE)
  ve <- res$block_tests$var_explained
  expect_true(all(!is.na(ve) & ve >= 0 & ve <= 1))
})

test_that("test_block_haplotypes: frequency in (0, 1)", {
  skip_if_not_installed("rrBLUP")
  res <- test_block_haplotypes(.haps, blues = .blues_vec,
                               blocks = ldx_blocks, n_pcs = 0L, verbose = FALSE)
  f <- res$allele_tests$frequency
  expect_true(all(!is.na(f) & f > 0 & f <= 1))
})

test_that("test_block_haplotypes: significant is logical", {
  skip_if_not_installed("rrBLUP")
  res <- test_block_haplotypes(.haps, blues = .blues_vec,
                               blocks = ldx_blocks, n_pcs = 0L, verbose = FALSE)
  expect_type(res$allele_tests$significant,       "logical")
  expect_type(res$block_tests$significant_omnibus, "logical")
})

test_that("test_block_haplotypes: p_wald_adj >= p_wald (Bonferroni inflation)", {
  skip_if_not_installed("rrBLUP")
  res <- test_block_haplotypes(.haps, blues = .blues_vec,
                               blocks = ldx_blocks, n_pcs = 0L, verbose = FALSE)
  at <- res$allele_tests
  expect_true(all(at$p_wald_adj >= at$p_wald - 1e-10))
})

test_that("test_block_haplotypes: n_pcs_used = 0 when n_pcs = 0", {
  skip_if_not_installed("rrBLUP")
  res <- test_block_haplotypes(.haps, blues = .blues_vec,
                               blocks = ldx_blocks, n_pcs = 0L, verbose = FALSE)
  expect_equal(res$n_pcs_used, 0L)
})

test_that("test_block_haplotypes: n_pcs = 3 uses 3 GRM PCs", {
  skip_if_not_installed("rrBLUP")
  res <- test_block_haplotypes(.haps, blues = .blues_vec,
                               blocks = ldx_blocks, n_pcs = 3L, verbose = FALSE)
  expect_equal(res$n_pcs_used, 3L)
  # Results should still be valid
  expect_true(nrow(res$allele_tests) > 0)
  expect_true(all(res$allele_tests$p_wald > 0 & res$allele_tests$p_wald <= 1))
})

test_that("test_block_haplotypes: n_pcs = NULL auto-selects (0 to 10 PCs)", {
  skip_if_not_installed("rrBLUP")
  res <- test_block_haplotypes(.haps, blues = .blues_vec,
                               blocks = ldx_blocks, n_pcs = NULL, verbose = FALSE)
  expect_true(res$n_pcs_used >= 0L && res$n_pcs_used <= 10L)
})

test_that("test_block_haplotypes: single-trait data.frame blues works", {
  skip_if_not_installed("rrBLUP")
  haps <- .haps
  df   <- data.frame(id = ldx_blues$id, blue = ldx_blues$YLD)
  res  <- test_block_haplotypes(haps, blues = df, blocks = ldx_blocks,
                                n_pcs = 0L, id_col = "id",
                                blue_col = "blue", verbose = FALSE)
  expect_equal(length(res$traits), 1L)    # one trait returned
  expect_equal(res$traits, "blue")         # trait name = blue_col name
  expect_true(nrow(res$allele_tests) > 0)
})

test_that("test_block_haplotypes: multi-trait blues returns results for each trait", {
  skip_if_not_installed("rrBLUP")
  res <- test_block_haplotypes(.haps, blues = ldx_blues,
                               blocks = ldx_blocks,
                               id_col = "id", blue_cols = c("YLD","RES"),
                               n_pcs = 0L, verbose = FALSE)
  expect_true(all(c("YLD","RES") %in% res$allele_tests$trait))
  expect_true(all(c("YLD","RES") %in% res$block_tests$trait))
})

test_that("test_block_haplotypes: named list blues works", {
  skip_if_not_installed("rrBLUP")
  bl <- list(YLD = .blues_vec)
  res <- test_block_haplotypes(.haps, blues = bl,
                               blocks = ldx_blocks, n_pcs = 0L, verbose = FALSE)
  expect_equal(res$traits, "YLD")
  expect_true(nrow(res$allele_tests) > 0)
})

test_that("test_block_haplotypes: alpha stored in result", {
  skip_if_not_installed("rrBLUP")
  res <- test_block_haplotypes(.haps, blues = .blues_vec,
                               blocks = ldx_blocks, alpha = 0.01, verbose = FALSE)
  expect_equal(res$alpha, 0.01)
})

test_that("test_block_haplotypes: Bonferroni alpha = 0.05 / n_tests", {
  skip_if_not_installed("rrBLUP")
  res <- test_block_haplotypes(.haps, blues = .blues_vec,
                               blocks = ldx_blocks, verbose = FALSE)
  expect_equal(res$alpha, 0.05 / res$n_tests, tolerance = 1e-12)
})

test_that("test_block_haplotypes: allele_tests sorted by CHR then start_bp", {
  skip_if_not_installed("rrBLUP")
  res <- test_block_haplotypes(.haps, blues = .blues_vec,
                               blocks = ldx_blocks, n_pcs = 0L, verbose = FALSE)
  at  <- res$allele_tests[res$allele_tests$trait == res$traits[1], ]
  if (nrow(at) > 1)
    expect_true(all(order(at$CHR, at$start_bp, at$p_wald) ==
                      seq_len(nrow(at))))
})

test_that("test_block_haplotypes: print method works without error", {
  skip_if_not_installed("rrBLUP")
  res <- test_block_haplotypes(.haps, blues = .blues_vec,
                               blocks = ldx_blocks, verbose = FALSE)
  expect_output(print(res), "LDxBlocks Haplotype Association")
})

test_that("test_block_haplotypes: full ldx dataset, 9 blocks, n_pcs=0", {
  skip_if_not_installed("rrBLUP")
  res <- test_block_haplotypes(.haps, blues = .blues_vec,
                               blocks = ldx_blocks, n_pcs = 0L, verbose = FALSE)
  # All 9 blocks should have omnibus tests
  expect_equal(nrow(res$block_tests), 9L)
  # block_tests and allele_tests reference same block IDs
  expect_true(all(res$block_tests$block_id %in% res$allele_tests$block_id))
})

# ══════════════════════════════════════════════════════════════════════════════
# 2. estimate_diplotype_effects
# ══════════════════════════════════════════════════════════════════════════════

test_that("estimate_diplotype_effects: returns LDxBlocks_diplotype", {
  skip_if_not_installed("rrBLUP")
  res <- estimate_diplotype_effects(.haps, blues = .blues_vec,
                                    blocks = ldx_blocks, verbose = FALSE)
  expect_s3_class(res, "LDxBlocks_diplotype")
  expect_true(is.list(res))
})

test_that("estimate_diplotype_effects: three required list elements present", {
  skip_if_not_installed("rrBLUP")
  res <- estimate_diplotype_effects(.haps, blues = .blues_vec,
                                    blocks = ldx_blocks, verbose = FALSE)
  expect_true(all(c("diplotype_means","dominance_table","omnibus_tests")
                  %in% names(res)))
})

test_that("estimate_diplotype_effects: diplotype_means columns correct", {
  skip_if_not_installed("rrBLUP")
  res <- estimate_diplotype_effects(.haps, blues = .blues_vec,
                                    blocks = ldx_blocks, verbose = FALSE)
  if (nrow(res$diplotype_means) > 0) {
    req <- c("block_id","CHR","start_bp","end_bp","trait",
             "diplotype","n","mean_blue","se_mean")
    expect_true(all(req %in% names(res$diplotype_means)),
                info = paste("Missing:", paste(setdiff(req, names(res$diplotype_means)),
                                               collapse=", ")))
  }
})

test_that("estimate_diplotype_effects: dominance_table columns correct", {
  skip_if_not_installed("rrBLUP")
  res <- estimate_diplotype_effects(.haps_dip, blues = .blues_dip,
                                    blocks = NULL, verbose = FALSE)
  req <- c("block_id","trait","allele_A","allele_B",
           "mean_AA","mean_AB","mean_BB","a","d","d_over_a","overdominance")
  expect_gt(nrow(res$dominance_table), 0L)
  expect_true(all(req %in% names(res$dominance_table)),
              info = paste("Missing:", paste(setdiff(req, names(res$dominance_table)),
                                             collapse=", ")))
})

test_that("estimate_diplotype_effects: omnibus_tests columns correct", {
  skip_if_not_installed("rrBLUP")
  res <- estimate_diplotype_effects(.haps_dip, blues = .blues_dip,
                                    blocks = NULL, verbose = FALSE)
  req <- c("block_id","trait","n_diplotypes","F_stat",
           "df1","df2","p_omnibus","p_omnibus_adj","significant")
  expect_gt(nrow(res$omnibus_tests), 0L)
  expect_true(all(req %in% names(res$omnibus_tests)),
              info = paste("Missing:", paste(setdiff(req, names(res$omnibus_tests)),
                                             collapse=", ")))
})

test_that("estimate_diplotype_effects: overdominance is logical", {
  skip_if_not_installed("rrBLUP")
  res <- estimate_diplotype_effects(.haps_dip, blues = .blues_dip,
                                    blocks = NULL, verbose = FALSE)
  expect_gt(nrow(res$dominance_table), 0L)
  expect_type(res$dominance_table$overdominance, "logical")
})

test_that("estimate_diplotype_effects: significant is logical", {
  skip_if_not_installed("rrBLUP")
  res <- estimate_diplotype_effects(.haps_dip, blues = .blues_dip,
                                    blocks = NULL, verbose = FALSE)
  expect_gt(nrow(res$omnibus_tests), 0L)
  expect_type(res$omnibus_tests$significant, "logical")
})

test_that("estimate_diplotype_effects: p_omnibus in (0, 1]", {
  skip_if_not_installed("rrBLUP")
  res <- estimate_diplotype_effects(.haps_dip, blues = .blues_dip,
                                    blocks = NULL, verbose = FALSE)
  expect_gt(nrow(res$omnibus_tests), 0L)
  p <- res$omnibus_tests$p_omnibus
  expect_true(all(!is.na(p) & p > 0 & p <= 1))
})

test_that("estimate_diplotype_effects: diplotype strings are canonical (sorted)", {
  skip_if_not_installed("rrBLUP")
  res <- estimate_diplotype_effects(.haps, blues = .blues_vec,
                                    blocks = ldx_blocks, verbose = FALSE)
  if (nrow(res$diplotype_means) > 0) {
    # Each diplotype string "A/B" should satisfy A <= B (alphabetically)
    dips <- res$diplotype_means$diplotype
    halves <- strsplit(dips, "/", fixed = TRUE)
    canonical <- vapply(halves, function(h) {
      length(h) == 2L && h[1] <= h[2]
    }, logical(1L))
    expect_true(all(canonical),
                label = "All diplotype strings should be sorted (A <= B)")
  }
})

test_that("estimate_diplotype_effects: n in diplotype_means >= min_n_diplotype", {
  skip_if_not_installed("rrBLUP")
  min_n <- 4L
  res <- estimate_diplotype_effects(.haps, blues = .blues_vec,
                                    blocks = ldx_blocks,
                                    min_n_diplotype = min_n, verbose = FALSE)
  if (nrow(res$diplotype_means) > 0)
    expect_true(all(res$diplotype_means$n >= min_n))
})

test_that("estimate_diplotype_effects: a = (mean_BB - mean_AA) / 2", {
  skip_if_not_installed("rrBLUP")
  res <- estimate_diplotype_effects(.haps_dip, blues = .blues_dip,
                                    blocks = NULL, verbose = FALSE)
  dt <- res$dominance_table
  expect_gt(nrow(dt), 0L)
  expect_equal(round(dt$a, 4), round((dt$mean_BB - dt$mean_AA) / 2, 4))
})

test_that("estimate_diplotype_effects: d = mean_AB - midpoint(AA,BB)", {
  skip_if_not_installed("rrBLUP")
  res <- estimate_diplotype_effects(.haps_dip, blues = .blues_dip,
                                    blocks = NULL, verbose = FALSE)
  dt         <- res$dominance_table
  expect_gt(nrow(dt), 0L)
  expected_d <- round(dt$mean_AB - (dt$mean_AA + dt$mean_BB) / 2, 4)
  expect_equal(round(dt$d, 4), expected_d)
})

test_that("estimate_diplotype_effects: d_over_a = d/a when |a| > 0", {
  skip_if_not_installed("rrBLUP")
  res <- estimate_diplotype_effects(.haps_dip, blues = .blues_dip,
                                    blocks = NULL, verbose = FALSE)
  dt_v <- res$dominance_table[!is.na(res$dominance_table$d_over_a), ]
  expect_gt(nrow(dt_v), 0L)
  expect_equal(round(dt_v$d_over_a, 3), round(dt_v$d / dt_v$a, 3))
})

test_that("estimate_diplotype_effects: overdominance TRUE iff |d/a| > 1", {
  skip_if_not_installed("rrBLUP")
  res <- estimate_diplotype_effects(.haps_dip, blues = .blues_dip,
                                    blocks = NULL, verbose = FALSE)
  dt <- res$dominance_table[!is.na(res$dominance_table$d_over_a), ]
  expect_gt(nrow(dt), 0L)
  # blues designed so AA~1, AB~3, BB~2 -> d/a~3 -> overdominance
  expect_equal(dt$overdominance, abs(dt$d_over_a) > 1)
  expect_true(any(dt$overdominance))  # fixture guarantees at least one
})

test_that("estimate_diplotype_effects: multi-trait works", {
  skip_if_not_installed("rrBLUP")
  # Use deterministic fixture with two traits
  bl2 <- list(trait1 = .blues_dip,
              trait2 = setNames(rev(.blues_dip), names(.blues_dip)))
  res <- estimate_diplotype_effects(.haps_dip, blues = bl2,
                                    blocks = NULL, verbose = FALSE)
  expect_gt(nrow(res$omnibus_tests), 0L)
  expect_true("trait1" %in% res$omnibus_tests$trait)
  expect_true("trait2" %in% res$omnibus_tests$trait)
})

test_that("estimate_diplotype_effects: print method works", {
  skip_if_not_installed("rrBLUP")
  res <- estimate_diplotype_effects(.haps, blues = .blues_vec,
                                    blocks = ldx_blocks, verbose = FALSE)
  expect_output(print(res), "LDxBlocks Diplotype Effect Results")
})

# ══════════════════════════════════════════════════════════════════════════════
# 3. score_favorable_haplotypes
# ══════════════════════════════════════════════════════════════════════════════

# Build a minimal allele_effects table from the haplotype feature matrix
.ae <- do.call(rbind, lapply(names(.haps), function(bn) {
  h     <- .haps[[bn]]
  valid <- h[!grepl("\\.", h)]
  if (!length(valid)) return(NULL)
  freq_tbl <- sort(table(valid) / length(valid), decreasing = TRUE)
  # Keep alleles with frequency >= 0.05 (at least 6/120 individuals)
  common   <- names(freq_tbl)[freq_tbl >= 0.05]
  if (!length(common)) return(NULL)
  alleles  <- common[seq_len(min(3L, length(common)))]
  data.frame(block_id = bn, allele = alleles,
             allele_effect = rnorm(length(alleles), 0, 0.5),
             stringsAsFactors = FALSE)
}))

test_that("score_favorable_haplotypes: returns data.frame with required columns", {
  res <- score_favorable_haplotypes(.haps, allele_effects = .ae)
  expect_s3_class(res, "data.frame")
  req <- c("id","stacking_index","n_blocks_scored","mean_block_score","rank")
  expect_true(all(req %in% names(res)),
              info = paste("Missing:", paste(setdiff(req, names(res)), collapse=", ")))
})

test_that("score_favorable_haplotypes: one row per individual", {
  res <- score_favorable_haplotypes(.haps, allele_effects = .ae)
  expect_equal(nrow(res), length(.haps[[1]]))
})

test_that("score_favorable_haplotypes: stacking_index in [0, 1] when normalize=TRUE", {
  res <- score_favorable_haplotypes(.haps, allele_effects = .ae,
                                    normalize = TRUE)
  si <- res$stacking_index
  expect_true(all(!is.na(si) & si >= 0 & si <= 1))
})

test_that("score_favorable_haplotypes: normalize=FALSE returns raw scores (unbounded)", {
  res_n <- score_favorable_haplotypes(.haps, allele_effects = .ae, normalize = TRUE)
  res_r <- score_favorable_haplotypes(.haps, allele_effects = .ae, normalize = FALSE)
  # Raw scores need not be in [0,1]
  expect_false(all(res_r$stacking_index >= 0 & res_r$stacking_index <= 1))
  # But ranks should be identical (same ordering)
  expect_equal(res_n$rank, res_r$rank)
})

test_that("score_favorable_haplotypes: rank 1 has highest stacking_index", {
  res <- score_favorable_haplotypes(.haps, allele_effects = .ae)
  top <- res[res$rank == 1L, ]
  expect_true(all(res$stacking_index <= max(res$stacking_index) + 1e-8))
  expect_equal(top$stacking_index[1], max(res$stacking_index))
})

test_that("score_favorable_haplotypes: ranks are positive integers, max = n", {
  res <- score_favorable_haplotypes(.haps, allele_effects = .ae)
  # With ties.method = "first", all ranks are unique integers from 1 to nrow
  expect_equal(length(unique(res$rank)), nrow(res))
  expect_equal(min(res$rank), 1L)
  expect_equal(max(res$rank), nrow(res))
  expect_type(res$rank, "integer")
})

test_that("score_favorable_haplotypes: n_blocks_scored >= 0 and <= n_blocks", {
  res <- score_favorable_haplotypes(.haps, allele_effects = .ae)
  n_blk <- length(unique(.ae$block_id))
  expect_true(all(res$n_blocks_scored >= 0 & res$n_blocks_scored <= n_blk))
})

test_that("score_favorable_haplotypes: per-block score columns present", {
  res <- score_favorable_haplotypes(.haps, allele_effects = .ae)
  score_cols <- grep("^score_", names(res), value = TRUE)
  expect_true(length(score_cols) > 0)
  expect_equal(length(score_cols), length(unique(.ae$block_id)))
})

test_that("score_favorable_haplotypes: error when allele_effects missing required column", {
  ae_bad <- .ae[, c("block_id","allele")]   # missing allele_effect
  expect_error(score_favorable_haplotypes(.haps, allele_effects = ae_bad),
               "allele_effects missing")
})

test_that("score_favorable_haplotypes: no overlap between haplotypes and effects -> zero scores", {
  ae_no_match <- data.frame(block_id = "nonexistent_block",
                            allele = "111",
                            allele_effect = 1.0,
                            stringsAsFactors = FALSE)
  res <- suppressWarnings(
    score_favorable_haplotypes(.haps, allele_effects = ae_no_match)
  )
  expect_true(all(res$stacking_index == 0 | is.na(res$stacking_index)))
})

test_that("score_favorable_haplotypes: id column matches haplotype individual names", {
  res <- score_favorable_haplotypes(.haps, allele_effects = .ae)
  expect_true(all(res$id %in% names(.haps[[1]])))
})

test_that("score_favorable_haplotypes: min_freq=1 excludes all alleles -> zero scores", {
  res <- score_favorable_haplotypes(.haps, allele_effects = .ae, min_freq = 1.0)
  expect_true(all(res$stacking_index == 0))
})

# ══════════════════════════════════════════════════════════════════════════════
# 4. summarize_parent_haplotypes
# ══════════════════════════════════════════════════════════════════════════════

test_that("summarize_parent_haplotypes: returns data.frame", {
  res <- summarize_parent_haplotypes(.haps)
  expect_s3_class(res, "data.frame")
})

test_that("summarize_parent_haplotypes: has all required columns", {
  res <- summarize_parent_haplotypes(.haps)
  req <- c("id","block_id","CHR","start_bp","end_bp",
           "allele","dosage","allele_freq","allele_effect","is_rare")
  expect_true(all(req %in% names(res)),
              info = paste("Missing:", paste(setdiff(req, names(res)), collapse=", ")))
})

test_that("summarize_parent_haplotypes: dosage values are 0, 1, or 2", {
  res <- summarize_parent_haplotypes(.haps)
  expect_true(all(res$dosage %in% c(0L, 1L, 2L)))
})

test_that("summarize_parent_haplotypes: allele_freq in (0, 1]", {
  res <- summarize_parent_haplotypes(.haps)
  f <- res$allele_freq
  expect_true(all(!is.na(f) & f > 0 & f <= 1))
})

test_that("summarize_parent_haplotypes: is_rare is logical", {
  res <- summarize_parent_haplotypes(.haps)
  expect_type(res$is_rare, "logical")
})

test_that("summarize_parent_haplotypes: is_rare TRUE iff allele_freq < 0.10", {
  res <- summarize_parent_haplotypes(.haps)
  expect_equal(res$is_rare, res$allele_freq < 0.10)
})

test_that("summarize_parent_haplotypes: allele_effect is NA without allele_effects arg", {
  res <- summarize_parent_haplotypes(.haps)
  expect_true(all(is.na(res$allele_effect)))
})

test_that("summarize_parent_haplotypes: allele_effect populated when allele_effects supplied", {
  res <- summarize_parent_haplotypes(.haps, allele_effects = .ae)
  # Some rows should have non-NA effects (those matching block_id + allele)
  expect_true(any(!is.na(res$allele_effect)))
})

test_that("summarize_parent_haplotypes: candidate_ids filters to subset", {
  all_ids <- names(.haps[[1]])
  cands   <- all_ids[1:5]
  res_all <- summarize_parent_haplotypes(.haps)
  res_sub <- summarize_parent_haplotypes(.haps, candidate_ids = cands)
  expect_true(all(unique(res_sub$id) %in% cands))
  expect_true(nrow(res_sub) < nrow(res_all))
})

test_that("summarize_parent_haplotypes: candidate_ids not in haplotypes are silently dropped", {
  res <- summarize_parent_haplotypes(.haps, candidate_ids = c("fake_ind_99", "fake_ind_100"))
  # Should return empty data.frame (no valid IDs)
  expect_equal(nrow(res), 0L)
})

test_that("summarize_parent_haplotypes: min_freq=0.5 reduces rows vs default", {
  res_lo <- summarize_parent_haplotypes(.haps, min_freq = 0.02)
  res_hi <- summarize_parent_haplotypes(.haps, min_freq = 0.50)
  expect_true(nrow(res_hi) <= nrow(res_lo))
})

test_that("summarize_parent_haplotypes: sorted by id, CHR, start_bp, desc dosage", {
  res <- summarize_parent_haplotypes(.haps)
  if (nrow(res) > 1) {
    # Within each id x block, higher dosage rows should come first
    for (uid in unique(res$id)) {
      for (bn in unique(res$block_id[res$id == uid])) {
        sub <- res[res$id == uid & res$block_id == bn, ]
        if (nrow(sub) > 1)
          expect_true(all(diff(sub$dosage) <= 0),
                      label = paste("dosage not descending for", uid, bn))
      }
    }
  }
})

test_that("summarize_parent_haplotypes: works on full ldx dataset", {
  top3 <- names(.haps[[1]])[1:3]
  res  <- summarize_parent_haplotypes(.haps, candidate_ids = top3)
  expect_true(nrow(res) > 0)
  expect_equal(sort(unique(res$id)), sort(top3))
})

test_that("summarize_parent_haplotypes: allele_effects with wrong columns errors", {
  ae_bad <- data.frame(block_id = "b", allele = "a", stringsAsFactors = FALSE)
  expect_error(
    summarize_parent_haplotypes(.haps, allele_effects = ae_bad),
    "allele_effects missing"
  )
})

# ══════════════════════════════════════════════════════════════════════════════
# 6. Gap-coverage: edge cases and alternative input paths
# ══════════════════════════════════════════════════════════════════════════════

# ── Gap 1: estimate_diplotype_effects with PHASED input ──────────────────────
# Phased haplotypes use "|" separator strings; heterozygous is determined by
# h1 != h2 (not by dosage containing "1"), so the code path through
# infer_block_haplotypes differs from the unphased path.

test_that("estimate_diplotype_effects: phased input produces non-empty dominance_table", {
  skip_if_not_installed("rrBLUP")
  # Build a phased version of the diplotype fixture:
  # 6 AA (hap1=hap2="00000"), 6 AB (hap1="00000", hap2="11111"),
  # 6 BB (hap1=hap2="11111")
  n_per  <- 6L
  n_snps <- 5L
  n_ind  <- 3L * n_per
  # Dosage matrix: individuals x SNPs
  G_p <- rbind(
    matrix(0L, nrow = n_per, ncol = n_snps),  # AA
    matrix(1L, nrow = n_per, ncol = n_snps),  # AB
    matrix(2L, nrow = n_per, ncol = n_snps)   # BB
  )
  rownames(G_p) <- paste0("p_ind", seq_len(n_ind))
  colnames(G_p) <- paste0("ps", seq_len(n_snps))
  # Explicit phased gametes: individuals x SNPs
  h1_mat <- rbind(
    matrix(0L, nrow = n_per, ncol = n_snps),  # AA: hap1 = all-0
    matrix(0L, nrow = n_per, ncol = n_snps),  # AB: hap1 = all-0
    matrix(1L, nrow = n_per, ncol = n_snps)   # BB: hap1 = all-1
  )
  h2_mat <- rbind(
    matrix(0L, nrow = n_per, ncol = n_snps),  # AA: hap2 = all-0
    matrix(1L, nrow = n_per, ncol = n_snps),  # AB: hap2 = all-1
    matrix(1L, nrow = n_per, ncol = n_snps)   # BB: hap2 = all-1
  )
  dimnames(h1_mat) <- dimnames(h2_mat) <- dimnames(G_p)
  phased_in <- list(
    hap1       = t(h1_mat),   # SNPs x individuals
    hap2       = t(h2_mat),   # SNPs x individuals
    dosage     = t(G_p),      # SNPs x individuals
    sample_ids = rownames(G_p),
    phased     = TRUE
  )
  si_p <- data.frame(
    SNP = paste0("ps", seq_len(n_snps)),
    CHR = "1",
    POS = seq(1000L, by = 2000L, length.out = n_snps),
    REF = "A",
    ALT = "T",
    stringsAsFactors = FALSE
  )
  blk_p <- data.frame(
    start      = 1L,
    end        = n_snps,
    start.rsID = "ps1",
    end.rsID   = paste0("ps", n_snps),
    start.bp   = 1000L,
    end.bp     = 1000L + (n_snps - 1L) * 2000L,
    CHR        = "1",
    length_bp  = (n_snps - 1L) * 2000L + 1L,
    stringsAsFactors = FALSE
  )
  haps_phased <- extract_haplotypes(phased_in, si_p, blk_p, min_snps = 3L)
  # Blues must use same IDs as phased fixture ("p_ind*"), not .blues_dip ("dip_ind*")
  set.seed(123L)
  blues_p <- setNames(
    c(rnorm(n_per, 1.0, 0.2), rnorm(n_per, 3.0, 0.2), rnorm(n_per, 2.0, 0.2)),
    rownames(G_p)
  )
  res <- estimate_diplotype_effects(
    haps_phased,
    blues   = blues_p,
    blocks  = NULL,
    verbose = FALSE
  )
  expect_gt(nrow(res$dominance_table), 0L)
  expect_gt(nrow(res$omnibus_tests), 0L)
  dt <- res$dominance_table
  expect_equal(round(dt$a, 4), round((dt$mean_BB - dt$mean_AA) / 2, 4))
  expect_equal(
    round(dt$d, 4),
    round(dt$mean_AB - (dt$mean_AA + dt$mean_BB) / 2, 4)
  )
})

# ── Gap 2: All-homozygous population ─────────────────────────────────────────
# When all individuals are homozygous (dosage in {0,2} only), there are no AB
# diplotypes. estimate_diplotype_effects must return empty dominance_table and
# test_block_haplotypes must still produce allele_tests without error.

test_that("estimate_diplotype_effects: all-homozygous population -> empty dominance_table", {
  skip_if_not_installed("rrBLUP")
  # 12 individuals: 6 homozygous-ref (0), 6 homozygous-alt (2)
  n_hom  <- 6L
  n_snps <- 5L
  n_ind  <- 2L * n_hom
  G_hom <- matrix(
    c(rep(0L, n_snps * n_hom), rep(2L, n_snps * n_hom)),
    nrow = n_ind, ncol = n_snps, byrow = FALSE
  )
  rownames(G_hom) <- paste0("hom_ind", seq_len(n_ind))
  colnames(G_hom) <- paste0("hs", seq_len(n_snps))
  si_hom <- data.frame(SNP = paste0("hs", seq_len(n_snps)), CHR = "1",
                       POS = seq(1000L, by = 2000L, length.out = n_snps),
                       REF = "A", ALT = "T", stringsAsFactors = FALSE)
  blk_hom <- data.frame(start = 1L, end = n_snps, start.rsID = "hs1",
                        end.rsID = paste0("hs", n_snps),
                        start.bp = 1000L, end.bp = 1000L + (n_snps-1L)*2000L,
                        CHR = "1", length_bp = (n_snps-1L)*2000L + 1L,
                        stringsAsFactors = FALSE)
  haps_hom <- extract_haplotypes(G_hom, si_hom, blk_hom, min_snps = 3L)
  blues_hom <- setNames(c(rnorm(n_hom, 1, 0.2), rnorm(n_hom, 2, 0.2)),
                        rownames(G_hom))

  res <- estimate_diplotype_effects(haps_hom, blues = blues_hom,
                                    blocks = NULL, verbose = FALSE)
  # No heterozygous individuals -> no AB class -> empty dominance_table
  expect_equal(nrow(res$dominance_table), 0L)
  # diplotype_means should still have entries (AA and BB classes)
  expect_gt(nrow(res$diplotype_means), 0L)
})

test_that("test_block_haplotypes: all-homozygous population runs without error", {
  skip_if_not_installed("rrBLUP")
  n_hom  <- 10L; n_snps <- 5L; n_ind <- 2L * n_hom
  G_hom <- matrix(
    c(rep(0L, n_snps * n_hom), rep(2L, n_snps * n_hom)),
    nrow = n_ind, ncol = n_snps, byrow = FALSE
  )
  rownames(G_hom) <- paste0("hom_ind", seq_len(n_ind))
  colnames(G_hom) <- paste0("ts", seq_len(n_snps))
  si_hom <- data.frame(SNP = paste0("ts", seq_len(n_snps)), CHR = "1",
                       POS = seq(1000L, by = 2000L, length.out = n_snps),
                       REF = "A", ALT = "T", stringsAsFactors = FALSE)
  blk_hom <- data.frame(start = 1L, end = n_snps, start.rsID = "ts1",
                        end.rsID = paste0("ts", n_snps),
                        start.bp = 1000L, end.bp = 1000L + (n_snps-1L)*2000L,
                        CHR = "1", length_bp = (n_snps-1L)*2000L + 1L,
                        stringsAsFactors = FALSE)
  haps_hom <- extract_haplotypes(G_hom, si_hom, blk_hom, min_snps = 3L)
  blues_hom <- setNames(c(rnorm(n_hom, 0, 1), rnorm(n_hom, 1, 1)),
                        rownames(G_hom))
  res <- expect_no_error(
    test_block_haplotypes(haps_hom, blues = blues_hom,
                          blocks = NULL, n_pcs = 0L, verbose = FALSE)
  )
  expect_s3_class(res, "LDxBlocks_haplotype_assoc")
})

# ── Gap 3: Degenerate GRM (all identical genotypes) ──────────────────────────
# When all individuals carry identical genotypes, the GRM is rank-deficient
# (all entries equal). Both test_block_haplotypes and estimate_diplotype_effects
# must handle this gracefully via their internal tryCatch on mixed.solve().

test_that("test_block_haplotypes: degenerate GRM handled gracefully", {
  skip_if_not_installed("rrBLUP")
  n_ind  <- 12L; n_snps <- 5L
  # All individuals identical -> GRM = matrix of ones (singular)
  G_deg <- matrix(1L, nrow = n_ind, ncol = n_snps)
  rownames(G_deg) <- paste0("deg_ind", seq_len(n_ind))
  colnames(G_deg) <- paste0("ds", seq_len(n_snps))
  si_deg <- data.frame(SNP = paste0("ds", seq_len(n_snps)), CHR = "1",
                       POS = seq(1000L, by = 2000L, length.out = n_snps),
                       REF = "A", ALT = "T", stringsAsFactors = FALSE)
  blk_deg <- data.frame(start = 1L, end = n_snps, start.rsID = "ds1",
                        end.rsID = paste0("ds", n_snps),
                        start.bp = 1000L, end.bp = 1000L + (n_snps-1L)*2000L,
                        CHR = "1", length_bp = (n_snps-1L)*2000L + 1L,
                        stringsAsFactors = FALSE)
  haps_deg <- extract_haplotypes(G_deg, si_deg, blk_deg, min_snps = 3L)
  blues_deg <- setNames(rnorm(n_ind), rownames(G_deg))
  # Should not throw — degenerate GRM causes mixed.solve to fail,
  # which is caught internally and the block is skipped
  expect_no_error(
    test_block_haplotypes(haps_deg, blues = blues_deg,
                          blocks = NULL, n_pcs = 0L, verbose = FALSE)
  )
})

test_that("estimate_diplotype_effects: degenerate GRM handled gracefully", {
  skip_if_not_installed("rrBLUP")
  # Reuse the diplotype fixture but make all blues identical -> zero variance
  blues_zero <- setNames(rep(1.0, length(.blues_dip)), names(.blues_dip))
  expect_no_error(
    estimate_diplotype_effects(.haps_dip, blues = blues_zero,
                               blocks = NULL, verbose = FALSE)
  )
})

# ── Gap 4: d/a edge cases ─────────────────────────────────────────────────────
# a=0 (mean_AA == mean_BB): d_over_a must be NA, overdominance must be FALSE.
# d=0 (no dominance, AB = midpoint): overdominance must be FALSE.

test_that("estimate_diplotype_effects: a=0 gives NA d_over_a and overdominance=FALSE", {
  skip_if_not_installed("rrBLUP")
  n_per <- 6L
  # Deterministic fixture: mean_AA == mean_BB exactly, so a = 0 exactly
  blues_a0 <- setNames(
    c(
      rep(2.0, n_per),  # AA
      rep(3.5, n_per),  # AB
      rep(2.0, n_per)   # BB
    ),
    names(.blues_dip)
  )
  res <- suppressWarnings(
    estimate_diplotype_effects(
      .haps_dip,
      blues   = blues_a0,
      blocks  = NULL,
      verbose = FALSE
    )
  )
  expect_gt(nrow(res$dominance_table), 0L)
  dt <- res$dominance_table
  # Match the implementation contract exactly:
  # d_over_a is NA only when |a| <= 1e-10
  expect_true(all(abs(dt$a) < 1e-10))
  expect_true(all(is.na(dt$d_over_a)))
  expect_true(all(!dt$overdominance))
})

test_that("estimate_diplotype_effects: d=0 (no dominance) gives overdominance=FALSE", {
  skip_if_not_installed("rrBLUP")
  # Blues at midpoint for AB: d = AB - (AA+BB)/2 = 0 -> overdominance FALSE
  n_per <- 6L
  blues_d0 <- setNames(
    c(rnorm(n_per, mean = 1.0, sd = 0.1),   # AA
      rnorm(n_per, mean = 1.5, sd = 0.1),   # AB = midpoint -> d~0
      rnorm(n_per, mean = 2.0, sd = 0.1)),  # BB
    names(.blues_dip)
  )
  res <- estimate_diplotype_effects(.haps_dip, blues = blues_d0,
                                    blocks = NULL, verbose = FALSE)
  if (nrow(res$dominance_table) > 0) {
    expect_true(all(!res$dominance_table$overdominance))
  }
})
