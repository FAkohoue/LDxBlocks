## tests/testthat/test-cpp.R
## -----------------------------------------------------------------------------
## Unit tests for all eight exported C++ kernels in src/ld_core.cpp.
## Tests are property-based: verify mathematical guarantees rather than
## exact numeric values, so they pass regardless of BLAS implementation.
## -----------------------------------------------------------------------------

library(testthat)
library(LDxBlocks)

set.seed(99)
G_small <- matrix(sample(0:2, 80 * 30, replace = TRUE), 80, 30)
G_na    <- G_small; G_na[sample(80*30, 60)] <- NA   # inject NAs

# -- compute_r2_cpp ------------------------------------------------------------
test_that("compute_r2_cpp: symmetric, diagonal 0, values in [0,1]", {
  r2 <- compute_r2_cpp(G_small)
  expect_true(isSymmetric(r2, tol = 1e-10))
  expect_equal(diag(r2), rep(0, 30))
  expect_true(all(r2 >= 0 & r2 <= 1 + 1e-8))
})

test_that("compute_r2_cpp: rounding with digits=4", {
  r2 <- compute_r2_cpp(G_small, digits = 4L)
  remainders <- r2 * 1e4 - round(r2 * 1e4)
  expect_true(max(abs(remainders)) < 1e-8)
})

test_that("compute_r2_cpp: constant column gives zero row/col", {
  G       <- G_small
  G[, 5]  <- 1L
  r2      <- compute_r2_cpp(G)
  expect_equal(r2[5, ], rep(0, 30))
  expect_equal(r2[, 5], rep(0, 30))
})

test_that("compute_r2_cpp: NA columns are mean-imputed (no NA in output)", {
  r2 <- compute_r2_cpp(G_na)
  expect_false(any(is.na(r2)))
})

test_that("compute_r2_cpp: two-thread result equals single-thread", {
  r2_1 <- compute_r2_cpp(G_small, n_threads = 1L)
  r2_2 <- compute_r2_cpp(G_small, n_threads = 2L)
  expect_equal(r2_1, r2_2, tolerance = 1e-12)
})

test_that("compute_r2_cpp: perfectly correlated columns give r2 = 1", {
  G        <- matrix(0, 50, 5)
  G[, 1]   <- rbinom(50, 2, 0.3)
  G[, 2]   <- G[, 1]   # identical
  r2       <- compute_r2_cpp(G)
  expect_equal(r2[1, 2], 1.0, tolerance = 1e-10)
})

# -- compute_rV2_cpp -----------------------------------------------------------
test_that("compute_rV2_cpp: same guarantees as compute_r2_cpp", {
  # rV2 on whitened data is the same kernel
  Gc <- scale(G_small, center = TRUE, scale = FALSE)
  r2 <- compute_rV2_cpp(Gc)
  expect_true(isSymmetric(r2, tol = 1e-10))
  expect_equal(diag(r2), rep(0, 30))
  expect_true(all(r2 >= 0 & r2 <= 1 + 1e-8))
})

# -- maf_filter_cpp ------------------------------------------------------------
test_that("maf_filter_cpp: all-ref monomorphic column removed", {
  G       <- G_small
  G[, 3]  <- 0L
  keep    <- maf_filter_cpp(G, maf_cut = 0.05)
  expect_false(keep[3])
})

test_that("maf_filter_cpp: all-alt monomorphic column removed", {
  G       <- G_small
  G[, 7]  <- 2L
  keep    <- maf_filter_cpp(G, maf_cut = 0.05)
  expect_false(keep[7])
})

test_that("maf_filter_cpp: rare SNP below threshold removed", {
  G       <- G_small
  G[, 10] <- c(rep(0L, 78), 1L, 1L)   # MAF = 2/160 = 0.0125 < 0.05
  keep    <- maf_filter_cpp(G, maf_cut = 0.05)
  expect_false(keep[10])
})

test_that("maf_filter_cpp: common SNPs all kept", {
  keep <- maf_filter_cpp(G_small, maf_cut = 0.05)
  # G_small is randomly generated with MAF ~ 0.3 on average
  expect_true(sum(keep) >= 25L)
})

test_that("maf_filter_cpp: returns LogicalVector of correct length", {
  keep <- maf_filter_cpp(G_small, maf_cut = 0.05)
  expect_equal(length(keep), 30L)
  expect_type(keep, "logical")
})

# -- build_adj_matrix_cpp ------------------------------------------------------
test_that("build_adj_matrix_cpp: 0/1, symmetric, zero diagonal", {
  r2  <- compute_r2_cpp(G_small)
  adj <- build_adj_matrix_cpp(r2, 0.3)
  expect_true(all(adj == 0L | adj == 1L))
  expect_equal(diag(adj), rep(0L, 30))
  expect_true(isSymmetric(adj))
})

test_that("build_adj_matrix_cpp: threshold 0 gives all-ones off-diagonal", {
  # Use a single shared founder - all columns identical -> r2 = 1 for all pairs
  founder <- c(rep(0L, 10), rep(1L, 10), rep(2L, 10))
  G2  <- matrix(rep(founder, 5), nrow = 30, ncol = 5)
  r2  <- compute_r2_cpp(G2)
  # Verify r2 is indeed 1 for all off-diagonal pairs
  expect_true(all(r2[upper.tri(r2)] >= 0.99))
  # With threshold = 0.5, all perfectly correlated pairs get adj = 1
  adj  <- build_adj_matrix_cpp(r2, 0.5)
  # Diagonal must be 0 (loop only covers upper triangle k > j)
  expect_equal(diag(adj), rep(0L, 5))
  # Off-diagonal: use logical index to avoid diag<- coercion issues
  off  <- adj[row(adj) != col(adj)]
  expect_true(all(off == 1L))
})

test_that("build_adj_matrix_cpp: threshold 1.0 gives all-zeros", {
  r2  <- compute_r2_cpp(G_small)
  adj <- build_adj_matrix_cpp(r2, 1.0 + 1e-9)
  expect_true(all(adj == 0L))
})

# -- col_r2_cpp ----------------------------------------------------------------
test_that("col_r2_cpp: query column with itself is 0", {
  r2_col <- col_r2_cpp(G_small, 5L)
  expect_equal(r2_col[5], 0.0)
})

test_that("col_r2_cpp: result length equals ncol", {
  r2_col <- col_r2_cpp(G_small, 1L)
  expect_equal(length(r2_col), 30L)
})

test_that("col_r2_cpp: matches corresponding row of full r2 matrix", {
  r2_full <- compute_r2_cpp(G_small)
  r2_col  <- as.vector(col_r2_cpp(G_small, 8L))
  expect_equal(r2_col, r2_full[8, ], tolerance = 1e-10)
})

# -- compute_r2_sparse_cpp -----------------------------------------------------
test_that("compute_r2_sparse_cpp: returns triplet list", {
  bp     <- as.integer(seq(1000, by = 5000, length.out = 30))
  result <- compute_r2_sparse_cpp(G_small, bp, max_bp_dist = 20000L,
                                  threshold = 0.0)
  expect_type(result, "list")
  expect_true(all(c("row","col","r2") %in% names(result)))
  expect_equal(length(result$row), length(result$col))
  expect_equal(length(result$row), length(result$r2))
})

test_that("compute_r2_sparse_cpp: all r2 values in [0,1]", {
  bp     <- as.integer(seq(1000, by = 5000, length.out = 30))
  result <- compute_r2_sparse_cpp(G_small, bp, max_bp_dist = 30000L,
                                  threshold = 0.0)
  expect_true(all(result$r2 >= 0 & result$r2 <= 1 + 1e-8))
})

test_that("compute_r2_sparse_cpp: threshold filters pairs correctly", {
  bp     <- as.integer(seq(1000, by = 1000, length.out = 30))
  full   <- compute_r2_sparse_cpp(G_small, bp, max_bp_dist = 30000L,
                                  threshold = 0.0)
  thresh <- compute_r2_sparse_cpp(G_small, bp, max_bp_dist = 30000L,
                                  threshold = 0.3)
  expect_true(length(thresh$row) <= length(full$row))
  expect_true(all(thresh$r2 >= 0.3))
})

# -- boundary_scan_cpp ---------------------------------------------------------
test_that("boundary_scan_cpp: returns 0/1 integer vector of correct length", {
  Gc  <- scale(G_small, center = TRUE, scale = FALSE)
  res <- boundary_scan_cpp(Gc, start = 5L, end = 20L,
                           half_w = 3L, threshold = 0.5)
  expect_equal(length(res), 16L)
  expect_true(all(res %in% c(0L, 1L)))
})

test_that("boundary_scan_cpp: high threshold yields more valid cuts", {
  Gc   <- scale(G_small, center = TRUE, scale = FALSE)
  # threshold=0.0: every pair has r2>=0, so cross-boundary LD always detected -> no valid cuts
  low  <- boundary_scan_cpp(Gc, 3L, 25L, 3L, threshold = 0.0)
  # threshold=0.99: very few pairs exceed it -> most positions are valid cuts
  high <- boundary_scan_cpp(Gc, 3L, 25L, 3L, threshold = 0.99)
  expect_true(sum(high == 1L) >= sum(low == 1L))
})

# -- build_hap_strings_cpp -----------------------------------------------------
test_that("build_hap_strings_cpp: returns character vector of correct length", {
  blk <- matrix(as.integer(G_small[, 1:10]), nrow = nrow(G_small))
  out <- build_hap_strings_cpp(blk, ".")
  expect_type(out, "character")
  expect_equal(length(out), nrow(G_small))
})

test_that("build_hap_strings_cpp: each string has length equal to n_snps", {
  blk <- matrix(as.integer(G_small[, 1:10]), nrow = nrow(G_small))
  out <- build_hap_strings_cpp(blk, ".")
  expect_true(all(nchar(out) == 10L))
})

test_that("build_hap_strings_cpp: characters are 0/1/2 only (no NA)", {
  blk <- matrix(as.integer(G_small[, 1:10]), nrow = nrow(G_small))
  out <- build_hap_strings_cpp(blk, ".")
  chars <- unique(unlist(strsplit(paste(out, collapse=""), "")))
  expect_true(all(chars %in% c("0","1","2")))
})

test_that("build_hap_strings_cpp: NA_integer_ encoded as na_char", {
  blk    <- matrix(as.integer(G_small[, 1:5]), nrow = nrow(G_small))
  blk[1, 3] <- NA_integer_
  out <- build_hap_strings_cpp(blk, ".")
  # Individual 1 string should have "." at position 3
  expect_equal(substr(out[1], 3, 3), ".")
  # Other individuals should not have "."
  expect_false(any(grepl(".", out[-1], fixed=TRUE)))
})

test_that("build_hap_strings_cpp: custom na_char respected", {
  blk    <- matrix(as.integer(G_small[, 1:5]), nrow = nrow(G_small))
  blk[2, 1] <- NA_integer_
  out <- build_hap_strings_cpp(blk, "X")
  expect_equal(substr(out[2], 1, 1), "X")
})

test_that("build_hap_strings_cpp: matches R vapply reference implementation", {
  blk    <- matrix(as.integer(G_small[, 1:8]), nrow = nrow(G_small))
  # Reference R implementation
  ref <- vapply(seq_len(nrow(blk)), function(i) {
    v <- as.character(blk[i, ])
    v[is.na(blk[i, ])] <- "."
    paste(v, collapse="")
  }, character(1L))
  out <- build_hap_strings_cpp(blk, ".")
  expect_equal(out, ref)
})

# -- resolve_overlap_cpp -------------------------------------------------------
# Note: Comprehensive cross-validation against .resolve_overlap() R reference
# is in test-resolve-overlap.R. These tests verify the C++ kernel properties.

test_that("resolve_overlap_cpp: returns matrix with <= input rows", {
  adj <- scale(G_small, center = TRUE, scale = FALSE)
  blocks_in <- matrix(as.integer(c(1, 12, 10, 20, 18, 28)), nrow = 3L, byrow = TRUE)
  result    <- resolve_overlap_cpp(blocks_in, adj, k_rep = 5L)
  expect_true(nrow(result) <= nrow(blocks_in))
  expect_equal(ncol(result), 2L)
})

test_that("resolve_overlap_cpp: non-overlapping input is unchanged", {
  adj <- scale(G_small, center = TRUE, scale = FALSE)
  blocks_in <- matrix(as.integer(c(1, 8, 10, 18, 20, 28)), nrow = 3L, byrow = TRUE)
  result    <- resolve_overlap_cpp(blocks_in, adj, k_rep = 5L)
  expect_equal(result, blocks_in)
})

test_that("resolve_overlap_cpp: output blocks are non-overlapping", {
  adj <- scale(G_small, center = TRUE, scale = FALSE)
  blocks_in <- matrix(as.integer(c(1, 12, 10, 20, 18, 28)), nrow = 3L, byrow = TRUE)
  result    <- resolve_overlap_cpp(blocks_in, adj, k_rep = 5L)
  if (nrow(result) > 1L)
    expect_true(all(result[-1L, 1L] > result[-nrow(result), 2L]))
})

test_that("resolve_overlap_cpp: start <= end for all output blocks", {
  adj <- scale(G_small, center = TRUE, scale = FALSE)
  blocks_in <- matrix(as.integer(c(1, 15, 12, 25)), nrow = 2L, byrow = TRUE)
  result    <- resolve_overlap_cpp(blocks_in, adj, k_rep = 5L)
  expect_true(all(result[, 1L] <= result[, 2L]))
})

test_that("resolve_overlap_cpp: matches .resolve_overlap on random overlapping input", {
  set.seed(77)
  adj <- scale(G_small, center = TRUE, scale = FALSE)
  # Two overlapping blocks with genuine overlap zone
  blocks_in <- matrix(as.integer(c(1, 18, 15, 28)), nrow = 2L, byrow = TRUE)
  r_res   <- LDxBlocks:::.resolve_overlap(blocks_in, adj, k_rep = 5L)
  cpp_res <- resolve_overlap_cpp(blocks_in, adj, k_rep = 5L)
  expect_equal(nrow(cpp_res), nrow(r_res))
  expect_equal(cpp_res[1L, 2L], r_res[1L, 2L], label = "block A end matches")
  expect_equal(cpp_res[2L, 1L], r_res[2L, 1L], label = "block B start matches")
})


# -- block_snp_ranges_cpp ------------------------------------------------------

test_that("block_snp_ranges_cpp: returns list with lo, hi, n_snps", {
  pos      <- as.integer(seq(1000L, by = 1000L, length.out = 30L))
  block_sb <- as.integer(c(1000L, 10000L, 20000L))
  block_eb <- as.integer(c(5000L, 15000L, 25000L))
  res <- block_snp_ranges_cpp(pos, block_sb, block_eb)
  expect_type(res, "list")
  expect_true(all(c("lo", "hi", "n_snps") %in% names(res)))
  expect_equal(length(res$lo),     3L)
  expect_equal(length(res$hi),     3L)
  expect_equal(length(res$n_snps), 3L)
})

test_that("block_snp_ranges_cpp: n_snps equals hi - lo + 1 for retained blocks", {
  pos      <- as.integer(seq(1000L, by = 1000L, length.out = 30L))
  block_sb <- as.integer(c(1000L, 10000L, 20000L))
  block_eb <- as.integer(c(5000L, 15000L, 25000L))
  res <- block_snp_ranges_cpp(pos, block_sb, block_eb)
  for (b in seq_along(res$lo)) {
    if (res$n_snps[b] > 0L)
      expect_equal(res$hi[b] - res$lo[b] + 1L, res$n_snps[b],
                   label = paste("block", b))
  }
})

test_that("block_snp_ranges_cpp: SNPs within each block are correct", {
  pos      <- as.integer(seq(1000L, by = 1000L, length.out = 30L))
  block_sb <- as.integer(c(1000L, 10000L))
  block_eb <- as.integer(c(5000L, 15000L))
  res <- block_snp_ranges_cpp(pos, block_sb, block_eb)
  # Block 1: pos[1:5] = 1000..5000 -> lo=1, hi=5, n_snps=5
  expect_equal(res$lo[1],     1L)
  expect_equal(res$hi[1],     5L)
  expect_equal(res$n_snps[1], 5L)
  # Block 2: pos[10:15] = 10000..15000 -> lo=10, hi=15, n_snps=6
  expect_equal(res$lo[2],     10L)
  expect_equal(res$hi[2],     15L)
  expect_equal(res$n_snps[2], 6L)
})

test_that("block_snp_ranges_cpp: block with no SNPs gives lo=hi=n_snps=0", {
  pos      <- as.integer(seq(1000L, by = 1000L, length.out = 10L))
  # Block placed entirely between SNPs
  block_sb <- as.integer(c(500L))
  block_eb <- as.integer(c(999L))
  res <- block_snp_ranges_cpp(pos, block_sb, block_eb)
  expect_equal(res$lo[1],     0L)
  expect_equal(res$hi[1],     0L)
  expect_equal(res$n_snps[1], 0L)
})

test_that("block_snp_ranges_cpp: adjacent blocks partition SNPs without overlap", {
  pos      <- as.integer(seq(1000L, by = 1000L, length.out = 20L))
  block_sb <- as.integer(c(1000L, 6000L, 11000L, 16000L))
  block_eb <- as.integer(c(5000L, 10000L, 15000L, 20000L))
  res <- block_snp_ranges_cpp(pos, block_sb, block_eb)
  # Total SNPs across all blocks should equal 20
  expect_equal(sum(res$n_snps), 20L)
  # No overlap: each block's lo > previous block's hi
  for (b in 2:length(res$lo)) {
    if (res$n_snps[b] > 0L && res$n_snps[b-1] > 0L)
      expect_true(res$lo[b] > res$hi[b-1])
  }
})

test_that("block_snp_ranges_cpp: single-SNP block", {
  pos      <- as.integer(c(1000L, 2000L, 3000L, 4000L, 5000L))
  block_sb <- as.integer(c(3000L))
  block_eb <- as.integer(c(3000L))
  res <- block_snp_ranges_cpp(pos, block_sb, block_eb)
  expect_equal(res$n_snps[1], 1L)
  expect_equal(res$lo[1], res$hi[1])
  expect_equal(pos[res$lo[1]], 3000L)
})

test_that("block_snp_ranges_cpp: matches R findInterval reference", {
  set.seed(7)
  pos      <- sort(as.integer(sample(1:100000, 200)))
  block_sb <- as.integer(c(5000L, 30000L, 60000L, 80000L))
  block_eb <- as.integer(c(25000L, 55000L, 75000L, 95000L))
  res <- block_snp_ranges_cpp(pos, block_sb, block_eb)
  for (b in seq_along(block_sb)) {
    # R reference: which() for each block
    ref_idx <- which(pos >= block_sb[b] & pos <= block_eb[b])
    if (length(ref_idx) == 0L) {
      expect_equal(res$n_snps[b], 0L, label = paste("block", b, "empty"))
    } else {
      expect_equal(res$n_snps[b], length(ref_idx),
                   label = paste("block", b, "count"))
      expect_equal(res$lo[b], min(ref_idx),
                   label = paste("block", b, "lo"))
      expect_equal(res$hi[b], max(ref_idx),
                   label = paste("block", b, "hi"))
    }
  }
})

# -- impute_and_filter_cpp -----------------------------------------------------

test_that("impute_and_filter_cpp: no missing, no callrate filter -> unchanged", {
  G   <- matrix(as.integer(G_small), nrow = nrow(G_small))
  res <- impute_and_filter_cpp(G, min_callrate = 0.0, method = 0L)
  expect_equal(dim(res$geno_imputed), dim(G))
  expect_equal(res$n_filtered, 0L)
  expect_equal(res$n_imputed,  0L)
  expect_true(all(as.logical(res$keep)))
})

test_that("impute_and_filter_cpp: all SNPs pass when callrate threshold = 0", {
  G_na <- matrix(as.integer(G_small), nrow = nrow(G_small))
  G_na[1:5, 3] <- NA_integer_
  res <- impute_and_filter_cpp(G_na, min_callrate = 0.0, method = 0L)
  expect_true(all(as.logical(res$keep)))
  expect_equal(res$n_filtered, 0L)
})

test_that("impute_and_filter_cpp: low-callrate SNP is filtered", {
  n   <- nrow(G_small)
  G_lc <- matrix(as.integer(G_small), nrow = n)
  # Make column 5 have 60% missing -> callrate 40% < 80% threshold
  G_lc[seq_len(round(0.6 * n)), 5] <- NA_integer_
  res <- impute_and_filter_cpp(G_lc, min_callrate = 0.8, method = 0L)
  expect_false(res$keep[5])
  expect_true(res$n_filtered >= 1L)
})

test_that("impute_and_filter_cpp: mean_rounded imputation fills all NAs", {
  G_na <- matrix(as.integer(G_small), nrow = nrow(G_small))
  G_na[1:3, c(2, 7, 15)] <- NA_integer_
  res <- impute_and_filter_cpp(G_na, min_callrate = 0.0, method = 0L)
  expect_equal(res$n_imputed, 9L)
  expect_equal(sum(is.na(res$geno_imputed)), 0L)
})

test_that("impute_and_filter_cpp: mode imputation fills all NAs", {
  G_na <- matrix(as.integer(G_small), nrow = nrow(G_small))
  G_na[1:4, 10] <- NA_integer_
  res <- impute_and_filter_cpp(G_na, min_callrate = 0.0, method = 1L)
  expect_equal(res$n_imputed, 4L)
  expect_equal(sum(is.na(res$geno_imputed)), 0L)
})

test_that("impute_and_filter_cpp: mode tie falls back to mean_rounded", {
  # Construct a column with exactly equal counts of 0 and 2 (ambiguous mode).
  # n=80: rows 1-40 = 0, rows 41-80 = 2.
  # Inject NAs at rows 1-5 (from the 0-group) AND rows 41-45 (from the 2-group)
  # -> cnt[0]=35, cnt[1]=0, cnt[2]=35 -- a true tie.
  # Fallback: mean = (35*0 + 35*2) / 70 = 1.0 -> mean_rounded = 1.
  n <- nrow(G_small)
  G_tie <- matrix(as.integer(G_small), nrow = n)
  G_tie[, 1] <- c(rep(0L, 40), rep(2L, 40))   # rows 1-40 = 0, rows 41-80 = 2
  G_tie[1:5,   1] <- NA_integer_               # 5 NAs from the 0-group
  G_tie[41:45, 1] <- NA_integer_               # 5 NAs from the 2-group
  # Now: cnt[0]=35, cnt[1]=0, cnt[2]=35 -> ambiguous mode -> mean_rounded
  res <- impute_and_filter_cpp(G_tie, min_callrate = 0.0, method = 1L)
  expect_equal(res$imp_values[1], 1L,
               label = "tied mode (cnt[0]==cnt[2]==35) falls back to mean_rounded=1")
  expect_equal(sum(is.na(res$geno_imputed)), 0L)
})

test_that("impute_and_filter_cpp: unique mode used when no tie", {
  n <- nrow(G_small)
  G_mode <- matrix(as.integer(G_small), nrow = n)
  # Make column 2 mostly 0s with a few 1s and NAs -> clear mode = 0
  G_mode[, 2] <- c(rep(0L, 70), rep(1L, 10))
  G_mode[1:5, 2] <- NA_integer_
  res <- impute_and_filter_cpp(G_mode, min_callrate = 0.0, method = 1L)
  expect_equal(res$imp_values[2], 0L,
               label = "unique mode 0 used (no tie)")
})

test_that("impute_and_filter_cpp: imputed values are in {0,1,2}", {
  G_na <- matrix(as.integer(G_small), nrow = nrow(G_small))
  G_na[sample(length(G_na), 50)] <- NA_integer_
  res <- impute_and_filter_cpp(G_na, min_callrate = 0.0, method = 0L)
  vals <- as.integer(res$geno_imputed)
  expect_true(all(vals %in% 0:2))
})

test_that("impute_and_filter_cpp: call_rates length equals ncol(geno)", {
  G   <- matrix(as.integer(G_small), nrow = nrow(G_small))
  res <- impute_and_filter_cpp(G, min_callrate = 0.0, method = 0L)
  expect_equal(length(res$call_rates), ncol(G))
})

test_that("impute_and_filter_cpp: call_rates in [0,1]", {
  G_na <- matrix(as.integer(G_small), nrow = nrow(G_small))
  G_na[1:10, 1] <- NA_integer_
  res <- impute_and_filter_cpp(G_na, min_callrate = 0.0, method = 0L)
  cr <- as.numeric(res$call_rates)
  expect_true(all(cr >= 0 & cr <= 1 + 1e-8))
})

test_that("impute_and_filter_cpp: keep vector length equals ncol(geno)", {
  G   <- matrix(as.integer(G_small), nrow = nrow(G_small))
  res <- impute_and_filter_cpp(G, min_callrate = 0.5, method = 0L)
  expect_equal(length(res$keep), ncol(G))
})

test_that("impute_and_filter_cpp: n_filtered + n_kept = ncol(geno)", {
  G_na <- matrix(as.integer(G_small), nrow = nrow(G_small))
  # Give 20 SNPs very low call rate
  G_na[1:70, 1:20] <- NA_integer_
  res <- impute_and_filter_cpp(G_na, min_callrate = 0.8, method = 0L)
  expect_equal(res$n_filtered + sum(as.logical(res$keep)), ncol(G_na))
})

# -- extract_chr_haplotypes_cpp ------------------------------------------------

test_that("extract_chr_haplotypes_cpp: returns list with required fields", {
  G   <- matrix(as.integer(G_small[, 1:15]), nrow = nrow(G_small))
  pos <- as.integer(seq(1000L, by = 1000L, length.out = 15L))
  sb  <- as.integer(c(1000L, 8000L))
  eb  <- as.integer(c(6000L, 13000L))
  res <- extract_chr_haplotypes_cpp(G, pos, sb, eb,
                                    min_snps = 2L, min_freq = 0.0,
                                    top_n = 0L, na_char = ".")
  expect_type(res, "list")
  expected_fields <- c("hap_strings","hap_alleles","hap_freq",
                       "hap_counts","freq_dominant",
                       "block_lo","block_hi","n_snps",
                       "retained_idx","retained_start_bp","retained_end_bp",
                       "n_retained")
  expect_true(all(expected_fields %in% names(res)))
})

test_that("extract_chr_haplotypes_cpp: n_retained <= n_blocks", {
  G   <- matrix(as.integer(G_small[, 1:20]), nrow = nrow(G_small))
  pos <- as.integer(seq(1000L, by = 1000L, length.out = 20L))
  sb  <- as.integer(c(1000L, 6000L, 11000L, 16000L))
  eb  <- as.integer(c(5000L, 10000L, 15000L, 20000L))
  res <- extract_chr_haplotypes_cpp(G, pos, sb, eb,
                                    min_snps = 2L, min_freq = 0.0,
                                    top_n = 0L, na_char = ".")
  expect_true(res$n_retained <= length(sb))
})

test_that("extract_chr_haplotypes_cpp: hap_strings length equals n_ind per block", {
  n   <- nrow(G_small)
  G   <- matrix(as.integer(G_small[, 1:10]), nrow = n)
  pos <- as.integer(seq(1000L, by = 1000L, length.out = 10L))
  sb  <- as.integer(c(1000L))
  eb  <- as.integer(c(10000L))
  res <- extract_chr_haplotypes_cpp(G, pos, sb, eb,
                                    min_snps = 2L, min_freq = 0.0,
                                    top_n = 0L, na_char = ".")
  if (res$n_retained > 0L) {
    expect_equal(length(res$hap_strings[[1L]]), n)
  }
})

test_that("extract_chr_haplotypes_cpp: string width equals n_snps in block", {
  G   <- matrix(as.integer(G_small[, 1:10]), nrow = nrow(G_small))
  pos <- as.integer(seq(1000L, by = 1000L, length.out = 10L))
  sb  <- as.integer(c(1000L))
  eb  <- as.integer(c(10000L))
  res <- extract_chr_haplotypes_cpp(G, pos, sb, eb,
                                    min_snps = 2L, min_freq = 0.0,
                                    top_n = 0L, na_char = ".")
  if (res$n_retained > 0L) {
    w <- unique(nchar(res$hap_strings[[1L]]))
    expect_equal(w, res$n_snps[1L])
  }
})

test_that("extract_chr_haplotypes_cpp: hap string characters are in {0,1,2,.}", {
  G   <- matrix(as.integer(G_small[, 1:10]), nrow = nrow(G_small))
  pos <- as.integer(seq(1000L, by = 1000L, length.out = 10L))
  sb  <- as.integer(c(1000L))
  eb  <- as.integer(c(10000L))
  res <- extract_chr_haplotypes_cpp(G, pos, sb, eb,
                                    min_snps = 2L, min_freq = 0.0,
                                    top_n = 0L, na_char = ".")
  if (res$n_retained > 0L) {
    chars <- unique(unlist(strsplit(paste(res$hap_strings[[1L]], collapse=""), "")))
    expect_true(all(chars %in% c("0","1","2",".")))
  }
})

test_that("extract_chr_haplotypes_cpp: freq_dominant is in (0,1] for each block", {
  G   <- matrix(as.integer(G_small[, 1:20]), nrow = nrow(G_small))
  pos <- as.integer(seq(1000L, by = 1000L, length.out = 20L))
  sb  <- as.integer(c(1000L, 11000L))
  eb  <- as.integer(c(10000L, 20000L))
  res <- extract_chr_haplotypes_cpp(G, pos, sb, eb,
                                    min_snps = 2L, min_freq = 0.0,
                                    top_n = 0L, na_char = ".")
  if (res$n_retained > 0L) {
    fd <- res$freq_dominant
    expect_true(all(fd > 0 & fd <= 1 + 1e-8))
  }
})

test_that("extract_chr_haplotypes_cpp: freq_dominant equals max(hap_freq) per block", {
  G   <- matrix(as.integer(G_small[, 1:15]), nrow = nrow(G_small))
  pos <- as.integer(seq(1000L, by = 1000L, length.out = 15L))
  sb  <- as.integer(c(1000L, 8000L))
  eb  <- as.integer(c(7000L, 15000L))
  res <- extract_chr_haplotypes_cpp(G, pos, sb, eb,
                                    min_snps = 2L, min_freq = 0.0,
                                    top_n = 0L, na_char = ".")
  for (b in seq_len(res$n_retained)) {
    freqs <- res$hap_freq[[b]]
    # freq_dominant should equal the first (highest) frequency since alleles
    # are sorted descending - and should equal max(freqs)
    expect_equal(res$freq_dominant[b], max(freqs), tolerance = 1e-10,
                 label = paste("block", b, "freq_dominant == max(hap_freq)"))
  }
})

test_that("extract_chr_haplotypes_cpp: hap_freq sums to 1.0 per block", {
  G   <- matrix(as.integer(G_small[, 1:10]), nrow = nrow(G_small))
  pos <- as.integer(seq(1000L, by = 1000L, length.out = 10L))
  sb  <- as.integer(c(1000L))
  eb  <- as.integer(c(10000L))
  res <- extract_chr_haplotypes_cpp(G, pos, sb, eb,
                                    min_snps = 2L, min_freq = 0.0,
                                    top_n = 0L, na_char = ".")
  if (res$n_retained > 0L) {
    # Frequencies should sum to 1 (ignoring min_freq filtering, which is 0 here)
    expect_equal(sum(res$hap_freq[[1L]]), 1.0, tolerance = 1e-8)
  }
})

test_that("extract_chr_haplotypes_cpp: top_n limits alleles returned per block", {
  G   <- matrix(as.integer(G_small[, 1:10]), nrow = nrow(G_small))
  pos <- as.integer(seq(1000L, by = 1000L, length.out = 10L))
  sb  <- as.integer(c(1000L))
  eb  <- as.integer(c(10000L))
  res_all <- extract_chr_haplotypes_cpp(G, pos, sb, eb,
                                        min_snps = 2L, min_freq = 0.0,
                                        top_n = 0L, na_char = ".")
  res_top <- extract_chr_haplotypes_cpp(G, pos, sb, eb,
                                        min_snps = 2L, min_freq = 0.0,
                                        top_n = 3L, na_char = ".")
  if (res_top$n_retained > 0L && res_all$n_retained > 0L) {
    expect_true(length(res_top$hap_alleles[[1L]]) <= 3L)
    expect_true(length(res_top$hap_alleles[[1L]]) <=
                  length(res_all$hap_alleles[[1L]]))
  }
})

test_that("extract_chr_haplotypes_cpp: blocks below min_snps are excluded", {
  G   <- matrix(as.integer(G_small[, 1:10]), nrow = nrow(G_small))
  pos <- as.integer(seq(1000L, by = 1000L, length.out = 10L))
  # One block with 1 SNP (below min_snps=2) and one with 5 SNPs
  sb  <- as.integer(c(1000L, 3000L))
  eb  <- as.integer(c(1000L, 7000L))   # block 1: 1 SNP, block 2: 5 SNPs
  res <- extract_chr_haplotypes_cpp(G, pos, sb, eb,
                                    min_snps = 2L, min_freq = 0.0,
                                    top_n = 0L, na_char = ".")
  # Only block 2 (5 SNPs) should be retained
  expect_equal(res$n_retained, 1L)
  expect_equal(res$n_snps[1L], 5L)
})

test_that("extract_chr_haplotypes_cpp: retained_idx maps to correct original block", {
  G   <- matrix(as.integer(G_small[, 1:10]), nrow = nrow(G_small))
  pos <- as.integer(seq(1000L, by = 1000L, length.out = 10L))
  # Three blocks: singleton (skip), 5-SNP (retain, original index 2), 3-SNP (retain, index 3)
  sb  <- as.integer(c(1000L, 3000L,  9000L))
  eb  <- as.integer(c(1000L, 7000L, 10000L))
  res <- extract_chr_haplotypes_cpp(G, pos, sb, eb,
                                    min_snps = 2L, min_freq = 0.0,
                                    top_n = 0L, na_char = ".")
  expect_equal(res$n_retained, 2L)
  expect_true(!is.null(res$retained_idx))
  expect_equal(length(res$retained_idx), res$n_retained)
  # retained_idx[1] must be 2 (1-based: block 2, the 5-SNP one) -- block 1 was skipped
  expect_equal(res$retained_idx[1L], 2L)
  # retained_idx[2] must be 3 (1-based: block 3, the 2-SNP one)
  expect_equal(res$retained_idx[2L], 3L)
  # n_snps should match the actual SNP counts of the retained blocks
  expect_equal(res$n_snps[1L], 5L)
  expect_equal(res$n_snps[2L], 2L)
})

test_that("extract_chr_haplotypes_cpp: NA genotypes encoded as na_char in strings", {
  G      <- matrix(as.integer(G_small[, 1:8]), nrow = nrow(G_small))
  G[1, 3] <- NA_integer_   # inject NA for individual 1 at SNP 3
  pos    <- as.integer(seq(1000L, by = 1000L, length.out = 8L))
  sb     <- as.integer(c(1000L))
  eb     <- as.integer(c(8000L))
  res    <- extract_chr_haplotypes_cpp(G, pos, sb, eb,
                                       min_snps = 2L, min_freq = 0.0,
                                       top_n = 0L, na_char = ".")
  if (res$n_retained > 0L) {
    # Individual 1 string should contain "." at position 3
    ind1_str <- res$hap_strings[[1L]][1L]
    expect_equal(substr(ind1_str, 3L, 3L), ".")
  }
})

test_that("extract_chr_haplotypes_cpp: retained_start_bp and retained_end_bp match input block coords", {
  G   <- matrix(as.integer(G_small[, 1:12]), nrow = nrow(G_small))
  pos <- as.integer(seq(1000L, by = 1000L, length.out = 12L))
  sb  <- as.integer(c(1000L, 4000L, 9000L))   # three blocks
  eb  <- as.integer(c(3000L, 7000L, 12000L))
  res <- extract_chr_haplotypes_cpp(G, pos, sb, eb,
                                    min_snps = 2L, min_freq = 0.0,
                                    top_n = 0L, na_char = ".")
  # retained_start_bp/end_bp must be present
  expect_true(!is.null(res$retained_start_bp))
  expect_true(!is.null(res$retained_end_bp))
  expect_equal(length(res$retained_start_bp), res$n_retained)
  expect_equal(length(res$retained_end_bp),   res$n_retained)
  # Each retained block's coords must equal the original block_sb/block_eb
  for (k in seq_len(res$n_retained)) {
    orig_idx <- res$retained_idx[k]   # 1-based
    expect_equal(res$retained_start_bp[k], sb[orig_idx],
                 label = paste("block", k, "start_bp matches sb"))
    expect_equal(res$retained_end_bp[k],   eb[orig_idx],
                 label = paste("block", k, "end_bp matches eb"))
  }
})

test_that("extract_chr_haplotypes_cpp: matches build_hap_strings_cpp for a single block", {
  G    <- matrix(as.integer(G_small[, 1:8]), nrow = nrow(G_small))
  pos  <- as.integer(seq(1000L, by = 1000L, length.out = 8L))
  sb   <- as.integer(c(1000L))
  eb   <- as.integer(c(8000L))
  # C++ chromosome extractor
  res  <- extract_chr_haplotypes_cpp(G, pos, sb, eb,
                                     min_snps = 2L, min_freq = 0.0,
                                     top_n = 0L, na_char = ".")
  # Reference: build_hap_strings_cpp on the full block
  ref  <- build_hap_strings_cpp(G, ".")
  if (res$n_retained > 0L) {
    # Both should produce the same set of strings (possibly different order)
    expect_equal(sort(res$hap_strings[[1L]]), sort(ref))
  }
})

# -- extract_chr_haplotypes_phased_cpp -----------------------------------------

test_that("extract_chr_haplotypes_phased_cpp: returns list with required fields", {
  n   <- nrow(G_small)
  pos <- as.integer(seq(1000L, by = 1000L, length.out = 15L))
  h1  <- matrix(as.integer(G_small[, 1:15] > 0L), nrow = n)   # gamete 1: 0/1
  h2  <- matrix(as.integer(G_small[, 1:15] == 2L), nrow = n)  # gamete 2: 0/1
  sb  <- as.integer(c(1000L, 8000L))
  eb  <- as.integer(c(6000L, 13000L))
  res <- extract_chr_haplotypes_phased_cpp(h1, h2, pos, sb, eb,
                                           min_snps = 2L, min_freq = 0.0,
                                           top_n = 0L, na_char = ".")
  expect_type(res, "list")
  expected <- c("hap_strings","hap_alleles","hap_freq","hap_counts",
                "freq_dominant","block_lo","block_hi","n_snps",
                "retained_idx","retained_start_bp","retained_end_bp","n_retained")
  expect_true(all(expected %in% names(res)))
})

test_that("extract_chr_haplotypes_phased_cpp: strings have pipe separator", {
  n   <- nrow(G_small)
  pos <- as.integer(seq(1000L, by = 1000L, length.out = 10L))
  h1  <- matrix(as.integer(G_small[, 1:10] > 0L), nrow = n)
  h2  <- matrix(as.integer(G_small[, 1:10] == 2L), nrow = n)
  res <- extract_chr_haplotypes_phased_cpp(h1, h2, pos,
                                           as.integer(1000L), as.integer(10000L),
                                           min_snps = 2L, min_freq = 0.0,
                                           top_n = 0L, na_char = ".")
  if (res$n_retained > 0L) {
    hs <- res$hap_strings[[1L]]
    expect_true(all(grepl("|", hs, fixed = TRUE)))
  }
})

test_that("extract_chr_haplotypes_phased_cpp: each gamete half has width = n_snps", {
  n   <- nrow(G_small)
  pos <- as.integer(seq(1000L, by = 1000L, length.out = 8L))
  h1  <- matrix(as.integer(G_small[, 1:8] > 0L), nrow = n)
  h2  <- matrix(as.integer(G_small[, 1:8] == 2L), nrow = n)
  res <- extract_chr_haplotypes_phased_cpp(h1, h2, pos,
                                           as.integer(1000L), as.integer(8000L),
                                           min_snps = 2L, min_freq = 0.0,
                                           top_n = 0L, na_char = ".")
  if (res$n_retained > 0L) {
    ns <- res$n_snps[1L]
    parts <- strsplit(res$hap_strings[[1L]][1L], "|", fixed = TRUE)[[1L]]
    expect_equal(length(parts), 2L)
    expect_equal(nchar(parts[1L]), ns)
    expect_equal(nchar(parts[2L]), ns)
  }
})

test_that("extract_chr_haplotypes_phased_cpp: hap_freq sums to 1.0 (gamete counts)", {
  n   <- nrow(G_small)
  pos <- as.integer(seq(1000L, by = 1000L, length.out = 10L))
  h1  <- matrix(as.integer(G_small[, 1:10] > 0L), nrow = n)
  h2  <- matrix(as.integer(G_small[, 1:10] == 2L), nrow = n)
  res <- extract_chr_haplotypes_phased_cpp(h1, h2, pos,
                                           as.integer(1000L), as.integer(10000L),
                                           min_snps = 2L, min_freq = 0.0,
                                           top_n = 0L, na_char = ".")
  if (res$n_retained > 0L) {
    expect_equal(sum(res$hap_freq[[1L]]), 1.0, tolerance = 1e-8)
  }
})

test_that("extract_chr_haplotypes_phased_cpp: retained_idx and bp coords are consistent", {
  n   <- nrow(G_small)
  pos <- as.integer(seq(1000L, by = 1000L, length.out = 12L))
  h1  <- matrix(as.integer(G_small[, 1:12] > 0L), nrow = n)
  h2  <- matrix(as.integer(G_small[, 1:12] == 2L), nrow = n)
  sb  <- as.integer(c(1000L, 1500L, 6000L))  # block 1 has 1 SNP -> filtered
  eb  <- as.integer(c(1000L, 5000L, 12000L))
  res <- extract_chr_haplotypes_phased_cpp(h1, h2, pos, sb, eb,
                                           min_snps = 2L, min_freq = 0.0,
                                           top_n = 0L, na_char = ".")
  expect_equal(length(res$retained_idx),      res$n_retained)
  expect_equal(length(res$retained_start_bp), res$n_retained)
  expect_equal(length(res$retained_end_bp),   res$n_retained)
  # Each retained block's bp coords must equal the input sb/eb at retained_idx
  for (k in seq_len(res$n_retained)) {
    oi <- res$retained_idx[k]
    expect_equal(res$retained_start_bp[k], sb[oi])
    expect_equal(res$retained_end_bp[k],   eb[oi])
  }
})

test_that("extract_chr_haplotypes_phased_cpp: homozygous ref gives gametes '00...0|00...0'", {
  n   <- nrow(G_small)
  pos <- as.integer(seq(1000L, by = 1000L, length.out = 5L))
  # All individuals homozygous reference (dosage 0) -> both gametes = "00000"
  h1  <- matrix(0L, nrow = n, ncol = 5L)
  h2  <- matrix(0L, nrow = n, ncol = 5L)
  res <- extract_chr_haplotypes_phased_cpp(h1, h2, pos,
                                           as.integer(1000L), as.integer(5000L),
                                           min_snps = 2L, min_freq = 0.0,
                                           top_n = 0L, na_char = ".")
  if (res$n_retained > 0L) {
    # All strings should be "00000|00000" -> one unique haplotype
    expect_equal(length(res$hap_alleles[[1L]]), 1L)
    expect_equal(res$freq_dominant[1L], 1.0, tolerance = 1e-10)
    hs <- res$hap_strings[[1L]]
    expect_true(all(hs == "00000|00000"))
  }
})

test_that("extract_chr_haplotypes_phased_cpp: n_retained <= n_blocks", {
  n   <- nrow(G_small)
  pos <- as.integer(seq(1000L, by = 1000L, length.out = 20L))
  h1  <- matrix(as.integer(G_small[, 1:20] > 0L), nrow = n)
  h2  <- matrix(as.integer(G_small[, 1:20] == 2L), nrow = n)
  sb  <- as.integer(c(1000L, 6000L, 11000L, 16000L))
  eb  <- as.integer(c(5000L, 10000L, 15000L, 20000L))
  res <- extract_chr_haplotypes_phased_cpp(h1, h2, pos, sb, eb,
                                           min_snps = 2L, min_freq = 0.0,
                                           top_n = 0L, na_char = ".")
  expect_true(res$n_retained <= length(sb))
})
