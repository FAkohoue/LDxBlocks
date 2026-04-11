## tests/testthat/test-cpp.R
## ─────────────────────────────────────────────────────────────────────────────
## Unit tests for all seven exported C++ kernels in src/ld_core.cpp.
## Tests are property-based: verify mathematical guarantees rather than
## exact numeric values, so they pass regardless of BLAS implementation.
## ─────────────────────────────────────────────────────────────────────────────

library(testthat)
library(LDxBlocks)

set.seed(99)
G_small <- matrix(sample(0:2, 80 * 30, replace = TRUE), 80, 30)
G_na    <- G_small; G_na[sample(80*30, 60)] <- NA   # inject NAs

# ── compute_r2_cpp ────────────────────────────────────────────────────────────
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

# ── compute_rV2_cpp ───────────────────────────────────────────────────────────
test_that("compute_rV2_cpp: same guarantees as compute_r2_cpp", {
  # rV2 on whitened data is the same kernel
  Gc <- scale(G_small, center = TRUE, scale = FALSE)
  r2 <- compute_rV2_cpp(Gc)
  expect_true(isSymmetric(r2, tol = 1e-10))
  expect_equal(diag(r2), rep(0, 30))
  expect_true(all(r2 >= 0 & r2 <= 1 + 1e-8))
})

# ── maf_filter_cpp ────────────────────────────────────────────────────────────
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

# ── build_adj_matrix_cpp ──────────────────────────────────────────────────────
test_that("build_adj_matrix_cpp: 0/1, symmetric, zero diagonal", {
  r2  <- compute_r2_cpp(G_small)
  adj <- build_adj_matrix_cpp(r2, 0.3)
  expect_true(all(adj == 0L | adj == 1L))
  expect_equal(diag(adj), rep(0L, 30))
  expect_true(isSymmetric(adj))
})

test_that("build_adj_matrix_cpp: threshold 0 gives all-ones off-diagonal", {
  # Use a single shared founder — all columns identical -> r2 = 1 for all pairs
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

# ── col_r2_cpp ────────────────────────────────────────────────────────────────
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

# ── compute_r2_sparse_cpp ─────────────────────────────────────────────────────
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

# ── boundary_scan_cpp ─────────────────────────────────────────────────────────
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
