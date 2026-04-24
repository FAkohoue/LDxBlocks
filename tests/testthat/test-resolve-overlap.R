## tests/testthat/test-resolve-overlap.R
## -----------------------------------------------------------------------------
## Unit tests for LD-informed overlap resolution:
##   - LDxBlocks:::.resolve_overlap()    R reference implementation
##   - LDxBlocks::resolve_overlap_cpp()  C++ accelerated implementation
##
## Parallel structure: every .resolve_overlap test has a matching
## resolve_overlap_cpp test, plus cross-validation tests that verify
## both implementations produce identical results on the same input.
##
## Biological guarantees tested:
##   1. Union merge when one core is empty (B inside A, or same start)
##   2. All-left path: disputed SNPs in stronger LD with block A
##   3. All-right path: disputed SNPs in stronger LD with block B
##   4. Mixed split: cumulative-score boundary rule
##   5. Tie-breaking: zero-variance overlap SNPs -> all-left (block A priority)
##   6. Non-overlapping input passes through unchanged
##   7. C++ and R produce identical boundaries on all cases above
## -----------------------------------------------------------------------------

library(testthat)
library(LDxBlocks)

data(ldx_geno,     package = "LDxBlocks")
data(ldx_snp_info, package = "LDxBlocks")

# -- Shared fixtures -----------------------------------------------------------

# Centred adjusted matrix for a 30-column window (used by most tests)
adj30 <- scale(ldx_geno[, 1:30], center = TRUE, scale = FALSE)

# Build a high-LD genotype block: n individuals sharing a haplotype
make_hap_block <- function(n, p, seed) {
  set.seed(seed)
  h <- sample(0:2, n, replace = TRUE)
  G <- matrix(h, n, p) + matrix(
    sample(-1:1, n * p, replace = TRUE, prob = c(.05, .9, .05)), n, p)
  pmax(pmin(G, 2L), 0L)
}


# -- 1. No overlap: input unchanged -------------------------------------------

test_that(".resolve_overlap: non-overlapping blocks pass through unchanged", {
  blocks_in <- matrix(c(1L, 10L, 12L, 20L), nrow = 2L, byrow = TRUE)
  result    <- LDxBlocks:::.resolve_overlap(blocks_in, adj30, k_rep = 5L)
  expect_equal(result, blocks_in)
})

test_that("resolve_overlap_cpp: non-overlapping blocks pass through unchanged", {
  blocks_in <- matrix(as.integer(c(1, 10, 12, 20)), nrow = 2L, byrow = TRUE)
  result    <- resolve_overlap_cpp(blocks_in, adj30, k_rep = 5L)
  expect_equal(result, blocks_in)
})

test_that("C++ vs R: non-overlapping input", {
  blocks_in <- matrix(as.integer(c(1, 8, 10, 18, 20, 28)), nrow = 3L, byrow = TRUE)
  r_res  <- LDxBlocks:::.resolve_overlap(blocks_in, adj30, k_rep = 5L)
  cpp_res <- resolve_overlap_cpp(blocks_in, adj30, k_rep = 5L)
  expect_equal(cpp_res, r_res)
})


# -- 2. Union merge: B fully inside A -----------------------------------------

test_that(".resolve_overlap: union merge when B fully inside A", {
  blocks_in <- matrix(c(1L, 30L, 10L, 20L), nrow = 2L, byrow = TRUE)
  result    <- LDxBlocks:::.resolve_overlap(blocks_in, adj30, k_rep = 5L)
  expect_equal(nrow(result), 1L)
  expect_equal(result[1L, 1L], 1L)
  expect_equal(result[1L, 2L], 30L)
})

test_that("resolve_overlap_cpp: union merge when B fully inside A", {
  blocks_in <- matrix(as.integer(c(1, 30, 10, 20)), nrow = 2L, byrow = TRUE)
  result    <- resolve_overlap_cpp(blocks_in, adj30, k_rep = 5L)
  expect_equal(nrow(result), 1L)
  expect_equal(result[1L, 1L], 1L)
  expect_equal(result[1L, 2L], 30L)
})

test_that("C++ vs R: union merge when B fully inside A", {
  blocks_in <- matrix(as.integer(c(1, 30, 10, 20)), nrow = 2L, byrow = TRUE)
  r_res   <- LDxBlocks:::.resolve_overlap(blocks_in, adj30, k_rep = 5L)
  cpp_res <- resolve_overlap_cpp(blocks_in, adj30, k_rep = 5L)
  expect_equal(nrow(cpp_res), nrow(r_res))
  expect_equal(cpp_res[1L, 1L], r_res[1L, 1L])
  expect_equal(cpp_res[1L, 2L], r_res[1L, 2L])
})


# -- 3. All-left path ----------------------------------------------------------

# One strong haplotype across all 25 cols -> overlap SNPs (16:20) all prefer left
.make_allleft_data <- function() {
  set.seed(5)
  G <- make_hap_block(80L, 25L, seed = 5L)
  list(adj = scale(G, center = TRUE, scale = FALSE),
       blocks = matrix(as.integer(c(1, 20, 16, 25)), nrow = 2L, byrow = TRUE))
}

test_that(".resolve_overlap: all-left produces two non-overlapping blocks", {
  d      <- .make_allleft_data()
  result <- LDxBlocks:::.resolve_overlap(d$blocks, d$adj, k_rep = 5L)
  expect_equal(nrow(result), 2L)
  expect_true(result[1L, 2L] < result[2L, 1L])
  expect_equal(result[1L, 1L], 1L)
  expect_equal(result[2L, 2L], 25L)
})

test_that("resolve_overlap_cpp: all-left produces two non-overlapping blocks", {
  d      <- .make_allleft_data()
  result <- resolve_overlap_cpp(d$blocks, d$adj, k_rep = 5L)
  expect_equal(nrow(result), 2L)
  expect_true(result[1L, 2L] < result[2L, 1L])
  expect_equal(result[1L, 1L], 1L)
  expect_equal(result[2L, 2L], 25L)
})

test_that("C++ vs R: all-left path produces identical boundaries", {
  d      <- .make_allleft_data()
  r_res  <- LDxBlocks:::.resolve_overlap(d$blocks, d$adj, k_rep = 5L)
  cpp_res <- resolve_overlap_cpp(d$blocks, d$adj, k_rep = 5L)
  expect_equal(nrow(cpp_res), nrow(r_res))
  expect_equal(cpp_res[1L, 2L], r_res[1L, 2L], label = "block A end matches")
  expect_equal(cpp_res[2L, 1L], r_res[2L, 1L], label = "block B start matches")
})


# -- 4. All-right path ---------------------------------------------------------

# Copy overlap columns to match right core exactly -> all scores negative
.make_allright_data <- function() {
  adj <- scale(ldx_geno[, 1:30], center = TRUE, scale = FALSE)
  adj[, 16:20] <- adj[, 21:25]  # overlap identical to right reps
  list(adj    = adj,
       blocks = matrix(as.integer(c(1, 20, 16, 30)), nrow = 2L, byrow = TRUE))
}

test_that(".resolve_overlap: all-right path -> block A ends at sB-1", {
  d      <- .make_allright_data()
  result <- LDxBlocks:::.resolve_overlap(d$blocks, d$adj, k_rep = 5L)
  expect_equal(nrow(result), 2L)
  expect_equal(result[1L, 2L], 15L)
  expect_equal(result[2L, 1L], 16L)
})

test_that("resolve_overlap_cpp: all-right path -> block A ends at sB-1", {
  d      <- .make_allright_data()
  result <- resolve_overlap_cpp(d$blocks, d$adj, k_rep = 5L)
  expect_equal(nrow(result), 2L)
  expect_equal(result[1L, 2L], 15L)
  expect_equal(result[2L, 1L], 16L)
})

test_that("C++ vs R: all-right path produces identical boundaries", {
  d       <- .make_allright_data()
  r_res   <- LDxBlocks:::.resolve_overlap(d$blocks, d$adj, k_rep = 5L)
  cpp_res <- resolve_overlap_cpp(d$blocks, d$adj, k_rep = 5L)
  expect_equal(cpp_res[1L, 2L], r_res[1L, 2L])
  expect_equal(cpp_res[2L, 1L], r_res[2L, 1L])
})


# -- 5. Tie-breaking: zero-variance overlap -> all-left ------------------------

.make_ties_data <- function() {
  adj <- scale(ldx_geno[, 1:30], center = TRUE, scale = FALSE)
  adj[, 16:20] <- 0  # constant -> r2 = 0 with everything -> score = 0 -> tie
  list(adj    = adj,
       blocks = matrix(as.integer(c(1, 20, 16, 30)), nrow = 2L, byrow = TRUE))
}

test_that(".resolve_overlap: ties -> all-left (block A priority)", {
  d      <- .make_ties_data()
  result <- LDxBlocks:::.resolve_overlap(d$blocks, d$adj, k_rep = 5L)
  expect_equal(nrow(result), 2L)
  expect_equal(result[1L, 2L], 20L)
  expect_equal(result[2L, 1L], 21L)
})

test_that("resolve_overlap_cpp: ties -> all-left (block A priority)", {
  d      <- .make_ties_data()
  result <- resolve_overlap_cpp(d$blocks, d$adj, k_rep = 5L)
  expect_equal(nrow(result), 2L)
  expect_equal(result[1L, 2L], 20L)
  expect_equal(result[2L, 1L], 21L)
})

test_that("C++ vs R: tie-breaking path produces identical boundaries", {
  d       <- .make_ties_data()
  r_res   <- LDxBlocks:::.resolve_overlap(d$blocks, d$adj, k_rep = 5L)
  cpp_res <- resolve_overlap_cpp(d$blocks, d$adj, k_rep = 5L)
  expect_equal(cpp_res[1L, 2L], r_res[1L, 2L])
  expect_equal(cpp_res[2L, 1L], r_res[2L, 1L])
})


# -- 6. Mixed split with cumulative-score rule ---------------------------------

.make_mixed_data <- function() {
  set.seed(42)
  n  <- 100L
  hA <- sample(0:2, n, replace = TRUE, prob = c(.2, .6, .2))
  hB <- sample(0:2, n, replace = TRUE, prob = c(.2, .6, .2))
  G_A  <- make_hap_block(n, 15L, seed = 42L)
  G_B  <- make_hap_block(n, 11L, seed = 43L)
  G_ov <- matrix(0L, n, 4L)
  for (j in 1:4) {
    w <- j / 5
    G_ov[, j] <- as.integer(round((1 - w) * hA + w * hB)) +
      sample(-1:1, n, replace = TRUE, prob = c(.05, .9, .05))
    G_ov[, j] <- pmax(pmin(G_ov[, j], 2L), 0L)
  }
  G <- cbind(G_A, G_ov, G_B)
  rownames(G) <- paste0("ind", seq_len(n))
  list(adj    = scale(G, center = TRUE, scale = FALSE),
       blocks = matrix(as.integer(c(1, 19, 16, 30)), nrow = 2L, byrow = TRUE))
}

test_that(".resolve_overlap: mixed split produces 2 non-overlapping blocks", {
  d      <- .make_mixed_data()
  result <- LDxBlocks:::.resolve_overlap(d$blocks, d$adj, k_rep = 10L)
  expect_equal(nrow(result), 2L)
  expect_true(result[1L, 2L] < result[2L, 1L])
  expect_equal(result[1L, 1L], 1L)
  expect_equal(result[2L, 2L], 30L)
})

test_that("resolve_overlap_cpp: mixed split produces 2 non-overlapping blocks", {
  d      <- .make_mixed_data()
  result <- resolve_overlap_cpp(d$blocks, d$adj, k_rep = 10L)
  expect_equal(nrow(result), 2L)
  expect_true(result[1L, 2L] < result[2L, 1L])
  expect_equal(result[1L, 1L], 1L)
  expect_equal(result[2L, 2L], 30L)
})

test_that("C++ vs R: mixed split produces identical boundaries", {
  d       <- .make_mixed_data()
  r_res   <- LDxBlocks:::.resolve_overlap(d$blocks, d$adj, k_rep = 10L)
  cpp_res <- resolve_overlap_cpp(d$blocks, d$adj, k_rep = 10L)
  expect_equal(nrow(cpp_res), nrow(r_res))
  expect_equal(cpp_res[1L, 2L], r_res[1L, 2L],
               label = "block A end: C++ matches R")
  expect_equal(cpp_res[2L, 1L], r_res[2L, 1L],
               label = "block B start: C++ matches R")
})


# -- 7. Output invariants ------------------------------------------------------

test_that("resolve_overlap_cpp: output has <= input rows (no new blocks created)", {
  blocks_in <- matrix(as.integer(c(1, 15, 12, 25, 23, 30)), nrow = 3L, byrow = TRUE)
  result    <- resolve_overlap_cpp(blocks_in, adj30, k_rep = 5L)
  expect_true(nrow(result) <= nrow(blocks_in))
})

test_that("resolve_overlap_cpp: all output blocks are non-overlapping", {
  blocks_in <- matrix(as.integer(c(1, 15, 12, 25, 23, 30)), nrow = 3L, byrow = TRUE)
  result    <- resolve_overlap_cpp(blocks_in, adj30, k_rep = 5L)
  if (nrow(result) > 1L) {
    expect_true(all(result[-1L, 1L] > result[-nrow(result), 2L]),
                label = "no overlaps remain in output")
  }
})

test_that("resolve_overlap_cpp: start <= end for all output blocks", {
  blocks_in <- matrix(as.integer(c(1, 15, 12, 25)), nrow = 2L, byrow = TRUE)
  result    <- resolve_overlap_cpp(blocks_in, adj30, k_rep = 5L)
  expect_true(all(result[, 1L] <= result[, 2L]))
})

test_that("resolve_overlap_cpp: global span preserved (min start, max end unchanged)", {
  blocks_in <- matrix(as.integer(c(1, 15, 12, 25)), nrow = 2L, byrow = TRUE)
  result    <- resolve_overlap_cpp(blocks_in, adj30, k_rep = 5L)
  expect_equal(min(result[, 1L]), 1L)
  expect_equal(max(result[, 2L]), 25L)
})


# -- 8. Big_LD integration: C++ resolver produces non-overlapping output -------

test_that("Big_LD with resolve_overlap_cpp: no overlapping blocks in output", {
  idx  <- which(ldx_snp_info$CHR == "1")
  blks <- LDxBlocks:::Big_LD(ldx_geno[, idx],
                             ldx_snp_info[idx, c("SNP", "POS")],
                             method = "r2", CLQcut = 0.5,
                             leng = 10L, subSegmSize = 70L,
                             verbose = FALSE)
  if (nrow(blks) > 1L) {
    expect_true(all(blks$start[-1L] > blks$end[-nrow(blks)]),
                label = "no index overlaps remain")
    expect_true(all(blks$start.bp[-1L] > blks$end.bp[-nrow(blks)]),
                label = "no bp overlaps remain")
  }
})

test_that("Big_LD with resolve_overlap_cpp: detects known block structure", {
  idx  <- which(ldx_snp_info$CHR == "1")
  blks <- LDxBlocks:::Big_LD(ldx_geno[, idx],
                             ldx_snp_info[idx, c("SNP", "POS")],
                             method = "r2", CLQcut = 0.5,
                             leng = 10L, subSegmSize = 70L,
                             verbose = FALSE)
  expect_true(nrow(blks) >= 3L)
})

# -- BUG 1 fix: non-adjacent contained blocks ---------------------------------
# Before the fix, blocks from different sub-segments arriving non-adjacent in
# the blocks matrix could produce overlapping output. The fix sorts by
# (start ASC, end DESC) and uses a max_end sweep to find ALL overlapping pairs.

test_that("resolve_overlap_cpp BUG1 fix: non-adjacent block B inside A is resolved", {
  # Block A = [1, 30], Block B = [8, 20] - B is fully inside A.
  # They are placed non-adjacently (rows 1 and 3) to trigger the old bug.
  # A third non-overlapping block C = [35, 45] separates them.
  data(ldx_geno,     package = "LDxBlocks")
  adj <- scale(ldx_geno[, 1:45], center = TRUE, scale = FALSE)
  # Row order: A [1,30], C [35,45], B [8,20] - B and A are NOT adjacent rows
  blocks_in <- matrix(
    as.integer(c(1, 30,   # row 1: block A (wide)
                 35, 45,  # row 2: block C (non-overlapping)
                 8, 20)), # row 3: block B (inside A, non-adjacent)
    nrow = 3L, byrow = TRUE
  )
  result <- resolve_overlap_cpp(blocks_in, adj, k_rep = 5L)
  # After fix: A and B must be merged (B inside A -> union [1,30])
  # or A extended; result must be non-overlapping
  if (nrow(result) > 1L) {
    expect_true(all(result[-1L, 1L] > result[-nrow(result), 2L]),
                label = "output must be non-overlapping after BUG1 fix")
  }
  # The block covering [1,30] should still exist
  covers_A <- any(result[, 1L] <= 1L & result[, 2L] >= 20L)
  expect_true(covers_A, label = "union of A and B should be retained")
})

test_that("resolve_overlap_cpp BUG2 fix: three consecutive overlapping blocks resolved", {
  # Three blocks with pairwise overlaps: A=[1,15], B=[12,25], C=[22,30].
  # A-B overlap [12,15] and B-C overlap [22,25].
  # BUG 2 (shed_row) caused index corruption when applying A-B merge then B-C merge.
  # Fix: keep[] vector + single-pass compaction.
  data(ldx_geno,     package = "LDxBlocks")
  adj <- scale(ldx_geno[, 1:30], center = TRUE, scale = FALSE)
  blocks_in <- matrix(
    as.integer(c(1,  15,   # A
                 12, 25,   # B (overlaps A)
                 22, 30)), # C (overlaps B)
    nrow = 3L, byrow = TRUE
  )
  result <- resolve_overlap_cpp(blocks_in, adj, k_rep = 5L)
  # Must produce non-overlapping output regardless of merge strategy
  expect_true(nrow(result) >= 1L, label = "at least one block must survive")
  if (nrow(result) > 1L) {
    expect_true(all(result[-1L, 1L] > result[-nrow(result), 2L]),
                label = "output must be non-overlapping after BUG2 fix")
  }
  # Global span must be preserved
  expect_equal(min(result[, 1L]), 1L,  label = "global start preserved")
  expect_equal(max(result[, 2L]), 30L, label = "global end preserved")
})
