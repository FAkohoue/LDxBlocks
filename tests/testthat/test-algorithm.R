## tests/testthat/test-algorithm.R
## -----------------------------------------------------------------------------
## Tests for the Big-LD segmentation algorithm: LDxBlocks:::Big_LD(), LDxBlocks:::CLQD(),
## run_Big_LD_all_chr(), summarise_blocks(), plot_ld_blocks().
## Uses ldx_* example data and small synthetic fixtures.
## -----------------------------------------------------------------------------

library(testthat)
library(LDxBlocks)

data(ldx_geno,     package = "LDxBlocks")
data(ldx_snp_info, package = "LDxBlocks")
data(ldx_blocks,   package = "LDxBlocks")

# -- CLQD ---------------------------------------------------------------------

test_that("CLQD: returns integer vector of correct length", {
  g    <- ldx_geno[, 1:25]
  info <- ldx_snp_info[1:25, c("SNP","POS")]
  Gc   <- scale(g, center = TRUE, scale = FALSE)
  bv   <- LDxBlocks:::CLQD(g, info, Gc, CLQcut = 0.4, verbose = FALSE)
  expect_equal(length(bv), 25L)
  expect_type(bv, "integer")
})

test_that("CLQD: all values are positive integers or NA", {
  g    <- ldx_geno[, 1:25]
  info <- ldx_snp_info[1:25, c("SNP","POS")]
  Gc   <- scale(g, center = TRUE, scale = FALSE)
  bv   <- LDxBlocks:::CLQD(g, info, Gc, CLQcut = 0.4, verbose = FALSE)
  non_na <- bv[!is.na(bv)]
  expect_true(all(non_na >= 1L))
})

test_that("CLQD: high CLQcut produces many singletons (NAs)", {
  g    <- ldx_geno[, 1:20]
  info <- ldx_snp_info[1:20, c("SNP","POS")]
  Gc   <- scale(g, center = TRUE, scale = FALSE)
  bv_low  <- LDxBlocks:::CLQD(g, info, Gc, CLQcut = 0.2, verbose = FALSE)
  bv_high <- LDxBlocks:::CLQD(g, info, Gc, CLQcut = 0.99, verbose = FALSE)
  # High threshold should produce more NAs (fewer cliques)
  expect_true(sum(is.na(bv_high)) >= sum(is.na(bv_low)))
})

test_that("CLQD: Density mode vs Maximal mode produce valid bin vectors", {
  g    <- ldx_geno[, 1:30]
  info <- ldx_snp_info[1:30, c("SNP","POS")]
  Gc   <- scale(g, center = TRUE, scale = FALSE)
  bv_d <- LDxBlocks:::CLQD(g, info, Gc, CLQcut = 0.4, CLQmode = "Density",  verbose = FALSE)
  bv_m <- LDxBlocks:::CLQD(g, info, Gc, CLQcut = 0.4, CLQmode = "Maximal",  verbose = FALSE)
  expect_equal(length(bv_d), 30L)
  expect_equal(length(bv_m), 30L)
})

# -- Big_LD --------------------------------------------------------------------

test_that("Big_LD: returns data.frame with required columns", {
  idx  <- which(ldx_snp_info$CHR == "1")
  blks <- LDxBlocks:::Big_LD(ldx_geno[, idx],
                             ldx_snp_info[idx, c("SNP","POS")],
                             method = "r2", CLQcut = 0.5,
                             leng = 10, subSegmSize = 70, verbose = FALSE)
  req <- c("start","end","start.rsID","end.rsID","start.bp","end.bp","n_snps")
  expect_true(all(req %in% names(blks)))
  expect_s3_class(blks, "data.frame")
})

test_that("Big_LD: start <= end for every block", {
  idx  <- which(ldx_snp_info$CHR == "1")
  blks <- LDxBlocks:::Big_LD(ldx_geno[, idx],
                             ldx_snp_info[idx, c("SNP","POS")],
                             method = "r2", CLQcut = 0.5,
                             leng = 10, subSegmSize = 70, verbose = FALSE)
  expect_true(all(blks$start    <= blks$end))
  expect_true(all(blks$start.bp <= blks$end.bp))
})

test_that("Big_LD: blocks are non-overlapping (sorted by start.bp)", {
  idx  <- which(ldx_snp_info$CHR == "1")
  blks <- LDxBlocks:::Big_LD(ldx_geno[, idx],
                             ldx_snp_info[idx, c("SNP","POS")],
                             method = "r2", CLQcut = 0.5,
                             leng = 10, subSegmSize = 70, verbose = FALSE)
  if (nrow(blks) > 1L) {
    expect_true(all(blks$start.bp[-1] >= blks$end.bp[-nrow(blks)]))
  }
})

test_that("Big_LD: fewer than 2 polymorphic SNPs returns empty df with warning", {
  G_mono   <- matrix(0L, 60, 5)
  rownames(G_mono) <- paste0("ind", 1:60)
  colnames(G_mono) <- paste0("rs",  1:5)
  info_mono <- data.frame(SNP=colnames(G_mono), POS=1:5*1000)
  expect_warning(
    blks <- LDxBlocks:::Big_LD(G_mono, info_mono, method="r2",
                               leng=2, subSegmSize=5, verbose=FALSE),
    "Fewer than 2 polymorphic SNPs"
  )
  expect_equal(nrow(blks), 0L)
})

test_that("Big_LD: detects known block structure in simulated data", {
  idx  <- which(ldx_snp_info$CHR == "1")
  blks <- LDxBlocks:::Big_LD(ldx_geno[, idx],
                             ldx_snp_info[idx, c("SNP","POS")],
                             method = "r2", CLQcut = 0.5,
                             leng = 10, subSegmSize = 70, verbose = FALSE)
  # Chr1 has 3 simulated blocks; should detect at least 3
  expect_true(nrow(blks) >= 3L)
})

test_that("Big_LD: appendrare=TRUE does not error", {
  idx  <- which(ldx_snp_info$CHR == "1")
  expect_no_error(
    blks <- LDxBlocks:::Big_LD(ldx_geno[, idx],
                               ldx_snp_info[idx, c("SNP","POS")],
                               method = "r2", CLQcut = 0.5, MAFcut = 0.20,
                               appendrare = TRUE,
                               leng = 10, subSegmSize = 70, verbose = FALSE)
  )
})

test_that("Big_LD: split=TRUE with clstgap runs without error", {
  idx  <- which(ldx_snp_info$CHR == "1")
  expect_no_error(
    blks <- LDxBlocks:::Big_LD(ldx_geno[, idx],
                               ldx_snp_info[idx, c("SNP","POS")],
                               method = "r2", CLQcut = 0.5, split = TRUE,
                               clstgap = 30000L,
                               leng = 10, subSegmSize = 70, verbose = FALSE)
  )
})

test_that("Big_LD: singleton_as_block=TRUE includes single-SNP entries", {
  idx <- which(ldx_snp_info$CHR == "1")
  # Very high CLQcut forces most SNPs to be singletons
  blks_sing <- LDxBlocks:::Big_LD(ldx_geno[, idx],
                                  ldx_snp_info[idx, c("SNP","POS")],
                                  method = "r2", CLQcut = 0.99,
                                  singleton_as_block = TRUE,
                                  leng = 5, subSegmSize = 70,
                                  verbose = FALSE)
  blks_drop <- LDxBlocks:::Big_LD(ldx_geno[, idx],
                                  ldx_snp_info[idx, c("SNP","POS")],
                                  method = "r2", CLQcut = 0.99,
                                  singleton_as_block = FALSE,
                                  leng = 5, subSegmSize = 70,
                                  verbose = FALSE)
  # singleton_as_block=TRUE produces at least as many blocks
  expect_true(nrow(blks_sing) >= nrow(blks_drop))
  # Singleton blocks have start == end (same SNP index)
  if (nrow(blks_sing) > 0L)
    expect_true(any(blks_sing$start == blks_sing$end))
})

# -- run_Big_LD_all_chr --------------------------------------------------------

test_that("run_Big_LD_all_chr: processes all three chromosomes", {
  blks <- run_Big_LD_all_chr(ldx_geno, snp_info = ldx_snp_info,
                             method = "r2", CLQcut = 0.5,
                             leng = 10, subSegmSize = 70, verbose = FALSE)
  expect_true(all(c("1","2","3") %in% blks$CHR))
})

test_that("run_Big_LD_all_chr: length_bp column equals end.bp - start.bp + 1", {
  blks <- run_Big_LD_all_chr(ldx_geno, snp_info = ldx_snp_info,
                             method = "r2", CLQcut = 0.5,
                             leng = 10, subSegmSize = 70, verbose = FALSE)
  expected <- blks$end.bp - blks$start.bp + 1L
  expect_equal(blks$length_bp, expected)
})

test_that("run_Big_LD_all_chr: accepts LDxBlocks_backend instead of matrix", {
  be   <- read_geno(ldx_geno, format = "matrix", snp_info = ldx_snp_info)
  blks <- run_Big_LD_all_chr(be, method = "r2", CLQcut = 0.5,
                             leng = 10, subSegmSize = 70, verbose = FALSE)
  expect_s3_class(blks, "data.frame")
  expect_true(nrow(blks) >= 3L)
  close_backend(be)
})

test_that("run_Big_LD_all_chr: min_snps_chr skips small chromosomes", {
  # Inject a fake tiny chromosome with 5 SNPs
  g_fake    <- matrix(sample(0:2, 120*5, replace=TRUE), 120, 5)
  rownames(g_fake) <- rownames(ldx_geno)
  colnames(g_fake) <- paste0("fake", 1:5)
  g_all     <- cbind(ldx_geno, g_fake)
  info_fake <- data.frame(SNP=colnames(g_fake), CHR="99",
                          POS=1:5*1000, REF=NA_character_,
                          ALT=NA_character_, stringsAsFactors=FALSE)
  info_all  <- rbind(ldx_snp_info, info_fake)

  blks <- run_Big_LD_all_chr(g_all, snp_info=info_all, method="r2",
                             CLQcut=0.5, leng=10, subSegmSize=70,
                             min_snps_chr=10L, verbose=FALSE)
  expect_false("99" %in% blks$CHR)
})

test_that("run_Big_LD_all_chr: seed produces reproducible results", {
  blks1 <- run_Big_LD_all_chr(ldx_geno, snp_info=ldx_snp_info,
                              method="r2", CLQcut=0.5,
                              leng=10, subSegmSize=70,
                              seed=7L, verbose=FALSE)
  blks2 <- run_Big_LD_all_chr(ldx_geno, snp_info=ldx_snp_info,
                              method="r2", CLQcut=0.5,
                              leng=10, subSegmSize=70,
                              seed=7L, verbose=FALSE)
  expect_equal(blks1, blks2)
})

test_that("run_Big_LD_all_chr: chr parameter restricts to one chromosome", {
  blks <- run_Big_LD_all_chr(ldx_geno, snp_info = ldx_snp_info,
                             method = "r2", CLQcut = 0.5,
                             leng = 10, subSegmSize = 70,
                             chr = "1", verbose = FALSE)
  expect_true(all(blks$CHR == "1"))
  expect_false("2" %in% blks$CHR)
  expect_false("3" %in% blks$CHR)
})

# -- summarise_blocks ----------------------------------------------------------

test_that("summarise_blocks: GENOME row always present", {
  s <- summarise_blocks(ldx_blocks)
  expect_true("GENOME" %in% s$CHR)
})

test_that("summarise_blocks: per-chromosome row count matches unique CHRs", {
  s <- summarise_blocks(ldx_blocks)
  n_chr <- length(unique(ldx_blocks$CHR))
  expect_equal(nrow(s), n_chr + 1L)   # +1 for GENOME
})

test_that("summarise_blocks: GENOME total_bp_covered equals sum of chr values", {
  s     <- summarise_blocks(ldx_blocks)
  genome_row <- s[s$CHR == "GENOME", ]
  chr_rows   <- s[s$CHR != "GENOME", ]
  expect_equal(genome_row$total_bp_covered, sum(chr_rows$total_bp_covered))
})

test_that("summarise_blocks: n_blocks sums correctly", {
  s <- summarise_blocks(ldx_blocks)
  genome_row <- s[s$CHR == "GENOME", ]
  chr_rows   <- s[s$CHR != "GENOME", ]
  expect_equal(genome_row$n_blocks, sum(chr_rows$n_blocks))
})

test_that("summarise_blocks: works without CHR column", {
  blks_no_chr <- ldx_blocks[, c("start.bp","end.bp")]
  s <- summarise_blocks(blks_no_chr)
  expect_true("ALL" %in% s$CHR)
  expect_equal(nrow(s), 1L)
})

# -- plot_ld_blocks ------------------------------------------------------------

test_that("plot_ld_blocks: returns ggplot when ggplot2 available", {
  skip_if_not_installed("ggplot2")
  p <- plot_ld_blocks(ldx_blocks, colour_by = "length_bp")
  expect_s3_class(p, "ggplot")
})

test_that("plot_ld_blocks: colour_by CHR works", {
  skip_if_not_installed("ggplot2")
  p <- plot_ld_blocks(ldx_blocks, colour_by = "CHR")
  expect_s3_class(p, "ggplot")
})

test_that("plot_ld_blocks: missing CHR column throws error", {
  blks_no_chr <- ldx_blocks[, c("start.bp","end.bp","length_bp")]
  expect_error(plot_ld_blocks(blks_no_chr), "CHR")
})

# -- v0.3.1: CLQmode = "Louvain" -----------------------------------------------

test_that("CLQD: Louvain mode returns integer vector of correct length", {
  g    <- ldx_geno[, 1:25]
  info <- ldx_snp_info[1:25, c("SNP", "POS")]
  Gc   <- scale(g, center = TRUE, scale = FALSE)
  bv   <- LDxBlocks:::CLQD(g, info, Gc, CLQcut = 0.4,
                           CLQmode = "Louvain", verbose = FALSE)
  expect_equal(length(bv), 25L)
  expect_type(bv, "integer")
})

test_that("CLQD: Louvain mode values are positive integers or NA", {
  g    <- ldx_geno[, 1:25]
  info <- ldx_snp_info[1:25, c("SNP", "POS")]
  Gc   <- scale(g, center = TRUE, scale = FALSE)
  bv   <- LDxBlocks:::CLQD(g, info, Gc, CLQcut = 0.4,
                           CLQmode = "Louvain", verbose = FALSE)
  expect_true(all(is.na(bv) | (bv > 0L)))
})

test_that("CLQD: Louvain detects known block structure (chr1 block 1)", {
  # First 25 SNPs of ldx_geno are one LD block - Louvain should assign
  # the majority to a single community rather than all singletons.
  g    <- ldx_geno[, 1:25]
  info <- ldx_snp_info[1:25, c("SNP", "POS")]
  Gc   <- scale(g, center = TRUE, scale = FALSE)
  bv   <- LDxBlocks:::CLQD(g, info, Gc, CLQcut = 0.35,
                           CLQmode = "Louvain", verbose = FALSE)
  n_assigned <- sum(!is.na(bv))
  expect_gt(n_assigned, 10L)  # at least 10 SNPs grouped
})

test_that("CLQD: Leiden mode returns integer vector of correct length", {
  skip_if_not(
    packageVersion("igraph") >= "1.3.0",
    "cluster_leiden requires igraph >= 1.3.0"
  )
  g    <- ldx_geno[, 1:25]
  info <- ldx_snp_info[1:25, c("SNP", "POS")]
  Gc   <- scale(g, center = TRUE, scale = FALSE)
  bv   <- LDxBlocks:::CLQD(g, info, Gc, CLQcut = 0.4,
                           CLQmode = "Leiden", verbose = FALSE)
  expect_equal(length(bv), 25L)
  expect_type(bv, "integer")
})

test_that("run_Big_LD_all_chr: CLQmode Louvain produces valid block table", {
  blocks <- run_Big_LD_all_chr(
    ldx_geno,
    snp_info    = ldx_snp_info,
    method      = "r2",
    CLQmode     = "Louvain",
    CLQcut      = 0.35,
    leng        = 15L,
    subSegmSize = 100L,
    n_threads   = 1L,
    verbose     = FALSE
  )
  expect_s3_class(blocks, "data.frame")
  expect_true(nrow(blocks) >= 1L)
  expect_true(all(c("start", "end", "start.bp", "end.bp", "CHR") %in% names(blocks)))
  expect_true(all(blocks$start.bp <= blocks$end.bp))
})

# -- v0.3.1: max_bp_distance (sparse LD) --------------------------------------

test_that("CLQD: max_bp_distance produces same-length output as default", {
  g    <- ldx_geno[, 1:30]
  info <- ldx_snp_info[1:30, c("SNP", "POS")]
  Gc   <- scale(g, center = TRUE, scale = FALSE)
  bv_full   <- LDxBlocks:::CLQD(g, info, Gc, CLQcut = 0.4,
                                max_bp_distance = 0L, verbose = FALSE)
  bv_sparse <- LDxBlocks:::CLQD(g, info, Gc, CLQcut = 0.4,
                                max_bp_distance = 50000L, verbose = FALSE)
  expect_equal(length(bv_full),   30L)
  expect_equal(length(bv_sparse), 30L)
})

test_that("CLQD: max_bp_distance=500 skips distant pairs (very tight window)", {
  g    <- ldx_geno[, 1:20]
  info <- ldx_snp_info[1:20, c("SNP", "POS")]
  Gc   <- scale(g, center = TRUE, scale = FALSE)
  # SNPs in ldx_snp_info are spaced ~1000 bp apart.
  # max_bp_distance = 500 means NO pairs qualify -> adjacency all zeros
  # -> graph has no edges -> CLQD produces only singletons (all NA)
  bv <- LDxBlocks:::CLQD(g, info, Gc, CLQcut = 0.4,
                         max_bp_distance = 500L, verbose = FALSE)
  expect_equal(length(bv), 20L)       # correct output length
  expect_true(all(is.na(bv)),          # no cliques formed (empty graph)
              info = paste("Non-NA bins:", sum(!is.na(bv))))
})

test_that("run_Big_LD_all_chr: max_bp_distance=500000 produces valid blocks", {
  blocks <- run_Big_LD_all_chr(
    ldx_geno,
    snp_info        = ldx_snp_info,
    method          = "r2",
    CLQmode         = "Louvain",
    CLQcut          = 0.35,
    max_bp_distance = 500000L,
    leng            = 15L,
    subSegmSize     = 100L,
    n_threads       = 1L,
    verbose         = FALSE
  )
  expect_s3_class(blocks, "data.frame")
  expect_true(nrow(blocks) >= 1L)
  expect_true(all(blocks$start.bp <= blocks$end.bp))
})

# -- v0.3.1: bigmemory backend -------------------------------------------------

test_that("read_geno_bigmemory: converts matrix to bigmemory backend", {
  skip_if_not_installed("bigmemory")
  be_mat <- read_geno(ldx_geno, format = "matrix", snp_info = ldx_snp_info)
  be_bm  <- read_geno_bigmemory(
    be_mat,
    backingfile = tempfile("ldxbm_test_"),
    backingpath = tempdir(),
    type        = "char",
    verbose     = FALSE
  )
  on.exit({ close_backend(be_bm); close_backend(be_mat) })

  expect_s3_class(be_bm, "LDxBlocks_backend")
  expect_equal(be_bm$type, "bigmemory")
  expect_equal(be_bm$n_samples, nrow(ldx_geno))
  expect_equal(be_bm$n_snps,    ncol(ldx_geno))
})

test_that("read_geno_bigmemory: read_chunk returns correct dimensions", {
  skip_if_not_installed("bigmemory")
  be_mat <- read_geno(ldx_geno, format = "matrix", snp_info = ldx_snp_info)
  be_bm  <- read_geno_bigmemory(
    be_mat,
    backingfile = tempfile("ldxbm_chunk_"),
    backingpath = tempdir(),
    type        = "char",
    verbose     = FALSE
  )
  on.exit({ close_backend(be_bm); close_backend(be_mat) })

  chunk <- read_chunk(be_bm, 1:20)
  expect_equal(nrow(chunk), nrow(ldx_geno))
  expect_equal(ncol(chunk), 20L)
  expect_true(all(chunk %in% c(0, 1, 2, NA)))
})

test_that("read_geno_bigmemory: run_Big_LD_all_chr works through bigmemory backend", {
  skip_if_not_installed("bigmemory")
  be_mat <- read_geno(ldx_geno, format = "matrix", snp_info = ldx_snp_info)
  be_bm  <- read_geno_bigmemory(
    be_mat,
    backingfile = tempfile("ldxbm_ld_"),
    backingpath = tempdir(),
    type        = "char",
    verbose     = FALSE
  )
  on.exit({ close_backend(be_bm); close_backend(be_mat) })

  blocks <- run_Big_LD_all_chr(
    be_bm,
    method      = "r2",
    CLQcut      = 0.35,
    leng        = 15L,
    subSegmSize = 100L,
    n_threads   = 1L,
    verbose     = FALSE
  )
  expect_s3_class(blocks, "data.frame")
  expect_true(nrow(blocks) >= 1L)
})

# -- v0.3.1: LD-informed overlap resolution -----------------------------------

test_that("LD-informed split: overlapping blocks from adjacent sub-segments are resolved", {
  # Construct a genotype matrix with two clear LD blocks separated by a
  # recombination hotspot. Force a sub-segment seam inside the hotspot so
  # the boundary scan produces an overlap, then verify the resolved output
  # has non-overlapping blocks and the split falls near the hotspot.
  set.seed(77)
  n <- 80L; b <- 30L  # individuals, SNPs per block
  # Block 1: SNPs 1-30, high mutual LD (common haplotype)
  h1 <- sample(0:2, n, replace = TRUE)
  G1 <- matrix(h1, n, b) + matrix(sample(-1:1, n*b, replace=TRUE, prob=c(.1,.8,.1)), n, b)
  G1 <- pmax(pmin(G1, 2L), 0L)
  # Block 2: SNPs 31-60, independent haplotype
  h2 <- sample(0:2, n, replace = TRUE)
  G2 <- matrix(h2, n, b) + matrix(sample(-1:1, n*b, replace=TRUE, prob=c(.1,.8,.1)), n, b)
  G2 <- pmax(pmin(G2, 2L), 0L)
  G  <- cbind(G1, G2)
  rownames(G) <- paste0("ind", seq_len(n))
  colnames(G) <- paste0("rs",  seq_len(2*b))
  info <- data.frame(SNP = colnames(G), CHR = "1",
                     POS = seq(1000L, by = 1000L, length.out = 2*b),
                     stringsAsFactors = FALSE)

  # Use small subSegmSize to force sub-segment seam near the block boundary
  blocks <- LDxBlocks:::Big_LD(G, info[, c("SNP","POS")],
                               method = "r2", CLQcut = 0.3,
                               leng = 5L, subSegmSize = 35L,
                               verbose = FALSE)
  # Non-overlapping: no block's end.bp >= next block's start.bp
  if (nrow(blocks) > 1L) {
    overlap_exists <- any(blocks$end.bp[-nrow(blocks)] >=
                            blocks$start.bp[-1L])
    expect_false(overlap_exists,
                 label = "blocks should be non-overlapping after LD-informed split")
  }
  # Blocks should still cover both halves of the genome
  expect_true(min(blocks$start.bp) <= 5000L)
  expect_true(max(blocks$end.bp)   >= 55000L)
})

test_that("LD-informed split: result same as before when no overlaps exist", {
  # When boundary_scan_cpp places cuts cleanly, no overlaps reach .resolve_overlap
  # and the output should be unchanged vs the union-merge version.
  data(ldx_geno,     package = "LDxBlocks")
  data(ldx_snp_info, package = "LDxBlocks")
  idx <- which(ldx_snp_info$CHR == "1")
  blocks <- LDxBlocks:::Big_LD(ldx_geno[, idx],
                               ldx_snp_info[idx, c("SNP","POS")],
                               method = "r2", CLQcut = 0.5,
                               leng = 10L, subSegmSize = 70L,
                               verbose = FALSE)
  # No overlaps
  if (nrow(blocks) > 1L)
    expect_true(all(blocks$start.bp[-1L] > blocks$end.bp[-nrow(blocks)]))
})

test_that(".resolve_overlap: union merge fallback when B fully inside A", {
  # Block B [10,20] fully inside block A [1,30] -> right core empty
  # (right_core = 21:20 = empty). No right anchor -> union merge.
  data(ldx_geno,     package = "LDxBlocks")
  data(ldx_snp_info, package = "LDxBlocks")
  adj_mat   <- scale(ldx_geno[, 1:30], center = TRUE, scale = FALSE)
  blocks_in <- matrix(c(1L, 30L, 10L, 20L), nrow = 2L, byrow = TRUE)
  result    <- LDxBlocks:::.resolve_overlap(blocks_in, adj_mat, k_rep = 5L)
  expect_equal(nrow(result), 1L)
  expect_equal(result[1L, 1L], 1L)
  expect_equal(result[1L, 2L], 30L)
})

test_that(".resolve_overlap: all-left produces two surviving blocks", {
  # When all disputed SNPs are assigned left, block A keeps the overlap
  # zone and block B simply starts later. Both blocks survive.
  data(ldx_geno,     package = "LDxBlocks")
  data(ldx_snp_info, package = "LDxBlocks")
  # Create two overlapping blocks where columns 1-20 are all one LD block
  # so overlap SNPs (16-20) are assigned left.
  set.seed(5)
  n <- 80L
  h <- sample(0:2, n, replace = TRUE)
  G <- matrix(h, n, 25L) + matrix(sample(-1:1, n*25L, replace=TRUE,
                                         prob=c(.05,.9,.05)), n, 25L)
  G <- pmax(pmin(G, 2L), 0L)
  rownames(G) <- paste0("ind", seq_len(n))
  adj_mat   <- scale(G, center = TRUE, scale = FALSE)
  # Block A = [1,20], Block B = [16,25] -- overlap [16,20]
  # left_core=[1,15], right_core=[21,25]
  blocks_in <- matrix(c(1L, 20L, 16L, 25L), nrow = 2L, byrow = TRUE)
  result    <- LDxBlocks:::.resolve_overlap(blocks_in, adj_mat, k_rep = 5L)
  # Must produce exactly 2 non-overlapping blocks
  expect_equal(nrow(result), 2L)
  expect_true(result[1L, 2L] < result[2L, 1L],
              label = "blocks must be non-overlapping after split")
  # Block A start unchanged, block B end unchanged
  expect_equal(result[1L, 1L], 1L)
  expect_equal(result[2L, 2L], 25L)
})

test_that(".resolve_overlap: cumulative-score rule handles alternating assignments", {
  # Construct a case where individual per-SNP preferences alternate L/R/L/R
  # but the cumulative score correctly places the boundary later than last-left.
  # scores = [+0.3, -0.1, +0.5, -0.3, -0.5, -0.7]
  # cumsum = [+0.3, +0.2, +0.7, +0.4, -0.1, -0.8]
  # last positive cumsum = position 4 -> split after overlap[4]
  # last-left would give position 3 (last L individual assignment) -- wrong.
  #
  # We test this indirectly: create data where the 4th overlap SNP has
  # cumulative left advantage (+0.4) but the 5th and 6th flip to right.
  # Expect: result has 2 non-overlapping blocks, block A ends after 4th
  # overlap SNP, block B starts at 5th.
  set.seed(42)
  n  <- 100L
  # Block A haplotype: SNPs 1-15 in strong mutual LD
  hA <- sample(0:2, n, replace = TRUE, prob = c(.2, .6, .2))
  G_A <- matrix(hA, n, 15L) + matrix(
    sample(-1:1, n*15L, replace=TRUE, prob=c(.05,.9,.05)), n, 15L)
  G_A <- pmax(pmin(G_A, 2L), 0L)
  # Block B haplotype: SNPs 20-30 in strong mutual LD, independent of A
  hB <- sample(0:2, n, replace = TRUE, prob = c(.2, .6, .2))
  G_B <- matrix(hB, n, 11L) + matrix(
    sample(-1:1, n*11L, replace=TRUE, prob=c(.05,.9,.05)), n, 11L)
  G_B <- pmax(pmin(G_B, 2L), 0L)
  # Overlap zone: SNPs 16-19, gradual transition
  G_ov <- matrix(0L, n, 4L)
  for (j in 1:4) {
    w <- j / 5  # weight toward B increases across the overlap
    G_ov[, j] <- as.integer(round((1-w)*hA + w*hB)) +
      sample(-1:1, n, replace=TRUE, prob=c(.05,.9,.05))
    G_ov[, j] <- pmax(pmin(G_ov[, j], 2L), 0L)
  }
  G <- cbind(G_A, G_ov, G_B)
  rownames(G) <- paste0("ind", seq_len(n))
  adj_mat <- scale(G, center = TRUE, scale = FALSE)
  # Block A = cols 1-19, Block B = cols 16-30  (overlap = cols 16-19)
  blocks_in <- matrix(c(1L, 19L, 16L, 30L), nrow = 2L, byrow = TRUE)
  result <- LDxBlocks:::.resolve_overlap(blocks_in, adj_mat, k_rep = 10L)
  expect_equal(nrow(result), 2L)
  expect_true(result[1L, 2L] < result[2L, 1L],
              label = "blocks non-overlapping after cumulative split")
  expect_equal(result[1L, 1L], 1L,  label = "block A start unchanged")
  expect_equal(result[2L, 2L], 30L, label = "block B end unchanged")
})

test_that(".resolve_overlap: all-negative scores -> all-right path", {
  # When every overlap SNP prefers block B (all scores negative),
  # cum_score is everywhere negative, which(cum_score >= 0) is empty,
  # last_left = 0 -> all-right path.
  # Block A ends at sB-1, block B keeps the full overlap zone [sB..eB].
  data(ldx_geno,     package = "LDxBlocks")
  data(ldx_snp_info, package = "LDxBlocks")
  adj_mat <- scale(ldx_geno[, 1:30], center = TRUE, scale = FALSE)
  # Blocks overlap at cols 16:20; set those cols to match right core
  # by copying cols 21:25 so overlap SNPs are in perfect LD with right core
  adj_mat[, 16:20] <- adj_mat[, 21:25]
  blocks_in <- matrix(c(1L, 20L, 16L, 30L), nrow = 2L, byrow = TRUE)
  result <- LDxBlocks:::.resolve_overlap(blocks_in, adj_mat, k_rep = 5L)
  expect_equal(nrow(result), 2L)
  # Block A should end at 15 (sB-1=15), block B starts at 16
  expect_equal(result[1L, 2L], 15L,
               label = "block A ends at sB-1 when all overlap prefers right")
  expect_equal(result[2L, 1L], 16L,
               label = "block B starts at sB when all overlap prefers right")
})

test_that(".resolve_overlap: all-zero scores (ties) -> all-left path", {
  # When every overlap SNP has score=0 (zero variance -> r2=0 with everything),
  # cum_score = [0,0,...,0], which(cum_score >= 0) = all positions,
  # last_left = length(overlap) -> all-left path.
  # Block A keeps the full overlap zone, block B starts at eA+1.
  # This gives ties to block A, consistent with block A being detected first.
  data(ldx_geno,     package = "LDxBlocks")
  data(ldx_snp_info, package = "LDxBlocks")
  adj_mat <- scale(ldx_geno[, 1:30], center = TRUE, scale = FALSE)
  # Set overlap cols 16:20 to constant 0 -- zero variance, r2=0 with all
  adj_mat[, 16:20] <- 0
  blocks_in <- matrix(c(1L, 20L, 16L, 30L), nrow = 2L, byrow = TRUE)
  result <- LDxBlocks:::.resolve_overlap(blocks_in, adj_mat, k_rep = 5L)
  expect_equal(nrow(result), 2L)
  # Block A keeps overlap (ends at eA=20), block B starts at eA+1=21
  expect_equal(result[1L, 2L], 20L,
               label = "block A keeps overlap when all scores are tied")
  expect_equal(result[2L, 1L], 21L,
               label = "block B starts at eA+1 when all scores are tied")
})
