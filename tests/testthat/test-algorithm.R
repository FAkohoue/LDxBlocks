## tests/testthat/test-algorithm.R
## ─────────────────────────────────────────────────────────────────────────────
## Tests for the Big-LD segmentation algorithm: LDxBlocks:::Big_LD(), LDxBlocks:::CLQD(),
## run_Big_LD_all_chr(), summarise_blocks(), plot_ld_blocks().
## Uses ldx_* example data and small synthetic fixtures.
## ─────────────────────────────────────────────────────────────────────────────

library(testthat)
library(LDxBlocks)

data(ldx_geno,     package = "LDxBlocks")
data(ldx_snp_info, package = "LDxBlocks")
data(ldx_blocks,   package = "LDxBlocks")

# ── CLQD ─────────────────────────────────────────────────────────────────────

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

# ── Big_LD ────────────────────────────────────────────────────────────────────

test_that("Big_LD: returns data.frame with required columns", {
  idx  <- which(ldx_snp_info$CHR == "1")
  blks <- LDxBlocks:::Big_LD(ldx_geno[, idx],
                             ldx_snp_info[idx, c("SNP","POS")],
                             method = "r2", CLQcut = 0.5,
                             leng = 10, subSegmSize = 70, verbose = FALSE)
  req <- c("start","end","start.rsID","end.rsID","start.bp","end.bp")
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

test_that("Big_LD: detectes known block structure in simulated data", {
  # Chr1 has 3 blocks. We should detect at least 3 blocks.
  idx  <- which(ldx_snp_info$CHR == "1")
  blks <- LDxBlocks:::Big_LD(ldx_geno[, idx],
                             ldx_snp_info[idx, c("SNP","POS")],
                             method = "r2", CLQcut = 0.5,
                             leng = 10, subSegmSize = 70, verbose = FALSE)
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

# ── run_Big_LD_all_chr ────────────────────────────────────────────────────────

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

  # With min_snps_chr = 10, chr 99 (5 SNPs) should be skipped
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

# ── summarise_blocks ──────────────────────────────────────────────────────────

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

# ── plot_ld_blocks ────────────────────────────────────────────────────────────

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
