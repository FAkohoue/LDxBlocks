## tests/testthat/test-basic.R
## ─────────────────────────────────────────────────────────────────────────────
## Basic smoke tests — verify that the package loads, example data is accessible,
## and high-level functions produce output of the correct type without crashing.
## These tests are intentionally fast (<10 s total) and use small subsets of
## ldx_geno so they pass on CRAN/CI with minimal resources.
## ─────────────────────────────────────────────────────────────────────────────

library(testthat)
library(LDxBlocks)

# ── Example data loads ────────────────────────────────────────────────────────
test_that("example datasets load and have correct dimensions", {
  data(ldx_geno,     package = "LDxBlocks")
  data(ldx_snp_info, package = "LDxBlocks")
  data(ldx_blocks,   package = "LDxBlocks")
  data(ldx_gwas,     package = "LDxBlocks")

  expect_true(is.matrix(ldx_geno))
  expect_equal(dim(ldx_geno), c(120L, 230L))
  expect_equal(nrow(ldx_snp_info), 230L)
  expect_true(all(c("SNP","CHR","POS","REF","ALT") %in% names(ldx_snp_info)))
  expect_equal(nrow(ldx_blocks), 9L)
  expect_true(all(c("start","end","start.rsID","end.rsID","start.bp","end.bp","CHR") %in%
                    names(ldx_blocks)))
  expect_equal(nrow(ldx_gwas), 20L)
  expect_true(all(c("Marker","CHR","POS") %in% names(ldx_gwas)))
})

test_that("example SNP IDs match between geno matrix and snp_info", {
  data(ldx_geno, package = "LDxBlocks")
  data(ldx_snp_info, package = "LDxBlocks")
  expect_equal(colnames(ldx_geno), ldx_snp_info$SNP)
})

test_that("example chromosomes are normalised (no chr prefix)", {
  data(ldx_snp_info, package = "LDxBlocks")
  expect_true(all(ldx_snp_info$CHR %in% c("1","2","3")))
})

# ── C++ kernels: basic contracts ──────────────────────────────────────────────
test_that("compute_r2_cpp returns symmetric matrix in [0,1] with zero diagonal", {
  data(ldx_geno, package = "LDxBlocks")
  r2 <- compute_r2_cpp(ldx_geno[, 1:20], digits = -1L, n_threads = 1L)
  expect_equal(dim(r2), c(20L, 20L))
  expect_true(isSymmetric(r2, tol = 1e-10))
  expect_true(all(r2 >= 0 & r2 <= 1 + 1e-8))
  expect_equal(diag(r2), rep(0, 20))
})

test_that("maf_filter_cpp removes low-MAF and monomorphic columns", {
  data(ldx_geno, package = "LDxBlocks")
  G        <- ldx_geno[, 1:30]
  G[, 1]   <- 0L      # monomorphic
  G[, 2]   <- c(rep(0L, 115), rep(1L, 5))  # MAF = 5/240 ~ 0.021 < 0.05
  keep     <- maf_filter_cpp(G, maf_cut = 0.05)
  expect_false(keep[1])
  expect_false(keep[2])
  expect_true(all(keep[3:30]))
})

test_that("build_adj_matrix_cpp produces valid 0/1 matrix", {
  data(ldx_geno, package = "LDxBlocks")
  r2  <- compute_r2_cpp(ldx_geno[, 1:15])
  adj <- build_adj_matrix_cpp(r2, 0.5)
  expect_true(all(adj == 0L | adj == 1L))
  expect_equal(diag(adj), rep(0L, 15))
  expect_true(isSymmetric(adj))
})

# ── read_geno: matrix backend ─────────────────────────────────────────────────
test_that("read_geno wraps matrix into valid LDxBlocks_backend", {
  data(ldx_geno, package = "LDxBlocks")
  data(ldx_snp_info, package = "LDxBlocks")
  be <- read_geno(ldx_geno, format = "matrix", snp_info = ldx_snp_info)
  expect_s3_class(be, "LDxBlocks_backend")
  expect_equal(be$type, "matrix")
  expect_equal(be$n_samples, 120L)
  expect_equal(be$n_snps, 230L)
  expect_equal(be$sample_ids, rownames(ldx_geno))
  close_backend(be)
})

test_that("read_chunk returns correct slice dimensions", {
  data(ldx_geno, package = "LDxBlocks")
  data(ldx_snp_info, package = "LDxBlocks")
  be    <- read_geno(ldx_geno, format = "matrix", snp_info = ldx_snp_info)
  chunk <- read_chunk(be, 10:25)
  expect_equal(dim(chunk), c(120L, 16L))
  close_backend(be)
})

test_that("close_backend is a no-op for matrix backend", {
  data(ldx_geno, package = "LDxBlocks")
  data(ldx_snp_info, package = "LDxBlocks")
  be <- read_geno(ldx_geno, format = "matrix", snp_info = ldx_snp_info)
  expect_null(close_backend(be))
})

# ── read_geno: flat-file formats ──────────────────────────────────────────────
test_that("read_geno reads numeric dosage CSV correctly", {
  f <- system.file("extdata", "example_genotypes_numeric.csv",
                   package = "LDxBlocks")
  skip_if(!file.exists(f), "extdata CSV not available")
  be <- read_geno(f)
  expect_s3_class(be, "LDxBlocks_backend")
  expect_equal(be$type, "numeric")
  expect_equal(be$n_snps, 230L)
  expect_equal(be$n_samples, 120L)
  chunk <- read_chunk(be, 1:10)
  expect_equal(dim(chunk), c(120L, 10L))
  expect_true(all(chunk %in% c(0, 1, 2, NA)))
  close_backend(be)
})

test_that("read_geno reads HapMap format correctly", {
  f <- system.file("extdata", "example_genotypes.hmp.txt",
                   package = "LDxBlocks")
  skip_if(!file.exists(f), "extdata HapMap not available")
  be <- read_geno(f)
  expect_s3_class(be, "LDxBlocks_backend")
  expect_equal(be$type, "hapmap")
  expect_equal(be$n_snps, 230L)
  close_backend(be)
})

test_that("read_geno reads VCF format correctly", {
  f <- system.file("extdata", "example_genotypes.vcf",
                   package = "LDxBlocks")
  skip_if(!file.exists(f), "extdata VCF not available")
  be <- read_geno(f)
  expect_s3_class(be, "LDxBlocks_backend")
  expect_equal(be$type, "vcf")
  expect_equal(be$n_snps, 230L)
  # Chromosome names should be normalised (no "chr" prefix)
  expect_true(all(be$snp_info$CHR %in% c("1","2","3")))
  chunk <- read_chunk(be, 1:5)
  expect_equal(dim(chunk), c(120L, 5L))
  close_backend(be)
})

# ── Big_LD: single chromosome ─────────────────────────────────────────────────
test_that("Big_LD (r2) returns valid block data.frame for chr1 subset", {
  data(ldx_geno, package = "LDxBlocks")
  data(ldx_snp_info, package = "LDxBlocks")
  chr1_idx <- which(ldx_snp_info$CHR == "1")
  g1       <- ldx_geno[, chr1_idx]
  s1       <- ldx_snp_info[chr1_idx, c("SNP","POS")]
  blocks   <- Big_LD(g1, s1, method = "r2", CLQcut = 0.5,
                     leng = 8, subSegmSize = 80, verbose = FALSE)
  expect_s3_class(blocks, "data.frame")
  expect_true(nrow(blocks) >= 1L)
  expect_true(all(c("start","end","start.rsID","end.rsID",
                    "start.bp","end.bp") %in% names(blocks)))
  expect_true(all(blocks$start <= blocks$end))
  expect_true(all(blocks$start.bp <= blocks$end.bp))
  # Blocks should not overlap
  if (nrow(blocks) > 1L) {
    expect_true(all(blocks$start.bp[-1] >= blocks$end.bp[-nrow(blocks)]))
  }
})

# ── run_Big_LD_all_chr ────────────────────────────────────────────────────────
test_that("run_Big_LD_all_chr returns blocks for all 3 chromosomes", {
  data(ldx_geno, package = "LDxBlocks")
  data(ldx_snp_info, package = "LDxBlocks")
  blocks <- run_Big_LD_all_chr(ldx_geno, snp_info = ldx_snp_info,
                               method = "r2", CLQcut = 0.4,
                               leng = 5, subSegmSize = 100, verbose = FALSE)
  expect_s3_class(blocks, "data.frame")
  expect_true(nrow(blocks) >= 3L)
  expect_true(all(c("CHR","start","end","start.bp","end.bp","length_bp") %in%
                    names(blocks)))
  expect_true(all(c("1","2","3") %in% blocks$CHR))
})

test_that("run_Big_LD_all_chr accepts LDxBlocks_backend", {
  data(ldx_geno, package = "LDxBlocks")
  data(ldx_snp_info, package = "LDxBlocks")
  be     <- read_geno(ldx_geno, format = "matrix", snp_info = ldx_snp_info)
  blocks <- run_Big_LD_all_chr(be, method = "r2", CLQcut = 0.5,
                               leng = 10, subSegmSize = 70, verbose = FALSE)
  expect_s3_class(blocks, "data.frame")
  expect_true(nrow(blocks) >= 3L)
  close_backend(be)
})

# ── summarise_blocks ──────────────────────────────────────────────────────────
test_that("summarise_blocks returns per-chromosome and genome-wide rows", {
  data(ldx_blocks, package = "LDxBlocks")
  s <- summarise_blocks(ldx_blocks)
  expect_true("GENOME" %in% s$CHR)
  expect_equal(nrow(s), 4L)   # chr1 + chr2 + chr3 + GENOME
  expect_true(all(s$n_blocks > 0))
  expect_true(all(s$median_bp > 0))
})

# ── compute_ld + prepare_geno ─────────────────────────────────────────────────
test_that("compute_r2 matches compute_r2_cpp for small window", {
  data(ldx_geno, package = "LDxBlocks")
  G  <- ldx_geno[, 1:15]
  r1 <- compute_r2(G, digits = 6L)
  r2 <- compute_r2_cpp(G, digits = 6L)
  expect_equal(r1, r2)
})

test_that("prepare_geno r2 returns centred matrix with NULL V_inv_sqrt", {
  data(ldx_geno, package = "LDxBlocks")
  prep <- prepare_geno(ldx_geno[, 1:20], method = "r2")
  expect_null(prep$V_inv_sqrt)
  expect_equal(dim(prep$adj_geno), c(120L, 20L))
  expect_true(max(abs(colMeans(prep$adj_geno))) < 1e-10)
})

# ── CLQD ─────────────────────────────────────────────────────────────────────
test_that("CLQD returns integer vector of correct length", {
  data(ldx_geno, package = "LDxBlocks")
  data(ldx_snp_info, package = "LDxBlocks")
  g    <- ldx_geno[, 1:25]
  info <- ldx_snp_info[1:25, c("SNP","POS")]
  Gc   <- scale(g, center = TRUE, scale = FALSE)
  bv   <- CLQD(g, info, Gc, CLQcut = 0.4, digits = -1L,
               n_threads = 1L, verbose = FALSE)
  expect_equal(length(bv), 25L)
  expect_type(bv, "integer")
})

# ── Haplotype pipeline ────────────────────────────────────────────────────────
test_that("extract_haplotypes returns list with one entry per block", {
  data(ldx_geno, package = "LDxBlocks")
  data(ldx_snp_info, package = "LDxBlocks")
  data(ldx_blocks, package = "LDxBlocks")
  haps <- extract_haplotypes(ldx_geno, ldx_snp_info, ldx_blocks, min_snps = 3)
  expect_type(haps, "list")
  expect_true(length(haps) >= 1L)
  expect_equal(length(haps[[1]]), 120L)
  expect_false(is.null(attr(haps, "block_info")))
})

test_that("compute_haplotype_diversity returns valid metrics", {
  data(ldx_geno, package = "LDxBlocks")
  data(ldx_snp_info, package = "LDxBlocks")
  data(ldx_blocks, package = "LDxBlocks")
  haps <- extract_haplotypes(ldx_geno, ldx_snp_info, ldx_blocks, min_snps = 3)
  div  <- compute_haplotype_diversity(haps)
  expect_s3_class(div, "data.frame")
  expect_true(all(c("block_id","n_ind","n_haplotypes","He","Shannon",
                    "freq_dominant") %in% names(div)))
  expect_true(all(div$He[!is.na(div$He)] >= 0 &
                    div$He[!is.na(div$He)] <= 1 + 1e-8))
  expect_true(all(div$Shannon[!is.na(div$Shannon)] >= 0))
  expect_true(all(div$freq_dominant[!is.na(div$freq_dominant)] > 0 &
                    div$freq_dominant[!is.na(div$freq_dominant)] <= 1))
})

test_that("build_haplotype_feature_matrix has correct dimensions", {
  data(ldx_geno, package = "LDxBlocks")
  data(ldx_snp_info, package = "LDxBlocks")
  data(ldx_blocks, package = "LDxBlocks")
  haps <- extract_haplotypes(ldx_geno, ldx_snp_info, ldx_blocks, min_snps = 3)
  feat <- build_haplotype_feature_matrix(haps, top_n = 3)
  expect_true(is.matrix(feat))
  expect_equal(nrow(feat), 120L)
  expect_equal(ncol(feat), length(haps) * 3L)
})

test_that("scaled haplotype feature matrix has near-zero column means", {
  data(ldx_geno, package = "LDxBlocks")
  data(ldx_snp_info, package = "LDxBlocks")
  data(ldx_blocks, package = "LDxBlocks")
  haps <- extract_haplotypes(ldx_geno, ldx_snp_info, ldx_blocks, min_snps = 3)
  feat <- build_haplotype_feature_matrix(haps, top_n = 3, scale_features = TRUE)
  col_means <- colMeans(feat, na.rm = TRUE)
  expect_true(all(abs(col_means) < 1e-10))
})

# ── tune_LD_params ────────────────────────────────────────────────────────────
test_that("tune_LD_params returns correct result structure", {
  data(ldx_geno, package = "LDxBlocks")
  data(ldx_snp_info, package = "LDxBlocks")
  data(ldx_gwas, package = "LDxBlocks")

  # Tiny grid to keep test fast
  grid <- expand.grid(
    CLQcut = c(0.5, 0.6), clstgap = 1e5, leng = 10,
    subSegmSize = 70, split = FALSE, checkLargest = FALSE,
    MAFcut = 0.05, CLQmode = "Density", kin_method = "chol",
    digits = -1, appendrare = FALSE, stringsAsFactors = FALSE
  )
  res <- tune_LD_params(ldx_geno, ldx_snp_info, ldx_gwas,
                        grid = grid, prefer_perfect = TRUE, seed = 1L)
  expect_type(res, "list")
  expect_true(all(c("best_params","score_table","final_blocks",
                    "gwas_assigned") %in% names(res)))
  expect_true(is.list(res$best_params))
  expect_s3_class(res$score_table, "data.frame")
  expect_true("LD_block" %in% names(res$gwas_assigned))
})
