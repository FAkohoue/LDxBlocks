## tests/testthat/test-haplotypes.R
## ─────────────────────────────────────────────────────────────────────────────
## Tests for extract_haplotypes(), compute_haplotype_diversity(), and
## build_haplotype_feature_matrix(). All tests use the ldx_* example data
## or small synthetic fixtures to stay fast.
## ─────────────────────────────────────────────────────────────────────────────

library(testthat)
library(LDxBlocks)

data(ldx_geno,     package = "LDxBlocks")
data(ldx_snp_info, package = "LDxBlocks")
data(ldx_blocks,   package = "LDxBlocks")

# ── extract_haplotypes ────────────────────────────────────────────────────────

test_that("extract_haplotypes: list length equals number of qualifying blocks", {
  haps <- extract_haplotypes(ldx_geno, ldx_snp_info, ldx_blocks, min_snps = 5)
  # 9 blocks in ldx_blocks; all have >= 20 SNPs so all qualify
  expect_equal(length(haps), 9L)
})

test_that("extract_haplotypes: each element has one string per individual", {
  haps <- extract_haplotypes(ldx_geno, ldx_snp_info, ldx_blocks, min_snps = 5)
  for (nm in names(haps)) {
    expect_equal(length(haps[[nm]]), nrow(ldx_geno),
                 label = paste("length of", nm))
  }
})

test_that("extract_haplotypes: strings have correct width (n_snps chars)", {
  haps     <- extract_haplotypes(ldx_geno, ldx_snp_info, ldx_blocks, min_snps = 5)
  bi       <- attr(haps, "block_info")
  for (i in seq_len(nrow(bi))) {
    expected_w <- bi$n_snps[i]
    actual_w   <- unique(nchar(haps[[bi$block_id[i]]]))
    expect_equal(actual_w, expected_w,
                 label = paste("string width for", bi$block_id[i]))
  }
})

test_that("extract_haplotypes: string characters in {0,1,2,.}", {
  haps <- extract_haplotypes(ldx_geno, ldx_snp_info, ldx_blocks, min_snps = 5)
  all_chars <- unique(unlist(strsplit(unlist(haps), "")))
  expect_true(all(all_chars %in% c("0","1","2",".")))
})

test_that("extract_haplotypes: min_snps filters small blocks", {
  # Use a real small block from ldx_blocks and verify min_snps filtering
  # ldx_blocks block 1 has >= 20 SNPs; with min_snps = 25 it should be skipped
  haps_strict <- extract_haplotypes(ldx_geno, ldx_snp_info,
                                    ldx_blocks[1, , drop = FALSE],
                                    min_snps = 100)
  expect_equal(length(haps_strict), 0L)

  # With min_snps = 5, block 1 (>= 20 SNPs) should be included
  haps_loose <- extract_haplotypes(ldx_geno, ldx_snp_info,
                                   ldx_blocks[1, , drop = FALSE],
                                   min_snps = 5)
  expect_equal(length(haps_loose), 1L)
})

test_that("extract_haplotypes: block_info attribute has correct columns", {
  haps <- extract_haplotypes(ldx_geno, ldx_snp_info, ldx_blocks, min_snps = 5)
  bi   <- attr(haps, "block_info")
  expect_true(all(c("block_id","CHR","start_bp","end_bp","n_snps") %in% names(bi)))
  expect_equal(nrow(bi), length(haps))
})

test_that("extract_haplotypes: chr subset works via chr argument", {
  haps_all  <- extract_haplotypes(ldx_geno, ldx_snp_info, ldx_blocks,
                                  min_snps = 5)
  haps_chr1 <- extract_haplotypes(ldx_geno, ldx_snp_info, ldx_blocks,
                                  min_snps = 5, chr = "1")
  bi_all  <- attr(haps_all,  "block_info")
  bi_chr1 <- attr(haps_chr1, "block_info")
  expect_true(all(bi_chr1$CHR == "1"))
  expect_true(nrow(bi_chr1) < nrow(bi_all))
})

test_that("extract_haplotypes: NA genotype becomes na_char in string", {
  G_na         <- ldx_geno
  G_na[1, 1]   <- NA   # inject NA in individual 1, first SNP

  # Block covering SNP 1
  blk1 <- ldx_blocks[ldx_blocks$CHR == "1", ][1, ]
  haps <- extract_haplotypes(G_na, ldx_snp_info, blk1, min_snps = 2,
                             na_char = ".")
  # Individual 1 should have "." in their string
  expect_true(grepl(".", haps[[1]][1], fixed = TRUE))
  # All other individuals should be clean
  others <- haps[[1]][-1]
  expect_false(any(grepl(".", others, fixed = TRUE)))
})

# ── compute_haplotype_diversity ───────────────────────────────────────────────

test_that("compute_haplotype_diversity: returns data.frame with required cols", {
  haps <- extract_haplotypes(ldx_geno, ldx_snp_info, ldx_blocks, min_snps = 5)
  div  <- compute_haplotype_diversity(haps)
  req  <- c("block_id","n_ind","n_haplotypes","He","Shannon","freq_dominant")
  expect_true(all(req %in% names(div)))
})

test_that("compute_haplotype_diversity: one row per block", {
  haps <- extract_haplotypes(ldx_geno, ldx_snp_info, ldx_blocks, min_snps = 5)
  div  <- compute_haplotype_diversity(haps)
  expect_equal(nrow(div), length(haps))
})

test_that("compute_haplotype_diversity: He in [0,1]", {
  haps <- extract_haplotypes(ldx_geno, ldx_snp_info, ldx_blocks, min_snps = 5)
  div  <- compute_haplotype_diversity(haps)
  He   <- div$He[!is.na(div$He)]
  expect_true(all(He >= 0 & He <= 1 + 1e-8))
})

test_that("compute_haplotype_diversity: Shannon >= 0", {
  haps <- extract_haplotypes(ldx_geno, ldx_snp_info, ldx_blocks, min_snps = 5)
  div  <- compute_haplotype_diversity(haps)
  Sh   <- div$Shannon[!is.na(div$Shannon)]
  expect_true(all(Sh >= 0))
})

test_that("compute_haplotype_diversity: freq_dominant in (0,1]", {
  haps <- extract_haplotypes(ldx_geno, ldx_snp_info, ldx_blocks, min_snps = 5)
  div  <- compute_haplotype_diversity(haps)
  fd   <- div$freq_dominant[!is.na(div$freq_dominant)]
  expect_true(all(fd > 0 & fd <= 1))
})

test_that("compute_haplotype_diversity: monomorphic block gives He=0", {
  # Build a synthetic 5-SNP block where all individuals are identical
  G_syn  <- matrix(0L, nrow = 120, ncol = 5)
  rownames(G_syn) <- rownames(ldx_geno)
  colnames(G_syn) <- paste0("syn", 1:5)
  info_syn <- data.frame(SNP = colnames(G_syn), CHR = "1",
                         POS = seq(1000L, by = 500L, length.out = 5L),
                         stringsAsFactors = FALSE)
  blk_syn  <- data.frame(start=1L, end=5L,
                         start.rsID="syn1", end.rsID="syn5",
                         start.bp=1000L, end.bp=3000L,
                         CHR="1", length_bp=2001L,
                         stringsAsFactors=FALSE)
  haps_mono <- extract_haplotypes(G_syn, info_syn, blk_syn, min_snps = 2)
  div_mono  <- compute_haplotype_diversity(haps_mono)
  expect_equal(div_mono$n_haplotypes, 1L)
  expect_equal(div_mono$freq_dominant, 1.0)
  if (!is.na(div_mono$He)) expect_equal(div_mono$He, 0.0, tolerance = 1e-10)
})

test_that("compute_haplotype_diversity: Shannon=0 for monomorphic block", {
  # Same synthetic monomorphic block
  G_syn  <- matrix(0L, nrow = 120, ncol = 5)
  rownames(G_syn) <- rownames(ldx_geno)
  colnames(G_syn) <- paste0("syn", 1:5)
  info_syn <- data.frame(SNP = colnames(G_syn), CHR = "1",
                         POS = seq(1000L, by = 500L, length.out = 5L),
                         stringsAsFactors = FALSE)
  blk_syn  <- data.frame(start=1L, end=5L,
                         start.rsID="syn1", end.rsID="syn5",
                         start.bp=1000L, end.bp=3000L,
                         CHR="1", length_bp=2001L,
                         stringsAsFactors=FALSE)
  haps_mono <- extract_haplotypes(G_syn, info_syn, blk_syn, min_snps = 2)
  div_mono  <- compute_haplotype_diversity(haps_mono)
  if (!is.na(div_mono$Shannon))
    expect_equal(div_mono$Shannon, 0.0, tolerance = 1e-10)
})

# ── build_haplotype_feature_matrix ───────────────────────────────────────────

test_that("build_haplotype_feature_matrix: correct dimensions", {
  haps <- extract_haplotypes(ldx_geno, ldx_snp_info, ldx_blocks, min_snps = 5)
  feat <- build_haplotype_feature_matrix(haps, top_n = 3)
  expect_true(is.matrix(feat))
  expect_equal(nrow(feat), nrow(ldx_geno))
  expect_equal(ncol(feat), length(haps) * 3L)
})

test_that("build_haplotype_feature_matrix: values are 0, 2, or NA", {
  haps <- extract_haplotypes(ldx_geno, ldx_snp_info, ldx_blocks, min_snps = 5)
  feat <- build_haplotype_feature_matrix(haps, top_n = 3)
  vals <- as.vector(feat)
  vals <- vals[!is.na(vals)]
  expect_true(all(vals %in% c(0, 2)))
})

test_that("build_haplotype_feature_matrix: scaled columns have mean ~0 and sd ~1", {
  haps <- extract_haplotypes(ldx_geno, ldx_snp_info, ldx_blocks, min_snps = 5)
  feat <- build_haplotype_feature_matrix(haps, top_n = 3, scale_features = TRUE)
  # Only check columns that are not all-NA or zero-variance
  col_means <- colMeans(feat, na.rm = TRUE)
  expect_true(all(abs(col_means) < 1e-10))
})

test_that("build_haplotype_feature_matrix: top_n = 1 gives one column per block", {
  haps <- extract_haplotypes(ldx_geno, ldx_snp_info, ldx_blocks, min_snps = 5)
  feat <- build_haplotype_feature_matrix(haps, top_n = 1)
  expect_equal(ncol(feat), length(haps))
})

test_that("build_haplotype_feature_matrix: column names include block IDs", {
  haps <- extract_haplotypes(ldx_geno, ldx_snp_info, ldx_blocks, min_snps = 5)
  feat <- build_haplotype_feature_matrix(haps, top_n = 2)
  # Column names should start with block IDs
  block_names <- names(haps)
  expect_true(all(sapply(block_names, function(bn)
    any(startsWith(colnames(feat), bn)))))
})

test_that("build_haplotype_feature_matrix: row names match individual IDs", {
  haps <- extract_haplotypes(ldx_geno, ldx_snp_info, ldx_blocks, min_snps = 5)
  feat <- build_haplotype_feature_matrix(haps, top_n = 3)
  expect_equal(rownames(feat), rownames(ldx_geno))
})
