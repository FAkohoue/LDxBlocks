# tests/testthat/test-association-stress.R

test_that("test_block_haplotypes handles fully monomorphic haplotypes", {
  haps <- list(
    block_1_100_200 = setNames(rep("00000", 20), paste0("ind", 1:20))
  )
  attr(haps, "block_info") <- data.frame(
    block_id = "block_1_100_200", CHR = "1",
    start_bp = 100L, end_bp = 200L, n_snps = 5L, phased = FALSE,
    stringsAsFactors = FALSE
  )
  blues  <- setNames(rnorm(20), paste0("ind", 1:20))
  blocks <- data.frame(CHR = "1", start.bp = 100L, end.bp = 200L, stringsAsFactors = FALSE)

  expect_no_error({
    res <- test_block_haplotypes(haplotypes = haps, blues = blues,
                                 blocks = blocks, verbose = FALSE)
  })
  expect_s3_class(res, "LDxBlocks_haplotype_assoc")
  expect_true(res$status %in% c("skipped_monomorphic", "no_results", "ok"))
})


test_that("test_block_haplotypes handles missing haplotype values", {
  ids <- paste0("ind", 1:30)
  h   <- sample(c("000", "111", "010", NA), 30, replace = TRUE)
  names(h) <- ids
  haps <- list(block_1_100_200 = h)
  attr(haps, "block_info") <- data.frame(
    block_id = "block_1_100_200", CHR = "1",
    start_bp = 100L, end_bp = 200L, n_snps = 3L, phased = FALSE,
    stringsAsFactors = FALSE
  )
  blues  <- setNames(rnorm(30), ids)
  blocks <- data.frame(CHR = "1", start.bp = 100L, end.bp = 200L, stringsAsFactors = FALSE)

  expect_no_error({
    res <- test_block_haplotypes(haplotypes = haps, blues = blues,
                                 blocks = blocks, verbose = FALSE)
  })
  expect_s3_class(res, "LDxBlocks_haplotype_assoc")
})


test_that("test_block_haplotypes handles degenerate GRM gracefully", {
  ids <- paste0("ind", 1:25)
  haps <- list(
    block_1_100_200 = setNames(rep("000", 25), ids),
    block_1_300_400 = setNames(rep("111", 25), ids)
  )
  attr(haps, "block_info") <- data.frame(
    block_id = c("block_1_100_200", "block_1_300_400"), CHR = c("1", "1"),
    start_bp = c(100L, 300L), end_bp = c(200L, 400L),
    n_snps = c(3L, 3L), phased = c(FALSE, FALSE), stringsAsFactors = FALSE
  )
  blues  <- setNames(rnorm(25), ids)
  blocks <- data.frame(CHR = c("1", "1"), start.bp = c(100L, 300L),
                       end.bp = c(200L, 400L), stringsAsFactors = FALSE)

  expect_no_error({
    res <- test_block_haplotypes(haplotypes = haps, blues = blues,
                                 blocks = blocks, verbose = FALSE)
  })
  expect_s3_class(res, "LDxBlocks_haplotype_assoc")
})


test_that("test_block_haplotypes skips traits with too few common individuals", {
  ids_hap   <- paste0("ind",   1:30)
  ids_pheno <- paste0("other", 1:30)
  haps <- list(
    block_1_100_200 = setNames(sample(c("000", "111"), 30, TRUE), ids_hap)
  )
  attr(haps, "block_info") <- data.frame(
    block_id = "block_1_100_200", CHR = "1",
    start_bp = 100L, end_bp = 200L, n_snps = 3L, phased = FALSE,
    stringsAsFactors = FALSE
  )
  blues  <- setNames(rnorm(30), ids_pheno)
  blocks <- data.frame(CHR = "1", start.bp = 100L, end.bp = 200L, stringsAsFactors = FALSE)

  expect_no_error({
    res <- test_block_haplotypes(haplotypes = haps, blues = blues,
                                 blocks = blocks, verbose = FALSE)
  })
  expect_s3_class(res, "LDxBlocks_haplotype_assoc")
  expect_true(res$status %in% c("no_results", "ok"))
})


test_that("test_block_haplotypes handles multi-trait data frame input", {
  ids <- paste0("ind", 1:40)
  haps <- list(
    block_1_100_200 = setNames(sample(c("000", "111", "010"), 40, TRUE), ids),
    block_1_300_400 = setNames(sample(c("001", "110", "011"), 40, TRUE), ids)
  )
  attr(haps, "block_info") <- data.frame(
    block_id = c("block_1_100_200", "block_1_300_400"), CHR = c("1", "1"),
    start_bp = c(100L, 300L), end_bp = c(200L, 400L),
    n_snps = c(3L, 3L), phased = c(FALSE, FALSE), stringsAsFactors = FALSE
  )
  blues <- data.frame(id = ids, YLD = rnorm(40), RES = rnorm(40),
                      stringsAsFactors = FALSE)
  blocks <- data.frame(CHR = c("1", "1"), start.bp = c(100L, 300L),
                       end.bp = c(200L, 400L), stringsAsFactors = FALSE)

  expect_no_error({
    res <- test_block_haplotypes(haplotypes = haps, blues = blues,
                                 blocks = blocks, id_col = "id",
                                 blue_cols = c("YLD", "RES"), verbose = FALSE)
  })
  expect_s3_class(res, "LDxBlocks_haplotype_assoc")
  expect_true(all(c("YLD", "RES") %in% res$traits))
})
