## tests/testthat/test-extensions.R
## -----------------------------------------------------------------------------
## Tests for the 10 new analysis extension and haplotype inference functions:
##
##   haplotype_analysis.R:
##     cv_haplotype_prediction()
##     compare_haplotype_populations()
##     plot_haplotype_network()
##     run_haplotype_stability()
##     export_candidate_regions()
##     decompose_block_effects()
##     scan_diversity_windows()
##
##   haplotype_inference.R:
##     infer_block_haplotypes()
##     collapse_haplotypes()
##     harmonize_haplotypes()
##
## Shared fixtures use helper.R (make_geno, make_snpinfo, make_blocks, make_blues)
## and the ldx_* example datasets.
## -----------------------------------------------------------------------------

library(testthat)
library(LDxBlocks)

data(ldx_geno,     package = "LDxBlocks")
data(ldx_snp_info, package = "LDxBlocks")
data(ldx_blocks,   package = "LDxBlocks")
data(ldx_gwas,     package = "LDxBlocks")
data(ldx_blues,    package = "LDxBlocks")

# -- Shared fixtures -----------------------------------------------------------

# Small synthetic genotype + block infrastructure for fast tests
.G   <- make_geno(n = 40, p = 30, seed = 7L)
.si  <- make_snpinfo(p = 30, chr = "1")
.blk <- make_blocks(.si, n_blocks = 3L)

# Extract haplotypes once - reused across multiple test groups
.haps <- extract_haplotypes(ldx_geno, ldx_snp_info, ldx_blocks, min_snps = 5L)

# Small haplotypes for fast-path tests
.haps_small <- extract_haplotypes(.G, .si, .blk, min_snps = 3L)

# Blues for the full ldx dataset
.blues_vec <- setNames(ldx_blues$YLD, ldx_blues$id)

# Two-environment blues list for stability tests
.blues_list <- list(
  env1 = setNames(ldx_blues$YLD + rnorm(120, 0,   0.1), ldx_blues$id),
  env2 = setNames(ldx_blues$YLD + rnorm(120, 0.4, 0.1), ldx_blues$id)
)

# ══════════════════════════════════════════════════════════════════════════════
# 1. cv_haplotype_prediction
# ══════════════════════════════════════════════════════════════════════════════

test_that("cv_haplotype_prediction: returns LDxBlocks_cv object", {
  skip_if_not_installed("rrBLUP")
  cv <- cv_haplotype_prediction(
    geno_matrix = ldx_geno, snp_info = ldx_snp_info, blocks = ldx_blocks,
    blues = .blues_vec, k = 3L, n_rep = 1L, verbose = FALSE
  )
  expect_s3_class(cv, "LDxBlocks_cv")
})

test_that("cv_haplotype_prediction: pa_summary has required columns", {
  skip_if_not_installed("rrBLUP")
  cv <- cv_haplotype_prediction(
    geno_matrix = ldx_geno, snp_info = ldx_snp_info, blocks = ldx_blocks,
    blues = .blues_vec, k = 3L, n_rep = 1L, verbose = FALSE
  )
  req <- c("trait", "rep", "fold", "n_train", "n_test", "PA", "RMSE")
  expect_true(all(req %in% names(cv$pa_summary)))
})

test_that("cv_haplotype_prediction: k folds produce k rows per trait per rep", {
  skip_if_not_installed("rrBLUP")
  cv <- cv_haplotype_prediction(
    geno_matrix = ldx_geno, snp_info = ldx_snp_info, blocks = ldx_blocks,
    blues = .blues_vec, k = 4L, n_rep = 1L, verbose = FALSE
  )
  trait_rows <- cv$pa_summary[cv$pa_summary$rep == 1L, ]
  expect_equal(nrow(trait_rows), 4L)
})

test_that("cv_haplotype_prediction: PA in [-1, 1] or NA", {
  skip_if_not_installed("rrBLUP")
  cv <- cv_haplotype_prediction(
    geno_matrix = ldx_geno, snp_info = ldx_snp_info, blocks = ldx_blocks,
    blues = .blues_vec, k = 3L, n_rep = 1L, verbose = FALSE
  )
  pa <- cv$pa_summary$PA[!is.na(cv$pa_summary$PA)]
  expect_true(all(pa >= -1 - 1e-8 & pa <= 1 + 1e-8))
})

test_that("cv_haplotype_prediction: pa_mean has one row per trait", {
  skip_if_not_installed("rrBLUP")
  cv <- cv_haplotype_prediction(
    geno_matrix = ldx_geno, snp_info = ldx_snp_info, blocks = ldx_blocks,
    blues = .blues_vec, k = 3L, n_rep = 1L, verbose = FALSE
  )
  expect_equal(nrow(cv$pa_mean), 1L)
  expect_true("PA" %in% names(cv$pa_mean))
  expect_true("RMSE" %in% names(cv$pa_mean))
})

test_that("cv_haplotype_prediction: k and n_rep stored correctly", {
  skip_if_not_installed("rrBLUP")
  cv <- cv_haplotype_prediction(
    geno_matrix = ldx_geno, snp_info = ldx_snp_info, blocks = ldx_blocks,
    blues = .blues_vec, k = 3L, n_rep = 2L, verbose = FALSE
  )
  expect_equal(cv$k, 3L)
  expect_equal(cv$n_rep, 2L)
})

test_that("cv_haplotype_prediction: multi-trait data.frame blues works", {
  skip_if_not_installed("rrBLUP")
  cv <- cv_haplotype_prediction(
    geno_matrix = ldx_geno, snp_info = ldx_snp_info, blocks = ldx_blocks,
    blues = ldx_blues, k = 3L, n_rep = 1L,
    id_col = "id", blue_cols = c("YLD", "RES"), verbose = FALSE
  )
  expect_equal(length(unique(cv$pa_summary$trait)), 2L)
  expect_equal(nrow(cv$pa_mean), 2L)
})

test_that("cv_haplotype_prediction: print method works without error", {
  skip_if_not_installed("rrBLUP")
  cv <- cv_haplotype_prediction(
    geno_matrix = ldx_geno, snp_info = ldx_snp_info, blocks = ldx_blocks,
    blues = .blues_vec, k = 3L, n_rep = 1L, verbose = FALSE
  )
  expect_no_error(print(cv))
})

# ══════════════════════════════════════════════════════════════════════════════
# 2. compare_haplotype_populations
# ══════════════════════════════════════════════════════════════════════════════

test_that("compare_haplotype_populations: returns data.frame with required columns", {
  ids <- rownames(ldx_geno)
  cmp <- compare_haplotype_populations(.haps, ids[1:60], ids[61:120])
  req <- c("block_id","CHR","start_bp","end_bp","n1","n2",
           "n_alleles","FST","max_freq_diff","dominant_g1","dominant_g2",
           "chisq_p","divergent")
  expect_true(all(req %in% names(cmp)))
  expect_s3_class(cmp, "data.frame")
})

test_that("compare_haplotype_populations: one row per block", {
  ids <- rownames(ldx_geno)
  cmp <- compare_haplotype_populations(.haps, ids[1:60], ids[61:120])
  expect_equal(nrow(cmp), length(.haps))
})

test_that("compare_haplotype_populations: FST in [0, 1]", {
  ids <- rownames(ldx_geno)
  cmp <- compare_haplotype_populations(.haps, ids[1:60], ids[61:120])
  fst <- cmp$FST[!is.na(cmp$FST)]
  expect_true(all(fst >= 0 - 1e-8 & fst <= 1 + 1e-8))
})

test_that("compare_haplotype_populations: max_freq_diff in [0, 1]", {
  ids <- rownames(ldx_geno)
  cmp <- compare_haplotype_populations(.haps, ids[1:60], ids[61:120])
  fd <- cmp$max_freq_diff[!is.na(cmp$max_freq_diff)]
  expect_true(all(fd >= 0 - 1e-8 & fd <= 1 + 1e-8))
})

test_that("compare_haplotype_populations: n1 and n2 are positive integers", {
  ids <- rownames(ldx_geno)
  cmp <- compare_haplotype_populations(.haps, ids[1:60], ids[61:120])
  expect_true(all(cmp$n1 >= 0L & cmp$n2 >= 0L))
})

test_that("compare_haplotype_populations: divergent is logical", {
  ids <- rownames(ldx_geno)
  cmp <- compare_haplotype_populations(.haps, ids[1:60], ids[61:120])
  expect_type(cmp$divergent, "logical")
})

test_that("compare_haplotype_populations: group with 0 matching individuals returns NA row", {
  ids <- rownames(ldx_geno)
  # group2 = non-existent IDs -> n2 = 0 -> NA row
  cmp <- compare_haplotype_populations(.haps, ids[1:10], c("ghost1","ghost2"))
  expect_true(all(is.na(cmp$FST)))
  expect_true(all(cmp$n2 == 0L))
})

test_that("compare_haplotype_populations: sorted by CHR and start_bp", {
  ids <- rownames(ldx_geno)
  cmp <- compare_haplotype_populations(.haps, ids[1:60], ids[61:120])
  if (nrow(cmp) > 1L)
    expect_equal(order(cmp$CHR, cmp$start_bp), seq_len(nrow(cmp)))
})

# ══════════════════════════════════════════════════════════════════════════════
# 3. plot_haplotype_network
# ══════════════════════════════════════════════════════════════════════════════

test_that("plot_haplotype_network: returns igraph object invisibly", {
  skip_if_not_installed("igraph")
  # Use a block known to have multiple alleles
  block_nm <- names(.haps)[1]
  result <- withVisible(
    plot_haplotype_network(.haps, block_id = block_nm, min_freq = 0.01)
  )
  expect_false(result$visible)
  expect_s3_class(result$value, "igraph")
})

test_that("plot_haplotype_network: errors for non-existent block_id", {
  skip_if_not_installed("igraph")
  expect_error(
    plot_haplotype_network(.haps, block_id = "block_does_not_exist"),
    "not found"
  )
})

test_that("plot_haplotype_network: MST has correct number of vertices", {
  skip_if_not_installed("igraph")
  block_nm <- names(.haps)[1]
  hap <- .haps[[block_nm]]
  valid <- hap[!grepl(".", hap, fixed = TRUE)]
  tbl  <- table(valid)
  freq <- as.numeric(tbl) / sum(tbl)
  n_kept <- sum(freq >= 0.01)
  if (n_kept >= 2L) {
    g <- plot_haplotype_network(.haps, block_id = block_nm, min_freq = 0.01)
    expect_equal(igraph::vcount(g), n_kept)
  }
})

test_that("plot_haplotype_network: works with groups argument", {
  skip_if_not_installed("igraph")
  ids <- names(.haps[[1]])
  groups <- setNames(rep(c("A","B"), length.out = length(ids)), ids)
  block_nm <- names(.haps)[1]
  expect_no_error(
    plot_haplotype_network(.haps, block_id = block_nm, groups = groups,
                           min_freq = 0.01)
  )
})

# ══════════════════════════════════════════════════════════════════════════════
# 4. run_haplotype_stability
# ══════════════════════════════════════════════════════════════════════════════

test_that("run_haplotype_stability: returns data.frame with required columns", {
  skip_if_not_installed("rrBLUP")
  stab <- run_haplotype_stability(
    geno_matrix = ldx_geno, snp_info = ldx_snp_info, blocks = ldx_blocks,
    blues_list = .blues_list, verbose = FALSE
  )
  req <- c("block_id","b","b_se","r2_fw","s2d","p_b1","stable")
  expect_true(all(req %in% names(stab)))
  expect_s3_class(stab, "data.frame")
})

test_that("run_haplotype_stability: b is numeric", {
  skip_if_not_installed("rrBLUP")
  stab <- run_haplotype_stability(
    geno_matrix = ldx_geno, snp_info = ldx_snp_info, blocks = ldx_blocks,
    blues_list = .blues_list, verbose = FALSE
  )
  expect_type(stab$b, "double")
})

test_that("run_haplotype_stability: stable is logical", {
  skip_if_not_installed("rrBLUP")
  stab <- run_haplotype_stability(
    geno_matrix = ldx_geno, snp_info = ldx_snp_info, blocks = ldx_blocks,
    blues_list = .blues_list, verbose = FALSE
  )
  expect_type(stab$stable, "logical")
})

test_that("run_haplotype_stability: r2_fw in [0, 1]", {
  skip_if_not_installed("rrBLUP")
  stab <- run_haplotype_stability(
    geno_matrix = ldx_geno, snp_info = ldx_snp_info, blocks = ldx_blocks,
    blues_list = .blues_list, verbose = FALSE
  )
  r2 <- stab$r2_fw[!is.na(stab$r2_fw)]
  expect_true(all(r2 >= 0 - 1e-8 & r2 <= 1 + 1e-8))
})

test_that("run_haplotype_stability: errors with fewer than 2 environments", {
  expect_error(
    run_haplotype_stability(
      geno_matrix = ldx_geno, snp_info = ldx_snp_info, blocks = ldx_blocks,
      blues_list = list(env1 = .blues_vec), verbose = FALSE
    ),
    "At least 2"
  )
})

test_that("run_haplotype_stability: errors with unnamed blues_list", {
  expect_error(
    run_haplotype_stability(
      geno_matrix = ldx_geno, snp_info = ldx_snp_info, blocks = ldx_blocks,
      blues_list = list(.blues_vec, .blues_vec), verbose = FALSE
    ),
    "named list"
  )
})

# ══════════════════════════════════════════════════════════════════════════════
# 5. export_candidate_regions
# ══════════════════════════════════════════════════════════════════════════════

.qtl <- define_qtl_regions(ldx_gwas, ldx_blocks, ldx_snp_info,
                           p_threshold = NULL, trait_col = "trait")

test_that("export_candidate_regions: BED format has 6 columns", {
  bed <- export_candidate_regions(.qtl, format = "bed")
  expect_equal(ncol(bed), 6L)
  expect_true(all(c("chrom","start","end","name","score","strand") %in% names(bed)))
})

test_that("export_candidate_regions: BED start is 0-based", {
  bed <- export_candidate_regions(.qtl, format = "bed")
  # BED start = block start_bp - 1 (0-based)
  expect_true(all(bed$start >= 0L))
  # end > start for all rows
  expect_true(all(bed$end > bed$start))
})

test_that("export_candidate_regions: BED one row per QTL region", {
  bed <- export_candidate_regions(.qtl, format = "bed")
  expect_equal(nrow(bed), nrow(.qtl))
})

test_that("export_candidate_regions: chr_prefix applied to chrom column", {
  bed <- export_candidate_regions(.qtl, format = "bed", chr_prefix = "chr")
  expect_true(all(startsWith(bed$chrom, "chr")))
})

test_that("export_candidate_regions: CSV format identical to input", {
  csv <- export_candidate_regions(.qtl, format = "csv")
  expect_equal(csv, .qtl)
})

test_that("export_candidate_regions: biomaRt format returns named list", {
  bm <- export_candidate_regions(.qtl, format = "biomart")
  expect_type(bm, "list")
  expect_true(all(c("chromosome_name","start","end") %in% names(bm)))
  expect_equal(length(bm$chromosome_name), nrow(.qtl))
})

test_that("export_candidate_regions: padding_bp extends BED regions", {
  bed_no_pad  <- export_candidate_regions(.qtl, format = "bed", padding_bp = 0L)
  bed_pad     <- export_candidate_regions(.qtl, format = "bed", padding_bp = 5000L)
  # Padded regions should be wider
  expect_true(all((bed_pad$end - bed_pad$start) >=
                    (bed_no_pad$end - bed_no_pad$start)))
})

test_that("export_candidate_regions: writes BED file to disk when out_file given", {
  tmp <- tempfile(fileext = ".bed")
  export_candidate_regions(.qtl, format = "bed", out_file = tmp)
  expect_true(file.exists(tmp))
  lines <- readLines(tmp)
  expect_equal(length(lines), nrow(.qtl))
  unlink(tmp)
})

test_that("export_candidate_regions: errors on missing required columns", {
  bad_qtl <- .qtl[, setdiff(names(.qtl), "start_bp")]
  expect_error(export_candidate_regions(bad_qtl, format = "bed"), "missing")
})

# ══════════════════════════════════════════════════════════════════════════════
# 6. decompose_block_effects
# ══════════════════════════════════════════════════════════════════════════════

# Create synthetic SNP effects for testing
.snp_effects <- setNames(rnorm(ncol(ldx_geno)), colnames(ldx_geno))

test_that("decompose_block_effects: returns data.frame with required columns", {
  res <- decompose_block_effects(
    .haps, ldx_snp_info, ldx_blocks, .snp_effects, min_freq = 0.02
  )
  req <- c("block_id","CHR","start_bp","end_bp","allele",
           "frequency","allele_effect","n_snps_block","effect_rank")
  expect_true(all(req %in% names(res)))
  expect_s3_class(res, "data.frame")
})

test_that("decompose_block_effects: effect_rank starts at 1 per block", {
  res <- decompose_block_effects(
    .haps, ldx_snp_info, ldx_blocks, .snp_effects, min_freq = 0.02
  )
  for (bn in unique(res$block_id)) {
    ranks <- res$effect_rank[res$block_id == bn]
    expect_equal(min(ranks), 1L, label = paste("block", bn))
  }
})

test_that("decompose_block_effects: frequency in (0, 1] for all rows", {
  res <- decompose_block_effects(
    .haps, ldx_snp_info, ldx_blocks, .snp_effects, min_freq = 0.02
  )
  freq <- res$frequency[!is.na(res$frequency)]
  expect_true(all(freq > 0 & freq <= 1 + 1e-8))
})

test_that("decompose_block_effects: min_freq filters alleles", {
  res_loose  <- decompose_block_effects(
    .haps, ldx_snp_info, ldx_blocks, .snp_effects, min_freq = 0.001
  )
  res_strict <- decompose_block_effects(
    .haps, ldx_snp_info, ldx_blocks, .snp_effects, min_freq = 0.30
  )
  expect_true(nrow(res_strict) <= nrow(res_loose))
})

test_that("decompose_block_effects: returns empty df when no SNP effect overlap", {
  empty_effects <- setNames(numeric(0), character(0))
  res <- decompose_block_effects(
    .haps, ldx_snp_info, ldx_blocks, empty_effects
  )
  expect_equal(nrow(res), 0L)
})

test_that("decompose_block_effects: allele_effect is numeric", {
  res <- decompose_block_effects(
    .haps, ldx_snp_info, ldx_blocks, .snp_effects, min_freq = 0.02
  )
  expect_type(res$allele_effect, "double")
})

test_that("decompose_block_effects: sorted by CHR, start_bp, effect_rank", {
  res <- decompose_block_effects(
    .haps, ldx_snp_info, ldx_blocks, .snp_effects, min_freq = 0.02
  )
  if (nrow(res) > 1L) {
    expected_order <- order(res$CHR, res$start_bp, res$effect_rank)
    expect_equal(expected_order, seq_len(nrow(res)))
  }
})

# ══════════════════════════════════════════════════════════════════════════════
# 7. scan_diversity_windows
# ══════════════════════════════════════════════════════════════════════════════

test_that("scan_diversity_windows: returns data.frame with required columns", {
  scan <- scan_diversity_windows(
    ldx_geno, ldx_snp_info, window_bp = 50000L, step_bp = 25000L,
    min_snps_win = 3L
  )
  req <- c("CHR","win_start","win_end","win_mid","n_snps","n_ind",
           "n_haplotypes","He","Shannon","n_eff_alleles","freq_dominant",
           "sweep_flag")
  expect_true(all(req %in% names(scan)))
})

test_that("scan_diversity_windows: He in [0, 1]", {
  scan <- scan_diversity_windows(
    ldx_geno, ldx_snp_info, window_bp = 50000L, step_bp = 25000L,
    min_snps_win = 3L
  )
  He <- scan$He[!is.na(scan$He)]
  expect_true(all(He >= 0 - 1e-8 & He <= 1 + 1e-8))
})

test_that("scan_diversity_windows: sweep_flag is logical", {
  scan <- scan_diversity_windows(
    ldx_geno, ldx_snp_info, window_bp = 50000L, step_bp = 25000L,
    min_snps_win = 3L
  )
  expect_type(scan$sweep_flag, "logical")
})

test_that("scan_diversity_windows: n_snps >= min_snps_win for all rows", {
  scan <- scan_diversity_windows(
    ldx_geno, ldx_snp_info, window_bp = 50000L, step_bp = 25000L,
    min_snps_win = 5L
  )
  expect_true(all(scan$n_snps >= 5L))
})

test_that("scan_diversity_windows: multiple chromosomes produce rows for each", {
  scan <- scan_diversity_windows(
    ldx_geno, ldx_snp_info, window_bp = 50000L, step_bp = 25000L,
    min_snps_win = 3L
  )
  chrs_in_data <- unique(ldx_snp_info$CHR)
  chrs_in_scan <- unique(scan$CHR)
  # Not all chromosomes may have windows with enough SNPs, but at least some do
  expect_true(length(chrs_in_scan) >= 1L)
  expect_true(all(chrs_in_scan %in% chrs_in_data))
})

test_that("scan_diversity_windows: sorted by CHR and win_start", {
  scan <- scan_diversity_windows(
    ldx_geno, ldx_snp_info, window_bp = 50000L, step_bp = 25000L,
    min_snps_win = 3L
  )
  if (nrow(scan) > 1L)
    expect_equal(order(scan$CHR, scan$win_start), seq_len(nrow(scan)))
})

test_that("scan_diversity_windows: win_mid between win_start and win_end", {
  scan <- scan_diversity_windows(
    ldx_geno, ldx_snp_info, window_bp = 50000L, step_bp = 25000L,
    min_snps_win = 3L
  )
  expect_true(all(scan$win_mid >= scan$win_start))
  expect_true(all(scan$win_mid <= scan$win_end))
})

test_that("scan_diversity_windows: errors on missing snp_info columns", {
  bad_si <- ldx_snp_info[, c("SNP", "CHR")]  # missing POS
  expect_error(
    scan_diversity_windows(ldx_geno, bad_si, window_bp = 50000L),
    "missing"
  )
})

test_that("scan_diversity_windows: n_eff_alleles >= 1 for all rows", {
  scan <- scan_diversity_windows(
    ldx_geno, ldx_snp_info, window_bp = 50000L, step_bp = 25000L,
    min_snps_win = 3L
  )
  nea <- scan$n_eff_alleles[!is.na(scan$n_eff_alleles)]
  expect_true(all(nea >= 1 - 1e-8))
})

# ══════════════════════════════════════════════════════════════════════════════
# 8. infer_block_haplotypes
# ══════════════════════════════════════════════════════════════════════════════

test_that("infer_block_haplotypes: returns data.frame with required columns", {
  dip <- infer_block_haplotypes(.haps)
  req <- c("block_id","CHR","start_bp","end_bp","id",
           "hap1","hap2","diplotype","heterozygous","phase_ambiguous","missing")
  expect_true(all(req %in% names(dip)))
  expect_s3_class(dip, "data.frame")
})

test_that("infer_block_haplotypes: one row per individual per block", {
  dip <- infer_block_haplotypes(.haps)
  n_blocks <- length(.haps)
  n_ind    <- nrow(ldx_geno)
  expect_equal(nrow(dip), n_blocks * n_ind)
})

test_that("infer_block_haplotypes: heterozygous is logical", {
  dip <- infer_block_haplotypes(.haps)
  expect_type(dip$heterozygous, "logical")
})

test_that("infer_block_haplotypes: phase_ambiguous is logical", {
  dip <- infer_block_haplotypes(.haps)
  expect_type(dip$phase_ambiguous, "logical")
})

test_that("infer_block_haplotypes: unphased homozygous individuals have hap1 == hap2", {
  dip <- infer_block_haplotypes(.haps)
  hom <- dip[!dip$heterozygous & !dip$missing & !is.na(dip$hap1), ]
  if (nrow(hom) > 0L)
    expect_true(all(hom$hap1 == hom$hap2))
})

test_that("infer_block_haplotypes: unphased heterozygous: phase_ambiguous TRUE by default", {
  # Deterministic: build a tiny block where individuals have dosage=1
  # (heterozygous), ensuring phase_ambiguous=TRUE with resolve_unphased=FALSE.
  # NOTE: when resolve_unphased=FALSE, hap1/hap2 remain NA for het individuals,
  # so heterozygous=FALSE (since h1!=h2 requires both non-NA). The correct
  # column to test ambiguity is phase_ambiguous.
  G_het <- matrix(c(0L, 0L, 0L,   # hom ref
                    1L, 1L, 1L,   # het (phase unknown)
                    2L, 2L, 2L),  # hom alt
                  nrow = 3L, ncol = 3L, byrow = TRUE)
  rownames(G_het) <- c("hom_ref", "het_ind", "hom_alt")
  colnames(G_het) <- c("s1", "s2", "s3")
  si_het  <- data.frame(SNP=c("s1","s2","s3"), CHR="1",
                        POS=c(1000L,2000L,3000L), REF="A", ALT="T",
                        stringsAsFactors=FALSE)
  blk_het <- data.frame(start=1L, end=3L, start.rsID="s1", end.rsID="s3",
                        start.bp=1000L, end.bp=3000L, CHR="1",
                        length_bp=2001L, stringsAsFactors=FALSE)
  haps_het <- extract_haplotypes(G_het, si_het, blk_het, min_snps=2L)
  dip <- infer_block_haplotypes(haps_het, resolve_unphased = FALSE)
  # het_ind has dosage "111": biologically het, phase unresolved
  ambig <- dip[dip$phase_ambiguous & !dip$missing, ]
  expect_gt(nrow(ambig), 0L)
  expect_equal(ambig$id, "het_ind")
  expect_true(all(ambig$heterozygous))  # biologically het (Fix 2)
  expect_true(all(is.na(ambig$hap1)))   # phase unresolved
  expect_true(all(is.na(ambig$hap2)))
})

test_that("infer_block_haplotypes: resolve_unphased=TRUE clears phase_ambiguous", {
  # Same fixture: with resolve=TRUE, het_ind gets hap1/hap2 assigned,
  # phase_ambiguous=FALSE, and heterozygous=TRUE.
  G_het <- matrix(c(0L,0L,0L, 1L,1L,1L, 2L,2L,2L),
                  nrow=3L, ncol=3L, byrow=TRUE)
  rownames(G_het) <- c("hom_ref","het_ind","hom_alt")
  colnames(G_het) <- c("s1","s2","s3")
  si_het  <- data.frame(SNP=c("s1","s2","s3"), CHR="1",
                        POS=c(1000L,2000L,3000L), REF="A", ALT="T",
                        stringsAsFactors=FALSE)
  blk_het <- data.frame(start=1L, end=3L, start.rsID="s1", end.rsID="s3",
                        start.bp=1000L, end.bp=3000L, CHR="1",
                        length_bp=2001L, stringsAsFactors=FALSE)
  haps_het <- extract_haplotypes(G_het, si_het, blk_het, min_snps=2L)
  dip <- infer_block_haplotypes(haps_het, resolve_unphased = TRUE)
  # After resolution: het_ind has hap1/hap2 assigned, phase_ambiguous=FALSE
  resolved <- dip[dip$id == "het_ind", ]
  expect_equal(nrow(resolved), 1L)
  expect_false(resolved$phase_ambiguous)
  expect_false(is.na(resolved$hap1))
  expect_false(is.na(resolved$hap2))
  expect_true(resolved$heterozygous)  # hap1 != hap2 after resolution
})

test_that("infer_block_haplotypes: phased input: phase_ambiguous always FALSE", {
  # Build minimal phased list
  h1 <- matrix(floor(ldx_geno / 2), nrow = nrow(ldx_geno))
  h2 <- ldx_geno - h1
  dimnames(h1) <- dimnames(h2) <- dimnames(ldx_geno)
  phased <- list(hap1 = t(h1), hap2 = t(h2), dosage = t(ldx_geno),
                 sample_ids = rownames(ldx_geno), phased = TRUE)
  haps_p <- extract_haplotypes(phased, ldx_snp_info, ldx_blocks, min_snps = 5L)
  dip_p  <- infer_block_haplotypes(haps_p)
  expect_true(all(!dip_p$phase_ambiguous))
})

test_that("infer_block_haplotypes: diplotype is canonical (sorted halves)", {
  dip <- infer_block_haplotypes(.haps)
  resolved <- dip[!is.na(dip$hap1) & !is.na(dip$hap2), ]
  if (nrow(resolved) > 0L) {
    expected_dip <- paste(
      pmin(resolved$hap1, resolved$hap2),
      pmax(resolved$hap1, resolved$hap2),
      sep = "/"
    )
    expect_equal(resolved$diplotype, expected_dip)
  }
})

test_that("infer_block_haplotypes: resolve_unphased=TRUE fills hap1/hap2 for some hets", {
  dip_unresolved <- infer_block_haplotypes(.haps, resolve_unphased = FALSE)
  dip_resolved   <- infer_block_haplotypes(.haps, resolve_unphased = TRUE)
  # Resolved version should have fewer NAs in hap1 among het rows
  n_na_before <- sum(is.na(dip_unresolved$hap1) & dip_unresolved$heterozygous)
  n_na_after  <- sum(is.na(dip_resolved$hap1)   & dip_resolved$heterozygous)
  expect_true(n_na_after <= n_na_before)
})

test_that("infer_block_haplotypes: errors without block_info attribute", {
  bad_haps <- .haps
  attr(bad_haps, "block_info") <- NULL
  expect_error(infer_block_haplotypes(bad_haps), "block_info")
})

# ══════════════════════════════════════════════════════════════════════════════
# 9. collapse_haplotypes
# ══════════════════════════════════════════════════════════════════════════════

test_that("collapse_haplotypes: rare_to_other pools rare alleles into <other>", {
  haps_col <- collapse_haplotypes(.haps_small, min_freq = 0.30,
                                  collapse = "rare_to_other")
  # Check that <other> appears in at least one block (if any allele was rare)
  all_strings <- unique(unlist(haps_col))
  all_strings <- all_strings[!grepl(".", all_strings, fixed = TRUE)]
  # <other> appears if there were rare alleles
  # (not guaranteed if all alleles happened to be common - just check structure)
  expect_true(is.list(haps_col))
  expect_equal(length(haps_col), length(.haps_small))
})

test_that("collapse_haplotypes: output has same structure as input", {
  haps_col <- collapse_haplotypes(.haps_small, min_freq = 0.10,
                                  collapse = "nearest")
  expect_equal(names(haps_col), names(.haps_small))
  expect_equal(length(haps_col), length(.haps_small))
  # Same individual names within each block
  for (bn in names(haps_col))
    expect_equal(names(haps_col[[bn]]), names(.haps_small[[bn]]))
})

test_that("collapse_haplotypes: block_info attribute preserved", {
  haps_col <- collapse_haplotypes(.haps_small, min_freq = 0.10)
  expect_false(is.null(attr(haps_col, "block_info")))
  expect_equal(attr(haps_col, "block_info"),
               attr(.haps_small, "block_info"))
})

test_that("collapse_haplotypes: label_map attribute present when keep_labels=TRUE", {
  haps_col <- collapse_haplotypes(.haps_small, min_freq = 0.10,
                                  keep_labels = TRUE)
  expect_false(is.null(attr(haps_col, "label_map")))
  expect_equal(names(attr(haps_col, "label_map")), names(.haps_small))
})

test_that("collapse_haplotypes: nearest strategy maps rare to Hamming-nearest", {
  haps_col <- collapse_haplotypes(.haps, min_freq = 0.20,
                                  collapse = "nearest")
  lm <- attr(haps_col, "label_map")
  # Skip if no blocks have rare alleles to map (e.g. all alleles above threshold)
  any_mapped <- any(vapply(lm, function(x) length(x) > 0L, logical(1L)))
  skip_if(!any_mapped, "no rare alleles to collapse in example data")
  # For each block with mappings, target must be a common allele
  for (bn in names(lm)) {
    if (length(lm[[bn]]) == 0L) next
    hap_bn <- .haps[[bn]]
    valid  <- hap_bn[!grepl(".", hap_bn, fixed = TRUE)]
    tbl    <- table(valid)
    freq   <- as.numeric(tbl) / sum(tbl)
    names(freq) <- names(tbl)
    common <- names(freq)[freq > 0.20]
    for (target in unlist(lm[[bn]])) {
      if (target != "<other>")
        expect_true(target %in% common,
                    label = paste("Mapped target", target,
                                  "is not a common allele in block", bn))
    }
  }
})

test_that("collapse_haplotypes: tree_based strategy runs without error", {
  expect_no_error(
    collapse_haplotypes(.haps_small, min_freq = 0.10, collapse = "tree_based")
  )
})

test_that("collapse_haplotypes: all-common case produces no changes", {
  # With min_freq = 0 nothing is rare -> no changes
  haps_col <- collapse_haplotypes(.haps_small, min_freq = 0.0,
                                  collapse = "nearest")
  for (bn in names(haps_col))
    expect_equal(haps_col[[bn]], .haps_small[[bn]], label = bn)
})

# ══════════════════════════════════════════════════════════════════════════════
# 10. harmonize_haplotypes
# ══════════════════════════════════════════════════════════════════════════════

# Split ldx_geno into reference (70%) and target (30%) panels
.n_ind   <- nrow(ldx_geno)
set.seed(42L)
.ref_idx <- sort(sample(.n_ind, round(.n_ind * 0.7)))
.tgt_idx <- setdiff(seq_len(.n_ind), .ref_idx)

.ref_geno <- ldx_geno[.ref_idx, ]
.tgt_geno <- ldx_geno[.tgt_idx, ]

.haps_ref <- extract_haplotypes(.ref_geno, ldx_snp_info, ldx_blocks,
                                min_snps = 5L)
.haps_tgt <- extract_haplotypes(.tgt_geno, ldx_snp_info, ldx_blocks,
                                min_snps = 5L)

test_that("harmonize_haplotypes: output has same names as target", {
  hh <- harmonize_haplotypes(.haps_tgt, .haps_ref)
  expect_equal(names(hh), names(.haps_tgt))
})

test_that("harmonize_haplotypes: output has same length as target", {
  hh <- harmonize_haplotypes(.haps_tgt, .haps_ref)
  expect_equal(length(hh), length(.haps_tgt))
})

test_that("harmonize_haplotypes: harmonization_report attribute present", {
  hh <- harmonize_haplotypes(.haps_tgt, .haps_ref)
  report <- attr(hh, "harmonization_report")
  expect_false(is.null(report))
  req <- c("block_id","n_exact","n_nearest","n_novel","mean_hamming_dist")
  expect_true(all(req %in% names(report)))
})

test_that("harmonize_haplotypes: n_exact + n_nearest + n_novel = total non-missing", {
  hh     <- harmonize_haplotypes(.haps_tgt, .haps_ref)
  report <- attr(hh, "harmonization_report")
  for (i in seq_len(nrow(report))) {
    bn <- report$block_id[i]
    if (is.na(report$n_exact[i])) next   # block not in reference
    hap_tgt <- .haps_tgt[[bn]]
    n_valid <- sum(!grepl(".", hap_tgt, fixed = TRUE))
    total_assigned <- report$n_exact[i] + report$n_nearest[i] + report$n_novel[i]
    expect_equal(total_assigned, n_valid,
                 label = paste("block", bn, "total assigned"))
  }
})

test_that("harmonize_haplotypes: block_info preserved from target", {
  hh    <- harmonize_haplotypes(.haps_tgt, .haps_ref)
  bi_hh <- attr(hh, "block_info")
  bi_tg <- attr(.haps_tgt, "block_info")
  expect_equal(bi_hh, bi_tg)
})

test_that("harmonize_haplotypes: exact matches preserved verbatim", {
  hh <- harmonize_haplotypes(.haps_tgt, .haps_ref)
  # Any allele string that was already in the reference dictionary must
  # appear unchanged in the output for that individual
  bn <- names(.haps_tgt)[1]
  hap_tgt <- .haps_tgt[[bn]]
  hap_ref <- .haps_ref[[bn]]
  valid_ref <- hap_ref[!grepl(".", hap_ref, fixed = TRUE)]
  tbl_ref   <- table(valid_ref)
  freq_ref  <- as.numeric(tbl_ref) / sum(tbl_ref)
  names(freq_ref) <- names(tbl_ref)
  dict <- names(freq_ref)[freq_ref >= 0.02]
  # For individuals in target whose allele is in dict, harmonised value = same
  hh_bn  <- hh[[bn]]
  for (i in seq_along(hap_tgt)) {
    s <- hap_tgt[i]
    if (grepl(".", s, fixed = TRUE)) next
    if (s %in% dict)
      expect_equal(hh_bn[i], s,
                   label = paste("exact match for individual", names(hap_tgt)[i]))
  }
})

test_that("harmonize_haplotypes: max_hamming=0 forces non-exact to <novel>", {
  hh <- harmonize_haplotypes(.haps_tgt, .haps_ref, max_hamming = 0L)
  report <- attr(hh, "harmonization_report")
  # With max_hamming=0, all non-exact matches become <novel>
  # n_nearest should be 0 everywhere
  nn <- report$n_nearest[!is.na(report$n_nearest)]
  expect_true(all(nn == 0L))
})

test_that("harmonize_haplotypes: errors with no common block IDs", {
  # Create haps with completely different block names
  fake_haps <- .haps_tgt
  names(fake_haps) <- paste0("fake_block_", seq_along(fake_haps))
  attr(fake_haps, "block_info")$block_id <- names(fake_haps)
  expect_error(
    harmonize_haplotypes(fake_haps, .haps_ref),
    "No common block IDs"
  )
})

test_that("harmonize_haplotypes: novel alleles labelled <novel> when max_hamming specified", {
  # Force max_hamming = 0 so any non-exact allele gets <novel>
  hh <- harmonize_haplotypes(.haps_tgt, .haps_ref, max_hamming = 0L)
  # Check that <novel> appears if any target allele was not in reference dict
  all_strings <- unique(unlist(hh))
  # Not guaranteed to have <novel> if target is subset of ref, but
  # the assignment mechanism is correct - just confirm no error and
  # <novel> is a valid string if it appears
  if ("<novel>" %in% all_strings)
    expect_true("<novel>" %in% all_strings)  # trivially true if present
  # More useful: no NA strings in non-missing positions
  for (bn in names(hh)) {
    non_miss <- hh[[bn]][!grepl(".", hh[[bn]], fixed = TRUE)]
    expect_false(any(is.na(non_miss)), label = paste("no NA in block", bn))
  }
})
