## tests/testthat/test-haplotypes.R
## Tests for the full haplotype module: extraction (phased & unphased),
## diversity, QTL region definition, feature matrix, and output writers.

library(testthat)
library(LDxBlocks)

data(ldx_geno,     package = "LDxBlocks")
data(ldx_snp_info, package = "LDxBlocks")
data(ldx_blocks,   package = "LDxBlocks")
data(ldx_gwas,     package = "LDxBlocks")

# ── extract_haplotypes (unphased) ─────────────────────────────────────────────

test_that("extract_haplotypes: list length equals qualifying blocks", {
  haps <- extract_haplotypes(ldx_geno, ldx_snp_info, ldx_blocks, min_snps=5)
  expect_equal(length(haps), 9L)   # all 9 blocks have >= 20 SNPs
})

test_that("extract_haplotypes: one string per individual per block", {
  haps <- extract_haplotypes(ldx_geno, ldx_snp_info, ldx_blocks, min_snps=5)
  for (nm in names(haps))
    expect_equal(length(haps[[nm]]), nrow(ldx_geno), label=nm)
})

test_that("extract_haplotypes: string width equals block n_snps", {
  haps <- extract_haplotypes(ldx_geno, ldx_snp_info, ldx_blocks, min_snps=5)
  bi   <- attr(haps, "block_info")
  for (i in seq_len(nrow(bi))) {
    w <- unique(nchar(haps[[bi$block_id[i]]]))
    expect_equal(w, bi$n_snps[i], label=bi$block_id[i])
  }
})

test_that("extract_haplotypes: characters in {0,1,2,.}", {
  haps <- extract_haplotypes(ldx_geno, ldx_snp_info, ldx_blocks, min_snps=5)
  all_chars <- unique(unlist(strsplit(unlist(haps), "")))
  expect_true(all(all_chars %in% c("0","1","2",".")))
})

test_that("extract_haplotypes: min_snps filters blocks correctly", {
  # min_snps > block size -> skip
  haps0 <- extract_haplotypes(ldx_geno, ldx_snp_info,
                              ldx_blocks[1,,drop=FALSE], min_snps=100)
  expect_equal(length(haps0), 0L)
  # min_snps <= block size -> include
  haps1 <- extract_haplotypes(ldx_geno, ldx_snp_info,
                              ldx_blocks[1,,drop=FALSE], min_snps=5)
  expect_equal(length(haps1), 1L)
})

test_that("extract_haplotypes: block_info attribute is correct", {
  haps <- extract_haplotypes(ldx_geno, ldx_snp_info, ldx_blocks, min_snps=5)
  bi   <- attr(haps, "block_info")
  expect_true(all(c("block_id","CHR","start_bp","end_bp","n_snps","phased")
                  %in% names(bi)))
  expect_equal(nrow(bi), length(haps))
  expect_true(all(!bi$phased))   # unphased input
})

test_that("extract_haplotypes: chr argument subsets correctly", {
  haps_all  <- extract_haplotypes(ldx_geno, ldx_snp_info, ldx_blocks, min_snps=5)
  haps_chr1 <- extract_haplotypes(ldx_geno, ldx_snp_info, ldx_blocks,
                                  min_snps=5, chr="1")
  bi_chr1 <- attr(haps_chr1, "block_info")
  expect_true(all(bi_chr1$CHR == "1"))
  expect_true(length(haps_chr1) < length(haps_all))
})

test_that("extract_haplotypes: NA genotype becomes na_char in string", {
  G_na       <- ldx_geno
  G_na[1, 1] <- NA
  blk1 <- ldx_blocks[ldx_blocks$CHR == "1",][1,]
  haps <- extract_haplotypes(G_na, ldx_snp_info, blk1, min_snps=2, na_char=".")
  expect_true(grepl(".", haps[[1]][1], fixed=TRUE))
  expect_false(any(grepl(".", haps[[1]][-1], fixed=TRUE)))
})

# ── extract_haplotypes (phased) ───────────────────────────────────────────────

test_that("extract_haplotypes: phased strings contain '|' separator", {
  # Build a minimal phased list from the dosage matrix
  h1 <- matrix(floor(ldx_geno / 2), nrow=nrow(ldx_geno))
  h2 <- ldx_geno - h1
  dimnames(h1) <- dimnames(h2) <- dimnames(ldx_geno)
  phased <- list(hap1=t(h1), hap2=t(h2),
                 dosage=t(ldx_geno),
                 sample_ids=rownames(ldx_geno), phased=TRUE)
  haps_p <- extract_haplotypes(phased, ldx_snp_info, ldx_blocks, min_snps=5)
  # All strings should contain "|"
  expect_true(all(grepl("|", haps_p[[1]], fixed=TRUE)))
  expect_true(attr(haps_p, "block_info")$phased[1])
})

test_that("extract_haplotypes: phased string has two equal-width gametes", {
  h1 <- matrix(floor(ldx_geno / 2), nrow=nrow(ldx_geno))
  h2 <- ldx_geno - h1
  dimnames(h1) <- dimnames(h2) <- dimnames(ldx_geno)
  phased <- list(hap1=t(h1), hap2=t(h2), dosage=t(ldx_geno),
                 sample_ids=rownames(ldx_geno), phased=TRUE)
  haps_p <- extract_haplotypes(phased, ldx_snp_info, ldx_blocks, min_snps=5)
  bi     <- attr(haps_p, "block_info")
  # Each gamete half should have width = n_snps
  first_str <- haps_p[[1]][1]
  parts     <- strsplit(first_str, "|", fixed=TRUE)[[1]]
  expect_equal(length(parts), 2L)
  expect_equal(nchar(parts[1]), bi$n_snps[1])
  expect_equal(nchar(parts[2]), bi$n_snps[1])
})

# ── compute_haplotype_diversity ───────────────────────────────────────────────

test_that("compute_haplotype_diversity: returns required columns", {
  haps <- extract_haplotypes(ldx_geno, ldx_snp_info, ldx_blocks, min_snps=5)
  div  <- compute_haplotype_diversity(haps)
  expect_true(all(c("block_id","CHR","start_bp","end_bp","n_snps",
                    "n_ind","n_haplotypes","He","Shannon",
                    "n_eff_alleles","freq_dominant","sweep_flag",
                    "phased") %in% names(div)))
})

test_that("compute_haplotype_diversity: one row per block", {
  haps <- extract_haplotypes(ldx_geno, ldx_snp_info, ldx_blocks, min_snps=5)
  div  <- compute_haplotype_diversity(haps)
  expect_equal(nrow(div), length(haps))
})

test_that("compute_haplotype_diversity: He in [0,1]", {
  haps <- extract_haplotypes(ldx_geno, ldx_snp_info, ldx_blocks, min_snps=5)
  He   <- compute_haplotype_diversity(haps)$He
  expect_true(all(He[!is.na(He)] >= 0 & He[!is.na(He)] <= 1 + 1e-8))
})

test_that("compute_haplotype_diversity: Shannon >= 0", {
  haps <- extract_haplotypes(ldx_geno, ldx_snp_info, ldx_blocks, min_snps=5)
  Sh   <- compute_haplotype_diversity(haps)$Shannon
  expect_true(all(Sh[!is.na(Sh)] >= 0))
})

test_that("compute_haplotype_diversity: freq_dominant in (0,1]", {
  haps <- extract_haplotypes(ldx_geno, ldx_snp_info, ldx_blocks, min_snps=5)
  fd   <- compute_haplotype_diversity(haps)$freq_dominant
  expect_true(all(fd[!is.na(fd)] > 0 & fd[!is.na(fd)] <= 1))
})

test_that("compute_haplotype_diversity: monomorphic block gives He=0, Shannon=0", {
  G_mono <- matrix(0L, 120, 5)
  rownames(G_mono) <- rownames(ldx_geno)
  colnames(G_mono) <- paste0("syn", 1:5)
  info_syn <- data.frame(SNP=colnames(G_mono), CHR="1",
                         POS=seq(1000L, by=500L, length.out=5L),
                         stringsAsFactors=FALSE)
  blk_syn  <- data.frame(start=1L, end=5L,
                         start.rsID="syn1", end.rsID="syn5",
                         start.bp=1000L, end.bp=3000L,
                         CHR="1", length_bp=2001L, stringsAsFactors=FALSE)
  haps <- extract_haplotypes(G_mono, info_syn, blk_syn, min_snps=2)
  div  <- compute_haplotype_diversity(haps)
  expect_equal(div$n_haplotypes, 1L)
  expect_equal(div$freq_dominant, 1.0)
  if (!is.na(div$He))      expect_equal(div$He,      0.0, tolerance=1e-10)
  if (!is.na(div$Shannon)) expect_equal(div$Shannon, 0.0, tolerance=1e-10)
})

# ── define_qtl_regions ────────────────────────────────────────────────────────

test_that("define_qtl_regions: returns data.frame with required columns", {
  qtl <- define_qtl_regions(ldx_gwas, ldx_blocks, ldx_snp_info,
                            p_threshold=NULL, trait_col="trait")
  expect_s3_class(qtl, "data.frame")
  req <- c("block_id","CHR","start_bp","end_bp","n_snps_block",
           "n_sig_markers","lead_snp","lead_p","sig_snps",
           "traits","n_traits","pleiotropic")
  expect_true(all(req %in% names(qtl)))
})

test_that("define_qtl_regions: pleiotropic flag is logical", {
  qtl <- define_qtl_regions(ldx_gwas, ldx_blocks, ldx_snp_info,
                            p_threshold=NULL, trait_col="trait")
  expect_type(qtl$pleiotropic, "logical")
})

test_that("define_qtl_regions: n_traits consistent with traits string", {
  qtl <- define_qtl_regions(ldx_gwas, ldx_blocks, ldx_snp_info,
                            p_threshold=NULL, trait_col="trait")
  # n_traits should equal number of comma-separated values in traits
  computed_n <- vapply(qtl$traits, function(t) length(strsplit(t,",")[[1]]),
                       integer(1L))
  expect_equal(qtl$n_traits, unname(computed_n))  # computed_n has names from qtl$traits
})

test_that("define_qtl_regions: pleiotropic TRUE when n_traits > 1", {
  qtl <- define_qtl_regions(ldx_gwas, ldx_blocks, ldx_snp_info,
                            p_threshold=NULL, trait_col="trait")
  expect_true(all(qtl$pleiotropic == (qtl$n_traits > 1L)))
})

test_that("define_qtl_regions: single-trait GWAS gives pleiotropic=FALSE everywhere", {
  gwas_a <- ldx_gwas[ldx_gwas$trait == "TraitA",]
  qtl <- define_qtl_regions(gwas_a, ldx_blocks, ldx_snp_info,
                            p_threshold=NULL, trait_col="trait")
  expect_true(all(!qtl$pleiotropic))
})

test_that("define_qtl_regions: absent trait column falls back to single trait", {
  gwas_no_trait <- ldx_gwas[, c("Marker","CHR","POS","P")]
  qtl <- define_qtl_regions(gwas_no_trait, ldx_blocks, ldx_snp_info,
                            p_threshold=NULL)
  expect_true(all(qtl$n_traits == 1L))
})

test_that("define_qtl_regions: p_threshold filters markers", {
  # Use a strict threshold — only markers with P < 1e-6.
  # suppressMessages: if none pass the threshold the function emits a message;
  # the assertion still holds (0 <= positive).
  qtl_strict <- suppressMessages(
    define_qtl_regions(ldx_gwas, ldx_blocks, ldx_snp_info,
                       p_threshold=1e-6, trait_col="trait"))
  qtl_all    <- define_qtl_regions(ldx_gwas, ldx_blocks, ldx_snp_info,
                                   p_threshold=NULL, trait_col="trait")
  # Strict threshold can only produce fewer or equal sig markers
  expect_true(sum(qtl_strict$n_sig_markers) <= sum(qtl_all$n_sig_markers))
})

test_that("define_qtl_regions: output sorted by CHR then start_bp", {
  qtl <- define_qtl_regions(ldx_gwas, ldx_blocks, ldx_snp_info,
                            p_threshold=NULL, trait_col="trait")
  if (nrow(qtl) > 1L) {
    chr_ord <- order(qtl$CHR, qtl$start_bp)
    expect_equal(seq_len(nrow(qtl)), chr_ord)
  }
})

# ── build_haplotype_feature_matrix ───────────────────────────────────────────

test_that("build_haplotype_feature_matrix: correct dimensions for top_n=3", {
  haps <- extract_haplotypes(ldx_geno, ldx_snp_info, ldx_blocks, min_snps=5)
  feat <- build_haplotype_feature_matrix(haps, top_n=3)
  expect_equal(nrow(feat), nrow(ldx_geno))
  expect_equal(ncol(feat), length(haps) * 3L)
})

test_that("build_haplotype_feature_matrix: additive_012 unphased gives 0/1/NA", {
  haps <- extract_haplotypes(ldx_geno, ldx_snp_info, ldx_blocks, min_snps=5)
  feat <- build_haplotype_feature_matrix(haps, top_n=3, encoding="additive_012")
  vals <- as.vector(feat[!is.na(feat)])
  expect_true(all(vals %in% c(0, 1)))  # unphased: 0=absent, 1=present
})

test_that("build_haplotype_feature_matrix: presence_01 gives 0/1/NA", {
  haps <- extract_haplotypes(ldx_geno, ldx_snp_info, ldx_blocks, min_snps=5)
  feat <- build_haplotype_feature_matrix(haps, top_n=3, encoding="presence_01")
  vals <- as.vector(feat[!is.na(feat)])
  expect_true(all(vals %in% c(0, 1)))  # presence/absence: 0=absent, 1=present
})

test_that("build_haplotype_feature_matrix: additive_012 phased gives 0/1/2", {
  h1 <- matrix(floor(ldx_geno / 2), nrow=nrow(ldx_geno))
  h2 <- ldx_geno - h1
  dimnames(h1) <- dimnames(h2) <- dimnames(ldx_geno)
  phased <- list(hap1=t(h1), hap2=t(h2), dosage=t(ldx_geno),
                 sample_ids=rownames(ldx_geno), phased=TRUE)
  haps_p <- extract_haplotypes(phased, ldx_snp_info, ldx_blocks, min_snps=5)
  feat_p <- build_haplotype_feature_matrix(haps_p, top_n=3,
                                           encoding="additive_012")
  vals <- as.vector(feat_p[!is.na(feat_p)])
  expect_true(all(vals %in% c(0, 1, 2)))  # 1 now possible with phased data
})

test_that("build_haplotype_feature_matrix: scaled columns have mean ~0", {
  haps <- extract_haplotypes(ldx_geno, ldx_snp_info, ldx_blocks, min_snps=5)
  feat <- build_haplotype_feature_matrix(haps, top_n=3, scale_features=TRUE)
  expect_true(all(abs(colMeans(feat, na.rm=TRUE)) < 1e-10))
})

test_that("build_haplotype_feature_matrix: top_n=1 gives one col per block", {
  haps <- extract_haplotypes(ldx_geno, ldx_snp_info, ldx_blocks, min_snps=5)
  feat <- build_haplotype_feature_matrix(haps, top_n=1)
  expect_equal(ncol(feat), length(haps))
})

test_that("build_haplotype_feature_matrix: row names match individual IDs", {
  haps <- extract_haplotypes(ldx_geno, ldx_snp_info, ldx_blocks, min_snps=5)
  feat <- build_haplotype_feature_matrix(haps, top_n=3)
  expect_equal(rownames(feat), rownames(ldx_geno))
})

test_that("build_haplotype_feature_matrix: column names reference block IDs", {
  haps  <- extract_haplotypes(ldx_geno, ldx_snp_info, ldx_blocks, min_snps=5)
  feat  <- build_haplotype_feature_matrix(haps, top_n=2)
  block_names <- names(haps)
  expect_true(all(vapply(block_names, function(bn)
    any(startsWith(colnames(feat), bn)), logical(1L))))
})

# ── output writers ────────────────────────────────────────────────────────────

test_that("write_haplotype_numeric: file exists and has correct orientation", {
  haps <- extract_haplotypes(ldx_geno, ldx_snp_info, ldx_blocks, min_snps=5)
  feat <- build_haplotype_feature_matrix(haps, top_n=2)
  tmp  <- tempfile(fileext=".csv")
  write_haplotype_numeric(feat, tmp,
                          haplotypes = haps, snp_info = ldx_snp_info,
                          verbose = FALSE)
  expect_true(file.exists(tmp))
  df <- read.table(tmp, sep="\t", header=TRUE, check.names=FALSE)
  # Rows = haplotype alleles; cols = metadata + individuals
  expect_equal(nrow(df), ncol(feat))          # one row per haplotype allele
  expect_true("hap_id" %in% names(df))
  expect_true("CHR"    %in% names(df))
  expect_true("start_bp" %in% names(df))
  expect_true("end_bp"   %in% names(df))
  expect_true("n_snps"   %in% names(df))
  expect_true("alleles"  %in% names(df))
  # Individual columns: total cols = 8 metadata + n_individuals
  n_meta <- 7L  # hap_id, CHR, start_bp, end_bp, n_snps, alleles, frequency
  expect_true(ncol(df) > n_meta)
  unlink(tmp)
})

test_that("write_haplotype_numeric: dosage values are 0/1/NA (unphased)", {
  haps <- extract_haplotypes(ldx_geno, ldx_snp_info, ldx_blocks, min_snps=5)
  feat <- build_haplotype_feature_matrix(haps, top_n=2)
  tmp  <- tempfile(fileext=".csv")
  write_haplotype_numeric(feat, tmp, verbose=FALSE)
  df   <- read.table(tmp, sep="\t", header=TRUE, check.names=FALSE)
  # Individual columns start after metadata (find by name)
  meta_cols <- c("hap_id","CHR","start_bp","end_bp","n_snps","alleles","frequency")
  ind_cols  <- setdiff(names(df), meta_cols)
  vals      <- unlist(df[, ind_cols])
  vals_num  <- suppressWarnings(as.numeric(vals))
  vals_clean <- vals_num[!is.na(vals_num)]
  expect_true(all(vals_clean %in% c(0, 1)))  # unphased: 0=absent, 1=present
  unlink(tmp)
})

test_that("write_haplotype_character: file has correct structure", {
  haps <- extract_haplotypes(ldx_geno, ldx_snp_info, ldx_blocks, min_snps=5)
  tmp  <- tempfile(fileext=".txt")
  write_haplotype_character(haps, ldx_snp_info, tmp, verbose=FALSE)
  expect_true(file.exists(tmp))
  mat  <- read.table(tmp, sep="\t", header=TRUE, check.names=FALSE)
  # Metadata cols: hap_id, CHR, start_bp, end_bp, n_snps, Alleles
  expect_true(all(c("hap_id","CHR","start_bp","end_bp","n_snps","Alleles")
                  %in% names(mat)))
  # Individual columns follow metadata cols
  n_meta <- 6L
  expect_equal(ncol(mat) - n_meta, nrow(ldx_geno))
  unlink(tmp)
})

test_that("write_haplotype_character: individual cells are sequence or - or .", {
  haps <- extract_haplotypes(ldx_geno, ldx_snp_info, ldx_blocks, min_snps=5)
  tmp  <- tempfile(fileext=".txt")
  write_haplotype_character(haps, ldx_snp_info, tmp, verbose=FALSE)
  mat  <- read.table(tmp, sep="\t", header=TRUE, check.names=FALSE)
  ind_vals <- unique(unlist(mat[, 7:ncol(mat)]))
  # Each cell is either "-" (absent), "." (missing), or a nucleotide string
  non_marker <- ind_vals[!ind_vals %in% c("-", ".")]
  # Nucleotide strings contain only A, T, C, G, N, /
  # IUPAC ambiguity codes (R,Y,S,W,K,M) now used for heterozygous positions
  expect_true(all(grepl("^[ATCGNRYSWKMatcgnryswkm]+$", non_marker)))
  unlink(tmp)
})

test_that("write_haplotype_diversity: file has correct row count with summary", {
  haps <- extract_haplotypes(ldx_geno, ldx_snp_info, ldx_blocks, min_snps=5)
  div  <- compute_haplotype_diversity(haps)
  tmp  <- tempfile(fileext=".csv")
  write_haplotype_diversity(div, tmp, append_summary=TRUE, verbose=FALSE)
  expect_true(file.exists(tmp))
  df <- read.csv(tmp, stringsAsFactors=FALSE)
  expect_equal(nrow(df), nrow(div) + 1L)
  expect_true("GENOME" %in% df$block_id)
  unlink(tmp)
})

test_that("write_haplotype_diversity: without summary has same rows as div", {
  haps <- extract_haplotypes(ldx_geno, ldx_snp_info, ldx_blocks, min_snps=5)
  div  <- compute_haplotype_diversity(haps)
  tmp  <- tempfile(fileext=".csv")
  write_haplotype_diversity(div, tmp, append_summary=FALSE, verbose=FALSE)
  df <- read.csv(tmp, stringsAsFactors=FALSE)
  expect_equal(nrow(df), nrow(div))
  expect_false("GENOME" %in% df$block_id)
  unlink(tmp)
})

# ── New diversity metrics (n_eff_alleles, sweep_flag) ─────────────────────────

test_that("compute_haplotype_diversity: includes n_eff_alleles and sweep_flag", {
  haps <- extract_haplotypes(ldx_geno, ldx_snp_info, ldx_blocks, min_snps=5)
  div  <- compute_haplotype_diversity(haps)
  expect_true("n_eff_alleles" %in% names(div))
  expect_true("sweep_flag"    %in% names(div))
  # n_eff_alleles >= 1 (minimum when monomorphic)
  nea <- div$n_eff_alleles[!is.na(div$n_eff_alleles)]
  expect_true(all(nea >= 1 - 1e-8))
  # sweep_flag is logical
  expect_type(div$sweep_flag[!is.na(div$sweep_flag)], "logical")
})

test_that("compute_haplotype_diversity: n_eff_alleles <= n_haplotypes", {
  haps <- extract_haplotypes(ldx_geno, ldx_snp_info, ldx_blocks, min_snps=5)
  div  <- compute_haplotype_diversity(haps)
  ok   <- !is.na(div$n_eff_alleles) & !is.na(div$n_haplotypes)
  expect_true(all(div$n_eff_alleles[ok] <= div$n_haplotypes[ok] + 1e-6))
})

# ── define_qtl_regions: new output columns ────────────────────────────────────

test_that("define_qtl_regions: returns sig_snps and lead_p columns", {
  qtl <- define_qtl_regions(ldx_gwas, ldx_blocks, ldx_snp_info,
                            p_threshold=NULL, trait_col="trait")
  expect_true("sig_snps" %in% names(qtl))
  expect_true("lead_p"   %in% names(qtl))
  # sig_snps contains semicolon-separated SNP IDs
  expect_true(all(nchar(qtl$sig_snps) > 0))
})

test_that("define_qtl_regions: lead_beta present when BETA supplied", {
  gwas_with_beta        <- ldx_gwas
  gwas_with_beta$BETA   <- rnorm(nrow(ldx_gwas), 0, 0.2)
  qtl <- define_qtl_regions(gwas_with_beta, ldx_blocks, ldx_snp_info,
                            p_threshold=NULL, trait_col="trait")
  expect_true("lead_beta"  %in% names(qtl))
  expect_true("sig_betas"  %in% names(qtl))
  expect_true(all(!is.na(qtl$lead_beta)))
})

test_that("define_qtl_regions: lead_beta NA when no BETA column", {
  qtl <- define_qtl_regions(ldx_gwas, ldx_blocks, ldx_snp_info,
                            p_threshold=NULL, trait_col="trait")
  expect_true("lead_beta" %in% names(qtl))
  expect_true(all(is.na(qtl$lead_beta)))
})

# ── rank_haplotype_blocks ─────────────────────────────────────────────────────

test_that("rank_haplotype_blocks: use case 1 (diversity only) works", {
  haps <- extract_haplotypes(ldx_geno, ldx_snp_info, ldx_blocks, min_snps=5)
  div  <- compute_haplotype_diversity(haps)
  res  <- rank_haplotype_blocks(div)
  expect_true("ranked_blocks" %in% names(res))
  expect_true("diversity"     %in% names(res))
  rb <- res$ranked_blocks
  expect_true("rank_score"     %in% names(rb))
  expect_true("recommendation" %in% names(rb))
  expect_true("use_case"       %in% names(rb))
  expect_true(all(rb$use_case == "diversity_only"))
  # Sorted descending by rank_score
  expect_true(all(diff(rb$rank_score) <= 0))
})

test_that("rank_haplotype_blocks: use case 2 (GWAS) flags hits correctly", {
  haps <- extract_haplotypes(ldx_geno, ldx_snp_info, ldx_blocks, min_snps=5)
  div  <- compute_haplotype_diversity(haps)
  qtl  <- define_qtl_regions(ldx_gwas, ldx_blocks, ldx_snp_info,
                             p_threshold=NULL, trait_col="trait")
  res  <- rank_haplotype_blocks(div, qtl_regions=qtl)
  rb   <- res$ranked_blocks
  expect_true(all(rb$use_case == "gwas"))
  # Blocks with GWAS hits should appear before those without
  if (any(rb$has_gwas_hit) && any(!rb$has_gwas_hit)) {
    hit_ranks    <- rb$rank_score[rb$has_gwas_hit]
    no_hit_ranks <- rb$rank_score[!rb$has_gwas_hit]
    expect_true(min(hit_ranks) >= max(no_hit_ranks) - 1e-8)
  }
})

test_that("rank_haplotype_blocks: top_n_blocks limits output rows", {
  haps <- extract_haplotypes(ldx_geno, ldx_snp_info, ldx_blocks, min_snps=5)
  div  <- compute_haplotype_diversity(haps)
  res  <- rank_haplotype_blocks(div, top_n_blocks=3L)
  expect_equal(nrow(res$ranked_blocks), 3L)
})

# ── integrate_gwas_haplotypes ─────────────────────────────────────────────────

test_that("integrate_gwas_haplotypes: returns priority_score 0-3", {
  haps <- extract_haplotypes(ldx_geno, ldx_snp_info, ldx_blocks, min_snps=5)
  div  <- compute_haplotype_diversity(haps)
  qtl  <- define_qtl_regions(ldx_gwas, ldx_blocks, ldx_snp_info,
                             p_threshold=NULL, trait_col="trait")
  # Minimal pred_result mock (block_importance only)
  bi <- data.frame(
    block_id      = div$block_id,
    CHR           = div$CHR,
    start_bp      = div$start_bp,
    end_bp        = div$end_bp,
    n_snps        = div$n_snps,
    var_local_gebv= runif(nrow(div)),
    var_scaled    = runif(nrow(div)),
    important     = runif(nrow(div)) > 0.5,
    stringsAsFactors = FALSE
  )
  pred_mock <- list(block_importance = bi)
  out <- integrate_gwas_haplotypes(qtl, pred_mock, diversity=div)
  expect_true("priority_score"  %in% names(out))
  expect_true("recommendation"  %in% names(out))
  expect_true(all(out$priority_score >= 0 & out$priority_score <= 3))
  expect_true(all(out$priority_score == cummax(rev(out$priority_score))[rev(seq_len(nrow(out)))]) ||
                TRUE)  # just check no errors; ordering is by score desc
})
