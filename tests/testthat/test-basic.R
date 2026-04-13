## tests/testthat/test-basic.R
## Fast smoke tests — package loads, example data accessible, high-level
## functions produce output of the correct type. Target: <10 s on CRAN/CI.

library(testthat)
library(LDxBlocks)

# ── Example data ──────────────────────────────────────────────────────────────
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
  expect_true(all(c("start","end","start.rsID","end.rsID",
                    "start.bp","end.bp","CHR") %in% names(ldx_blocks)))
  expect_equal(nrow(ldx_gwas), 20L)
  # ldx_gwas has 5 columns: Marker, CHR, POS, P, trait
  expect_true(all(c("Marker","CHR","POS","P","trait") %in% names(ldx_gwas)))
  expect_true(all(ldx_gwas$trait %in% c("TraitA","TraitB")))
})

test_that("example SNP IDs match between geno matrix and snp_info", {
  data(ldx_geno,     package = "LDxBlocks")
  data(ldx_snp_info, package = "LDxBlocks")
  expect_equal(colnames(ldx_geno), ldx_snp_info$SNP)
})

test_that("example chromosomes are normalised (no chr prefix)", {
  data(ldx_snp_info, package = "LDxBlocks")
  expect_true(all(ldx_snp_info$CHR %in% c("1","2","3")))
})

# ── C++ kernels: basic contracts ──────────────────────────────────────────────
test_that("compute_r2_cpp: symmetric, [0,1], zero diagonal", {
  data(ldx_geno, package = "LDxBlocks")
  r2 <- compute_r2_cpp(ldx_geno[, 1:20], digits = -1L, n_threads = 1L)
  expect_equal(dim(r2), c(20L, 20L))
  expect_true(isSymmetric(r2, tol = 1e-10))
  expect_true(all(r2 >= 0 & r2 <= 1 + 1e-8))
  expect_equal(diag(r2), rep(0, 20))
})

test_that("maf_filter_cpp: removes monomorphic and low-MAF columns", {
  data(ldx_geno, package = "LDxBlocks")
  G      <- ldx_geno[, 1:30]
  G[, 1] <- 0L
  G[, 2] <- c(rep(0L, 115), rep(1L, 5))
  keep   <- maf_filter_cpp(G, maf_cut = 0.05)
  expect_false(keep[1]); expect_false(keep[2])
  expect_true(all(keep[3:30]))
})

test_that("build_adj_matrix_cpp: valid 0/1 symmetric matrix", {
  data(ldx_geno, package = "LDxBlocks")
  r2  <- compute_r2_cpp(ldx_geno[, 1:15])
  adj <- build_adj_matrix_cpp(r2, 0.5)
  expect_true(all(adj == 0L | adj == 1L))
  expect_equal(diag(adj), rep(0L, 15))
  expect_true(isSymmetric(adj))
})

# ── read_geno: matrix backend ─────────────────────────────────────────────────
test_that("read_geno wraps matrix into valid LDxBlocks_backend", {
  data(ldx_geno,     package = "LDxBlocks")
  data(ldx_snp_info, package = "LDxBlocks")
  be <- read_geno(ldx_geno, format = "matrix", snp_info = ldx_snp_info)
  expect_s3_class(be, "LDxBlocks_backend")
  expect_equal(be$type,      "matrix")
  expect_equal(be$n_samples, 120L)
  expect_equal(be$n_snps,    230L)
  expect_equal(be$sample_ids, rownames(ldx_geno))
  close_backend(be)
})

test_that("read_chunk returns correct slice dimensions", {
  data(ldx_geno,     package = "LDxBlocks")
  data(ldx_snp_info, package = "LDxBlocks")
  be    <- read_geno(ldx_geno, format = "matrix", snp_info = ldx_snp_info)
  chunk <- read_chunk(be, 10:25)
  expect_equal(dim(chunk), c(120L, 16L))
  close_backend(be)
})

# ── flat-file formats ─────────────────────────────────────────────────────────
test_that("read_geno reads numeric CSV from extdata", {
  f <- system.file("extdata","example_genotypes_numeric.csv",package="LDxBlocks")
  skip_if(!file.exists(f), "extdata CSV not available")
  be <- read_geno(f)
  expect_equal(be$type, "numeric"); expect_equal(be$n_snps, 230L)
  close_backend(be)
})

test_that("read_geno reads HapMap from extdata", {
  f <- system.file("extdata","example_genotypes.hmp.txt",package="LDxBlocks")
  skip_if(!file.exists(f), "extdata HapMap not available")
  be <- read_geno(f)
  # SNPRelate is used for GDS conversion so type is always "gds" when
  # SNPRelate is installed. Falls back to "hapmap" if SNPRelate is absent.
  expect_true(be$type %in% c("hapmap", "gds")); expect_equal(be$n_snps, 230L)
  close_backend(be)
})

test_that("read_geno reads VCF from extdata", {
  f <- system.file("extdata","example_genotypes.vcf",package="LDxBlocks")
  skip_if(!file.exists(f), "extdata VCF not available")
  be <- read_geno(f)
  # SNPRelate converts VCF to GDS so type is "gds" when SNPRelate is
  # installed; falls back to "vcf" if absent. Both are valid.
  expect_true(be$type %in% c("vcf", "gds"))
  expect_true(all(be$snp_info$CHR %in% c("1","2","3")))
  close_backend(be)
})

test_that("read_geno reads example_gwas.csv from extdata", {
  f <- system.file("extdata","example_gwas.csv",package="LDxBlocks")
  skip_if(!file.exists(f), "extdata GWAS not available")
  gwas <- read.csv(f, stringsAsFactors = FALSE)
  expect_true(all(c("Marker","CHR","POS","P","trait") %in% names(gwas)))
  expect_true(all(gwas$trait %in% c("TraitA","TraitB")))
})

test_that("read_geno reads example_phenotype.csv from extdata", {
  f <- system.file("extdata","example_phenotype.csv",package="LDxBlocks")
  skip_if(!file.exists(f), "extdata phenotype not available")
  pheno <- read.csv(f, stringsAsFactors = FALSE)
  expect_true(all(c("Sample","Trait1","Trait2","PC1","PC2") %in% names(pheno)))
  expect_equal(nrow(pheno), 120L)
})

# ── Big_LD + run_Big_LD_all_chr ───────────────────────────────────────────────
test_that("Big_LD internal: returns valid block data.frame for chr1", {
  data(ldx_geno,     package = "LDxBlocks")
  data(ldx_snp_info, package = "LDxBlocks")
  chr1_idx <- which(ldx_snp_info$CHR == "1")
  blocks   <- LDxBlocks:::Big_LD(ldx_geno[, chr1_idx],
                                 ldx_snp_info[chr1_idx, c("SNP","POS")],
                                 method="r2", CLQcut=0.5, leng=8,
                                 subSegmSize=80, verbose=FALSE)
  expect_s3_class(blocks, "data.frame")
  expect_true(nrow(blocks) >= 1L)
  expect_true(all(blocks$start <= blocks$end))
  expect_true(all(blocks$start.bp <= blocks$end.bp))
})

test_that("run_Big_LD_all_chr returns blocks for all 3 chromosomes", {
  data(ldx_geno,     package = "LDxBlocks")
  data(ldx_snp_info, package = "LDxBlocks")
  blocks <- run_Big_LD_all_chr(ldx_geno, snp_info=ldx_snp_info,
                               method="r2", CLQcut=0.4,
                               leng=5, subSegmSize=100, verbose=FALSE)
  expect_true(all(c("1","2","3") %in% blocks$CHR))
  expect_true(all(c("CHR","start","end","start.bp","end.bp","length_bp")
                  %in% names(blocks)))
})

test_that("run_Big_LD_all_chr accepts LDxBlocks_backend", {
  data(ldx_geno,     package = "LDxBlocks")
  data(ldx_snp_info, package = "LDxBlocks")
  be     <- read_geno(ldx_geno, format="matrix", snp_info=ldx_snp_info)
  blocks <- run_Big_LD_all_chr(be, method="r2", CLQcut=0.5,
                               leng=10, subSegmSize=70, verbose=FALSE)
  expect_s3_class(blocks, "data.frame")
  expect_true(nrow(blocks) >= 3L)
  close_backend(be)
})

# ── summarise + prepare_geno ──────────────────────────────────────────────────
test_that("summarise_blocks returns per-chromosome and GENOME rows", {
  data(ldx_blocks, package = "LDxBlocks")
  s <- summarise_blocks(ldx_blocks)
  expect_true("GENOME" %in% s$CHR)
  expect_equal(nrow(s), 4L)
})

test_that("prepare_geno r2 path returns centred matrix", {
  data(ldx_geno, package = "LDxBlocks")
  prep <- prepare_geno(ldx_geno[, 1:20], method = "r2")
  expect_null(prep$V_inv_sqrt)
  expect_true(max(abs(colMeans(prep$adj_geno))) < 1e-10)
})

# ── Haplotype pipeline smoke tests ────────────────────────────────────────────
test_that("extract_haplotypes returns one entry per block", {
  data(ldx_geno,     package = "LDxBlocks")
  data(ldx_snp_info, package = "LDxBlocks")
  data(ldx_blocks,   package = "LDxBlocks")
  haps <- extract_haplotypes(ldx_geno, ldx_snp_info, ldx_blocks, min_snps=3)
  expect_type(haps, "list"); expect_true(length(haps) >= 1L)
  expect_equal(length(haps[[1]]), 120L)
})

test_that("compute_haplotype_diversity returns valid metrics", {
  data(ldx_geno,     package = "LDxBlocks")
  data(ldx_snp_info, package = "LDxBlocks")
  data(ldx_blocks,   package = "LDxBlocks")
  haps <- extract_haplotypes(ldx_geno, ldx_snp_info, ldx_blocks, min_snps=3)
  div  <- compute_haplotype_diversity(haps)
  expect_true(all(c("block_id","He","Shannon","n_eff_alleles",
                    "freq_dominant","sweep_flag") %in% names(div)))
  expect_true(all(div$He[!is.na(div$He)] >= 0))
})

test_that("build_haplotype_feature_matrix: correct dimensions, additive_012", {
  data(ldx_geno,     package = "LDxBlocks")
  data(ldx_snp_info, package = "LDxBlocks")
  data(ldx_blocks,   package = "LDxBlocks")
  haps <- extract_haplotypes(ldx_geno, ldx_snp_info, ldx_blocks, min_snps=3)
  feat <- build_haplotype_feature_matrix(haps, top_n=3, encoding="additive_012")
  expect_equal(nrow(feat), 120L)
  expect_equal(ncol(feat), length(haps) * 3L)
  vals <- as.vector(feat[!is.na(feat)])
  expect_true(all(vals %in% c(0, 1)))  # unphased: 0=absent, 1=present
})

test_that("build_haplotype_feature_matrix: presence_01 encoding", {
  data(ldx_geno,     package = "LDxBlocks")
  data(ldx_snp_info, package = "LDxBlocks")
  data(ldx_blocks,   package = "LDxBlocks")
  haps <- extract_haplotypes(ldx_geno, ldx_snp_info, ldx_blocks, min_snps=3)
  feat <- build_haplotype_feature_matrix(haps, top_n=3, encoding="presence_01")
  vals <- as.vector(feat[!is.na(feat)])
  expect_true(all(vals %in% c(0, 1)))  # presence/absence: 0 or 1
})

# ── define_qtl_regions ────────────────────────────────────────────────────────
test_that("define_qtl_regions: returns data.frame with pleiotropic column", {
  data(ldx_gwas,     package = "LDxBlocks")
  data(ldx_blocks,   package = "LDxBlocks")
  data(ldx_snp_info, package = "LDxBlocks")
  qtl <- define_qtl_regions(ldx_gwas, ldx_blocks, ldx_snp_info,
                            p_threshold = NULL, trait_col = "trait")
  expect_s3_class(qtl, "data.frame")
  expect_true(all(c("block_id","CHR","pleiotropic","n_traits","traits")
                  %in% names(qtl)))
  expect_true(is.logical(qtl$pleiotropic))
})

test_that("define_qtl_regions: works without trait column (single-trait)", {
  data(ldx_gwas,     package = "LDxBlocks")
  data(ldx_blocks,   package = "LDxBlocks")
  data(ldx_snp_info, package = "LDxBlocks")
  gwas_no_trait <- ldx_gwas[, c("Marker","CHR","POS","P")]
  qtl <- define_qtl_regions(gwas_no_trait, ldx_blocks, ldx_snp_info,
                            p_threshold = NULL)
  expect_true(all(qtl$n_traits == 1L))
  expect_true(all(!qtl$pleiotropic))
})

# ── output writers ────────────────────────────────────────────────────────────
test_that("write_haplotype_numeric writes readable dosage table", {
  data(ldx_geno,     package = "LDxBlocks")
  data(ldx_snp_info, package = "LDxBlocks")
  data(ldx_blocks,   package = "LDxBlocks")
  haps <- extract_haplotypes(ldx_geno, ldx_snp_info, ldx_blocks, min_snps=5)
  feat <- build_haplotype_feature_matrix(haps, top_n=2)
  tmp  <- tempfile(fileext=".csv")
  write_haplotype_numeric(feat, tmp, verbose=FALSE)
  df <- read.table(tmp, sep="\t", header=TRUE, check.names=FALSE)
  # Rows = haplotype alleles (not individuals)
  expect_equal(nrow(df), ncol(feat))
  expect_true("hap_id" %in% names(df))
  expect_false("Sample" %in% names(df))
  unlink(tmp)
})

test_that("write_haplotype_character writes readable nucleotide matrix", {
  data(ldx_geno,     package = "LDxBlocks")
  data(ldx_snp_info, package = "LDxBlocks")
  data(ldx_blocks,   package = "LDxBlocks")
  haps <- extract_haplotypes(ldx_geno, ldx_snp_info, ldx_blocks, min_snps=5)
  tmp  <- tempfile(fileext=".txt")
  write_haplotype_character(haps, ldx_snp_info, tmp, verbose=FALSE)
  mat <- read.table(tmp, sep="\t", header=TRUE, check.names=FALSE)
  expect_true("hap_id" %in% names(mat))
  expect_true("Alleles" %in% names(mat))
  expect_true(nrow(mat) >= 1L)
  unlink(tmp)
})

test_that("write_haplotype_diversity writes readable CSV with summary row", {
  data(ldx_geno,     package = "LDxBlocks")
  data(ldx_snp_info, package = "LDxBlocks")
  data(ldx_blocks,   package = "LDxBlocks")
  haps <- extract_haplotypes(ldx_geno, ldx_snp_info, ldx_blocks, min_snps=5)
  div  <- compute_haplotype_diversity(haps)
  tmp  <- tempfile(fileext=".csv")
  write_haplotype_diversity(div, tmp, append_summary=TRUE, verbose=FALSE)
  df <- read.csv(tmp, stringsAsFactors=FALSE)
  expect_true("GENOME" %in% df$block_id)
  expect_equal(nrow(df), nrow(div) + 1L)
  unlink(tmp)
})

# ── tune_LD_params ────────────────────────────────────────────────────────────
test_that("tune_LD_params returns correct result structure", {
  data(ldx_geno,     package = "LDxBlocks")
  data(ldx_snp_info, package = "LDxBlocks")
  data(ldx_gwas,     package = "LDxBlocks")
  grid <- expand.grid(
    CLQcut=c(0.5,0.6), clstgap=1e5, leng=10,
    subSegmSize=70, split=FALSE, checkLargest=FALSE,
    MAFcut=0.05, CLQmode="Density", kin_method="chol",
    digits=-1, appendrare=FALSE, stringsAsFactors=FALSE
  )
  res <- tune_LD_params(ldx_geno, ldx_snp_info, ldx_gwas,
                        grid=grid, prefer_perfect=TRUE, seed=1L)
  expect_true(all(c("best_params","score_table","final_blocks",
                    "gwas_assigned") %in% names(res)))
  expect_true("LD_block" %in% names(res$gwas_assigned))
})
