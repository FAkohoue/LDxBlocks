## tests/testthat/test-basic.R
## Fast smoke tests - package loads, example data accessible, high-level
## functions produce output of the correct type. Target: <10 s on CRAN/CI.

library(testthat)
library(LDxBlocks)

# -- Example data --------------------------------------------------------------
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
                    "start.bp","end.bp","CHR","n_snps") %in% names(ldx_blocks)))
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

# -- C++ kernels: basic contracts ----------------------------------------------
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

# -- read_geno: matrix backend -------------------------------------------------
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

# -- flat-file formats ---------------------------------------------------------
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

# -- Big_LD + run_Big_LD_all_chr -----------------------------------------------
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
  expect_true(all(c("CHR","start","end","start.bp","end.bp","length_bp","n_snps")
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

# -- summarise + prepare_geno --------------------------------------------------
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

# -- Haplotype pipeline smoke tests --------------------------------------------
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
  feat <- build_haplotype_feature_matrix(haps, top_n=3, encoding="additive_012")$matrix
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
  feat <- build_haplotype_feature_matrix(haps, top_n=3, encoding="presence_01")$matrix
  vals <- as.vector(feat[!is.na(feat)])
  expect_true(all(vals %in% c(0, 1)))  # presence/absence: 0 or 1
})

# -- define_qtl_regions --------------------------------------------------------
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

# -- output writers ------------------------------------------------------------
test_that("write_haplotype_numeric writes readable dosage table", {
  data(ldx_geno,     package = "LDxBlocks")
  data(ldx_snp_info, package = "LDxBlocks")
  data(ldx_blocks,   package = "LDxBlocks")
  haps <- extract_haplotypes(ldx_geno, ldx_snp_info, ldx_blocks, min_snps=5)
  feat <- build_haplotype_feature_matrix(haps, top_n=2)$matrix
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

# -- tune_LD_params ------------------------------------------------------------
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

# -- ldx_blues dataset ---------------------------------------------------------

test_that("ldx_blues loads and has correct structure", {
  data(ldx_blues, package = "LDxBlocks")
  expect_s3_class(ldx_blues, "data.frame")
  expect_equal(nrow(ldx_blues), 120L)
  expect_true(all(c("id", "YLD", "RES") %in% names(ldx_blues)))
  expect_type(ldx_blues$YLD, "double")
  expect_type(ldx_blues$RES, "double")
  expect_false(any(is.na(ldx_blues$YLD)))
  expect_false(any(is.na(ldx_blues$RES)))
})

test_that("ldx_blues IDs match rownames of ldx_geno", {
  data(ldx_blues, package = "LDxBlocks")
  data(ldx_geno,  package = "LDxBlocks")
  expect_true(all(ldx_blues$id %in% rownames(ldx_geno)))
  expect_equal(length(unique(ldx_blues$id)), 120L)
})

test_that("example_blues.csv extdata file exists and is readable", {
  f <- system.file("extdata", "example_blues.csv", package = "LDxBlocks")
  skip_if(!file.exists(f), "example_blues.csv not found in extdata")
  blues <- read.csv(f, stringsAsFactors = FALSE)
  expect_true(all(c("id", "YLD", "RES") %in% names(blues)))
  expect_equal(nrow(blues), 120L)
  expect_type(blues$YLD, "double")
})

# -- ldx_blues_list dataset ----------------------------------------------------

test_that("ldx_blues_list loads and has correct structure", {
  data(ldx_blues_list, package = "LDxBlocks")
  expect_type(ldx_blues_list, "list")
  expect_equal(names(ldx_blues_list), c("env1", "env2"))
  expect_equal(length(ldx_blues_list$env1), 120L)
  expect_equal(length(ldx_blues_list$env2), 120L)
  expect_type(ldx_blues_list$env1, "double")
  expect_type(ldx_blues_list$env2, "double")
})

test_that("ldx_blues_list IDs match ldx_geno rownames", {
  data(ldx_blues_list, package = "LDxBlocks")
  data(ldx_geno,       package = "LDxBlocks")
  expect_equal(names(ldx_blues_list$env1), rownames(ldx_geno))
  expect_equal(names(ldx_blues_list$env2), rownames(ldx_geno))
})

test_that("ldx_blues_list: env2 has higher mean than env1 (simulated offset)", {
  data(ldx_blues_list, package = "LDxBlocks")
  # env2 was generated with +0.4 mean shift above env1
  diff_means <- mean(ldx_blues_list$env2) - mean(ldx_blues_list$env1)
  expect_true(diff_means > 0.1,
              label = paste("mean(env2) - mean(env1) =", round(diff_means, 3)))
})

test_that("example_blues_env.csv extdata file exists and is readable", {
  f <- system.file("extdata", "example_blues_env.csv", package = "LDxBlocks")
  skip_if(!file.exists(f), "example_blues_env.csv not found in extdata")
  blues_env <- read.csv(f, stringsAsFactors = FALSE)
  expect_true(all(c("id", "YLD_env1", "YLD_env2") %in% names(blues_env)))
  expect_equal(nrow(blues_env), 120L)
})

# -- run_ldx_pipeline ----------------------------------------------------------
test_that("run_ldx_pipeline: Path A (file path) returns correct result list", {
  f <- system.file("extdata", "example_genotypes_numeric.csv",
                   package = "LDxBlocks")
  skip_if(!file.exists(f), "extdata CSV not available")
  res <- run_ldx_pipeline(
    geno_source    = f,
    out_blocks     = tempfile(fileext = ".csv"),
    out_diversity  = tempfile(fileext = ".csv"),
    out_hap_matrix = tempfile(fileext = ".csv"),
    CLQcut         = 0.5, leng = 10L, subSegmSize = 80L,
    min_snps_block = 3L, verbose = FALSE
  )
  expect_true(all(c("blocks","diversity","hap_matrix","haplotypes",
                    "geno_matrix","snp_info_filtered") %in% names(res)))
  expect_s3_class(res$blocks, "data.frame")
  expect_true(nrow(res$blocks) >= 1L)
  expect_true(all(c("He","Shannon","sweep_flag") %in% names(res$diversity)))
  expect_equal(nrow(res$hap_matrix), 120L)
})

test_that("run_ldx_pipeline: Path A accepts LDxBlocks_backend as geno_source", {
  data(ldx_geno,     package = "LDxBlocks")
  data(ldx_snp_info, package = "LDxBlocks")
  be <- read_geno(ldx_geno, format = "matrix", snp_info = ldx_snp_info)
  on.exit(close_backend(be))
  res <- run_ldx_pipeline(
    geno_source    = be,
    out_blocks     = tempfile(fileext = ".csv"),
    out_diversity  = tempfile(fileext = ".csv"),
    out_hap_matrix = tempfile(fileext = ".csv"),
    CLQcut         = 0.5, leng = 10L, subSegmSize = 70L,
    min_snps_block = 3L, verbose = FALSE
  )
  expect_s3_class(res$blocks, "data.frame")
  expect_true(nrow(res$blocks) >= 1L)
})

test_that("run_ldx_pipeline: use_bigmemory=TRUE builds bigmemory backend internally", {
  skip_if_not_installed("bigmemory")
  f <- system.file("extdata", "example_genotypes_numeric.csv",
                   package = "LDxBlocks")
  skip_if(!file.exists(f), "extdata CSV not available")
  # Use a unique directory per test - avoids cross-test contamination
  # from shared tempdir() and Windows memory-mapped file locks.
  bm_path <- file.path(tempdir(), paste0("ldxbm_test_", as.integer(Sys.time())))
  dir.create(bm_path, recursive = TRUE, showWarnings = FALSE)
  # Ensure directory was actually created (path may be invalid on some systems)
  skip_if(!dir.exists(bm_path), "Could not create bigmemory temp directory")
  on.exit(unlink(bm_path, recursive = TRUE), add = TRUE)

  res <- run_ldx_pipeline(
    geno_source    = f,
    out_blocks     = tempfile(fileext = ".csv"),
    out_diversity  = tempfile(fileext = ".csv"),
    out_hap_matrix = tempfile(fileext = ".csv"),
    CLQcut         = 0.5, leng = 10L, subSegmSize = 80L,
    min_snps_block = 3L, verbose = FALSE,
    use_bigmemory  = TRUE,
    bigmemory_path = bm_path,
    bigmemory_type = "char"
  )
  expect_s3_class(res$blocks, "data.frame")
  expect_true(nrow(res$blocks) >= 1L)
  # Check that backing files were created. On Windows, bigmemory uses
  # forward slashes internally so we check both separator variants.
  # We also use list.files() as a robust alternative to file.exists().
  bm_path_norm <- normalizePath(bm_path, mustWork = FALSE)
  bm_path_fwd  <- gsub("\\\\", "/", bm_path_norm, fixed = TRUE)
  created_files <- list.files(bm_path_norm, recursive = FALSE)
  # bigmemory >= 1.4.7 names the file "ldxblocks_bm.bin";
  # older versions omit the .bin extension ("ldxblocks_bm").
  # Accept either form.
  bin_exists  <- any(c("ldxblocks_bm.bin", "ldxblocks_bm") %in% created_files) ||
    file.exists(file.path(bm_path_fwd, "ldxblocks_bm.bin")) ||
    file.exists(file.path(bm_path_fwd, "ldxblocks_bm"))
  desc_exists <- "ldxblocks_bm.desc"      %in% created_files ||
    file.exists(file.path(bm_path_fwd, "ldxblocks_bm.desc"))
  si_exists   <- "ldxblocks_bm_snpinfo.rds" %in% created_files ||
    file.exists(file.path(bm_path_fwd, "ldxblocks_bm_snpinfo.rds"))
  expect_true(bin_exists,
              label = paste0("ldxblocks_bm(.bin) in '", bm_path_norm, "' (found: ",
                             paste(created_files, collapse=", "), ")"))
  expect_true(desc_exists, label = "ldxblocks_bm.desc should exist")
  expect_true(si_exists,   label = "ldxblocks_bm_snpinfo.rds should exist")
})

test_that("run_ldx_pipeline: use_bigmemory reattaches on second call", {
  skip_if_not_installed("bigmemory")
  f <- system.file("extdata", "example_genotypes_numeric.csv",
                   package = "LDxBlocks")
  skip_if(!file.exists(f), "extdata CSV not available")
  bm_path <- file.path(tempdir(), paste0("ldxbm_reattach_", as.integer(Sys.time())))
  dir.create(bm_path, showWarnings = FALSE)
  on.exit(unlink(bm_path, recursive = TRUE), add = TRUE)

  # First run: builds the backing files
  run_ldx_pipeline(
    geno_source    = f,
    out_blocks     = tempfile(fileext = ".csv"),
    out_diversity  = tempfile(fileext = ".csv"),
    out_hap_matrix = tempfile(fileext = ".csv"),
    CLQcut = 0.5, leng = 10L, subSegmSize = 80L,
    min_snps_block = 3L, verbose = FALSE,
    use_bigmemory = TRUE, bigmemory_path = bm_path
  )
  # Normalize path to match what pipeline.R writes to (forward slashes).
  # Accept both ldxblocks_bm.bin (bigmemory >= 1.4.7) and ldxblocks_bm (older).
  bm_path_r <- normalizePath(bm_path, mustWork = FALSE)
  bm_path_r <- gsub("\\\\", "/", bm_path_r, fixed = TRUE)
  created_r <- list.files(bm_path_r, recursive = FALSE)
  skip_if(!any(c("ldxblocks_bm.bin", "ldxblocks_bm") %in% created_r),
          "First run did not create backing files -- skipping reattach test")
  # Second run: should reattach without rebuilding
  expect_no_error(
    run_ldx_pipeline(
      geno_source    = f,
      out_blocks     = tempfile(fileext = ".csv"),
      out_diversity  = tempfile(fileext = ".csv"),
      out_hap_matrix = tempfile(fileext = ".csv"),
      CLQcut = 0.5, leng = 10L, subSegmSize = 80L,
      min_snps_block = 3L, verbose = FALSE,
      use_bigmemory = TRUE, bigmemory_path = bm_path
    )
  )
})

test_that("run_ldx_pipeline: stale imputed cache rebuilt when maf_cut changes", {
  skip_if_not_installed("bigmemory")
  f <- system.file("extdata", "example_genotypes_numeric.csv",
                   package = "LDxBlocks")
  skip_if(!file.exists(f), "extdata CSV not available")
  bm_path <- file.path(tempdir(), paste0("ldxbm_fp_", as.integer(Sys.time())))
  dir.create(bm_path, showWarnings = FALSE)
  on.exit(unlink(bm_path, recursive = TRUE), add = TRUE)

  # First run with maf_cut = 0.05 -> builds imputed backend + fingerprint
  run_ldx_pipeline(
    geno_source    = f,
    out_blocks     = tempfile(fileext = ".csv"),
    out_diversity  = tempfile(fileext = ".csv"),
    out_hap_matrix = tempfile(fileext = ".csv"),
    CLQcut = 0.5, leng = 10L, subSegmSize = 80L,
    min_snps_block = 3L, verbose = FALSE,
    use_bigmemory = TRUE, bigmemory_path = bm_path,
    maf_cut = 0.05
  )
  fp_file <- file.path(bm_path, "ldxblocks_bm_imputed_params.rds")
  # Normalise to forward slashes (Windows bigmemory behaviour)
  fp_file <- gsub("\\", "/", normalizePath(fp_file, mustWork = FALSE), fixed = TRUE)
  # Check fingerprint was written
  expect_true(file.exists(fp_file) ||
                file.exists(file.path(bm_path, "ldxblocks_bm_imputed_params.rds")),
              label = "fingerprint file must exist after first run")

  # Second run with different maf_cut -> pipeline must detect stale cache and rebuild
  # No error is the key assertion; if stale detection is broken, bigmemory raises
  # "Backing file already exists" instead.
  expect_no_error(
    run_ldx_pipeline(
      geno_source    = f,
      out_blocks     = tempfile(fileext = ".csv"),
      out_diversity  = tempfile(fileext = ".csv"),
      out_hap_matrix = tempfile(fileext = ".csv"),
      CLQcut = 0.5, leng = 10L, subSegmSize = 80L,
      min_snps_block = 3L, verbose = FALSE,
      use_bigmemory = TRUE, bigmemory_path = bm_path,
      maf_cut = 0.10   # different -> stale -> rebuild
    )
  )
})


test_that("run_ldx_pipeline: partial backing files cleaned automatically", {
  skip_if_not_installed("bigmemory")
  f <- system.file("extdata", "example_genotypes_numeric.csv",
                   package = "LDxBlocks")
  skip_if(!file.exists(f), "extdata CSV not available")
  bm_path <- file.path(tempdir(), paste0("ldxbm_partial_", as.integer(Sys.time())))
  dir.create(bm_path, showWarnings = FALSE)
  on.exit(unlink(bm_path, recursive = TRUE), add = TRUE)

  # Write only the .bin file (simulate interrupted previous run)
  # Use a small text file - bigmemory will remove and replace it
  writeLines("stale content", file.path(bm_path, "ldxblocks_bm.bin"))
  expect_true(file.exists(file.path(bm_path, "ldxblocks_bm.bin")))
  expect_false(file.exists(file.path(bm_path, "ldxblocks_bm.desc")))

  # Pipeline should detect partial state, clean it, and rebuild without error
  expect_no_error(
    run_ldx_pipeline(
      geno_source    = f,
      out_blocks     = tempfile(fileext = ".csv"),
      out_diversity  = tempfile(fileext = ".csv"),
      out_hap_matrix = tempfile(fileext = ".csv"),
      CLQcut = 0.5, leng = 10L, subSegmSize = 80L,
      min_snps_block = 3L, verbose = FALSE,
      use_bigmemory = TRUE, bigmemory_path = bm_path
    )
  )
})
