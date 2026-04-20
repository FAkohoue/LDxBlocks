## data-raw/generate_example_data.R
## -----------------------------------------------------------------------------
## Generates all example datasets shipped with LDxBlocks.
## Run once with: source("data-raw/generate_example_data.R")
## Outputs land in data/ (as .rda) and inst/extdata/ (as flat files).
##
## Datasets produced
## -----------------
##   ldx_geno      numeric matrix 120 ind x 230 SNPs  (three-block structure per chromosome)
##   ldx_snp_info  data.frame SNP / CHR / POS / REF / ALT
##   ldx_blocks    data.frame  - pre-computed block table for examples/tests
##   ldx_gwas      data.frame  - 20 toy GWAS markers for tune_LD_params demos
##   ldx_blues     data.frame  - 120 ind x 3 col BLUEs (id, YLD, RES) for run_haplotype_prediction demos
##   ldx_blues_list named list - 2 environments (env1, env2) of named numeric vectors for
##                              run_haplotype_stability() and cv_haplotype_prediction() demos
##
## Flat-file copies
## ----------------
##   inst/extdata/example_genotypes_numeric.csv   (numeric dosage)
##   inst/extdata/example_genotypes.hmp.txt       (HapMap)
##   inst/extdata/example_genotypes.vcf           (VCF)
##   inst/extdata/example_gwas.csv                (GWAS marker table)
##   inst/extdata/example_phenotype.csv           (simulated phenotypes)
##   inst/extdata/example_blues.csv               (pre-adjusted BLUEs, id/YLD/RES)
##   inst/extdata/example_blues_env.csv           (BLUEs by environment: id/YLD_env1/YLD_env2)
## -----------------------------------------------------------------------------

set.seed(20250407)

# -- Parameters ----------------------------------------------------------------
N_IND   <- 120      # individuals
N_CHR   <- 3        # chromosomes
BLOCK_SIZES <- list(
  chr1 = c(25, 20, 25),   # three blocks on chr1 (sums to 70)
  chr2 = c(30, 20, 20),   # three blocks on chr2
  chr3 = c(20, 20, 20)    # three blocks on chr3
)
INTER_BLOCK_SNP <- 5   # low-LD SNPs between blocks (singleton stretches)

ALLELES <- matrix(c(
  "A","T",  "G","C",  "C","G",  "T","A",  "A","G",
  "G","A",  "C","T",  "T","C",  "A","C",  "G","T"
), ncol = 2, byrow = TRUE)

# -- Helper: simulate one LD block using founder haplotypes -------------------
# Creates realistic block structure with K=4 distinct founder haplotypes.
# Each individual is a diploid combination of two founders.
# - MAF enforced >= 0.10 per SNP so all SNPs pass maf_filter_cpp(maf_cut=0.05).
# - Low flip rate (1%) preserves founder dosage strings so diplotype classes
#   (AA/AB/BB combinations) contain >= 3 individuals each, enabling
#   estimate_diplotype_effects() to run without data-driven skips.
make_block <- function(n_ind, n_snp, p = NULL, n_founders = 4L,
                       flip = 0.01, min_maf = 0.10) {
  # Retry until all SNPs have MAF >= min_maf
  for (attempt in seq_len(50L)) {
    # Step 1: K founder haplotypes (binary {0,1} vectors)
    founders <- matrix(0L, nrow = n_founders, ncol = n_snp)
    for (k in seq_len(n_founders)) {
      p_k <- runif(1, min_maf + 0.05, 1 - min_maf - 0.05)
      founders[k, ] <- rbinom(n_snp, 1L, p_k)
    }

    # Step 2: Sample founder pairs for each individual
    founder_freq <- runif(n_founders, 0.5, 1.5)
    founder_freq <- founder_freq / sum(founder_freq)
    hap1_idx <- sample(n_founders, n_ind, replace = TRUE, prob = founder_freq)
    hap2_idx <- sample(n_founders, n_ind, replace = TRUE, prob = founder_freq)

    # Step 3: Dosage = sum of the two haplotypes (values in {0, 1, 2})
    M <- founders[hap1_idx, , drop = FALSE] +
      founders[hap2_idx, , drop = FALSE]

    # Check MAF; retry if any SNP is near-monomorphic
    alt_freq <- colMeans(M) / 2
    maf      <- pmin(alt_freq, 1 - alt_freq)
    if (all(maf >= min_maf)) break
  }

  # Step 4: Very low flip rate - preserves founder haplotype identity
  # so the same founder pair produces the same dosage string across individuals,
  # creating repeatable diplotype classes.
  for (j in seq_len(n_snp)) {
    idx   <- sample.int(n_ind, max(1L, floor(flip * n_ind)))
    delta <- sample(c(-1L, 1L), length(idx), replace = TRUE)
    M[idx, j] <- pmin(2L, pmax(0L, M[idx, j] + delta))
  }
  M
}

# -- Helper: simulate low-LD singleton SNPs ------------------------------------
make_singleton <- function(n_ind, n_snp) {
  sapply(seq_len(n_snp), function(i) {
    p <- runif(1, 0.1, 0.4)
    rbinom(n_ind, 2, p)
  })
}

# -- Build genome --------------------------------------------------------------
geno_list    <- list()
snpinfo_list <- list()
global_pos   <- 1L

for (ci in seq_len(N_CHR)) {
  chr_label <- as.character(ci)
  bs        <- BLOCK_SIZES[[ci]]
  cur_pos   <- 1000L
  chr_geno  <- NULL
  snp_ids   <- character(0)
  snp_pos   <- integer(0)
  snp_ref   <- character(0)
  snp_alt   <- character(0)
  snp_count <- 0L

  for (bi in seq_along(bs)) {
    n_b <- bs[bi]
    blk <- make_block(N_IND, n_b)
    chr_geno <- cbind(chr_geno, blk)

    for (j in seq_len(n_b)) {
      snp_count <- snp_count + 1L
      al <- ALLELES[(snp_count %% nrow(ALLELES)) + 1L, ]
      snp_ids <- c(snp_ids, sprintf("rs%s%03d", chr_label, snp_count))
      snp_pos <- c(snp_pos, cur_pos)
      snp_ref <- c(snp_ref, al[1])
      snp_alt <- c(snp_alt, al[2])
      cur_pos <- cur_pos + sample(800L:1200L, 1L)
    }

    # Inter-block singletons (except after last block)
    if (bi < length(bs)) {
      sing <- make_singleton(N_IND, INTER_BLOCK_SNP)
      chr_geno <- cbind(chr_geno, sing)
      cur_pos <- cur_pos + 50000L   # genomic gap
      for (j in seq_len(INTER_BLOCK_SNP)) {
        snp_count <- snp_count + 1L
        al <- ALLELES[(snp_count %% nrow(ALLELES)) + 1L, ]
        snp_ids <- c(snp_ids, sprintf("rs%s%03d", chr_label, snp_count))
        snp_pos <- c(snp_pos, cur_pos)
        snp_ref <- c(snp_ref, al[1])
        snp_alt <- c(snp_alt, al[2])
        cur_pos <- cur_pos + sample(800L:1200L, 1L)
      }
    }
  }

  geno_list[[ci]] <- chr_geno
  snpinfo_list[[ci]] <- data.frame(
    SNP = snp_ids, CHR = chr_label,
    POS = snp_pos, REF = snp_ref, ALT = snp_alt,
    stringsAsFactors = FALSE
  )
}

ldx_geno     <- do.call(cbind, geno_list)
ldx_snp_info <- do.call(rbind, snpinfo_list)
rownames(ldx_geno) <- paste0("ind", sprintf("%03d", seq_len(N_IND)))
colnames(ldx_geno) <- ldx_snp_info$SNP

stopifnot(ncol(ldx_geno) == nrow(ldx_snp_info))
message("Simulated genotype matrix: ", nrow(ldx_geno), " ind x ", ncol(ldx_geno), " SNPs")

# -- Pre-compute block table (stored as example output) ------------------------
# Use in-memory path (no external deps at data-raw stage)
# We derive a simple block table by chromosome for illustration.
# Full run: blocks <- run_Big_LD_all_chr(ldx_geno, ldx_snp_info, CLQcut=0.6)
# We build a hand-crafted reference table matching the simulated structure.

make_block_row <- function(snp_info, snp_ids) {
  idx <- match(snp_ids, snp_info$SNP)
  data.frame(
    start      = min(idx), end = max(idx),
    start.rsID = snp_info$SNP[min(idx)],
    end.rsID   = snp_info$SNP[max(idx)],
    start.bp   = snp_info$POS[min(idx)],
    end.bp     = snp_info$POS[max(idx)],
    CHR        = snp_info$CHR[min(idx)],
    length_bp  = snp_info$POS[max(idx)] - snp_info$POS[min(idx)] + 1L,
    n_snps     = length(snp_ids),
    stringsAsFactors = FALSE
  )
}

# Derive expected block SNP ranges from simulation structure
chr1_snps <- ldx_snp_info$SNP[ldx_snp_info$CHR == "1"]
chr2_snps <- ldx_snp_info$SNP[ldx_snp_info$CHR == "2"]
chr3_snps <- ldx_snp_info$SNP[ldx_snp_info$CHR == "3"]

ldx_blocks <- rbind(
  # chr1: 3 blocks of sizes 25, 20, 25 with 5-SNP gaps between
  make_block_row(ldx_snp_info, chr1_snps[1:25]),
  make_block_row(ldx_snp_info, chr1_snps[31:50]),
  make_block_row(ldx_snp_info, chr1_snps[56:80]),
  # chr2: 3 blocks of sizes 30, 20, 20
  make_block_row(ldx_snp_info, chr2_snps[1:30]),
  make_block_row(ldx_snp_info, chr2_snps[36:55]),
  make_block_row(ldx_snp_info, chr2_snps[61:80]),
  # chr3: 3 blocks of sizes 20, 20, 20
  make_block_row(ldx_snp_info, chr3_snps[1:20]),
  make_block_row(ldx_snp_info, chr3_snps[26:45]),
  make_block_row(ldx_snp_info, chr3_snps[51:70])
)
rownames(ldx_blocks) <- NULL

message("Pre-computed block table: ", nrow(ldx_blocks), " blocks")

# -- Toy GWAS marker table for tune_LD_params demos ---------------------------
set.seed(42)
gwas_snp_idx <- c(
  sample(1:25, 3),           # from chr1 block1
  sample(31:50, 3),          # from chr1 block2
  sample(56:80, 3),          # from chr1 block3
  sample(which(ldx_snp_info$CHR=="2")[1:30], 4),
  sample(which(ldx_snp_info$CHR=="2")[36:55], 3),
  sample(which(ldx_snp_info$CHR=="3")[1:20], 4)
)
gwas_snp_idx <- sort(unique(gwas_snp_idx))[1:20]

ldx_gwas <- data.frame(
  Marker = ldx_snp_info$SNP[gwas_snp_idx],
  CHR    = ldx_snp_info$CHR[gwas_snp_idx],
  POS    = ldx_snp_info$POS[gwas_snp_idx],
  P      = c(runif(8, 1e-10, 9e-7),   # 8 markers clearly below 1e-6
             runif(7, 1e-5, 1e-3),   # 7 suggestive
             runif(5, 1e-3, 0.05)),  # 5 sub-threshold
  trait  = sample(c("TraitA", "TraitB"), 20, replace = TRUE, prob = c(0.6, 0.4)),
  stringsAsFactors = FALSE
)

message("GWAS marker table: ", nrow(ldx_gwas), " markers")



# -- Write flat-file copies to inst/extdata/ ----------------------------------

## 1. Numeric dosage CSV (SNPs x samples: SNP CHR POS REF ALT + sample cols)
numeric_df <- cbind(
  ldx_snp_info,
  as.data.frame(t(ldx_geno))
)
write.csv(numeric_df,
          file = "inst/extdata/example_genotypes_numeric.csv",
          row.names = FALSE, quote = FALSE
)

## 2. HapMap format
# Iterate over SNPs (rows of ldx_geno after transposing back to SNPs x samples)
# ldx_geno is individuals x SNPs, so columns = SNPs
hmp_decode_snp <- function(snp_dosages, ref, alt) {
  # snp_dosages: numeric vector length n_ind
  vapply(snp_dosages, function(x) {
    if (is.na(x)) return("NN")
    switch(as.character(as.integer(x)),
           "0" = paste0(ref, ref),
           "1" = paste0(ref, alt),
           "2" = paste0(alt, alt),
           "NN"
    )
  }, character(1))
}

# Build sample call matrix: n_snp rows x n_ind cols
n_snp <- nrow(ldx_snp_info)
n_ind <- nrow(ldx_geno)
hmp_calls <- matrix("NN", nrow = n_snp, ncol = n_ind)
for (si in seq_len(n_snp)) {
  hmp_calls[si, ] <- hmp_decode_snp(
    ldx_geno[, si],
    ldx_snp_info$REF[si],
    ldx_snp_info$ALT[si]
  )
}

hmp_header <- data.frame(
  "rs#"       = ldx_snp_info$SNP,
  alleles     = paste0(ldx_snp_info$REF, "/", ldx_snp_info$ALT),
  chrom       = ldx_snp_info$CHR,
  pos         = ldx_snp_info$POS,
  strand      = "+",
  "assembly#" = "NA",
  center      = "NA",
  protLSID    = "NA",
  assayLSID   = "NA",
  panelLSID   = "NA",
  QCcode      = "NA",
  check.names = FALSE,
  stringsAsFactors = FALSE
)
hmp_out <- cbind(hmp_header, as.data.frame(hmp_calls, stringsAsFactors = FALSE))
colnames(hmp_out)[12:ncol(hmp_out)] <- rownames(ldx_geno)
write.table(hmp_out,
            file = "inst/extdata/example_genotypes.hmp.txt",
            sep = "\t", row.names = FALSE, quote = FALSE
)

## 3. VCF
vcf_lines <- c(
  "##fileformat=VCFv4.2",
  paste0("##source=LDxBlocks example data, generated ", Sys.Date()),
  "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">",
  paste(c("#CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT",
          rownames(ldx_geno)), collapse = "\t")
)
gt_encode <- function(g) {
  vapply(g, function(x) {
    if (is.na(x)) return("./.")
    switch(as.character(x), "0"="0/0", "1"="0/1", "2"="1/1", "./.")
  }, character(1))
}
for (si in seq_len(nrow(ldx_snp_info))) {
  gt_vec <- gt_encode(ldx_geno[, si])
  vcf_lines <- c(vcf_lines, paste(c(
    ldx_snp_info$CHR[si], ldx_snp_info$POS[si], ldx_snp_info$SNP[si],
    ldx_snp_info$REF[si], ldx_snp_info$ALT[si], ".", "PASS", ".", "GT",
    gt_vec
  ), collapse = "\t"))
}
writeLines(vcf_lines, "inst/extdata/example_genotypes.vcf")

## 4. GWAS flat file (inst/extdata/example_gwas.csv)
write.csv(ldx_gwas,
          file = "inst/extdata/example_gwas.csv",
          row.names = FALSE, quote = FALSE
)

## 5. Simulated phenotype file (raw, unadjusted)
## Two quantitative traits (polygenic architecture from the LD blocks)
## plus two principal components (population structure surrogates).
## Trait1 / Trait2 are linear combinations of 10 random SNP effects plus noise.
## This file shows the raw phenotype format a user would start with before
## running a mixed model to obtain pre-adjusted means (BLUEs). The BLUEs
## for direct use with run_haplotype_prediction() are in example_blues.csv.
set.seed(20250408)
n_qtl <- 10L

# Random SNP effects for each trait
qtl_idx1  <- sample(ncol(ldx_geno), n_qtl)
qtl_eff1  <- rnorm(n_qtl, 0, 0.3)
trait1_val <- as.numeric(ldx_geno[, qtl_idx1] %*% qtl_eff1) + rnorm(N_IND, 0, 0.5)
trait1_val <- round((trait1_val - mean(trait1_val)) / sd(trait1_val), 4)

qtl_idx2  <- sample(ncol(ldx_geno), n_qtl)
qtl_eff2  <- rnorm(n_qtl, 0, 0.3)
trait2_val <- as.numeric(ldx_geno[, qtl_idx2] %*% qtl_eff2) + rnorm(N_IND, 0, 0.5)
trait2_val <- round((trait2_val - mean(trait2_val)) / sd(trait2_val), 4)

## -- Build ldx_blues: must come after trait1_val and trait2_val are defined --
## This data.frame is saved as .rda (usethis::use_data above) and as
## inst/extdata/example_blues.csv (section 5b below).
ldx_blues <- data.frame(
  id   = rownames(ldx_geno),
  YLD  = trait1_val,   # standardised yield-like BLUE
  RES  = trait2_val,   # standardised resistance-like BLUE
  stringsAsFactors = FALSE
)

# -- Save .rda to data/ --------------------------------------------------------
## -- ldx_blues_list: per-environment BLUEs for run_haplotype_stability() -----
## Two environments are simulated by adding small environment-specific offsets
## and noise to the YLD BLUEs. This mimics a multi-environment trial where
## the same genotypes were evaluated in two locations.
set.seed(20250409)
env1_blues <- setNames(
  trait1_val + rnorm(N_IND, mean = 0.0, sd = 0.15),
  rownames(ldx_geno)
)
env2_blues <- setNames(
  trait1_val + rnorm(N_IND, mean = 0.4, sd = 0.15),
  rownames(ldx_geno)
)
ldx_blues_list <- list(
  env1 = round(env1_blues, 4),
  env2 = round(env2_blues, 4)
)

usethis::use_data(ldx_geno,       overwrite = TRUE, compress = "xz")
usethis::use_data(ldx_snp_info,   overwrite = TRUE, compress = "xz")
usethis::use_data(ldx_blocks,     overwrite = TRUE, compress = "xz")
usethis::use_data(ldx_gwas,       overwrite = TRUE, compress = "xz")
usethis::use_data(ldx_blues,      overwrite = TRUE, compress = "xz")
usethis::use_data(ldx_blues_list, overwrite = TRUE, compress = "xz")

message("Saved all .rda files to data/")  # 6 datasets: geno, snp_info, blocks, gwas, blues, blues_list

# PCA on genotype matrix for population structure covariates
Gc_pca  <- scale(ldx_geno, center = TRUE, scale = FALSE)
pca_out <- prcomp(Gc_pca, rank. = 2L)
pc_scores <- round(pca_out$x[, 1:2], 4)

ldx_pheno <- data.frame(
  Sample = rownames(ldx_geno),
  Trait1 = trait1_val,
  Trait2 = trait2_val,
  PC1    = pc_scores[, 1],
  PC2    = pc_scores[, 2],
  stringsAsFactors = FALSE
)

write.csv(ldx_pheno,
          file = "inst/extdata/example_phenotype.csv",
          row.names = FALSE, quote = FALSE
)

## 5b. Pre-adjusted phenotype means (BLUEs) for run_haplotype_prediction demos
## ldx_blues was defined and saved as .rda before this section.
## Here we write the flat-file copy.
write.csv(ldx_blues,
          file = "inst/extdata/example_blues.csv",
          row.names = FALSE, quote = FALSE
)

## 5c. Per-environment BLUEs for run_haplotype_stability() demos
## ldx_blues_list contains two named numeric vectors (env1, env2).
## The flat-file version is a wide data frame for readability.
blues_env_df <- data.frame(
  id      = names(ldx_blues_list$env1),
  YLD_env1 = ldx_blues_list$env1,
  YLD_env2 = ldx_blues_list$env2,
  row.names = NULL, stringsAsFactors = FALSE
)
write.csv(blues_env_df,
          file = "inst/extdata/example_blues_env.csv",
          row.names = FALSE, quote = FALSE
)

## 6. Remove any stale GDS cache files
## These are auto-generated by read_geno() during tests or previous sessions.
## Deleting them here forces clean re-conversion via SNPRelate on the next
## test run, ensuring the GDS schema matches what snpgdsGetGeno() expects.
stale_gds <- list.files("inst/extdata", pattern = "\\.gds$", full.names = TRUE)
if (length(stale_gds)) {
  file.remove(stale_gds)
  message("Removed ", length(stale_gds), " stale GDS cache file(s) from inst/extdata/")
}

message("Written flat-file examples to inst/extdata/")
message("Data generation complete.")
