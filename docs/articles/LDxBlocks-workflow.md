# LDxBlocks Workflow: From Genotype File to Haplotype Features

## Overview

This vignette walks through the complete LDxBlocks workflow, from
reading a genotype file to obtaining haplotype feature matrices and
diversity statistics. The pipeline proceeds in four stages:

1.  **Read** genotype data via
    [`read_geno()`](https://FAkohoue.github.io/LDxBlocks/reference/read_geno.md)
    or the all-in-one
    [`run_ldx_pipeline()`](https://FAkohoue.github.io/LDxBlocks/reference/run_ldx_pipeline.md)
2.  **Detect** LD blocks chromosome-by-chromosome with
    [`run_Big_LD_all_chr()`](https://FAkohoue.github.io/LDxBlocks/reference/run_Big_LD_all_chr.md)
3.  **Extract** haplotypes and compute diversity with
    [`extract_haplotypes()`](https://FAkohoue.github.io/LDxBlocks/reference/extract_haplotypes.md)
    and
    [`compute_haplotype_diversity()`](https://FAkohoue.github.io/LDxBlocks/reference/compute_haplotype_diversity.md)
4.  **Predict** or **associate**: build a feature matrix for genomic
    prediction or run haplotype association tests

All examples use `eval = FALSE` so the vignette builds without requiring
the large input files. Adapt paths and parameters to your dataset.

------------------------------------------------------------------------

## Reading genotype data

### Supported formats

[`read_geno()`](https://FAkohoue.github.io/LDxBlocks/reference/read_geno.md)
auto-detects format from the file extension:

``` r
library(LDxBlocks)

# VCF (recommended for WGS data; auto-converted to GDS cache on first call)
be_vcf  <- read_geno("mydata.vcf.gz")

# PLINK BED (requires BEDMatrix)
be_bed  <- read_geno("mydata.bed")

# SNPRelate GDS (requires SNPRelate; fastest for large panels)
be_gds  <- read_geno("mydata.gds")

# Numeric dosage CSV (SNP × individuals, values in {0,1,2,NA})
be_csv  <- read_geno("mydata.csv", format = "numeric")

# Plain R matrix (backward-compatible)
be_mat  <- read_geno(my_matrix, snp_info = my_snp_df)
```

The `LDxBlocks_backend` object stores metadata without loading the full
genotype matrix:

``` r
be_vcf
# LDxBlocks backend
#   Type    : vcf
#   Samples : 204
#   SNPs    : 2960965
#   Chr     : 1, 2, 3, ..., 12
```

Always call
[`close_backend()`](https://FAkohoue.github.io/LDxBlocks/reference/close_backend.md)
when finished:

``` r
close_backend(be_vcf)
```

### Chromosome label normalisation

At read time, the leading `chr` / `Chr` / `CHR` prefix is stripped from
chromosome labels. Polyploid sub-genome labels (`1A`, `2D`) are
preserved verbatim:

| Input label            | Stored as |
|------------------------|-----------|
| `chr1`, `Chr1`, `CHR1` | `1`       |
| `Chr1A`                | `1A`      |
| `Chr2D`                | `2D`      |

------------------------------------------------------------------------

## LD block detection

### Recommended settings

**For chip-density panels (\< 100k SNPs per chromosome):**

``` r
blocks <- run_Big_LD_all_chr(
  geno_matrix = be,
  method      = "r2",         # standard r² — no kinship correction
  CLQcut      = 0.70,
  CLQmode     = "Density",    # exact Bron-Kerbosch clique enumeration
  subSegmSize = 1500L,
  leng        = 200L,
  MAFcut      = 0.05,
  n_threads   = 8L,
  verbose     = TRUE
)
```

**For WGS panels (2M–10M+ SNPs):**

``` r
blocks <- run_Big_LD_all_chr(
  geno_matrix     = be,
  method          = "r2",
  CLQcut          = 0.80,
  CLQmode         = "Leiden",       # guaranteed connected communities, O(n log n)
  max_bp_distance = 500000L,        # skip pairs > 500 kb (70-90% of pairs)
  subSegmSize     = 500L,
  leng            = 50L,
  checkLargest    = TRUE,
  n_threads       = 8L,
  verbose         = TRUE
)
```

**For related populations (livestock, inbred lines):**

``` r
blocks <- run_Big_LD_all_chr(
  geno_matrix = be,
  method      = "rV2",         # kinship-adjusted LD — requires AGHmatrix + ASRgenomics
  kin_method  = "chol",        # Cholesky whitening (default); "eigen" for near-singular GRM
  CLQcut      = 0.70,
  n_threads   = 8L
)
```

### Block table structure

``` r
head(blocks)
#   start  end start.rsID end.rsID start.bp   end.bp CHR length_bp
# 1     1   18      rs001    rs018     1000   103000   1    102001

summarise_blocks(blocks)
#    CHR n_blocks min_bp median_bp mean_bp max_bp total_bp_covered
# 1    1     1214   2000     82000   88940 512000        108034316
```

### Parameter auto-tuning with GWAS markers

When GWAS-significant markers are available,
[`tune_LD_params()`](https://FAkohoue.github.io/LDxBlocks/reference/tune_LD_params.md)
finds the `CLQcut` value that best captures those markers within blocks:

``` r
result <- tune_LD_params(
  geno_matrix    = my_geno,
  snp_info       = my_snp_info,
  gwas_df        = my_gwas,           # data.frame: Marker, CHR, POS
  prefer_perfect = TRUE,              # zero unassigned + zero forced first
  target_bp_band = c(5e4, 5e5),
  parallel       = FALSE,
  seed           = 42L
)

result$best_params       # selected parameters
result$gwas_assigned     # GWAS markers with LD_block column
result$final_blocks      # block table from best_params
```

------------------------------------------------------------------------

## Haplotype extraction

### Phase-free (diploid allele strings)

[`extract_haplotypes()`](https://FAkohoue.github.io/LDxBlocks/reference/extract_haplotypes.md)
concatenates each individual’s allele codes (0, 1, or 2) for all SNPs
within a block:

``` r
haps <- extract_haplotypes(
  geno   = be,           # LDxBlocks_backend or plain matrix
  snp_info = be$snp_info,
  blocks = blocks,
  min_snps = 3L          # skip blocks with fewer than 3 SNPs
)

# Named list; one character vector (length n_individuals) per block
length(haps)      # number of haplotype blocks
haps[[1]][1:5]    # first 5 individuals, first block
# [1] "0120" "1100" "0210" "2100" "0021"
```

The `block_info` attribute contains block metadata:

``` r
bi <- attr(haps, "block_info")
head(bi)
#        block_id CHR start_bp  end_bp n_snps
# block_1_1000_103000   1     1000  103000     25
```

### Diversity metrics

``` r
div <- compute_haplotype_diversity(haps)

head(div[, c("block_id","CHR","n_snps","n_ind","n_haplotypes","He","Shannon","freq_dominant","sweep_flag")])

# Blocks under potential selective sweep
div[div$sweep_flag, c("block_id","CHR","start_bp","He","freq_dominant")]
```

### Haplotype feature matrix for genomic prediction

``` r
feat_obj <- build_haplotype_feature_matrix(
  haplotypes     = haps,
  top_n          = 5L,         # keep top 5 most frequent alleles per block
  min_freq       = 0.01,       # discard alleles < 1% frequency
  scale_features = TRUE        # standardise columns to zero mean, unit variance
)

feat <- feat_obj$matrix        # individuals × (n_blocks × top_n)
dim(feat)                      # e.g. 204 × 52105

# Haplotype GRM for GBLUP
G_hap <- tcrossprod(feat) / ncol(feat)
```

------------------------------------------------------------------------

## End-to-end pipeline: `run_ldx_pipeline()`

For the most common use case — VCF or BED input → blocks → haplotypes →
feature matrix —
[`run_ldx_pipeline()`](https://FAkohoue.github.io/LDxBlocks/reference/run_ldx_pipeline.md)
wraps all steps in a single call:

``` r
result <- run_ldx_pipeline(
  geno_source    = "mydata.vcf.gz",
  out_dir        = "ldx_results",
  out_blocks     = "ldx_results/blocks.csv",
  out_diversity  = "ldx_results/diversity.csv",
  out_hap_matrix = "ldx_results/hap_matrix.csv",
  hap_format     = "numeric",

  maf_cut        = 0.05,
  impute         = "mean_rounded",

  CLQcut         = 0.70,
  method         = "r2",
  CLQmode        = "Leiden",
  max_bp_distance = 500000L,
  subSegmSize    = 500L,
  leng           = 50L,
  n_threads      = 8L,

  min_snps_block = 3L,
  top_n          = 5L,

  verbose = TRUE
)

result$blocks          # LD block table
result$diversity       # per-block diversity
result$haplotypes      # raw haplotype strings
result$geno_matrix     # cleaned + imputed genotype matrix
result$phase_method    # "unphased"
```

### With bigmemory caching (fast restart)

For WGS datasets, enable bigmemory caching so subsequent runs reattach
the imputed genotype matrix from disk without re-reading the raw VCF:

``` r
result <- run_ldx_pipeline(
  geno_source    = "mydata.vcf.gz",
  out_dir        = "ldx_results",
  out_blocks     = "ldx_results/blocks.csv",
  out_hap_matrix = "ldx_results/hap_matrix.csv",
  CLQcut         = 0.70,
  CLQmode        = "Leiden",
  max_bp_distance = 500000L,
  n_threads      = 8L,

  use_bigmemory  = TRUE,
  bigmemory_path = "ldx_results/bm_cache",
  bigmemory_type = "char",     # 1 byte per cell; "double" for full precision

  verbose = TRUE
)
# Second run: reattaches from bm_cache instantly
```

------------------------------------------------------------------------

## Genomic prediction

### Single-trait GBLUP

``` r
pred <- run_haplotype_prediction(
  geno_matrix = geno,
  snp_info    = snp_info,
  blocks      = blocks,
  blues       = setNames(my_blues_df$YLD, my_blues_df$id),
  id_col      = "id",
  blue_col    = "YLD"
)

# GEBVs
head(pred$gebv[, c("id", "YLD")])

# Most important blocks
important <- pred$block_importance[pred$block_importance$important, ]
head(important[order(-important$var_scaled), c("block_id","CHR","var_scaled")])
```

### Multi-trait GBLUP

``` r
pred_mt <- run_haplotype_prediction(
  geno_matrix     = geno,
  snp_info        = snp_info,
  blocks          = blocks,
  blues           = my_blues_df,          # wide data frame
  id_col          = "id",
  blue_cols       = c("YLD", "DIS", "PHT"),
  importance_rule = "any"                 # important in at least one trait
)

# Blocks important for all three traits
all_imp <- pred_mt$block_importance[pred_mt$block_importance$important_all, ]
```

### Cross-validation

``` r
cv <- cv_haplotype_prediction(
  geno_matrix = geno,
  snp_info    = snp_info,
  blocks      = blocks,
  blues       = my_blues_vec,
  k           = 5L,
  n_rep       = 3L,
  seed        = 42L
)
cv$pa_mean     # mean PA ± SD per trait across folds
```

------------------------------------------------------------------------

## GWAS integration

### Define QTL regions

``` r
qtl <- define_qtl_regions(
  gwas_df       = my_gwas,          # data.frame: Marker, CHR, POS, P
  blocks        = blocks,
  snp_info      = snp_info,
  p_threshold   = 5e-8,
  verbose       = TRUE
)

# Pleiotropic blocks (QTL for multiple traits)
qtl[qtl$n_traits > 1, ]
```

### Integrate evidence for breeding priority

``` r
priority <- integrate_gwas_haplotypes(
  qtl_regions = qtl,
  pred_result = pred_mt,
  diversity   = div
)

# Top candidates: supported by GWAS + variance + diversity
priority[priority$priority_score == 3,
         c("block_id","CHR","start_bp","priority_score")]
```

### Export candidate regions for annotation

``` r
export_candidate_regions(
  qtl_regions = qtl,
  format      = "bed",
  chr_prefix  = "chr",
  out_file    = "candidate_regions.bed"
)
```

------------------------------------------------------------------------

## Haplotype association testing

### EMMAX model (pure GRM, recommended for continuous kinship)

``` r
blues_vec <- setNames(my_blues_df$YLD, my_blues_df$id)

# EMMAX: GRM random effect only, FDR-based discovery
assoc <- test_block_haplotypes(
  haplotypes = haps,
  blues      = blues_vec,
  blocks     = blocks,
  n_pcs      = 0L,          # EMMAX: GRM random effect only
  sig_metric = "p_fdr",     # Benjamini-Hochberg FDR (recommended for discovery)
  verbose    = FALSE
)

# Top blocks by omnibus significance
head(assoc$block_tests[order(assoc$block_tests$p_omnibus), ])

# Significant per-allele associations (FDR <= 0.05)
assoc$allele_tests[assoc$allele_tests$significant, ]
```

### Q+K model with simpleM correction (recommended for family-wise control)

Use `sig_metric = "p_simplem_sidak"` and `meff_scope = "chromosome"`
when you want family-wise error control that accounts for the
correlation among haplotype alleles (less conservative than raw
Bonferroni). All four p-value columns are always present in the output
regardless of `sig_metric`.

``` r
assoc_sm <- test_block_haplotypes(
  haplotypes       = haps,
  blues            = blues_vec,
  blocks           = blocks,
  n_pcs            = 3L,                 # Q+K: 3 GRM-derived PCs as fixed effects
  sig_metric       = "p_simplem_sidak",  # simpleM Šidák (recommended)
  meff_scope       = "chromosome",       # separate Meff per chromosome (recommended)
  meff_percent_cut = 0.995,              # 99.5% variance threshold (Gao default)
  verbose          = FALSE
)

# Per-chromosome effective test counts
assoc_sm$meff$YLD$allele$chromosome

# All four p-value columns always present
names(assoc_sm$allele_tests)
# ... "p_wald" "p_fdr" "p_simplem" "p_simplem_sidak" "Meff" ...
```

### Automatic PC model selection (optional, recommended when lambda != 1)

When genomic control lambda is deflated (\< 1, over-correction) or
inflated (\> 1, under-correction), `optimize_pcs = TRUE` lets the data
choose n_pcs:

``` r
assoc_opt <- test_block_haplotypes(
  haplotypes       = haps,
  blues            = blues_vec,
  blocks           = blocks,
  optimize_pcs     = TRUE,           # triggers automatic PC selection
  optimize_pcs_max = 10L,            # test k = 0..10
  optimize_method  = "bic_lambda",   # recommended: |lambda-1| + 0.01*BIC
  sig_metric       = "p_simplem_sidak",
  meff_scope       = "chromosome",
  plot             = TRUE,           # PDF plots: Manhattan, QQ, PCA, scree
  out_dir          = "my_results",
  verbose          = TRUE
)

# Model selection table: one row per k tested
assoc_opt$pc_model_selection
#  n_pcs     BIC lambda_gc   score selected
#      0  1234.1     1.021  0.0210    FALSE
#      1  1231.5     1.015  0.0150    FALSE
#      2  1230.8     1.003  0.0033    FALSE
#      3  1232.1     0.999  0.0011     TRUE

# n_pcs chosen
assoc_opt$n_pcs_used
```

### Diplotype effects (additive + dominance)

[`estimate_diplotype_effects()`](https://FAkohoue.github.io/LDxBlocks/reference/estimate_diplotype_effects.md)
now supports the full correction set with a `sig_metric` parameter,
mirroring
[`test_block_haplotypes()`](https://FAkohoue.github.io/LDxBlocks/reference/test_block_haplotypes.md):

``` r
dip <- estimate_diplotype_effects(
  haplotypes       = haps,
  blues            = blues_vec,
  blocks           = blocks,
  min_n_diplotype  = 3L,
  sig_metric       = "p_omnibus_simplem_sidak",  # recommended
  meff_percent_cut = 0.995,
  verbose          = FALSE
)

# All five correction columns always present in $omnibus_tests:
# p_omnibus_adj | p_omnibus_fdr | p_omnibus_simplem | p_omnibus_simplem_sidak | Meff
head(dip$omnibus_tests[, c("block_id","trait","p_omnibus","p_omnibus_fdr",
                            "p_omnibus_simplem_sidak","Meff","significant")], 5)

# Blocks showing overdominance (heterosis candidates)
dip$dominance_table[dip$dominance_table$overdominance,
                    c("block_id","allele_A","allele_B","a","d","d_over_a")]
```

### Cross-population effect concordance

[`compare_block_effects()`](https://FAkohoue.github.io/LDxBlocks/reference/compare_block_effects.md)
validates whether haplotype effects discovered in one population
replicate in another, using IVW meta-analysis, Cochran Q heterogeneity,
and direction agreement:

``` r
# Assumes assoc_popA and assoc_popB are test_block_haplotypes() results
# from two independent populations, both run using the same block boundaries

conc <- compare_block_effects(
  assoc_pop1            = assoc_popA,
  assoc_pop2            = assoc_popB,
  pop1_name             = "PopA",
  pop2_name             = "PopB",
  blocks_pop1           = blocks_popA,
  blocks_pop2           = blocks_popB, # Pop B own block table (may differ)
  block_match           = "position",  # match by genomic IoU; match_type column reports result
  overlap_min           = 0.50,   # min IoU for 'position' match (else 'pop1_only')
  direction_threshold   = 0.75,
  boundary_overlap_warn = 0.80
)

# boundary_overlap_ratio is an OUTPUT column (automatically computed from block
# tables); boundary_overlap_warn is the INPUT threshold (default 0.80) that
# controls when boundary_warning = TRUE is set.

# Replicated blocks: directionally concordant AND Q_p > 0.05
conc$concordance[conc$concordance$replicated, ]

# Per-allele IVW detail
head(conc$shared_alleles)

print(conc)   # summary: n compared, n replicated, median I²
```

------------------------------------------------------------------------

## Breeding decision tools

### Score individuals for haplotype stacking

``` r
# Get allele effects from prediction pipeline
ae <- decompose_block_effects(
  haplotypes  = haps,
  snp_info    = snp_info,
  blocks      = blocks,
  snp_effects = pred$snp_effects[["YLD"]]
)

scores <- score_favorable_haplotypes(haps, allele_effects = ae, normalize = TRUE)

# Top 10 selection candidates
head(scores[order(scores$rank), c("id","stacking_index","n_blocks_scored","rank")], 10)
```

### Parent haplotype inventory

``` r
top10 <- scores$id[scores$rank <= 10]

inv <- summarize_parent_haplotypes(
  haplotypes    = haps,
  candidate_ids = top10,
  allele_effects = ae
)

# Parents carrying rare alleles with positive effects
inv[inv$dosage > 0 & inv$is_rare & inv$allele_effect > 0,
    c("id","block_id","allele","dosage","allele_freq","allele_effect")]
```

------------------------------------------------------------------------

## Between-population comparison

``` r
ids <- rownames(geno)

cmp <- compare_haplotype_populations(
  haplotypes  = haps,
  group1      = ids[ids %in% wild_ids],
  group2      = ids[ids %in% elite_ids],
  group1_name = "wild",
  group2_name = "elite"
)

# Divergent blocks (FST > 0.1, chi-squared p < 0.05)
cmp[cmp$divergent, c("block_id","FST","max_freq_diff","chisq_p")]
```

------------------------------------------------------------------------

## Cross-population validation from external GWAS

When GWAS was run in an external tool (GAPIT, TASSEL, FarmCPU, PLINK),
use
[`compare_gwas_effects()`](https://FAkohoue.github.io/LDxBlocks/reference/compare_gwas_effects.md)
rather than
[`compare_block_effects()`](https://FAkohoue.github.io/LDxBlocks/reference/compare_block_effects.md).

### Pre-mapped path (recommended)

``` r
# Map each population's markers to blocks first
qtl_popA <- define_qtl_regions(
  gwas_df     = my_gwas_A,
  blocks      = blocks,
  snp_info    = snp_info,
  p_threshold = 5e-8
)
qtl_popB <- define_qtl_regions(
  gwas_df     = my_gwas_B,
  blocks      = blocks,
  snp_info    = snp_info,
  p_threshold = 5e-8
)

conc_gwas <- compare_gwas_effects(
  qtl_pop1    = qtl_popA,
  qtl_pop2    = qtl_popB,
  blocks_pop1 = blocks_popA,
  blocks_pop2 = blocks_popB,
  block_match = "position",
  overlap_min = 0.50,
  pop1_name   = "PopA",
  pop2_name   = "PopB"
)

# Replicated: direction concordant AND meta_p <= 0.05
conc_gwas$concordance[conc_gwas$concordance$replicated, ]
print(conc_gwas)
```

### Raw GWAS path (convenience)

SE is derived automatically from BETA and P when absent. `Marker` is
accepted as an alias for `SNP`.

``` r
conc_raw <- compare_gwas_effects(
  gwas_pop1     = my_gwas_A,
  gwas_pop2     = my_gwas_B,
  blocks_pop1   = blocks,
  blocks_pop2   = blocks,
  snp_info_pop1 = snp_info,     # reused for pop2 when snp_info_pop2 = NULL
  p_threshold   = 5e-8,
  beta_col      = "BETA",       # adjust to match your GWAS output column names
  se_col        = "SE",         # NULL or absent -> derived from z-score
  p_col         = "P",
  pop1_name     = "PopA",
  pop2_name     = "PopB",
  verbose       = FALSE
)

# Extra output columns vs compare_block_effects():
#   lead_snp_pop1/pop2  — which SNP tagged the block in each population
#   lead_p_pop1/pop2    — lead SNP p-values
#   se_derived_pop1/pop2 — TRUE when SE was derived from z-score
#   both_pleiotropic    — TRUE when block is pleiotropic in both populations
head(conc_raw$concordance[, c("block_id","lead_snp_pop1","lead_snp_pop2",
                                "meta_effect","meta_p","direction_agreement",
                                "se_derived_pop1","replicated")])
```

## Epistasis detection

Within-block and between-block epistasis functions operate on
GRM-corrected REML residuals from the same null model as
[`test_block_haplotypes()`](https://FAkohoue.github.io/LDxBlocks/reference/test_block_haplotypes.md),
ensuring population-structure-corrected tests throughout.

### Within-block pairwise SNP interaction scan

``` r
# Scan all C(p,2) SNP pairs within significant blocks
# Tests H0: aa_ij = 0  in  y = mu + ai*xi + aj*xj + aa_ij*(xi*xj) + e
# Corrected by Bonferroni AND simpleM Sidak within each block
epi_within <- scan_block_epistasis(
  assoc              = assoc,          # test_block_haplotypes() result
  geno_matrix        = geno_matrix,
  snp_info           = snp_info,
  blocks             = blocks,
  blues              = blues,
  haplotypes         = haps,
  trait              = "YLD",
  sig_blocks         = NULL,           # uses significant_omnibus blocks
  max_snps_per_block = 300L,           # exhaustive up to 300 SNPs/block
  sig_metric         = "p_simplem_sidak",
  sig_threshold      = 0.05
)

# Output: class LDxBlocks_epistasis
print(epi_within)

# Significant interacting pairs
epi_within$results[epi_within$results$significant, ]

# Per-block summary: n_snps, n_pairs, n_significant, min_p
epi_within$scan_summary

# All correction columns always present
# p_wald | p_bonf | p_simplem | p_simplem_sidak | Meff
# significant | significant_bonf | significant_simplem | significant_simplem_sidak
```

### Between-block trans-haplotype epistasis scan

``` r
# Tests significant haplotype alleles x all other blocks
# O(n_sig x n_total_alleles) tests -- tractable at WGS scale
# Identifies genetic background dependence: a haplotype at block A
# only functions in the presence of a specific haplotype at block B
epi_between <- scan_block_by_block_epistasis(
  assoc      = assoc,
  haplotypes = haps,
  blues      = blues,
  blocks     = blocks,
  trait      = "YLD",
  sig_alleles = NULL,    # NULL = uses significant alleles from assoc
  sig_threshold = 0.05
)

print(epi_between)
epi_between$results[epi_between$results$significant, ]

# Results columns:
# block_i, allele_i, block_j, allele_j, CHR_i, CHR_j, same_chr
# aa_effect, SE, t_stat, p_wald, p_bonf, significant
```

### Single-block epistasis fine-mapping

``` r
# Fine-map specific interacting SNP pairs within one block.
# Requires pre-computed REML residuals from the null model.
# method = "auto": pairwise scan for p <= 200 SNPs, LASSO for larger blocks.

# Pre-compute residuals (example using the null model from test_block_haplotypes)
feat <- build_haplotype_feature_matrix(haps, min_freq = 0.05,
                                        encoding = "additive_012")$matrix
G_grm <- compute_haplotype_grm(feat)
common <- intersect(rownames(G_grm), names(blues))
y  <- blues[common]
fit <- rrBLUP::mixed.solve(y = y, K = G_grm[common, common], method = "REML")
y_resid <- y - as.numeric(fit$beta) - fit$u[common]

# Exhaustive pairwise scan for blocks with <= 200 SNPs
fine <- fine_map_epistasis_block(
  block_id    = "block_12_1054210_1086071",
  geno_matrix = geno_matrix,
  snp_info    = snp_info,
  blocks      = blocks,
  y_resid     = y_resid,
  method      = "auto"   # pairwise <= 200 SNPs, lasso otherwise
)
head(fine)

# Columns (pairwise): SNP_i, SNP_j, POS_i, POS_j, dist_bp,
#                     aa_effect, SE, t_stat, p_wald, p_bonf, significant
# Columns (lasso):   SNP_i, SNP_j, POS_i, POS_j, dist_bp, lasso_coef, selected
```

> **Statistical note.** With typical WGS sample sizes (n \< 300), the
> detection threshold for within-block Bonferroni correction over C(p,2)
> pairs is stringent. The between-block trans-haplotype scan is better
> powered because each haplotype dosage column aggregates multi-SNP
> variation, providing a stronger signal-to-noise ratio per test.

------------------------------------------------------------------------

## Session information

``` r
sessionInfo()
```

    #> R version 4.5.0 (2025-04-11 ucrt)
    #> Platform: x86_64-w64-mingw32/x64
    #> Running under: Windows 11 x64 (build 26200)
    #> 
    #> Matrix products: default
    #>   LAPACK version 3.12.1
    #> 
    #> locale:
    #> [1] LC_COLLATE=English_United States.utf8 
    #> [2] LC_CTYPE=English_United States.utf8   
    #> [3] LC_MONETARY=English_United States.utf8
    #> [4] LC_NUMERIC=C                          
    #> [5] LC_TIME=English_United States.utf8    
    #> 
    #> time zone: America/Bogota
    #> tzcode source: internal
    #> 
    #> attached base packages:
    #> [1] stats     graphics  grDevices utils     datasets  methods   base     
    #> 
    #> loaded via a namespace (and not attached):
    #>  [1] digest_0.6.39     desc_1.4.3        R6_2.6.1          fastmap_1.2.0    
    #>  [5] xfun_0.57         cachem_1.1.0      knitr_1.51        htmltools_0.5.9  
    #>  [9] rmarkdown_2.31    lifecycle_1.0.5   cli_3.6.6         sass_0.4.10      
    #> [13] pkgdown_2.2.0     textshaping_1.0.5 jquerylib_0.1.4   systemfonts_1.3.2
    #> [17] compiler_4.5.0    rstudioapi_0.18.0 tools_4.5.0       ragg_1.5.2       
    #> [21] bslib_0.10.0      evaluate_1.0.5    yaml_2.3.12       otel_0.2.0       
    #> [25] jsonlite_2.0.0    htmlwidgets_1.6.4 rlang_1.2.0       fs_2.0.1
