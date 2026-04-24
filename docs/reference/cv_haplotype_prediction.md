# K-Fold Cross-Validation for Haplotype-Based Genomic Prediction

Estimates the predictive ability of the haplotype GBLUP model via k-fold
cross-validation. In each fold, a subset of individuals is masked from
the phenotype and predicted from the haplotype GRM; Pearson correlation
between predicted and observed BLUEs is returned as the predictive
ability (PA). Runs per trait when multiple traits are supplied.

## Usage

``` r
cv_haplotype_prediction(
  geno_matrix,
  snp_info,
  blocks,
  blues,
  k = 5L,
  n_rep = 1L,
  top_n = NULL,
  min_freq = 0.05,
  min_snps = 3L,
  id_col = "id",
  blue_col = "blue",
  blue_cols = NULL,
  seed = 42L,
  verbose = TRUE
)
```

## Arguments

- geno_matrix:

  Numeric matrix (individuals x SNPs), MAF-filtered dosage.

- snp_info:

  Data frame with columns `SNP`, `CHR`, `POS`.

- blocks:

  LD block table from
  [`run_Big_LD_all_chr`](https://FAkohoue.github.io/LDxBlocks/reference/run_Big_LD_all_chr.md).

- blues:

  Pre-adjusted phenotype means. Accepts the same four formats as
  [`run_haplotype_prediction`](https://FAkohoue.github.io/LDxBlocks/reference/run_haplotype_prediction.md):
  named numeric vector, single-trait data frame, multi-trait data frame,
  or named list.

- k:

  Integer. Number of folds. Default `5L`.

- n_rep:

  Integer. Number of CV replications (each with a different random fold
  assignment). Default `1L`.

- top_n:

  Integer or `NULL`. Maximum haplotype alleles per block passed to
  [`build_haplotype_feature_matrix`](https://FAkohoue.github.io/LDxBlocks/reference/build_haplotype_feature_matrix.md).
  Default `NULL` (all alleles above `min_freq`).

- min_freq:

  Numeric. Minimum haplotype allele frequency. Default `0.05`.

- min_snps:

  Integer. Minimum SNPs per block for haplotype extraction. Default
  `3L`.

- id_col:

  Character. Name of the individual ID column when `blues` is a data
  frame. Default `"id"`.

- blue_col:

  Character. Name of the BLUE column for single-trait data frames.
  Default `"blue"`.

- blue_cols:

  Character vector. Trait column names for multi-trait data frames.
  Default `NULL` (auto-detect all numeric non-ID columns).

- seed:

  Integer. RNG seed for reproducible fold assignment. Default `42L`.

- verbose:

  Logical. Print progress. Default `TRUE`.

## Value

A named list of class `LDxBlocks_cv`:

- `pa_summary`:

  Data frame: `trait`, `rep`, `fold`, `n_train`, `n_test`, `PA` (Pearson
  r), `RMSE`.

- `pa_mean`:

  Data frame: mean PA and RMSE per trait across all folds and
  replications.

- `gebv_all`:

  Data frame of out-of-fold GEBVs for all individuals and traits (one
  row per individual x trait).

- `k`:

  Number of folds used.

- `n_rep`:

  Number of replications.

## See also

[`run_haplotype_prediction`](https://FAkohoue.github.io/LDxBlocks/reference/run_haplotype_prediction.md),
[`build_haplotype_feature_matrix`](https://FAkohoue.github.io/LDxBlocks/reference/build_haplotype_feature_matrix.md)

## Examples

``` r
# \donttest{
data(ldx_geno, ldx_snp_info, ldx_blocks, ldx_blues, package = "LDxBlocks")
cv <- cv_haplotype_prediction(
  geno_matrix = ldx_geno,
  snp_info    = ldx_snp_info,
  blocks      = ldx_blocks,
  blues       = ldx_blues,
  k           = 5L,
  id_col      = "id",
  verbose     = FALSE
)
cv$pa_mean
#>   trait        PA      RMSE     PA_sd    RMSE_sd
#> 1   RES 0.3208779 0.9430814 0.1538289 0.07182249
#> 2   YLD 0.1742259 1.0030557 0.1953338 0.09088870
# }
```
