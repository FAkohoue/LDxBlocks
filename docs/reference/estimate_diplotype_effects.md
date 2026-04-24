# Estimate Diplotype Effects and Dominance Deviations Per LD Block

Fits per-block diplotype models to quantify additive and dominance
genetic effects at LD blocks using haplotype-derived diplotypes.

This function integrates haplotype structure with mixed-model correction
for relatedness to provide biologically interpretable measures of gene
action, including dominance and overdominance.

**Statistical framework:**

The analysis proceeds in three stages:

1.  A haplotype feature matrix is constructed using
    [`build_haplotype_feature_matrix`](https://FAkohoue.github.io/LDxBlocks/reference/build_haplotype_feature_matrix.md).

2.  A haplotype genomic relationship matrix (GRM) is computed via
    [`compute_haplotype_grm`](https://FAkohoue.github.io/LDxBlocks/reference/compute_haplotype_grm.md),
    capturing genome-wide similarity based on haplotype composition.

3.  For each trait:

    - A null mixed model is fitted: \$\$y = \mu + u + e\$\$ where \\u
      \sim N(0, G\sigma^2_u)\\ using
      [`rrBLUP::mixed.solve()`](https://rdrr.io/pkg/rrBLUP/man/mixed.solve.html).

    - Residuals (de-regressed phenotypes) are computed: \\y^\* = y -
      \hat{\mu} - \hat{u}\\

    - Diplotype effects are estimated from these residuals.

**Diplotype representation:**

Diplotypes are canonical strings of two haplotypes sorted alphabetically
and joined by `"/"`:

- Homozygote: `"010/010"`

- Heterozygote: `"010/110"`

This ensures phase invariance (`"010/110"` == `"110/010"`).

**Genetic effects:**

For each allele pair \\A, B\\, the following quantities are estimated:

- Additive effect: \$\$a = (\bar{y}\_{BB} - \bar{y}\_{AA}) / 2\$\$

- Dominance deviation: \$\$d = \bar{y}\_{AB} - (\bar{y}\_{AA} +
  \bar{y}\_{BB}) / 2\$\$

- Dominance ratio: \$\$d/a\$\$

Interpretation:

- \\d/a = 0\\ -\> additive

- \\\|d/a\| \< 1\\ -\> partial dominance

- \\\|d/a\| = 1\\ -\> complete dominance

- \\\|d/a\| \> 1\\ -\> overdominance (heterosis)

**Important:**

- Phased haplotypes are strongly recommended.

- Unphased heterozygotes are resolved arbitrarily
  (`resolve_unphased = TRUE`), which does not bias mean-based dominance
  estimates.

- GRM scaling is automatically adjusted based on whether haplotypes are
  phased.

## Usage

``` r
estimate_diplotype_effects(
  haplotypes,
  blues,
  blocks,
  min_freq = 0.05,
  min_n_diplotype = 3L,
  id_col = "id",
  blue_col = "blue",
  blue_cols = NULL,
  verbose = TRUE
)
```

## Arguments

- haplotypes:

  Named list from
  [`extract_haplotypes`](https://FAkohoue.github.io/LDxBlocks/reference/extract_haplotypes.md).
  Must include a `block_info` attribute with: `block_id`, `CHR`,
  `start_bp`, `end_bp`, and `phased`.

- blues:

  Phenotypic values (BLUEs or adjusted means). Supported formats:

  - Named numeric vector

  - Data frame (single or multi-trait)

  - Named list of numeric vectors

  Individual IDs must match haplotype IDs.

- blocks:

  LD block table from
  [`run_Big_LD_all_chr`](https://FAkohoue.github.io/LDxBlocks/reference/run_Big_LD_all_chr.md)
  (used for metadata).

- min_freq:

  Numeric in \[0, 1\]. Minimum haplotype allele frequency for inclusion
  in dominance decomposition. Default `0.05`.

- min_n_diplotype:

  Integer \>= 2. Minimum number of individuals per diplotype class.
  Default `3L`.

- id_col:

  Character. ID column name for data frame input. Default `"id"`.

- blue_col:

  Character. Single-trait column name. Default `"blue"`.

- blue_cols:

  Character vector for multi-trait input. Default `NULL`.

- verbose:

  Logical. Print progress messages. Default `TRUE`.

## Value

Object of class `LDxBlocks_diplotype` containing:

- diplotype_means:

  One row per diplotype class per block per trait.

  Columns:

  - block_id, CHR, start_bp, end_bp

  - trait

  - diplotype

  - n (class size)

  - n_total

  - mean_blue

  - se_mean

- dominance_table:

  One row per allele pair per block per trait.

  Columns:

  - allele_A, allele_B

  - mean_AA, mean_AB, mean_BB

  - a, d

  - d_over_a

  - overdominance (logical)

- omnibus_tests:

  Per-block ANOVA results:

  - F_stat, df1, df2

  - p_omnibus

  - p_omnibus_adj (Bonferroni per trait)

  - significant

## See also

[`compute_haplotype_grm`](https://FAkohoue.github.io/LDxBlocks/reference/compute_haplotype_grm.md),
[`infer_block_haplotypes`](https://FAkohoue.github.io/LDxBlocks/reference/infer_block_haplotypes.md),
[`test_block_haplotypes`](https://FAkohoue.github.io/LDxBlocks/reference/test_block_haplotypes.md)

## Examples

``` r
# \donttest{
if (requireNamespace("LDxBlocks", quietly = TRUE)) {
  data(ldx_geno, ldx_snp_info, ldx_blocks, ldx_blues, package = "LDxBlocks")

haps <- extract_haplotypes(
  geno     = ldx_geno,
  snp_info = ldx_snp_info,
  blocks   = ldx_blocks,
  min_snps = 5L
)

res <- estimate_diplotype_effects(
  haplotypes = haps,
  blues      = setNames(ldx_blues$YLD, ldx_blues$id),
  blocks     = ldx_blocks,
  verbose    = FALSE
)

subset(res$dominance_table, overdominance)
subset(res$omnibus_tests, significant)
}
#> [1] block_id      trait         n_diplotypes  F_stat        df1          
#> [6] df2           p_omnibus     p_omnibus_adj significant  
#> <0 rows> (or 0-length row.names)
# }
```
