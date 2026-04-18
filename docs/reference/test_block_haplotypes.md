# Block-Level Haplotype Association Testing (Q+K Mixed Linear Model)

Performs genome-wide haplotype block association tests for one or more
quantitative traits. Each LD block is tested as a unit: per-allele Wald
tests identify which specific haplotype alleles drive association, and
an omnibus F-test evaluates the block as a whole. Population structure
and kinship are corrected jointly through a unified mixed linear model.

**Statistical model (Q+K / EMMAX formulation):**

\$\$y = \mu + \alpha \cdot x\_{\mathrm{hap}} + \sum\_{k=1}^{K} \beta_k
\\ PC_k + g + \varepsilon\$\$

- \\x\_{\mathrm{hap}}\\:

  Haplotype allele dosage (0, 1, or 2 copies for phased data; 0 or 1 for
  unphased) - the quantity being tested for association.

- \\PC_k\\:

  The k-th eigenvector of the haplotype GRM, included as a fixed-effect
  covariate to capture discrete population structure explicitly. Derived
  from `eigen(G_hap)`, so both fixed-effect PCs and the random-effect
  GRM use the same kinship model.

- \\g \sim MVN(0,\\\sigma_g^2 G)\\:

  Polygenic background - captures residual continuous kinship
  (within-family, cryptic relatedness) as a random effect after PC
  removal.

- \\\varepsilon \sim MVN(0,\\\sigma_e^2 I)\\:

  Residual error.

**Implementation and scaling:** The GRM is inverted once per trait via
[`rrBLUP::mixed.solve()`](https://rdrr.io/pkg/rrBLUP/man/mixed.solve.html)
(dominant cost, O(n³)). Per-allele tests on the de-regressed residuals
are then fully vectorised across all blocks simultaneously using a
single
[`crossprod()`](https://rdrr.io/pkg/Matrix/man/matmult-methods.html)
call (O(n x p), analogous to the BLAS DGEMV trick in marginal SNP
screening). For 5,000 individuals and 51,000 haplotype allele columns
across 17,000 blocks, the scan step takes seconds after the one-time GRM
inversion (~30 s for n = 5,000).

## Usage

``` r
test_block_haplotypes(
  haplotypes,
  blues,
  blocks,
  n_pcs = 0L,
  top_n = NULL,
  min_freq = 0.05,
  id_col = "id",
  blue_col = "blue",
  blue_cols = NULL,
  alpha = NULL,
  verbose = TRUE
)
```

## Arguments

- haplotypes:

  Named list produced by
  [`extract_haplotypes`](https://FAkohoue.github.io/LDxBlocks/reference/extract_haplotypes.md).
  Each element is a named character vector (one haplotype dosage string
  per individual), and the list must carry a `block_info` attribute
  (added automatically by `extract_haplotypes`). The number of elements
  equals the number of qualifying LD blocks.

- blues:

  Pre-adjusted phenotype means (BLUEs or BLUPs from a field trial mixed
  model). Accepted in four formats:

  - Named numeric vector: `c(ind1 = 2.3, ind2 = 1.8, ...)`. Names must
    match individual IDs (row names of the genotype matrix).

  - Single-trait data frame: columns `id_col` (individual IDs) and
    `blue_col` (numeric phenotype values).

  - Multi-trait data frame: columns `id_col` plus one column per trait
    specified in `blue_cols`. All named traits are tested in a single
    call sharing the same GRM.

  - Named list of named numeric vectors: one element per trait, e.g.
    `list(YLD = c(ind1=2.3,...), RES = c(ind1=0.8,...))`.

  Individuals in `blues` not present in `haplotypes` are silently
  dropped. At least 10 common individuals are required per trait.

- blocks:

  LD block table returned by
  [`run_Big_LD_all_chr`](https://FAkohoue.github.io/LDxBlocks/reference/run_Big_LD_all_chr.md)
  or
  [`run_ldx_pipeline`](https://FAkohoue.github.io/LDxBlocks/reference/run_ldx_pipeline.md).
  Required columns: `CHR`, `start.bp`, `end.bp`. Used only for block
  metadata annotation in the output; the actual haplotype allele columns
  come from `haplotypes`.

- n_pcs:

  Integer (`>= 0`) or `NULL`. Number of haplotype-GRM eigenvectors to
  include as fixed-effect population structure covariates in the null
  model. Controls the trade-off between correction strategies:

  - `0L` (default): Pure GRM correction - EMMAX / P3D approximation. The
    GRM random effect absorbs all structure and kinship. Appropriate for
    populations with diffuse continuous kinship (e.g. livestock half-sib
    families, advanced inbred lines) where there are no sharp discrete
    subpopulation boundaries.

  - `1` to `10`: Q+K model (Yu et al. 2006) - top-k GRM eigenvectors as
    fixed effects plus the GRM random effect. The fixed- effect PCs
    capture discrete subpopulation membership (e.g. breeds, ecotypes,
    geographic clusters); the GRM random effect captures
    within-subpopulation continuous kinship. Use when strong cluster
    structure inflates the Q-Q plot under `n_pcs = 0`. Typically 3-5 PCs
    are sufficient; using more than 10 risks over-correction.

  - `NULL`: Auto-select via the elbow of the GRM eigenvalue scree plot
    (first position where the marginal gain in variance explained drops
    below 1%), capped at 10.

  PCs are derived from `eigen(G_hap)` - the same GRM that enters as the
  random effect - ensuring mathematical consistency.

- top_n:

  Integer or `NULL`. Maximum number of haplotype alleles per block to
  include in the feature matrix before testing. `NULL` (default) retains
  all alleles above `min_freq`. Supply an integer (e.g. `5L`) only to
  cap column count for very large panels where memory is constrained.
  Alleles are ranked by frequency; the most common `top_n` are retained.

- min_freq:

  Numeric in (0, 1). Minimum haplotype allele frequency in the panel.
  Alleles with frequency below this threshold are excluded from both the
  feature matrix and all tests. Default `0.05`. Lower values (e.g.
  `0.02`) include rarer alleles at the cost of reduced power and
  increased multiple testing burden. Values below `0.02` are not
  recommended for panels smaller than 200 individuals.

- id_col:

  Character. Name of the individual-ID column when `blues` is a data
  frame. Must exactly match a column name in the data frame. Default
  `"id"`.

- blue_col:

  Character. Name of the phenotype column when `blues` is a single-trait
  data frame. Must be numeric. Default `"blue"`.

- blue_cols:

  Character vector. Names of phenotype columns when `blues` is a
  multi-trait wide data frame. Each named column is treated as a
  separate trait and tested independently (but using the same GRM).
  Default `NULL` (ignored unless `blues` is a data frame with more than
  one numeric column beyond `id_col`).

- alpha:

  Numeric or `NULL`. Significance threshold for the `significant` flag
  in `allele_tests` and `significant_omnibus` in `block_tests`. `NULL`
  (default) applies genome-wide Bonferroni correction: \\\alpha = 0.05 /
  n\_{\mathrm{tests}}\\ where \\n\_{\mathrm{tests}}\\ is the total
  number of allele-level tests across all blocks and traits. Supply a
  fixed value (e.g. `0.05`) to use a less conservative threshold, or
  `0.05 / nrow(blocks)` to apply Bonferroni at the block level rather
  than the allele level.

- verbose:

  Logical. If `TRUE` (default), prints timestamped progress messages:
  model type, trait name, number of alleles scanned, and a final
  summary. Set `FALSE` for batch use or inside loops.

## Value

A named list of class `c("LDxBlocks_haplotype_assoc", "list")` with the
following elements:

- `allele_tests`:

  Data frame. One row per haplotype allele per LD block per trait.
  Contains 13 columns:

  - `block_id` (character) - Block identifier string matching
    `names(haplotypes)`, e.g. `"block_1_1000_103000"`.

  - `CHR` (character) - Chromosome label.

  - `start_bp` (integer) - Block start coordinate (base pairs).

  - `end_bp` (integer) - Block end coordinate (base pairs).

  - `trait` (character) - Trait name.

  - `allele` (character) - Haplotype allele identifier string.

  - `frequency` (numeric, (0,1\]) - Allele frequency in the panel
    (proportion of individuals carrying the allele).

  - `effect` (numeric) - Estimated additive effect: mean phenotype
    difference per unit increase in allele dosage on the de-regressed
    residual scale. Positive = favourable if higher trait values are
    desirable.

  - `SE` (numeric) - Standard error of the effect estimate.

  - `t_stat` (numeric) - t-statistic (`effect / SE`).

  - `p_wald` (numeric, (0,1\]) - Two-sided Wald p-value (raw,
    uncorrected).

  - `p_wald_adj` (numeric, (0,1\]) - Bonferroni-adjusted p-value:
    `min(p_wald * n_tests, 1)`.

  - `significant` (logical) - `TRUE` when `p_wald <= alpha`.

  Rows are sorted ascending by `CHR`, `start_bp`, then `p_wald` within
  each trait.

- `block_tests`:

  Data frame. One row per LD block per trait. Contains 12 columns:

  - `block_id`, `CHR`, `start_bp`, `end_bp`, `trait` - as above.

  - `n_alleles_tested` (integer) - Number of alleles above `min_freq` in
    this block that were included in the omnibus test.

  - `F_stat` (numeric) - Omnibus F-statistic testing all
    `n_alleles_tested` allele columns jointly against the de- regressed
    residual.

  - `df_LRT` (integer) - Numerator degrees of freedom (=
    `n_alleles_tested`).

  - `p_omnibus` (numeric, (0,1\]) - Raw omnibus p-value from the
    F-distribution.

  - `p_omnibus_adj` (numeric, (0,1\]) - Bonferroni-adjusted across all
    blocks per trait: `min(p_omnibus * n_blocks, 1)`.

  - `var_explained` (numeric, \[0,1\]) - Proportion of de-regressed
    phenotypic variance explained by all haplotype alleles of this block
    jointly: `1 - RSS_full / RSS_null`.

  - `significant_omnibus` (logical) - `TRUE` when
    `p_omnibus_adj < 0.05`.

  Rows sorted ascending by `CHR`, `start_bp`, `p_omnibus`.

- `traits`:

  Character vector of trait names that were tested.

- `n_pcs_used`:

  Integer. Number of GRM eigenvectors included as fixed-effect
  covariates in the null model. `0L` when `n_pcs = 0L` (pure EMMAX).

- `alpha`:

  Numeric. The significance threshold actually used (either supplied or
  computed by Bonferroni).

- `n_tests`:

  Integer. Total number of allele-level Wald tests performed across all
  blocks and traits (the Bonferroni denominator).

## References

Endelman JB (2011). Ridge regression and other kernels for genomic
selection with R package rrBLUP. *Plant Genome* **4**:250-255.
[doi:10.3835/plantgenome2011.08.0024](https://doi.org/10.3835/plantgenome2011.08.0024)

Yu J, Pressoir G, Briggs WH, et al. (2006). A unified mixed-model method
for association mapping that accounts for multiple levels of
relatedness. *Nature Genetics* **38**(2):203-208.
[doi:10.1038/ng1702](https://doi.org/10.1038/ng1702)

Kang HM, Sul JH, Service SK, et al. (2010). Variance component model to
account for sample structure in genome-wide association studies. *Nature
Genetics* **42**(4):348-354.
[doi:10.1038/ng.548](https://doi.org/10.1038/ng.548)

## See also

[`extract_haplotypes`](https://FAkohoue.github.io/LDxBlocks/reference/extract_haplotypes.md),
[`compute_haplotype_grm`](https://FAkohoue.github.io/LDxBlocks/reference/compute_haplotype_grm.md),
[`run_haplotype_prediction`](https://FAkohoue.github.io/LDxBlocks/reference/run_haplotype_prediction.md),
[`estimate_diplotype_effects`](https://FAkohoue.github.io/LDxBlocks/reference/estimate_diplotype_effects.md)

## Examples

``` r
# \donttest{
data(ldx_geno, ldx_snp_info, ldx_blocks, ldx_blues, package = "LDxBlocks")
haps <- extract_haplotypes(ldx_geno, ldx_snp_info, ldx_blocks)

# Pure GRM correction (EMMAX, default - n_pcs = 0)
assoc <- test_block_haplotypes(
  haplotypes = haps,
  blues      = setNames(ldx_blues$YLD, ldx_blues$id),
  blocks     = ldx_blocks,
  verbose    = FALSE
)
head(assoc$block_tests[order(assoc$block_tests$p_omnibus), ])
#>                block_id CHR start_bp end_bp trait n_alleles_tested F_stat
#> 7    block_3_1000_19068   3     1000  19068 trait                8 1.8388
#> 4    block_2_1000_30023   2     1000  30023 trait                7 1.3505
#> 8   block_3_74532_93854   3    74532  93854 trait                8 0.9660
#> 2   block_1_81064_99022   1    81064  99022 trait                7 0.9341
#> 1    block_1_1000_25027   1     1000  25027 trait                7 0.8244
#> 3 block_1_155368_179371   1   155368 179371 trait                9 0.7209
#>   df_LRT  p_omnibus var_explained p_omnibus_adj significant_omnibus
#> 7      8 0.07717218        0.1170     0.6945496               FALSE
#> 4      7 0.23347078        0.0778     1.0000000               FALSE
#> 8      8 0.46630057        0.0651     1.0000000               FALSE
#> 2      7 0.48322056        0.0552     1.0000000               FALSE
#> 1      7 0.56903648        0.0490     1.0000000               FALSE
#> 3      9 0.68877228        0.0557     1.0000000               FALSE

# Q+K model (3 GRM-derived PCs + GRM, for structured populations)
assoc_qk <- test_block_haplotypes(
  haplotypes = haps,
  blues      = setNames(ldx_blues$YLD, ldx_blues$id),
  blocks     = ldx_blocks,
  n_pcs      = 3L,
  verbose    = FALSE
)

# Multi-trait
assoc_mt <- test_block_haplotypes(
  haplotypes = haps,
  blues      = ldx_blues,
  blocks     = ldx_blocks,
  id_col     = "id",
  blue_cols  = c("YLD", "RES"),
  n_pcs      = 3L,
  verbose    = FALSE
)
# }
```
