# Cross-Population Haplotype Effect Concordance

Given two `allele_tests` data frames produced by
[`test_block_haplotypes`](https://FAkohoue.github.io/LDxBlocks/reference/test_block_haplotypes.md)
on two independent populations, computes per-block (and optionally
per-trait) statistics that quantify how consistently haplotype allele
effects replicate across populations.

**Statistics computed per block:**

- Effect correlation (\\r\\):

  Pearson correlation of per-allele effect sizes between populations
  across all shared alleles. Requires at least 3 shared alleles; `NA`
  otherwise.

- Direction agreement:

  Proportion of shared alleles where both populations assign the same
  sign to the effect. A value \>= 0.75 (i.e. at least 3 out of 4 alleles
  agree in direction) is considered strong directional replication.

- Inverse-variance weighted (IVW) meta-analytic effect:

  The weighted mean effect \\\hat\beta\_{\mathrm{IVW}} = \sum w_i
  \beta_i / \sum w_i\\ where \\w_i = 1/\mathrm{SE}\_i^2\\. Computed
  separately for each shared allele and summarised as the mean of
  per-allele IVW estimates. The IVW meta-analysis is the same framework
  used in two-sample Mendelian randomisation and cross-population GWAS
  meta-analysis (see Borenstein et al. 2009).

- Cochran's Q heterogeneity statistic:

  Tests whether effect sizes differ significantly between the two
  populations: \\Q = \sum w_i (\beta_i - \hat\beta\_{\mathrm{IVW}})^2\\.
  Under the null of no heterogeneity, \\Q \sim \chi^2\_{n-1}\\.
  Significant Q (large \\Q_p\\) indicates that effect sizes differ
  between populations - a sign of GxE interaction, LD structure
  differences, or population-specific allelic action.

- \\I^2\\ inconsistency:

  \\I^2 = 100 \times \max(0,\\ (Q - df)/Q)\\. Values \> 50% indicate
  substantial between-population heterogeneity.

- Block boundary concordance:

  When `blocks_pop1` and `blocks_pop2` are both supplied, the bp overlap
  ratio of the two populations' block definitions is computed for each
  block. A ratio \< 0.8 flags blocks where LD structure likely differs
  between populations and effect comparisons should be interpreted
  cautiously.

**Addressing the two principal limitations of cross-population
validation:**

1.  *Population structure confounding*: The Q+K mixed model in
    [`test_block_haplotypes()`](https://FAkohoue.github.io/LDxBlocks/reference/test_block_haplotypes.md)
    already corrects for within-population structure via the haplotype
    GRM. Cross-population validation inherently controls
    between-population confounding by design: the same haplotype allele
    must associate with the phenotype in a genetically distinct
    background, making false-positive carry-over from Pop A's
    stratification implausible. The meta-analytic \\I^2\\ statistic
    additionally flags cases where the effect sizes are so heterogeneous
    that a shared biological mechanism is unlikely.

2.  *LD block boundary differences*: Block boundaries can differ between
    populations due to different historical recombination rates and
    ancestral haplotype structure. The `boundary_overlap_ratio` column
    quantifies this directly for every block. Low overlap (\< 0.8) means
    the two populations carve the same genomic region into blocks of
    different sizes; in that case a haplotype allele tested over a 50 kb
    window in Pop A is compared to an allele over a 30 kb window in Pop
    B and the strings will not match well. The `n_shared_alleles` column
    is the empirical consequence: very low shared-allele counts despite
    adequate allele frequencies in both populations are a direct symptom
    of block boundary mismatch. The recommended remedy - passing Pop A's
    block table as the `blocks` argument to
    [`test_block_haplotypes()`](https://FAkohoue.github.io/LDxBlocks/reference/test_block_haplotypes.md)
    in Pop B - eliminates this issue entirely by forcing both runs to
    use identical coordinates.

## Usage

``` r
compare_block_effects(
  assoc_pop1,
  assoc_pop2,
  pop1_name = "pop1",
  pop2_name = "pop2",
  traits = NULL,
  min_shared_alleles = 2L,
  blocks_pop1 = NULL,
  blocks_pop2 = NULL,
  block_match = c("id", "position"),
  overlap_min = 0.5,
  direction_threshold = 0.75,
  boundary_overlap_warn = 0.8,
  verbose = TRUE
)
```

## Arguments

- assoc_pop1:

  Output of
  [`test_block_haplotypes`](https://FAkohoue.github.io/LDxBlocks/reference/test_block_haplotypes.md)
  for population 1 (discovery). Must contain an `allele_tests` data
  frame with columns `block_id`, `CHR`, `start_bp`, `end_bp`, `trait`,
  `allele`, `effect`, `SE`, `p_wald`.

- assoc_pop2:

  Output of
  [`test_block_haplotypes`](https://FAkohoue.github.io/LDxBlocks/reference/test_block_haplotypes.md)
  for population 2 (validation). Same structure as `assoc_pop1`.

- pop1_name:

  Character. Label for population 1 in the output. Default `"pop1"`.

- pop2_name:

  Character. Label for population 2 in the output. Default `"pop2"`.

- traits:

  Character vector or `NULL`. Traits to compare. If `NULL` (default),
  all traits present in both result objects are included.

- min_shared_alleles:

  Integer. Minimum number of alleles shared between the two populations
  for a block to be included in the output. Blocks with fewer shared
  alleles are retained but marked with `enough_shared = FALSE`. Default
  `2L`.

- blocks_pop1:

  Optional data frame of LD blocks from population 1 (output of
  [`run_Big_LD_all_chr`](https://FAkohoue.github.io/LDxBlocks/reference/run_Big_LD_all_chr.md)).
  When both `blocks_pop1` and `blocks_pop2` are supplied, a
  `boundary_overlap_ratio` column is computed for every block. Required
  columns: `block_id` (or constructible from `CHR`/`start.bp`/`end.bp`),
  `start.bp`, `end.bp`.

- blocks_pop2:

  Optional data frame of LD blocks from population 2. Same structure as
  `blocks_pop1`.

- block_match:

  Character. How to match blocks between populations.

  - `"id"` (default) - match by exact `block_id` string. Fast and
    backward-compatible. Use when both populations were analysed with
    the same block table (recommended workflow: run
    [`test_block_haplotypes()`](https://FAkohoue.github.io/LDxBlocks/reference/test_block_haplotypes.md)
    on both populations using Pop A's block boundaries for Pop B as
    well).

  - `"position"` - match by genomic interval overlap
    (Intersection-over-Union, IoU, in base pairs). For each Pop1 block,
    the best-matching Pop2 block on the same chromosome is found. Blocks
    with IoU \\\geq\\ `overlap_min` are matched; those below are
    labelled `"pop1_only"`. Use when block boundaries genuinely differ
    between populations (different ancestral LD structures, different
    `CLQcut` used, or independent block detection runs).

- overlap_min:

  Numeric in (0, 1\]. Minimum Intersection-over-Union (IoU) in base
  pairs for two blocks to be considered the same region when
  `block_match = "position"`. Default `0.50`. Blocks below this
  threshold are assigned `match_type = "pop1_only"` and not compared.
  Lower values (e.g. `0.30`) tolerate more boundary discordance; higher
  values (e.g. `0.80`) require tighter boundary agreement. Ignored when
  `block_match = "id"`.

- direction_threshold:

  Numeric in (0.5, 1\]. Minimum direction-agreement proportion to
  consider a block directionally concordant. Default `0.75`.

- boundary_overlap_warn:

  Numeric in (0, 1). Threshold used to raise a warning flag in the
  output. When the automatically computed `boundary_overlap_ratio` (an
  output column, not a user-set value) is below this threshold, the
  corresponding `boundary_warning` output column is set to `TRUE`.
  Default `0.80`. Raise this value (e.g. `0.90`) to flag more
  conservatively, or lower it (e.g. `0.60`) to flag only severely
  mismatched blocks.

- verbose:

  Logical. Print progress. Default `TRUE`.

## Value

A named list of class `c("LDxBlocks_effect_concordance", "list")`:

- `concordance`:

  Data frame with one row per block per trait. Columns:

  - `block_id`, `CHR`, `start_bp`, `end_bp`, `trait` - block coordinates
    and trait name.

  - `n_alleles_pop1`, `n_alleles_pop2` - alleles tested in each
    population (before intersection).

  - `n_shared_alleles` - alleles present in both populations.

  - `enough_shared` - logical; `TRUE` when
    `n_shared_alleles >= min_shared_alleles`.

  - `effect_correlation` - Pearson r of per-allele effects across
    populations (NA when n_shared \< 3).

  - `direction_agreement` - proportion of shared alleles with concordant
    effect signs.

  - `directionally_concordant` - logical; direction_agreement \>=
    `direction_threshold`.

  - `meta_effect` - IVW meta-analytic effect (mean over shared alleles,
    each weighted by combined inverse-variance).

  - `meta_SE` - SE of the IVW estimate.

  - `meta_z` - meta-analytic z-score.

  - `meta_p` - two-sided p-value of meta-analytic effect.

  - `Q_stat` - Cochran Q heterogeneity statistic.

  - `Q_df` - degrees of freedom of Q (n_shared - 1).

  - `Q_p` - p-value of Q under chi-squared distribution.

  - `I2` - \\I^2\\ inconsistency (0-100%).

  - `replicated` - logical; `TRUE` when `directionally_concordant` AND
    `Q_p > 0.05` AND `enough_shared`.

  - `boundary_overlap_ratio` - bp intersection / bp union of the two
    populations' block boundaries. `NA` when block tables not supplied.

  - `boundary_warning` - logical; `TRUE` when
    `boundary_overlap_ratio < boundary_overlap_warn`. These blocks
    should be interpreted cautiously because different LD structures
    likely produce non-comparable haplotype strings.

  - `match_type` - character; how the block was matched between
    populations: `"exact"` (same `block_id` string), `"position"`
    (matched by genomic IoU \\\geq\\ `overlap_min` when
    `block_match = "position"`), `"pop1_only"` (no Pop2 block overlapped
    at the threshold), or `NA` when no block tables were supplied.

  Sorted by `CHR`, `start_bp`, `meta_p` (ascending).

- `shared_alleles`:

  Data frame. One row per shared allele per block per trait. Contains
  `effect_pop1`, `SE_pop1`, `p_wald_pop1`, `effect_pop2`, `SE_pop2`,
  `p_wald_pop2`, `direction_agree`, `ivw_effect`, `ivw_SE`, for detailed
  per-allele inspection.

- `pop1_name`, `pop2_name`:

  Character. Population labels.

- `block_match`:

  Character. Matching strategy used.

- `overlap_min`:

  Numeric. IoU threshold used (relevant when
  `block_match = "position"`).

- `traits`:

  Character vector of traits compared.

- `direction_threshold`:

  Numeric. Threshold used.

- `boundary_overlap_warn`:

  Numeric. Boundary warning threshold.

## References

Borenstein M, Hedges LV, Higgins JPT, Rothstein HR (2009). *Introduction
to Meta-Analysis*. Wiley.

Higgins JPT, Thompson SG (2002). Quantifying heterogeneity in a
meta-analysis. *Statistics in Medicine* **21**(11):1539-1558.
[doi:10.1002/sim.1186](https://doi.org/10.1002/sim.1186)

## See also

[`test_block_haplotypes`](https://FAkohoue.github.io/LDxBlocks/reference/test_block_haplotypes.md),
[`harmonize_haplotypes`](https://FAkohoue.github.io/LDxBlocks/reference/harmonize_haplotypes.md),
[`score_favorable_haplotypes`](https://FAkohoue.github.io/LDxBlocks/reference/score_favorable_haplotypes.md)

## Examples

``` r
# \donttest{
data(ldx_geno, ldx_snp_info, ldx_blocks, ldx_blues, package = "LDxBlocks")

# Simulate two populations by splitting the example dataset
set.seed(1L)
n     <- nrow(ldx_geno)
idx_1 <- sample(n, round(n * 0.6))
idx_2 <- setdiff(seq_len(n), idx_1)

haps_1 <- extract_haplotypes(ldx_geno[idx_1, ], ldx_snp_info, ldx_blocks)
haps_2 <- extract_haplotypes(ldx_geno[idx_2, ], ldx_snp_info, ldx_blocks)

blues_1 <- setNames(ldx_blues$YLD[idx_1], ldx_blues$id[idx_1])
blues_2 <- setNames(ldx_blues$YLD[idx_2], ldx_blues$id[idx_2])

res_1 <- test_block_haplotypes(haps_1, blues = blues_1,
                               blocks = ldx_blocks, verbose = FALSE)
res_2 <- test_block_haplotypes(haps_2, blues = blues_2,
                               blocks = ldx_blocks, verbose = FALSE)

conc <- compare_block_effects(res_1, res_2,
                              pop1_name = "pop1", pop2_name = "pop2",
                              blocks_pop1 = ldx_blocks,
                              blocks_pop2 = ldx_blocks)
#> [compare_block_effects] Comparing traits: trait
#> [compare_block_effects] Done. Blocks compared: 9 | Replicated (concordant, Q_p > 0.05): 1

# Replicated blocks
conc$concordance[conc$concordance$replicated, ]
#>               block_id CHR start_bp end_bp trait n_alleles_pop1 n_alleles_pop2
#> 5 block_2_86236_105290   2    86236 105290 trait              8              5
#>   n_shared_alleles enough_shared effect_correlation direction_agreement
#> 5                5          TRUE             0.6684                 0.8
#>   directionally_concordant meta_effect  meta_SE  meta_z    meta_p Q_stat Q_df
#> 5                     TRUE   -0.054497 0.102792 -0.5302 0.5959947 1.4741    4
#>         Q_p I2 replicated boundary_overlap_ratio boundary_warning match_type
#> 5 0.8312194  0       TRUE                      1            FALSE      exact

# Full per-allele details
head(conc$shared_alleles)
#>             block_id CHR start_bp end_bp trait                  allele
#> 1 block_1_1000_25027   1     1000  25027 trait block_1_1000_25027_hap1
#> 2 block_1_1000_25027   1     1000  25027 trait block_1_1000_25027_hap2
#> 3 block_1_1000_25027   1     1000  25027 trait block_1_1000_25027_hap3
#> 4 block_1_1000_25027   1     1000  25027 trait block_1_1000_25027_hap5
#> 5 block_1_1000_25027   1     1000  25027 trait block_1_1000_25027_hap4
#> 6 block_1_1000_25027   1     1000  25027 trait block_1_1000_25027_hap6
#>   effect_pop1  SE_pop1 p_wald_pop1 effect_pop2  SE_pop2 p_wald_pop2
#> 1   -0.397472 0.251976   0.1192082    0.211897 0.347946   0.5455243
#> 2   -0.310596 0.273817   0.2605300   -0.174085 0.401162   0.6663514
#> 3    0.237420 0.302744   0.4355541    0.757311 0.501033   0.1375008
#> 4    0.226424 0.321419   0.4834891   -0.729625 0.501925   0.1528313
#> 5    0.143343 0.303588   0.6382789   -0.338970 0.510885   0.5103244
#> 6   -0.181619 0.416618   0.6642237    0.562911 0.506570   0.2722476
#>   direction_agree ivw_effect   ivw_SE
#> 1           FALSE  -0.187836 0.204082
#> 2            TRUE  -0.267210 0.226157
#> 3            TRUE   0.376468 0.259115
#> 4           FALSE  -0.051613 0.270676
#> 5           FALSE   0.017475 0.260985
#> 6           FALSE   0.118784 0.321774

# Summary
print(conc)
#> LDxBlocks Cross-Population Effect Concordance
#>   Populations: pop1 vs pop2
#>   Traits:       trait 
#>   Blocks compared:           9 
#>   With enough shared alleles: 9 
#>   Directionally concordant:   1 
#>   Replicated (dir + Q_p>0.05): 1 
#>   Boundary warnings:          0 (overlap ratio < 0.8 )
#>   Median I2 (heterogeneity):  14.1 %
#>   Shared allele comparisons:  61 
# }
```
