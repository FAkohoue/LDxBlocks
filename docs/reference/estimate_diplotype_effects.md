# Estimate Diplotype Effects and Dominance Deviations Per LD Block

Fits a per-block diplotype model to decompose phenotypic variation into
additive and dominance components at each LD block, after correcting for
population structure and kinship via the haplotype GRM. This function
provides direct evidence for non-additive gene action: dominance,
overdominance, and heterosis.

**Method:** For each LD block, unique diplotypes (canonical
haplotype-pair strings inferred by
[`infer_block_haplotypes`](https://FAkohoue.github.io/LDxBlocks/reference/infer_block_haplotypes.md),
e.g. `"010/110"`) are treated as levels of a fixed factor. A mixed model
with the haplotype GRM as a kinship random effect is fitted via
[`rrBLUP::mixed.solve()`](https://rdrr.io/pkg/rrBLUP/man/mixed.solve.html)
to obtain de-regressed phenotype residuals. Diplotype class means are
computed on these residuals. For each pair of common alleles (A, B),
classical quantitative genetics parameters are derived from the three
diplotype class means:

- Additive effect: \\a = (\bar{y}\_{BB} - \bar{y}\_{AA}) / 2\\

- Dominance deviation: \\d = \bar{y}\_{AB} - (\bar{y}\_{AA} +
  \bar{y}\_{BB}) / 2\\

- Dominance ratio: \\d/a\\ - the key quantity for interpreting gene
  action (see `d_over_a` column description).

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

  Named list produced by
  [`extract_haplotypes`](https://FAkohoue.github.io/LDxBlocks/reference/extract_haplotypes.md).
  Phase-ambiguous diplotypes (unphased heterozygotes) are excluded from
  the dominance analysis but included in the omnibus F-test if their
  diplotype string is unambiguous. For accurate dominance decomposition,
  phased input (from `read_phased_vcf`, `phase_with_beagle`, or via
  [`phase_with_beagle`](https://FAkohoue.github.io/LDxBlocks/reference/phase_with_beagle.md))
  is strongly recommended.

- blues:

  Pre-adjusted phenotype means. Accepts the same four formats as
  [`test_block_haplotypes`](https://FAkohoue.github.io/LDxBlocks/reference/test_block_haplotypes.md):
  named numeric vector, single-trait data frame, multi-trait data frame,
  or named list. All formats require individual IDs matching the names
  of haplotype strings.

- blocks:

  LD block table from
  [`run_Big_LD_all_chr`](https://FAkohoue.github.io/LDxBlocks/reference/run_Big_LD_all_chr.md).
  Used only to supply block coordinate metadata (`CHR`, `start_bp`,
  `end_bp`) to the output tables.

- min_freq:

  Numeric in (0, 1). Minimum allele frequency. Alleles below this
  threshold in the full panel are excluded from the dominance
  decomposition (`dominance_table`) but their diplotypes still
  contribute to `diplotype_means` and `omnibus_tests` if sufficient
  individuals carry them. Default `0.05`.

- min_n_diplotype:

  Integer. Minimum number of individuals that must carry a given
  diplotype class for it to be included in `diplotype_means` and the
  dominance decomposition. Diplotype classes with fewer individuals are
  silently excluded. The omnibus F-test uses only the retained classes.
  Default `3L`. Increasing to 5-10 improves mean estimate stability for
  small panels; decreasing below 3 is not recommended.

- id_col:

  Character. Name of the individual-ID column when `blues` is a data
  frame. Default `"id"`.

- blue_col:

  Character. Name of the phenotype column for single-trait data frames.
  Default `"blue"`.

- blue_cols:

  Character vector. Phenotype column names for multi-trait data frames.
  Default `NULL`.

- verbose:

  Logical. `TRUE` (default) prints progress per trait.

## Value

A named list of class `c("LDxBlocks_diplotype", "list")` with three
elements:

- `diplotype_means`:

  Data frame. One row per diplotype class per LD block per trait, for
  diplotype classes with `>= min_n_diplotype` individuals. Contains 9
  columns:

  - `block_id` (character) - Block identifier.

  - `CHR` (character) - Chromosome.

  - `start_bp`, `end_bp` (integer) - Block coordinates.

  - `trait` (character) - Trait name.

  - `diplotype` (character) - Canonical diplotype string: two haplotype
    allele strings sorted alphabetically and joined by `"/"`, e.g.
    `"010/110"`. Homozygotes have identical halves: `"010/010"`. This
    format ensures that `"010/110"` and `"110/010"` are always
    represented identically.

  - `n_class` (integer) - Number of individuals carrying this specific
    diplotype combination.

  - `n_total` (integer) - Total individuals with non-missing data in
    this block (sum of n_class across all diplotype classes).

  - `mean_blue` (numeric) - Mean de-regressed phenotype value for this
    diplotype class (on the residual scale after GRM correction).

  - `se_mean` (numeric) - Standard error of the mean
    (`sd / sqrt(n_class)`).

- `dominance_table`:

  Data frame. One row per ordered allele pair (A, B) per block per
  trait, for pairs where all three diplotype classes (AA, AB, BB) each
  have `>= min_n_diplotype` individuals. Contains 14 columns:

  - `block_id`, `CHR`, `start_bp`, `end_bp`, `trait` - as above.

  - `allele_A`, `allele_B` (character) - The two alleles being compared
    (alphabetically ordered: A comes before B).

  - `mean_AA`, `mean_AB`, `mean_BB` (numeric) - Diplotype class means on
    the de-regressed phenotype scale.

  - `a` (numeric) - Additive effect: \\(\bar{y}\_{BB} - \bar{y}\_{AA}) /
    2\\. Positive means allele B increases trait value relative to A.

  - `d` (numeric) - Dominance deviation: \\\bar{y}\_{AB} -
    (\bar{y}\_{AA} + \bar{y}\_{BB}) / 2\\. Positive: heterozygote
    advantage (overdominance when \\\|d\| \> \|a\|\\). Negative:
    heterozygote disadvantage (underdominance).

  - `d_over_a` (numeric or `NA`) - Dominance ratio \\d/a\\.
    Interpretation: 0 = purely additive gene action; \\\pm 0.5\\ =
    partial dominance; \\\pm 1\\ = complete dominance (heterozygote
    equal to one homozygote); \\\|d/a\| \> 1\\ = overdominance
    (heterozygote exceeds both homozygotes - evidence for heterosis at
    this locus). `NA` when \\\|a\| \< 10^{-10}\\ (no additive variation
    between homozygotes).

  - `overdominance` (logical) - `TRUE` when \\\|d/a\| \> 1\\ and
    `d_over_a` is not `NA`.

- `omnibus_tests`:

  Data frame. One row per LD block per trait, for blocks where at least
  `min_n_diplotype` individuals carry at least 2 distinct diplotype
  classes. Contains 9 columns:

  - `block_id`, `trait` - as above.

  - `n_diplotypes` (integer) - Number of diplotype classes with
    `>= min_n_diplotype` individuals, used as factor levels.

  - `F_stat` (numeric) - F-statistic from one-way ANOVA (diplotype
    factor) on de-regressed phenotype residuals.

  - `df1` (integer) - Numerator df = `n_diplotypes - 1`.

  - `df2` (integer) - Denominator df = `n_obs - n_diplotypes`.

  - `p_omnibus` (numeric, (0,1\]) - Raw p-value.

  - `p_omnibus_adj` (numeric, (0,1\]) - Bonferroni-adjusted across all
    tested blocks per trait.

  - `significant` (logical) - `TRUE` when `p_omnibus_adj < 0.05`.

  Sorted ascending by `p_omnibus`.

## See also

[`infer_block_haplotypes`](https://FAkohoue.github.io/LDxBlocks/reference/infer_block_haplotypes.md),
[`test_block_haplotypes`](https://FAkohoue.github.io/LDxBlocks/reference/test_block_haplotypes.md)

## Examples

``` r
# \donttest{
data(ldx_geno, ldx_snp_info, ldx_blocks, ldx_blues, package = "LDxBlocks")
haps <- extract_haplotypes(ldx_geno, ldx_snp_info, ldx_blocks, min_snps = 5)
dip_res <- estimate_diplotype_effects(
  haplotypes = haps,
  blues      = setNames(ldx_blues$YLD, ldx_blues$id),
  blocks     = ldx_blocks,
  verbose    = FALSE
)
# Overdominant blocks
dip_res$dominance_table[dip_res$dominance_table$overdominance, ]
#>                block_id CHR start_bp end_bp trait                  allele_A
#> 1 block_1_155368_179371   1   155368 179371 trait 0100010000000000001000100
#> 2 block_2_161515_180473   2   161515 180473 trait      00100000100000000000
#>                    allele_B   mean_AA   mean_AB   mean_BB         a         d
#> 1 0111011011101110001000101 -0.438808  0.135218 -0.125740  0.156534  0.417491
#> 2      00100011100000110110  0.695056 -0.011832  0.214904 -0.240076 -0.466812
#>   d_over_a overdominance
#> 1   2.6671          TRUE
#> 2   1.9444          TRUE
# Significant diplotype effects
dip_res$omnibus_tests[dip_res$omnibus_tests$significant, ]
#> [1] block_id      trait         n_diplotypes  F_stat        df1          
#> [6] df2           p_omnibus     p_omnibus_adj significant  
#> <0 rows> (or 0-length row.names)
# }
```
