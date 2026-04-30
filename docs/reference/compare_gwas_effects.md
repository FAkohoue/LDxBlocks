# Cross-Population GWAS Effect Concordance from External Results

Compares haplotype block-level association signals between two
populations when GWAS was run externally (e.g. in GAPIT, TASSEL,
FarmCPU, PLINK, or any other tool) rather than through
[`test_block_haplotypes`](https://FAkohoue.github.io/LDxBlocks/reference/test_block_haplotypes.md).

The function accepts external GWAS results in two ways:

1.  **Pre-mapped** (recommended): supply the output of
    [`define_qtl_regions`](https://FAkohoue.github.io/LDxBlocks/reference/define_qtl_regions.md)
    for each population. The QTL table already maps each marker to its
    LD block, providing the lead SNP, its effect size (BETA), and its
    p-value per block. This is the most reliable input because the block
    assignment is explicit.

2.  **Raw GWAS + blocks**: supply raw GWAS data frames plus block
    tables; the function calls
    [`define_qtl_regions`](https://FAkohoue.github.io/LDxBlocks/reference/define_qtl_regions.md)
    internally to map markers to blocks before comparing.

**What is compared and how:**

External GWAS produces marker-level effects, not haplotype-allele-level
effects. At the block level each population contributes one observation:
the lead SNP (smallest p-value in the block). The comparison is
therefore analogous to two-sample Mendelian randomisation - one effect
estimate per exposure (block) per population, combined via IVW
meta-analysis.

Because there is only one "allele" per block (the lead SNP), several
statistics differ from
[`compare_block_effects`](https://FAkohoue.github.io/LDxBlocks/reference/compare_block_effects.md):

- `n_shared_alleles` is always 1 (the lead SNP concept).

- `effect_correlation` is always `NA` (needs \>= 3 alleles).

- `direction_agreement` is 0 or 1 (sign of lead SNP effect agrees or
  not).

- `Q_df` is always 0 (one observation per population; Cochran Q is
  undefined with two populations and one allele per block).

- `meta_p` tests whether the IVW combined effect differs from zero - the
  primary replication test.

**Standard error derivation:**

Many GWAS tools report BETA and P but not SE. When SE is absent,
`compare_gwas_effects()` derives it from the z-score: \$\$SE =
\frac{\|\beta\|}{\|z\|},\quad z = \Phi^{-1}(P/2)\$\$ where \\\Phi^{-1}\\
is the standard normal quantile function. This is exact for large-sample
Wald tests (the dominant GWAS framework) and provides a valid SE for IVW
weighting.

**Interpreting pleiotropic blocks:**

When `trait_col` is supplied and GWAS results contain multiple traits,
the function compares trait by trait. Blocks significant for multiple
traits in both populations (`pleiotropic = TRUE` in both QTL tables) are
flagged with `both_pleiotropic = TRUE` in the output, identifying the
most robust cross-population replication targets.

## Usage

``` r
compare_gwas_effects(
  qtl_pop1 = NULL,
  qtl_pop2 = NULL,
  gwas_pop1 = NULL,
  gwas_pop2 = NULL,
  blocks_pop1 = NULL,
  blocks_pop2 = NULL,
  snp_info_pop1 = NULL,
  snp_info_pop2 = NULL,
  pop1_name = "pop1",
  pop2_name = "pop2",
  beta_col = "BETA",
  se_col = "SE",
  p_col = "P",
  trait_col = "trait",
  p_threshold = 5e-08,
  min_snps = 3L,
  block_match = c("id", "position"),
  overlap_min = 0.5,
  direction_threshold = 0.75,
  boundary_overlap_warn = 0.8,
  verbose = TRUE
)
```

## Arguments

- qtl_pop1:

  Data frame. Output of
  [`define_qtl_regions`](https://FAkohoue.github.io/LDxBlocks/reference/define_qtl_regions.md)
  for population 1 (discovery). Required columns: `block_id`, `CHR`,
  `start_bp`, `end_bp`, `lead_marker`, `lead_p`, `lead_beta`. Optional:
  `traits`, `pleiotropic`. Ignored when `gwas_pop1` is supplied instead.

- qtl_pop2:

  Data frame. Output of
  [`define_qtl_regions`](https://FAkohoue.github.io/LDxBlocks/reference/define_qtl_regions.md)
  for population 2 (validation). Same structure as `qtl_pop1`. Ignored
  when `gwas_pop2` is supplied instead.

- gwas_pop1:

  Data frame. Raw GWAS results for population 1. Required when
  `qtl_pop1` is `NULL`. Must contain `SNP` (or `Marker`), `CHR`, `POS`,
  and the columns named by `beta_col`, `se_col` (optional), `p_col`.

- gwas_pop2:

  Data frame. Raw GWAS results for population 2. Same structure as
  `gwas_pop1`.

- blocks_pop1:

  Data frame. LD block table for population 1 (output of
  [`run_Big_LD_all_chr`](https://FAkohoue.github.io/LDxBlocks/reference/run_Big_LD_all_chr.md)).
  Required when `gwas_pop1` is supplied. Also used for
  `boundary_overlap_ratio` computation when `qtl_pop1` is the input.

- blocks_pop2:

  Data frame. LD block table for population 2. Same structure as
  `blocks_pop1`.

- snp_info_pop1:

  Data frame. SNP metadata for population 1 (columns: `SNP`, `CHR`,
  `POS`). Required when `gwas_pop1` is supplied.

- snp_info_pop2:

  Data frame. SNP metadata for population 2. Same structure as
  `snp_info_pop1`. If `NULL` (default), `snp_info_pop1` is reused
  (appropriate when both populations share the same marker panel).

- pop1_name:

  Character. Label for population 1. Default `"pop1"`.

- pop2_name:

  Character. Label for population 2. Default `"pop2"`.

- beta_col:

  Character. Name of the effect-size column in raw GWAS data frames.
  Default `"BETA"`.

- se_col:

  Character or `NULL`. Name of the standard-error column. When `NULL` or
  absent, SE is derived from `beta_col` and `p_col` via the z-score
  formula. Default `"SE"`.

- p_col:

  Character. Name of the p-value column. Default `"P"`.

- trait_col:

  Character. Name of the trait column when GWAS results contain multiple
  traits. Default `"trait"`.

- p_threshold:

  Numeric. Significance threshold applied when mapping raw GWAS results
  to blocks via
  [`define_qtl_regions`](https://FAkohoue.github.io/LDxBlocks/reference/define_qtl_regions.md).
  Ignored when `qtl_pop1` / `qtl_pop2` are supplied directly. Default
  `5e-8`.

- min_snps:

  Integer. Minimum number of SNPs in a block for it to be included in
  the QTL mapping step. Default `3L`.

- block_match:

  Character. How to match blocks between populations. `"id"` (default)
  matches by exact `block_id` string - fast and backward-compatible,
  appropriate when both populations were mapped against the same block
  table. `"position"` matches by genomic interval overlap
  (Intersection-over-Union \\\geq\\ `overlap_min`) - recommended when
  block boundaries differ between populations due to different ancestral
  LD structures or independent block-detection runs.

- overlap_min:

  Numeric in (0, 1\]. Minimum Intersection-over-Union (IoU) in base
  pairs for two blocks to be considered the same region when
  `block_match = "position"`. Default `0.50`. Ignored when
  `block_match = "id"`.

- direction_threshold:

  Numeric in (0.5, 1\]. Minimum direction agreement to consider a block
  directionally concordant. Default `0.75`. Note: with one lead SNP per
  block this is effectively a 0/1 flag at this threshold, but the
  parameter is kept for consistency with
  [`compare_block_effects`](https://FAkohoue.github.io/LDxBlocks/reference/compare_block_effects.md).

- boundary_overlap_warn:

  Numeric in (0, 1). `boundary_overlap_ratio` below this value triggers
  `boundary_warning = TRUE`. The `boundary_overlap_ratio` is an
  \*\*output column\*\* automatically computed from the block tables -
  it cannot be set by the user. Default `0.80`.

- verbose:

  Logical. Print progress. Default `TRUE`.

## Value

A named list of class `c("LDxBlocks_effect_concordance", "list")` with
the same structure as
[`compare_block_effects`](https://FAkohoue.github.io/LDxBlocks/reference/compare_block_effects.md).
Additional columns in `$concordance` specific to GWAS input:

- `lead_marker_pop1`, `lead_marker_pop2` - lead SNP ID from each
  population (same SNP = same tag; different SNP = different LD tag,
  same region).

- `lead_p_pop1`, `lead_p_pop2` - lead SNP p-values.

- `se_derived_pop1`, `se_derived_pop2` - logical; `TRUE` when SE was
  derived from BETA and P rather than read directly.

- `both_pleiotropic` - logical; `TRUE` when the block is pleiotropic in
  both populations (requires `pleiotropic` column in QTL tables).

The `$shared_alleles` data frame contains one row per block per trait
with `lead_marker` instead of `allele`, and `effect_pop1`, `SE_pop1`,
`effect_pop2`, `SE_pop2`, `direction_agree`, `ivw_effect`, `ivw_SE`.

## References

Borenstein M, Hedges LV, Higgins JPT, Rothstein HR (2009). *Introduction
to Meta-Analysis*. Wiley.

Higgins JPT, Thompson SG (2002). Quantifying heterogeneity in a
meta-analysis. *Statistics in Medicine* **21**(11):1539-1558.
[doi:10.1002/sim.1186](https://doi.org/10.1002/sim.1186)

## See also

[`define_qtl_regions`](https://FAkohoue.github.io/LDxBlocks/reference/define_qtl_regions.md),
[`compare_block_effects`](https://FAkohoue.github.io/LDxBlocks/reference/compare_block_effects.md),
[`test_block_haplotypes`](https://FAkohoue.github.io/LDxBlocks/reference/test_block_haplotypes.md)

## Examples

``` r
# \donttest{
data(ldx_geno, ldx_snp_info, ldx_blocks, ldx_gwas, package = "LDxBlocks")

# Simulate two populations
set.seed(1L)
n    <- nrow(ldx_geno)
idx1 <- sample(n, 70L)
idx2 <- setdiff(seq_len(n), idx1)

# Suppose we have external GWAS results for both populations
# (here we use ldx_gwas as a stand-in)
gwas_A <- ldx_gwas; gwas_A$BETA <- rnorm(nrow(gwas_A), 0.5, 0.1)
gwas_B <- ldx_gwas; gwas_B$BETA <- rnorm(nrow(gwas_B), 0.4, 0.15)

# Path 1: raw GWAS + blocks (convenience)
conc <- compare_gwas_effects(
  gwas_pop1     = gwas_A,
  gwas_pop2     = gwas_B,
  blocks_pop1   = ldx_blocks,
  blocks_pop2   = ldx_blocks,
  snp_info_pop1 = ldx_snp_info,
  pop1_name     = "PopA",
  pop2_name     = "PopB",
  p_threshold   = NULL   # keep all markers for this demo
)
#> [compare_gwas_effects] Mapping PopA GWAS results to blocks ...
#> [compare_gwas_effects] Mapping PopB GWAS results to blocks ...
#> [compare_gwas_effects] Comparing effects for blocks with hits in both populations ...
#> [compare_gwas_effects] Traits: trait
#> [compare_gwas_effects] Done. Blocks in both pops: 6 | Replicated (dir concordant + meta_p <= 0.05): 6
print(conc)
#> LDxBlocks Cross-Population Effect Concordance
#>   Populations: PopA vs PopB
#>   Traits:       trait 
#>   Blocks compared:           6 
#>   With enough shared alleles: 6 
#>   Directionally concordant:   6 
#>   Replicated (dir + Q_p>0.05): 6 
#>   Boundary warnings:          0 (overlap ratio < 0.8 )
#>   Shared allele comparisons:  6 

# Path 2: pre-mapped QTL tables (recommended)
qtl_A <- define_qtl_regions(gwas_A, ldx_blocks, ldx_snp_info,
                             p_threshold = NULL)
qtl_B <- define_qtl_regions(gwas_B, ldx_blocks, ldx_snp_info,
                             p_threshold = NULL)
conc2 <- compare_gwas_effects(
  qtl_pop1      = qtl_A,
  qtl_pop2      = qtl_B,
  blocks_pop1   = ldx_blocks,
  blocks_pop2   = ldx_blocks,
  pop1_name     = "PopA",
  pop2_name     = "PopB"
)
#> [compare_gwas_effects]   PopA: lead_se absent - deriving from lead_beta + lead_p
#> [compare_gwas_effects]   PopB: lead_se absent - deriving from lead_beta + lead_p
#> [compare_gwas_effects] Comparing effects for blocks with hits in both populations ...
#> [compare_gwas_effects] Traits: trait
#> [compare_gwas_effects] Done. Blocks in both pops: 6 | Replicated (dir concordant + meta_p <= 0.05): 6
conc2$concordance[conc2$concordance$replicated, ]
#>                block_id CHR start_bp end_bp trait lead_marker_pop1
#> 1    block_1_1000_25027   1     1000  25027 trait           rs1005
#> 2   block_1_81064_99022   1    81064  99022 trait           rs1048
#> 3 block_1_155368_179371   1   155368 179371 trait           rs1070
#> 4    block_2_1000_30023   2     1000  30023 trait           rs2004
#> 5  block_2_86236_105290   2    86236 105290 trait           rs2050
#> 6    block_3_1000_19068   3     1000  19068 trait           rs3004
#>   lead_marker_pop2  lead_p_pop1  lead_p_pop2 se_derived_pop1 se_derived_pop2
#> 1           rs1005 4.023280e-07 4.023280e-07            TRUE            TRUE
#> 2           rs1048 3.493586e-07 3.493586e-07            TRUE            TRUE
#> 3           rs1070 3.653110e-09 3.653110e-09            TRUE            TRUE
#> 4           rs2004 1.726081e-05 1.726081e-05            TRUE            TRUE
#> 5           rs2050 3.857636e-04 3.857636e-04            TRUE            TRUE
#> 6           rs3004 2.215581e-02 2.215581e-02            TRUE            TRUE
#>   n_alleles_pop1 n_alleles_pop2 n_shared_alleles enough_shared
#> 1              1              1                1          TRUE
#> 2              1              1                1          TRUE
#> 3              1              1                1          TRUE
#> 4              1              1                1          TRUE
#> 5              1              1                1          TRUE
#> 6              1              1                1          TRUE
#>   effect_correlation direction_agreement directionally_concordant meta_effect
#> 1                 NA                   1                     TRUE    0.511102
#> 2                 NA                   1                     TRUE    0.377802
#> 3                 NA                   1                     TRUE    0.462585
#> 4                 NA                   1                     TRUE    0.379837
#> 5                 NA                   1                     TRUE    0.453007
#> 6                 NA                   1                     TRUE    0.481742
#>    meta_SE meta_z       meta_p Q_stat Q_df Q_p I2 replicated both_pleiotropic
#> 1 0.071374 7.1609 8.016005e-13     NA   NA  NA NA       TRUE            FALSE
#> 2 0.053127 7.1112 1.150019e-12     NA   NA  NA NA       TRUE             TRUE
#> 3 0.056073 8.2496 1.588641e-16     NA   NA  NA NA       TRUE            FALSE
#> 4 0.067388 5.6366 1.734364e-08     NA   NA  NA NA       TRUE             TRUE
#> 5 0.090458 5.0079 5.502050e-07     NA   NA  NA NA       TRUE             TRUE
#> 6 0.149185 3.2292 1.241580e-03     NA   NA  NA NA       TRUE             TRUE
#>   boundary_overlap_ratio boundary_warning match_type
#> 1                      1            FALSE      exact
#> 2                      1            FALSE      exact
#> 3                      1            FALSE      exact
#> 4                      1            FALSE      exact
#> 5                      1            FALSE      exact
#> 6                      1            FALSE      exact
# }
```
