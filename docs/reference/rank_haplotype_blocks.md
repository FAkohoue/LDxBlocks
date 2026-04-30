# Rank Haplotype Blocks by Evidence Strength

Provides a unified block ranking that works across three use cases,
depending on what data the user has available. In all cases, `min_freq`
filtering inside
[`build_haplotype_feature_matrix`](https://FAkohoue.github.io/LDxBlocks/reference/build_haplotype_feature_matrix.md)
is applied first as a hard population-level filter – equivalent to MAF
filtering for single SNPs. Blocks are ranked only among those that
survive this filter.

## Usage

``` r
rank_haplotype_blocks(
  diversity,
  qtl_regions = NULL,
  pred_result = NULL,
  He_threshold = 0.3,
  top_n_blocks = NULL
)
```

## Arguments

- diversity:

  Data frame from
  [`compute_haplotype_diversity`](https://FAkohoue.github.io/LDxBlocks/reference/compute_haplotype_diversity.md).
  Required for all three use cases.

- qtl_regions:

  Optional. Data frame from
  [`define_qtl_regions`](https://FAkohoue.github.io/LDxBlocks/reference/define_qtl_regions.md).
  When supplied, blocks are binary-flagged as containing a GWAS hit or
  not. The p-value is not used for ranking. Default `NULL`.

- pred_result:

  Optional. List from
  [`run_haplotype_prediction`](https://FAkohoue.github.io/LDxBlocks/reference/run_haplotype_prediction.md).
  When supplied, blocks are ranked by `Var(local GEBV)`. Default `NULL`.

- He_threshold:

  Minimum He to consider a block diverse enough for haplotype stacking.
  Default `0.3`.

- top_n_blocks:

  Return only the top n blocks. Default `NULL` (return all).

## Value

Data frame with one row per block, sorted by evidence strength, with
columns: `block_id`, `CHR`, `start_bp`, `end_bp`, `n_snps`, `He`,
`n_eff_alleles`, `freq_dominant`, `sweep_flag`, `is_diverse`,
`has_gwas_hit` (if `qtl_regions` supplied), `lead_marker`, `lead_beta`,
`n_sig_markers` (if `qtl_regions` supplied), `var_scaled`,
`is_important` (if `pred_result` supplied), `use_case`, `rank_score`,
`recommendation`.

## The three use cases

1.  **Genotype only** (no GWAS, no phenotype): blocks ranked by
    haplotype diversity (He, effective number of alleles). Blocks with
    high diversity have the most potential for haplotype stacking.

2.  **Genotype + GWAS** (no phenotype): blocks are binary-flagged by
    whether they contain a GWAS-significant marker. Within the GWAS-hit
    group and within the non-hit group, blocks are further ordered by
    He. The p-value is not used for ranking – a marker either crosses
    the significance threshold or it does not.

3.  **Genotype + phenotype** (? GWAS): blocks ranked by scaled
    `Var(local GEBV)` from
    [`run_haplotype_prediction`](https://FAkohoue.github.io/LDxBlocks/reference/run_haplotype_prediction.md).
    When GWAS results are also available, the binary GWAS flag is added
    as a secondary layer via
    [`integrate_gwas_haplotypes`](https://FAkohoue.github.io/LDxBlocks/reference/integrate_gwas_haplotypes.md).

## Link to `min_freq`

`min_freq` is a hard population-level pre-filter identical in purpose to
MAF filtering for single SNPs: haplotype alleles observed at frequency
below this threshold cannot have their effects reliably estimated
regardless of trait association, and are dropped before the dosage
matrix is built. `rank_haplotype_blocks` operates entirely downstream of
this filter – it ranks blocks that have already passed `min_freq`, not
individual alleles. A block survives as long as at least one of its
alleles passes `min_freq`.

## References

Difabachew YF et al. (2023). Genomic prediction with haplotype blocks in
wheat. *Frontiers in Plant Science* **14**:1168547.
[doi:10.3389/fpls.2023.1168547](https://doi.org/10.3389/fpls.2023.1168547)

Weber SE, Frisch M, Snowdon RJ, Voss-Fels KP (2023). Haplotype blocks
for genomic prediction: a comparative evaluation in multiple crop
datasets. *Frontiers in Plant Science* **14**:1217589.
[doi:10.3389/fpls.2023.1217589](https://doi.org/10.3389/fpls.2023.1217589)

Tong J et al. (2024). Stacking beneficial haplotypes from the Vavilov
wheat collection. *Theoretical and Applied Genetics* **137**:274.
[doi:10.1007/s00122-024-04784-w](https://doi.org/10.1007/s00122-024-04784-w)

Tong J et al. (2025). Haplotype stacking to improve stability of stripe
rust resistance in wheat. *Theoretical and Applied Genetics*
**138**:267.
[doi:10.1007/s00122-025-05045-0](https://doi.org/10.1007/s00122-025-05045-0)

## Examples

``` r
if (FALSE) { # \dontrun{
haps <- extract_haplotypes(geno, snp_info, blocks)
div  <- compute_haplotype_diversity(haps)

# Use case 1: genotype only -- rank by diversity
ranked <- rank_haplotype_blocks(div)

# Use case 2: genotype + GWAS -- binary flag, then diversity within groups
qtl    <- define_qtl_regions(gwas_df, blocks, snp_info)
ranked <- rank_haplotype_blocks(div, qtl_regions = qtl)

# Use case 3: genotype + phenotype -- rank by Var(local GEBV)
pred   <- run_haplotype_prediction(geno, snp_info, blocks, blues = blues)
ranked <- rank_haplotype_blocks(div, qtl_regions = qtl, pred_result = pred)

# Top 20 blocks
ranked <- rank_haplotype_blocks(div, pred_result = pred, top_n_blocks = 20)
} # }
```
