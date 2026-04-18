# Compare Haplotype Allele Frequencies Between Two Population Groups

For each LD block, computes allele frequencies in two sample groups and
returns FST (Weir-Cockerham 1984), frequency differences, and a
chi-squared test of independence. Useful for detecting blocks under
divergent selection or monitoring diversity changes between breeding
cycles.

## Usage

``` r
compare_haplotype_populations(
  haplotypes,
  group1,
  group2,
  group1_name = "group1",
  group2_name = "group2",
  min_freq = 0.02,
  missing_string = "."
)
```

## Arguments

- haplotypes:

  Named list from
  [`extract_haplotypes`](https://FAkohoue.github.io/LDxBlocks/reference/extract_haplotypes.md).

- group1:

  Character vector of individual IDs for group 1 (e.g. wild/landrace
  accessions).

- group2:

  Character vector of individual IDs for group 2 (e.g. elite breeding
  lines).

- group1_name:

  Character. Label for group 1. Default `"group1"`.

- group2_name:

  Character. Label for group 2. Default `"group2"`.

- min_freq:

  Numeric. Alleles below this frequency in both groups are pooled into
  an "other" category before testing. Default `0.02`.

- missing_string:

  Character. Missing haplotype placeholder. Default `"."`.

## Value

Data frame with one row per block, sorted by `CHR` and `start_bp`.

- `block_id`, `CHR`, `start_bp`, `end_bp`:

  Block coordinates.

- `n1`, `n2`:

  Sample sizes for group 1 and group 2.

- `n_alleles`:

  Number of distinct alleles in this block.

- `FST`:

  Weir-Cockerham FST, clamped to \[0,1\].

- `max_freq_diff`:

  Maximum absolute allele frequency difference.

- `dominant_g1`, `dominant_g2`:

  Most frequent allele in each group.

- `chisq_p`:

  Chi-squared p-value (Monte Carlo). `NA` if \< 2 alleles.

- `divergent`:

  Logical; TRUE when FST \> 0.1 and chisq_p \< 0.05.

## References

Weir BS, Cockerham CC (1984). Estimating F-statistics for the analysis
of population structure. *Evolution* **38**(6):1358-1370.

## See also

[`extract_haplotypes`](https://FAkohoue.github.io/LDxBlocks/reference/extract_haplotypes.md),
[`compute_haplotype_diversity`](https://FAkohoue.github.io/LDxBlocks/reference/compute_haplotype_diversity.md)

## Examples

``` r
# \donttest{
data(ldx_geno, ldx_snp_info, ldx_blocks, package = "LDxBlocks")
haps <- extract_haplotypes(ldx_geno, ldx_snp_info, ldx_blocks)
ids  <- rownames(ldx_geno)
cmp  <- compare_haplotype_populations(
  haplotypes   = haps,
  group1       = ids[1:60],
  group2       = ids[61:120],
  group1_name  = "cycle1",
  group2_name  = "cycle2"
)
cmp[cmp$divergent, c("block_id", "FST", "max_freq_diff", "chisq_p")]
#> [1] block_id      FST           max_freq_diff chisq_p      
#> <0 rows> (or 0-length row.names)
# }
```
