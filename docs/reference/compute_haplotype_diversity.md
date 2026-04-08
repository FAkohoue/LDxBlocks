# Compute Haplotype Diversity Metrics Per LD Block

For each block produced by
[`extract_haplotypes`](https://FAkohoue.github.io/LDxBlocks/reference/extract_haplotypes.md),
computes:

- Richness (\\k\\):

  Number of unique haplotype strings.

- Expected heterozygosity (\\H_e\\):

  Nei's (1973) gene diversity: \\H_e = \frac{n}{n-1}\left(1 - \sum_i
  p_i^2\right)\\ where \\p_i\\ is the frequency of haplotype \\i\\ and
  \\n\\ is the number of individuals (with non-missing haplotypes).

- Shannon entropy (\\H'\\):

  \\H' = -\sum_i p_i \log_2(p_i)\\. Sensitive to both richness and
  evenness.

- Dominant haplotype frequency (\\f\_{\max}\\):

  Frequency of the most common haplotype; high values indicate a
  selective sweep or strong founder effect.

## Usage

``` r
compute_haplotype_diversity(haplotypes, missing_string = ".")
```

## Arguments

- haplotypes:

  List as returned by
  [`extract_haplotypes`](https://FAkohoue.github.io/LDxBlocks/reference/extract_haplotypes.md).

- missing_string:

  Character. Haplotype strings containing this pattern are treated as
  missing and excluded from frequency calculations. Default `"."`.

## Value

A `data.frame` with one row per block and columns: `block_id`, `n_ind`
(non-missing individuals), `n_haplotypes` (richness), `He` (expected
heterozygosity), `Shannon` (entropy in bits), `freq_dominant` (frequency
of most common haplotype).

## References

Nei M (1973) Analysis of gene diversity in subdivided populations.
*PNAS* **70**(12):3321–3323.

## See also

[`extract_haplotypes`](https://FAkohoue.github.io/LDxBlocks/reference/extract_haplotypes.md),
[`build_haplotype_feature_matrix`](https://FAkohoue.github.io/LDxBlocks/reference/build_haplotype_feature_matrix.md)

## Examples

``` r
# \donttest{
# (Continuing extract_haplotypes example)
set.seed(1)
geno <- matrix(sample(0:2, 50 * 40, replace = TRUE), 50, 40)
rownames(geno) <- paste0("ind", 1:50)
colnames(geno) <- paste0("rs", 1:40)
snp_info <- data.frame(SNP = colnames(geno), CHR = "chr1",
  POS = seq(1000, by = 5000, length.out = 40))
blocks <- data.frame(start.bp = c(1000, 80000), end.bp = c(70000, 200000), CHR = "chr1")
haps <- extract_haplotypes(geno, snp_info, blocks)
diversity <- compute_haplotype_diversity(haps)
print(diversity)
#>            block_id n_ind n_haplotypes He  Shannon freq_dominant
#> 1  block_1000_70000    50           50  1 5.643856          0.02
#> 2 block_80000_2e+05    50           50  1 5.643856          0.02
# }
```
