# Build a Haplotype Dosage Feature Matrix for Genomic Prediction

Converts the phase-free haplotype strings from
[`extract_haplotypes`](https://FAkohoue.github.io/LDxBlocks/reference/extract_haplotypes.md)
into a numeric dosage matrix suitable for genomic prediction (GBLUP,
BayesB, random forests, etc.).

For each block, the `top_n` most frequent haplotypes are selected as
\*reference haplotypes\*. Each individual then receives an integer
dosage count (0, 1, or 2 for diploids) for each reference haplotype,
analogous to allele dosage in single-SNP analysis. Rare haplotypes (not
in the top `top_n`) are collapsed into an `"OTHER"` category.

## Usage

``` r
build_haplotype_feature_matrix(
  haplotypes,
  top_n = 5L,
  missing_string = ".",
  scale_features = FALSE
)
```

## Arguments

- haplotypes:

  List as returned by
  [`extract_haplotypes`](https://FAkohoue.github.io/LDxBlocks/reference/extract_haplotypes.md).

- top_n:

  Integer. Number of most frequent haplotypes per block to include as
  features. Remaining haplotypes are pooled as `"hap_OTHER"`. Default
  `5`.

- missing_string:

  Character. Haplotype strings matching this pattern are coded as `NA`
  in the output matrix. Default `"."`.

- scale_features:

  Logical. If `TRUE`, each haplotype dosage column is centred and scaled
  (mean 0, sd 1). Recommended when combining blocks in a single model.
  Default `FALSE`.

## Value

A numeric matrix of dimension `n_individuals × n_features`, where
features are haplotype dosages named `"<block_id>__hap<rank>"` (rank 1 =
most frequent). Row names are individual IDs (from the input haplotype
list). Columns for the `"OTHER"` category are omitted (it is the
implicit reference level).

## Why haplotype features?

Single-SNP models assume additive effects and miss multi-locus
interactions implicit in haplotype structure. Haplotype dosages capture
these interactions without requiring explicit interaction terms, often
improving prediction accuracy in structured populations (Calus et al.
2008; de Roos et al. 2009).

## References

Calus MPL et al. (2008) Accuracy of genomic selection using different
methods to define haplotypes. *Genetics* **178**(1):553–561.  
de Roos APW et al. (2009) Linkage disequilibrium and persistence of
phase in Holstein–Friesian, Jersey and Angus cattle. *Genetics*
**179**(3):1503–1512.

## See also

[`extract_haplotypes`](https://FAkohoue.github.io/LDxBlocks/reference/extract_haplotypes.md),
[`compute_haplotype_diversity`](https://FAkohoue.github.io/LDxBlocks/reference/compute_haplotype_diversity.md)

## Examples

``` r
# \donttest{
set.seed(1)
geno <- matrix(sample(0:2, 50 * 40, replace = TRUE), 50, 40)
rownames(geno) <- paste0("ind", 1:50)
colnames(geno) <- paste0("rs", 1:40)
snp_info <- data.frame(SNP = colnames(geno), CHR = "chr1",
  POS = seq(1000, by = 5000, length.out = 40))
blocks <- data.frame(start.bp = c(1000, 80000), end.bp = c(70000, 200000), CHR = "chr1")
haps <- extract_haplotypes(geno, snp_info, blocks)
feat <- build_haplotype_feature_matrix(haps, top_n = 3)
dim(feat)
#> [1] 50  6
feat[1:5, 1:min(6, ncol(feat))]
#>      block_1000_70000__hap1 block_1000_70000__hap2 block_1000_70000__hap3
#> ind1                      0                      0                      0
#> ind2                      0                      0                      0
#> ind3                      0                      0                      0
#> ind4                      0                      0                      0
#> ind5                      0                      0                      0
#>      block_80000_2e+05__hap1 block_80000_2e+05__hap2 block_80000_2e+05__hap3
#> ind1                       0                       0                       0
#> ind2                       0                       0                       0
#> ind3                       2                       0                       0
#> ind4                       0                       0                       0
#> ind5                       0                       0                       0
# }
```
