# Collapse Rare Haplotype Alleles Into Biologically Meaningful Groups

Rather than dropping rare haplotype alleles below a frequency threshold
(as
[`build_haplotype_feature_matrix`](https://FAkohoue.github.io/LDxBlocks/reference/build_haplotype_feature_matrix.md)
does), this function merges them into the most appropriate existing
allele. Three strategies are supported:

- `"rare_to_other"`:

  All alleles below `min_freq` are pooled into a single catch-all
  `"<other>"` category. Simple and lossless for total frequency; best
  when rare alleles are truly heterogeneous.

- `"nearest"`:

  Each rare allele is merged with the most similar common allele
  (minimum Hamming distance). Preserves biological similarity; best when
  rare alleles are likely sequencing errors or very recent recombinants
  of existing haplotypes.

- `"tree_based"`:

  Builds a UPGMA dendrogram from pairwise Hamming distances and merges
  rare alleles by cutting the tree at the level that produces the fewest
  groups while keeping all groups above `min_freq`. Computationally
  heavier but most biologically principled.

## Usage

``` r
collapse_haplotypes(
  haplotypes,
  min_freq = 0.05,
  collapse = c("nearest", "rare_to_other", "tree_based"),
  missing_string = ".",
  keep_labels = TRUE
)
```

## Arguments

- haplotypes:

  Named list from
  [`extract_haplotypes`](https://FAkohoue.github.io/LDxBlocks/reference/extract_haplotypes.md).

- min_freq:

  Numeric. Alleles at or below this frequency are considered rare.
  Default `0.05`.

- collapse:

  Character. Collapsing strategy: `"rare_to_other"`, `"nearest"`, or
  `"tree_based"`. Default `"nearest"`.

- missing_string:

  Character. Missing data placeholder. Default `"."`.

- keep_labels:

  Logical. If `TRUE` (default), a `"label_map"` attribute is attached to
  the output list: a named list per block giving the original →
  collapsed allele mapping.

## Value

Named list of the same structure as the input `haplotypes`, with rare
allele strings replaced by their collapsed equivalents. The `block_info`
attribute is preserved. If `keep_labels = TRUE`, a `label_map` attribute
is also attached (used by
[`harmonize_haplotypes`](https://FAkohoue.github.io/LDxBlocks/reference/harmonize_haplotypes.md)).

## See also

[`extract_haplotypes`](https://FAkohoue.github.io/LDxBlocks/reference/extract_haplotypes.md),
[`harmonize_haplotypes`](https://FAkohoue.github.io/LDxBlocks/reference/harmonize_haplotypes.md),
[`build_haplotype_feature_matrix`](https://FAkohoue.github.io/LDxBlocks/reference/build_haplotype_feature_matrix.md)

## Examples

``` r
# \donttest{
data(ldx_geno, ldx_snp_info, ldx_blocks, package = "LDxBlocks")
haps      <- extract_haplotypes(ldx_geno, ldx_snp_info, ldx_blocks)
haps_col  <- collapse_haplotypes(haps, min_freq = 0.05, collapse = "nearest")
# Compare diversity before/after
div_before <- compute_haplotype_diversity(haps)
div_after  <- compute_haplotype_diversity(haps_col)
summary(div_after$He - div_before$He)  # small reduction expected
#>     Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
#> -0.07549 -0.06681 -0.04398 -0.05121 -0.03950 -0.03319 
# }
```
