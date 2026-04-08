# Compute Standard r² LD Matrix

Fast C++/Armadillo implementation of the standard pairwise squared
Pearson correlation (r²) for a window of SNP columns. Missing genotypes
(NA) are mean-imputed per column before computation. This is the default
LD metric in LDxBlocks and is 10–50× faster than
[`stats::cor()`](https://rdrr.io/r/stats/cor.html).

## Usage

``` r
compute_r2(X, digits = -1L, n_threads = 1L)
```

## Arguments

- X:

  NumericMatrix (individuals × SNPs). Values 0/1/2. NA allowed.

- digits:

  Integer. Rounding decimal places. `-1` (default) = none.

- n_threads:

  Integer. OpenMP threads. Default 1.

## Value

Symmetric p × p NumericMatrix, diagonal 0, values in \[0, 1\].

## When to use r² vs rV²

- r²:

  Use for large unstructured datasets (\> 500 k markers), random mating
  populations, or whenever speed matters. The standard estimator is
  inflated in related populations (i.e. will over-estimate LD) but this
  usually leads to slightly more conservative (larger) blocks rather
  than catastrophically wrong ones.

- rV²:

  Use for highly structured / related populations (livestock, inbred
  lines, family-based human cohorts) where kinship inflation would
  meaningfully distort block boundaries. Requires computing and
  inverting the GRM — prohibitive beyond ~5 k individuals.

## See also

[`compute_rV2`](https://FAkohoue.github.io/LDxBlocks/reference/compute_rV2.md),
[`compute_ld`](https://FAkohoue.github.io/LDxBlocks/reference/compute_ld.md)

## Examples

``` r
set.seed(1)
G <- matrix(sample(0:2, 80 * 30, replace = TRUE), 80, 30)
r2 <- compute_r2(G, digits = 6L)
range(r2)   # [0, 1]
#> [1] 0.000000 0.117788
```
