# Compute LD Matrix: Standard r^2 or Kinship-Adjusted rV^2

Unified dispatcher. With `method = "r2"` (default and recommended for
large datasets) it calls the C++ Armadillo back-end for fast standard
squared Pearson correlations. With `method = "rV2"` it applies kinship
whitening first and then the same C++ correlation kernel.

## Usage

``` r
compute_ld(X, method = c("r2", "rV2"), digits = -1L, n_threads = 1L)
```

## Arguments

- X:

  NumericMatrix (individuals x SNPs). For `method = "r2"`: raw or
  mean-centred genotypes (0/1/2). For `method = "rV2"`: the pre-whitened
  matrix \\V^{-1/2} \tilde{G}\\ produced by
  [`prepare_geno()`](https://FAkohoue.github.io/LDxBlocks/reference/prepare_geno.md).

- method:

  Character. `"r2"` (default) or `"rV2"`.

- digits:

  Integer. Rounding precision. `-1` skips rounding.

- n_threads:

  Integer. OpenMP threads for the C++ kernel. Default 1. Use
  [`parallel::detectCores()`](https://rdrr.io/r/parallel/detectCores.html)
  to choose automatically.

## Value

Symmetric p x p numeric matrix, diagonal 0, values in \[0, 1\].

## See also

[`compute_r2`](https://FAkohoue.github.io/LDxBlocks/reference/compute_r2.md),
[`compute_rV2`](https://FAkohoue.github.io/LDxBlocks/reference/compute_rV2.md),
[`prepare_geno`](https://FAkohoue.github.io/LDxBlocks/reference/prepare_geno.md),
[`Big_LD`](https://FAkohoue.github.io/LDxBlocks/reference/Big_LD.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# Internal function — use compute_r2() or compute_rV2() instead
set.seed(1)
G <- matrix(sample(0:2, 60 * 20, replace = TRUE), 60, 20)
ld_r2  <- LDxBlocks:::compute_ld(G, method = "r2")
range(ld_r2)
} # }
```
