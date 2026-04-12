# Compute Kinship-Adjusted rV^2 LD Matrix

Computes the kinship-adjusted squared correlation (rV^2) for a
pre-whitened genotype matrix \\X = V^{-1/2} \tilde{G}\\. The whitening
step must be performed first via
[`prepare_geno`](https://FAkohoue.github.io/LDxBlocks/reference/prepare_geno.md).

Mathematically identical to
[`compute_r2`](https://FAkohoue.github.io/LDxBlocks/reference/compute_r2.md)
applied to \\X\\ because the whitening already removes kinship
structure. The C++ kernel is the same - the distinction is purely in the
preparation step.

## Usage

``` r
compute_rV2(X, digits = -1L, n_threads = 1L)
```

## Arguments

- X:

  NumericMatrix (n x p). Pre-whitened, mean-centred genotype matrix.

- digits:

  Integer. Rounding. `-1` = none.

- n_threads:

  Integer. OpenMP threads.

## Value

Symmetric p x p NumericMatrix, diagonal 0.

## References

Kim S-A et al. (2018) GENETICS 209(3):855-868.

## See also

[`compute_r2`](https://FAkohoue.github.io/LDxBlocks/reference/compute_r2.md),
[`prepare_geno`](https://FAkohoue.github.io/LDxBlocks/reference/prepare_geno.md)
