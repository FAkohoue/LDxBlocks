# Compute the Inverse Square Root (Whitening Factor) of a Kinship Matrix

Returns \\A\\ such that \\A V A^\top = I\\. Used by
[`prepare_geno`](https://FAkohoue.github.io/LDxBlocks/reference/prepare_geno.md)
when `method = "rV2"`.

## Usage

``` r
get_V_inv_sqrt(V, method = c("chol", "eigen"))
```

## Arguments

- V:

  Symmetric positive-definite matrix (n × n). Typically a VanRaden GRM
  after bending/tuning.

- method:

  Character. `"chol"` (default, faster) or `"eigen"` (more robust for
  near-singular matrices).

## Value

Numeric matrix A (n × n).
