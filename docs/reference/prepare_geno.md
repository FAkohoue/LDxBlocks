# Prepare Genotype Matrix for LD Computation

Central preparation function called inside
[`Big_LD`](https://FAkohoue.github.io/LDxBlocks/reference/Big_LD.md) and
[`run_Big_LD_all_chr`](https://FAkohoue.github.io/LDxBlocks/reference/run_Big_LD_all_chr.md).
Depending on `method`:

- `"r2"`:

  Mean-centres each SNP column. No kinship needed. Returns the centred
  matrix directly. Fast - O(np).

- `"rV2"`:

  Computes VanRaden GRM via AGHmatrix, tunes it via ASRgenomics, inverts
  via Cholesky/eigen, and left-multiplies the centred genotype matrix.
  Returns the whitened matrix and the whitening factor for use in
  `appendSGTs`.

## Usage

``` r
prepare_geno(
  geno,
  method = c("r2", "rV2"),
  kin_method = c("chol", "eigen"),
  verbose = FALSE
)
```

## Arguments

- geno:

  Numeric matrix (individuals x SNPs), 0/1/2.

- method:

  Character. `"r2"` or `"rV2"`.

- kin_method:

  Character. Whitening decomposition for rV^2: `"chol"` or `"eigen"`.

- verbose:

  Logical. Print progress.

## Value

Named list:

- adj_geno:

  n x p numeric matrix ready for
  [`compute_ld()`](https://FAkohoue.github.io/LDxBlocks/reference/compute_ld.md).

- V_inv_sqrt:

  n x n whitening matrix, or `NULL` for r^2.
