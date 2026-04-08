# Extract a Genotype Slice from an LDxBlocks Backend

Returns an `n_samples × length(col_idx)` numeric matrix of dosage values
(0/1/2/NA) for the selected SNP columns. Works identically for all
supported backend types.

## Usage

``` r
read_chunk(backend, col_idx)
```

## Arguments

- backend:

  An object of class `"LDxBlocks_backend"` as returned by
  [`read_geno`](https://FAkohoue.github.io/LDxBlocks/reference/read_geno.md).

- col_idx:

  Integer vector of column indices (1-based SNP positions).

## Value

Numeric matrix (n_samples × length(col_idx)).
