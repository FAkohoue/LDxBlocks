# Write Haplotype Matrix as Numeric CSV

Rows=individuals, columns=haplotype alleles (0/1/2 or 0/2). Compatible
with rrBLUP, BGLR, ASReml-R.

## Usage

``` r
write_haplotype_numeric(
  hap_matrix,
  out_file,
  sep = ",",
  na_str = "NA",
  verbose = TRUE
)
```

## Arguments

- hap_matrix:

  Numeric matrix (individuals x haplotype alleles).

- out_file:

  Output CSV path.

- sep:

  Separator. Default ",".

- na_str:

  NA string. Default "NA".

- verbose:

  Logical. Default TRUE.

## Value

Invisibly returns out_file.
