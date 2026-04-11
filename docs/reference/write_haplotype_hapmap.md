# Write Haplotype Matrix in HapMap Format

Rows=haplotype alleles, columns=individuals. Encoding: 0-\>HH, 1-\>HA,
2-\>AA, NA-\>NN (H=ref token, A=alt token). Compatible with TASSEL and
GAPIT.

## Usage

``` r
write_haplotype_hapmap(
  hap_matrix,
  out_file,
  ref_token = "H",
  alt_token = "A",
  verbose = TRUE
)
```

## Arguments

- hap_matrix:

  Numeric matrix (individuals x haplotype alleles).

- out_file:

  Output path (.hmp.txt recommended).

- ref_token:

  Reference allele symbol. Default "H".

- alt_token:

  Alternate allele symbol. Default "A".

- verbose:

  Logical. Default TRUE.

## Value

Invisibly returns out_file.
