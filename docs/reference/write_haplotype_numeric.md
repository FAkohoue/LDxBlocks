# Write Haplotype Feature Matrix as Numeric Dosage Table

Writes the haplotype dosage matrix in a tab-delimited format with
haplotype alleles as rows and individuals as columns. Metadata columns
(`hap_id`, `CHR`, `start_bp`, `end_bp`, `n_snps`, `alleles`,
`frequency`) precede the individual columns. Individual cells contain
0/1/2/NA dosage values.

## Usage

``` r
write_haplotype_numeric(
  hap_matrix,
  out_file,
  haplotypes = NULL,
  snp_info = NULL,
  hap_info = NULL,
  sep = "\t",
  na_str = "NA",
  min_freq = 0.01,
  missing_string = ".",
  verbose = TRUE
)
```

## Arguments

- hap_matrix:

  Numeric matrix (individuals x haplotype alleles) from
  [`build_haplotype_feature_matrix`](https://FAkohoue.github.io/LDxBlocks/reference/build_haplotype_feature_matrix.md).

- out_file:

  Output file path.

- haplotypes:

  List from
  [`extract_haplotypes`](https://FAkohoue.github.io/LDxBlocks/reference/extract_haplotypes.md).
  When supplied together with `snp_info`, the `alleles` and `frequency`
  metadata columns are populated.

- snp_info:

  Data frame with `CHR`, `POS`, `REF`, `ALT`. Required for `alleles`
  column.

- hap_info:

  Data frame of exact per-column metadata from
  [`build_haplotype_feature_matrix`](https://FAkohoue.github.io/LDxBlocks/reference/build_haplotype_feature_matrix.md)`()\$info`.
  When supplied, the `alleles`, `frequency`, `n_snps`, `CHR`,
  `start_bp`, and `end_bp` columns are written directly from this object
  without any reconstruction. Recommended - pass
  `hap_info = feat_out\$info` where `feat_out` is the return value of
  [`build_haplotype_feature_matrix()`](https://FAkohoue.github.io/LDxBlocks/reference/build_haplotype_feature_matrix.md).
  Default `NULL` (falls back to legacy reconstruction from `haplotypes`
  and `snp_info`).

- sep:

  Field separator. Default `","`.

- na_str:

  NA string. Default `"NA"`.

- min_freq:

  Minimum frequency used when computing `alleles`. Default `0.01`.

- missing_string:

  Missing genotype marker. Default `"."`.

- verbose:

  Logical. Default `TRUE`.

## Value

Invisibly returns `out_file`.
