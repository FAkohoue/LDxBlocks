# Build Haplotype Dosage Matrix for Genomic Prediction

Converts haplotype strings to a numeric matrix for genomic prediction.
Supports phased and unphased input with two encoding schemes.

encoding="additive_012" (default, recommended for GBLUP/rrBLUP/BGLR):
Phased: 0=0 copies, 1=1 copy (het), 2=2 copies (hom) Unphased: 0=no
match, 2=match (1 not identifiable without phase)

encoding="presence_02" (kernel methods, random forest): Phased: 2=either
gamete matches, 0=neither, NA=missing Unphased: 2=match, 0=no match,
NA=missing

## Usage

``` r
build_haplotype_feature_matrix(
  haplotypes,
  top_n = 5L,
  encoding = c("additive_012", "presence_02"),
  missing_string = ".",
  scale_features = FALSE,
  min_freq = 0.01
)
```

## Arguments

- haplotypes:

  List from extract_haplotypes().

- top_n:

  Top haplotype alleles per block. Default 5L.

- encoding:

  "additive_012" (default) or "presence_02".

- missing_string:

  Missing data marker. Default ".".

- scale_features:

  Center and scale columns. Default FALSE.

- min_freq:

  Minimum allele frequency to include. Default 0.01.

## Value

Numeric matrix (individuals x haplotype allele columns).
