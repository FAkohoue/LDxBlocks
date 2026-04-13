# Build Haplotype Dosage Matrix for Genomic Prediction

Converts haplotype strings to a numeric matrix for genomic prediction.
Supports phased and unphased input with two encoding schemes.

`encoding="additive_012"` (default, recommended for GBLUP/rrBLUP/BGLR):
Phased: 0=0 copies, 1=1 copy (het), 2=2 copies (hom) Unphased: 0=no
match, 2=match (1 not identifiable without phase)

encoding="presence_02" (kernel methods, random forest): Phased: 2=either
gamete matches, 0=neither, NA=missing Unphased: 2=match, 0=no match,
NA=missing

## Usage

``` r
build_haplotype_feature_matrix(
  haplotypes,
  top_n = NULL,
  encoding = c("additive_012", "presence_01"),
  missing_string = ".",
  scale_features = FALSE,
  min_freq = 0.01
)
```

## Arguments

- haplotypes:

  List from extract_haplotypes().

- top_n:

  Integer or `NULL`. Maximum number of haplotype alleles to retain per
  block, ranked by frequency. `NULL` (default) retains all alleles that
  pass `min_freq` – recommended for most analyses. Set an integer cap
  (e.g. `top_n = 5L`) only when you need to limit matrix width for
  memory reasons on panels with thousands of blocks and highly diverse
  haplotypes (many rare alleles above `min_freq`).

- encoding:

  Dosage encoding for the feature matrix:

  - `"additive_012"` (default): values are 0, 1, or 2 for **phased**
    data (0 = neither gamete, 1 = one gamete, 2 = both gametes carry
    this allele); values are 0 or 2 for **unphased** data (0 = absent, 2
    = homozygous for this allele). For **unphased** data the value is 0
    or 1/NA: 1 = the individual's haplotype string matches this allele
    exactly, 0 = it does not. The value 2 is not used for unphased data
    because we cannot confirm both chromosomes carry the same allele
    from unphased dosage strings. Compatible with rrBLUP, BGLR, sommer,
    ASReml-R.

  - `"presence_01"`: values are 0 or 1 — clean presence/absence
    encoding. For phased data: 1 if either gamete carries the allele.
    For unphased data: 1 if the individual's allele string matches.
    Loses copy-number information compared to `"additive_012"` but may
    be preferable for Bayesian variable selection models (BayesB,
    BayesC) where the prior expects binary indicators.

- missing_string:

  Missing data marker. Default ".".

- scale_features:

  Center and scale columns. Default FALSE.

- min_freq:

  Minimum allele frequency to include. Default 0.01.

## Value

Numeric matrix (individuals x haplotype allele columns).

## References

Difabachew YF et al. (2023). Genomic prediction with haplotype blocks in
wheat. *Frontiers in Plant Science* **14**:1168547.
[doi:10.3389/fpls.2023.1168547](https://doi.org/10.3389/fpls.2023.1168547)

Weber SE, Frisch M, Snowdon RJ, Voss-Fels KP (2023). Haplotype blocks
for genomic prediction: a comparative evaluation in multiple crop
datasets. *Frontiers in Plant Science* **14**:1217589.
[doi:10.3389/fpls.2023.1217589](https://doi.org/10.3389/fpls.2023.1217589)
