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

  - `"additive_012"` (default): **Phased data**: values are 0 (neither
    gamete carries the allele), 1 (one gamete carries it -
    heterozygous), or 2 (both gametes - homozygous). This gives true
    allele dosage on the standard 0/1/2 scale. **Unphased data**: values
    are 0 or 1 only (presence/absence of the dosage-pattern haplotype).
    The value 2 is never produced for unphased data because it is
    impossible to confirm that both chromosomes carry the same haplotype
    string without phase information. Compatible with rrBLUP, BGLR,
    sommer, ASReml-R.

  - `"presence_01"`: values are 0 or 1 for both phased and unphased
    data. For phased data: 1 if either gamete carries the allele (loses
    copy-number information vs `"additive_012"`). For unphased data:
    identical to `"additive_012"` since that already gives 0/1. May be
    preferable for Bayesian variable selection models (BayesB, BayesC)
    where the prior expects binary indicators.

- missing_string:

  Missing data marker. Default ".".

- scale_features:

  Center and scale columns. Default FALSE.

- min_freq:

  Minimum allele frequency to include. Default 0.01.

## Value

Numeric matrix (individuals x haplotype allele columns).

## Phased vs unphased haplotypes

**Phased data** (from
[`read_phased_vcf`](https://FAkohoue.github.io/LDxBlocks/reference/read_phased_vcf.md),
[`phase_with_beagle`](https://FAkohoue.github.io/LDxBlocks/reference/phase_with_beagle.md)):
each individual's block string is `"g1|g2"` where `g1` and `g2` are the
two gametic sequences. Haplotype alleles are identified at the *gamete*
level - true haplotypes in the biological sense. Frequencies are gamete
frequencies (each individual contributes two observations). Dosage
values 0, 1, 2 measure actual allele copy number.

**Unphased data** (from an unphased VCF or dosage matrix): each
individual's block string is a single multi-SNP dosage string (e.g.
`"021002"`). Haplotype *alleles* are distinct dosage patterns, not true
gametic haplotypes. Frequencies measure the proportion of *individuals*
carrying each dosage pattern. This is biologically meaningful (distinct
genotypic patterns at the block level) but should not be equated with
true haplotype allele frequency without phasing. The recommended
workflow for true haplotype analysis is: phase first with Beagle, then
call
[`extract_haplotypes()`](https://FAkohoue.github.io/LDxBlocks/reference/extract_haplotypes.md)
on the phased output.

## References

Difabachew YF et al. (2023). Genomic prediction with haplotype blocks in
wheat. *Frontiers in Plant Science* **14**:1168547.
[doi:10.3389/fpls.2023.1168547](https://doi.org/10.3389/fpls.2023.1168547)

Weber SE, Frisch M, Snowdon RJ, Voss-Fels KP (2023). Haplotype blocks
for genomic prediction: a comparative evaluation in multiple crop
datasets. *Frontiers in Plant Science* **14**:1217589.
[doi:10.3389/fpls.2023.1217589](https://doi.org/10.3389/fpls.2023.1217589)
