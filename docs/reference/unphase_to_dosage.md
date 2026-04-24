# Collapse Phased Gametes to 0/1/2 Dosage

Collapse Phased Gametes to 0/1/2 Dosage

## Usage

``` r
unphase_to_dosage(phased_list)
```

## Arguments

- phased_list:

  List from
  [`read_phased_vcf`](https://FAkohoue.github.io/LDxBlocks/reference/read_phased_vcf.md).

## Value

Numeric matrix (individuals x SNPs, 0/1/2/NA).
