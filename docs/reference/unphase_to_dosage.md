# Collapse Phased Gametes to 0/1/2 Dosage

Collapse Phased Gametes to 0/1/2 Dosage

## Usage

``` r
unphase_to_dosage(phased_list)
```

## Arguments

- phased_list:

  List from read_phased_vcf() or phase_with_pedigree().

## Value

Numeric matrix (individuals x SNPs, 0/1/2/NA).
