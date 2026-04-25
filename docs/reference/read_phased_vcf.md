# Read Pre-Phased VCF

Read Pre-Phased VCF

## Usage

``` r
read_phased_vcf(vcf_file, min_maf = 0, verbose = TRUE)
```

## Arguments

- vcf_file:

  Path to phased VCF or VCF.gz with 0\|1 GT fields.

- min_maf:

  Minimum MAF. Default 0.0.

- verbose:

  Logical. Default TRUE.

## Value

List: hap1, hap2 (SNPs x individuals, 0/1), dosage (0/1/2), snp_info,
sample_ids.
