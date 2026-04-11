# Compute Haplotype Diversity Per Block

Calculates n_haplotypes, He, Shannon entropy, and freq_dominant per
block. Phased data contributes two gamete observations per individual.

## Usage

``` r
compute_haplotype_diversity(haplotypes, missing_string = ".")
```

## Arguments

- haplotypes:

  List from extract_haplotypes().

- missing_string:

  Missing data marker. Default ".".

## Value

Data frame with block_id, CHR, start_bp, end_bp, n_snps, n_ind,
n_haplotypes, He, Shannon, freq_dominant, phased.
