# Compute Haplotype Diversity Per Block

Calculates per-block haplotype diversity metrics: richness
(n_haplotypes), expected heterozygosity (He, Nei 1973 sample-size
corrected), Shannon entropy, effective number of alleles (1/\\\sum
p_i^2\\), dominant haplotype frequency, and a sweep flag (TRUE when
freq_dominant \\\geq\\ 0.90). These metrics directly correspond to those
used to characterise block diversity and identify selection signatures
in Difabachew et al. (2023) and Tong et al. (2024). Phased data
contributes two gamete observations per individual, doubling the
effective sample size.

## Usage

``` r
compute_haplotype_diversity(haplotypes, missing_string = ".")
```

## Arguments

- haplotypes:

  Named list from
  [`extract_haplotypes`](https://FAkohoue.github.io/LDxBlocks/reference/extract_haplotypes.md).

- missing_string:

  Missing data marker. Default `"."`.

## Value

Data frame with one row per block: `block_id`, `CHR`, `start_bp`,
`end_bp`, `n_snps`, `n_ind`, `n_haplotypes`, `He` (corrected),
`Shannon`, `n_eff_alleles`, `freq_dominant`, `sweep_flag`, `phased`.

## References

Nei M (1973). Analysis of gene diversity in subdivided populations.
*Proceedings of the National Academy of Sciences* **70**(12):3321-3323.
[doi:10.1073/pnas.70.12.3321](https://doi.org/10.1073/pnas.70.12.3321)

Difabachew YF et al. (2023). Genomic prediction with haplotype blocks in
wheat. *Frontiers in Plant Science* **14**:1168547.
[doi:10.3389/fpls.2023.1168547](https://doi.org/10.3389/fpls.2023.1168547)

Tong J et al. (2024). Stacking beneficial haplotypes from the Vavilov
wheat collection to accelerate breeding for multiple disease resistance.
*Theoretical and Applied Genetics* **137**:274.
[doi:10.1007/s00122-024-04784-w](https://doi.org/10.1007/s00122-024-04784-w)
