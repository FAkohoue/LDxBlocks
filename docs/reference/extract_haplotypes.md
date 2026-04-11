# Extract Haplotype Strings Within LD Blocks

Builds haplotype strings per individual per LD block. Unphased mode:
"012201" (diploid). Phased mode: "011\|100" (gametes). Blocks are always
defined within a single chromosome.

## Usage

``` r
extract_haplotypes(
  geno,
  snp_info,
  blocks,
  chr = NULL,
  min_snps = 3L,
  na_char = "."
)
```

## Arguments

- geno:

  Dosage matrix (individuals x SNPs) OR phased list from
  read_phased_vcf()/phase_with_pedigree().

- snp_info:

  Data frame with SNP, CHR, POS.

- blocks:

  LD block data frame from run_Big_LD_all_chr() with CHR, start.bp,
  end.bp.

- chr:

  Optional chromosome filter. Default NULL (all).

- min_snps:

  Minimum SNPs per block. Default 3L.

- na_char:

  Missing allele character. Default ".".

## Value

Named list of character vectors (length=n_individuals), one per block.
Carries block_info attribute.
