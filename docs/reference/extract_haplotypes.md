# Extract Haplotype Dosage Strings from LD Blocks

Builds per-block haplotype dosage strings for all individuals across the
LD blocks in `blocks`. Each block is processed by the C++ engine
`extract_chr_haplotypes_cpp()` (unphased) or
`extract_chr_haplotypes_phased_cpp()` (phased VCF input), which assigns
each individual a dosage string of 0/1/2 characters (one per SNP in the
block) and identifies the top haplotype alleles by frequency.

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

  One of:

  - An `LDxBlocks_backend` from
    [`read_geno`](https://FAkohoue.github.io/LDxBlocks/reference/read_geno.md)
    or
    [`read_geno_bigmemory`](https://FAkohoue.github.io/LDxBlocks/reference/read_geno_bigmemory.md)
    (streaming, one chromosome at a time).

  - A named list with elements `hap1` and `hap2` (phased SNPs x
    individuals matrices from
    [`read_phased_vcf`](https://FAkohoue.github.io/LDxBlocks/reference/read_phased_vcf.md)).

  - A numeric matrix (individuals x SNPs, values 0/1/2/NA).

- snp_info:

  Data frame with columns `SNP`, `CHR`, `POS`.

- blocks:

  Data frame of LD blocks from
  [`run_Big_LD_all_chr`](https://FAkohoue.github.io/LDxBlocks/reference/run_Big_LD_all_chr.md),
  with columns `CHR`, `start.bp`, `end.bp`, `n_snps`.

- chr:

  Character vector of chromosomes to process. `NULL` (default) processes
  all chromosomes present in `blocks`.

- min_snps:

  Integer. Minimum number of SNPs a block must contain to be included.
  Default `3L`.

- na_char:

  Character. Symbol used to denote missing genotype in the dosage
  string. Default `"."`.

## Value

A named list of per-block haplotype dosage matrices (individuals x
haplotype alleles, values 0/1/2 for phased data or 0/1 for unphased).
The list carries a `block_info` attribute (data frame with one row per
block: `block_id`, `CHR`, `start_bp`, `end_bp`, `n_snps`,
`n_haplotypes`, `phased`).

## See also

[`run_Big_LD_all_chr`](https://FAkohoue.github.io/LDxBlocks/reference/run_Big_LD_all_chr.md),
[`build_haplotype_feature_matrix`](https://FAkohoue.github.io/LDxBlocks/reference/build_haplotype_feature_matrix.md),
[`compute_haplotype_diversity`](https://FAkohoue.github.io/LDxBlocks/reference/compute_haplotype_diversity.md),
[`decode_haplotype_strings`](https://FAkohoue.github.io/LDxBlocks/reference/decode_haplotype_strings.md)

## Examples

``` r
data(ldx_geno, ldx_snp_info, ldx_blocks)
haps <- extract_haplotypes(ldx_geno, ldx_snp_info, ldx_blocks, min_snps = 3L)
length(haps)                     # one element per block
#> [1] 9
names(haps)[1]                   # e.g. "block_1_1000_25000"
#> [1] "block_1_1000_25027"
dim(haps[[1]])                   # individuals x haplotype alleles
#> NULL
```
