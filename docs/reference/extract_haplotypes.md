# Extract Phase-Free Haplotypes Within LD Blocks

For each LD block in `blocks`, extracts the SNP columns from `geno` that
fall within the block's bp interval and constructs a \*phase-free
haplotype\* for each individual by concatenating their allele codes (0,
1, or 2) into a character string. Missing genotypes are represented as
`"."`.

These haplotype strings serve as multi-locus genotype identifiers and
can be used directly with
[`compute_haplotype_diversity`](https://FAkohoue.github.io/LDxBlocks/reference/compute_haplotype_diversity.md)
or
[`build_haplotype_feature_matrix`](https://FAkohoue.github.io/LDxBlocks/reference/build_haplotype_feature_matrix.md).

## Usage

``` r
extract_haplotypes(
  geno,
  snp_info,
  blocks,
  chr = NULL,
  min_snps = 2L,
  na_char = "."
)
```

## Arguments

- geno:

  Numeric matrix (individuals × SNPs; 0/1/2). Column names must be SNP
  IDs matching `snp_info$SNP`. Row names are used as individual IDs.

- snp_info:

  Data frame with columns `SNP`, `CHR`, `POS`.

- blocks:

  Data frame of LD blocks as returned by
  [`run_Big_LD_all_chr`](https://FAkohoue.github.io/LDxBlocks/reference/run_Big_LD_all_chr.md)
  or
  [`Big_LD`](https://FAkohoue.github.io/LDxBlocks/reference/Big_LD.md).
  Must include columns `start.bp`, `end.bp`, and (if `chr` is given)
  `CHR`.

- chr:

  Character. If supplied, only the blocks on this chromosome and the
  SNPs on this chromosome (from `snp_info`) are processed. `NULL`
  processes all blocks jointly.

- min_snps:

  Integer. Minimum number of SNPs required for a block to be included.
  Blocks with fewer SNPs are skipped. Default `2`.

- na_char:

  Character. Replacement for `NA` alleles in strings. Default `"."`.

## Value

A named list with one element per block (named
`"block_<start.bp>_<end.bp>"`). Each element is a character vector of
length `nrow(geno)`, giving the haplotype string for each individual.
The list has an attribute `"block_info"` — a data frame summarising
block ID, chromosome, positions, and SNP count.

## Note on phasing

True gametic haplotypes require statistical or read-based phasing (e.g.
SHAPEIT, Beagle). The strings produced here are \*diploid allele
strings\*, not gametic phases. They are nonetheless highly informative
for diversity metrics and genomic prediction because within a high-LD
block, unphased multi-locus codes are nearly 1:1 with true haplotype
classes.

## See also

[`compute_haplotype_diversity`](https://FAkohoue.github.io/LDxBlocks/reference/compute_haplotype_diversity.md),
[`build_haplotype_feature_matrix`](https://FAkohoue.github.io/LDxBlocks/reference/build_haplotype_feature_matrix.md)

## Examples

``` r
# \donttest{
set.seed(1)
geno <- matrix(sample(0:2, 50 * 40, replace = TRUE), 50, 40)
rownames(geno) <- paste0("ind", 1:50)
colnames(geno) <- paste0("rs", 1:40)
snp_info <- data.frame(
  SNP = colnames(geno),
  CHR = "chr1",
  POS = seq(1000, by = 5000, length.out = 40)
)
# Toy blocks
blocks <- data.frame(
  start.bp = c(1000, 80000),
  end.bp   = c(70000, 200000),
  CHR      = "chr1"
)
haps <- extract_haplotypes(geno, snp_info, blocks)
head(haps[[1]])
#>             ind1             ind2             ind3             ind4 
#> "02100002201212" "21120200122200" "01010110100211" "11200000111122" 
#>             ind5             ind6 
#> "01220101001002" "20100220002201" 
# }
```
