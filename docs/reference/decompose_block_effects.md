# Decompose Per-SNP Effects into Per-Haplotype-Allele Effect Table

Aggregates the per-SNP additive effects (from
[`backsolve_snp_effects`](https://FAkohoue.github.io/LDxBlocks/reference/backsolve_snp_effects.md))
into a per-allele effect for each LD block: the effect of carrying
allele h is the sum of SNP effects weighted by the allele's dosage at
each SNP position. Returns an interpretable table of "which allele of
which block is worth how many units of the trait?"

## Usage

``` r
decompose_block_effects(
  haplotypes,
  snp_info,
  blocks,
  snp_effects,
  min_freq = 0.02,
  missing_string = "."
)
```

## Arguments

- haplotypes:

  Named list from
  [`extract_haplotypes`](https://FAkohoue.github.io/LDxBlocks/reference/extract_haplotypes.md).

- snp_info:

  Data frame with `SNP`, `CHR`, `POS`.

- blocks:

  LD block table from
  [`run_Big_LD_all_chr`](https://FAkohoue.github.io/LDxBlocks/reference/run_Big_LD_all_chr.md).

- snp_effects:

  Named numeric vector of per-SNP additive effects, as returned by
  [`backsolve_snp_effects`](https://FAkohoue.github.io/LDxBlocks/reference/backsolve_snp_effects.md).

- min_freq:

  Numeric. Alleles below this frequency are excluded. Default `0.02`.

- missing_string:

  Character. Missing haplotype placeholder. Default `"."`.

## Value

Data frame with one row per allele per block, sorted by `CHR`,
`start_bp`, `effect_rank`.

- `block_id`, `CHR`, `start_bp`, `end_bp`:

  Block coordinates.

- `allele`:

  Haplotype allele string.

- `frequency`:

  Allele frequency in the panel.

- `allele_effect`:

  Sum of SNP effects weighted by allele dosage.

- `effect_rank`:

  Rank within block; 1 = most positive effect.

- `n_snps_block`:

  Number of SNPs in this block.

## See also

[`backsolve_snp_effects`](https://FAkohoue.github.io/LDxBlocks/reference/backsolve_snp_effects.md),
[`rank_haplotype_blocks`](https://FAkohoue.github.io/LDxBlocks/reference/rank_haplotype_blocks.md)

## Examples

``` r
# \donttest{
data(ldx_geno, ldx_snp_info, ldx_blocks, ldx_blues, package = "LDxBlocks")
haps    <- extract_haplotypes(ldx_geno, ldx_snp_info, ldx_blocks)
res     <- run_haplotype_prediction(ldx_geno, ldx_snp_info, ldx_blocks,
                                     blues = ldx_blues, id_col = "id",
                                     blue_col = "YLD", verbose = FALSE)
# Single-trait: res$snp_effects is the named numeric vector directly
snp_fx  <- res$snp_effects
allele_tbl <- decompose_block_effects(haps, ldx_snp_info, ldx_blocks,
                                       snp_effects = snp_fx)
head(allele_tbl[order(-allele_tbl$allele_effect), ])
#>                                          block_id CHR start_bp end_bp
#> 222222222222222222222022202222 block_2_1000_30023   2     1000  30023
#> 122121212112211211122011201112 block_2_1000_30023   2     1000  30023
#> 111122112222221222211122211112 block_2_1000_30023   2     1000  30023
#> 022020202002200200022000200002 block_2_1000_30023   2     1000  30023
#> 011021102112210211111111210002 block_2_1000_30023   2     1000  30023
#> 000022002222220222200222220002 block_2_1000_30023   2     1000  30023
#>                                                        allele frequency
#> 222222222222222222222022202222 222222222222222222222022202222    0.0250
#> 122121212112211211122011201112 122121212112211211122011201112    0.0417
#> 111122112222221222211122211112 111122112222221222211122211112    0.0917
#> 022020202002200200022000200002 022020202002200200022000200002    0.0333
#> 011021102112210211111111210002 011021102112210211111111210002    0.0917
#> 000022002222220222200222220002 000022002222220222200222220002    0.1000
#>                                allele_effect n_snps_block effect_rank
#> 222222222222222222222022202222      0.790777           30           1
#> 122121212112211211122011201112      0.704303           30           2
#> 111122112222221222211122211112      0.672643           30           3
#> 022020202002200200022000200002      0.617828           30           4
#> 011021102112210211111111210002      0.586168           30           5
#> 000022002222220222200222220002      0.554508           30           6
# }
```
