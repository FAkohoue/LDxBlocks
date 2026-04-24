# Example LD Block Table

A reference LD block table for
[`ldx_geno`](https://FAkohoue.github.io/LDxBlocks/reference/ldx_geno.md),
matching the nine blocks (three per chromosome) embedded in the
simulation. This object is used in examples and tests as the expected
output of
[`run_Big_LD_all_chr`](https://FAkohoue.github.io/LDxBlocks/reference/run_Big_LD_all_chr.md)
on `ldx_geno`.

## Usage

``` r
ldx_blocks
```

## Format

A `data.frame` with 9 rows and 9 columns:

- `start`:

  Integer. Index of the first SNP in the block (1-based, over the full
  230-SNP set including monomorphics).

- `end`:

  Integer. Index of the last SNP.

- `start.rsID`:

  Character. SNP ID at the block start.

- `end.rsID`:

  Character. SNP ID at the block end.

- `start.bp`:

  Integer. Base-pair position of the block start.

- `end.bp`:

  Integer. Base-pair position of the block end.

- `CHR`:

  Character. Chromosome label.

- `length_bp`:

  Integer. `end.bp - start.bp + 1`. Physical span of the block in base
  pairs.

- `n_snps`:

  Integer. Number of polymorphic SNPs in the block (those that passed
  the MAF filter and were assigned to a clique). For singletons
  `n_snps = 1`. Note that `end - start + 1` may exceed `n_snps` when
  monomorphic SNPs fall within the block boundaries.

## Source

Derived analytically from the founder-haplotype simulation in
`data-raw/generate_example_data.R`; matches expected
`run_Big_LD_all_chr(ldx_geno, ldx_snp_info, CLQcut = 0.6)` output.

## See also

[`ldx_geno`](https://FAkohoue.github.io/LDxBlocks/reference/ldx_geno.md),
[`run_Big_LD_all_chr`](https://FAkohoue.github.io/LDxBlocks/reference/run_Big_LD_all_chr.md),
[`summarise_blocks`](https://FAkohoue.github.io/LDxBlocks/reference/summarise_blocks.md),
[`compare_block_effects`](https://FAkohoue.github.io/LDxBlocks/reference/compare_block_effects.md),
[`compare_gwas_effects`](https://FAkohoue.github.io/LDxBlocks/reference/compare_gwas_effects.md)

## Examples

``` r
data(ldx_blocks)
ldx_blocks
#>   start end start.rsID end.rsID start.bp end.bp CHR length_bp n_snps
#> 1     1  25     rs1001   rs1025     1000  25027   1     24028     25
#> 2    31  50     rs1031   rs1050    81064  99022   1     17959     20
#> 3    56  80     rs1056   rs1080   155368 179371   1     24004     25
#> 4    81 110     rs2001   rs2030     1000  30023   2     29024     30
#> 5   116 135     rs2036   rs2055    86236 105290   2     19055     20
#> 6   141 160     rs2061   rs2080   161515 180473   2     18959     20
#> 7   161 180     rs3001   rs3020     1000  19068   3     18069     20
#> 8   186 205     rs3026   rs3045    74532  93854   3     19323     20
#> 9   211 230     rs3051   rs3070   149647 168376   3     18730     20
summarise_blocks(ldx_blocks)
#>      CHR n_blocks min_bp median_bp  mean_bp max_bp total_bp_covered
#> 1      1        3  17959     24004 21997.00  24028            65991
#> 2      2        3  18959     19055 22346.00  29024            67038
#> 3      3        3  18069     18730 18707.33  19323            56122
#> 4 GENOME        9  17959     19055 21016.78  29024           189151

# \donttest{
# Compare haplotype allele frequencies between two sample groups
data(ldx_geno, ldx_snp_info)
haps <- extract_haplotypes(ldx_geno, ldx_snp_info, ldx_blocks)
ids  <- rownames(ldx_geno)
cmp  <- compare_haplotype_populations(
  haplotypes  = haps,
  group1      = ids[1:60],
  group2      = ids[61:120],
  group1_name = "cycle1",
  group2_name = "cycle2"
)
head(cmp[, c("block_id","n1","n2","FST","divergent")])
#>                block_id n1 n2    FST divergent
#> 1    block_1_1000_25027 60 60 0.0096     FALSE
#> 2   block_1_81064_99022 60 60 0.0068     FALSE
#> 3 block_1_155368_179371 60 60 0.0072     FALSE
#> 4    block_2_1000_30023 60 60 0.0046     FALSE
#> 5  block_2_86236_105290 60 60 0.0034     FALSE
#> 6 block_2_161515_180473 60 60 0.0099     FALSE

# Cross-population haplotype effect concordance.
# test_block_haplotypes() produces LDxBlocks_haplotype_assoc objects;
# compare_block_effects() compares their allele effects block by block.
data(ldx_blues)
blues_vec <- setNames(ldx_blues$YLD, ldx_blues$id)
# Simulate two populations by splitting the 120 individuals
set.seed(2L)
idx1 <- sample(ids, 70L);  idx2 <- setdiff(ids, idx1)
haps1 <- extract_haplotypes(ldx_geno[idx1, ], ldx_snp_info, ldx_blocks)
haps2 <- extract_haplotypes(ldx_geno[idx2, ], ldx_snp_info, ldx_blocks)
assoc1 <- test_block_haplotypes(haps1, blues = blues_vec[idx1],
                                 blocks = ldx_blocks, n_pcs = 0L,
                                 verbose = FALSE)
assoc2 <- test_block_haplotypes(haps2, blues = blues_vec[idx2],
                                 blocks = ldx_blocks, n_pcs = 0L,
                                 verbose = FALSE)
# block_match = "position" handles different block boundaries;
# use "id" (default) when both populations share the same block table.
conc <- compare_block_effects(
  assoc1, assoc2,
  pop1_name   = "Pop1",
  pop2_name   = "Pop2",
  blocks_pop1 = ldx_blocks,
  blocks_pop2 = ldx_blocks,
  block_match = "id",
  verbose     = FALSE
)
conc$concordance[, c("block_id","n_shared_alleles",
                     "direction_agreement","meta_p","replicated")]
#>                block_id n_shared_alleles direction_agreement    meta_p
#> 1    block_1_1000_25027                7              0.4286 0.7964556
#> 2   block_1_81064_99022                5              0.4000 0.7000451
#> 3 block_1_155368_179371                7              0.4286 0.8408251
#> 4    block_2_1000_30023                7              0.5714 0.6899425
#> 5  block_2_86236_105290                7              0.4286 0.7257983
#> 6 block_2_161515_180473                6              0.6667 0.7767358
#> 7    block_3_1000_19068                7              0.8571 0.5867187
#> 8   block_3_74532_93854                6              0.5000 0.8366144
#> 9 block_3_149647_168376                8              0.5000 0.9410410
#>   replicated
#> 1      FALSE
#> 2      FALSE
#> 3      FALSE
#> 4      FALSE
#> 5      FALSE
#> 6      FALSE
#> 7       TRUE
#> 8      FALSE
#> 9      FALSE
# }
```
