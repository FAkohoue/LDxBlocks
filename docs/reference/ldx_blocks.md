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

A `data.frame` with 9 rows and 8 columns:

- `start`:

  Integer. Index of the first SNP in the block (1-based, over the full
  230-SNP set).

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

  Integer. `end.bp - start.bp + 1`.

## Source

Derived analytically from the founder-haplotype simulation in
`data-raw/generate_example_data.R`; matches expected
`run_Big_LD_all_chr(ldx_geno, ldx_snp_info, CLQcut = 0.6)` output.

## See also

[`ldx_geno`](https://FAkohoue.github.io/LDxBlocks/reference/ldx_geno.md),
[`run_Big_LD_all_chr`](https://FAkohoue.github.io/LDxBlocks/reference/run_Big_LD_all_chr.md),
[`summarise_blocks`](https://FAkohoue.github.io/LDxBlocks/reference/summarise_blocks.md)

## Examples

``` r
data(ldx_blocks)
ldx_blocks
#>   start end start.rsID end.rsID start.bp end.bp CHR length_bp
#> 1     1  25     rs1001   rs1025     1000  25027   1     24028
#> 2    31  50     rs1031   rs1050    81064  99022   1     17959
#> 3    56  80     rs1056   rs1080   155368 179371   1     24004
#> 4    81 110     rs2001   rs2030     1000  30023   2     29024
#> 5   116 135     rs2036   rs2055    86236 105290   2     19055
#> 6   141 160     rs2061   rs2080   161515 180473   2     18959
#> 7   161 180     rs3001   rs3020     1000  19068   3     18069
#> 8   186 205     rs3026   rs3045    74532  93854   3     19323
#> 9   211 230     rs3051   rs3070   149647 168376   3     18730
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
# }
```
