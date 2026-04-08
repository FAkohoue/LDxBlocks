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
  200-SNP set).

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

Derived analytically from the simulation structure in
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
#> 1     1  25     rs1001   rs1025     1000  25535   1     24536
#> 2    31  50     rs1031   rs1050    81986 100878   1     18893
#> 3    56  80     rs1056   rs1080   156776 181114   1     24339
#> 4    81 110     rs2001   rs2030     1000  29445   2     28446
#> 5   116 135     rs2036   rs2055    85463 104532   2     19070
#> 6   141 160     rs2061   rs2080   160237 178996   2     18760
#> 7   161 180     rs3001   rs3020     1000  19994   3     18995
#> 8   186 205     rs3026   rs3045    76186  95838   3     19653
#> 9   211 230     rs3051   rs3070   151654 171449   3     19796
summarise_blocks(ldx_blocks)
#>      CHR n_blocks min_bp median_bp  mean_bp max_bp total_bp_covered
#> 1      1        3  18893     24339 22589.33  24536            67768
#> 2      2        3  18760     19070 22092.00  28446            66276
#> 3      3        3  18995     19653 19481.33  19796            58444
#> 4 GENOME        9  18760     19653 21387.56  28446           192488
```
