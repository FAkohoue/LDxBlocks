# Example SNP Information Table

Genomic metadata for the 230 SNPs in
[`ldx_geno`](https://FAkohoue.github.io/LDxBlocks/reference/ldx_geno.md).

## Usage

``` r
ldx_snp_info
```

## Format

A `data.frame` with 230 rows and 5 columns:

- `SNP`:

  Character. SNP identifier (e.g. `"rs1001"`).

- `CHR`:

  Character. Chromosome label (`"1"`, `"2"`, or `"3"`). No `chr`
  prefix - consistent with the normalisation applied by
  [`read_geno`](https://FAkohoue.github.io/LDxBlocks/reference/read_geno.md).

- `POS`:

  Integer. Base-pair position on the chromosome. Positions increase
  monotonically within each chromosome; inter-block gaps of
  approximately 50,000 bp separate LD blocks.

- `REF`:

  Character. Reference allele (A/C/G/T).

- `ALT`:

  Character. Alternate allele.

## Source

Simulated with `data-raw/generate_example_data.R`.

## See also

[`ldx_geno`](https://FAkohoue.github.io/LDxBlocks/reference/ldx_geno.md),
[`ldx_blocks`](https://FAkohoue.github.io/LDxBlocks/reference/ldx_blocks.md)

## Examples

``` r
data(ldx_snp_info)
head(ldx_snp_info)
#>      SNP CHR  POS REF ALT
#> 1 rs1001   1 1000   G   C
#> 2 rs1002   1 2192   C   G
#> 3 rs1003   1 3253   T   A
#> 4 rs1004   1 4132   A   G
#> 5 rs1005   1 5188   G   A
#> 6 rs1006   1 6314   C   T
table(ldx_snp_info$CHR)   # 80 80 70
#> 
#>  1  2  3 
#> 80 80 70 

# \donttest{
# Compute LD decay and chromosome-specific decay distances
data(ldx_geno)
decay <- compute_ld_decay(
  geno         = ldx_geno,
  snp_info     = ldx_snp_info,
  sampling     = "random",
  r2_threshold = "both",
  n_pairs      = 3000L,
  verbose      = FALSE
)
decay$critical_r2_param   # background LD level
#>        95% 
#> 0.03718883 
decay$decay_dist          # per-chromosome decay distances (kb)
#>   CHR decay_dist_bp decay_dist_kb threshold_used r2_col_used censored
#> 1   1         39475         39.48     0.03718883    r2_loess    FALSE
#> 2   2         50920         50.92     0.03718883    r2_loess    FALSE
#> 3   3         29560         29.56     0.03718883    r2_loess    FALSE
# }
```
