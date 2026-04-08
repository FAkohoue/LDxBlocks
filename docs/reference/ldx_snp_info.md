# Example SNP Information Table

Genomic metadata for the 200 SNPs in
[`ldx_geno`](https://FAkohoue.github.io/LDxBlocks/reference/ldx_geno.md).

## Usage

``` r
ldx_snp_info
```

## Format

A `data.frame` with 200 rows and 5 columns:

- `SNP`:

  Character. SNP identifier (e.g. `"rs1001"`).

- `CHR`:

  Character. Chromosome label (`"1"`, `"2"`, or `"3"`). No `chr` prefix
  — consistent with the normalisation applied by
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
#> 2 rs1002   1 2134   C   G
#> 3 rs1003   1 3089   T   A
#> 4 rs1004   1 4067   A   G
#> 5 rs1005   1 5246   G   A
#> 6 rs1006   1 6312   C   T
table(ldx_snp_info$CHR)   # 70 70 60
#> 
#>  1  2  3 
#> 80 80 70 
```
