# Example GWAS Marker Table

Twenty toy GWAS-significant markers drawn from the simulated LD blocks,
used to demonstrate
[`tune_LD_params`](https://FAkohoue.github.io/LDxBlocks/reference/tune_LD_params.md)
parameter auto-tuning. All 20 markers fall within one of the nine
simulated blocks in
[`ldx_blocks`](https://FAkohoue.github.io/LDxBlocks/reference/ldx_blocks.md),
so a correctly parameterised run of `tune_LD_params` should return zero
unassigned markers.

## Usage

``` r
ldx_gwas
```

## Format

A `data.frame` with 20 rows and 4 columns:

- `Marker`:

  Character. SNP identifier, matching entries in `ldx_snp_info$SNP`.

- `CHR`:

  Character. Chromosome label.

- `POS`:

  Integer. Base-pair position.

- `P`:

  Numeric. Toy p-value in (1e-8, 0.001); not used by `tune_LD_params`
  itself but useful for demonstration.

## Source

Simulated with `data-raw/generate_example_data.R`. Seed: `set.seed(42)`.

## See also

[`tune_LD_params`](https://FAkohoue.github.io/LDxBlocks/reference/tune_LD_params.md),
[`ldx_geno`](https://FAkohoue.github.io/LDxBlocks/reference/ldx_geno.md),
[`ldx_snp_info`](https://FAkohoue.github.io/LDxBlocks/reference/ldx_snp_info.md)

## Examples

``` r
data(ldx_gwas)
head(ldx_gwas)
#>   Marker CHR   POS            P
#> 1 rs1001   1  1000 0.0009057391
#> 2 rs1005   1  5246 0.0004469752
#> 3 rs1017   1 17510 0.0008360059
#> 4 rs1034   1 84716 0.0007375982
#> 5 rs1040   1 90503 0.0008110570
#> 6 rs1048   1 98806 0.0003881144
# All markers should fall within ldx_blocks:
data(ldx_blocks)
data(ldx_snp_info)
all(ldx_gwas$Marker %in% ldx_snp_info$SNP)
#> [1] TRUE
```
