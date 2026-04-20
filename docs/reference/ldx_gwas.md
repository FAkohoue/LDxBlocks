# Example GWAS Marker Table

Twenty toy GWAS markers drawn from simulated LD blocks for use with
[`tune_LD_params`](https://FAkohoue.github.io/LDxBlocks/reference/tune_LD_params.md)
and
[`define_qtl_regions`](https://FAkohoue.github.io/LDxBlocks/reference/define_qtl_regions.md).
All 20 markers fall within LD blocks, but not all nine blocks are
represented: chromosome 3 contributes markers only from its first block;
blocks 2 and 3 of chromosome 3 have no GWAS markers. A correctly
parameterised run of `tune_LD_params` should return zero unassigned
markers.

## Usage

``` r
ldx_gwas
```

## Format

A `data.frame` with 20 rows and 5 columns:

- `Marker`:

  Character. SNP identifier, matching entries in `ldx_snp_info$SNP`.

- `CHR`:

  Character. Chromosome label.

- `POS`:

  Integer. Base-pair position.

- `P`:

  Numeric. Toy p-value. Eight markers have P \< 1e-6 (clearly
  significant), seven have P in (1e-5, 1e-3) (suggestive), and five have
  P in (1e-3, 0.05) (sub-threshold).

- `trait`:

  Character. Toy trait label (`"TraitA"` or `"TraitB"`), assigned
  randomly at 60/40 probability. Used to demonstrate multi-trait
  pleiotropic block detection in
  [`define_qtl_regions`](https://FAkohoue.github.io/LDxBlocks/reference/define_qtl_regions.md).
  When a single trait is analysed, this column is simply ignored or the
  `trait_col` argument omitted.

## Source

Simulated with `data-raw/generate_example_data.R`. Seed: `set.seed(42)`.

## See also

[`tune_LD_params`](https://FAkohoue.github.io/LDxBlocks/reference/tune_LD_params.md),
[`define_qtl_regions`](https://FAkohoue.github.io/LDxBlocks/reference/define_qtl_regions.md),
[`ldx_geno`](https://FAkohoue.github.io/LDxBlocks/reference/ldx_geno.md),
[`ldx_snp_info`](https://FAkohoue.github.io/LDxBlocks/reference/ldx_snp_info.md)

## Examples

``` r
data(ldx_gwas)
head(ldx_gwas)
#>   Marker CHR   POS            P  trait
#> 1 rs1001   1  1000 8.151737e-07 TraitB
#> 2 rs1005   1  5188 4.023280e-07 TraitB
#> 3 rs1017   1 16833 7.524202e-07 TraitB
#> 4 rs1034   1 83635 6.638623e-07 TraitA
#> 5 rs1040   1 89846 7.299685e-07 TraitA
#> 6 rs1048   1 97212 3.493586e-07 TraitA
# All markers should fall within ldx_blocks:
data(ldx_blocks)
data(ldx_snp_info)
all(is.element(ldx_gwas$Marker, ldx_snp_info$SNP))
#> [1] TRUE
```
