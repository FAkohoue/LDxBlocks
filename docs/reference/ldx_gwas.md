# Example GWAS Marker Table

Twenty toy GWAS markers drawn from simulated LD blocks for use with
[`tune_LD_params`](https://FAkohoue.github.io/LDxBlocks/reference/tune_LD_params.md),
[`define_qtl_regions`](https://FAkohoue.github.io/LDxBlocks/reference/define_qtl_regions.md),
and
[`compare_gwas_effects`](https://FAkohoue.github.io/LDxBlocks/reference/compare_gwas_effects.md).
All 20 markers fall within LD blocks, but not all nine blocks are
represented: chromosome 3 contributes markers only from its first block;
blocks 2 and 3 of chromosome 3 have no GWAS markers. A correctly
parameterised run of `tune_LD_params` should return zero unassigned
markers.

The dataset includes effect sizes (`BETA`, `SE`) so that it can serve as
a minimal external GWAS result for `compare_gwas_effects` demonstrations
without requiring a full association analysis.

## Usage

``` r
ldx_gwas
```

## Format

A `data.frame` with 20 rows and 7 columns:

- `Marker`:

  Character. SNP identifier, matching entries in `ldx_snp_info$SNP`.

- `CHR`:

  Character. Chromosome label.

- `POS`:

  Integer. Base-pair position.

- `BETA`:

  Numeric. Simulated effect size. Derived from the z-score of `P` with
  5% random noise. Positive for the first 10 markers, negative for the
  last 10, providing a testable directional signal for
  [`compare_gwas_effects`](https://FAkohoue.github.io/LDxBlocks/reference/compare_gwas_effects.md)
  demonstrations.

- `SE`:

  Numeric. Standard error of `BETA`. Derived as `|BETA|/|z|` inflated by
  10%, consistent with standard mixed-model GWAS output.

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
[`compare_gwas_effects`](https://FAkohoue.github.io/LDxBlocks/reference/compare_gwas_effects.md),
[`ldx_geno`](https://FAkohoue.github.io/LDxBlocks/reference/ldx_geno.md),
[`ldx_snp_info`](https://FAkohoue.github.io/LDxBlocks/reference/ldx_snp_info.md)

## Examples

``` r
data(ldx_gwas)
head(ldx_gwas)
#>   Marker CHR   POS   BETA     SE            P  trait
#> 1 rs1001   1  1000 0.2510 0.0560 8.151737e-07 TraitA
#> 2 rs1005   1  5188 0.2572 0.0558 4.023280e-07 TraitA
#> 3 rs1017   1 16833 0.2425 0.0539 7.524202e-07 TraitA
#> 4 rs1034   1 83635 0.2584 0.0572 6.638623e-07 TraitB
#> 5 rs1040   1 89846 0.2560 0.0569 7.299685e-07 TraitA
#> 6 rs1048   1 97212 0.2426 0.0524 3.493586e-07 TraitA
# All markers should fall within ldx_blocks:
data(ldx_blocks)
data(ldx_snp_info)
all(is.element(ldx_gwas$Marker, ldx_snp_info$SNP))
#> [1] TRUE

# \donttest{
# Map markers to LD blocks and identify lead SNPs per block
qtl <- define_qtl_regions(
  gwas_results = ldx_gwas,
  blocks       = ldx_blocks,
  snp_info     = ldx_snp_info,
  p_threshold  = 1e-5,
  verbose      = FALSE
)
qtl[, c("block_id", "CHR", "lead_snp", "lead_p", "lead_beta", "n_sig_markers")]
#>                block_id CHR lead_snp       lead_p lead_beta n_sig_markers
#> 1    block_1_1000_25027   1   rs1005 4.023280e-07    0.2572             3
#> 2   block_1_81064_99022   1   rs1048 3.493586e-07    0.2426             3
#> 3 block_1_155368_179371   1   rs1070 3.653110e-09    0.3262             2

# Cross-population effect concordance from external GWAS results.
# Here ldx_gwas is used for both populations; in practice these
# would be independent GWAS runs on separate cohorts.
set.seed(1L)
gwas_pop2        <- ldx_gwas
gwas_pop2$BETA   <- ldx_gwas$BETA * 0.8 + rnorm(20, 0, 0.02)
gwas_pop2$SE     <- ldx_gwas$SE   * 0.9 + runif(20, 0, 0.005)
gwas_pop2$P      <- 2 * pnorm(-abs(gwas_pop2$BETA / gwas_pop2$SE))

# Path 1: pre-mapped (recommended for auditability)
qtl1 <- define_qtl_regions(ldx_gwas,  ldx_blocks, ldx_snp_info,
                           p_threshold = NULL, verbose = FALSE)
qtl2 <- define_qtl_regions(gwas_pop2, ldx_blocks, ldx_snp_info,
                           p_threshold = NULL, verbose = FALSE)
conc <- compare_gwas_effects(
  qtl_pop1    = qtl1,
  qtl_pop2    = qtl2,
  blocks_pop1 = ldx_blocks,
  blocks_pop2 = ldx_blocks,
  pop1_name   = "DiscoveryPop",
  pop2_name   = "ReplicationPop",
  verbose     = FALSE
)
conc$concordance[, c("block_id","direction_agreement","meta_p","replicated")]
#>                block_id direction_agreement       meta_p replicated
#> 1    block_1_1000_25027                   1 1.855521e-10       TRUE
#> 2   block_1_81064_99022                   1 1.687084e-11       TRUE
#> 3 block_1_155368_179371                   1 3.037882e-14       TRUE
#> 4    block_2_1000_30023                   1 1.409021e-07       TRUE
#> 5  block_2_86236_105290                   1 3.304721e-07       TRUE
#> 6    block_3_1000_19068                   1 7.263598e-03       TRUE

# Path 2: raw GWAS + blocks (calls define_qtl_regions internally)
conc2 <- compare_gwas_effects(
  gwas_pop1     = ldx_gwas,
  gwas_pop2     = gwas_pop2,
  blocks_pop1   = ldx_blocks,
  blocks_pop2   = ldx_blocks,
  snp_info_pop1 = ldx_snp_info,
  p_threshold   = NULL,
  verbose       = FALSE
)
print(conc2)
#> LDxBlocks Cross-Population Effect Concordance
#>   Populations: pop1 vs pop2
#>   Traits:       trait 
#>   Blocks compared:           6 
#>   With enough shared alleles: 6 
#>   Directionally concordant:   6 
#>   Replicated (dir + Q_p>0.05): 6 
#>   Boundary warnings:          0 (overlap ratio < 0.8 )
#>   Shared allele comparisons:  6 
# }
```
