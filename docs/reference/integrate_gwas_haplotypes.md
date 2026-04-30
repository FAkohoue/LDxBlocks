# Integrate GWAS QTL Regions with Haplotype Prediction Results

Links the output of
[`define_qtl_regions`](https://FAkohoue.github.io/LDxBlocks/reference/define_qtl_regions.md)
(biological evidence from GWAS) with the output of
[`run_haplotype_prediction`](https://FAkohoue.github.io/LDxBlocks/reference/run_haplotype_prediction.md)
(statistical evidence from haplotype variance) to identify blocks that
are supported by both lines of evidence. These are the priority
candidates for haplotype stacking in breeding.

## Usage

``` r
integrate_gwas_haplotypes(
  qtl_regions,
  pred_result,
  diversity = NULL,
  He_threshold = 0.3
)
```

## Arguments

- qtl_regions:

  Data frame from
  [`define_qtl_regions`](https://FAkohoue.github.io/LDxBlocks/reference/define_qtl_regions.md).

- pred_result:

  List from
  [`run_haplotype_prediction`](https://FAkohoue.github.io/LDxBlocks/reference/run_haplotype_prediction.md).

- diversity:

  Data frame from
  [`compute_haplotype_diversity`](https://FAkohoue.github.io/LDxBlocks/reference/compute_haplotype_diversity.md).
  Optional. If supplied, adds He and sweep_flag to the output.

- He_threshold:

  Minimum expected heterozygosity to flag a block as sufficiently
  diverse for haplotype stacking. Default `0.3`.

## Value

Data frame with one row per block that appears in at least one of the
input sources, with columns: `block_id`, `CHR`, `start_bp`, `end_bp`,
`n_snps`, `has_gwas_hit` (logical), `lead_marker`, `lead_p`,
`lead_beta`, `n_sig_markers`, `is_important` (logical, scaled var \>=
0.9), `var_scaled`, `He`, `sweep_flag`, `is_diverse` (logical, He \>=
He_threshold), `priority_score` (0-3), `recommendation` (character
label). Sorted by `priority_score` descending, then `var_scaled`
descending.

## Details

Three complementary evidence layers are combined per block:

- Biological (GWAS):

  Does the block contain a genome-wide significant marker? Sourced from
  [`define_qtl_regions()`](https://FAkohoue.github.io/LDxBlocks/reference/define_qtl_regions.md).

- Statistical (variance):

  Does the block explain substantial variance in the trait? Blocks with
  scaled `Var(local GEBV)` \>= 0.9 are flagged `important` by
  [`run_haplotype_prediction()`](https://FAkohoue.github.io/LDxBlocks/reference/run_haplotype_prediction.md).

- Diversity (He):

  Does the block have enough haplotype diversity to stack favourable
  alleles? Sourced from
  [`compute_haplotype_diversity()`](https://FAkohoue.github.io/LDxBlocks/reference/compute_haplotype_diversity.md).

The output `priority_score` is the sum of binary flags for each layer
(0-3). Blocks scoring 3 are supported by all three lines of evidence and
are the strongest candidates for haplotype stacking. Blocks scoring 2
are worth investigating. Blocks scoring 1 or 0 require caution.

## Interpretation guide

|  |  |  |
|----|----|----|
| **Score** | **Meaning** | **Action** |
| 3 | GWAS hit + high variance + diverse | Top priority for stacking |
| 2 (GWAS + var) | Real effect, low diversity | Select across populations |
| 2 (GWAS + div) | Real locus, small effect | Include if trait is oligogenic |
| 2 (var + div) | Variance explained, no GWAS | May be pop. structure – verify |
| 1 | Single evidence only | Use with caution |
| 0 | No evidence | Exclude from stacking |

## References

Difabachew YF et al. (2023). Genomic prediction with haplotype blocks in
wheat. *Frontiers in Plant Science* **14**:1168547.
[doi:10.3389/fpls.2023.1168547](https://doi.org/10.3389/fpls.2023.1168547)

Weber SE, Frisch M, Snowdon RJ, Voss-Fels KP (2023). Haplotype blocks
for genomic prediction: a comparative evaluation in multiple crop
datasets. *Frontiers in Plant Science* **14**:1217589.
[doi:10.3389/fpls.2023.1217589](https://doi.org/10.3389/fpls.2023.1217589)

Tong J et al. (2024). Stacking beneficial haplotypes from the Vavilov
wheat collection. *Theoretical and Applied Genetics* **137**:274.
[doi:10.1007/s00122-024-04784-w](https://doi.org/10.1007/s00122-024-04784-w)

Tong J et al. (2025). Haplotype stacking to improve stability of stripe
rust resistance in wheat. *Theoretical and Applied Genetics*
**138**:267.
[doi:10.1007/s00122-025-05045-0](https://doi.org/10.1007/s00122-025-05045-0)

## Examples

``` r
if (FALSE) { # \dontrun{
# 1. GWAS integration
gwas   <- read.csv("gwas_results.csv")   # SNP, CHR, POS, P, BETA
qtl    <- define_qtl_regions(gwas, blocks, snp_info, p_threshold = 5e-8)

# 2. Haplotype prediction
blues  <- read.csv("blues.csv")
pred   <- run_haplotype_prediction(geno, snp_info, blocks,
                                    blues    = blues,
                                    id_col   = "id",
                                    blue_col = "YLD")

# 3. Haplotype diversity
haps   <- extract_haplotypes(geno, snp_info, blocks, min_snps = 3)
div    <- compute_haplotype_diversity(haps)

# 4. Integrate all three
priority <- integrate_gwas_haplotypes(qtl, pred, diversity = div)

# Top priority blocks for haplotype stacking
priority[priority$priority_score == 3, ]
} # }
```
