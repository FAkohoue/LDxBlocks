# Map GWAS Hits to LD Blocks (Post-GWAS QTL Region Definition)

Maps significant GWAS markers onto LD blocks to define QTL regions.
Blocks with significant markers from multiple traits are flagged
pleiotropic. Implements the approach of Tong et al. (2024).

## Usage

``` r
define_qtl_regions(
  gwas_results,
  blocks,
  snp_info,
  p_threshold = 5e-08,
  trait_col = "trait",
  min_snps = 3L,
  ld_decay = NULL,
  verbose = TRUE
)
```

## Arguments

- gwas_results:

  Data frame with columns `SNP` (or `Marker`), `CHR`, `POS`. Optional
  columns: `P` (p-value), `BETA` (additive effect estimate), `trait`.

- blocks:

  LD block data frame from
  [`run_Big_LD_all_chr`](https://FAkohoue.github.io/LDxBlocks/reference/run_Big_LD_all_chr.md).

- snp_info:

  Full SNP metadata data frame.

- p_threshold:

  Significance threshold. Default `5e-8`. `NULL` uses all markers
  regardless of p-value.

- trait_col:

  Trait column name in `gwas_results`. Default `"trait"`.

- min_snps:

  Minimum SNPs per block. Default `3L`.

- ld_decay:

  Optional `LDxBlocks_decay` object from
  [`compute_ld_decay`](https://FAkohoue.github.io/LDxBlocks/reference/compute_ld_decay.md),
  or a data frame with columns `CHR` and `decay_dist_bp`. When supplied,
  candidate gene windows are extended by the chromosome-specific decay
  distance on both sides of each lead SNP position, adding columns
  `candidate_region_start`, `candidate_region_end`, and
  `candidate_region_size_kb` to the output. Default `NULL`.

- verbose:

  Logical. Print progress messages. Default `TRUE`.

## Value

Data frame with one row per block containing significant markers:
`block_id`, `CHR`, `start_bp`, `end_bp`, `n_snps_block`,
`n_sig_markers`, `lead_marker`, `lead_p`, `lead_beta` (if BETA
supplied), `sig_markers`, `sig_betas` (if BETA supplied), `traits`,
`n_traits`, `pleiotropic`.

## Block effect estimation

When multiple SNPs within a block are GWAS-significant, their marginal
BETA values are correlated (due to LD) and cannot be summed directly.
LDxBlocks returns three complementary columns when `BETA` is present in
`gwas_results`:

- `lead_beta`:

  BETA of the lead (lowest-p) SNP. The simplest proxy for the block
  effect – assumes the lead SNP fully tags the causal variant.
  Underestimates the block effect when multiple independent signals
  exist.

- `sig_markers`:

  Semicolon-separated IDs of all significant SNPs in the block. Pass
  these to a conditional/joint analysis tool (e.g. COJO, SuSiE, finemap)
  to obtain independent within-block effects.

- `sig_betas`:

  Semicolon-separated marginal BETA values for all significant SNPs
  (same order as `sig_markers`). Their absolute values are an upper
  bound on the true block effect because they include LD-induced
  inflation; use `lead_beta` or joint analysis for calibrated estimates.

For the biologically correct block-level effect, fit a haplotype model:
[`build_haplotype_feature_matrix()`](https://FAkohoue.github.io/LDxBlocks/reference/build_haplotype_feature_matrix.md)
produces the 0/1/2 dosage columns whose regression coefficients
(alpha_h) are the true per-allele effects.

## References

Tong J, Tarekegn ZT, Jambuthenne D, Alahmad S, Periyannan S, Hickey L,
Dinglasan E, Hayes B (2024). Stacking beneficial haplotypes from the
Vavilov wheat collection to accelerate breeding for multiple disease
resistance. *Theoretical and Applied Genetics* **137**:274.
[doi:10.1007/s00122-024-04784-w](https://doi.org/10.1007/s00122-024-04784-w)

Yang J et al. (2012). Conditional and joint multiple-SNP analysis of
GWAS summary statistics identifies additional variants influencing
complex traits. *Nature Genetics* **44**(4):369-375.
[doi:10.1038/ng.2213](https://doi.org/10.1038/ng.2213)
