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
  min_snps = 3L
)
```

## Arguments

- gwas_results:

  Data frame: SNP, CHR, POS. Optional: P, trait.

- blocks:

  LD block data frame from run_Big_LD_all_chr().

- snp_info:

  Full SNP metadata data frame.

- p_threshold:

  Significance threshold. Default 5e-8. NULL = use all markers.

- trait_col:

  Trait column name. Default "trait".

- min_snps:

  Minimum block SNP count. Default 3L.

## Value

Data frame: block_id, CHR, start_bp, end_bp, n_snps_block,
n_sig_markers, lead_snp, lead_p, traits, n_traits, pleiotropic.

## References

Tong et al. (2024) Theor Appl Genet 137:274.
