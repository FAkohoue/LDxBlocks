# Fine-Map Epistatic SNP Pairs Within a Single Block

Exhaustively or adaptively identifies the specific SNP pairs within a
single LD block that drive the block-level epistatic signal. Three
methods are available:

- `"pairwise"`:

  Exhaustive C(p,2) scan. Tests all pairs using the same model as
  [`scan_block_epistasis`](https://FAkohoue.github.io/LDxBlocks/reference/scan_block_epistasis.md).
  Recommended for blocks with p \<= 200 SNPs.

- `"lasso"`:

  Fits a LASSO model with main effects and all pairwise interaction
  terms using `glmnet`. Identifies the interaction terms with non-zero
  coefficients at the cross-validated lambda. Recommended for blocks
  with p \> 200 SNPs.

- `"auto"`:

  Dispatches to `"pairwise"` when p \<= 200, `"lasso"` otherwise.

## Usage

``` r
fine_map_epistasis_block(
  block_id,
  geno_matrix,
  snp_info,
  blocks,
  y_resid,
  method = c("auto", "pairwise", "lasso"),
  min_freq = 0.05,
  sig_threshold = 0.05,
  sig_metric = c("p_simplem_sidak", "p_simplem", "p_bonf", "p_fdr"),
  meff_percent_cut = 0.995,
  lasso_nfolds = 5L,
  verbose = TRUE
)
```

## Arguments

- block_id:

  Character. The block to fine-map.

- geno_matrix:

  Numeric matrix (individuals x SNPs) or `LDxBlocks_backend`.

- snp_info:

  Data frame with columns `SNP`, `CHR`, `POS`.

- blocks:

  LD block table.

- y_resid:

  Named numeric vector of GRM-corrected REML residuals (from
  `test_block_haplotypes` null model). Names must match individual IDs.

- method:

  Character. One of `"pairwise"`, `"lasso"`, `"auto"`. Default `"auto"`.

- min_freq:

  Numeric. Minimum MAF. Default `0.05`.

- sig_threshold:

  Numeric. Significance threshold applied to the p-value chosen by
  `sig_metric`. Default `0.05`.

- sig_metric:

  Character. Which p-value drives the `significant` flag in pairwise
  output. One of `"p_simplem_sidak"` (default, recommended),
  `"p_simplem"`, `"p_bonf"`, or `"p_fdr"`. All four p-value columns are
  always present regardless of this choice. Ignored for the `"lasso"`
  method (which returns `lasso_coef` and `selected` instead of
  p-values).

- meff_percent_cut:

  Numeric in (0, 1). Variance threshold for simpleM eigendecomposition.
  Default `0.995`.

- lasso_nfolds:

  Integer. CV folds for glmnet lambda selection. Default `5L`.

- verbose:

  Logical. Default `TRUE`.

## Value

A data frame sorted by p_wald (pairwise) or \|coefficient\| (lasso).
Columns: `SNP_i`, `SNP_j`, `POS_i`, `POS_j`, `dist_bp`, `aa_effect` or
`lasso_coef`, `SE`, `t_stat`, `p_wald`, `p_bonf`, `significant`
(pairwise) or `selected` (lasso).

## See also

[`scan_block_epistasis`](https://FAkohoue.github.io/LDxBlocks/reference/scan_block_epistasis.md),
[`scan_block_by_block_epistasis`](https://FAkohoue.github.io/LDxBlocks/reference/scan_block_by_block_epistasis.md)
