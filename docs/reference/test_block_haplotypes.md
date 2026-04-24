# Block-Level Haplotype Association Testing (Q+K Mixed Linear Model with simpleM Multiple-Testing Correction)

Performs genome-wide haplotype block association tests for one or more
quantitative traits. Each LD block is tested as a unit: per-allele Wald
tests identify which specific haplotype alleles drive association, and
an omnibus F-test evaluates the block as a whole. Population structure
and kinship are corrected jointly through a unified mixed linear model.

In addition to raw and FDR-adjusted p-values, this function applies
**simpleM**-style multiple-testing correction to account for correlation
among tested haplotype alleles. simpleM estimates the effective number
of independent tests (\\M\_{\mathrm{eff}}\\) from the eigenvalues of the
haplotype allele correlation matrix and derives either a
Bonferroni-style or Sidak-style adjusted threshold.

**Statistical model (Q+K / EMMAX formulation):**

\$\$y = \mu + \alpha \cdot x\_{\mathrm{hap}} + \sum\_{k=1}^{K} \beta_k
\\ PC_k + g + \varepsilon\$\$

- \\x\_{\mathrm{hap}}\\:

  Haplotype allele dosage (0, 1, or 2 copies for phased data; 0 or 1 for
  unphased).

- \\PC_k\\:

  k-th eigenvector of the haplotype GRM, as a fixed-effect covariate for
  population structure.

- \\g \sim MVN(0,\\\sigma_g^2 G)\\:

  Polygenic background as a random effect.

- \\\varepsilon \sim MVN(0,\\\sigma_e^2 I)\\:

  Residual error.

**simpleM extension (Gao et al. 2008, 2010, 2011):**
\\M\_{\mathrm{eff}}\\ is estimated from the eigenspectrum of the
haplotype allele dosage correlation matrix - the number of principal
components needed to explain `meff_percent_cut` of variance (default
99.5%). Two simpleM adjusted p-value paths are always computed:

- `p_simplem` - Bonferroni-style: \\\min(p \times
  M\_{\mathrm{eff}},\\1)\\.

- `p_simplem_sidak` - Sidak-style: \\1 - (1 - p)^{M\_{\mathrm{eff}}}\\.

Because haplotype alleles are correlated within and across nearby
blocks, simpleM is less conservative than raw Bonferroni while still
providing family-wise error control.

## Usage

``` r
test_block_haplotypes(
  haplotypes,
  blues,
  blocks,
  n_pcs = 0L,
  top_n = NULL,
  min_freq = 0.05,
  id_col = "id",
  blue_col = "blue",
  blue_cols = NULL,
  sig_threshold = 0.05,
  sig_metric = c("p_wald", "p_fdr", "p_simplem", "p_simplem_sidak"),
  meff_scope = c("chromosome", "global", "block"),
  meff_percent_cut = 0.995,
  meff_max_cols = 1000L,
  plot = FALSE,
  out_dir = ".",
  verbose = TRUE
)
```

## Arguments

- haplotypes:

  Named list produced by
  [`extract_haplotypes`](https://FAkohoue.github.io/LDxBlocks/reference/extract_haplotypes.md).
  Each element is a named character vector (one haplotype dosage string
  per individual), and the list must carry a `block_info` attribute.

- blues:

  Pre-adjusted phenotype means. Accepted in four formats: named numeric
  vector, single-trait data frame, multi-trait data frame, or named list
  of named numeric vectors. At least 10 common individuals are required
  per trait.

- blocks:

  LD block table from
  [`run_Big_LD_all_chr`](https://FAkohoue.github.io/LDxBlocks/reference/run_Big_LD_all_chr.md).
  Used for block metadata annotation in the output only.

- n_pcs:

  Integer (`>= 0`) or `NULL`. Number of haplotype-GRM eigenvectors to
  include as fixed-effect population structure covariates. `0L`
  (default): pure GRM / EMMAX. `1-10`: Q+K model. `NULL`: auto-selected
  from the GRM scree-plot elbow, capped at 10.

- top_n:

  Integer or `NULL`. Maximum haplotype alleles per block in the feature
  matrix. `NULL` retains all alleles above `min_freq`.

- min_freq:

  Numeric in (0, 1). Minimum haplotype allele frequency. Default `0.05`.

- id_col:

  Character. Individual-ID column name when `blues` is a data frame.
  Default `"id"`.

- blue_col:

  Character. Phenotype column name for single-trait data frames. Default
  `"blue"`.

- blue_cols:

  Character vector. Phenotype column names for multi-trait data frames.
  Default `NULL`.

- sig_threshold:

  Numeric in (0, 1\]. Significance cutoff applied to the p-value chosen
  by `sig_metric`. Default `0.05`. Common choices:

  - `0.05` - recommended with `sig_metric = "p_fdr"`, `"p_simplem"`, or
    `"p_simplem_sidak"`.

  - `0.05 / n_tests` - raw Bonferroni with `sig_metric = "p_wald"`.

- sig_metric:

  Character. Which p-value to use for the `significant` and
  `significant_omnibus` flags. One of:

  - `"p_wald"` - raw Wald p-value. Use with a pre-corrected
    `sig_threshold` (e.g. `0.05 / n_tests`).

  - `"p_fdr"` - Benjamini-Hochberg FDR. Recommended for
    discovery-oriented analyses.

  - `"p_simplem"` - simpleM Bonferroni-style adjusted p-value.

  - `"p_simplem_sidak"` - simpleM Sidak-style adjusted p-value.
    Recommended default when family-wise error control is desired with
    correlated haplotype predictors.

  All four p-value columns are always present in the output regardless
  of this choice.

- meff_scope:

  Character. Scope at which \\M\_{\mathrm{eff}}\\ is estimated. One of:

  - `"chromosome"` (recommended) - separate \\M\_{\mathrm{eff}}\\ per
    chromosome; best balance of scalability and LD awareness.

  - `"global"` - one genome-wide \\M\_{\mathrm{eff}}\\; suitable for
    moderate-sized scans.

  - `"block"` - one \\M\_{\mathrm{eff}}\\ per LD block at the allele
    level; for block omnibus tests \\M\_{\mathrm{eff}} = 1\\.

- meff_percent_cut:

  Numeric in (0, 1). Proportion of variance explained by the retained
  PCs in the simpleM eigendecomposition. Default `0.995` (99.5%),
  following the original simpleM recommendation.

- meff_max_cols:

  Integer. Maximum columns per eigendecomposition block when computing
  \\M\_{\mathrm{eff}}\\. Larger groups are chunked and summed. Default
  `1000L`.

- plot:

  Logical. If `TRUE`, save Manhattan and Q-Q plots. Default `FALSE`.
  Significance markers reflect `sig_metric`.

- out_dir:

  Character. Directory for optional plots. Default `"."`.

- verbose:

  Logical. Print progress messages. Default `TRUE`.

## Value

A named list of class `c("LDxBlocks_haplotype_assoc", "list")`:

- `allele_tests`:

  Data frame with one row per allele per block per trait. Always
  includes: `p_wald`, `p_fdr`, `Meff`, `alpha_simplem`,
  `alpha_simplem_sidak`, `p_simplem`, `p_simplem_sidak`, `significant`.

- `block_tests`:

  Data frame with one row per block per trait. Always includes:
  `p_omnibus`, `p_omnibus_fdr`, `p_omnibus_adj` (backward-compat
  Bonferroni), `Meff`, `alpha_simplem`, `alpha_simplem_sidak`,
  `p_omnibus_simplem`, `p_omnibus_simplem_sidak`, `significant_omnibus`.

- `traits`:

  Character vector of trait names tested.

- `n_pcs_used`:

  Integer. GRM PCs used in the null model.

- `sig_threshold`:

  Numeric. Significance cutoff used.

- `sig_metric`:

  Character. P-value used for significance flags.

- `n_tests`:

  Integer. Total allele-level tests performed.

- `meff_scope`:

  Character. simpleM estimation scope.

- `meff_percent_cut`:

  Numeric. Variance cutoff for simpleM.

- `meff`:

  Named list of \\M\_{\mathrm{eff}}\\ summaries per trait, each with
  `$allele` (global/chromosome/block) and `$block` (global/chromosome)
  components.

## References

Endelman JB (2011). Ridge regression and other kernels for genomic
selection with R package rrBLUP. *Plant Genome* **4**:250-255.
[doi:10.3835/plantgenome2011.08.0024](https://doi.org/10.3835/plantgenome2011.08.0024)

Yu J et al. (2006). A unified mixed-model method for association
mapping. *Nature Genetics* **38**(2):203-208.
[doi:10.1038/ng1702](https://doi.org/10.1038/ng1702)

Kang HM et al. (2010). Variance component model to account for sample
structure in genome-wide association studies. *Nature Genetics*
**42**(4):348-354. [doi:10.1038/ng.548](https://doi.org/10.1038/ng.548)

Gao X et al. (2008). A multiple testing correction method for genetic
association studies using correlated single nucleotide polymorphisms.
*Genetic Epidemiology* **32**:361-369.
[doi:10.1002/gepi.20310](https://doi.org/10.1002/gepi.20310)

Gao X et al. (2010). Avoiding the high Bonferroni penalty in genome-wide
association studies. *Genetic Epidemiology* **34**:100-105.
[doi:10.1002/gepi.20430](https://doi.org/10.1002/gepi.20430)

Gao X (2011). Multiple testing corrections for imputed SNPs. *Genetic
Epidemiology* **35**:154-158.
[doi:10.1002/gepi.20563](https://doi.org/10.1002/gepi.20563)

## See also

[`extract_haplotypes`](https://FAkohoue.github.io/LDxBlocks/reference/extract_haplotypes.md),
[`compute_haplotype_grm`](https://FAkohoue.github.io/LDxBlocks/reference/compute_haplotype_grm.md),
[`run_haplotype_prediction`](https://FAkohoue.github.io/LDxBlocks/reference/run_haplotype_prediction.md),
[`estimate_diplotype_effects`](https://FAkohoue.github.io/LDxBlocks/reference/estimate_diplotype_effects.md)

## Examples

``` r
# \donttest{
data(ldx_geno, ldx_snp_info, ldx_blocks, ldx_blues, package = "LDxBlocks")
haps <- extract_haplotypes(ldx_geno, ldx_snp_info, ldx_blocks)

# FDR-based discovery (default EMMAX model)
assoc_fdr <- test_block_haplotypes(
  haplotypes = haps, blues = setNames(ldx_blues$YLD, ldx_blues$id),
  blocks = ldx_blocks, sig_metric = "p_fdr", verbose = FALSE
)

# simpleM Sidak with 3 GRM PCs (Q+K), chromosome-wise Meff
assoc_sm <- test_block_haplotypes(
  haplotypes = haps, blues = setNames(ldx_blues$YLD, ldx_blues$id),
  blocks = ldx_blocks, n_pcs = 3L,
  sig_metric = "p_simplem_sidak", meff_scope = "chromosome",
  verbose = FALSE
)

# Check per-trait Meff
assoc_sm$meff$trait$allele$chromosome
#>  1  2  3 
#> 23 23 24 

# Multi-trait with simpleM
assoc_mt <- test_block_haplotypes(
  haplotypes = haps, blues = ldx_blues,
  blocks = ldx_blocks, id_col = "id", blue_cols = c("YLD","RES"),
  sig_metric = "p_simplem", meff_scope = "chromosome", verbose = FALSE
)
# }
```
