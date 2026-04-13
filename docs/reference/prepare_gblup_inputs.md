# Prepare Genomic Prediction Inputs for External GBLUP Software

Assembles the inputs required to fit a GBLUP model in external software
(rrBLUP, sommer, ASReml-R, BGLR) and subsequently run the Tong et al.
(2025) haplotype stacking pipeline. Returns the VanRaden GRM computed
from the haplotype feature matrix, aligned with a user-supplied
phenotype table.

## Usage

``` r
prepare_gblup_inputs(
  hap_matrix,
  pheno_df,
  id_col = "id",
  trait_col = NULL,
  bend = TRUE
)
```

## Arguments

- hap_matrix:

  Numeric matrix (individuals x haplotype alleles) from
  [`build_haplotype_feature_matrix`](https://FAkohoue.github.io/LDxBlocks/reference/build_haplotype_feature_matrix.md).

- pheno_df:

  Data frame of phenotypes. Must contain:

  - An ID column (set via `id_col`, default `"id"`). Values must match
    `rownames(hap_matrix)` exactly (case-sensitive, no leading/trailing
    spaces).

  - One or more numeric trait columns (referenced via `trait_col`).
    Column names are arbitrary.

  - `NA` values in trait columns are allowed.

  Minimal single-trait format:

  |      |      |
  |------|------|
  | id   | YLD  |
  | G001 | 4.21 |
  | G002 | 3.87 |
  | G003 | NA   |

- id_col:

  Name of the individual ID column in `pheno_df`. Default `"id"`.

- trait_col:

  Name of the trait column to extract as a numeric vector. Default
  `NULL` – no `y_vec` is returned, only the aligned data frame.

- bend:

  Logical. Add 0.001 to diagonal of G for positive-definiteness. Default
  `TRUE` (recommended for mixed model solvers).

## Value

A named list:

- `G`:

  VanRaden GRM (n x n), aligned to individuals present in both
  `hap_matrix` and `pheno_df`.

- `pheno_df`:

  Phenotype data frame subset and reordered to match rows of `G`.

- `y_vec`:

  Named numeric vector of the requested trait (only if `trait_col` is
  supplied). `NA` values are preserved so the user can decide how to
  handle them (e.g. set to `NA` for prediction-only individuals in a
  training/validation split).

- `n_train`:

  Number of individuals with non-missing trait values.

- `n_predict`:

  Number of individuals with missing trait values (prediction
  candidates).

## Workflow

LDxBlocks handles genotype processing and block detection. Phenotype
handling and GBLUP fitting are intentionally left to dedicated R
packages because phenotype data requires preprocessing
(multi-environment adjustment, outlier removal, covariate inclusion)
that is dataset-specific. The handoff is:

1.  **LDxBlocks** (this function): produce aligned G matrix and
    phenotype vector.

2.  **External GBLUP** (rrBLUP / sommer / ASReml-R / BGLR): fit the
    model, obtain GEBV.

3.  **LDxBlocks**
    ([`backsolve_snp_effects`](https://FAkohoue.github.io/LDxBlocks/reference/backsolve_snp_effects.md) +
    [`compute_local_gebv`](https://FAkohoue.github.io/LDxBlocks/reference/compute_local_gebv.md)):
    derive block-level haplotype effects from the GEBV.

## Example GBLUP calls after this function

    # rrBLUP
    library(rrBLUP)
    fit  <- kin.blup(data = inp$pheno_df, geno = "id",
                     pheno = "trait", K = inp$G)
    gebv <- fit$g

    # sommer
    library(sommer)
    fit  <- sommer::mmes(trait ~ 1, random = ~vsm(ism(id), Gu = inp$G),
                 data = inp$pheno_df)
    gebv <- fit$U$`u:id`$trait

    # BGLR
    library(BGLR)
    fit  <- BGLR(y = inp$y_vec, ETA = list(list(K = inp$G, model = "RKHS")))
    gebv <- fit$yHat

## References

Tong J et al. (2025). Haplotype stacking to improve stability of stripe
rust resistance in wheat. *Theoretical and Applied Genetics*
**138**:267.
[doi:10.1007/s00122-025-05045-0](https://doi.org/10.1007/s00122-025-05045-0)

VanRaden PM (2008). Efficient methods to compute genomic predictions.
*Journal of Dairy Science* **91**(11):4414-4423.
[doi:10.3168/jds.2007-0980](https://doi.org/10.3168/jds.2007-0980)

## Examples

``` r
if (FALSE) { # \dontrun{
# After building haplotype feature matrix:
feat  <- build_haplotype_feature_matrix(haps, top_n = 5)
pheno <- read.csv("phenotypes.csv")   # columns: id, YLD, PHT, ...

inp <- prepare_gblup_inputs(feat, pheno, id_col = "id",
                             trait_col = "YLD")

# Fit GBLUP with rrBLUP
library(rrBLUP)
fit  <- kin.blup(data = inp$pheno_df, geno = "id",
                 pheno = "YLD", K = inp$G)
gebv <- fit$g

# Then derive block-level haplotype effects
snp_fx <- backsolve_snp_effects(geno_matrix, gebv, G = inp$G)
loc    <- compute_local_gebv(geno_matrix, snp_info, blocks, snp_fx)
} # }
```
