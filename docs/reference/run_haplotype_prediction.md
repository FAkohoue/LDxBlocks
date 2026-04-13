# Haplotype Prediction and Block Importance from Pre-Adjusted Phenotypes

Runs the complete Tong et al. (2025) haplotype stacking pipeline using
pre-adjusted phenotype values (BLUEs, BLUPs, or adjusted entry means).
Accepts either a single trait or multiple traits simultaneously.

When a single trait is supplied, GBLUP is fitted via
[`rrBLUP::kin.blup()`](https://rdrr.io/pkg/rrBLUP/man/kin.blup.html) and
block importance is ranked by `Var(local GEBV)` for that trait.

When multiple traits are supplied, a single trait-agnostic GRM is
computed once and shared across all traits.
[`rrBLUP::kin.blup()`](https://rdrr.io/pkg/rrBLUP/man/kin.blup.html) is
fitted per trait using this shared GRM. Block importance is summarised
across traits with per-trait columns plus cross-trait aggregates.

## Usage

``` r
run_haplotype_prediction(
  geno_matrix,
  snp_info,
  blocks,
  blues,
  id_col = "id",
  blue_col = "blue",
  blue_cols = NULL,
  importance_rule = c("any", "all", "mean"),
  top_n = NULL,
  min_freq = 0.01,
  min_snps = 3L,
  bend = TRUE,
  verbose = TRUE
)
```

## Arguments

- geno_matrix:

  Numeric matrix (individuals x SNPs), values 0/1/2/NA. Row names must
  be genotype IDs.

- snp_info:

  Data frame with columns `SNP`, `CHR`, `POS`.

- blocks:

  Block table from
  [`run_Big_LD_all_chr`](https://FAkohoue.github.io/LDxBlocks/reference/run_Big_LD_all_chr.md).

- blues:

  Pre-adjusted phenotype values. See section above for accepted formats.

- id_col:

  Name of the ID column when `blues` is a data frame. Default `"id"`.

- blue_col:

  Name of the BLUE column for single-trait data frames. Default
  `"blue"`.

- blue_cols:

  Character vector of trait column names for multi-trait data frames.
  Default `NULL` – all numeric non-ID columns are used.

- importance_rule:

  How to set the combined `important` flag for multi-trait results:

  - `"any"` (default): TRUE if important for \>= 1 trait.

  - `"all"`: TRUE only if important for all traits.

  - `"mean"`: TRUE if `var_scaled_mean` \>= 0.9.

  Ignored for single-trait runs (single-trait `important` is always
  scaled variance \>= 0.9).

- top_n:

  Integer or `NULL`. Max haplotype alleles per block. Default `NULL`.

- min_freq:

  Minimum haplotype allele frequency. Default `0.01`.

- min_snps:

  Minimum SNPs per block. Default `3L`.

- bend:

  Logical. Add ridge to GRM diagonal. Default `TRUE`.

- verbose:

  Logical. Print progress. Default `TRUE`.

## Value

Named list. For single-trait runs, contains:

- `blocks`, `diversity`, `hap_matrix`, `haplotypes`, `G`:

  Core pipeline outputs.

- `n_blocks`, `n_hap_columns`:

  Summary counts.

- `n_traits`:

  Integer: 1 for single-trait.

- `traits`:

  Character vector of trait names.

- `solver_used`:

  Character: always `"rrBLUP"`.

- `gebv`:

  Named numeric vector of GEBV.

- `snp_effects`:

  Per-SNP additive effects.

- `local_gebv`:

  Matrix (individuals x blocks).

- `block_importance`:

  Data frame with `block_id`, coordinates, `var_local_gebv`,
  `var_scaled`, `important`.

- `n_train`, `n_predict`:

  Training/prediction counts.

For multi-trait runs, the same structure but `gebv`, `snp_effects`,
`local_gebv` are named lists (one per trait), and `block_importance`
contains additional per-trait and cross-trait columns as described in
the *Block importance* section. `block_importance_list` is added as a
named list of single-trait block importance data frames.

## Input format for `blues`

Accepts any of three formats:

- A **named numeric vector** – single trait, names are genotype IDs:
  `c(G001 = 4.2, G002 = 3.8)`. `id_col`/`blue_col` are ignored.

- A **data frame with one trait column** – single trait, `id_col` and
  `blue_col` specify the columns.

- A **data frame with multiple trait columns** – multi-trait, `id_col`
  names the ID column, `blue_cols` names the trait columns (or `NULL` to
  use all numeric non-ID columns).

- A **named list of named numeric vectors** – multi-trait with
  potentially different individuals per trait:
  `list(YLD = c(G001=4.2, ...), DIS = c(G001=0.3, ...))`.

## Multi-trait GBLUP solver strategy

All traits are fitted using
[`rrBLUP::kin.blup()`](https://rdrr.io/pkg/rrBLUP/man/kin.blup.html) in
a per-trait loop sharing the same GRM. This produces numerically
identical results to
[`sommer::mmes()`](https://rdrr.io/pkg/sommer/man/mmes.html) for
single-trait models (verified empirically to 4 decimal places), while
avoiding the multi-trait [`cbind()`](https://rdrr.io/r/base/cbind.html)
formula which fails in sommer \<= 4.4.5 (Armadillo fixed-size matrix
error in the C++ coefficient matrix construction). Because all traits
use the same GRM, cross-trait block importance values are directly
comparable. `solver_used` is always `"rrBLUP"`.

## Block importance – single trait

`Var(local GEBV)` is computed per block after backsolving per-SNP
effects from GEBV. Blocks are scaled to \[0,1\]; those with scaled
variance \>= 0.9 are flagged `important`. This follows Tong et al.
(2024).

## Block importance – multi-trait

Per-trait columns `var_scaled_<trait>` and `important_<trait>` are added
for every trait. Cross-trait aggregates:

- `var_scaled_mean`:

  Mean scaled variance across all traits. Primary ranking criterion –
  blocks consistently important across traits are more robust stacking
  candidates.

- `n_traits_important`:

  Count of traits for which this block is flagged important.

- `important_any`:

  TRUE if important for at least one trait.

- `important_all`:

  TRUE if important for all traits.

- `important`:

  Controlled by `importance_rule`.

## References

Tong J et al. (2025). Haplotype stacking to improve stability of stripe
rust resistance in wheat. *Theoretical and Applied Genetics*
**138**:267.
[doi:10.1007/s00122-025-05045-0](https://doi.org/10.1007/s00122-025-05045-0)

Tong J et al. (2024). Stacking beneficial haplotypes from the Vavilov
wheat collection. *Theoretical and Applied Genetics* **137**:274.
[doi:10.1007/s00122-024-04784-w](https://doi.org/10.1007/s00122-024-04784-w)

Endelman JB (2011). Ridge regression and other kernels for genomic
selection with R package rrBLUP. *Plant Genome* **4**:250-255.
[doi:10.3835/plantgenome2011.08.0024](https://doi.org/10.3835/plantgenome2011.08.0024)

Covarrubias-Pazaran G (2016). Genome-assisted prediction of quantitative
traits using the R package sommer. *PLOS ONE* **11**:e0156744.
[doi:10.1371/journal.pone.0156744](https://doi.org/10.1371/journal.pone.0156744)

## Examples

``` r
if (FALSE) { # \dontrun{
library(LDxBlocks)
be     <- read_geno("mydata.vcf.gz")
blocks <- run_Big_LD_all_chr(be, CLQcut = 0.70)
geno   <- read_chunk(be, seq_len(be$n_snps))
rownames(geno) <- be$sample_ids
colnames(geno) <- be$snp_info$SNP

# -- Single trait: named vector ---------------------------------------------
blues_vec <- c(G001 = 4.21, G002 = 3.87, G003 = 5.14)
res <- run_haplotype_prediction(geno, be$snp_info, blocks, blues = blues_vec)
res$block_importance[res$block_importance$important, ]
sort(res$gebv, decreasing = TRUE)

# -- Single trait: data frame -----------------------------------------------
blues_df <- read.csv("blues.csv")   # columns: Genotype, YLD_BLUE
res <- run_haplotype_prediction(geno, be$snp_info, blocks,
                                 blues    = blues_df,
                                 id_col   = "Genotype",
                                 blue_col = "YLD_BLUE")

# -- Multi-trait: data frame ------------------------------------------------
blues_mt <- read.csv("blues_mt.csv")  # columns: id, YLD, DIS, PHT
res_mt <- run_haplotype_prediction(geno, be$snp_info, blocks,
                                    blues           = blues_mt,
                                    id_col          = "id",
                                    blue_cols       = c("YLD","DIS","PHT"),
                                    importance_rule = "any")
res_mt$n_traits        # 3
res_mt$solver_used     # "rrBLUP"
res_mt$block_importance[
  res_mt$block_importance$important_any,
  c("block_id","var_scaled_YLD","var_scaled_DIS","var_scaled_mean",
    "n_traits_important")]

# -- Multi-trait: named list (different individuals per trait) --------------
res_mt2 <- run_haplotype_prediction(geno, be$snp_info, blocks,
  blues = list(
    YLD = c(G001 = 4.2, G002 = 3.8),
    DIS = c(G001 = 0.3, G003 = 0.7)
  ))
} # }
```
