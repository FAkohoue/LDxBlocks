# Haplotype Prediction and Block Importance from Pre-Adjusted Phenotypes

Runs the complete Tong et al. (2025) haplotype stacking pipeline using
pre-adjusted phenotype values (BLUEs, BLUPs, or adjusted entry means).
Requires the rrBLUP package for REML-based GBLUP fitting.

## Usage

``` r
run_haplotype_prediction(
  geno_matrix,
  snp_info,
  blocks,
  blues,
  id_col = "id",
  blue_col = "blue",
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

  Pre-adjusted phenotype values. Either:

  - A **named numeric vector** (names = genotype IDs), or

  - A **data frame** with ID and BLUE columns specified by `id_col` and
    `blue_col`.

  Genotypes in `geno_matrix` but absent from `blues` are treated as
  prediction candidates.

- id_col:

  Name of the ID column when `blues` is a data frame. Default `"id"`.

- blue_col:

  Name of the BLUE column when `blues` is a data frame. Default
  `"blue"`.

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

Named list containing all outputs needed for downstream analysis,
matching the structure of
[`run_ldx_pipeline`](https://FAkohoue.github.io/LDxBlocks/reference/run_ldx_pipeline.md)
plus prediction outputs:

- `blocks`:

  LD block table from `run_Big_LD_all_chr`.

- `diversity`:

  Per-block diversity metrics (He, Shannon, n_eff_alleles,
  freq_dominant, sweep_flag).

- `hap_matrix`:

  Numeric matrix (individuals x haplotype allele columns) for genomic
  prediction — the primary output for dimensionality-reduced GP.

- `haplotypes`:

  Raw dosage strings from `extract_haplotypes` — pass to
  `decode_haplotype_strings` for nucleotide sequences.

- `snp_info_filtered`:

  SNP metadata after MAF filter.

- `n_blocks`, `n_hap_columns`:

  Summary counts.

- `gebv`:

  Named numeric vector of GEBV for all genotyped individuals (training +
  prediction candidates).

- `snp_effects`:

  Per-SNP additive effects backsolved from GEBV.

- `local_gebv`:

  Matrix (individuals x blocks) of per-block local GEBV values —
  quantifies each block's contribution to the trait.

- `block_importance`:

  Data frame ranking blocks by `Var(local GEBV)`: block_id, CHR,
  start_bp, end_bp, n_snps, var_local_gebv, var_scaled, important
  (scaled \>= 0.9).

- `G`:

  VanRaden GRM used in the GBLUP model.

- `n_train`, `n_predict`:

  Training and prediction individual counts.

## Input format for `blues`

Accepts either:

- A **named numeric vector** where names are genotype IDs:
  `c(G001 = 4.2, G002 = 3.8, ...)`.

- A **data frame** with one column of genotype IDs and one column of
  BLUE values. Column names are specified via `id_col` and `blue_col`.
  This is the natural format of field trial output from ASReml-R, lme4,
  or SpATS.

## Pipeline steps

1.  Extract haplotypes from LD blocks.

2.  Build haplotype feature matrix.

3.  Compute VanRaden GRM from haplotype features.

4.  Fit GBLUP via
    [`rrBLUP::kin.blup()`](https://rdrr.io/pkg/rrBLUP/man/kin.blup.html)
    (REML-based).

5.  Backsolve per-SNP additive effects from GEBV.

6.  Compute local haplotype GEBV per block.

7.  Rank blocks by `Var(local GEBV)`.

## References

Tong J et al. (2025). Haplotype stacking to improve stability of stripe
rust resistance in wheat. *Theoretical and Applied Genetics*
**138**:267.
[doi:10.1007/s00122-025-05045-0](https://doi.org/10.1007/s00122-025-05045-0)

Endelman JB (2011). Ridge regression and other kernels for genomic
selection with R package rrBLUP. *Plant Genome* **4**:250-255.
[doi:10.3835/plantgenome2011.08.0024](https://doi.org/10.3835/plantgenome2011.08.0024)

## Examples

``` r
if (FALSE) { # \dontrun{
library(LDxBlocks)

be     <- read_geno("mydata.vcf.gz")
geno   <- read_chunk(be, seq_len(be$n_snps))
rownames(geno) <- be$sample_ids
colnames(geno) <- be$snp_info$SNP
blocks <- run_Big_LD_all_chr(be, CLQcut = 0.5)

# Option 1: named numeric vector (same format as GWAS input)
blues_vec <- c(G001 = 4.21, G002 = 3.87, G003 = 5.14)
res <- run_haplotype_prediction(geno, be$snp_info, blocks,
                                 blues = blues_vec)

# Option 2: data frame from field trial analysis
blues_df <- read.csv("blues.csv")   # columns: Genotype, YLD_BLUE
res <- run_haplotype_prediction(geno, be$snp_info, blocks,
                                 blues    = blues_df,
                                 id_col   = "Genotype",
                                 blue_col = "YLD_BLUE")

res$block_importance[res$block_importance$important, ]
sort(res$gebv, decreasing = TRUE)
} # }
```
