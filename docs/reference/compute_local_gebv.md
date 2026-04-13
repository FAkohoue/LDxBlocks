# Compute Local Haplotype GEBV per Block (Tong et al. 2025)

For each LD block, computes the local GEBV (haplotype effect) for every
individual by summing the per-SNP additive effects of the alleles they
carry within the block:

\$\$\text{local GEBV}\_f = \sum\_{t \in f} \left\[ h_t \alpha\_{1t} +
(1-h_t)\alpha\_{0t} \right\]\$\$

where \\h_t \in \\0, 0.5, 1\\\\ is the allele dosage at SNP \\t\\ scaled
to \[0,1\], \\\alpha\_{1t}\\ is the ALT allele effect, and
\\\alpha\_{0t} = -\alpha\_{1t}\\ is the REF allele effect (by symmetry
of the additive model).

Blocks are ranked by `Var(local GEBV)` – blocks with high variance
contribute strongly to trait differences among individuals and likely
harbour causal loci (Tong et al. 2025).

## Usage

``` r
compute_local_gebv(geno_matrix, snp_info, blocks, snp_effects, scale = TRUE)
```

## Arguments

- geno_matrix:

  Numeric matrix (individuals x SNPs), values 0/1/2/NA.

- snp_info:

  Data frame with columns `SNP`, `CHR`, `POS`.

- blocks:

  Block table from
  [`run_Big_LD_all_chr`](https://FAkohoue.github.io/LDxBlocks/reference/run_Big_LD_all_chr.md).

- snp_effects:

  Named numeric vector of per-SNP additive effects from
  [`backsolve_snp_effects`](https://FAkohoue.github.io/LDxBlocks/reference/backsolve_snp_effects.md)
  or from a marker model directly.

- scale:

  Logical. If `TRUE` (default), scale `Var(local GEBV)` to \[0,1\] so
  blocks are comparable across traits and datasets.

## Value

A list with two elements:

- `local_gebv`:

  Numeric matrix (individuals x blocks) of per-block local GEBV values.

- `block_importance`:

  Data frame with one row per block: block_id, CHR, start_bp, end_bp,
  n_snps, var_local_gebv, var_scaled, important (logical: scaled
  variance \>= 0.9).

## References

Tong J et al. (2025). Haplotype stacking to improve stability of stripe
rust resistance in wheat. *Theoretical and Applied Genetics*
**138**:267.
[doi:10.1007/s00122-025-05045-0](https://doi.org/10.1007/s00122-025-05045-0)

## Examples

``` r
if (FALSE) { # \dontrun{
snp_fx  <- backsolve_snp_effects(my_geno, gebv)
loc     <- compute_local_gebv(my_geno, snp_info, blocks, snp_fx)
# Blocks with scaled variance >= 0.9 are most important
head(loc$block_importance[loc$block_importance$important, ])
} # }
```
