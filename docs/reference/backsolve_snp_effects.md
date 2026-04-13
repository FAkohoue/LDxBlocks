# Backsolve SNP Effects from GEBV (Tong et al. 2025)

Derives per-SNP additive effect estimates from genome-wide genomic
estimated breeding values (GEBV) without re-fitting a marker model. This
implements Step 2 of the Tong et al. (2025) haplotype stacking pipeline:

\$\$\hat{\alpha} = \frac{M^\top G^{-1} \hat{g}}{2 \sum_t p_t(1-p_t)}\$\$

where \\M\\ is the centred genotype matrix, \\G\\ is the VanRaden GRM,
\\\hat{g}\\ are the GEBV, and \\p_t\\ are allele frequencies.

## Usage

``` r
backsolve_snp_effects(geno_matrix, gebv, G = NULL)
```

## Arguments

- geno_matrix:

  Numeric matrix (individuals x SNPs), values 0/1/2/NA. Row names must
  match names of `gebv`.

- gebv:

  Named numeric vector of GEBV, one per individual. Names must match
  `rownames(geno_matrix)`.

- G:

  Optional pre-computed VanRaden GRM (n x n). If `NULL` (default),
  computed internally via
  [`compute_haplotype_grm()`](https://FAkohoue.github.io/LDxBlocks/reference/compute_haplotype_grm.md)-style
  logic. Supply your own if you used a bended or tuned G in the GBLUP.

## Value

Named numeric vector of length p (SNPs), one effect per SNP. Names match
`colnames(geno_matrix)`.

## Details

This approach is preferred over direct marker-effect estimation when
GEBV are already available from a GBLUP run (e.g. from ASReml-R, sommer,
or rrBLUP). It avoids refitting the marker model and produces marker
effects on the same scale as the original GEBV.

Missing genotype values are mean-imputed per column before computation.

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
# After fitting GBLUP with rrBLUP:
# fit  <- rrBLUP::kin.blup(data, geno="id", pheno="trait", K=G)
# gebv <- fit$g
snp_fx <- backsolve_snp_effects(geno_matrix = my_geno, gebv = gebv)
head(sort(abs(snp_fx), decreasing = TRUE))
} # }
```
