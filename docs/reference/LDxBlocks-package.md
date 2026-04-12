# LDxBlocks: Genome-Wide LD Block Detection and Haplotype Analysis

Fast, scalable linkage disequilibrium (LD) block detection for
genome-wide SNP data using a C++/Armadillo computational core and two LD
metrics: standard r^2 and kinship-adjusted rV^2. Includes haplotype
reconstruction, diversity metrics, prediction feature matrices, and
multi-format genotype I/O.

## Choosing between r^2 and rV^2

- r^2 (`method = "r2"`):

  Standard squared Pearson correlation. No kinship matrix required. C++
  kernel via RcppArmadillo + OpenMP. Suitable for unstructured
  populations and large datasets (10 M+ markers).

- rV^2 (`method = "rV2"`):

  Kinship-whitened squared correlation (Mangin et al. 2012, Heredity
  108:285-291). Corrects LD inflation from relatedness. For livestock,
  inbred lines, or family-based cohorts with \< 200 k markers per
  chromosome.

## Pipeline overview

1.  Read genotype data:
    [`read_geno`](https://FAkohoue.github.io/LDxBlocks/reference/read_geno.md) -
    auto-detects CSV, HapMap, VCF, GDS, BED, or plain matrix.

2.  Detect LD blocks chromosome-wise:
    [`run_Big_LD_all_chr`](https://FAkohoue.github.io/LDxBlocks/reference/run_Big_LD_all_chr.md).

3.  Optionally auto-tune parameters:
    [`tune_LD_params`](https://FAkohoue.github.io/LDxBlocks/reference/tune_LD_params.md).

4.  Reconstruct haplotypes:
    [`extract_haplotypes`](https://FAkohoue.github.io/LDxBlocks/reference/extract_haplotypes.md).

5.  Compute diversity metrics:
    [`compute_haplotype_diversity`](https://FAkohoue.github.io/LDxBlocks/reference/compute_haplotype_diversity.md).

6.  Build prediction feature matrix:
    [`build_haplotype_feature_matrix`](https://FAkohoue.github.io/LDxBlocks/reference/build_haplotype_feature_matrix.md).

## Example data

- [`ldx_geno`](https://FAkohoue.github.io/LDxBlocks/reference/ldx_geno.md):

  120 x 200 genotype matrix, 3 chromosomes, 9 simulated LD blocks.

- [`ldx_snp_info`](https://FAkohoue.github.io/LDxBlocks/reference/ldx_snp_info.md):

  SNP metadata (SNP, CHR, POS, REF, ALT).

- [`ldx_blocks`](https://FAkohoue.github.io/LDxBlocks/reference/ldx_blocks.md):

  Reference block table for `ldx_geno`.

- [`ldx_gwas`](https://FAkohoue.github.io/LDxBlocks/reference/ldx_gwas.md):

  20 toy GWAS markers for tuning demos.

## References

Kim S-A et al. (2018) Bioinformatics 34(4):588-596.
[doi:10.1093/bioinformatics/btx609](https://doi.org/10.1093/bioinformatics/btx609)

Mangin B et al. (2012) Heredity 108(3):285-291.
[doi:10.1038/hdy.2011.73](https://doi.org/10.1038/hdy.2011.73)

VanRaden PM (2008) J. Dairy Sci. 91(11):4414-4423.
[doi:10.3168/jds.2007-0980](https://doi.org/10.3168/jds.2007-0980)

## See also

Useful links:

- <https://github.com/FAkohoue/LDxBlocks>

- <https://FAkohoue.github.io/LDxBlocks/>

- Report bugs at <https://github.com/FAkohoue/LDxBlocks/issues>

## Author

**Maintainer**: Félicien Akohoue <akohoue.f@gmail.com>
([ORCID](https://orcid.org/0000-0002-2160-0182))
