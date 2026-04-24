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
    [`read_geno`](https://FAkohoue.github.io/LDxBlocks/reference/read_geno.md)
    – auto-detects CSV, HapMap, VCF, GDS, BED, or plain matrix. For WGS
    panels where peak RAM is a concern, use
    [`read_geno_bigmemory`](https://FAkohoue.github.io/LDxBlocks/reference/read_geno_bigmemory.md)
    to build a memory-mapped store (requires bigmemory).

2.  Detect LD blocks chromosome-wise:
    [`run_Big_LD_all_chr`](https://FAkohoue.github.io/LDxBlocks/reference/run_Big_LD_all_chr.md).

3.  Optionally auto-tune parameters:
    [`tune_LD_params`](https://FAkohoue.github.io/LDxBlocks/reference/tune_LD_params.md).

4.  Reconstruct haplotypes:
    [`extract_haplotypes`](https://FAkohoue.github.io/LDxBlocks/reference/extract_haplotypes.md).

5.  Decode haplotype strings to nucleotides:
    [`decode_haplotype_strings`](https://FAkohoue.github.io/LDxBlocks/reference/decode_haplotype_strings.md).

6.  Compute diversity metrics:
    [`compute_haplotype_diversity`](https://FAkohoue.github.io/LDxBlocks/reference/compute_haplotype_diversity.md).

7.  Build prediction feature matrix:
    [`build_haplotype_feature_matrix`](https://FAkohoue.github.io/LDxBlocks/reference/build_haplotype_feature_matrix.md).

8.  Compute LD decay and chromosome-specific decay distances:
    [`compute_ld_decay`](https://FAkohoue.github.io/LDxBlocks/reference/compute_ld_decay.md).
    Provides the critical r\\^2\\ threshold (parametric: 95th percentile
    of unlinked-marker r\\^2\\) and per-chromosome decay distances for
    candidate gene windows.

9.  Map GWAS hits to QTL regions (LD-aware windows when `ld_decay` is
    supplied):
    [`define_qtl_regions`](https://FAkohoue.github.io/LDxBlocks/reference/define_qtl_regions.md).

10. Genomic prediction (Tong et al. 2024/2025):
    [`run_haplotype_prediction`](https://FAkohoue.github.io/LDxBlocks/reference/run_haplotype_prediction.md)
    – full pipeline from BLUEs to block importance; or step-by-step via
    [`compute_haplotype_grm`](https://FAkohoue.github.io/LDxBlocks/reference/compute_haplotype_grm.md),
    [`backsolve_snp_effects`](https://FAkohoue.github.io/LDxBlocks/reference/backsolve_snp_effects.md),
    [`compute_local_gebv`](https://FAkohoue.github.io/LDxBlocks/reference/compute_local_gebv.md),
    [`rank_haplotype_blocks`](https://FAkohoue.github.io/LDxBlocks/reference/rank_haplotype_blocks.md).

11. Haplotype association testing with simpleM correction (Gao et al.
    2008, 2010, 2011):
    [`test_block_haplotypes`](https://FAkohoue.github.io/LDxBlocks/reference/test_block_haplotypes.md).
    `sig_metric` selects among `p_wald`, `p_fdr`, `p_simplem`
    (Bonferroni-style), and `p_simplem_sidak` (Sidak-style,
    recommended). All four p-value columns always present. `meff_scope`
    controls Meff estimation scope.

12. Cross-population GWAS validation - two functions depending on how
    the association analysis was performed:

    - [`compare_block_effects`](https://FAkohoue.github.io/LDxBlocks/reference/compare_block_effects.md) -
      for results from
      [`test_block_haplotypes`](https://FAkohoue.github.io/LDxBlocks/reference/test_block_haplotypes.md).
      IVW meta-analysis, Cochran Q, I-squared, direction agreement per
      block.

    - [`compare_gwas_effects`](https://FAkohoue.github.io/LDxBlocks/reference/compare_gwas_effects.md) -
      for results from external GWAS tools (GAPIT, TASSEL, FarmCPU,
      PLINK, etc.). Accepts raw GWAS data frames or
      [`define_qtl_regions`](https://FAkohoue.github.io/LDxBlocks/reference/define_qtl_regions.md)
      output. Derives SE from z-score when absent.

    - Both functions support `block_match = "position"` to match LD
      blocks by genomic interval overlap (Intersection-over-Union)
      rather than `block_id` string equality. This handles the common
      case where block boundaries differ between populations due to
      different ancestral LD structures. `overlap_min` (default 0.50)
      sets the minimum IoU for a valid positional match. The
      `match_type` output column reports `"exact"`, `"position"`, or
      `"pop1_only"` for every block.

13. Write outputs:
    [`write_haplotype_numeric`](https://FAkohoue.github.io/LDxBlocks/reference/write_haplotype_numeric.md),
    [`write_haplotype_character`](https://FAkohoue.github.io/LDxBlocks/reference/write_haplotype_character.md),
    [`write_haplotype_diversity`](https://FAkohoue.github.io/LDxBlocks/reference/write_haplotype_diversity.md).

14. Or run everything at once:
    [`run_ldx_pipeline`](https://FAkohoue.github.io/LDxBlocks/reference/run_ldx_pipeline.md).

## Example data

- [`ldx_geno`](https://FAkohoue.github.io/LDxBlocks/reference/ldx_geno.md):

  120 x 230 genotype matrix, 3 chromosomes, 9 simulated LD blocks.

- [`ldx_snp_info`](https://FAkohoue.github.io/LDxBlocks/reference/ldx_snp_info.md):

  SNP metadata (SNP, CHR, POS, REF, ALT).

- [`ldx_blocks`](https://FAkohoue.github.io/LDxBlocks/reference/ldx_blocks.md):

  Reference block table for `ldx_geno`.

- [`ldx_gwas`](https://FAkohoue.github.io/LDxBlocks/reference/ldx_gwas.md):

  20 toy GWAS markers for tuning demos.

- [`ldx_blues`](https://FAkohoue.github.io/LDxBlocks/reference/ldx_blues.md):

  Pre-adjusted BLUEs (id, YLD, RES) for
  [`run_haplotype_prediction`](https://FAkohoue.github.io/LDxBlocks/reference/run_haplotype_prediction.md)
  demos. Available as an .rda dataset and
  `inst/extdata/example_blues.csv`.

## References

Kim S-A et al. (2018). A new haplotype block detection method for dense
genome sequencing data based on interval graph modeling and dynamic
programming. *Bioinformatics* **34**(4):588-596.
[doi:10.1093/bioinformatics/btx609](https://doi.org/10.1093/bioinformatics/btx609)

Mangin B et al. (2012). Novel measures of linkage disequilibrium that
correct the bias due to population structure and relatedness. *Heredity*
**108**(3):285-291.
[doi:10.1038/hdy.2011.73](https://doi.org/10.1038/hdy.2011.73)

VanRaden PM (2008). Efficient methods to compute genomic predictions.
*Journal of Dairy Science* **91**(11):4414-4423.
[doi:10.3168/jds.2007-0980](https://doi.org/10.3168/jds.2007-0980)

Difabachew YF et al. (2023). Genomic prediction with haplotype blocks in
wheat. *Frontiers in Plant Science* **14**:1168547.
[doi:10.3389/fpls.2023.1168547](https://doi.org/10.3389/fpls.2023.1168547)

Weber SE, Frisch M, Snowdon RJ, Voss-Fels KP (2023). Haplotype blocks
for genomic prediction: a comparative evaluation in multiple crop
datasets. *Frontiers in Plant Science* **14**:1217589.
[doi:10.3389/fpls.2023.1217589](https://doi.org/10.3389/fpls.2023.1217589)

Pook T et al. (2019). HaploBlocker: Creation of subgroup-specific
haplotype blocks and libraries. *Genetics* **212**(4):1045-1061.
[doi:10.1534/genetics.119.302283](https://doi.org/10.1534/genetics.119.302283)

Tong J et al. (2024). Stacking beneficial haplotypes from the Vavilov
wheat collection to accelerate breeding for multiple disease resistance.
*Theoretical and Applied Genetics* **137**:274.
[doi:10.1007/s00122-024-04784-w](https://doi.org/10.1007/s00122-024-04784-w)

Tong J et al. (2025). Haplotype stacking to improve stability of stripe
rust resistance in wheat. *Theoretical and Applied Genetics*
**138**:267.
[doi:10.1007/s00122-025-05045-0](https://doi.org/10.1007/s00122-025-05045-0)

Blondel VD, Guillaume J-L, Lambiotte R, Lefebvre E (2008). Fast
unfolding of communities in large networks. *Journal of Statistical
Mechanics: Theory and Experiment* **2008**:P10008.
[doi:10.1088/1742-5468/2008/10/P10008](https://doi.org/10.1088/1742-5468/2008/10/P10008)

Traag VA, Waltman L, van Eck NJ (2019). From Louvain to Leiden:
guaranteeing well-connected communities. *Scientific Reports*
**9**:5233.
[doi:10.1038/s41598-019-41695-z](https://doi.org/10.1038/s41598-019-41695-z)

Hill WG, Weir BS (1988). Variances and covariances of squared linkage
disequilibria in finite populations. *Theoretical Population Biology*
**33**(1):54-78.
[doi:10.1016/0040-5809(88)90004-4](https://doi.org/10.1016/0040-5809%2888%2990004-4)

Remington DL et al. (2001). Structure of linkage disequilibrium and
phenotypic associations in the maize genome. *PNAS*
**98**(20):11479-11484.
[doi:10.1073/pnas.201394398](https://doi.org/10.1073/pnas.201394398)

Gao X, Starmer J, Martin ER (2008). A multiple testing correction method
for genetic association studies using correlated single nucleotide
polymorphisms. *Genetic Epidemiology* **32**:361-369.
[doi:10.1002/gepi.20310](https://doi.org/10.1002/gepi.20310)

Gao X, Becker LC, Becker DM, Starmer JD, Province MA (2010). Avoiding
the high Bonferroni penalty in genome-wide association studies. *Genetic
Epidemiology* **34**:100-105.
[doi:10.1002/gepi.20430](https://doi.org/10.1002/gepi.20430)

Gao X (2011). Multiple testing corrections for imputed SNPs. *Genetic
Epidemiology* **35**:154-158.
[doi:10.1002/gepi.20563](https://doi.org/10.1002/gepi.20563)

Borenstein M, Hedges LV, Higgins JPT, Rothstein HR (2009). *Introduction
to Meta-Analysis*. Wiley.

Higgins JPT, Thompson SG (2002). Quantifying heterogeneity in a
meta-analysis. *Statistics in Medicine* **21**(11):1539-1558.
[doi:10.1002/sim.1186](https://doi.org/10.1002/sim.1186)

## See also

Useful links:

- <https://github.com/FAkohoue/LDxBlocks>

- <https://FAkohoue.github.io/LDxBlocks/>

- Report bugs at <https://github.com/FAkohoue/LDxBlocks/issues>

## Author

**Maintainer**: Félicien Akohoue <akohoue.f@gmail.com>
([ORCID](https://orcid.org/0000-0002-2160-0182))
