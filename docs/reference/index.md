# Package index

## Package

Package-level documentation and overview.

- [`LDxBlocks-package`](https://FAkohoue.github.io/LDxBlocks/reference/LDxBlocks-package.md)
  [`LDxBlocks`](https://FAkohoue.github.io/LDxBlocks/reference/LDxBlocks-package.md)
  : LDxBlocks: Genome-Wide LD Block Detection and Haplotype Analysis

## Main pipeline

High-level functions for genome-wide LD block detection.
run_Big_LD_all_chr() is the recommended entry point. run_ldx_pipeline()
runs the complete end-to-end workflow from a single file path.

- [`run_Big_LD_all_chr()`](https://FAkohoue.github.io/LDxBlocks/reference/run_Big_LD_all_chr.md)
  : Genome-Wide LD Block Detection by Chromosome
- [`tune_LD_params()`](https://FAkohoue.github.io/LDxBlocks/reference/tune_LD_params.md)
  : Auto-Tune LD Block Detection Parameters
- [`run_ldx_pipeline()`](https://FAkohoue.github.io/LDxBlocks/reference/run_ldx_pipeline.md)
  : End-to-End Haplotype Block Pipeline

## Genotype I/O

Read genotype data from six formats (numeric dosage CSV, HapMap, VCF,
GDS, PLINK BED, R matrix) through a unified backend interface. For
WGS-scale panels where peak RAM is a constraint, read_geno_bigmemory()
creates a file-backed memory-mapped store (requires bigmemory) with
OS-page-fault column access.

- [`read_geno()`](https://FAkohoue.github.io/LDxBlocks/reference/read_geno.md)
  : Read Genotype Data into an LDxBlocks Backend
- [`read_chunk()`](https://FAkohoue.github.io/LDxBlocks/reference/read_chunk.md)
  : Extract a Genotype Slice from an LDxBlocks Backend
- [`close_backend()`](https://FAkohoue.github.io/LDxBlocks/reference/close_backend.md)
  : Close an LDxBlocks Backend and Release File Handles
- [`print(`*`<LDxBlocks_backend>`*`)`](https://FAkohoue.github.io/LDxBlocks/reference/print.LDxBlocks_backend.md)
  : Print Method for LDxBlocks Backend
- [`summary(`*`<LDxBlocks_backend>`*`)`](https://FAkohoue.github.io/LDxBlocks/reference/summary.LDxBlocks_backend.md)
  : Summary Method for LDxBlocks Backend
- [`read_geno_bigmemory()`](https://FAkohoue.github.io/LDxBlocks/reference/read_geno_bigmemory.md)
  : Open a bigmemory-backed Genotype Store

## LD computation

LD matrix computation, decay analysis, and threshold estimation.

- [`compute_r2()`](https://FAkohoue.github.io/LDxBlocks/reference/compute_r2.md)
  : Compute Standard r^2 LD Matrix
- [`compute_rV2()`](https://FAkohoue.github.io/LDxBlocks/reference/compute_rV2.md)
  : Compute Kinship-Adjusted rV^2 LD Matrix
- [`compute_ld_decay()`](https://FAkohoue.github.io/LDxBlocks/reference/compute_ld_decay.md)
  : Compute LD Decay and Chromosome-Specific Decay Distances
- [`plot_ld_decay()`](https://FAkohoue.github.io/LDxBlocks/reference/plot_ld_decay.md)
  : Plot LD Decay Curve
- [`prepare_geno()`](https://FAkohoue.github.io/LDxBlocks/reference/prepare_geno.md)
  : Prepare Genotype Matrix for LD Computation
- [`get_V_inv_sqrt()`](https://FAkohoue.github.io/LDxBlocks/reference/get_V_inv_sqrt.md)
  : Compute the Inverse Square Root (Whitening Factor) of a Kinship
  Matrix

## Phasing

Functions for obtaining gametic phase from WGS data (Beagle), pedigree
records, or pre-phased VCF files.

- [`read_phased_vcf()`](https://FAkohoue.github.io/LDxBlocks/reference/read_phased_vcf.md)
  : Read Pre-Phased VCF
- [`phase_with_beagle()`](https://FAkohoue.github.io/LDxBlocks/reference/phase_with_beagle.md)
  : Statistical Phasing via Beagle 5.x
- [`phase_with_pedigree()`](https://FAkohoue.github.io/LDxBlocks/reference/phase_with_pedigree.md)
  : Pedigree-Based Allele Transmission Phasing

## Haplotype analysis

Extract haplotypes within LD blocks (phased or unphased), decode
haplotype strings to nucleotide sequences, compute diversity metrics,
map GWAS hits to QTL regions, and build feature matrices for genomic
prediction.

- [`extract_haplotypes()`](https://FAkohoue.github.io/LDxBlocks/reference/extract_haplotypes.md)
  : Extract Haplotype Strings Within LD Blocks
- [`decode_haplotype_strings()`](https://FAkohoue.github.io/LDxBlocks/reference/decode_haplotype_strings.md)
  : Decode Haplotype Strings to Nucleotide Sequences
- [`compute_haplotype_diversity()`](https://FAkohoue.github.io/LDxBlocks/reference/compute_haplotype_diversity.md)
  : Compute Haplotype Diversity Per Block
- [`define_qtl_regions()`](https://FAkohoue.github.io/LDxBlocks/reference/define_qtl_regions.md)
  : Map GWAS Hits to LD Blocks (Post-GWAS QTL Region Definition)
- [`build_haplotype_feature_matrix()`](https://FAkohoue.github.io/LDxBlocks/reference/build_haplotype_feature_matrix.md)
  : Build Haplotype Dosage Matrix for Genomic Prediction
- [`compute_haplotype_grm()`](https://FAkohoue.github.io/LDxBlocks/reference/compute_haplotype_grm.md)
  : Compute Haplotype-Based Genomic Relationship Matrix

## Haplotype output writers

Write haplotype matrices and diversity tables to disk. Numeric format:
phased=0/1/2/NA, unphased=0/1/NA dosage; rows=haplotypes,
cols=individuals. Character format: nucleotide sequences with IUPAC
ambiguity codes for heterozygous positions; rows=haplotypes,
cols=individuals.

- [`write_haplotype_numeric()`](https://FAkohoue.github.io/LDxBlocks/reference/write_haplotype_numeric.md)
  : Write Haplotype Feature Matrix as Numeric Dosage Table
- [`write_haplotype_character()`](https://FAkohoue.github.io/LDxBlocks/reference/write_haplotype_character.md)
  : Write Haplotype Character (Nucleotide) Matrix
- [`write_haplotype_diversity()`](https://FAkohoue.github.io/LDxBlocks/reference/write_haplotype_diversity.md)
  : Write Haplotype Diversity Table

## Genomic prediction pipeline (Tong et al. 2024/2025)

Haplotype-based genomic prediction following Tong et al. (2024/2025).
run_haplotype_prediction() accepts a single trait (named vector or
one-column data frame) or multiple traits (wide data frame or named
list). All traits are fitted via rrBLUP::kin.blup() per trait using a
shared GRM. Block importance is aggregated across traits for robust
haplotype stacking candidate identification.

- [`prepare_gblup_inputs()`](https://FAkohoue.github.io/LDxBlocks/reference/prepare_gblup_inputs.md)
  : Prepare Genomic Prediction Inputs for External GBLUP Software
- [`run_haplotype_prediction()`](https://FAkohoue.github.io/LDxBlocks/reference/run_haplotype_prediction.md)
  : Haplotype Prediction and Block Importance from Pre-Adjusted
  Phenotypes
- [`backsolve_snp_effects()`](https://FAkohoue.github.io/LDxBlocks/reference/backsolve_snp_effects.md)
  : Backsolve SNP Effects from GEBV (Tong et al. 2025)
- [`compute_local_gebv()`](https://FAkohoue.github.io/LDxBlocks/reference/compute_local_gebv.md)
  : Compute Local Haplotype GEBV per Block (Tong et al. 2025)
- [`integrate_gwas_haplotypes()`](https://FAkohoue.github.io/LDxBlocks/reference/integrate_gwas_haplotypes.md)
  : Integrate GWAS QTL Regions with Haplotype Prediction Results
- [`rank_haplotype_blocks()`](https://FAkohoue.github.io/LDxBlocks/reference/rank_haplotype_blocks.md)
  : Rank Haplotype Blocks by Evidence Strength

## Analysis extensions

Cross-validation, population comparison, stability, export, and effect
decomposition

- [`cv_haplotype_prediction()`](https://FAkohoue.github.io/LDxBlocks/reference/cv_haplotype_prediction.md)
  : K-Fold Cross-Validation for Haplotype-Based Genomic Prediction
- [`compare_haplotype_populations()`](https://FAkohoue.github.io/LDxBlocks/reference/compare_haplotype_populations.md)
  : Compare Haplotype Allele Frequencies Between Two Population Groups
- [`plot_haplotype_network()`](https://FAkohoue.github.io/LDxBlocks/reference/plot_haplotype_network.md)
  : Plot a Minimum-Spanning Haplotype Network for One LD Block
- [`run_haplotype_stability()`](https://FAkohoue.github.io/LDxBlocks/reference/run_haplotype_stability.md)
  : Finlay-Wilkinson Stability Analysis of Haplotype Effects Across
  Environments
- [`export_candidate_regions()`](https://FAkohoue.github.io/LDxBlocks/reference/export_candidate_regions.md)
  : Export Candidate Gene Regions to BED, CSV, or biomaRt Format
- [`decompose_block_effects()`](https://FAkohoue.github.io/LDxBlocks/reference/decompose_block_effects.md)
  : Decompose Per-SNP Effects into Per-Haplotype-Allele Effect Table
- [`scan_diversity_windows()`](https://FAkohoue.github.io/LDxBlocks/reference/scan_diversity_windows.md)
  : Sliding-Window Genome-Wide Diversity Scan

## Haplotype inference and harmonisation

True diplotype inference, rare-allele collapsing, and cross-panel label
harmonisation

- [`infer_block_haplotypes()`](https://FAkohoue.github.io/LDxBlocks/reference/infer_block_haplotypes.md)
  : Infer Structured Block-Level Diplotypes Per Individual
- [`collapse_haplotypes()`](https://FAkohoue.github.io/LDxBlocks/reference/collapse_haplotypes.md)
  : Collapse Rare Haplotype Alleles Into Biologically Meaningful Groups
- [`harmonize_haplotypes()`](https://FAkohoue.github.io/LDxBlocks/reference/harmonize_haplotypes.md)
  : Harmonize Haplotype Allele Labels Across Panels or Analysis Runs

## Haplotype association testing

Block-level association tests and diplotype effect estimation (Q+K mixed
model)

- [`test_block_haplotypes()`](https://FAkohoue.github.io/LDxBlocks/reference/test_block_haplotypes.md)
  : Block-Level Haplotype Association Testing (Q+K Mixed Linear Model)
- [`estimate_diplotype_effects()`](https://FAkohoue.github.io/LDxBlocks/reference/estimate_diplotype_effects.md)
  : Estimate Diplotype Effects and Dominance Deviations Per LD Block

## Breeding decision tools

Haplotype portfolio scoring and parent allele inventory for stacking
decisions

- [`score_favorable_haplotypes()`](https://FAkohoue.github.io/LDxBlocks/reference/score_favorable_haplotypes.md)
  : Score Individual Haplotype Portfolios Against Known Allele Effects
- [`summarize_parent_haplotypes()`](https://FAkohoue.github.io/LDxBlocks/reference/summarize_parent_haplotypes.md)
  : Summarise Haplotype Allele Inventory Per Candidate Parent

## Utilities

Summary statistics and visualisation.

- [`summarise_blocks()`](https://FAkohoue.github.io/LDxBlocks/reference/summarise_blocks.md)
  : Summarise LD Block Characteristics
- [`plot_ld_blocks()`](https://FAkohoue.github.io/LDxBlocks/reference/plot_ld_blocks.md)
  : Plot LD Block Structure Across Chromosomes

## Example data

Simulated datasets (120 individuals, 230 SNPs, 3 chromosomes, 9 LD
blocks) for examples and tests.

- [`ldx_geno`](https://FAkohoue.github.io/LDxBlocks/reference/ldx_geno.md)
  : Example Genotype Matrix
- [`ldx_snp_info`](https://FAkohoue.github.io/LDxBlocks/reference/ldx_snp_info.md)
  : Example SNP Information Table
- [`ldx_blocks`](https://FAkohoue.github.io/LDxBlocks/reference/ldx_blocks.md)
  : Example LD Block Table
- [`ldx_gwas`](https://FAkohoue.github.io/LDxBlocks/reference/ldx_gwas.md)
  : Example GWAS Marker Table
- [`ldx_blues`](https://FAkohoue.github.io/LDxBlocks/reference/ldx_blues.md)
  : Pre-Adjusted Phenotype Means (BLUEs) for Genomic Prediction Examples
- [`ldx_blues_list`](https://FAkohoue.github.io/LDxBlocks/reference/ldx_blues_list.md)
  : Per-Environment BLUEs for Stability Analysis Examples
