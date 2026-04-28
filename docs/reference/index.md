# Package index

## Package overview

Main documentation entry point for LDxBlocks.

- [`LDxBlocks-package`](https://FAkohoue.github.io/LDxBlocks/reference/LDxBlocks-package.md)
  [`LDxBlocks`](https://FAkohoue.github.io/LDxBlocks/reference/LDxBlocks-package.md)
  : LDxBlocks: Genome-Wide LD Block Detection and Haplotype Analysis

## Main workflows

High-level functions covering the full LDxBlocks pipeline, from genotype
input to LD block detection, haplotype construction, and genomic
prediction.

- [`run_ldx_pipeline()`](https://FAkohoue.github.io/LDxBlocks/reference/run_ldx_pipeline.md)
  : End-to-End Haplotype Block Pipeline
- [`run_Big_LD_all_chr()`](https://FAkohoue.github.io/LDxBlocks/reference/run_Big_LD_all_chr.md)
  : Genome-Wide LD Block Detection by Chromosome
- [`tune_LD_params()`](https://FAkohoue.github.io/LDxBlocks/reference/tune_LD_params.md)
  : Auto-Tune LD Block Detection Parameters

## Genotype input and backend

Import genotype data from multiple formats and manage the LDxBlocks
memory-efficient backend for large genomic datasets.

- [`read_geno()`](https://FAkohoue.github.io/LDxBlocks/reference/read_geno.md)
  : Read Genotype Data into an LDxBlocks Backend
- [`read_geno_bigmemory()`](https://FAkohoue.github.io/LDxBlocks/reference/read_geno_bigmemory.md)
  : Open a bigmemory-backed Genotype Store
- [`read_chunk()`](https://FAkohoue.github.io/LDxBlocks/reference/read_chunk.md)
  : Extract a Genotype Slice from an LDxBlocks Backend
- [`prepare_geno()`](https://FAkohoue.github.io/LDxBlocks/reference/prepare_geno.md)
  : Prepare Genotype Matrix for LD Computation
- [`close_backend()`](https://FAkohoue.github.io/LDxBlocks/reference/close_backend.md)
  : Close an LDxBlocks Backend and Release File Handles
- [`print(`*`<LDxBlocks_backend>`*`)`](https://FAkohoue.github.io/LDxBlocks/reference/print.LDxBlocks_backend.md)
  : Print Method for LDxBlocks Backend
- [`summary(`*`<LDxBlocks_backend>`*`)`](https://FAkohoue.github.io/LDxBlocks/reference/summary.LDxBlocks_backend.md)
  : Summary Method for LDxBlocks Backend

## LD computation and block detection

Core algorithms for LD estimation, LD decay analysis, and genome-wide LD
block segmentation using clique-based and community-detection
approaches.

- [`compute_r2()`](https://FAkohoue.github.io/LDxBlocks/reference/compute_r2.md)
  : Compute Standard r^2 LD Matrix
- [`compute_rV2()`](https://FAkohoue.github.io/LDxBlocks/reference/compute_rV2.md)
  : Compute Kinship-Adjusted rV^2 LD Matrix
- [`get_V_inv_sqrt()`](https://FAkohoue.github.io/LDxBlocks/reference/get_V_inv_sqrt.md)
  : Compute the Inverse Square Root (Whitening Factor) of a Kinship
  Matrix
- [`compute_ld_decay()`](https://FAkohoue.github.io/LDxBlocks/reference/compute_ld_decay.md)
  : Compute LD Decay and Chromosome-Specific Decay Distances
- [`plot_ld_decay()`](https://FAkohoue.github.io/LDxBlocks/reference/plot_ld_decay.md)
  : Plot LD Decay Curve
- [`summarise_blocks()`](https://FAkohoue.github.io/LDxBlocks/reference/summarise_blocks.md)
  : Summarise LD Block Characteristics
- [`plot_ld_blocks()`](https://FAkohoue.github.io/LDxBlocks/reference/plot_ld_blocks.md)
  : Plot LD Block Structure Across Chromosomes

## Phasing tools

Functions for obtaining gametic phase from WGS data (Beagle) or
pre-phased VCF files.

- [`read_phased_vcf()`](https://FAkohoue.github.io/LDxBlocks/reference/read_phased_vcf.md)
  : Read Pre-Phased VCF
- [`phase_with_beagle()`](https://FAkohoue.github.io/LDxBlocks/reference/phase_with_beagle.md)
  : Statistical Phasing via Beagle 5.x

## Haplotype construction and diversity

Extract haplotypes within LD blocks, decode nucleotide sequences,
compute haplotype diversity statistics, and construct genomic prediction
features.

- [`extract_haplotypes()`](https://FAkohoue.github.io/LDxBlocks/reference/extract_haplotypes.md)
  : Extract Haplotype Dosage Strings from LD Blocks
- [`decode_haplotype_strings()`](https://FAkohoue.github.io/LDxBlocks/reference/decode_haplotype_strings.md)
  : Decode Haplotype Strings to Nucleotide Sequences
- [`compute_haplotype_diversity()`](https://FAkohoue.github.io/LDxBlocks/reference/compute_haplotype_diversity.md)
  : Compute Haplotype Diversity Per Block
- [`build_haplotype_feature_matrix()`](https://FAkohoue.github.io/LDxBlocks/reference/build_haplotype_feature_matrix.md)
  : Build Haplotype Dosage Matrix for Genomic Prediction
- [`compute_haplotype_grm()`](https://FAkohoue.github.io/LDxBlocks/reference/compute_haplotype_grm.md)
  : Compute Haplotype-Based Genomic Relationship Matrix
- [`define_qtl_regions()`](https://FAkohoue.github.io/LDxBlocks/reference/define_qtl_regions.md)
  : Map GWAS Hits to LD Blocks (Post-GWAS QTL Region Definition)

## Haplotype harmonisation

Harmonise haplotypes across populations, infer diplotypes, and collapse
rare alleles into biologically meaningful groups.

- [`infer_block_haplotypes()`](https://FAkohoue.github.io/LDxBlocks/reference/infer_block_haplotypes.md)
  : Infer Structured Block-Level Diplotypes Per Individual
- [`collapse_haplotypes()`](https://FAkohoue.github.io/LDxBlocks/reference/collapse_haplotypes.md)
  : Collapse Rare Haplotype Alleles Into Biologically Meaningful Groups
- [`harmonize_haplotypes()`](https://FAkohoue.github.io/LDxBlocks/reference/harmonize_haplotypes.md)
  : Harmonize Haplotype Allele Labels Across Panels or Analysis Runs

## Haplotype export

Export haplotype dosage matrices and diversity summaries for downstream
statistical analyses.

- [`write_haplotype_numeric()`](https://FAkohoue.github.io/LDxBlocks/reference/write_haplotype_numeric.md)
  : Write Haplotype Feature Matrix as Numeric Dosage Table
- [`write_haplotype_character()`](https://FAkohoue.github.io/LDxBlocks/reference/write_haplotype_character.md)
  : Write Haplotype Character (Nucleotide) Matrix
- [`write_haplotype_diversity()`](https://FAkohoue.github.io/LDxBlocks/reference/write_haplotype_diversity.md)
  : Write Haplotype Diversity Table

## Association and QTL interpretation

Block-level association tests (Q+K mixed model with simpleM
multiple-testing correction), diplotype effect decomposition,
cross-population effect concordance from haplotype association results
(compare_block_effects) or external GWAS tools (compare_gwas_effects),
GWAS integration, and epistasis detection.

- [`test_block_haplotypes()`](https://FAkohoue.github.io/LDxBlocks/reference/test_block_haplotypes.md)
  : Block-Level Haplotype Association Testing (Q+K Mixed Linear Model
  with simpleM Multiple-Testing Correction)
- [`estimate_diplotype_effects()`](https://FAkohoue.github.io/LDxBlocks/reference/estimate_diplotype_effects.md)
  : Estimate Diplotype Effects and Dominance Deviations Per LD Block
- [`compare_block_effects()`](https://FAkohoue.github.io/LDxBlocks/reference/compare_block_effects.md)
  : Cross-Population Haplotype Effect Concordance
- [`compare_gwas_effects()`](https://FAkohoue.github.io/LDxBlocks/reference/compare_gwas_effects.md)
  : Cross-Population GWAS Effect Concordance from External Results
- [`integrate_gwas_haplotypes()`](https://FAkohoue.github.io/LDxBlocks/reference/integrate_gwas_haplotypes.md)
  : Integrate GWAS QTL Regions with Haplotype Prediction Results

## Genomic prediction

Haplotype-based genomic prediction following Tong et al. (2024-2025),
including GBLUP integration and local GEBV estimation.

- [`prepare_gblup_inputs()`](https://FAkohoue.github.io/LDxBlocks/reference/prepare_gblup_inputs.md)
  : Prepare Genomic Prediction Inputs for External GBLUP Software
- [`run_haplotype_prediction()`](https://FAkohoue.github.io/LDxBlocks/reference/run_haplotype_prediction.md)
  : Haplotype Prediction and Block Importance from Pre-Adjusted
  Phenotypes
- [`compute_local_gebv()`](https://FAkohoue.github.io/LDxBlocks/reference/compute_local_gebv.md)
  : Compute Local Haplotype GEBV per Block (Tong et al. 2025)
- [`backsolve_snp_effects()`](https://FAkohoue.github.io/LDxBlocks/reference/backsolve_snp_effects.md)
  : Backsolve SNP Effects from GEBV (Tong et al. 2025)

## Epistasis detection

Within-block pairwise SNP interaction scan (scan_block_epistasis),
trans-haplotype between-block epistasis scan
(scan_block_by_block_epistasis), and single-block fine-mapping with
pairwise or LASSO dispatch (fine_map_epistasis_block). All functions
operate on GRM-corrected REML residuals from test_block_haplotypes().

- [`scan_block_epistasis()`](https://FAkohoue.github.io/LDxBlocks/reference/scan_block_epistasis.md)
  : Within-Block Pairwise SNP Epistasis Scan
- [`scan_block_by_block_epistasis()`](https://FAkohoue.github.io/LDxBlocks/reference/scan_block_by_block_epistasis.md)
  : Between-Block Haplotype Allele Epistasis Scan
- [`fine_map_epistasis_block()`](https://FAkohoue.github.io/LDxBlocks/reference/fine_map_epistasis_block.md)
  : Fine-Map Epistatic SNP Pairs Within a Single Block

## Breeding decision support

Identify favourable haplotypes and rank LD blocks for haplotype stacking
strategies.

- [`rank_haplotype_blocks()`](https://FAkohoue.github.io/LDxBlocks/reference/rank_haplotype_blocks.md)
  : Rank Haplotype Blocks by Evidence Strength
- [`score_favorable_haplotypes()`](https://FAkohoue.github.io/LDxBlocks/reference/score_favorable_haplotypes.md)
  : Score Individual Haplotype Portfolios Against Known Allele Effects
- [`summarize_parent_haplotypes()`](https://FAkohoue.github.io/LDxBlocks/reference/summarize_parent_haplotypes.md)
  : Summarise Haplotype Allele Inventory Per Candidate Parent

## Population and stability analysis

Compare haplotype distributions across populations and evaluate
stability across environments.

- [`compare_haplotype_populations()`](https://FAkohoue.github.io/LDxBlocks/reference/compare_haplotype_populations.md)
  : Compare Haplotype Allele Frequencies Between Two Population Groups
- [`run_haplotype_stability()`](https://FAkohoue.github.io/LDxBlocks/reference/run_haplotype_stability.md)
  : Finlay-Wilkinson Stability Analysis of Haplotype Effects Across
  Environments
- [`scan_diversity_windows()`](https://FAkohoue.github.io/LDxBlocks/reference/scan_diversity_windows.md)
  : Sliding-Window Genome-Wide Diversity Scan
- [`decompose_block_effects()`](https://FAkohoue.github.io/LDxBlocks/reference/decompose_block_effects.md)
  : Decompose Per-SNP Effects into Per-Haplotype-Allele Effect Table
- [`export_candidate_regions()`](https://FAkohoue.github.io/LDxBlocks/reference/export_candidate_regions.md)
  : Export Candidate Gene Regions to BED, CSV, or biomaRt Format

## Cross-validation

Evaluate predictive ability of haplotype-based genomic models.

- [`cv_haplotype_prediction()`](https://FAkohoue.github.io/LDxBlocks/reference/cv_haplotype_prediction.md)
  : K-Fold Cross-Validation for Haplotype-Based Genomic Prediction

## Visualisation

Plot LD structure and haplotype relationships.

- [`plot_haplotype_network()`](https://FAkohoue.github.io/LDxBlocks/reference/plot_haplotype_network.md)
  : Plot a Minimum-Spanning Haplotype Network for One LD Block

## Example datasets

Simulated datasets for tutorials and reproducible examples.

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
