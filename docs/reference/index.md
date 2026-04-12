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
GDS, PLINK BED, R matrix) through a unified backend interface.

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

## LD computation

LD matrix computation: standard r2 and kinship-adjusted rV2.

- [`compute_r2()`](https://FAkohoue.github.io/LDxBlocks/reference/compute_r2.md)
  : Compute Standard r^2 LD Matrix
- [`compute_rV2()`](https://FAkohoue.github.io/LDxBlocks/reference/compute_rV2.md)
  : Compute Kinship-Adjusted rV^2 LD Matrix
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

## Haplotype output writers

Write haplotype matrices and diversity tables to disk. Numeric format:
0/1/2/NA dosage, rows=haplotypes, cols=individuals. Character format:
nucleotide sequences, rows=haplotypes, cols=individuals.

- [`write_haplotype_numeric()`](https://FAkohoue.github.io/LDxBlocks/reference/write_haplotype_numeric.md)
  : Write Haplotype Feature Matrix as Numeric Dosage Table
- [`write_haplotype_character()`](https://FAkohoue.github.io/LDxBlocks/reference/write_haplotype_character.md)
  : Write Haplotype Character (Nucleotide) Matrix
- [`write_haplotype_diversity()`](https://FAkohoue.github.io/LDxBlocks/reference/write_haplotype_diversity.md)
  : Write Haplotype Diversity Table

## Utilities

Summary statistics and visualisation for block tables.

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
