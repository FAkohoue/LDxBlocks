# Changelog

## LDxBlocks 0.3.0

### Breaking changes

- [`read_geno()`](https://FAkohoue.github.io/LDxBlocks/reference/read_geno.md):
  removed `pheno_ids` parameter. Sample subsetting is the caller’s
  responsibility – use `be$sample_ids` to match externally.
- Removed `read_pheno()` and `align_geno_pheno()` entirely. Phenotype
  handling is not part of LD block detection and these functions had no
  algorithmic role.
- [`run_Big_LD_all_chr()`](https://FAkohoue.github.io/LDxBlocks/reference/run_Big_LD_all_chr.md):
  parameter `rV2method` renamed to `kin_method`.
- `digits` default changed from `6` to `-1` (no rounding) throughout.
- SeqArray backend replaced by SNPRelate throughout. GDS files are now
  created via
  [`SNPRelate::snpgdsVCF2GDS()`](https://rdrr.io/pkg/SNPRelate/man/snpgdsVCF2GDS.html)
  and read via
  [`SNPRelate::snpgdsGetGeno()`](https://rdrr.io/pkg/SNPRelate/man/snpgdsGetGeno.html).
  Existing `.gds` files created by SeqArray are not compatible; delete
  and re-convert with `read_geno(..., verbose = TRUE)`.

### New features – haplotype module (comprehensive rewrite)

- **Phasing functions**:
  - [`read_phased_vcf()`](https://FAkohoue.github.io/LDxBlocks/reference/read_phased_vcf.md):
    read pre-phased VCF with `0|1` GT fields; returns `hap1`/`hap2`
    gamete matrices plus combined dosage.
  - [`phase_with_beagle()`](https://FAkohoue.github.io/LDxBlocks/reference/phase_with_beagle.md):
    call Beagle 5.x for statistical phasing of WGS data; returns path to
    phased VCF.gz.
  - [`phase_with_pedigree()`](https://FAkohoue.github.io/LDxBlocks/reference/phase_with_pedigree.md):
    Mendelian allele transmission phasing within parent-offspring trios;
    exact when parents are homozygous.
  - [`unphase_to_dosage()`](https://FAkohoue.github.io/LDxBlocks/reference/unphase_to_dosage.md):
    collapse phased gamete matrices back to 0/1/2 (internal helper; not
    exported).
- **Haplotype extraction**
  ([`extract_haplotypes()`](https://FAkohoue.github.io/LDxBlocks/reference/extract_haplotypes.md)):
  - Auto-detects phased vs unphased input from list structure.
  - Unphased mode: diploid string `"012201"` per individual per block.
  - Phased mode: two-gamete string `"011|100"` per individual per block.
  - Blocks are strictly per-chromosome; cross-chromosome blocks are
    architecturally impossible.
- **Haplotype string decoder**
  ([`decode_haplotype_strings()`](https://FAkohoue.github.io/LDxBlocks/reference/decode_haplotype_strings.md)):
  - Converts raw dosage strings (e.g. `"02110"`) to nucleotide sequences
    (e.g. `"AGTTA"`) using REF/ALT from `snp_info`.
  - Returns a data frame: block_id, CHR, start_bp, end_bp, hap_rank,
    hap_id, dosage_string, nucleotide_sequence, frequency, n_carriers,
    snp_positions, snp_alleles.
- **Feature matrix**
  ([`build_haplotype_feature_matrix()`](https://FAkohoue.github.io/LDxBlocks/reference/build_haplotype_feature_matrix.md)):
  - New `encoding =` parameter: `"additive_012"` (default) or
    `"presence_02"`.
  - Phased + `"additive_012"`: true 0/1/2 allele counts per gamete.
  - Unphased + `"additive_012"`: 0 or 2 (1 not identifiable without
    phase).
  - `"presence_02"`: 0/2 presence/absence for kernel methods or random
    forests.
  - New `min_freq =` parameter to drop rare haplotype allele columns.
  - `top_n` default changed from `5L` to `NULL` (retain all alleles
    above `min_freq`).
- **QTL region definition**
  ([`define_qtl_regions()`](https://FAkohoue.github.io/LDxBlocks/reference/define_qtl_regions.md)):
  - Maps significant GWAS markers onto LD blocks post-GWAS.
  - Flags pleiotropic blocks (significant hits from multiple traits).
  - Implements the haploblock-based QTL cataloguing approach of Tong et
    al. (2024) *Theoretical and Applied Genetics* 137:274.
- **Output writers** (all produce rows = haplotype alleles, cols =
  individuals, with metadata columns hap_id, CHR, start_bp, end_bp,
  n_snps before individual columns):
  - [`write_haplotype_numeric()`](https://FAkohoue.github.io/LDxBlocks/reference/write_haplotype_numeric.md):
    tab-delimited dosage matrix (0/1/2/NA). Metadata column `alleles`
    shows the nucleotide sequence of this allele. Compatible with
    rrBLUP, BGLR, ASReml-R.
  - [`write_haplotype_character()`](https://FAkohoue.github.io/LDxBlocks/reference/write_haplotype_character.md):
    tab-delimited nucleotide matrix. Each cell shows the nucleotide
    sequence if the individual carries that allele, `"-"` if absent,
    `"."` if missing. `Alleles` column shows all REF/ALT at every block
    SNP (semicolon-separated).
  - [`write_haplotype_diversity()`](https://FAkohoue.github.io/LDxBlocks/reference/write_haplotype_diversity.md):
    CSV of per-block diversity metrics with optional genome-wide mean
    summary row.
  - `write_haplotype_hapmap()` is **removed**. The diploid AA/AT/TT
    HapMap encoding was not meaningful for multi-SNP haplotype alleles.
- **Haplotype string decoder**
  ([`decode_haplotype_strings()`](https://FAkohoue.github.io/LDxBlocks/reference/decode_haplotype_strings.md)):
  exported function mapping dosage strings to nucleotide sequences with
  full SNP-level metadata.
- **End-to-end pipeline**
  ([`run_ldx_pipeline()`](https://FAkohoue.github.io/LDxBlocks/reference/run_ldx_pipeline.md)):
  single-call wrapper from genotype file to all outputs.
  `hap_format = "hapmap"` renamed to `hap_format = "character"`.

### New features – I/O and backend

- **Multi-format I/O**:
  [`read_geno()`](https://FAkohoue.github.io/LDxBlocks/reference/read_geno.md)
  supports numeric dosage CSV, HapMap, VCF / bgzipped VCF, SNPRelate
  GDS, PLINK BED/BIM/FAM, and plain R matrices.
- **SNPRelate GDS streaming backend**: peak RAM proportional to one
  window. Replaces the previous SeqArray backend entirely.
- **PLINK BED backend**: `BEDMatrix`-backed memory-mapped access.
- **[`read_chunk()`](https://FAkohoue.github.io/LDxBlocks/reference/read_chunk.md)**,
  **[`close_backend()`](https://FAkohoue.github.io/LDxBlocks/reference/close_backend.md)**:
  unified accessor interface.
- **[`run_Big_LD_all_chr()`](https://FAkohoue.github.io/LDxBlocks/reference/run_Big_LD_all_chr.md)**
  accepts `LDxBlocks_backend` directly.
- **`min_snps_chr`** parameter: skip scaffolds with fewer SNPs (default
  10).
- **`clean_malformed`** parameter in
  [`read_geno()`](https://FAkohoue.github.io/LDxBlocks/reference/read_geno.md),
  [`run_Big_LD_all_chr()`](https://FAkohoue.github.io/LDxBlocks/reference/run_Big_LD_all_chr.md),
  and
  [`run_ldx_pipeline()`](https://FAkohoue.github.io/LDxBlocks/reference/run_ldx_pipeline.md):
  streams input files to remove malformed lines before GDS conversion
  (for NGSEP and some GATK-origin VCFs).
- **`chr`** parameter in
  [`run_Big_LD_all_chr()`](https://FAkohoue.github.io/LDxBlocks/reference/run_Big_LD_all_chr.md)
  and
  [`run_ldx_pipeline()`](https://FAkohoue.github.io/LDxBlocks/reference/run_ldx_pipeline.md):
  restrict processing to specific chromosomes.
- **`ldx_gwas`** example dataset now includes a `trait` column
  (`"TraitA"` / `"TraitB"`) to demonstrate multi-trait pleiotropic block
  detection via `define_qtl_regions(trait_col = "trait")`.

### Performance improvements – C++ core

- `boundary_scan_cpp()`: C++ boundary scan replaces R inner loop (~20x
  faster).
- `maf_filter_cpp()`: combined MAF + monomorphic filter, single O(np)
  pass (~10x faster for panels \> 100 k SNPs).
- `build_adj_matrix_cpp()`: in-place adjacency construction, eliminates
  intermediate allocation.
- `compute_r2_sparse_cpp()`: sparse r² within bp window, avoids O(p²)
  cost for large sub-segments.

### Performance improvements – never-full-genome memory model

- **Chunked pre-allocated numeric CSV reader**: two-pass strategy –
  header scan only in Pass 1, then fixed 50,000-row chunks fill a single
  pre-allocated matrix in Pass 2. Peak RAM = one chunk, not 2x the file.
- **VCF and HapMap mandatory SNPRelate GDS auto-conversion**:
  transparent streaming GDS backend on first call; cache reused on
  subsequent calls.
- **Explicit per-chromosome gc()** in
  [`run_Big_LD_all_chr()`](https://FAkohoue.github.io/LDxBlocks/reference/run_Big_LD_all_chr.md)
  and
  [`extract_haplotypes()`](https://FAkohoue.github.io/LDxBlocks/reference/extract_haplotypes.md):
  prevents heap fragmentation across 20-30 chromosome passes.

### Tests

- `test-basic.R`: smoke tests via example data (230 SNPs, 3
  chromosomes).
- `test-cpp.R`: property-based C++ kernel tests.
- `test-io.R`: all I/O format round-trip tests via tempfiles.
- `test-haplotypes.R`: haplotype extraction, diversity, QTL, feature
  matrix, and output writers (`write_haplotype_numeric`,
  `write_haplotype_character`, `write_haplotype_diversity`).

------------------------------------------------------------------------

## LDxBlocks 0.2.0 (2025-04-07)

### New features

- C++/Armadillo core: `compute_r2_cpp()`, `compute_rV2_cpp()`,
  `maf_filter_cpp()`, `build_adj_matrix_cpp()`, `col_r2_cpp()`,
  `compute_r2_sparse_cpp()`, `boundary_scan_cpp()`.
- OpenMP parallelism via `n_threads` parameter.
- Dual LD metric: `method = "r2"` or `"rV2"`.
- [`prepare_geno()`](https://FAkohoue.github.io/LDxBlocks/reference/prepare_geno.md),
  [`compute_ld()`](https://FAkohoue.github.io/LDxBlocks/reference/compute_ld.md).

### Bug fixes

- [`CLQD()`](https://FAkohoue.github.io/LDxBlocks/reference/CLQD.md):
  fixed `re_idx` tracking after dense-core pre-pass.
- [`Big_LD()`](https://FAkohoue.github.io/LDxBlocks/reference/Big_LD.md):
  [`vapply()`](https://rdrr.io/r/base/lapply.html) used throughout for
  type-stability.

------------------------------------------------------------------------

## LDxBlocks 0.1.0 (2025-04-07)

Initial release.

- [`Big_LD()`](https://FAkohoue.github.io/LDxBlocks/reference/Big_LD.md),
  [`CLQD()`](https://FAkohoue.github.io/LDxBlocks/reference/CLQD.md),
  [`run_Big_LD_all_chr()`](https://FAkohoue.github.io/LDxBlocks/reference/run_Big_LD_all_chr.md),
  [`tune_LD_params()`](https://FAkohoue.github.io/LDxBlocks/reference/tune_LD_params.md).
- [`extract_haplotypes()`](https://FAkohoue.github.io/LDxBlocks/reference/extract_haplotypes.md),
  [`compute_haplotype_diversity()`](https://FAkohoue.github.io/LDxBlocks/reference/compute_haplotype_diversity.md),
  [`build_haplotype_feature_matrix()`](https://FAkohoue.github.io/LDxBlocks/reference/build_haplotype_feature_matrix.md).
- [`summarise_blocks()`](https://FAkohoue.github.io/LDxBlocks/reference/summarise_blocks.md),
  [`plot_ld_blocks()`](https://FAkohoue.github.io/LDxBlocks/reference/plot_ld_blocks.md).
