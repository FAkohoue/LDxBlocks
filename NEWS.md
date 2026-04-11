## LDxBlocks 0.3.0

### Breaking changes

- `read_geno()`: removed `pheno_ids` parameter. Sample subsetting is the
  caller's responsibility — use `be$sample_ids` to match externally.
- Removed `read_pheno()` and `align_geno_pheno()` entirely. Phenotype handling
  is not part of LD block detection and these functions had no algorithmic role.
- `run_Big_LD_all_chr()`: parameter `rV2method` renamed to `kin_method`.
- `digits` default changed from `6` to `-1` (no rounding) throughout.

### New features — haplotype module (comprehensive rewrite)

- **Phasing functions**:
  - `read_phased_vcf()`: read pre-phased VCF with `0|1` GT fields; returns
    `hap1`/`hap2` gamete matrices plus combined dosage.
  - `phase_with_beagle()`: call Beagle 5.x for statistical phasing of WGS
    data; returns path to phased VCF.gz.
  - `phase_with_pedigree()`: Mendelian allele transmission phasing within
    parent-offspring trios; exact when parents are homozygous.
  - `unphase_to_dosage()`: collapse phased gamete matrices back to 0/1/2.
- **Haplotype extraction** (`extract_haplotypes()`):
  - Auto-detects phased vs unphased input from list structure.
  - Unphased mode: diploid string `"012201"` per individual per block.
  - Phased mode: two-gamete string `"011|100"` per individual per block.
  - Blocks are strictly per-chromosome; cross-chromosome blocks are
    architecturally impossible.
- **Feature matrix** (`build_haplotype_feature_matrix()`):
  - New `encoding =` parameter: `"additive_012"` (default) or `"presence_02"`.
  - Phased + `"additive_012"`: true 0/1/2 allele counts per gamete.
  - Unphased + `"additive_012"`: 0 or 2 (1 not identifiable without phase).
  - `"presence_02"`: 0/2 presence/absence for kernel methods or random forests.
  - New `min_freq =` parameter to drop rare haplotype allele columns.
- **QTL region definition** (`define_qtl_regions()`):
  - Maps significant GWAS markers onto LD blocks post-GWAS.
  - Flags pleiotropic blocks (significant hits from multiple traits).
  - Implements the haploblock-based QTL cataloguing approach of
    Tong et al. (2024) *Theor Appl Genet* 137:274.
- **Output writers**:
  - `write_haplotype_numeric()`: CSV with rows = individuals,
    columns = haplotype alleles (0/1/2 or 0/2); compatible with rrBLUP,
    BGLR, ASReml-R.
  - `write_haplotype_hapmap()`: HapMap tab-delimited format with rows =
    haplotype alleles, columns = individuals; nucleotide encoding
    0→HH, 1→HA, 2→AA, NA→NN; compatible with TASSEL and GAPIT.
  - `write_haplotype_diversity()`: CSV of per-block diversity metrics with
    optional genome-wide mean summary row.
- **End-to-end pipeline** (`run_ldx_pipeline()`): single-call wrapper from
  genotype file to all outputs (blocks, diversity, haplotype matrix).

### New features — I/O and backend

- **Multi-format I/O**: `read_geno()` supports numeric dosage CSV, HapMap,
  VCF / bgzipped VCF, SeqArray GDS, PLINK BED/BIM/FAM, and plain R matrices.
- **GDS streaming backend**: peak RAM proportional to one window.
- **PLINK BED backend**: `BEDMatrix`-backed memory-mapped access.
- **`read_chunk()`**, **`close_backend()`**: unified accessor interface.
- **`run_Big_LD_all_chr()`** accepts `LDxBlocks_backend` directly.
- **`min_snps_chr`** parameter: skip scaffolds with fewer SNPs (default 10).
- **`ldx_gwas`** example dataset now includes a `trait` column
  (`"TraitA"` / `"TraitB"`) to demonstrate multi-trait pleiotropic block
  detection via `define_qtl_regions(trait_col = "trait")`. For single-trait
  GWAS the argument is omitted entirely.

### Performance improvements — C++ core

- `boundary_scan_cpp()`: C++ boundary scan replaces R inner loop.
- `maf_filter_cpp()`: combined MAF + monomorphic filter, single O(np) pass.
- `build_adj_matrix_cpp()`: in-place adjacency construction.
- `compute_r2_sparse_cpp()`: sparse r² within bp window.

### Performance improvements — never-full-genome memory model (OptSLDP-inspired)

- **Chunked pre-allocated numeric CSV reader** (`read_geno()` numeric path):
  two-pass strategy — header scan only in Pass 1, then fixed 50,000-row chunks
  fill a single pre-allocated matrix in Pass 2. Peak RAM = one chunk, not 2×
  the file. `gc(FALSE)` called after each chunk. Directly follows the
  `.read_dosage_chunked()` pattern of Akohoue et al. (2026) *OptSLDP*.
- **VCF and HapMap mandatory GDS auto-conversion**: `read_geno()` for VCF and
  HapMap auto-converts to a SeqArray GDS cache (placed next to the source file)
  when SeqArray is available, transparently switching to the streaming GDS
  backend. Cache reused on subsequent calls.
- **Explicit per-chromosome gc()** in `run_Big_LD_all_chr()`: `rm(geno_chr,
  info_chr)` and `gc(FALSE)` called after each chromosome's block detection
  completes, preventing heap fragmentation across 20–30 chromosome passes.
- **Backend streaming path in `extract_haplotypes()`**: when an
  `LDxBlocks_backend` object is passed (instead of a matrix), haplotypes are
  extracted one chromosome at a time using `read_chunk()`. Each chromosome is
  freed with `rm()` and `gc(FALSE)` before the next chromosome is loaded. The
  full genome is never in RAM simultaneously.

### Tests

- `test-basic.R`: smoke tests via example data (230 SNPs, 3 chromosomes).
- `test-cpp.R`: property-based C++ kernel tests.
- `test-io.R`: all I/O format round-trip tests via tempfiles.
- `test-haplotypes.R`: haplotype extraction, diversity, QTL, feature matrix,
  and all three output writers.

---

## LDxBlocks 0.2.0 (2025-04-07)

### New features

- C++/Armadillo core: `compute_r2_cpp()`, `compute_rV2_cpp()`,
  `maf_filter_cpp()`, `build_adj_matrix_cpp()`, `col_r2_cpp()`,
  `compute_r2_sparse_cpp()`, `boundary_scan_cpp()`.
- OpenMP parallelism via `n_threads` parameter.
- Dual LD metric: `method = "r2"` or `"rV2"`.
- `prepare_geno()`, `compute_ld()`.

### Bug fixes

- `CLQD()`: fixed `re_idx` tracking after dense-core pre-pass.
- `Big_LD()`: `vapply()` used throughout for type-stability.

---

## LDxBlocks 0.1.0 (2025-04-07)

Initial release.

- `Big_LD()`, `CLQD()`, `run_Big_LD_all_chr()`, `tune_LD_params()`.
- `extract_haplotypes()`, `compute_haplotype_diversity()`,
  `build_haplotype_feature_matrix()`.
- `summarise_blocks()`, `plot_ld_blocks()`.
