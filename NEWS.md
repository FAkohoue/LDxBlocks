# LDxBlocks News

## LDxBlocks 0.3.0 (2025-04-07)

### Breaking changes

- `read_geno()`: removed `pheno_ids` parameter. Sample subsetting is the
  caller's responsibility — use `be$sample_ids` to match externally.
- Removed `read_pheno()` and `align_geno_pheno()` entirely. Phenotype handling
  is not part of LD block detection and these functions had no algorithmic
  role.
- `run_Big_LD_all_chr()`: parameter `rV2method` renamed to `kin_method` for
  clarity.
- `digits` default changed from `6` to `-1` (no rounding) throughout. Rounding
  was a source of subtle threshold-comparison instability. Set `digits = 6`
  explicitly to restore previous behaviour.

### New features

- **Multi-format I/O**: `read_geno()` now supports numeric dosage CSV, HapMap,
  VCF / bgzipped VCF, SeqArray GDS, PLINK BED/BIM/FAM, and plain R matrices
  through a unified `LDxBlocks_backend` S3 interface.
- **GDS streaming backend**: SeqArray GDS files are accessed via
  `seqSetFilter()` + `seqGetData("$dosage")` per window, keeping peak RAM
  proportional to one chromosome window rather than the full genome.
- **PLINK BED backend**: `BEDMatrix`-backed streaming access for binary PLINK
  files.
- **`read_chunk()`**: unified genotype slice accessor — works identically
  for all six backend types.
- **`close_backend()`**: releases GDS file handles; no-op for in-memory
  backends.
- **`run_Big_LD_all_chr()` now accepts `LDxBlocks_backend`** directly, removing
  the need to load the full genome matrix before calling it.
- **`min_snps_chr`** parameter in `run_Big_LD_all_chr()`: skip chromosomes or
  scaffolds with fewer than this many SNPs (default 10).
- **Example datasets**: `ldx_geno`, `ldx_snp_info`, `ldx_blocks`, `ldx_gwas`
  — simulated with clear block structure for examples and tests.
- **Flat-file examples** in `inst/extdata/`: numeric CSV, HapMap, and VCF
  copies of the example data for format-reader tests.

### Performance improvements

- **`boundary_scan_cpp()`**: C++ replacement for the R inner loop in
  `cutsequence.modi()`. For a 50,000-SNP chromosome with `leng = 200`,
  eliminates ~150,000 small R-level matrix operations per chromosome.
- **`maf_filter_cpp()`**: combined MAF + monomorphic filter in one O(np) C++
  pass; ~10× faster than the `apply()` equivalent.
- **`build_adj_matrix_cpp()`**: in-place C++ adjacency construction, avoiding
  the intermediate logical matrix from `ifelse()`.
- **`compute_r2_sparse_cpp()`**: sparse r² within a bp distance window,
  returning only above-threshold pairs as a triplet (row, col, r²). Avoids
  the O(p²) dense matrix for large sub-segments where distant pairs are always
  below threshold.

### Tests

- Three test files: `test-basic.R` (smoke tests via example data),
  `test-cpp.R` (property-based C++ kernel tests), `test-io.R` (all I/O format
  tests via tempfiles).

---

## LDxBlocks 0.2.0 (2025-04-07)

### New features

- **C++/Armadillo core**: `compute_r2_cpp()`, `compute_rV2_cpp()`,
  `maf_filter_cpp()`, `build_adj_matrix_cpp()`, `col_r2_cpp()`,
  `compute_r2_sparse_cpp()`, `boundary_scan_cpp()` via RcppArmadillo.
- **OpenMP parallelism**: `n_threads` parameter for `compute_r2_cpp()`.
- **Dual LD metric**: `method = "r2"` (default) or `"rV2"` in `Big_LD()` and
  `run_Big_LD_all_chr()`.
- **`prepare_geno()`**: unified preparation step — centres (r²) or centres +
  whitens (rV²).
- **`compute_ld()`**: dispatcher routing to `compute_r2_cpp()` or
  `compute_rV2_cpp()`.

### Bug fixes

- `CLQD()`: fixed latent bug in `re_idx` tracking after dense-core pre-pass;
  original indices were not correctly mapped back after removal.
- `Big_LD()`: `vapply()` used throughout instead of `sapply()` for
  type-stability.

---

## LDxBlocks 0.1.0 (2025-04-07)

Initial release.

- `Big_LD()`: core per-chromosome LD block segmentation (Kim et al. 2018).
- `CLQD()`: clique detection on rV² adjacency graph.
- `run_Big_LD_all_chr()`: chromosome-wise wrapper.
- `tune_LD_params()`: grid-search parameter auto-tuner.
- `extract_haplotypes()`, `compute_haplotype_diversity()`,
  `build_haplotype_feature_matrix()`: haplotype analysis module.
- `summarise_blocks()`, `plot_ld_blocks()`: utilities.
