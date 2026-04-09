# Changelog

## LDxBlocks 0.3.0

#### Breaking changes

- [`read_geno()`](https://FAkohoue.github.io/LDxBlocks/reference/read_geno.md):
  removed `pheno_ids` parameter. Sample subsetting is the caller’s
  responsibility — use `be$sample_ids` to match externally.
- Removed `read_pheno()` and `align_geno_pheno()` entirely. Phenotype
  handling is not part of LD block detection and these functions had no
  algorithmic role.
- [`run_Big_LD_all_chr()`](https://FAkohoue.github.io/LDxBlocks/reference/run_Big_LD_all_chr.md):
  parameter `rV2method` renamed to `kin_method` for clarity.
- `digits` default changed from `6` to `-1` (no rounding) throughout.
  Rounding was a source of subtle threshold-comparison instability. Set
  `digits = 6` explicitly to restore previous behaviour.

#### New features

- **Multi-format I/O**:
  [`read_geno()`](https://FAkohoue.github.io/LDxBlocks/reference/read_geno.md)
  now supports numeric dosage CSV, HapMap, VCF / bgzipped VCF, SeqArray
  GDS, PLINK BED/BIM/FAM, and plain R matrices through a unified
  `LDxBlocks_backend` S3 interface.
- **GDS streaming backend**: SeqArray GDS files are accessed via
  [`seqSetFilter()`](https://rdrr.io/pkg/SeqArray/man/seqSetFilter.html) +
  `seqGetData("$dosage")` per window, keeping peak RAM proportional to
  one chromosome window rather than the full genome.
- **PLINK BED backend**: `BEDMatrix`-backed streaming access for binary
  PLINK files.
- **[`read_chunk()`](https://FAkohoue.github.io/LDxBlocks/reference/read_chunk.md)**:
  unified genotype slice accessor — works identically for all six
  backend types.
- **[`close_backend()`](https://FAkohoue.github.io/LDxBlocks/reference/close_backend.md)**:
  releases GDS file handles; no-op for in-memory backends.
- **[`run_Big_LD_all_chr()`](https://FAkohoue.github.io/LDxBlocks/reference/run_Big_LD_all_chr.md)
  now accepts `LDxBlocks_backend`** directly, removing the need to load
  the full genome matrix before calling it.
- **`min_snps_chr`** parameter in
  [`run_Big_LD_all_chr()`](https://FAkohoue.github.io/LDxBlocks/reference/run_Big_LD_all_chr.md):
  skip chromosomes or scaffolds with fewer than this many SNPs (default
  10).
- **Example datasets**: `ldx_geno`, `ldx_snp_info`, `ldx_blocks`,
  `ldx_gwas` — simulated with clear block structure for examples and
  tests.
- **Flat-file examples** in `inst/extdata/`: numeric CSV, HapMap, and
  VCF copies of the example data for format-reader tests.

#### Performance improvements

- **`boundary_scan_cpp()`**: C++ replacement for the R inner loop in
  `cutsequence.modi()`. For a 50,000-SNP chromosome with `leng = 200`,
  eliminates ~150,000 small R-level matrix operations per chromosome.
- **`maf_filter_cpp()`**: combined MAF + monomorphic filter in one O(np)
  C++ pass; ~10× faster than the
  [`apply()`](https://rdrr.io/r/base/apply.html) equivalent.
- **`build_adj_matrix_cpp()`**: in-place C++ adjacency construction,
  avoiding the intermediate logical matrix from
  [`ifelse()`](https://rdrr.io/r/base/ifelse.html).
- **`compute_r2_sparse_cpp()`**: sparse r² within a bp distance window,
  returning only above-threshold pairs as a triplet (row, col, r²).
  Avoids the O(p²) dense matrix for large sub-segments where distant
  pairs are always below threshold.

#### Tests

- Three test files: `test-basic.R` (smoke tests via example data),
  `test-cpp.R` (property-based C++ kernel tests), `test-io.R` (all I/O
  format tests via tempfiles).

------------------------------------------------------------------------

### LDxBlocks 0.2.0 (2025-04-07)

#### New features

- **C++/Armadillo core**: `compute_r2_cpp()`, `compute_rV2_cpp()`,
  `maf_filter_cpp()`, `build_adj_matrix_cpp()`, `col_r2_cpp()`,
  `compute_r2_sparse_cpp()`, `boundary_scan_cpp()` via RcppArmadillo.
- **OpenMP parallelism**: `n_threads` parameter for `compute_r2_cpp()`.
- **Dual LD metric**: `method = "r2"` (default) or `"rV2"` in
  [`Big_LD()`](https://FAkohoue.github.io/LDxBlocks/reference/Big_LD.md)
  and
  [`run_Big_LD_all_chr()`](https://FAkohoue.github.io/LDxBlocks/reference/run_Big_LD_all_chr.md).
- **[`prepare_geno()`](https://FAkohoue.github.io/LDxBlocks/reference/prepare_geno.md)**:
  unified preparation step — centres (r²) or centres + whitens (rV²).
- **[`compute_ld()`](https://FAkohoue.github.io/LDxBlocks/reference/compute_ld.md)**:
  dispatcher routing to `compute_r2_cpp()` or `compute_rV2_cpp()`.

#### Bug fixes

- [`CLQD()`](https://FAkohoue.github.io/LDxBlocks/reference/CLQD.md):
  fixed latent bug in `re_idx` tracking after dense-core pre-pass;
  original indices were not correctly mapped back after removal.
- [`Big_LD()`](https://FAkohoue.github.io/LDxBlocks/reference/Big_LD.md):
  [`vapply()`](https://rdrr.io/r/base/lapply.html) used throughout
  instead of [`sapply()`](https://rdrr.io/r/base/lapply.html) for
  type-stability.

------------------------------------------------------------------------

### LDxBlocks 0.1.0 (2025-04-07)

Initial release.

- [`Big_LD()`](https://FAkohoue.github.io/LDxBlocks/reference/Big_LD.md):
  core per-chromosome LD block segmentation (Kim et al. 2018).
- [`CLQD()`](https://FAkohoue.github.io/LDxBlocks/reference/CLQD.md):
  clique detection on rV² adjacency graph.
- [`run_Big_LD_all_chr()`](https://FAkohoue.github.io/LDxBlocks/reference/run_Big_LD_all_chr.md):
  chromosome-wise wrapper.
- [`tune_LD_params()`](https://FAkohoue.github.io/LDxBlocks/reference/tune_LD_params.md):
  grid-search parameter auto-tuner.
- [`extract_haplotypes()`](https://FAkohoue.github.io/LDxBlocks/reference/extract_haplotypes.md),
  [`compute_haplotype_diversity()`](https://FAkohoue.github.io/LDxBlocks/reference/compute_haplotype_diversity.md),
  [`build_haplotype_feature_matrix()`](https://FAkohoue.github.io/LDxBlocks/reference/build_haplotype_feature_matrix.md):
  haplotype analysis module.
- [`summarise_blocks()`](https://FAkohoue.github.io/LDxBlocks/reference/summarise_blocks.md),
  [`plot_ld_blocks()`](https://FAkohoue.github.io/LDxBlocks/reference/plot_ld_blocks.md):
  utilities.
