## R CMD check results
0 errors ✔ | 0 warnings ✔ | 0 notes ✔

* This is a new submission.

---
## Test environments
* Local: Windows 11, R 4.5.0, Rtools45
* GitHub Actions:
  - ubuntu-latest, R release
  - ubuntu-latest, R devel
  - ubuntu-latest, R oldrel-1
  - macos-latest,  R release
  - windows-latest, R release

---
## Notes on compiled code (C++ / OpenMP)
The package compiles C++ source code via Rcpp and RcppArmadillo.

**Linux / macOS** (`src/Makevars`):
    PKG_CXXFLAGS = $(SHLIB_OPENMP_CXXFLAGS)
    PKG_LIBS     = $(SHLIB_OPENMP_CXXFLAGS) $(LAPACK_LIBS) $(BLAS_LIBS) $(FLIBS)

**Windows** (`src/Makevars.win`):
    PKG_CXXFLAGS = $(SHLIB_OPENMP_CXXFLAGS)
    PKG_LIBS     = $(SHLIB_OPENMP_CXXFLAGS) $(LAPACK_LIBS) $(BLAS_LIBS) $(FLIBS)

Notes:
- No non-portable flags are used. `-march=native`, `-O2`, and
  `-funroll-loops` are absent from both Makevars files.
- OpenMP is enabled via `$(SHLIB_OPENMP_CXXFLAGS)`, which expands to the
  correct platform flag (`-fopenmp` on Linux/Windows, empty on macOS when
  OpenMP is unavailable). The package degrades gracefully to single-threaded
  execution -- every OpenMP call is guarded with `#ifdef _OPENMP`.
- LAPACK and BLAS are linked via `$(LAPACK_LIBS)` and `$(BLAS_LIBS)`,
  which resolve to whatever BLAS/LAPACK R itself was built against.

**Exported C++ functions (10 total, all in `src/ld_core.cpp`, 1,053 lines):**
- `compute_r2_cpp()` — standard r² matrix, OpenMP outer loop
- `compute_rV2_cpp()` — kinship-adjusted rV² (same kernel as compute_r2_cpp)
- `maf_filter_cpp()` — MAF + monomorphic filter, single O(np) pass
- `build_adj_matrix_cpp()` — threshold LD matrix to 0/1 adjacency
- `col_r2_cpp()` — r² of one query column against all others
- `compute_r2_sparse_cpp()` — sparse r² within bp distance window, OpenMP
- `boundary_scan_cpp()` — weak-LD cut position scan for subsegmentation
- `build_hap_strings_cpp()` — C++ haplotype string builder (replaces R vapply loop)
- `resolve_overlap_cpp()` — LD-informed block overlap resolution,
    BLAS DGEMM scoring + lazy column cache + OpenMP over pairs
- `block_snp_ranges_cpp()` — O(p + n_blocks) single linear sweep mapping
    all chromosome blocks to SNP index ranges; replaces per-block findInterval()
- `extract_chr_haplotypes_cpp()` — full chromosome haplotype extractor
    combining B7 (interval lookup + string building with OpenMP parallel for),
    B9 (C++ unordered_map frequency tabulation with min_freq/top_n filtering
    and freq_dominant = sorted_cnt[0].first / total), and B10 (allele strings,
    frequencies, and counts returned directly without R-side table() calls)

Two additional `static` (non-exported) helpers:
- `score_overlap_cpp()` — BLAS DGEMM overlap scoring kernel
- `resolve_seam_cpp()` — seam-local resolver for future pangenome-scale use

---
## Notes on optional dependencies
The following `Suggests` packages are used conditionally at runtime:

- `AGHmatrix`, `ASRgenomics` -- required only when `method = "rV2"` is
  selected in `prepare_geno()` and `run_Big_LD_all_chr()`. All default
  behaviour uses `method = "r2"`. Both are guarded with
  `requireNamespace(..., quietly = TRUE)` at the entry of the rV2 code path.

- `SNPRelate`, `gdsfmt` -- required only for the "gds" backend in
  `read_geno()`, and for automatic VCF/HapMap-to-GDS cache conversion.
  SNPRelate is guarded with `requireNamespace("SNPRelate", quietly = TRUE)`.
  gdsfmt is always available when SNPRelate is installed (SNPRelate depends
  on gdsfmt); listed separately to make the dependency explicit. Both are
  accessed only inside the GDS reader path via the internal `.require_pkg()`
  helper, which calls `requireNamespace()` and stops with an informative
  install message if the package is absent.

- `BEDMatrix` -- required only for the "bed" backend in `read_geno()`.
  Guarded via `.require_pkg("BEDMatrix", ...)` at the BED reader entry.

- `bigmemory` -- required only for `read_geno_bigmemory()` and when
  `use_bigmemory = TRUE` in `run_ldx_pipeline()`. Guarded with
  `requireNamespace("bigmemory", quietly = TRUE)`.

- `future.apply` -- required only when `parallel = TRUE` in
  `tune_LD_params()`. Guarded with
  `requireNamespace("future.apply", quietly = TRUE)`.

- `ggplot2` -- required only for `plot_ld_blocks()` and `plot_ld_decay()`.
  Guarded with `requireNamespace("ggplot2", quietly = TRUE)`.

- `usethis` -- used only in `data-raw/generate_example_data.R` (the script
  that generates bundled example datasets). This script is excluded from the
  tarball via `.Rbuildignore` and is never executed during installation or
  package checks.

- `knitr`, `rmarkdown` -- required to build the vignette. Not used at runtime.

- `testthat` -- required to run the test suite. Not used at runtime.

All `Suggests` packages are either guarded with
`requireNamespace(..., quietly = TRUE)` before use, or are needed only for
non-runtime tasks (vignette building, test suite, dataset generation script).
No `Suggests` package is loaded unconditionally at runtime.

---
## Notes on example data
Six datasets are bundled in `data/`: `ldx_geno`, `ldx_snp_info`, `ldx_blocks`,
`ldx_gwas`, `ldx_blues`, and `ldx_blues_list`. All are small simulated objects
(120 individuals, 230 SNPs, 3 chromosomes) stored as compressed .rda files
(`compress = "xz"`), occupying less than 200 KB combined. The `data-raw/`
directory (dataset generation scripts) is excluded from the tarball via
`.Rbuildignore`.

- `ldx_geno` / `ldx_snp_info` / `ldx_blocks` -- genotype matrix (120 x 230),
  SNP metadata, and pre-computed block table (9 blocks, 3 chromosomes) for
  block detection examples and tests.
- `ldx_gwas` -- 20 toy GWAS markers with p-values and a `trait` column
  ("TraitA" / "TraitB") for `tune_LD_params()` and `define_qtl_regions()`
  examples including pleiotropic block detection.
- `ldx_blues` -- pre-adjusted phenotype means (BLUEs) for two simulated traits
  (YLD, RES) across all 120 individuals, for `run_haplotype_prediction()`,
  `test_block_haplotypes()`, and `estimate_diplotype_effects()` examples.
- `ldx_blues_list` -- named list of two environments (env1, env2) of named
  numeric BLUEs (120 individuals each) for `run_haplotype_stability()` and
  multi-environment cross-validation examples.

---
## Notes on LD decay analysis
`compute_ld_decay()` estimates chromosome-specific r-squared decay distances
using a memory-efficient position-first sampling strategy: pair indices are
drawn from SNP positions (no genotype I/O), then `read_chunk()` loads only the
unique columns needed. Two threshold approaches are supported:
- Fixed numeric threshold (e.g. 0.1, standard GWAS practice).
- Parametric threshold: 95th percentile of r-squared between unlinked markers
  on different chromosomes, measuring background kinship-induced LD.
Optional curve fitting: LOESS smoothing or Hill-Weir (1988) nonlinear model.
Decay distances integrate directly with `define_qtl_regions(ld_decay = decay)`
for LD-aware candidate gene windows around GWAS hits. Censored distances
(threshold not crossed within `max_dist`) are flagged with `censored = TRUE`
and emitted as warnings.

---
## Notes on haplotype association testing
`test_block_haplotypes()` implements a unified Q+K mixed linear model
(EMMAX/GAPIT3 formulation). Population structure PCs are derived from the
eigendecomposition of the haplotype GRM, ensuring mathematical consistency
between fixed-effect covariates and the kinship random effect. The GRM is
inverted once per trait via `rrBLUP::mixed.solve()` (dominant cost, O(n^3));
per-allele tests across all blocks are then vectorised in a single BLAS
crossprod() call, analogous to marginal SNP screening. Setting `n_pcs = 0`
(default) gives pure GRM correction (EMMAX); `n_pcs > 0` enables the Q+K model.

---
## Downstream dependencies
None. This is a new submission.
