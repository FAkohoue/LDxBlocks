## R CMD check results
0 errors | 0 warnings | 1 note

* checking for future file timestamps: unable to verify current time
  This is a transient network issue on the check machine and is not
  reproducible. No file in the package has a future timestamp.

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

---
## Notes on optional dependencies
The following `Suggests` packages are optional:
- `AGHmatrix`, `ASRgenomics` -- required only when `method = "rV2"` is
  selected. All default behaviour uses `method = "r2"`.
- `SNPRelate` -- required only for the `"gds"` backend in `read_geno()`,
  and VCF/HapMap auto-conversion to GDS cache.
- `gdsfmt` -- required alongside `SNPRelate` for direct GDS node operations.
  SNPRelate depends on gdsfmt, so it is always available when SNPRelate is
  installed. Listed separately to make the dependency explicit.
- `BEDMatrix` -- required only for the `"bed"` backend in `read_geno()`.
- `future.apply` -- required only when `parallel = TRUE` in `tune_LD_params()`.
- `rrBLUP` -- required by `run_haplotype_prediction()` for REML-based GBLUP.
  Guarded with `requireNamespace("rrBLUP", quietly = TRUE)`.
- `sommer` -- optional multi-trait solver in `run_haplotype_prediction()`.
  When absent or when the multi-trait model fails to converge, the function
  silently falls back to `rrBLUP::kin.blup()` per trait. Guarded with
  `requireNamespace("sommer", quietly = TRUE)`.
- `ggplot2` -- required only for `plot_ld_blocks()`.
- `usethis` -- used only in `data-raw/generate_example_data.R`, which is
  excluded from the tarball via `.Rbuildignore`.
- `knitr`, `rmarkdown` -- required to build vignettes.
- `testthat` -- required to run the test suite.

All `Suggests` packages are guarded with `requireNamespace(..., quietly = TRUE)`
before use. No `Suggests` package is loaded unconditionally.

---
## Notes on example data
The four datasets in `data/` (`ldx_geno`, `ldx_snp_info`, `ldx_blocks`,
`ldx_gwas`) are small simulated objects (120 individuals, 230 SNPs, 3
chromosomes) stored as compressed `.rda` files (`compress = "xz"`) and
occupy less than 100 KB combined. The `data-raw/` directory is excluded
from the tarball via `.Rbuildignore`.

---
## Downstream dependencies
None. This is a new package.
