# LDxBlocks — Genome-Wide LD Block Detection, Haplotype Analysis, and Genomic Prediction Features

<p align="center">
  <img src="man/figures/logo.png" alt="LDxBlocks logo" width="190px">
</p>

<!-- badges: start -->
[![R-CMD-check](https://github.com/FAkohoue/LDxBlocks/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/FAkohoue/LDxBlocks/actions/workflows/R-CMD-check.yaml)
[![pkgdown](https://github.com/FAkohoue/LDxBlocks/actions/workflows/pkgdown.yaml/badge.svg)](https://github.com/FAkohoue/LDxBlocks/actions/workflows/pkgdown.yaml)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
<!-- badges: end -->

---

## Motivation

Linkage disequilibrium (LD) block detection is a foundational step in modern
genomic analyses. Knowing which SNPs co-segregate as a unit determines how GWAS
results are interpreted, how haplotypes are defined for population genetics, and
how genomic prediction models should be structured. Despite its importance, most
implementations of LD block detection suffer from three problems that limit
their usefulness in practice.

**The kinship problem.** Classical LD estimators (r²) assume independence
between individuals. In livestock, crop, or family-based human cohorts this
assumption is systematically violated. Cryptic relatedness inflates pairwise
correlations, causing LD to appear stronger than it is and blocks to be drawn
too broadly or in the wrong places. The kinship-adjusted squared correlation
rV² (Mangin et al. 2012, *Heredity* 108:285-291) corrects for this by
whitening the genotype matrix with the inverse square root of the genomic
relationship matrix (GRM):

$$rV^2_{jk} = \left[\mathrm{Cor}(V^{-1/2} G_j,\; V^{-1/2} G_k)\right]^2$$

where $V$ is the VanRaden (2008) GRM and $V^{-1/2}$ is its Cholesky whitening
factor. Every pairwise correlation in the Big-LD algorithm is replaced with
rV², giving block boundaries that reflect true recombination structure rather
than population-structure artefacts.

**The scale problem.** The original Big-LD implementation (Kim et al. 2018) contains
no compiled code -- every matrix operation runs through the R interpreter. For modern whole-genome sequencing panels with 2-10 million markers the
inner loops are prohibitively slow, and loading the full genotype matrix
before detection is impossible on most workstations. LDxBlocks addresses this
with a C++/Armadillo computational core compiled via Rcpp, OpenMP-parallelised
LD computation, a unified multi-format I/O layer, and a strict never-full-genome
memory model: the full genotype matrix is never loaded into RAM simultaneously
regardless of dataset size or format.

**The pipeline gap.** The original Big-LD stops at block boundaries. In
breeding programmes and GWAS follow-up studies the immediate next questions
are: What haplotypes exist within each block? How diverse are they? Which
blocks harbour QTLs? How do I build a genomic prediction model that captures
multi-locus block effects? These questions require a coherent downstream
pipeline that the original implementation does not provide.

**Why LDxBlocks and not other tools?**

- Unlike PLINK `--blocks`, LDxBlocks uses the clique-based Big-LD segmentation
  algorithm which does not require a pre-specified number of blocks and
  naturally handles complex LD patterns. It also supports kinship correction
  (rV²) and six input formats natively.
- Unlike the original Big-LD R package, LDxBlocks has a compiled C++ core
  (approximately 40x faster for typical window sizes), streams genotype data
  from disk without loading the full genome into RAM, adds rV² kinship
  correction, handles singleton SNPs explicitly, and extends the pipeline into
  haplotype analysis and genomic prediction.
- Unlike LDstore2 or LDpred2, LDxBlocks produces interpretable genomic
  intervals (start/end position, rsID) as output — not just summary statistics
  — making blocks immediately usable for annotation, diversity analysis, and
  region-based modelling.
- Unlike gpart (the Bioconductor successor to Big-LD), LDxBlocks integrates
  directly with R-based genomic prediction workflows (rrBLUP, BGLR, ASReml-R)
  through its haplotype feature matrix module, supports automatic GDS conversion
  for WGS-scale datasets, and provides GWAS-driven parameter auto-tuning.

---

## Summary

`LDxBlocks` is a complete pipeline for genome-wide LD block detection in
related or structured populations, with downstream haplotype analysis and
genomic prediction utilities. It extends the Big-LD clique-based segmentation
algorithm of Kim et al. (2018) with 12 concrete improvements:

1. **Dual LD metric** — standard r² (default, fast, no kinship correction,
   suitable for 10 M+ markers) and kinship-adjusted rV² (for structured or
   related populations) selectable via a single `method =` argument.
2. **C++/Armadillo core** — seven compiled functions handle r² computation,
   adjacency matrix construction, MAF filtering, boundary scanning, sparse LD,
   and column-wise correlation. The full r² matrix for a 1,500-SNP window
   computes in milliseconds rather than seconds.
3. **OpenMP parallelism** — the outer loop of `compute_r2_cpp()` is
   parallelised across threads; count controlled with `n_threads =`.
4. **Unified multi-format I/O** — `read_geno()` auto-detects and reads numeric
   dosage CSV, HapMap, VCF/VCF.gz, SNPRelate GDS, PLINK BED/BIM/FAM, and
   in-memory R matrices through a single entry point with a common backend
   interface (`LDxBlocks_backend`).
5. **Never-full-genome memory model** — the full genotype matrix is never held
   in RAM at once for any format. Numeric dosage CSV is read in pre-allocated
   50,000-row chunks (peak RAM = one chunk, not 2× the file). VCF and HapMap
   auto-convert to a streaming GDS cache. GDS and PLINK BED backends load only
   the SNP window per CLQD call. `gc(FALSE)` is called after each chromosome
   to prevent heap fragmentation across 20–30 chromosome passes.
6. **MAF filter in C++** — `maf_filter_cpp()` runs in a single O(np) pass,
   handling NA imputation and monomorphic detection simultaneously.
7. **C++ boundary scan** — `boundary_scan_cpp()` replaces the R inner loop in
   the subsegmentation step that runs hundreds of times per chromosome,
   eliminating interpreter overhead for the most-called function in the pipeline.
8. **Sparse r² computation** — `compute_r2_sparse_cpp()` computes pairwise r²
   only for SNP pairs within a user-specified bp distance, avoiding O(p²) cost
   for large sub-segments where distant SNPs will always be below threshold.
9. **Automatic parameter tuning** — `tune_LD_params()` performs a grid search
   over CLQcut and other parameters, selecting the combination that minimises
   unassigned and forced GWAS marker placements.
10. **Haplotype reconstruction** — `extract_haplotypes()` constructs phase-free
    diploid allele strings per LD block and individual, which are nearly 1:1
    with true haplotype classes within high-LD blocks.
11. **Diversity metrics** — `compute_haplotype_diversity()` computes per-block
    richness, Nei's expected heterozygosity, Shannon entropy, and dominant
    haplotype frequency.
12. **Prediction feature matrix** — `build_haplotype_feature_matrix()` converts
    haplotype strings to numeric dosage columns for GBLUP, BayesB, or machine
    learning models, capturing multi-locus epistatic effects implicitly.

---

## Relationship to the original Big-LD algorithm

LDxBlocks is built on the clique-based segmentation algorithm of Kim et al.
(2018), published as the `BigLD` R package and later updated as the `gpart`
Bioconductor package. The mathematical core — interval graph modelling of LD
bins, maximum-weight independent set block construction, and Bron-Kerbosch
clique enumeration via igraph — is preserved exactly. LDxBlocks extends that
foundation in the following concrete ways.

### Computational core: R loops replaced by C++

The original `Big_LD()` calls `cor()` inside the boundary-scan loop — up to
three separate R-level matrix operations per candidate cut position. For a
50,000-SNP chromosome with `leng = 200` this is approximately 150,000 small
matrix multiplications running through the R interpreter. LDxBlocks replaces
these with seven compiled functions:

| Original R operation | LDxBlocks C++ function | Speedup |
|---|---|---|
| `cor(subgeno)` per CLQD call | `compute_r2_cpp()` + OpenMP | ~40x for 1,500-SNP window |
| `apply(Ogeno, 2, ...)` MAF filter | `maf_filter_cpp()` single pass | ~10x for 100k+ SNPs |
| `cor()` inside boundary-scan loop | `boundary_scan_cpp()` compiled | ~20x per chromosome |
| `r2Mat[r2Mat >= CLQcut^2] <- 1` | `build_adj_matrix_cpp()` | eliminates intermediate allocation |
| Single-column correlation | `col_r2_cpp()` | used in boundary scan helper |
| Sparse within-window r² | `compute_r2_sparse_cpp()` | avoids O(p^2) for large segments |

The outer loop of `compute_r2_cpp()` is parallelised with OpenMP, controlled
by `n_threads =`. Thread scaling is efficient up to 8-16 threads for typical
`subSegmSize = 1500` windows.

### Memory model: never-full-genome

The original `Big_LD()` accepts a plain R matrix — the entire genotype dataset
must fit in RAM simultaneously. For a 10M-marker, 5,000-individual WGS panel
this requires approximately 400 GB as an R `double` matrix. LDxBlocks enforces
a strict memory contract regardless of format:

- **Numeric dosage CSV**: two-pass chunked pre-allocated reader. Pass 1 counts
  rows with zero data loaded. A single pre-allocated matrix is filled in
  50,000-row chunks via successive `fread()` calls. Peak RAM = one chunk, never
  2x the file.
- **VCF and HapMap**: auto-converted to a SNPRelate GDS cache on first call.
  All subsequent access streams per chromosome window via `read_chunk()`.
- **GDS and PLINK BED**: `read_chunk(backend, col_idx)` is called once per
  sub-segment per chromosome. Only those columns are loaded; the rest of the
  genome remains on disk.
- **Chromosome loop**: `rm()` and `gc(FALSE)` are called after each chromosome
  completes, preventing heap fragmentation across 20-30 chromosome passes.

### Kinship correction: rV²

The original implementation uses Pearson r as the LD metric. For related
populations (livestock half-sib families, inbred plant lines, family-based
human cohorts) kinship-induced allele sharing inflates r between all SNP pairs,
producing blocks that are too broad and incorrectly delimited.

LDxBlocks adds `method = "rV2"` which replaces every pairwise correlation with
the kinship-whitened equivalent. The whitening factor A is computed once per
chromosome from the VanRaden (2008) GRM via `get_V_inv_sqrt()` (Cholesky or
eigendecomposition), then applied to the centred genotype matrix before passing
to the same `compute_r2_cpp()` kernel. In related populations, rV² blocks are
typically 10-30% smaller and more precisely delimited than r² blocks.

### Singleton SNP handling

The original `Big_LD()` silently discards SNPs that pass MAF filtering but
receive `NA` from `CLQD()` (no clique partner above `CLQcut`). These
singletons — which mark recombination hotspots and rapidly-evolving loci — are
never returned to the user and are invisible in the block table.

LDxBlocks adds `singleton_as_block = TRUE` which collects singleton indices
during the sub-segment loop and appends them to the final block table as
single-SNP entries (`start == end`, `length_bp == 1`). With the default
`singleton_as_block = FALSE` the original behaviour is preserved for backward
compatibility.

### Bug fix: zero-row assignment

The original main loop body contains:

```r
LDblocks[(preleng1 + 1):(preleng1 + dim(nowLDblocks)[1]), ] <- nowLDblocks
```

When a sub-segment contains no valid cliques (all singletons above the MAF
threshold) `nowLDblocks` has zero rows and R evaluates `(preleng1+1):preleng1`
as a backwards sequence, causing `replacement has length zero`. LDxBlocks wraps
this with `if (nrow(nowLD) > 0L)`, making the function robust to sparse
chromosomal regions and small input datasets.

### Downstream pipeline

The original Big-LD stops at the block table. LDxBlocks adds a complete
downstream module:

| Capability | Original Big-LD / gpart | LDxBlocks |
|---|---|---|
| Block detection | Yes (core algorithm) | Yes (same algorithm + C++ + rV²) |
| Statistical phasing | No | `phase_with_beagle()`, `phase_with_pedigree()`, `read_phased_vcf()` |
| Haplotype extraction | No | `extract_haplotypes()` — phased and unphased, backend streaming |
| Diversity metrics | No | `compute_haplotype_diversity()` — He, Shannon, richness, f_max |
| Post-GWAS QTL mapping | No | `define_qtl_regions()` — pleiotropic block detection |
| Genomic prediction features | No | `build_haplotype_feature_matrix()` — additive 0/1/2 or presence/absence |
| Output writers | No | Numeric CSV, HapMap, diversity CSV — compatible with TASSEL, GAPIT, rrBLUP |
| Parameter auto-tuning | No | `tune_LD_params()` — grid search against GWAS marker coverage |
| Multi-format I/O | PLINK, VCF (gpart) | Numeric CSV, HapMap, VCF, GDS, BED, R matrix via unified backend |
| WGS-scale streaming | Partial (gpart GDS) | Full never-full-genome model for all formats |

### What is kept exactly

- The interval graph modelling of LD bins (cliques of strong pairwise LD SNPs)
- `CLQD()`: bin vector assignment via maximal clique enumeration and greedy
  density-priority selection
- `constructLDblock()`: maximum-weight independent set via dynamic programming
  on sorted interval sequences
- `appendSGTs()`: rare-SNP appending logic
- `cutsequence.modi()`: boundary-scan logic and forced-split fall-back
- All `CLQmode = "Density"` and `CLQmode = "Maximal"` clique scoring
- The `clstgap` physical distance splitting within cliques
- Block table column format (`start`, `end`, `start.rsID`, `end.rsID`,
  `start.bp`, `end.bp`) for drop-in compatibility with downstream tools that
  accept original Big-LD output

---

## Table of contents

1. [Installation](#installation)
2. [Quick start](#quick-start)
3. [Input formats](#input-formats)
4. [Statistical background](#statistical-background)
   - [MAF filtering](#1-maf-filtering)
   - [Genotype preparation](#2-genotype-preparation)
   - [Subsegmentation](#3-subsegmentation)
   - [Clique detection](#4-clique-detection-clqd)
   - [Block construction](#5-block-construction)
5. [LD metrics: r² versus rV²](#ld-metrics-r-versus-rv)
6. [Haplotype analysis](#haplotype-analysis)
7. [Parameter auto-tuning](#parameter-auto-tuning)
8. [Scale strategies and backends](#scale-strategies-and-backends)
9. [Full pipeline walkthrough](#full-pipeline-walkthrough)
10. [Function reference](#function-reference)
11. [Output objects](#output-objects)
12. [Memory and performance notes](#memory-and-performance-notes)
13. [Documentation](#documentation)
14. [Citation](#citation)
15. [Contributing](#contributing)
16. [License](#license)
17. [References](#references)

---

## Installation

Install from GitHub with vignettes (recommended):

```r
install.packages("remotes")
remotes::install_github("FAkohoue/LDxBlocks",
  build_vignettes = TRUE,
  dependencies    = TRUE
)
```

Install without vignettes for a faster build:

```r
remotes::install_github("FAkohoue/LDxBlocks",
  build_vignettes = FALSE,
  dependencies    = TRUE
)
```

**Required dependencies** (installed automatically):

```r
install.packages(c(
  "Rcpp",            # C++ interface
  "RcppArmadillo",   # linear algebra C++ library
  "igraph",          # maximal clique detection
  "data.table",      # fast file I/O and chromosome operations
  "dplyr"            # pipeline utilities
))
```

**Optional dependencies** for additional backends and features:

```r
# GDS backend — required for .gds files; recommended for panels > 2 M SNPs
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("SNPRelate")

# PLINK BED backend — required for .bed/.bim/.fam input
install.packages("BEDMatrix")

# Kinship-adjusted rV² (method = "rV2") — structured/related populations
install.packages(c("AGHmatrix", "ASRgenomics"))

# Parallel parameter tuning
install.packages("future.apply")

# LD block visualisation
install.packages("ggplot2")
```

> **C++ compilation note.** LDxBlocks compiles C++ with OpenMP enabled via
> `$(SHLIB_OPENMP_CXXFLAGS)` in `src/Makevars` and `src/Makevars.win`. R's
> default `CXXFLAGS` already include `-O2` optimisation. No non-portable flags
> are used, so the package passes `R CMD check --as-cran` without notes.

---

## Quick start

```r
library(LDxBlocks)

# ── Option A: from files (all formats accepted) ───────────────────────────────
be <- read_geno("mydata.vcf.gz")    # auto-detected as VCF
be
# LDxBlocks backend
#   Type    : vcf
#   Samples : 500
#   SNPs    : 850000
#   Chr     : 1, 2, 3, ..., 12

blocks <- run_Big_LD_all_chr(
  be,
  method    = "r2",    # standard r² — fast, no kinship needed
  CLQcut    = 0.70,
  n_threads = 8L
)

head(blocks)
#   start  end start.rsID end.rsID start.bp   end.bp CHR length_bp
# 1     1   18      rs001    rs018     1000   103000   1    102001
# 2    22   41      rs022    rs041   125000   228000   1    103001

summarise_blocks(blocks)
#    CHR n_blocks min_bp median_bp  mean_bp  max_bp total_bp_covered
# 1    1      124   2000     85000    91432  420000        11337568
# ...
# 13 GENOME  1502   2000     82000    88940  512000       133548316

close_backend(be)
```

```r
# ── Option B: plain R matrix (backward-compatible, unchanged API) ─────────────
blocks <- run_Big_LD_all_chr(
  my_geno_matrix,
  snp_info  = my_snp_info,   # data.frame: SNP, CHR, POS
  method    = "r2",
  CLQcut    = 0.70
)
```

```r
# ── Haplotype analysis — pass the backend directly for chromosome streaming ───
# extract_haplotypes() detects LDxBlocks_backend input and processes one
# chromosome at a time, freeing RAM before moving to the next.
haps <- extract_haplotypes(be, be$snp_info, blocks, min_snps = 3)

div <- compute_haplotype_diversity(haps)
head(div)
#            block_id n_ind n_haplotypes    He Shannon freq_dominant
# block_1000_103000     500           12 0.891   3.142         0.182

feat <- build_haplotype_feature_matrix(haps, top_n = 5, scale_features = TRUE)
dim(feat)   # 500 x (n_blocks * 5)
```

```r
# ── Parameter auto-tuning against GWAS markers ────────────────────────────────
result <- tune_LD_params(
  geno_matrix    = my_geno,
  snp_info       = my_snp_info,
  gwas_df        = my_gwas,      # data.frame: Marker, CHR, POS
  prefer_perfect = TRUE
)

result$best_params
result$score_table
result$gwas_assigned   # every GWAS marker assigned to a block
```

---

## Input formats

`read_geno()` accepts a path to a genotype file (or an in-memory R matrix)
and an optional `format =` override. The format is auto-detected from the
file extension when `format` is not supplied.

Phenotype data is not an input to block detection. If you need to align sample
IDs between a genotype file and a phenotype table — for example, to
subset to phenotyped individuals before running the pipeline — see the
`read_pheno()` and `align_geno_pheno()` utilities documented under
[Utility functions](#utilities).

### Genotype file

| Format | Extension | `format =` |
|--------|-----------|------------|
| Numeric dosage | `.csv`, `.txt` | `"numeric"` |
| HapMap | `.hmp.txt` | `"hapmap"` |
| VCF / bgzipped VCF | `.vcf`, `.vcf.gz` | `"vcf"` |
| SNPRelate GDS | `.gds` | `"gds"` |
| PLINK binary | `.bed` (+ `.bim`, `.fam`) | `"bed"` |
| R matrix | (in-memory) | `"matrix"` |

**Numeric dosage format** — one row per SNP; columns `SNP`, `CHR`, `POS`,
`REF`, `ALT` followed by one column per sample, values in `{0, 1, 2, NA}`.
`REF` and `ALT` are optional and filled with `NA` if absent:

| SNP | CHR | POS | REF | ALT | Line01 | Line02 | … |
|-----|-----|-----|-----|-----|--------|--------|---|
| SNP001 | 1 | 10000 | A | T | 0 | 1 | … |
| SNP002 | 1 | 20000 | G | C | 2 | 0 | … |

**HapMap format** — standard 11-column header (`rs#`, `alleles`, `chrom`,
`pos`, `strand`, `assembly#`, `center`, `protLSID`, `assayLSID`, `panelLSID`,
`QCcode`) followed by sample columns with two-character nucleotide calls.
`AA`, `AT`, `TT`, `NN` (missing) are all accepted. Dosage is decoded from the
alleles column — e.g. `alleles = A/T` and call `AT` → heterozygous ALT → 1:

| rs# | alleles | chrom | pos | strand | … | QCcode | Line01 | Line02 | … |
|-----|---------|-------|-----|--------|---|--------|--------|--------|---|
| SNP001 | A/T | 1 | 10000 | + | … | NA | AA | AT | … |
| SNP002 | G/C | 1 | 20000 | + | … | NA | CC | GG | … |
| SNP003 | C/G | 1 | 30000 | + | … | NA | CG | GG | … |

**VCF** — standard VCF v4.2. Both phased (`0|1`) and unphased (`0/1`) GT
fields are accepted. Multi-allelic sites use the first ALT allele. Missing
calls (`./.`) become `NA`. The `##` meta-information lines are skipped
automatically:

| #CHROM | POS | ID | REF | ALT | QUAL | FILTER | INFO | FORMAT | Line01 | Line02 | … |
|--------|-----|----|-----|-----|------|--------|------|--------|--------|--------|---|
| 1 | 10000 | SNP001 | A | T | . | PASS | . | GT | 0/0 | 0/1 | … |
| 1 | 20000 | SNP002 | G | C | . | PASS | . | GT | 1/1 | 0/0 | … |
| 1 | 30000 | SNP003 | C | G | . | PASS | . | GT | 0/1 | 1/1 | … |

**GDS** — SNPRelate GDS file. Requires `BiocManager::install("SNPRelate")`.
Chunk access is via `SNPRelate::snpgdsGetGeno()` with explicit `snp.id` and
`sample.id` vectors — the full genome matrix is never held in RAM simultaneously.

**PLINK BED** — binary PLINK format. The companion `.bim` and `.fam` files must
exist at the same path stem. Requires `install.packages("BEDMatrix")`. Row and
column indexing are memory-mapped via `BEDMatrix`.

**Chromosome normalisation.** At read time, all chromosome labels are stripped
of leading `chr`, `Chr`, or `CHR` so that `chr1`, `Chr01`, and `1` are all
stored as `"1"` in `snp_info$CHR`. This applies to all six formats and ensures
correct cross-format matching.


---

## Statistical background

### 1. MAF filtering

The ALT allele frequency is estimated from the dosage matrix as:

$$\mathrm{AF}_i = \frac{\sum_j g_{ij}}{2\, n_j}$$

where $g_{ij} \in \{0, 1, 2, \mathrm{NA}\}$ is the dosage for SNP $i$ in
sample $j$ and $n_j$ is the number of non-missing observations. The minor
allele frequency is:

$$\mathrm{MAF}_i = \min\!\left(\mathrm{AF}_i,\; 1 - \mathrm{AF}_i\right)$$

SNPs with $\mathrm{MAF}_i < \tau_{\mathrm{maf}}$ (default 0.05) are removed
before any LD computation. Monomorphic SNPs (all observations identical) are
removed unconditionally regardless of the MAF threshold. Both operations run in
a single O(np) C++ pass by `maf_filter_cpp()`, handling NA imputation in the
same loop.

### 2. Genotype preparation

Before computing LD, the genotype matrix $G$ (individuals × SNPs) is centred
and optionally whitened depending on the chosen LD metric. This step is
performed by `prepare_geno()`.

**Standard r² path (`method = "r2"`):**

$$\tilde{G}_j = G_j - \bar{G}_j$$

where $\bar{G}_j$ is the column mean for SNP $j$. No kinship matrix is
required. Runs in O(np) inside the same C++ kernel as `compute_r2_cpp()`.

**Kinship-adjusted rV² path (`method = "rV2"`):**

1. Compute the VanRaden (2008) GRM: $V = \frac{\tilde{G}\tilde{G}^\top}{2 \sum_j p_j(1-p_j)}$ via `AGHmatrix::Gmatrix()`.
2. Bend and condition-number tune $V$ via `ASRgenomics::G.tuneup(bend = TRUE, rcn = TRUE)` to ensure positive-definiteness.
3. Compute the whitening factor:
   - `kin_method = "chol"` (default): $A = R^{-1}$ where $V = R^\top R$ via Cholesky. Fast and numerically stable.
   - `kin_method = "eigen"`: $A = Q\Lambda^{-1/2}Q^\top$. Symmetric; eigenvalues floored at $10^{-6}$ for stability. Preferred for near-singular matrices.
4. Apply whitening: $X = A\tilde{G}$ — used in all downstream LD computations.

Both paths expose the same `compute_r2_cpp()` C++ kernel; the distinction is
only in the preparation step.

### 3. Subsegmentation

For chromosomes with more SNPs than `subSegmSize` (default 1,500), the
algorithm first identifies weak-LD boundary positions to divide the chromosome
into manageable sub-segments, avoiding the O(p²) cost of computing the full
chromosome LD matrix.

For each candidate cut position $i$, the maximum cross-boundary squared
correlation between a left window $L$ and right window $R$ of half-width
`leng` (default 200 SNPs) is evaluated:

$$\text{cut}(i) = 1 \iff \max_{j \in L,\; k \in R} r^2(j, k) < \tau_{\mathrm{CLQ}}$$

This scan is implemented by `boundary_scan_cpp()`: a C++ function that
pre-standardises all columns once, then iterates over every candidate cut
position in compiled code. For a chromosome with 50,000 SNPs and `leng = 200`,
this eliminates approximately 150,000 small matrix operations that would
otherwise run through the R interpreter.

Three candidate window sizes (1, 10, `leng` SNPs each side) are tested at each
position in order of increasing cost. The first window size that shows
cross-boundary LD moves on immediately; only if all three show no LD is a
cut-point declared.

If no weak-LD boundaries are found within $5 \times$ `subSegmSize` SNPs, the
algorithm switches to forced equal-size splits placed at the minimum-LD
positions within each oversized segment, identified by a secondary scan with
a narrow tick window of `leng / 5` SNPs each side.

### 4. Clique detection (CLQD)

Within each sub-segment, an undirected graph is constructed where SNPs are
nodes and edges connect pairs with $r^2 \geq \tau_{\mathrm{CLQ}}$ (`CLQcut`,
default 0.5). The adjacency matrix is built by `build_adj_matrix_cpp()` in a
single O(p²) compiled pass.

Maximal cliques in this graph are enumerated by `igraph::max_cliques()`, which
wraps the Bron-Kerbosch algorithm. Each maximal clique represents a set of SNPs
in mutual high LD. The cliques are prioritised by a score:

$$\text{score}(\mathcal{K}) = \begin{cases} |\mathcal{K}| / (\mathrm{span_{kb}}(\mathcal{K}) + 1) & \text{Density mode (default)} \\ |\mathcal{K}| & \text{Maximal mode} \end{cases}$$

**Density mode** (`CLQmode = "Density"`) prefers compact, high-density cliques
that correspond to biologically meaningful LD blocks. **Maximal mode**
(`CLQmode = "Maximal"`) prefers the largest cliques regardless of span —
useful when the goal is to retain as many co-inherited SNPs as possible.

A greedy assignment iteratively selects the highest-scoring clique and removes
its SNPs from subsequent rounds until all SNPs are assigned or no cliques
remain. The result is an integer bin vector assigning each SNP to a clique-bin
or `NA` (singleton).

**Large window optimisation.** When `checkLargest = TRUE` and a sub-segment
has ≥ 500 SNPs, a dense-core decomposition pre-pass runs before `max_cliques()`.
The graph coreness is computed; if the largest connected component has median
coreness > 80 and maximum coreness > 100 (indicative of a near-complete
subgraph that would cause exponential blowup in clique enumeration), the
hub SNP and its immediate neighbourhood are extracted as a bin and removed
from the graph before clique detection proceeds.

**Clique splitting.** When `split = TRUE`, cliques whose SNP positions span
more than `clstgap` base pairs (default 40,000) are split at the largest
internal gap. This prevents a single bin from containing biologically unrelated
SNPs that happen to share high r² due to long-range LD.

### 5. Block construction

Bin assignments from CLQD are converted to genomic intervals by a
maximum-weight independent set (MWIS) procedure operating on the **interval
graph** of overlapping clique ranges:

1. For each bin, compute its positional range [start_bp, end_bp].
2. Build an interval graph where bins sharing any positional overlap are
   adjacent, weighted by size (number of SNPs covered).
3. Find the MWIS — the largest collection of non-overlapping bins by a
   dynamic programming algorithm on the sorted interval sequence.
4. The selected bins become block intervals. Remaining unassigned SNPs are
   fed back as a new CLQD round until no bins remain.

After MWIS, overlapping blocks are merged greedily. Blocks are then re-indexed
against the full SNP set including monomorphic SNPs removed at step 4, so
that `start` and `end` columns count from 1 over all SNPs in the original input.

---

## LD metrics: r² versus rV²

The choice between `method = "r2"` and `method = "rV2"` is the most important
modelling decision in LDxBlocks. The table below summarises when each is
appropriate.

| | r² | rV² |
|---|---|---|
| **Population type** | Random mating, unrelated, weakly structured | Livestock, inbred lines, family-based cohorts, diverse panels |
| **Computational cost** | O(np) prep + O(p²) per window via C++ | O(n²p) GRM + O(n³) Cholesky + O(np) whitening |
| **RAM** | Proportional to one window | n×n GRM + n×n whitening factor held per chromosome |
| **Markers** | Scales to 10 M+ | Practical to ~200 k per chromosome |
| **Block accuracy** | Slightly inflated in related populations | Correct for structured populations |
| **External dependencies** | None (Rcpp + RcppArmadillo always installed) | AGHmatrix, ASRgenomics (optional Suggests) |

**When r² produces wrong blocks.** In a livestock panel where sires appear many
times through progeny, r² between two SNPs in the same family cluster can reach
0.9 even when they are in different LD blocks. The whitened matrix removes this
familial signal. In practice, rV² blocks are typically 10–30% smaller and more
precisely delimited in related populations.

**When r² is preferable.** For large human biobank panels (n > 10,000,
p > 1,000,000) computing and inverting the GRM is computationally infeasible.
The standard r² inflates LD modestly but the effect on block boundaries is
small in approximately random-mating populations, and the 50× speed advantage
of the pure C++ path dominates.

```r
# Switching LD metric requires only a single argument change
blocks_r2  <- run_Big_LD_all_chr(be, method = "r2",  CLQcut = 0.70)
blocks_rv2 <- run_Big_LD_all_chr(be, method = "rV2", CLQcut = 0.70,
                                  kin_method = "chol")
```

---

## Haplotype analysis

Within each LD block, the SNPs co-segregate as a unit. The multiallelic
haplotype defined by the joint allele configuration at those SNPs is more
informative than any single SNP in the block. LDxBlocks provides three
functions for haplotype-level analyses downstream of block detection.

### Phase-free haplotype extraction

`extract_haplotypes()` concatenates each individual's allele codes (0, 1, or
2) for all SNPs within a block into a single character string:

```
Individual i, block b covering SNPs j1, j2, j3, j4:
  haplotype = paste0(g[i,j1], g[i,j2], g[i,j3], g[i,j4]) = "0120"
```

Individuals with any missing genotype in the block receive a haplotype string
containing the `na_char` placeholder (default `"."`); these are excluded from
frequency calculations in `compute_haplotype_diversity()`.

> **A note on phasing.** True gametic haplotypes require statistical phasing
> (SHAPEIT, Beagle) or read-backed phasing. The strings produced here are
> *diploid allele strings*, not gametic phases. Within a high-LD block these
> strings are nonetheless nearly 1:1 with true haplotype classes (Calus et al.
> 2008) and are fully sufficient for diversity analysis and genomic prediction.
> If gametic phases are required, phase first and convert the phased VCF to
> 0/1/2 dosages before reading with `read_geno()`.

### Haplotype diversity metrics

`compute_haplotype_diversity()` returns four metrics per block from the
haplotype string frequency distribution:

**Richness ($k$):** number of unique haplotype strings. A highly recombined
region will show many distinct haplotypes; a selective sweep will have one
dominant haplotype with very high frequency.

**Expected heterozygosity ($H_e$):** Nei's (1973) gene diversity, corrected
for sample size:

$$H_e = \frac{n}{n-1}\left(1 - \sum_i p_i^2\right)$$

where $p_i$ is the frequency of haplotype $i$ and $n$ is the number of
individuals with non-missing haplotypes. $H_e = 0$ for a monomorphic block;
$H_e \to 1$ for many equally frequent haplotypes.

**Shannon entropy ($H'$):**

$$H' = -\sum_i p_i \log_2 p_i$$

Sensitive to both richness and evenness. Measured in bits; equals $\log_2 k$
for equal-frequency haplotypes.

**Dominant haplotype frequency ($f_{\max}$):** frequency of the most common
haplotype string. Values near 1.0 indicate a selective sweep or a strong
founder effect in the region.

### Haplotype feature matrix for genomic prediction

`build_haplotype_feature_matrix()` converts haplotype strings to a numeric
dosage matrix suitable for GBLUP, BayesB, random forests, or any other genomic
prediction framework. For each block, the `top_n` most frequent haplotypes are
selected as reference haplotypes and each individual receives a dosage of 2
(homozygous match) or 0 (absent) for each reference.

The resulting matrix has dimension n_individuals × (n_blocks × top_n):

```r
# Build a haplotype-GRM for use in GBLUP
feat  <- build_haplotype_feature_matrix(haps, top_n = 5, scale_features = TRUE)
G_hap <- tcrossprod(feat) / ncol(feat)   # haplotype GRM
# ... feed G_hap to rrBLUP::kinship.BLUP or ASReml-R
```

**Why haplotype features outperform single-SNP models.** Single-SNP additive
models miss multi-locus epistatic effects and are sensitive to the choice of
which tag SNP represents a block. Haplotype dosages implicitly encode the joint
configuration of all SNPs in the block, capturing dominance and inter-SNP
interactions without explicit interaction terms. Calus et al. (2008) and de
Roos et al. (2009) demonstrated consistent accuracy improvements of 2–8% over
SNP-based models in livestock panels; the advantage is greatest for traits with
strong dominance or in populations where blocks are well-defined by long-range
LD.

---

## Parameter auto-tuning

When GWAS-significant markers are available, `tune_LD_params()` automatically
selects the `CLQcut` (and optionally other parameters) that best captures those
markers within LD blocks, minimising in priority order:

1. **Unassigned markers** — markers not falling within any block boundary.
2. **Forced assignments** — nearest-block placements, flagged with `*` in the
   output.
3. **Number of blocks** — parsimony; fewer blocks preferred among equal-score
   combinations.
4. **Deviation from target median block size** — biological plausibility via
   `target_bp_band` (default 50 kb – 500 kb).

When `prefer_perfect = TRUE` (default), combinations achieving zero unassigned
and zero forced assignments are isolated first and the winner is selected from
that subset by criteria 3 and 4.

```r
result <- tune_LD_params(
  geno_matrix    = my_geno,
  snp_info       = my_snp_info,
  gwas_df        = my_gwas,          # data.frame: Marker, CHR, POS
  grid           = NULL,             # use built-in 4-point CLQcut grid
  prefer_perfect = TRUE,
  target_bp_band = c(5e4, 5e5),
  parallel       = FALSE,            # TRUE + future::plan() for speed
  seed           = 42L
)

result$best_params       # named list of selected parameters
result$score_table       # all combinations and scores
result$perfect_table     # zero-zero combinations (if any found)
result$final_blocks      # block table from best_params, all chromosomes
result$gwas_assigned     # gwas_df with LD_block column added
```

The default grid searches `CLQcut` ∈ {0.65, 0.70, 0.75, 0.80} with all other
parameters held at their defaults. A custom grid can be supplied as a
`data.frame` where each row is one parameter combination.

For parallel execution across grid combinations:

```r
library(future)
plan(multisession, workers = 4L)
result <- tune_LD_params(my_geno, my_snp_info, my_gwas, parallel = TRUE)
```

---

## Scale strategies and backends

### The LDxBlocks_backend interface

All functions that require genotype data operate through the
`LDxBlocks_backend` S3 object rather than directly on a matrix. This single
abstraction makes all six input formats interchangeable with zero code change
in the algorithm:

```r
be    <- read_geno("mydata.bed")      # opens PLINK BED, reads .bim/.fam
chunk <- read_chunk(be, col_idx)      # returns n_samples x length(col_idx) matrix
close_backend(be)                     # releases memory-mapped file handle
```

`read_chunk()` dispatches to the appropriate low-level reader based on
`be$type`:

| Backend type | `read_chunk()` implementation |
|---|---|
| `"numeric"`, `"hapmap"`, `"vcf"`, `"matrix"` | Direct column slice of in-memory matrix |
| `"gds"` | `snpgdsGetGeno(snp.id=..., sample.id=...)` — one call, no filter cycle |
| `"bed"` | `BEDMatrix` row × column index (OS-level memory mapping) |

For GDS and BED backends, only the requested column slice is loaded per call.
The genome-wide matrix never exists in RAM simultaneously.

### Memory requirements by configuration

The table below gives approximate peak RAM for a 500-individual, 10 M-SNP
whole-genome dataset under different backends and methods.

| Configuration | Peak RAM | Notes |
|---|---|---|
| Plain matrix, r², all in RAM | ~40 GB | Impractical |
| GDS or BED backend, r², `subSegmSize = 1500` | ~300 MB | One chromosome window at a time |
| GDS or BED backend, r², `subSegmSize = 5000` | ~1 GB | Faster; larger windows |
| Any backend, rV², `method = "rV2"` | ~4 GB | GRM n×n + whitening factor always in RAM |

### Recommended configurations by dataset size

| Markers | Individuals | Recommended configuration |
|---|---|---|
| < 100 k | any | `method = "r2"`, `format = "matrix"` |
| 100 k – 500 k | < 5,000 | `method = "rV2"` for structured populations |
| 100 k – 2 M | any | `method = "r2"`, `format = "vcf"` or `"bed"` |
| 2 M – 10 M | any | `method = "r2"`, `format = "gds"`, `n_threads = 8+` |
| > 10 M | any | `method = "r2"`, `format = "gds"`, increase `subSegmSize` |

---

## Full pipeline walkthrough

The numbered steps below describe what `run_Big_LD_all_chr()` executes for
each chromosome in order.

| Step | Action | Key parameter(s) |
|------|--------|-----------------|
| 1 | Accept genotype backend or wrap plain matrix | `read_geno()`, `format =` |
| 2 | Extract per-chromosome genotype slice | `read_chunk(backend, chr_idx)` |
| 4 | MAF filter + monomorphic removal in C++ | `MAFcut`, `maf_filter_cpp()` |
| 5 | Centre (r²) or centre + whiten (rV²) | `method`, `kin_method`, `prepare_geno()` |
| 6 | C++ boundary scan — find weak-LD cut points | `leng`, `subSegmSize`, `boundary_scan_cpp()` |
| 7 | Divide chromosome into sub-segments | `subSegmSize` |
| 8 | Per sub-segment: compute r² or rV² matrix in C++ | `CLQcut`, `compute_r2_cpp()` |
| 9 | Build binary adjacency matrix in C++ | `CLQcut`, `build_adj_matrix_cpp()` |
| 10 | Find maximal cliques (igraph Bron-Kerbosch) | `checkLargest`, `CLQmode` |
| 11 | Greedy clique assignment → bin vector | `split`, `clstgap`, `CLQD()` |
| 12 | MWIS block construction | internal |
| 13 | Re-merge across forced cut-points | automatic |
| 14 | Merge overlapping blocks | automatic |
| 15 | Map indices to bp position and rsID | `SNPinfo` |
| 16 | Re-index over full SNP set including monomorphics | automatic |
| 17 | Optionally append rare SNPs | `appendrare` |

**Detailed example with all parameters:**

```r
blocks <- run_Big_LD_all_chr(
  # ── Genotype input ──────────────────────────────────────────────────────────
  geno_matrix  = be,                  # LDxBlocks_backend or plain matrix

  # ── LD metric ───────────────────────────────────────────────────────────────
  method       = "r2",                # "r2" (default) or "rV2"
  kin_method   = "chol",              # "chol" (default) or "eigen" for rV2

  # ── Clique detection ────────────────────────────────────────────────────────
  CLQcut       = 0.70,                # r2 threshold for graph edges
  CLQmode      = "Density",           # "Density" (default) or "Maximal"
  clstgap      = 40000L,              # max bp gap within clique (split=TRUE)
  split        = FALSE,               # split cliques at genomic gaps

  # ── Subsegmentation ─────────────────────────────────────────────────────────
  leng         = 200L,                # boundary-scan half-window (SNPs)
  subSegmSize  = 1500L,               # max SNPs per CLQD call

  # ── Filtering ───────────────────────────────────────────────────────────────
  MAFcut       = 0.05,                # minor allele frequency minimum
  appendrare   = FALSE,               # append rare SNPs after detection

  # ── Large window heuristics ─────────────────────────────────────────────────
  checkLargest = FALSE,               # dense-core pre-pass for >= 500 SNPs

  # ── Parallelism and precision ────────────────────────────────────────────────
  n_threads    = 8L,                  # OpenMP threads for C++ LD kernel
  digits       = -1L,                 # -1 = no rounding (default)

  # ── Chromosome minimum ──────────────────────────────────────────────────────
  min_snps_chr = 10L,                 # skip chromosomes with fewer SNPs

  # ── Reproducibility ─────────────────────────────────────────────────────────
  seed         = 42L,
  verbose      = TRUE
)
```

---

## Function reference

### Main pipeline

| Function | Description |
|----------|-------------|
| `run_Big_LD_all_chr()` | Chromosome-wise LD block detection. Accepts both plain matrices and `LDxBlocks_backend` objects. **Recommended entry point.** |
| `Big_LD()` | Core per-chromosome segmentation. Exposed for fine-grained control or single-chromosome analyses. |
| `CLQD()` | Clique detection within one sub-segment. Returns integer bin vector. |
| `tune_LD_params()` | Grid-search auto-tuner minimising unassigned GWAS marker placements. |

### I/O

| Function | Description |
|----------|-------------|
| `read_geno()` | Auto-dispatch genotype reader. Returns an `LDxBlocks_backend` object. |
| `read_chunk()` | Extract a genotype slice (n_samples × width) from any backend type. |
| `close_backend()` | Release file handles. No-op for in-memory backends. |

### LD computation

| Function | Description |
|----------|-------------|
| `compute_ld()` | Unified dispatcher: routes to `compute_r2_cpp()` or `compute_rV2_cpp()` based on `method =`. |
| `compute_r2()` | Standard r² matrix via C++ Armadillo + optional OpenMP. |
| `compute_rV2()` | rV² on a pre-whitened matrix (same C++ kernel as `compute_r2()`). |
| `prepare_geno()` | Centre (r²) or centre + whiten (rV²). Returns list of `adj_geno` and `V_inv_sqrt`. |
| `get_V_inv_sqrt()` | Whitening factor A such that AVA' = I (Cholesky or eigendecomposition). |

### C++ kernels (direct access)

| Function | Description |
|----------|-------------|
| `compute_r2_cpp()` | Full r² matrix. OpenMP outer loop. NA mean-imputed per column before computation. |
| `compute_rV2_cpp()` | rV² on pre-whitened matrix (identical kernel to `compute_r2_cpp()`). |
| `maf_filter_cpp()` | MAF + monomorphic filter in one O(np) C++ pass. Returns logical keep vector. |
| `build_adj_matrix_cpp()` | LD threshold → 0/1 integer adjacency matrix. |
| `col_r2_cpp()` | r² of one query column against all others. Used in boundary scan helper. |
| `compute_r2_sparse_cpp()` | Sparse r² for pairs within a bp distance window. Returns triplet (row, col, r²). |
| `boundary_scan_cpp()` | Cross-boundary LD scan. Returns 0/1 vector of valid cut positions. |

### Haplotype analysis

| Function | Description |
|----------|-------------|
| `extract_haplotypes()` | Phase-free diploid allele strings per block × individual. |
| `compute_haplotype_diversity()` | Per-block richness, He, Shannon entropy, dominant haplotype frequency. |
| `build_haplotype_feature_matrix()` | Haplotype dosage matrix (individuals × haplotype features) for genomic prediction. |

### Utilities

| Function | Description |
|----------|-------------|
| `summarise_blocks()` | Per-chromosome and genome-wide block size summary statistics. |
| `plot_ld_blocks()` | ggplot2 block diagram coloured by block size or chromosome. |

**Sample-alignment helpers** (not required for block detection; useful when
subsetting to phenotyped individuals before running the pipeline):

| Function | Description |
|----------|-------------|
| `read_pheno()` | Read a delimited phenotype/covariate file; auto-detect sample ID column; optionally align rows to a backend. |
| `align_geno_pheno()` | Intersect and reorder a backend and phenotype data frame to their common sample IDs. |

---

## Output objects

### `run_Big_LD_all_chr()` — block table

A `data.frame` with one row per detected LD block:

| Column | Type | Description |
|--------|------|-------------|
| `start` | integer | Index of the first SNP in the block (full SNP set, including monomorphics). |
| `end` | integer | Index of the last SNP. |
| `start.rsID` | character | SNP identifier at the block start. |
| `end.rsID` | character | SNP identifier at the block end. |
| `start.bp` | numeric | Base-pair position of the block start. |
| `end.bp` | numeric | Base-pair position of the block end. |
| `CHR` | character | Chromosome label (normalised, no `chr` prefix). |
| `length_bp` | integer | `end.bp - start.bp + 1`. |

### `tune_LD_params()` — named list

| Element | Type | Description |
|---------|------|-------------|
| `best_params` | named list | Selected parameter values. |
| `score_table` | data.frame | All grid combinations with: n_unassigned, n_forced, n_blocks, median_block_bp, penalty_bp. |
| `perfect_table` | data.frame or NULL | Combinations with n_unassigned = 0 and n_forced = 0. |
| `final_blocks` | data.table | Block table from `best_params`, all chromosomes. |
| `gwas_assigned` | data.frame | Input GWAS data with `LD_block` column added. Entries ending in `*` denote forced (nearest-block) assignments. |

### `extract_haplotypes()` — named list

| Element | Type | Description |
|---------|------|-------------|
| `block_<start>_<end>` | character vector | One haplotype string per individual (length n_samples). One element per block. |
| `attr(., "block_info")` | data.frame | Block metadata: block_id, CHR, start_bp, end_bp, n_snps. |

### `compute_haplotype_diversity()` — data.frame

| Column | Description |
|--------|-------------|
| `block_id` | Block name matching `names(haplotypes)`. |
| `n_ind` | Individuals with non-missing haplotypes. |
| `n_haplotypes` | Richness: number of unique haplotype strings. |
| `He` | Expected heterozygosity (Nei 1973), corrected for sample size. |
| `Shannon` | Shannon entropy in bits. |
| `freq_dominant` | Frequency of the most common haplotype. |

### `read_geno()` — LDxBlocks_backend

| Element | Type | Description |
|---------|------|-------------|
| `type` | character | Backend type: `"numeric"`, `"hapmap"`, `"vcf"`, `"gds"`, `"bed"`, or `"matrix"`. |
| `n_samples` | integer | Number of individuals. |
| `n_snps` | integer | Number of SNPs. |
| `sample_ids` | character | Individual identifiers. |
| `snp_info` | data.frame | SNP metadata: SNP, CHR, POS, REF, ALT. CHR is normalised (no `chr` prefix). |

---

## Memory and performance notes

### C++ core

The seven compiled functions in `src/ld_core.cpp` (428 lines,
RcppArmadillo + OpenMP) replace the most expensive R operations:

**`compute_r2_cpp()`** replaces `stats::cov()` + R arithmetic for the full r²
matrix. Speedup over pure R is approximately 40× for a 1,500 × 1,500 window
with 500 individuals. The outer loop is OpenMP-parallelised; all threads share
the pre-standardised matrix Z (zero-mean, unit-variance columns) and write to
disjoint (j, k) element pairs.

**`boundary_scan_cpp()`** replaces the triple-nested R loop in
`cutsequence.modi()`. For a chromosome with 50,000 SNPs and `leng = 200`, this
eliminates approximately 150,000 small R-level matrix operations that the
subsegmentation step would otherwise perform.

**`maf_filter_cpp()`** replaces `apply(G, 2, ...)` with a single compiled pass.
Approximately 10× faster for panels with > 100,000 SNPs; handles NA mean
imputation and monomorphic detection in the same loop.

**`build_adj_matrix_cpp()`** replaces `ifelse(LD >= cut, 1L, 0L)` with a C++
write in place, avoiding the allocation of the intermediate logical matrix.

### Never-full-genome memory model

LDxBlocks enforces a strict memory contract: **the full genotype matrix is never
held in RAM at once for any dataset size or format.**

**Numeric dosage CSV** — Two-pass chunked reading following the OptSLDP pattern
(Akohoue et al. 2026): Pass 1 scans the header and counts rows with zero data
loading. A single pre-allocated matrix is then filled in 50,000-row chunks via
successive `data.table::fread()` calls. Peak RAM = one chunk (not 2× the file).
`gc(FALSE)` is called after each chunk.

**VCF and HapMap** — Auto-converted to a SNPRelate GDS cache on first call
(placed next to the source file). Subsequent calls reuse the cache. All access
is streaming via `read_chunk()`.

**GDS and PLINK BED** — `read_chunk(backend, col_idx)` is called once per
sub-segment per chromosome. With `subSegmSize = 1500` and 50,000 SNPs per
chromosome, this is approximately 33 disk reads per chromosome. Each read loads
only a 1,500-column slice; the rest of the genome remains on disk.

**Chromosome loop** — In `run_Big_LD_all_chr()` and `extract_haplotypes()`,
each chromosome is extracted, processed, freed with `rm()`, and `gc(FALSE)` is
called before the next chromosome is touched, preventing heap fragmentation
from accumulating across 20–30 chromosome passes.

### OpenMP thread count

The `n_threads` parameter controls the number of OpenMP threads in
`compute_r2_cpp()`. Thread scaling is efficient up to approximately 8–16
threads for typical window sizes (1,500 SNPs). For very large windows
(`subSegmSize = 5000+`), higher thread counts remain beneficial:

```r
# Auto-detect physical cores
n_thr  <- parallel::detectCores(logical = FALSE)
blocks <- run_Big_LD_all_chr(be, n_threads = n_thr)
```

---

## Documentation

Full documentation, function reference, and tutorials are available at:

<https://FAkohoue.github.io/LDxBlocks/>

To read the vignette after installation:

```r
vignette("LDxBlocks-workflow", package = "LDxBlocks")
```

---

## Citation

If you use `LDxBlocks` in published research, please cite:

```
LDxBlocks Development Team (2025).
LDxBlocks: Genome-Wide LD Block Detection, Haplotype Analysis, and Genomic
Prediction Features with Kinship-Adjusted Correlations.
R package version 0.3.0.
https://github.com/FAkohoue/LDxBlocks
```

Please also cite the underlying methodological papers:

```
Kim S-A, Cho C-S, Kim S-R, Bull SB, Yoo Y-J (2018).
A new haplotype block detection method for dense genome sequencing data based
on interval graph modeling and dynamic programming.
Bioinformatics, 34(4), 588-596.
https://doi.org/10.1093/bioinformatics/btx609

VanRaden PM (2008).
Efficient methods to compute genomic predictions.
Journal of Dairy Science, 91(11), 4414-4423.
https://doi.org/10.3168/jds.2007-0980

Calus MPL, Meuwissen THE, de Roos APW, Veerkamp RF (2008).
Accuracy of genomic selection using different methods to define haplotypes.
Genetics, 178(1), 553-561.
https://doi.org/10.1534/genetics.107.080838

de Roos APW, Hayes BJ, Goddard ME (2009).
Reliability of genomic predictions across multiple populations.
Genetics, 183(4), 1545-1553.
https://doi.org/10.1534/genetics.109.104935

Nei M (1973).
Analysis of gene diversity in subdivided populations.
Proceedings of the National Academy of Sciences, 70(12), 3321-3323.
https://doi.org/10.1073/pnas.70.12.3321
```

---

## Contributing

Bug reports, feature requests, and pull requests are welcome:

<https://github.com/FAkohoue/LDxBlocks/issues>

Before opening a pull request, please:

1. Run `devtools::check()` with zero errors and zero warnings.
2. Add or update `tests/testthat/test-core.R` for all new functionality.
3. Rebuild documentation with `devtools::document()`.

---

## License

MIT + file LICENSE © Félicien Akohoue

---

## References

Kim S-A, Cho C-S, Kim S-R, Bull SB, Yoo Y-J (2018). A new haplotype block
detection method for dense genome sequencing data based on interval graph
modeling and dynamic programming. *Bioinformatics* **34**(4):588-596.
<https://doi.org/10.1093/bioinformatics/btx609>

Mangin B, Siberchicot A, Nicolas S, Doligez A, This P, Cierco-Ayrolles C (2012).
Novel measures of linkage disequilibrium that correct the bias due to population
structure and relatedness. *Heredity* **108**(3):285-291.
<https://doi.org/10.1038/hdy.2011.73>

VanRaden PM (2008). Efficient methods to compute genomic predictions. *Journal
of Dairy Science* **91**(11):4414–4423.
<https://doi.org/10.3168/jds.2007-0980>

Gabriel SB, Schaffner SF, Nguyen H, et al. (2002). The structure of haplotype
blocks in the human genome. *Science* **296**(5576):2225–2229.
<https://doi.org/10.1126/science.1069424>

Calus MPL, Meuwissen THE, de Roos APW, Veerkamp RF (2008). Accuracy of genomic
selection using different methods to define haplotypes. *Genetics*
**178**(1):553–561. <https://doi.org/10.1534/genetics.107.080838>

de Roos APW, Hayes BJ, Goddard ME (2009). Reliability of genomic predictions
across multiple populations. *Genetics* **183**(4):1545–1553.
<https://doi.org/10.1534/genetics.109.104935>

Nei M (1973). Analysis of gene diversity in subdivided populations.
*Proceedings of the National Academy of Sciences* **70**(12):3321–3323.
<https://doi.org/10.1073/pnas.70.12.3321>
