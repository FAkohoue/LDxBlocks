# LDxBlocks — Genome-Wide LD Block Detection, Haplotype Analysis, and Genomic Prediction Features

<p align="center">
  <img src="man/figures/logo.png" alt="LDxBlocks logo" width="190%">
</p>

<!-- badges: start -->
[![R-CMD-check](https://github.com/FAkohoue/LDxBlocks/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/FAkohoue/LDxBlocks/actions/workflows/R-CMD-check.yaml)
[![pkgdown](https://github.com/FAkohoue/LDxBlocks/actions/workflows/pkgdown.yaml/badge.svg)](https://github.com/FAkohoue/LDxBlocks/actions/workflows/pkgdown.yaml)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
<!-- badges: end -->

---

## Table of contents
1. [Motivation](#1-motivation)
2. [Summary](#2-summary)
3. [Relationship to the original Big-LD algorithm](#3-relationship-to-the-original-big-ld-algorithm)
   - 3.1. [Computational core: R loops replaced by C++](#31-computational-core-r-loops-replaced-by-c)
   - 3.2. [Memory model: never-full-genome](#32-memory-model-never-full-genome)
   - 3.3. [Kinship correction: rV²](#33-kinship-correction-rv)
   - 3.4. [Singleton SNP handling](#34-singleton-snp-handling)
   - 3.5. [Bug fix: zero-row assignment](#35-bug-fix-zero-row-assignment)
   - 3.6. [Downstream pipeline](#36-downstream-pipeline)
   - 3.7. [What is kept exactly](#37-what-is-kept-exactly)
4. [Installation](#4-installation)
5. [Quick start](#5-quick-start)
6. [Input formats](#6-input-formats)
   - 6.1. [Genotype input format](#61-genotype-input-format)
   - 6.2. [Phenotype input format](#62-phenotype-input-format)
     - 6.2.1. [Format 1 — Named numeric vector](#621-format-1-named-numeric-vector-single-trait-simplest)
     - 6.2.2. [Format 2 — Data frame, single trait](#622-format-2-data-frame-single-trait)
     - 6.2.3. [Format 3 — Data frame, multiple traits](#623-format-3-data-frame-multiple-traits)
     - 6.2.4. [Format 4 — Named list](#624-format-4-named-list-different-individuals-per-trait)
   - 6.3. [ID matching rules](#63-id-matching-rules)
   - 6.4. [Preparing BLUEs from raw phenotype data](#64-preparing-blues-from-raw-phenotype-data)
7. [Statistical background](#7-statistical-background)
   - 7.1. [MAF filtering](#71-maf-filtering)
   - 7.2. [Genotype preparation](#72-genotype-preparation)
   - 7.3. [Subsegmentation](#73-subsegmentation)
   - 7.4. [Clique detection (CLQD)](#74-clique-detection-clqd)
   - 7.5. [Block construction](#75-block-construction)
8. [LD metrics: r² versus rV²](#8-ld-metrics-r-versus-rv)
9. [Clique detection mode (CLQmode)](#9-clique-detection-mode-clqmode)
   - 9.1. [Mode 1 — Density](#91-mode-1-density-clqmode-density-default)
   - 9.2. [Mode 2 — Maximal](#92-mode-2-maximal-clqmode-maximal)
   - 9.3. [Mode 3 — Louvain](#93-mode-3-louvain-clqmode-louvain)
   - 9.4. [Mode 4 — Leiden](#94-mode-4-leiden-clqmode-leiden)
   - 9.5. [Summary comparison](#95-summary-comparison)
   - 9.6. [Recommended configurations](#96-recommended-configurations)
10. [Haplotype analysis](#10-haplotype-analysis)
    - 10.1. [Phase-free haplotype extraction](#101-phase-free-haplotype-extraction)
    - 10.2. [Haplotype diversity metrics](#102-haplotype-diversity-metrics)
    - 10.3. [Haplotype feature matrix for genomic prediction](#103-haplotype-feature-matrix-for-genomic-prediction)
    - 10.4. [Haplotype-based genomic prediction pipeline](#104-haplotype-based-genomic-prediction-pipeline)
    - 10.5. [Cross-validation and prediction accuracy](#105-cross-validation-and-prediction-accuracy)
    - 10.6. [Between-population comparison](#106-between-population-comparison)
    - 10.7. [Haplotype network visualisation](#107-haplotype-network-visualisation)
    - 10.8. [Multi-environment stability](#108-multi-environment-stability)
    - 10.9. [Candidate region export](#109-candidate-region-export)
    - 10.10. [Per-allele effect decomposition](#1010-per-allele-effect-decomposition)
    - 10.11. [Sliding-window diversity scan](#1011-sliding-window-diversity-scan)
    - 10.12. [True diplotype inference](#1012-true-diplotype-inference)
    - 10.13. [Rare-allele collapsing](#1013-rare-allele-collapsing)
    - 10.14. [Cross-panel harmonisation](#1014-cross-panel-harmonisation)
    - 10.15. [Haplotype association testing](#1015-haplotype-association-testing)
    - 10.16. [Breeding decision tools](#1016-breeding-decision-tools)
    - 10.17. [Cross-population effect concordance (haplotype)](#1017-cross-population-effect-concordance)
    - 10.18. [Cross-population effect concordance (external GWAS)](#1018-cross-population-effect-concordance-external-gwas)
    - 10.19. [Within-block and between-block epistasis detection](#1019-within-block-and-between-block-epistasis-detection)
11. [Parameter auto-tuning](#11-parameter-auto-tuning)
12. [Scale strategies and backends](#12-scale-strategies-and-backends)
    - 12.1. [The LDxBlocks_backend interface](#121-the-ldxblocksbackend-interface)
    - 12.2. [Memory requirements by configuration](#122-memory-requirements-by-configuration)
    - 12.3. [Recommended configurations by dataset size](#123-recommended-configurations-by-dataset-size)
13. [Full pipeline walkthrough](#13-full-pipeline-walkthrough)
14. [Function reference](#14-function-reference)
    - 14.1. [Main pipeline](#141-main-pipeline)
    - 14.2. [I/O](#142-io)
    - 14.3. [LD computation](#143-ld-computation)
    - 14.4. [C++ kernels (direct access)](#144-c-kernels-direct-access)
    - 14.5. [Haplotype analysis](#145-haplotype-analysis)
    - 14.6. [Analysis extensions](#146-analysis-extensions)
    - 14.7. [True haplotype inference and harmonisation](#147-true-haplotype-inference-and-harmonisation)
    - 14.8. [Haplotype association testing](#148-haplotype-association-testing)
    - 14.9. [Breeding decision tools](#149-breeding-decision-tools)
    - 14.10. [Utilities](#1410-utilities)
15. [Output objects](#15-output-objects)
    - 15.1. [`run_Big_LD_all_chr()` — block table](#151-runbigldallchr-block-table)
    - 15.2. [`run_ldx_pipeline()` — named list](#152-runldxpipeline-named-list)
    - 15.3. [`tune_LD_params()` — named list](#153-tuneldparams-named-list)
    - 15.4. [`extract_haplotypes()` — named list](#154-extracthaplotypes-named-list)
    - 15.5. [`compute_haplotype_diversity()` — data.frame](#155-computehaplotypediversity-dataframe)
    - 15.6. [`read_geno()` — LDxBlocks_backend](#156-readgeno-ldxblocksbackend)
    - 15.7. [`test_block_haplotypes()` — LDxBlocks_haplotype_assoc](#157-test_block_haplotypes--ldxblocks_haplotype_assoc)
    - 15.8. [`estimate_diplotype_effects()` — LDxBlocks_diplotype](#158-estimate_diplotype_effects--ldxblocks_diplotype)
    - 15.9. [`score_favorable_haplotypes()` — data frame](#159-score_favorable_haplotypes--data-frame)
    - 15.10. [`summarize_parent_haplotypes()` — data frame](#1510-summarize_parent_haplotypes--data-frame)
    - 15.11. [`compare_block_effects()` — LDxBlocks_effect_concordance](#1511-compare_block_effects--ldxblocks_effect_concordance)
    - 15.12. [`compare_gwas_effects()` — LDxBlocks_effect_concordance](#1512-compare_gwas_effects--ldxblocks_effect_concordance)
16. [Memory and performance notes](#16-memory-and-performance-notes)
    - 16.1. [C++ core](#161-c-core)
    - 16.2. [Never-full-genome memory model](#162-never-full-genome-memory-model)
    - 16.3. [OpenMP thread count](#163-openmp-thread-count)
17. [Documentation](#17-documentation)
18. [Citation](#18-citation)
19. [Contributing](#19-contributing)
20. [License](#20-license)
21. [References](#21-references)

---
---

## 1. Motivation

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

    rV²ᵢⱼ = Cov(Xᵛᵢ, Xᵛⱼ)² / [Var(Xᵛᵢ) · Var(Xᵛⱼ)]

where **X**ᵛ = **V**⁻¹⁄² **G̃** is the kinship-whitened, mean-centred
genotype matrix; **V** = **ZZ**ᵀ / (2 Σⱼ pⱼqⱼ) is the VanRaden
(2008) GRM; and **V**⁻¹⁄² is its inverse square root.

**The scale problem.** The original Big-LD implementation (Kim et al. 2018) contains
no compiled code. For modern whole-genome sequencing panels with 2-10 million markers
the inner loops are prohibitively slow, and loading the full genotype matrix
before detection is impossible on most workstations. LDxBlocks addresses this
with a C++/Armadillo computational core compiled via Rcpp, OpenMP-parallelised
LD computation, a unified multi-format I/O layer, and a strict never-full-genome
memory model.

**The pipeline gap.** The original Big-LD stops at block boundaries. LDxBlocks
adds a complete downstream pipeline: statistical phasing via Beagle 5.x with
bigmemory-backed caching, haplotype extraction, diversity metrics, genomic
prediction features, association testing, and breeding decision tools.

---

## 2. Summary

`LDxBlocks` is a complete pipeline for genome-wide LD block detection in
related or structured populations, with downstream haplotype analysis and
genomic prediction utilities. It extends the Big-LD algorithm of Kim et al.
(2018) with 15 core improvements:

1. **Dual LD metric** — standard r² (default) and kinship-adjusted rV²
2. **C++/Armadillo core** — eleven compiled functions handle all expensive operations
3. **OpenMP parallelism** — outer loop of `compute_r2_cpp()` parallelised
4. **Unified multi-format I/O** — numeric CSV, HapMap, VCF/VCF.gz, GDS, PLINK BED, R matrix
5. **Never-full-genome memory model** — genome never held in RAM at once for any format
6. **MAF filter in C++** — single O(np) pass with NA imputation
7. **C++ boundary scan** — replaces R inner loop in subsegmentation
8. **Sparse r² computation** — O(p) cost for large sub-segments
9. **Automatic parameter tuning** — `tune_LD_params()` grid search
10. **Haplotype reconstruction** — phase-free and Beagle-phased pathways
11. **Diversity metrics** — richness, He, Shannon entropy, dominant haplotype frequency
12. **Prediction feature matrix** — multi-locus dosage columns for GBLUP/BayesB/ML
13. **Polynomial community detection** — Louvain and Leiden for WGS panels
14. **Sparse LD computation** — `max_bp_distance` restricts pairs to a physical window
15. **Memory-mapped genotype store** — `read_geno_bigmemory()` with OS-page access

---

## 3. Relationship to the original Big-LD algorithm

### 3.1. Computational core: R loops replaced by C++

| Original R operation | LDxBlocks C++ function | Speedup |
|---|---|---|
| `cor(subgeno)` per CLQD call | `compute_r2_cpp()` + OpenMP | ~40× for 1,500-SNP window |
| `apply(Ogeno, 2, ...)` MAF filter | `maf_filter_cpp()` single pass | ~10× for 100k+ SNPs |
| `cor()` inside boundary-scan loop | `boundary_scan_cpp()` compiled | ~20× per chromosome |
| `r2Mat[r2Mat >= CLQcut^2] <- 1` | `build_adj_matrix_cpp()` | eliminates intermediate allocation |
| Single-column correlation | `col_r2_cpp()` | used in boundary scan helper |
| Sparse within-window r² | `compute_r2_sparse_cpp()` | avoids O(p²) for large segments |
| LD-informed overlap resolution | `resolve_overlap_cpp()` | 15,700× per-SNP reduction on chr1 |
| Haplotype string building | `build_hap_strings_cpp()` | ~20-50× per block |
| Block-to-SNP interval lookup | `block_snp_ranges_cpp()` | O(p + n_blocks) single sweep |
| Chromosome haplotype extraction | `extract_chr_haplotypes_cpp()` | strings + freq tabulation in one OpenMP pass |
| Call-rate filter + imputation | `impute_and_filter_cpp()` | single O(n×p) pass |

### 3.2. Memory model: never-full-genome

The full genotype matrix is never held in RAM at once for any format. Numeric
dosage CSV is read in pre-allocated 50,000-row chunks. VCF and HapMap
auto-convert to a streaming GDS cache. GDS and PLINK BED backends load only
the SNP window per CLQD call. `gc(FALSE)` is called after each chromosome.

### 3.3. Kinship correction: rV²

`method = "rV2"` replaces every pairwise correlation with the
kinship-whitened equivalent. The whitening factor **A** is computed once per
chromosome from the VanRaden (2008) GRM via `get_V_inv_sqrt()` (Cholesky or
eigendecomposition), then applied to the centred genotype matrix. In related
populations, rV² blocks are typically 10-30% smaller and more precisely
delimited than r² blocks.

### 3.4. Singleton SNP handling

`singleton_as_block = TRUE` collects SNPs that receive `NA` from `CLQD()` and
appends them to the block table as single-SNP entries (`start == end`). Default
`FALSE` preserves original behaviour.

### 3.5. Bug fix: zero-row assignment

When a sub-segment contains no valid cliques, `nowLDblocks` has zero rows.
LDxBlocks wraps the assignment with `if (nrow(nowLD) > 0L)`, preventing the
backwards-sequence error in the original code.

### 3.6. Downstream pipeline

| Capability | Original Big-LD / gpart | LDxBlocks |
|---|---|---|
| Block detection | Yes (core algorithm) | Yes (same algorithm + C++ + rV²) |
| Statistical phasing | No | `run_ldx_pipeline(phase = TRUE)` calls `phase_with_beagle()` internally. Phased hap1/hap2/dosage cached as bigmemory backends for fast restart. `read_phased_vcf()` reads pre-phased VCF output. |
| Haplotype extraction | No | `extract_haplotypes()` — phased and unphased, backend streaming |
| Diversity metrics | No | `compute_haplotype_diversity()` — He, Shannon, richness, f_max |
| Post-GWAS QTL mapping | No | `define_qtl_regions()` — pleiotropic block detection |
| Genomic prediction features | No | `build_haplotype_feature_matrix()` |
| Output writers | No | Numeric dosage matrix, nucleotide character matrix, diversity CSV |
| Parameter auto-tuning | No | `tune_LD_params()` — grid search against GWAS marker coverage |
| Multi-format I/O | PLINK, VCF (gpart) | Numeric CSV, HapMap, VCF, GDS, BED, R matrix via unified backend |
| WGS-scale streaming | Partial (gpart GDS) | Full never-full-genome model for all formats |

### 3.7. What is kept exactly

- The interval graph modelling of LD bins (cliques of strong pairwise LD SNPs)
- `CLQD()`: bin vector assignment via maximal clique enumeration and greedy density-priority selection
- `constructLDblock()`: maximum-weight independent set via dynamic programming
- `appendSGTs()`: rare-SNP appending logic
- `cutsequence.modi()`: boundary-scan logic and forced-split fall-back
- All `CLQmode = "Density"` and `CLQmode = "Maximal"` clique scoring
- The `clstgap` physical distance splitting within cliques
- Block table column format for drop-in compatibility with downstream tools

---

## 4. Installation

Install from GitHub with vignettes (recommended):

```r
install.packages("remotes")
remotes::install_github("FAkohoue/LDxBlocks",
  build_vignettes = TRUE,
  dependencies    = TRUE
)
```

**Required dependencies** (installed automatically):

```r
install.packages(c("Rcpp", "RcppArmadillo", "igraph", "data.table", "dplyr"))
```

**Optional dependencies:**

```r
# GDS backend — required for .gds files; recommended for panels > 2 M SNPs
BiocManager::install("SNPRelate")

# PLINK BED backend
install.packages("BEDMatrix")

# Kinship-adjusted rV² (method = "rV2")
install.packages(c("AGHmatrix", "ASRgenomics"))

# Parallel parameter tuning
install.packages("future.apply")

# LD block visualisation
install.packages("ggplot2")

# Bigmemory-backed genotype store
install.packages("bigmemory")
```

---

## 5. Quick start

```r
library(LDxBlocks)

# ── Option A: end-to-end pipeline ────────────────────────────────────────────
result <- run_ldx_pipeline(
  geno_source    = "mydata.vcf.gz",
  out_dir        = "ldx_results",
  out_blocks     = "ldx_results/blocks.csv",
  out_diversity  = "ldx_results/diversity.csv",
  out_hap_matrix = "ldx_results/hap_matrix.csv",
  CLQcut         = 0.70,
  n_threads      = 8L,
  verbose        = TRUE
)

result$blocks          # data.frame of LD blocks
result$diversity       # per-block diversity table
result$haplotypes      # named list of haplotype strings

# ── Option B: block detection only ───────────────────────────────────────────
be <- read_geno("mydata.vcf.gz")
blocks <- run_Big_LD_all_chr(be, method = "r2", CLQcut = 0.70, n_threads = 8L)
close_backend(be)

# ── Option C: with Beagle phasing ────────────────────────────────────────────
# Place beagle.jar in out_dir first, then:
result <- run_ldx_pipeline(
  geno_source        = "mydata.vcf.gz",
  out_dir            = "ldx_results",
  out_blocks         = "ldx_results/blocks.csv",
  out_hap_matrix     = "ldx_results/hap_matrix.csv",
  phase              = TRUE,
  beagle_jar         = "ldx_results/beagle.jar",
  beagle_threads     = 8L,
  beagle_java_mem_gb = 16L,
  beagle_seed        = 42L,
  CLQcut             = 0.70,
  n_threads          = 8L
)

result$phase_method   # "beagle"
result$phased_vcf     # path to Beagle-phased VCF.gz
```

---

## 6. Input formats

`read_geno()` accepts a path to a genotype file (or an in-memory R matrix).
The format is auto-detected from the file extension when `format` is not
supplied.

### 6.1. Genotype input format

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

**VCF** — standard VCF v4.2. Both phased (`0|1`) and unphased (`0/1`) GT
fields are accepted. Multi-allelic sites use the first ALT allele. Missing
calls (`./.`) become `NA`.

**Chromosome normalisation.** The leading prefix `chr`, `Chr`, or `CHR` is
stripped from chromosome labels at read time. Polyploid sub-genome labels
(`1A`, `2D`) are preserved verbatim.

### 6.2. Phenotype input format

The `blues` argument of `run_haplotype_prediction()` and
`prepare_gblup_inputs()` accepts four formats:

**Format 1 — Named numeric vector:**
```r
blues <- c(G001 = 4.21, G002 = 3.87, G003 = 5.14)
```

**Format 2 — Data frame, single trait:**
```r
blues <- read.csv("blues.csv")   # any column names
res <- run_haplotype_prediction(geno, snp_info, blocks,
                                 blues = blues, id_col = "id", blue_col = "YLD")
```

**Format 3 — Data frame, multiple traits:**
```r
res <- run_haplotype_prediction(geno, snp_info, blocks,
                                 blues     = blues_df,
                                 id_col    = "id",
                                 blue_cols = c("YLD", "DIS", "PHT"))
```

**Format 4 — Named list (different individuals per trait):**
```r
blues <- list(
  YLD = c(G001 = 4.21, G002 = 3.87),
  DIS = c(G001 = 0.32, G003 = 0.28)
)
```

### 6.3. ID matching rules

Genotype IDs in phenotype data must match `rownames(geno_matrix)` exactly
(case-sensitive). The function takes the intersection, issues a message for
individuals present in only one source, and errors if no common individuals
are found.

### 6.4. Preparing BLUEs from raw phenotype data

```r
# Example with lme4 (single environment)
library(lme4)
m   <- lmer(YLD ~ (1|id) + (1|rep), data = field_data)
blues_lme4 <- data.frame(id  = rownames(coef(m)$id),
                          YLD = coef(m)$id[, 1])
```

---

## 7. Statistical background

### 7.1. MAF filtering

    AF_i = (Σⱼ g_ij) / (2 n_i),    MAF_i = min(AF_i, 1 − AF_i)

SNPs with MAF_i < τ_maf (default 0.05) are removed. Both operations run in a
single O(np) C++ pass by `maf_filter_cpp()`.

### 7.2. Genotype preparation

**Standard r² path:** Pearson r² via column standardisation followed by
BLAS-level matrix multiplication in `compute_r2_cpp()`.

**Kinship-adjusted rV² path:** VanRaden (2008) GRM computed via
`AGHmatrix::Gmatrix()`, bent and conditioned via `ASRgenomics::G.tuneup()`,
then **A** = **R**⁻¹ (Cholesky) or **Q** Λ⁻¹⁄² **Q**ᵀ (eigen) applied to
mean-centred genotypes before the same `compute_r2_cpp()` kernel.

### 7.3–7.5. Subsegmentation, clique detection, and block construction

See the [full pipeline walkthrough](#13-full-pipeline-walkthrough) and the
original Kim et al. (2018) paper for algorithmic details.

---

## 8. LD metrics: r² versus rV²

| | r² | rV² |
|---|---|---|
| **Population type** | Random mating, unrelated | Livestock, inbred lines, family-based cohorts |
| **Computational cost** | O(np) prep + O(p²) per window | O(n²p) GRM + O(n³) Cholesky + O(np) whitening |
| **RAM** | Proportional to one window | n×n GRM + whitening factor held per chromosome |
| **Block accuracy** | Slightly inflated in related populations | Correct for structured populations |
| **External dependencies** | None | AGHmatrix, ASRgenomics |

```r
blocks_r2  <- run_Big_LD_all_chr(be, method = "r2",  CLQcut = 0.70)
blocks_rv2 <- run_Big_LD_all_chr(be, method = "rV2", CLQcut = 0.70, kin_method = "chol")
```

---

## 9. Clique detection mode (CLQmode)

| | Density | Maximal | Louvain | Leiden |
|---|---|---|---|---|
| **Algorithm** | Bron-Kerbosch + density score | Bron-Kerbosch + size score | Community detection | Community detection |
| **Complexity** | Exponential (worst case) | Exponential (worst case) | O(n log n) | O(n log n) |
| **Connected communities?** | Yes | Yes | Post-processed in LDxBlocks | **Yes (formal guarantee)** |
| **WGS feasible?** | Only with `max_bp_distance` | Only with `max_bp_distance` | Yes | **Yes** |
| **Recommended for** | Chip panels (< 100k SNPs/chr) | Sparse chip panels | Not recommended (use Leiden) | **WGS panels (2M+ SNPs)** |

```r
# WGS recommended configuration
blocks <- run_Big_LD_all_chr(be, CLQmode = "Leiden", CLQcut = 0.80,
                              max_bp_distance = 500000L,
                              subSegmSize = 500L, leng = 50L,
                              checkLargest = TRUE, n_threads = n_threads)
```

---

## 10. Haplotype analysis

### 10.1. Phase-free haplotype extraction

`extract_haplotypes()` concatenates each individual's allele codes (0, 1, or
2) for all SNPs within a block into a single character string:

```
Individual i, block b covering SNPs j1–j4:
  haplotype = paste0(g[i,j1], g[i,j2], g[i,j3], g[i,j4]) = "0120"
```

> **A note on phasing.** True gametic haplotypes require statistical phasing.
> The strings produced by `extract_haplotypes()` are *diploid allele strings*,
> not gametic phases. Within a high-LD block these strings are nearly 1:1 with
> true haplotype classes (Calus et al. 2008) and are sufficient for diversity
> analysis and genomic prediction.
>
> If gametic phases are required, set `phase = TRUE` in `run_ldx_pipeline()`,
> which calls Beagle 5.x after LD block detection and caches phased
> hap1/hap2/dosage as bigmemory backends for fast restart. Place `beagle.jar`
> in `out_dir` before running.


### 10.2. Haplotype diversity metrics

`compute_haplotype_diversity()` returns four metrics per block:

**Richness (k):** number of unique haplotype strings.

**Expected heterozygosity (Hₑ):**

    Hₑ = n/(n−1) · (1 − Σᵢ pᵢ²)

**Shannon entropy (H'):**

    H' = −Σᵢ pᵢ log₂(pᵢ)

**Dominant haplotype frequency (f_max):** frequency of the most common
haplotype. Values near 1.0 indicate a selective sweep or strong founder effect.

### 10.3. Haplotype feature matrix for genomic prediction

`build_haplotype_feature_matrix()` converts haplotype strings to a numeric
dosage matrix. For each block, the `top_n` most frequent haplotypes are
selected and each individual receives:

- **Phased data**: 0 (neither gamete), 1 (one gamete — heterozygous), or 2
  (both gametes — homozygous). True allele copy number.
- **Unphased data**: 0 (absent) or 1 (present). The value 2 is never produced
  because the two chromosomes cannot be distinguished — an individual
  homozygous for an allele and one heterozygous for it produce different dosage
  strings and are treated as distinct allele classes.

```r
feat  <- build_haplotype_feature_matrix(haps, scale_features = TRUE)$matrix
G_hap <- tcrossprod(feat) / ncol(feat)   # haplotype GRM
```

### 10.4. Haplotype-based genomic prediction pipeline

`run_haplotype_prediction()` produces:

1. VanRaden GRM from haplotype features
2. GEBV for all genotyped individuals via REML-based GBLUP
3. Per-SNP additive effects backsolved from GEBV (Tong et al. 2025)
4. Local haplotype GEBV per block per individual
5. Block importance table ranked by scaled Var(local GEBV)

```r
# Single trait
res <- run_haplotype_prediction(geno, snp_info, blocks,
                                 blues = my_blues, id_col = "id", blue_col = "YLD")

# Multiple traits
res_mt <- run_haplotype_prediction(geno, snp_info, blocks,
                                    blues     = my_blues_df,
                                    id_col    = "id",
                                    blue_cols = c("YLD", "DIS"))

# Integrate GWAS evidence
qtl      <- define_qtl_regions(gwas, blocks, snp_info, p_threshold = 5e-8)
priority <- integrate_gwas_haplotypes(qtl, res, diversity = div)
priority[priority$priority_score == 3, ]
```

### 10.5. Cross-validation and prediction accuracy

`cv_haplotype_prediction()` estimates prediction accuracy via k-fold
cross-validation. Returns a `LDxBlocks_cv` object with `pa_mean` (PA ± SD
per trait), `pa_summary` (per-fold), `k`, and `n_rep`.

### 10.6. Between-population comparison

`compare_haplotype_populations()` computes per-block Weir-Cockerham FST,
maximum allele frequency difference, and a chi-squared test (Monte Carlo,
B = 2000). `divergent` flag: FST > 0.1 AND p < 0.05.

### 10.7. Haplotype network visualisation

`plot_haplotype_network()` draws a minimum-spanning network using
`igraph::mst()`. Nodes sized by frequency; edges weighted by Hamming distance.
Optional group colouring.

### 10.8. Multi-environment stability

`run_haplotype_stability()` runs Finlay-Wilkinson (1963) regression of
per-block local GEBV contributions against the environmental index. Returns
`b` (stability coefficient), `b_se`, `R²`, `s²d`, and `stable` flag.

### 10.9. Candidate region export

`export_candidate_regions()` converts `define_qtl_regions()` output to
`"bed"` (UCSC-compatible, 0-based), `"csv"`, or `"biomart"` format.

### 10.10. Per-allele effect decomposition

`decompose_block_effects()` aggregates per-SNP additive effects into a
per-haplotype-allele effect table: effect = sum(SNP_effect × allele_dosage).

### 10.11. Sliding-window diversity scan

`scan_diversity_windows()` computes He, Shannon entropy, n_eff_alleles, and
dominant haplotype frequency in sliding windows across the genome,
independently of LD block boundaries.

### 10.12. True diplotype inference

`infer_block_haplotypes()` converts raw haplotype strings into a structured
per-individual, per-block diplotype table. For phased input (from
`run_ldx_pipeline(phase = TRUE)`), `phase_ambiguous` is always `FALSE`.

### 10.13. Rare-allele collapsing

`collapse_haplotypes()` merges alleles below `min_freq` using
`"rare_to_other"`, `"nearest"` (Hamming-based), or `"tree_based"` (UPGMA)
strategies. A `label_map` attribute records every original → collapsed mapping.

### 10.14. Cross-panel harmonisation

`harmonize_haplotypes()` anchors allele identity to a reference dictionary,
matching by exact string or nearest Hamming neighbour within `max_hamming`.
Novel alleles receive label `"<novel>"`. A `harmonization_report` attribute
reports match quality per block.

### 10.15. Haplotype association testing

`test_block_haplotypes()` uses a unified Q+K mixed linear model with
**simpleM multiple-testing correction** (Gao et al. 2008, 2010, 2011).
The model jointly corrects for population structure and polygenic kinship:

> y = μ + α · x_hap + Σ β_k · PC_k + g + ε,   g ~ MVN(0, σ²_g G)

Key parameters:

- `n_pcs = 0L` (default, EMMAX): GRM random effect absorbs all structure
- `n_pcs = k`: top-k GRM eigenvectors as additional fixed effects (Q+K model)
- `n_pcs = NULL`: auto-selected from the GRM scree-plot elbow
- `optimize_pcs = FALSE`: when `TRUE`, auto-selects n_pcs by fitting null
  models for k = 0..`optimize_pcs_max` and minimising the criterion set by
  `optimize_method`. More principled than the elbow heuristic for GWAS.
- `optimize_pcs_max = 10L`: maximum PCs evaluated when `optimize_pcs = TRUE`
- `optimize_method = c("bic_lambda", "bic", "lambda")`: criterion for PC
  selection (only used when `optimize_pcs = TRUE`):
  - `"bic"` — minimise BIC of the null REML model
  - `"lambda"` — minimise |λ_GC − 1| (genomic control calibration)
  - `"bic_lambda"` (**default, recommended**) — hybrid: |λ−1| + 0.01·scaled_BIC.
    Minimises inflation/deflation while BIC breaks ties toward fewer PCs
- `sig_metric`: which p-value drives the `significant` flag:
  - `"p_wald"` — raw Wald p-value (use with a pre-corrected threshold)
  - `"p_fdr"` — Benjamini-Hochberg FDR (recommended for discovery)
  - `"p_simplem"` — simpleM Bonferroni-style: min(p × Meff, 1)
  - `"p_simplem_sidak"` — simpleM Šidák-style: 1−(1−p)^Meff (**recommended** for correlated haplotype predictors)
- `meff_scope`: scope for estimating the effective number of independent tests (Meff):
  - `"chromosome"` (**recommended**) — separate Meff per chromosome
  - `"global"` — one genome-wide Meff
  - `"block"` — one Meff per LD block (omnibus tests get Meff = 1)
- `meff_percent_cut`: variance threshold for simpleM eigendecomposition (default `0.995`)
- `plot = TRUE`: saves **PDF** plots (not PNG). Three files per run:
  `manhattan_<trait>.pdf`, `qq_<trait>.pdf`, `pca_grm.pdf` (GRM PCA coloured
  by phenotype), `grm_scree.pdf` (eigenvalue scree with selected PC in red)

All four p-value columns (`p_wald`, `p_fdr`, `p_simplem`, `p_simplem_sidak`)
are **always present** in every output regardless of `sig_metric`.

`estimate_diplotype_effects()` also now supports the full correction set with
a new `sig_metric` parameter. All four p-value columns (`p_omnibus_adj`,
`p_omnibus_fdr`, `p_omnibus_simplem`, `p_omnibus_simplem_sidak`) are always
present in `$omnibus_tests` regardless of the chosen `sig_metric`.

```r
# FDR-based discovery (EMMAX, default)
assoc_fdr <- test_block_haplotypes(
  haps, blues = blues_vec, blocks = blocks,
  sig_metric = "p_fdr", verbose = FALSE
)

# simpleM Šidák with Q+K correction, chromosome-wise Meff
assoc_sm <- test_block_haplotypes(
  haps, blues = blues_vec, blocks = blocks,
  n_pcs            = 3L,
  sig_metric       = "p_simplem_sidak",
  meff_scope       = "chromosome",
  meff_percent_cut = 0.995,
  verbose          = FALSE
)
assoc_sm$meff$trait$allele$chromosome   # per-chromosome effective test counts
head(assoc_sm$block_tests[order(assoc_sm$block_tests$p_omnibus), ])
```

`estimate_diplotype_effects()` decomposes variation at each block into
additive (a) and dominance (d) components. Dominance ratio d/a classifies
gene action: 0 = additive, ±1 = complete dominance, |d/a| > 1 = overdominance.

### 10.16. Breeding decision tools

`score_favorable_haplotypes()` produces a genome-wide stacking index for each
individual (sum of allele_effect × dosage across blocks, normalised to [0, 1]).

`summarize_parent_haplotypes()` produces a long-format allele inventory for
candidate parents: which alleles they carry, at what dosage, and with what
population frequency. `is_rare = TRUE` flags alleles with frequency < 10%.

### 10.17. Cross-population effect concordance

`compare_block_effects()` takes two `test_block_haplotypes()` result objects
from independent populations and computes per-block statistics quantifying
how consistently haplotype allele effects replicate across panels.

**Key output columns (one row per block per trait):**

| Column | What it answers |
|--------|----------------|
| `n_shared_alleles` | Alleles tested in both populations |
| `effect_correlation` | Pearson r of per-allele effects across populations (NA when < 3 shared alleles) |
| `direction_agreement` | Fraction of shared alleles with the same effect sign |
| `directionally_concordant` | `direction_agreement >= direction_threshold` (default 0.75) |
| `meta_effect` / `meta_SE` / `meta_p` | IVW meta-analytic effect (inverse-variance weighted, same as two-sample MR) |
| `Q_stat` / `Q_p` | Cochran Q: significant Q means effect sizes differ between populations |
| `I2` | I² inconsistency (0–100%); > 50% = substantial heterogeneity |
| `replicated` | `TRUE` when directionally concordant AND Q_p > 0.05 AND enough shared alleles |
| `boundary_overlap_ratio` | **Automatically computed output** (not a user-set value): bp(intersection) / bp(union) of the two populations' block boundaries. Requires `blocks_pop1` / `blocks_pop2` to be supplied; `NA` otherwise |
| `boundary_warning` | `TRUE` when `boundary_overlap_ratio < boundary_overlap_warn` (the **input parameter**, default `0.80`) |
| `match_type` | How the block was matched: `"exact"` (same `block_id` string), `"position"` (matched by genomic overlap), `"pop1_only"` (no Pop2 block overlaps at `overlap_min`), or `NA` when no block tables were supplied |

```r
# Use Pop A's block boundaries for Pop B (maximises shared alleles)
haps_B_harm <- harmonize_haplotypes(
  extract_haplotypes(geno_B, snp_info, blocks_A),   # same block coords
  reference = haps_A
)
assoc_A <- test_block_haplotypes(haps_A, blues = blues_A, blocks = blocks_A,
                                  sig_metric = "p_simplem_sidak")
assoc_B <- test_block_haplotypes(haps_B_harm, blues = blues_B, blocks = blocks_A,
                                  sig_metric = "p_simplem_sidak")

# block_match = "id" (default): match by block_id string — fast, backward-compatible
# block_match = "position":      match by genomic interval overlap — handles different
#                                 LD block boundaries between populations
conc <- compare_block_effects(
  assoc_A, assoc_B,
  pop1_name             = "PopA",
  pop2_name             = "PopB",
  blocks_pop1           = blocks_A,   # required for boundary_overlap_ratio + position matching
  blocks_pop2           = blocks_B,   # supply Pop B's OWN block table (may differ from A)
  block_match           = "position", # recommended when blocks differ between populations
  overlap_min           = 0.50,       # min IoU for two blocks to be considered the same region
  direction_threshold   = 0.75,
  boundary_overlap_warn = 0.80
)
# $concordance$match_type: "exact" | "position" | "pop1_only"
# "position" rows: same QTL region, different boundary definitions
conc$concordance[conc$concordance$replicated, ]   # replicated blocks
head(conc$shared_alleles)                          # per-allele IVW detail
print(conc)
```

---

### 10.18. Cross-population effect concordance (external GWAS)

When GWAS was run externally (GAPIT, TASSEL, FarmCPU, PLINK, or any other
tool), `compare_gwas_effects()` compares block-level lead-SNP effects between
two populations. It accepts either raw GWAS data frames or the pre-mapped
output of `define_qtl_regions()`.

**Key differences from `compare_block_effects()`:**

| | `compare_block_effects()` | `compare_gwas_effects()` |
|---|---|---|
| Input | `test_block_haplotypes()` results | GWAS tables or `define_qtl_regions()` output |
| Unit per block | Multiple haplotype alleles | One lead SNP |
| `effect_correlation` | Pearson r across alleles | Always NA |
| `direction_agreement` | Fraction of alleles with same sign | 0 or 1 |
| `Q_stat` / `I2` | Cochran Q, df = n_alleles − 1 | Always NA (df = 0) |
| `replicated` | dir_concordant AND Q_p > 0.05 | dir_concordant AND meta_p ≤ 0.05 |

**SE derivation.** Many GWAS tools omit the standard error. When `SE` is absent,
`compare_gwas_effects()` derives it from the z-score:
`SE = |BETA| / |z|`, `z = Φ⁻¹(P/2)`. The `se_derived_pop1` / `se_derived_pop2`
output columns flag which population required this derivation.

```r
# Path 1: pre-mapped (recommended) — most auditable
qtl_A <- define_qtl_regions(gwas_A, blocks, snp_info, p_threshold = 5e-8)
qtl_B <- define_qtl_regions(gwas_B, blocks, snp_info, p_threshold = 5e-8)

conc <- compare_gwas_effects(
  qtl_pop1      = qtl_A,
  qtl_pop2      = qtl_B,
  blocks_pop1   = blocks_A,       # Pop A's block table
  blocks_pop2   = blocks_B,       # Pop B's block table (may differ from A)
  block_match   = "position",     # match by genomic overlap — handles different boundaries
  overlap_min   = 0.50,           # min IoU for a valid match
  pop1_name     = "PopA",
  pop2_name     = "PopB"
)

# Path 2: raw GWAS + blocks (convenience — calls define_qtl_regions internally)
conc2 <- compare_gwas_effects(
  gwas_pop1     = gwas_A,
  gwas_pop2     = gwas_B,
  blocks_pop1   = blocks_A,
  blocks_pop2   = blocks_B,
  snp_info_pop1 = snp_info,
  # snp_info_pop2 = NULL: reuses snp_info_pop1 when marker panels are shared
  pop1_name     = "PopA",
  pop2_name     = "PopB",
  block_match   = "position",  # match blocks by genomic interval overlap
  overlap_min   = 0.50,        # blocks with IoU < 0.5 become "pop1_only"
  p_threshold   = 5e-8,
  beta_col      = "BETA",      # column name — change to match your GWAS tool
  se_col        = "SE",        # NULL or absent: derived from z-score automatically
  p_col         = "P"
)

# Output uses the same LDxBlocks_effect_concordance class
# GWAS-specific extra columns in $concordance:
#   lead_snp_pop1, lead_snp_pop2  — which SNP tagged the block in each pop
#   lead_p_pop1, lead_p_pop2      — lead SNP p-values
#   se_derived_pop1/pop2          — TRUE when SE was derived from z-score
#   both_pleiotropic              — TRUE when block is pleiotropic in both pops
conc$concordance[conc$concordance$replicated, ]
print(conc)
```

> **Shared marker panel.** When both populations were genotyped on the same
> array or sequenced with the same reference, set `snp_info_pop2 = NULL`
> (default) to reuse `snp_info_pop1`. When populations have different marker
> sets, supply `snp_info_pop2` explicitly.

### 10.19. Within-block and between-block epistasis detection

LDxBlocks provides three functions for epistasis detection that extend the
haplotype association framework. All operate on GRM-corrected REML residuals
from the same null model as `test_block_haplotypes()`, ensuring
population-structure-corrected tests throughout.

**`scan_block_epistasis()`** tests all C(p,2) SNP pairs within each
significant block for pairwise interaction on GRM-corrected REML residuals.
The model for each pair is:

> y = μ + aᵢxᵢ + aⱼxⱼ + aaᵢⱼ(xᵢ × xⱼ) + ε

Restricting to significant blocks avoids the genome-wide explosion: for 15
significant blocks with ~200 SNPs each, the total number of tests is ~300,000
rather than ~4.4 billion. Multiple testing is corrected by Bonferroni and
simpleM Sidak within each block, where Meff is estimated from the eigenspectrum
of the pairwise interaction column matrix.

**`scan_block_by_block_epistasis()`** is a trans-haplotype epistasis scan
that tests each significant haplotype allele against every allele at all
other blocks. This is conceptually analogous to a trans-eQTL scan but for
phenotypic haplotype interactions — it detects genetic background dependence
where a resistance haplotype at one locus only functions in the presence of a
specific background at another locus. For 25 significant alleles × 17,943
total alleles, the scan involves ~450,000 tests corrected by Bonferroni.

**`fine_map_epistasis_block()`** fine-maps a single block by identifying the
specific interacting SNP pairs. For blocks with p ≤ 200 SNPs it runs an
exhaustive pairwise scan. For larger blocks it uses LASSO with pairwise
interaction terms (`glmnet::cv.glmnet()`, `lambda.1se`), which avoids the
multiple-testing burden while identifying the most influential pairs.

```r
# Within-block epistasis scan (significant blocks only)
epi_within <- scan_block_epistasis(
  assoc              = assoc,          # test_block_haplotypes() result
  geno_matrix        = res$geno_matrix,
  snp_info           = snp_info,
  blocks             = blocks,
  blues              = blues_list,
  haplotypes         = haps,
  trait              = "BL",
  sig_blocks         = NULL,           # NULL = use significant_omnibus blocks
  max_snps_per_block = 300L,
  sig_metric         = "p_simplem_sidak",
  sig_threshold      = 0.05
)
print(epi_within)
epi_within$results[epi_within$results$significant, ]
epi_within$scan_summary

# Between-block trans-haplotype epistasis scan
epi_between <- scan_block_by_block_epistasis(
  assoc         = assoc,
  haplotypes    = haps,
  blues         = blues_list,
  blocks        = blocks,
  trait         = "BL",
  sig_alleles   = NULL,    # NULL = significant alleles from assoc
  sig_threshold = 0.05
)
print(epi_between)
epi_between$results[epi_between$results$significant, ]

# Fine-map a single block (auto-dispatches: pairwise or LASSO)
fine <- fine_map_epistasis_block(
  block_id   = "block_12_1054210_1086071",
  geno_matrix = res$geno_matrix,
  snp_info    = snp_info,
  blocks      = blocks,
  y_resid     = my_reml_residuals,   # from .fit_null_reml() or null model
  method      = "auto",              # pairwise <= 200 SNPs, lasso otherwise
  sig_threshold = 0.05
)
head(fine)
```

---

## 11. Parameter auto-tuning

`tune_LD_params()` selects `CLQcut` (and optionally other parameters)
minimising in priority order: unassigned GWAS markers, forced assignments,
number of blocks, deviation from target median block size.

```r
result <- tune_LD_params(
  geno_matrix    = my_geno,
  snp_info       = my_snp_info,
  gwas_df        = my_gwas,
  prefer_perfect = TRUE,
  target_bp_band = c(5e4, 5e5),
  parallel       = FALSE,
  seed           = 42L
)
result$best_params
result$gwas_assigned
```

---

## 12. Scale strategies and backends

### 12.1. The LDxBlocks_backend interface

```r
be    <- read_geno("mydata.bed")       # opens PLINK BED
chunk <- read_chunk(be, col_idx)       # n_samples × length(col_idx) matrix
close_backend(be)                      # release file handle
```

### 12.2. Memory requirements by configuration

| Configuration | Peak RAM |
|---|---|
| Plain matrix, r², all in RAM | ~40 GB (500 ind, 10M SNPs) |
| GDS or BED backend, r², `subSegmSize = 1500` | ~300 MB |
| Any backend, rV², `method = "rV2"` | ~4 GB |
| bigmemory, r², `subSegmSize = 500` | ~0.8 MB/window |

### 12.3. Recommended configurations by dataset size

| Markers | Individuals | Recommended configuration |
|---|---|---|
| < 100 k | any | `method = "r2"`, `format = "matrix"` |
| 100 k – 500 k | < 5,000 | `method = "rV2"` for structured populations |
| 100 k – 2 M | any | `method = "r2"`, `format = "vcf"` or `"bed"` |
| 2 M – 10 M | any | `CLQmode = "Leiden"`, `max_bp_distance = 500000L`, `format = "gds"` |
| > 10 M | any | `CLQmode = "Leiden"`, `max_bp_distance = 500000L`, `read_geno_bigmemory()` |

---

## 13. Full pipeline walkthrough

| Step | Action | Key parameter(s) |
|------|--------|-----------------|
| 1 | Accept genotype backend or wrap plain matrix | `read_geno()`, `format =` |
| 2 | Extract per-chromosome genotype slice | `read_chunk(backend, chr_idx)` |
| 3 | MAF filter + monomorphic removal in C++ | `MAFcut`, `maf_filter_cpp()` |
| 4 | Centre (r²) or centre + whiten (rV²) | `method`, `kin_method`, `prepare_geno()` |
| 4b | **Optional: Beagle phasing** (`phase = TRUE`). Calls `phase_with_beagle()` on the original VCF after imputation/filtering. Reads the phased VCF via `read_phased_vcf()`, aligns samples and SNPs via CHR+POS+REF+ALT+SNP composite key, and caches hap1/hap2/dosage as bigmemory backends (`ldxblocks_bm_phased_hap1`, `_hap2`, `_dos`) in `bigmemory_path`. On restart the VCF is **not re-read** — backends reattach from disk via fingerprint match (VCF path + mtime + n_snps). Haplotype extraction then uses gametic strings (`g1\|g2`). | `phase`, `beagle_jar`, `beagle_threads`, `beagle_java_mem_gb`, `beagle_seed`, `beagle_ref_panel`, `beagle_map_file` |
| 5 | C++ boundary scan — find weak-LD cut points | `leng`, `subSegmSize`, `boundary_scan_cpp()` |
| 6 | Divide chromosome into sub-segments | `subSegmSize` |
| 7 | Per sub-segment: compute r² or rV² matrix in C++ | `CLQcut`, `compute_r2_cpp()` |
| 8 | Build binary adjacency matrix in C++ | `CLQcut`, `build_adj_matrix_cpp()` |
| 9 | Find communities/cliques | `CLQmode`, `checkLargest`, `max_bp_distance` |
| 10 | Greedy clique assignment → bin vector | `split`, `clstgap`, `CLQD()` |
| 11 | MWIS block construction | internal |
| 12 | Re-merge across forced cut-points | automatic |
| 13 | Merge overlapping blocks | automatic |
| 14 | Map indices to bp position and rsID | `SNPinfo` |
| 15 | Re-index over full SNP set including monomorphics | automatic |
| 16 | Optionally append rare SNPs | `appendrare` |

**Detailed example with all parameters:**

```r
blocks <- run_Big_LD_all_chr(
  # ── Genotype input ──────────────────────────────────────────────────────────
  geno_matrix  = be,

  # ── LD metric ───────────────────────────────────────────────────────────────
  method       = "r2",          # "r2" (default) or "rV2"
  kin_method   = "chol",        # "chol" (default) or "eigen" for rV2

  # ── Clique detection ────────────────────────────────────────────────────────
  CLQcut       = 0.70,
  CLQmode      = "Density",     # "Density"/"Maximal"/"Louvain"/"Leiden"
  clstgap      = 40000L,
  split        = FALSE,

  # ── Subsegmentation ─────────────────────────────────────────────────────────
  leng         = 200L,
  subSegmSize  = 1500L,

  # ── Filtering ───────────────────────────────────────────────────────────────
  MAFcut       = 0.05,
  appendrare   = FALSE,

  # ── Large window heuristics ─────────────────────────────────────────────────
  checkLargest = FALSE,

  # ── Parallelism ─────────────────────────────────────────────────────────────
  n_threads    = 8L,
  digits       = -1L,

  # ── Chromosome minimum ──────────────────────────────────────────────────────
  min_snps_chr = 10L,

  # ── Reproducibility ─────────────────────────────────────────────────────────
  seed         = 42L,
  verbose      = TRUE
)
```

**Full `run_ldx_pipeline()` example with phasing:**

```r
result <- run_ldx_pipeline(
  # ── Genotype input ──────────────────────────────────────────────────────────
  geno_source    = "mydata.vcf.gz",   # VCF required for phase = TRUE
  out_dir        = "ldx_results",

  # ── Output paths ────────────────────────────────────────────────────────────
  out_blocks     = "ldx_results/blocks.csv",
  out_diversity  = "ldx_results/diversity.csv",
  out_hap_matrix = "ldx_results/hap_matrix.csv",
  hap_format     = "numeric",

  # ── Phasing (optional — requires VCF input and beagle.jar in out_dir) ───────
  phase              = TRUE,
  beagle_jar         = "ldx_results/beagle.jar",
  beagle_threads     = 8L,
  beagle_java_mem_gb = 16L,          # -Xmx16g; NULL = JVM default
  beagle_seed        = 42L,          # integer seed for reproducibility
  beagle_ref_panel   = NULL,         # phased reference VCF (optional)
  beagle_map_file    = NULL,         # genetic map for improved accuracy

  # ── Filtering & imputation ───────────────────────────────────────────────────
  maf_cut        = 0.05,
  impute         = "mean_rounded",

  # ── LD block detection ───────────────────────────────────────────────────────
  CLQcut         = 0.70,
  method         = "r2",
  CLQmode        = "Leiden",
  subSegmSize    = 500L,
  leng           = 50L,
  max_bp_distance = 500000L,
  n_threads      = 8L,

  # ── Haplotype extraction ─────────────────────────────────────────────────────
  min_snps_block = 3L,
  top_n          = 5L,

  # ── Bigmemory caching (restart-safe) ─────────────────────────────────────────
  use_bigmemory  = TRUE,
  bigmemory_path = "ldx_results/bm_cache",
  bigmemory_type = "char",

  # ── General ──────────────────────────────────────────────────────────────────
  verbose = TRUE
)
```

---

## 14. Function reference

### 14.1. Main pipeline

| Function | Description |
|----------|-------------|
| `run_ldx_pipeline()` | **Recommended end-to-end pipeline.** One call from genotype source to haplotype matrix, diversity table, and block CSV. Accepts `phase = TRUE` for Beagle phasing (place `beagle.jar` in `out_dir`). Phased data cached as bigmemory backends for fast restart. |
| `run_Big_LD_all_chr()` | Chromosome-wise LD block detection. Accepts both plain matrices and `LDxBlocks_backend` objects. |
| `Big_LD()` | Core per-chromosome segmentation. Called internally by `run_Big_LD_all_chr()`. Not exported. |
| `tune_LD_params()` | Grid-search auto-tuner minimising unassigned GWAS marker placements. |

### 14.2. I/O

| Function | Description |
|----------|-------------|
| `run_ldx_pipeline()` | See Main pipeline. Also handles all I/O for the end-to-end run. |
| `phase_with_beagle()` | Statistical phasing via Beagle 5.x. Called internally by `run_ldx_pipeline(phase = TRUE)`. Exposes all Beagle parameters: `java_mem_gb`, `map_file`, `chrom`, `seed`, `burnin`, `iterations`, `window`, `overlap`. Log written to `out_prefix.log`. Requires `beagle.jar`. |
| `read_phased_vcf()` | Read a phased VCF (pipe-separated `0\|1` GT fields) into `list(hap1, hap2, dosage, snp_info, sample_ids, phased=TRUE)`. Called internally after Beagle phasing. |
| `read_geno()` | Auto-dispatch genotype reader. Returns an `LDxBlocks_backend` object. |
| `read_chunk()` | Extract a genotype slice (n_samples × width) from any backend type. |
| `close_backend()` | Release file handles. No-op for in-memory backends. |
| `read_geno_bigmemory()` | Create a file-backed memory-mapped store from any source. |

### 14.3. LD computation

| Function | Description |
|----------|-------------|
| `compute_r2()` | Standard r² matrix via C++ Armadillo + optional OpenMP. |
| `compute_rV2()` | rV² on a pre-whitened matrix. |
| `prepare_geno()` | Centre (r²) or centre + whiten (rV²). |
| `get_V_inv_sqrt()` | Whitening factor **A** such that AVA' = I. |

### 14.4. C++ kernels (direct access)

| Function | Description |
|----------|-------------|
| `compute_r2_cpp()` | Full r² matrix. OpenMP outer loop. |
| `maf_filter_cpp()` | MAF + monomorphic filter in one O(np) C++ pass. |
| `build_adj_matrix_cpp()` | LD threshold → 0/1 integer adjacency matrix. |
| `col_r2_cpp()` | r² of one query column against all others. |
| `compute_r2_sparse_cpp()` | Sparse r² for pairs within a bp distance window. |
| `boundary_scan_cpp()` | Cross-boundary LD scan. Returns 0/1 vector of valid cut positions. |
| `resolve_overlap_cpp()` | BLAS DGEMM scoring of overlapping block boundaries. |

### 14.5. Haplotype analysis

| Function | Description |
|----------|-------------|
| `extract_haplotypes()` | Phase-free diploid allele strings per block × individual. For phased input (from `run_ldx_pipeline(phase = TRUE)`), produces gametic strings `g1\|g2`. |
| `compute_haplotype_diversity()` | Per-block richness, He, n_eff_alleles, Shannon, sweep_flag. |
| `build_haplotype_feature_matrix()` | Haplotype dosage matrix for genomic prediction. Phased: 0/1/2. Unphased: 0/1. |
| `compute_haplotype_grm()` | VanRaden GRM from haplotype feature matrix. |
| `decode_haplotype_strings()` | Decode dosage strings to nucleotide sequences. |
| `write_haplotype_numeric()` | Write haplotype dosage matrix to file. |
| `write_haplotype_character()` | Write nucleotide character matrix to file. |
| `write_haplotype_diversity()` | Write diversity table to CSV. |
| `define_qtl_regions()` | Map GWAS hits to LD blocks; detect pleiotropic blocks. |
| `backsolve_snp_effects()` | Derive per-SNP effects from GEBV (Tong et al. 2025). |
| `compute_local_gebv()` | Local haplotype GEBV per block per individual. |
| `prepare_gblup_inputs()` | Align phenotype data and haplotype GRM for external GBLUP solvers. |
| `run_haplotype_prediction()` | Single or multi-trait Tong et al. (2025) haplotype stacking pipeline. |
| `integrate_gwas_haplotypes()` | Combine GWAS, variance, and diversity evidence per block. |
| `rank_haplotype_blocks()` | Unified block ranking across 3 use cases. |

### 14.6. Analysis extensions

| Function | Description |
|----------|-------------|
| `cv_haplotype_prediction()` | K-fold cross-validation for the haplotype GBLUP model. |
| `compare_haplotype_populations()` | Per-block Weir-Cockerham FST between two sample groups. |
| `plot_haplotype_network()` | Minimum-spanning network of haplotype alleles for one block. |
| `run_haplotype_stability()` | Finlay-Wilkinson stability regression across environments. |
| `export_candidate_regions()` | Convert QTL regions to BED, CSV, or biomaRt format. |
| `decompose_block_effects()` | Per-haplotype-allele effect table from per-SNP additive effects. |
| `scan_diversity_windows()` | Sliding-window diversity scan across the genome. |

### 14.7. True haplotype inference and harmonisation

| Function | Description |
|----------|-------------|
| `infer_block_haplotypes()` | Structured diplotype table from raw haplotype strings. Phased input: `phase_ambiguous = FALSE`. |
| `collapse_haplotypes()` | Merge rare alleles via `"rare_to_other"`, `"nearest"`, or `"tree_based"`. |
| `harmonize_haplotypes()` | Cross-panel allele label harmonisation using reference dictionary. |

### 14.8. Haplotype association testing

| Function | Description |
|----------|-------------|
| `test_block_haplotypes()` | Block-level association tests via Q+K mixed model with simpleM multiple-testing correction (Gao et al. 2008, 2010, 2011). Parameters `sig_metric`, `meff_scope`, `meff_percent_cut`, `meff_max_cols` control correction method and Meff estimation scope. All four p-value flavours always present in output. Returns `LDxBlocks_haplotype_assoc`. |
| `estimate_diplotype_effects()` | Additive (a) and dominance (d) effects per block. Returns `LDxBlocks_diplotype`. |
| `compare_block_effects()` | Cross-population haplotype effect concordance. Computes IVW meta-analytic effects, Cochran Q heterogeneity, I² inconsistency, direction agreement, and replication flag per block. New `block_match` parameter: `"id"` (default, backward-compatible) matches by `block_id` string; `"position"` matches by genomic interval overlap (IoU ≥ `overlap_min`, default 0.50) — handles different LD block boundaries between populations. `boundary_overlap_ratio` is automatically computed (not user-set); `boundary_overlap_warn` (default 0.80) controls `boundary_warning`. Output `match_type` column: `"exact"` / `"position"` / `"pop1_only"`. Returns `LDxBlocks_effect_concordance`. |
| `compare_gwas_effects()` | Cross-population concordance from **external GWAS results** (GAPIT, TASSEL, FarmCPU, PLINK, etc.). Accepts raw GWAS data frames or pre-mapped `define_qtl_regions()` output. Derives SE from z-score when absent. `block_match = "position"` handles different block boundaries between populations (same as `compare_block_effects()`). One lead SNP per block; `effect_correlation` and Cochran Q are NA. Same output class as `compare_block_effects()`. |
| `scan_block_epistasis()` | Within-block pairwise SNP epistasis scan. Tests all C(p,2) SNP pairs within significant blocks on GRM-corrected REML residuals. Model: y = mu + a_i*x_i + a_j*x_j + aa_ij*(x_i*x_j) + e. Corrected via Bonferroni and simpleM Sidak within each block (Meff from interaction column eigenspectrum). `sig_metric` controls which correction drives the `significant` flag. Returns `LDxBlocks_epistasis`. |
| `scan_block_by_block_epistasis()` | Trans-haplotype between-block epistasis scan. Tests significant haplotype alleles against all other block alleles: O(n_sig x n_total_alleles) tests corrected by Bonferroni. Identifies genetic background dependence that single-block analyses cannot detect. Returns `LDxBlocks_block_epistasis`. |
| `fine_map_epistasis_block()` | Single-block epistasis fine-mapping. Dispatches to exhaustive pairwise scan (p <= 200 SNPs) or LASSO interaction search via `glmnet` (p > 200 SNPs). Requires pre-computed REML residuals (`y_resid`). |

### 14.9. Breeding decision tools

| Function | Description |
|----------|-------------|
| `score_favorable_haplotypes()` | Genome-wide haplotype stacking index per individual. |
| `summarize_parent_haplotypes()` | Long-format allele inventory for candidate parents. |

### 14.10. Utilities

| Function | Description |
|----------|-------------|
| `summarise_blocks()` | Per-chromosome and genome-wide block size summary statistics. |
| `plot_ld_blocks()` | ggplot2 block diagram coloured by block size or chromosome. |

---

## 15. Output objects

### 15.1. `run_Big_LD_all_chr()` — block table

| Column | Type | Description |
|--------|------|-------------|
| `start` | integer | Index of the first SNP in the block (full SNP set). |
| `end` | integer | Index of the last SNP. |
| `start.rsID` | character | SNP identifier at the block start. |
| `end.rsID` | character | SNP identifier at the block end. |
| `start.bp` | numeric | Base-pair position of the block start. |
| `end.bp` | numeric | Base-pair position of the block end. |
| `CHR` | character | Chromosome label (normalised, no `chr` prefix). |
| `length_bp` | integer | `end.bp - start.bp + 1`. |

### 15.2. `run_ldx_pipeline()` — named list

| Element | Type | Description |
|---------|------|-------------|
| `blocks` | data.frame | LD block table (same columns as Section 15.1). |
| `diversity` | data.frame | Per-block diversity metrics (same columns as Section 15.4). |
| `hap_matrix` | matrix | Haplotype feature matrix (individuals × haplotype columns). |
| `hap_matrix_info` | data.frame | Column metadata for `hap_matrix`. |
| `haplotypes` | named list | Raw haplotype strings per block (same structure as Section 15.3). |
| `geno_matrix` | matrix | Cleaned, imputed genotype matrix (individuals × SNPs). |
| `snp_info_filtered` | data.frame | SNP metadata after MAF filtering. |
| `phased_vcf` | character | Path to Beagle-phased VCF.gz in `out_dir`. `NULL` when `phase = FALSE`. |
| `phased_backend_desc` | character | Path to `ldxblocks_bm_phased_dos.desc` for reattachment. `NULL` when `use_bigmemory = FALSE` or `phase = FALSE`. |
| `phase_method` | character | `"beagle"` when Beagle was run; `"unphased"` otherwise. |
| `n_blocks` | integer | Total number of LD blocks detected. |
| `n_hap_columns` | integer | Number of columns in `hap_matrix`. |

**Bigmemory phased backends** (written to `bigmemory_path` when `use_bigmemory = TRUE` and `phase = TRUE`):

| File | Content | Type |
|------|---------|------|
| `ldxblocks_bm_phased_hap1.bin/.desc` | Gamete 1 alleles (individuals × SNPs) | `char`, 0/1 |
| `ldxblocks_bm_phased_hap2.bin/.desc` | Gamete 2 alleles (individuals × SNPs) | `char`, 0/1 |
| `ldxblocks_bm_phased_dos.bin/.desc` | Dosage = hap1 + hap2 (individuals × SNPs) | `char`, 0/1/2 |
| `ldxblocks_bm_phased_snpinfo.rds` | Aligned SNP metadata | — |
| `ldxblocks_bm_phased_sampleids.rds` | Sample IDs in aligned order | — |
| `ldxblocks_bm_phased_params.rds` | Fingerprint for cache validation | — |

On subsequent runs with the same `bigmemory_path` and unchanged phased VCF,
backends are reattached from disk without re-reading the VCF.

### 15.3. `tune_LD_params()` — named list

| Element | Type | Description |
|---------|------|-------------|
| `best_params` | named list | Selected parameter values. |
| `score_table` | data.frame | All grid combinations and scores. |
| `perfect_table` | data.frame or NULL | Combinations with n_unassigned = 0 and n_forced = 0. |
| `final_blocks` | data.table | Block table from `best_params`, all chromosomes. |
| `gwas_assigned` | data.frame | Input GWAS data with `LD_block` column added. Entries ending in `*` denote forced (nearest-block) assignments. |

### 15.4. `extract_haplotypes()` — named list

| Element | Type | Description |
|---------|------|-------------|
| `block_<start>_<end>` | character vector | One haplotype string per individual. |
| `attr(., "block_info")` | data.frame | block_id, CHR, start_bp, end_bp, n_snps. |

### 15.5. `compute_haplotype_diversity()` — data.frame

| Column | Description |
|--------|-------------|
| `block_id` | Block name matching `names(haplotypes)`. |
| `CHR` | Chromosome (normalised). |
| `start_bp`, `end_bp` | Block coordinates in base pairs. |
| `n_snps` | Number of SNPs in the block. |
| `n_ind` | Individuals with non-missing haplotypes. |
| `n_haplotypes` | Richness: number of unique haplotype strings. |
| `He` | Expected heterozygosity (Nei 1973), sample-size corrected. |
| `Shannon` | Shannon entropy in bits. |
| `n_eff_alleles` | Effective number of alleles = 1/Σpᵢ². |
| `freq_dominant` | Frequency of the most common haplotype. |
| `sweep_flag` | TRUE when freq_dominant ≥ 0.90. |
| `phased` | Logical: was phased input used? |

### 15.6. `read_geno()` — LDxBlocks_backend

| Element | Type | Description |
|---------|------|-------------|
| `type` | character | `"numeric"`, `"hapmap"`, `"vcf"`, `"gds"`, `"bed"`, or `"matrix"`. |
| `n_samples` | integer | Number of individuals. |
| `n_snps` | integer | Number of SNPs. |
| `sample_ids` | character | Individual identifiers. |
| `snp_info` | data.frame | SNP metadata: SNP, CHR, POS, REF, ALT. |

### 15.7. `test_block_haplotypes()` — LDxBlocks_haplotype_assoc

**`$allele_tests`** — one row per allele per block per trait:

| Column | Type | Description |
|--------|------|-------------|
| `block_id`, `CHR`, `start_bp`, `end_bp`, `trait` | — | Identifiers |
| `allele` | character | Haplotype allele string identifier |
| `allele_freq_tested` | numeric [0,1] | Allele frequency among tested individuals |
| `effect` | numeric | Additive effect on de-regressed scale |
| `SE` | numeric | Standard error |
| `t_stat` | numeric | t-statistic |
| `p_wald` | numeric | Two-sided Wald p-value (raw) |
| `p_fdr` | numeric | Benjamini-Hochberg FDR-adjusted p-value |
| `Meff` | numeric | Effective number of independent tests (simpleM, scope = `meff_scope`) |
| `alpha_simplem` | numeric | simpleM Bonferroni-style significance threshold: α / Meff |
| `alpha_simplem_sidak` | numeric | simpleM Šidák-style threshold: 1 − (1−α)^(1/Meff) |
| `p_simplem` | numeric | simpleM Bonferroni-style adjusted p-value: min(p × Meff, 1) |
| `p_simplem_sidak` | numeric | simpleM Šidák-style adjusted p-value: 1 − (1−p)^Meff |
| `significant` | logical | TRUE when the p-value chosen by `sig_metric` ≤ `sig_threshold` |

All p-value columns are always present regardless of `sig_metric`. The `significant`
flag is the only column that changes based on `sig_metric`.

**`$block_tests`** — one row per block per trait:

| Column | Type | Description |
|--------|------|-------------|
| `n_alleles_tested` | integer | Alleles above `min_freq` tested jointly |
| `F_stat` | numeric | Omnibus F-statistic |
| `df_LRT` | integer | Numerator df = n_alleles_tested |
| `p_omnibus` | numeric | Raw omnibus p-value |
| `p_omnibus_fdr` | numeric | BH FDR-adjusted omnibus p-value |
| `p_omnibus_adj` | numeric | Bonferroni-style per-trait adjustment (backward-compat) |
| `var_explained` | numeric [0,1] | Proportion of de-regressed variance explained |
| `Meff` | numeric | Block-level effective test count (from block-summary PC1 matrix) |
| `alpha_simplem` | numeric | Block-level simpleM Bonferroni threshold |
| `alpha_simplem_sidak` | numeric | Block-level simpleM Šidák threshold |
| `p_omnibus_simplem` | numeric | simpleM Bonferroni-style adjusted omnibus p-value |
| `p_omnibus_simplem_sidak` | numeric | simpleM Šidák-style adjusted omnibus p-value |
| `significant_omnibus` | logical | TRUE when the p-value chosen by `sig_metric` ≤ `sig_threshold` |

**Additional return list elements:**

| Element | Description |
|---------|-------------|
| `sig_metric` | Which p-value was used for `significant` / `significant_omnibus` |
| `meff_scope` | Scope used for Meff estimation (`"chromosome"`, `"global"`, or `"block"`) |
| `meff_percent_cut` | Variance threshold used for simpleM eigendecomposition |
| `meff` | Named list of Meff summaries per trait: `$allele$global`, `$allele$chromosome`, `$allele$block`, `$block$global`, `$block$chromosome` |
| `pc_model_selection` | `data.frame` with one row per k tested when `optimize_pcs = TRUE`: `n_pcs`, `BIC`, `lambda_gc`, `score`, `selected`. `NULL` when `optimize_pcs = FALSE` |

**Plots produced when `plot = TRUE`** (all saved as PDF):

| File | Description |
|------|-------------|
| `manhattan_<trait>.pdf` | Manhattan plot — one per trait |
| `qq_<trait>.pdf` | QQ plot — one per trait |
| `pca_grm.pdf` | PCA of GRM eigenvectors in PC1 × PC2, coloured by first trait phenotype |
| `grm_scree.pdf` | Scree plot of GRM eigenvalues (up to PC30); red bar and vertical line at selected k |

### 15.8. `estimate_diplotype_effects()` — LDxBlocks_diplotype

**`$omnibus_tests`** — one row per LD block per trait (updated columns):

| Column | Type | Description |
|--------|------|-------------|
| `block_id`, `trait` | — | Identifiers |
| `n_diplotypes` | integer | Diplotype classes with ≥ `min_n_diplotype` individuals |
| `F_stat` | numeric | F-statistic (one-way ANOVA on GRM-corrected residuals) |
| `df1`, `df2` | integer | Numerator and denominator degrees of freedom |
| `p_omnibus` | numeric | Raw p-value |
| `p_omnibus_adj` | numeric | Plain Bonferroni × n_blocks_per_trait (backward-compat) |
| `p_omnibus_fdr` | numeric | BH-FDR per trait |
| `p_omnibus_simplem` | numeric | simpleM Bonferroni-style: min(p × Meff, 1) |
| `p_omnibus_simplem_sidak` | numeric | simpleM Šidák-style: 1 − (1−p)^Meff (**recommended**) |
| `Meff` | numeric | Effective number of blocks (block-summary PC1 eigenspectrum) |
| `significant` | logical | TRUE when the p-value chosen by `sig_metric` < `sig_threshold` |

All five p-value columns are always present regardless of `sig_metric`.

**`$dominance_table`** — one row per allele pair per block per trait:

| Column | Type | Description |
|--------|------|-------------|
| `allele_A`, `allele_B` | character | The two alleles (A ≤ B alphabetically) |
| `mean_AA`, `mean_AB`, `mean_BB` | numeric | Diplotype class means |
| `a` | numeric | Additive effect: (mean_BB − mean_AA) / 2 |
| `d` | numeric | Dominance deviation: mean_AB − midpoint |
| `d_over_a` | numeric or NA | Dominance ratio. 0 = additive; ±1 = complete dominance; \|d/a\| > 1 = overdominance |
| `overdominance` | logical | TRUE when \|d/a\| > 1 |

### 15.9. `score_favorable_haplotypes()` — data frame

| Column | Type | Description |
|--------|------|-------------|
| `id` | character | Individual identifier |
| `stacking_index` | numeric [0,1] | Genome-wide normalised score |
| `n_blocks_scored` | integer | Blocks with at least one matched allele effect |
| `mean_block_score` | numeric | Raw mean per-block score |
| `rank` | integer | Rank by stacking_index (1 = highest) |
| `score_<block_id>` | numeric | Per-block score columns |

### 15.10. `summarize_parent_haplotypes()` — data frame

| Column | Type | Description |
|--------|------|-------------|
| `id`, `block_id`, `CHR`, `start_bp`, `end_bp` | — | Identifiers |
| `allele` | character | Haplotype allele string |
| `dosage` | integer {0,1,2} | Copies carried. Phased: 0/1/2. Unphased: 0/1. |
| `allele_freq` | numeric [0,1] | Population frequency in full panel |
| `allele_effect` | numeric or NA | Effect from allele_effects argument |
| `is_rare` | logical | TRUE when allele_freq < 0.10 |

### 15.11. `compare_block_effects()` — LDxBlocks_effect_concordance

**`$concordance`** — one row per block per trait:

| Column | Type | Description |
|--------|------|-------------|
| `block_id`, `CHR`, `start_bp`, `end_bp`, `trait` | — | Identifiers |
| `n_alleles_pop1`, `n_alleles_pop2` | integer | Alleles tested in each population before intersection |
| `n_shared_alleles` | integer | Alleles present (by string match) in both populations |
| `enough_shared` | logical | `n_shared_alleles >= min_shared_alleles` |
| `effect_correlation` | numeric | Pearson r of per-allele effects across populations (NA when n_shared < 3) |
| `direction_agreement` | numeric [0,1] | Fraction of shared alleles with the same effect sign |
| `directionally_concordant` | logical | `direction_agreement >= direction_threshold` |
| `meta_effect` | numeric | IVW meta-analytic effect (weighted mean of per-allele IVW estimates) |
| `meta_SE` | numeric | SE of IVW estimate |
| `meta_z` | numeric | meta-analytic z-score |
| `meta_p` | numeric | Two-sided p-value of meta-analytic effect |
| `Q_stat` | numeric | Cochran Q heterogeneity statistic |
| `Q_df` | integer | Degrees of freedom (n_shared_alleles − 1) |
| `Q_p` | numeric | p-value of Q under chi-squared; significant = heterogeneity between populations |
| `I2` | numeric [0,100] | I² inconsistency: max(0, (Q − df) / Q × 100). > 50% = substantial heterogeneity |
| `replicated` | logical | `enough_shared AND directionally_concordant AND Q_p > 0.05` |
| `boundary_overlap_ratio` | numeric [0,1] | **Automatically computed output** from `blocks_pop1` / `blocks_pop2`: bp(intersection) / bp(union). `NA` when block tables not supplied |
| `boundary_warning` | logical | `TRUE` when `boundary_overlap_ratio < boundary_overlap_warn` (the **input parameter**, default 0.80); flags blocks where differing LD structure may make haplotype strings non-comparable |

**`$shared_alleles`** — one row per shared allele per block per trait:

| Column | Type | Description |
|--------|------|-------------|
| `allele` | character | Haplotype allele string |
| `effect_pop1`, `SE_pop1`, `p_wald_pop1` | numeric | Effect, SE, and p-value from population 1 |
| `effect_pop2`, `SE_pop2`, `p_wald_pop2` | numeric | Effect, SE, and p-value from population 2 |
| `direction_agree` | logical | Same sign in both populations? |
| `ivw_effect` | numeric | Per-allele IVW combined effect (weighted mean of pop1 and pop2 by 1/SE²) |
| `ivw_SE` | numeric | SE of IVW combined effect |
| `match_type` | character | How this block was matched: `"exact"` (same block_id), `"position"` (matched by genomic IoU ≥ `overlap_min`), `"pop1_only"` (no Pop2 block overlapped), or `NA` when no block tables supplied. Set when `block_match = "position"`. |

---

### 15.12. `compare_gwas_effects()` — LDxBlocks_effect_concordance

Same output class as `compare_block_effects()` (Section 15.11). Additional
columns in `$concordance` specific to external GWAS input:

| Column | Type | Description |
|--------|------|-------------|
| `lead_snp_pop1`, `lead_snp_pop2` | character | Lead SNP ID from each population (same SNP = same LD tag; different = different proxies of the same QTL region) |
| `lead_p_pop1`, `lead_p_pop2` | numeric | Lead SNP p-values |
| `se_derived_pop1`, `se_derived_pop2` | logical | `TRUE` when SE was derived from `BETA` and `P` via z-score rather than read directly |
| `both_pleiotropic` | logical | `TRUE` when the block is pleiotropic in both populations (requires `pleiotropic` column from `define_qtl_regions()`) |

Columns that differ in meaning from `compare_block_effects()`:

| Column | In `compare_block_effects()` | In `compare_gwas_effects()` |
|--------|------------------------------|------------------------------|
| `n_shared_alleles` | Haplotype alleles in both pops | Always 1 (one lead SNP per block) |
| `effect_correlation` | Pearson r across alleles | Always `NA` (needs ≥ 3 alleles) |
| `direction_agreement` | Fraction of alleles with same sign | 0 or 1 only |
| `Q_stat`, `Q_p`, `I2` | Cochran Q heterogeneity test | Always `NA` (df = 0) |
| `replicated` | dir_concordant AND Q_p > 0.05 | dir_concordant AND meta_p ≤ 0.05 |
| `match_type` | `"exact"` / `"position"` / `"pop1_only"` / `NA` | Same semantics as in `compare_block_effects()`. `block_match = "position"` applies to both functions. |

The `$shared_alleles` data frame contains one row per block per trait with
`lead_snp_pop1`, `lead_snp_pop2` instead of `allele`.

## 16. Memory and performance notes

### 16.1. C++ core

The eleven compiled functions in `src/ld_core.cpp` replace the most expensive
R operations. Key speedups:

- **`compute_r2_cpp()`**: ~40× over pure R for a 1,500 × 1,500 window with 500 individuals.
- **`boundary_scan_cpp()`**: eliminates ~150,000 small R-level matrix operations for a 50,000-SNP chromosome.
- **`maf_filter_cpp()`**: ~10× faster for panels > 100,000 SNPs.
- **`resolve_overlap_cpp()`**: 15,700× per-SNP reduction on chr1 via BLAS DGEMM scoring with a lazy column cache.

### 16.2. Never-full-genome memory model

- **Numeric dosage CSV**: two-pass chunked reading. Peak RAM = one 50,000-row chunk.
- **VCF and HapMap**: auto-converted to GDS cache; all access streaming via `read_chunk()`.
- **GDS and PLINK BED**: `read_chunk()` called once per sub-segment per chromosome.
- **Chromosome loop**: `rm()` and `gc(FALSE)` called after each chromosome.

### 16.3. OpenMP thread count

```r
n_thr  <- parallel::detectCores(logical = FALSE)
blocks <- run_Big_LD_all_chr(be, n_threads = n_thr)
```

Efficient scaling up to approximately 8-16 threads for `subSegmSize = 1500`.

---

## 17. Documentation

Full documentation, function reference, and tutorials:

<https://FAkohoue.github.io/LDxBlocks/>

```r
vignette("LDxBlocks-workflow", package = "LDxBlocks")
vignette("LDxBlocks-phasing",  package = "LDxBlocks")
```

---

## 18. Citation

```
LDxBlocks Development Team (2025).
LDxBlocks: Genome-Wide LD Block Detection, Haplotype Analysis, and Genomic
Prediction Features with Kinship-Adjusted Correlations.
R package version 0.3.1.
https://github.com/FAkohoue/LDxBlocks
```

Please also cite:

```
Kim S-A et al. (2018). Bioinformatics 34(4):588-596.
VanRaden PM (2008). Journal of Dairy Science 91(11):4414-4423.
Calus MPL et al. (2008). Genetics 178(1):553-561.
Nei M (1973). PNAS 70(12):3321-3323.
Blondel VD et al. (2008). J. Stat. Mech. P10008.
Traag VA et al. (2019). Scientific Reports 9:5233.
```

---

## 19. Contributing

<https://github.com/FAkohoue/LDxBlocks/issues>

Before opening a pull request: run `devtools::check()` with zero errors/warnings,
add tests in `tests/testthat/test-core.R`, rebuild docs with `devtools::document()`.

---

## 20. License

MIT + file LICENSE © Félicien Akohoue

---

## 21. References

Gao X, Starmer J, Martin ER (2008). A multiple testing correction method for
genetic association studies using correlated single nucleotide polymorphisms.
*Genetic Epidemiology* **32**:361-369. doi:10.1002/gepi.20310

Gao X, Becker LC, Becker DM, Starmer JD, Province MA (2010). Avoiding the
high Bonferroni penalty in genome-wide association studies.
*Genetic Epidemiology* **34**:100-105. doi:10.1002/gepi.20430

Gao X (2011). Multiple testing corrections for imputed SNPs.
*Genetic Epidemiology* **35**:154-158. doi:10.1002/gepi.20563

Borenstein M, Hedges LV, Higgins JPT, Rothstein HR (2009).
*Introduction to Meta-Analysis*. Wiley.

Higgins JPT, Thompson SG (2002). Quantifying heterogeneity in a meta-analysis.
*Statistics in Medicine* **21**(11):1539-1558. doi:10.1002/sim.1186

Kim S-A, Cho C-S, Kim S-R, Bull SB, Yoo Y-J (2018). A new haplotype block
detection method for dense genome sequencing data based on interval graph
modeling and dynamic programming. *Bioinformatics* **34**(4):588-596.

VanRaden PM (2008). Efficient methods to compute genomic predictions.
*Journal of Dairy Science* **91**(11):4414-4423.

Mangin B et al. (2012). Novel measures of linkage disequilibrium that correct
the bias due to population structure and relatedness. *Heredity* **108**(3):285-291.

Calus MPL et al. (2008). Accuracy of genomic selection using different methods
to define haplotypes. *Genetics* **178**(1):553-561.

de Roos APW, Hayes BJ, Goddard ME (2009). Reliability of genomic predictions
across multiple populations. *Genetics* **183**(4):1545-1553.

Nei M (1973). Analysis of gene diversity in subdivided populations.
*PNAS* **70**(12):3321-3323.

Tong J et al. (2024). Stacking beneficial haplotypes from the Vavilov wheat
collection. *Theor. Appl. Genet.* **137**:274.

Tong J et al. (2025). Haplotype stacking to improve stability of stripe rust
resistance in wheat. *Theor. Appl. Genet.* **138**:267.

Blondel VD et al. (2008). Fast unfolding of communities in large networks.
*J. Stat. Mech.* P10008.

Traag VA, Waltman L, van Eck NJ (2019). From Louvain to Leiden: guaranteeing
well-connected communities. *Scientific Reports* **9**:5233.
