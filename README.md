# LDxBlocks

[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)

**LDxBlocks** is an R package for genome-wide detection of linkage disequilibrium (LD) blocks using a kinship-adjusted squared correlation metric (rV²).It implements an improved, kinship-aware version of the Big-LD algorithm by Kim et al. (2018), making it well-suited for structured or related populations (e.g., MAGIC, NAM).

## Overview

Traditional LD block detection can be biased in structured populations due to kinship. **LDxBlocks** addresses this by:

- computing a kinship-adjusted LD matrix rV² (via V⁻¹ᐟ²G),
- detecting cliques in a graph built from rV²,
- merging interval cliques into contiguous LD blocks,
- optionally appending rare SNPs,
- and auto-tuning parameters to minimize unassigned GWAS markers.

## Features

- Kinship-aware LD estimation
- Per-chromosome processing
- Clique-based block detection and gap-aware merging
- Optional rare SNP handling
- Parameter tuning to reduce unassigned/forced GWAS mappings

## Core Functions

| Function                    | Description                                                              |
|-----------------------------|--------------------------------------------------------------------------|
| `run_Big_LD_all_chr()` | Main wrapper to run LD block detection chromosome by chromosome          |
| `Big_LD()`             | Core LD segmentation function with kinship adjustment                    |
| `CLQD()`               | Clique detection based on rV² matrix                                     |
| `compute_rV2()`        | Computes kinship-adjusted squared correlation matrix (rV²)               |
| `get_V_inv_sqrt()`     | Computes the inverse square root of the kinship matrix V⁻¹ᐟ²             |
| `tune_LD_params()`     | Auto-tunes parameters to minimize unassigned/forced GWAS marker mappings |

## Installation

To install the development version from GitHub:

```r

# Install devtools if needed
install.packages("devtools")

# Install LDxBlocks from GitHub
devtools::install_github("FAkohoue/LDxBlocks")

```

```r

library(LDxBlocks)

```

## Quick Start

This section simulates a small genotype matrix, SNP metadata, and runs `run_Big_LD_all_chr()`.

```r

# Simulate a small genotype matrix: n individuals x p SNPs
n  <- 80
p  <- 600
chr_sizes <- c(200, 200, 200)   # 3 chromosomes, equal SNPs
stopifnot(sum(chr_sizes) == p)

# Random minor allele frequencies between 0.05 and 0.5
mafs <- runif(p, min = 0.05, max = 0.5)

# Genotypes: Binomial(2, maf) per SNP (independent across SNPs for demo)
geno_matrix <- sapply(mafs, function(m) rbinom(n, size = 2, prob = m))
mode(geno_matrix) <- "numeric"
rownames(geno_matrix) <- paste0("ind", seq_len(n))
colnames(geno_matrix) <- paste0("snp", seq_len(p))

# SNP metadata: CHR, SNP, POS
CHR <- rep(paste0("chr", 1:3), times = chr_sizes)
POS <- ave(seq_len(p), CHR, FUN = function(i) sort(round(cumsum(runif(length(i), 500L, 1500L)))))
SNP <- colnames(geno_matrix)

snp_info <- data.frame(CHR = CHR, SNP = SNP, POS = POS, row.names = NULL)
# Ensure types are clean
snp_info$CHR <- as.character(snp_info$CHR)
snp_info$POS <- as.numeric(snp_info$POS)

head(snp_info)

```
```r

# Run LD block detection
blocks <- run_Big_LD_all_chr(
  geno_matrix = geno_matrix,
  snp_info    = snp_info,
  CLQcut      = 0.5,
  clstgap     = 40000,
  leng        = 200,
  subSegmSize = 1500,
  MAFcut      = 0.05,
  appendrare  = FALSE,
  checkLargest= FALSE,
  CLQmode     = "Density",
  rV2method   = "chol",
  split       = TRUE
)

head(blocks)

```
```r

# Basic summaries
if (!is.null(blocks) && nrow(blocks) > 0) {
  blocks$length_bp <- blocks$end.bp - blocks$start.bp + 1
  by_chr <- aggregate(length_bp ~ CHR, data = blocks, FUN = function(x) c(n = length(x), median = median(x)))
  by_chr
}
```
```r

if (!is.null(blocks) && nrow(blocks) > 0) {
  hist(blocks$end.bp - blocks$start.bp + 1,
       breaks = 30,
       main = "LD block length (bp)",
       xlab = "bp")
}

```
## Optional: Parameter Tuning

`tune_LD_params()` grid-searches over Big-LD parameters, selecting the combo that (in order) minimizes unassigned GWAS markers, forced assignments, block count, and deviation from a target median block size band.

**Note:** Tuning can be compute-intensive. Below we use a tiny grid so it runs quickly. Consider enabling parallelization for larger grids.

```r

# Create a toy GWAS marker table by sampling SNPs
set.seed(42)
gwas_subset <- snp_info[sample.int(nrow(snp_info), size = 40), c("SNP","CHR","POS")]
gwas_df <- data.frame(
  Marker = paste0("mk_", seq_len(nrow(gwas_subset))),
  CHR    = gwas_subset$CHR,
  POS    = gwas_subset$POS,
  stringsAsFactors = FALSE
)

head(gwas_df)

```
```r

# A tiny grid to keep runtime modest
grid <- expand.grid(
  CLQcut       = c(0.50, 0.60),
  clstgap      = c(2e6),
  leng         = c(500),
  subSegmSize  = c(5000),
  split        = c(TRUE),
  checkLargest = c(FALSE),
  MAFcut       = c(0.05),
  CLQmode      = c("Density"),
  rV2method    = c("chol"),
  digits       = c(6),
  appendrare   = c(FALSE),
  stringsAsFactors = FALSE
)
grid

```

```r

# Run tuner (serial for the demo). For speed on real data:
# library(future); library(future.apply); plan(multisession, workers = 8)
tune_out <- tune_LD_params(
  geno_matrix        = geno_matrix,
  snp_info           = snp_info[, c("SNP","CHR","POS")],
  gwas_df            = gwas_df[,  c("Marker","CHR","POS")],
  grid               = grid,
  chromosomes        = NULL,            # or restrict e.g., c("chr1","chr2")
  target_bp_band     = c(5e4, 5e5),
  parallel           = FALSE,
  seed               = 123,
  prefer_perfect     = TRUE,
  return_all_perfect = TRUE
)

str(tune_out$best_params)
head(tune_out$final_blocks)

```

```r

# Inspect the score table
head(tune_out$score_table[order(tune_out$score_table$n_unassigned,
                                tune_out$score_table$n_forced,
                                tune_out$score_table$n_blocks,
                                tune_out$score_table$penalty_bp), ])
                                
```

## Reproducibility Knobs

- `digits`: rounding precision used inside computations to reduce floating-point jitter (typical 6–8).

- `seed`: if set, makes any stochastic operations/tie-breakers reproducible.

- `verbose`: print progress and diagnostics.

## Session Info

```r

sessionInfo()

```

## Notes on Performance

- Runtime scales with SNP count per chromosome; consider pre-filtering on MAF and/or restricting chromosomes while tuning.

- For large cohorts, set up parallelization (`future` + `future.apply`) and adjust `grid` size.

## License

This project is licensed under the MIT License. See the LICENSE file for details.

## Citation

Kim SA, Cho CS, Kim SR, Bull SB, Yoo YJ (2018). A new haplotype block detection method for dense genome sequencing data based on interval graph modeling of clusters of highly correlated SNPs. 
**Bioinformatics** 34(3):388–397. https://doi.org/10.1093/bioinformatics/btx609.

If you use LDxBlocks in a publication, please cite the package and the above reference.

## References

Kim SA, Cho CS, Kim SR, Bull SB, Yoo YJ (2018). A new haplotype block detection method for dense genome sequencing data based on interval graph modeling of clusters of highly correlated SNPs.
**Bioinformatics** 34(3):388–397. https://doi.org/10.1093/bioinformatics/btx609.

Mangin B, Siberchicot A, Nicolas S, Doligez A, This P, and Cierco-Ayrolles C (2012) Novel measures of linkage disequilibrium that correct the bias due to population structure and relatedness. 
**Heredity** 108(3), 285-291. doi: https://doi.org/10.1038/hdy.2011.73.
