# LDxBlocks

[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)

**LDxBlocks** is an R package for genome-wide detection of linkage disequilibrium (LD) blocks using a kinship-adjusted squared correlation metric \( rV^2 \). It implements an improved version of the Big-LD algorithm by Kim et al. (2017), making it more suitable for structured or related populations such as MAGIC or NAM designs.

## Overview

Traditional LD block detection methods can be biased in structured populations due to the presence of kinship. **LDxBlocks** addresses this by using a kinship-adjusted squared correlation metric \( rV^2 \) to compute LD, followed by graph-based clique detection and interval merging.

The package supports:
- Kinship-aware LD estimation
- Per-chromosome scalability
- Clique-based block detection
- Rare SNP handling
- Block merging across gap intervals

## Installation

To install the development version from GitHub:

```r
# Install devtools if not already installed
install.packages("devtools")

# Install LDxBlocks from GitHub
devtools::install_github("fakohoue/LDxBlocks")

## Quick Start
library(LDxBlocks)

# Example genotype matrix and SNP info
geno_matrix <- your_genotype_matrix  # matrix: individuals × SNPs
snp_info <- data.frame(
  CHR = your_chr_vector,     # Chromosome ID
  SNP = your_snp_ids,        # SNP ID (e.g., rsID)
  POS = your_snp_positions   # Base-pair position
)

# Run LD block detection
ld_blocks <- run_Big_LD_all_chr(
  geno_matrix = geno_matrix,
  snp_info = snp_info,
  CLQcut = 0.5,
  split = TRUE
)

# View result
head(ld_blocks)

## Core Functions
Function | Description
run_Big_LD_all_chr() | Main wrapper to run LD block detection chromosome by chromosome
Big_LD() | Core LD segmentation function with kinship adjustment
CLQD() | Clique detection based on rV² matrix
compute_rV2() | Computes kinship-adjusted squared correlation matrix (rV²)
get_V_inv_sqrt() | Computes the inverse square root of the kinship matrix (V⁻¹ᐟ²)

## License
This project is licensed under the MIT License. See the LICENSE file for details.
