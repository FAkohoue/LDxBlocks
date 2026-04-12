# Introduction to LDxBlocks

## Overview

`LDxBlocks` detects linkage disequilibrium (LD) blocks from genome-wide
SNP data and extends the pipeline into haplotype reconstruction,
diversity analysis, and genomic prediction feature engineering. Two LD
metrics are available:

- **Standard r²** – fast, no kinship correction, suitable for 10 M+
  markers and unstructured populations.
- **Kinship-adjusted rV²** – corrects for relatedness-induced LD
  inflation, for livestock, inbred lines, or family-based cohorts. See
  the *LD Metrics* vignette for the statistical detail and decision
  table.

The C++/Armadillo computational core makes the algorithm approximately
40x faster than the original Big-LD implementation (no compiled code)
for typical window sizes.

This vignette covers the complete pipeline – reading genotype data,
block detection, haplotype analysis, and parameter tuning – using the
`ldx_geno` example dataset (120 individuals, 230 SNPs, 3 chromosomes, 9
simulated LD blocks).

------------------------------------------------------------------------

## Why LDxBlocks? Relationship to the original Big-LD

LDxBlocks is built on the clique-based segmentation algorithm of Kim et
al. (2018), which introduced interval graph modelling of LD bins as a
principled alternative to sliding-window approaches. The mathematical
core – CLQD bin assignment, maximum-weight independent set block
construction, and Bron-Kerbosch clique enumeration – is preserved
exactly. LDxBlocks extends that foundation to address three limitations
of the original implementation that become critical for modern breeding
and genomics programmes.

### 1. Computational bottleneck

The original
[`Big_LD()`](https://FAkohoue.github.io/LDxBlocks/reference/Big_LD.md)
calls [`cor()`](https://rdrr.io/r/stats/cor.html) inside the
boundary-scan loop and once per
[`CLQD()`](https://FAkohoue.github.io/LDxBlocks/reference/CLQD.md) call,
with all matrix operations running through the R interpreter. For a
50,000-SNP chromosome with the default `leng = 200` this amounts to
approximately 150,000 small R-level matrix multiplications per
chromosome. For a WGS rice panel with 3 million SNPs across 12
chromosomes, the original implementation would take days on a
workstation.

LDxBlocks replaces the critical paths with seven compiled C++ functions:

``` r
# These calls are now C++/Armadillo + OpenMP -- not R loops
compute_r2_cpp(geno, digits = -1L, n_threads = 8L)   # ~40x faster
maf_filter_cpp(geno, maf_cut = 0.05)                  # ~10x faster
boundary_scan_cpp(geno, start, end, half_w, threshold) # ~20x faster
```

The outer loop of `compute_r2_cpp()` is parallelised with OpenMP. For a
1,500-SNP window with 500 individuals, the r² matrix computes in
milliseconds rather than seconds.

### 2. Memory wall

The original
[`Big_LD()`](https://FAkohoue.github.io/LDxBlocks/reference/Big_LD.md)
requires the complete genotype matrix in RAM before detection begins.
For a 10M-marker, 5,000-individual WGS dataset this is approximately 400
GB as an R `double` matrix – impossible on any workstation.

LDxBlocks enforces a strict never-full-genome memory model. All six
supported formats stream genotypes one chromosome window at a time:

``` r
# VCF auto-converts to GDS cache; detection streams one window at a time
be     <- read_geno("wgs_panel.vcf.gz")
blocks <- run_Big_LD_all_chr(be, method = "r2", n_threads = 8L)
# Peak RAM ~ n_samples x subSegmSize x 8 bytes = 60 MB for n=5000, w=1500
```

### 3. Pipeline gap

The original Big-LD stops at the block table. For breeding programmes
the immediate next questions are: What haplotypes exist? How diverse are
they? Which blocks harbour QTLs? How do I build a multi-locus genomic
prediction model? LDxBlocks provides a complete answer:

``` r
# From blocks to genomic prediction features in four lines
haps <- extract_haplotypes(be, be$snp_info, blocks, min_snps = 5)
div  <- compute_haplotype_diversity(haps)
qtl  <- define_qtl_regions(gwas_results, blocks, be$snp_info)
feat <- build_haplotype_feature_matrix(haps, top_n = 5, encoding = "additive_012")
```

### What is preserved from the original

The following components are unchanged from Kim et al. (2018):

- [`CLQD()`](https://FAkohoue.github.io/LDxBlocks/reference/CLQD.md)
  (internal): bin vector assignment via maximal clique enumeration and
  greedy density/maximal priority selection
- `constructLDblock()` (internal): maximum-weight independent set via
  dynamic programming on sorted interval sequences
- All `CLQmode`, `clstgap`, and `split` logic
- Block table column format (`start`, `end`, `start.rsID`, `end.rsID`,
  `start.bp`, `end.bp`) – drop-in compatible with tools that accept
  original Big-LD output

### Comparison with related tools

| Feature | Original Big-LD | gpart (Bioconductor) | PLINK –blocks | LDxBlocks |
|----|----|----|----|----|
| Algorithm | Clique/MWIS | Clique/MWIS + GPART | Gabriel et al. | Clique/MWIS |
| LD metric | r | r, D’ | r² | r², rV² |
| C++ core | No | Partial | Yes | Yes (7 functions) |
| WGS streaming | No | Partial | Yes | Yes (all formats) |
| Kinship correction | No | No | No | Yes (rV²) |
| Haplotype analysis | No | No | No | Yes |
| Genomic prediction | No | No | No | Yes |
| GWAS-driven tuning | No | No | No | Yes |
| Input formats | R matrix | PLINK, VCF | PLINK | 6 formats |
| Output | Block table | Block table + gene regions | Block table | Block table + 8 downstream functions |

------------------------------------------------------------------------

## Example data

Four datasets are shipped with the package, all generated by
`data-raw/generate_example_data.R`:

``` r
dim(ldx_geno)             # 120 individuals x 230 SNPs
#> [1] 120 230
head(ldx_snp_info)        # SNP, CHR, POS, REF, ALT
#>      SNP CHR  POS REF ALT
#> 1 rs1001   1 1000   G   C
#> 2 rs1002   1 2134   C   G
#> 3 rs1003   1 3089   T   A
#> 4 rs1004   1 4067   A   G
#> 5 rs1005   1 5246   G   A
#> 6 rs1006   1 6312   C   T
table(ldx_snp_info$CHR)   # chr1, chr2, chr3
#> 
#>  1  2  3 
#> 80 80 70
nrow(ldx_gwas)            # toy GWAS markers
#> [1] 20
```

Each chromosome has three LD blocks separated by inter-block gaps and
short stretches of low-LD singleton SNPs.

------------------------------------------------------------------------

## Step 1 — Reading genotype data

[`read_geno()`](https://FAkohoue.github.io/LDxBlocks/reference/read_geno.md)
auto-detects the file format from the extension and returns an
`LDxBlocks_backend` object. Six formats are supported:

| Format             | Extension                 | `format =`  |
|--------------------|---------------------------|-------------|
| Numeric dosage     | `.csv`, `.txt`            | `"numeric"` |
| HapMap             | `.hmp.txt`                | `"hapmap"`  |
| VCF / bgzipped VCF | `.vcf`, `.vcf.gz`         | `"vcf"`     |
| SNPRelate GDS      | `.gds`                    | `"gds"`     |
| PLINK binary       | `.bed` (+ `.bim`, `.fam`) | `"bed"`     |
| R matrix           | (in-memory)               | `"matrix"`  |

``` r
be_vcf <- read_geno("mydata.vcf.gz")
be_csv <- read_geno("mydata.csv")
be_hmp <- read_geno("mydata.hmp.txt")
be_gds <- read_geno("mydata.gds")    # SNPRelate GDS; streams SNP windows
be_bed <- read_geno("mydata.bed")    # PLINK BED; memory-mapped
```

``` r
be <- read_geno(ldx_geno, format = "matrix", snp_info = ldx_snp_info)
be
#> LDxBlocks backend
#>   Type      : matrix 
#>   Samples   : 120 
#>   SNPs      : 230 
#>   Chr       : 1, 2, 3
```

``` r
be$n_samples
#> [1] 120
be$n_snps
#> [1] 230
be$sample_ids[1:4]
#> [1] "ind001" "ind002" "ind003" "ind004"
head(be$snp_info)
#>      SNP CHR  POS REF ALT
#> 1 rs1001   1 1000   G   C
#> 2 rs1002   1 2134   C   G
#> 3 rs1003   1 3089   T   A
#> 4 rs1004   1 4067   A   G
#> 5 rs1005   1 5246   G   A
#> 6 rs1006   1 6312   C   T
chunk <- read_chunk(be, 1:30)
dim(chunk)
#> [1] 120  30
close_backend(be)
```

------------------------------------------------------------------------

## Step 2 — LD block detection

### Genome-wide (recommended)

``` r
blocks <- run_Big_LD_all_chr(
  ldx_geno,
  snp_info    = ldx_snp_info,
  method      = "r2",
  CLQcut      = 0.55,
  leng        = 15L,
  subSegmSize = 100L,
  n_threads   = 1L,
  verbose     = FALSE
)
nrow(blocks)
#> [1] 9
head(blocks)
#>   start end start.rsID end.rsID start.bp end.bp CHR length_bp
#> 1     1  25     rs1001   rs1025     1000  25535   1     24536
#> 2    31  50     rs1031   rs1050    81986 100878   1     18893
#> 3    56  80     rs1056   rs1080   156776 181114   1     24339
#> 4     1  30     rs2001   rs2030     1000  29445   2     28446
#> 5    36  55     rs2036   rs2055    85463 104532   2     19070
#> 6    61  80     rs2061   rs2080   160237 178996   2     18760
```

### Single chromosome

Use the `chr` parameter of
[`run_Big_LD_all_chr()`](https://FAkohoue.github.io/LDxBlocks/reference/run_Big_LD_all_chr.md)
to restrict detection to one or more specific chromosomes. This uses the
same streaming memory model and produces the same output columns as a
genome-wide run.

``` r
blocks_chr1 <- run_Big_LD_all_chr(
  ldx_geno,
  snp_info    = ldx_snp_info,
  method      = "r2",
  CLQcut      = 0.55,
  leng        = 15L,
  subSegmSize = 100L,
  chr         = "1",         # restrict to chromosome 1 only
  n_threads   = 1L,
  verbose     = FALSE
)
nrow(blocks_chr1)
#> [1] 3
unique(blocks_chr1$CHR)      # confirms only chr 1 processed
#> [1] "1"
```

------------------------------------------------------------------------

## Step 3 — Summarising and visualising blocks

``` r
summarise_blocks(blocks)
#>      CHR n_blocks min_bp median_bp  mean_bp max_bp total_bp_covered
#> 1      1        3  18893     24339 22589.33  24536            67768
#> 2      2        3  18760     19070 22092.00  28446            66276
#> 3      3        3  18995     19653 19481.33  19796            58444
#> 4 GENOME        9  18760     19653 21387.56  28446           192488
```

``` r
if (requireNamespace("ggplot2", quietly = TRUE)) {
  print(plot_ld_blocks(blocks, colour_by = "length_bp", mb_scale = TRUE))
}
```

![LD block structure coloured by block
size](LDxBlocks-intro_files/figure-html/plot-1.png)

LD block structure coloured by block size

------------------------------------------------------------------------

## Step 4 — Haplotype extraction

### Unphased mode (default)

[`extract_haplotypes()`](https://FAkohoue.github.io/LDxBlocks/reference/extract_haplotypes.md)
accepts either a numeric dosage matrix (0/1/2/NA) or a phased list from
[`read_phased_vcf()`](https://FAkohoue.github.io/LDxBlocks/reference/read_phased_vcf.md)
/
[`phase_with_pedigree()`](https://FAkohoue.github.io/LDxBlocks/reference/phase_with_pedigree.md).
In unphased mode each individual receives one diploid allele string per
block:

``` r
# For large real datasets, pass the backend object directly:
#   haps <- extract_haplotypes(be, be$snp_info, blocks, min_snps = 5L)
# This streams one chromosome at a time — the full genome is never in RAM.
#
# For the small example dataset we use the in-memory matrix path:
haps <- extract_haplotypes(
  geno     = ldx_geno,
  snp_info = ldx_snp_info,
  blocks   = blocks,
  min_snps = 5L,
  na_char  = "."
)
length(haps)
#> [1] 9
attr(haps, "block_info")
#>                block_id CHR start_bp end_bp n_snps phased
#> 1    block_1_1000_25535   1     1000  25535     25  FALSE
#> 2  block_1_81986_100878   1    81986 100878     20  FALSE
#> 3 block_1_156776_181114   1   156776 181114     25  FALSE
#> 4    block_2_1000_29445   2     1000  29445     30  FALSE
#> 5  block_2_85463_104532   2    85463 104532     20  FALSE
#> 6 block_2_160237_178996   2   160237 178996     20  FALSE
#> 7    block_3_1000_19994   3     1000  19994     20  FALSE
#> 8   block_3_76186_95838   3    76186  95838     20  FALSE
#> 9 block_3_151654_171449   3   151654 171449     20  FALSE
head(haps[[1]])
#>                      ind001                      ind002 
#> "0000000000000000000000000" "1111111111111111111111111" 
#>                      ind003                      ind004 
#> "1111111111111111111111111" "1111111111111111111111111" 
#>                      ind005                      ind006 
#> "0000000000000000000000000" "1111101111211111111111111"
```

### Phased mode

When WGS data have been phased externally (Beagle, SHAPEIT, pedigree),
load the phased VCF first and pass the resulting list to
[`extract_haplotypes()`](https://FAkohoue.github.io/LDxBlocks/reference/extract_haplotypes.md).
Phased mode is detected automatically from the `hap1`/`hap2` structure:

``` r
# Statistical phasing via Beagle 5.x
phase_with_beagle(input_vcf = "raw.vcf.gz", out_prefix = "phased",
                  nthreads = 8L)
phased <- read_phased_vcf("phased.vcf.gz", min_maf = 0.05)

# Pedigree-based phasing
ped    <- data.frame(id = c("off1","sire1","dam1"),
                     sire = c("sire1", NA, NA),
                     dam  = c("dam1",  NA, NA))
phased <- phase_with_pedigree(ldx_geno, pedigree = ped)

# Extract haplotypes — phased strings: "011|100"
haps_p <- extract_haplotypes(phased, ldx_snp_info, blocks, min_snps = 5)
head(haps_p[[1]])   # "011|100" format
```

------------------------------------------------------------------------

## Step 5 — Haplotype diversity

[`compute_haplotype_diversity()`](https://FAkohoue.github.io/LDxBlocks/reference/compute_haplotype_diversity.md)
works with both phased and unphased strings. For phased data each
individual contributes two gamete observations, doubling the effective
sample size for frequency estimates:

``` r
div <- compute_haplotype_diversity(haps)
head(div)
#>                block_id CHR start_bp end_bp n_snps n_ind n_haplotypes        He
#> 1    block_1_1000_25535   1     1000  25535     25   120           44 0.8418056
#> 2  block_1_81986_100878   1    81986 100878     20   120           43 0.8390278
#> 3 block_1_156776_181114   1   156776 181114     25   120           50 0.8752778
#> 4    block_2_1000_29445   2     1000  29445     30   120           59 0.8793056
#> 5  block_2_85463_104532   2    85463 104532     20   120           44 0.8463889
#> 6 block_2_160237_178996   2   160237 178996     20   120           39 0.8413889
#>    Shannon freq_dominant phased
#> 1 2.768290     0.3583333  FALSE
#> 2 2.654193     0.3000000  FALSE
#> 3 2.966925     0.2916667  FALSE
#> 4 3.108972     0.2833333  FALSE
#> 5 2.649231     0.2833333  FALSE
#> 6 2.604423     0.3083333  FALSE
```

| Metric           | Formula                | Interpretation              |
|------------------|------------------------|-----------------------------|
| Richness ($`k`$) | unique strings         | High = diverse              |
| $`H_e`$          | $`1 - \sum p_i^2`$     | Nei expected heterozygosity |
| Shannon ($`H'`$) | $`-\sum p_i \log p_i`$ | Entropy                     |
| $`f_{\max}`$     | $`\max_i p_i`$         | Near 1.0 = sweep            |

``` r
if (requireNamespace("ggplot2", quietly = TRUE)) {
  ggplot2::ggplot(div, ggplot2::aes(x = block_id, y = He)) +
    ggplot2::geom_col(fill = "#1D9E75") +
    ggplot2::theme_minimal() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle=90, size=8)) +
    ggplot2::labs(x=NULL, y="He", title="Haplotype diversity per LD block")
}
```

![Expected heterozygosity per LD
block](LDxBlocks-intro_files/figure-html/diversity-plot-1.png)

Expected heterozygosity per LD block

Write diversity results to disk:

``` r
write_haplotype_diversity(div, "haplotype_diversity.csv", append_summary = TRUE)
```

------------------------------------------------------------------------

## Step 6 — Post-GWAS QTL region definition

When GWAS results are available, map significant markers onto LD blocks
to define QTL regions and flag pleiotropic blocks:

``` r
# ldx_gwas has a `trait` column (TraitA / TraitB) enabling pleiotropic detection.
# For single-trait GWAS simply omit trait_col — the function still works,
# setting n_traits = 1 and pleiotropic = FALSE for all blocks.
qtl <- define_qtl_regions(
  gwas_results = ldx_gwas,
  blocks       = blocks,
  snp_info     = ldx_snp_info,
  p_threshold  = 5e-8,
  trait_col    = "trait"      # omit this argument for single-trait GWAS
)
head(qtl[, c("block_id","CHR","n_snps_block","n_sig_markers",
             "lead_snp","traits","pleiotropic")])
#>                block_id CHR n_snps_block n_sig_markers lead_snp traits
#> 1 block_1_156776_181114   1           25             1   rs1070 TraitA
#>   pleiotropic
#> 1       FALSE
subset(qtl, pleiotropic)      # blocks with hits from both TraitA and TraitB
#>  [1] block_id      CHR           start_bp      end_bp        n_snps_block 
#>  [6] n_sig_markers lead_snp      lead_p        traits        n_traits     
#> [11] pleiotropic  
#> <0 rows> (or 0-length row.names)
```

------------------------------------------------------------------------

## Step 7 — Feature matrix for genomic prediction

[`build_haplotype_feature_matrix()`](https://FAkohoue.github.io/LDxBlocks/reference/build_haplotype_feature_matrix.md)
supports two encoding schemes and both phased and unphased input:

``` r
# Additive 0/1/2 (phased) or 0/2 (unphased) — recommended for GBLUP/rrBLUP
feat_add <- build_haplotype_feature_matrix(
  haplotypes = haps,
  top_n      = 5L,
  encoding   = "additive_012",
  scale_features = TRUE
)
dim(feat_add)
#> [1] 120  45

# Presence/absence 0/2 — for kernel methods or random forest
feat_pa <- build_haplotype_feature_matrix(
  haplotypes = haps,
  top_n      = 5L,
  encoding   = "presence_02"
)
dim(feat_pa)
#> [1] 120  45
```

Write the feature matrix in numeric or HapMap format:

``` r
write_haplotype_numeric(feat_add, "hap_matrix_numeric.csv")
write_haplotype_hapmap(feat_add, "hap_matrix.hmp.txt")
```

### Building a haplotype GRM for GBLUP

``` r
G_hap <- tcrossprod(feat_add) / ncol(feat_add)
dim(G_hap)
#> [1] 120 120
round(range(diag(G_hap)), 3)
#> [1] 0.300 2.969
```

------------------------------------------------------------------------

## Step 8 — Parameter auto-tuning

``` r
grid <- expand.grid(
  CLQcut = c(0.45, 0.55, 0.65), clstgap = 1e5, leng = 15,
  subSegmSize = 100, split = FALSE, checkLargest = FALSE,
  MAFcut = 0.05, CLQmode = "Density", kin_method = "chol",
  digits = -1, appendrare = FALSE, stringsAsFactors = FALSE
)

result <- tune_LD_params(
  geno_matrix    = ldx_geno,
  snp_info       = ldx_snp_info,
  gwas_df        = ldx_gwas,
  grid           = grid,
  prefer_perfect = TRUE,
  target_bp_band = c(5e4, 5e5),
  seed           = 42L
)
result$best_params
#> $CLQcut
#> [1] 0.45
#> 
#> $clstgap
#> [1] 1e+05
#> 
#> $leng
#> [1] 15
#> 
#> $subSegmSize
#> [1] 100
#> 
#> $split
#> [1] FALSE
#> 
#> $checkLargest
#> [1] FALSE
#> 
#> $MAFcut
#> [1] 0.05
#> 
#> $CLQmode
#> [1] "Density"
#> 
#> $kin_method
#> [1] "chol"
#> 
#> $digits
#> [1] -1
#> 
#> $appendrare
#> [1] FALSE
result$score_table[, c("CLQcut","n_unassigned","n_forced","n_blocks")]
#>   CLQcut n_unassigned n_forced n_blocks
#> 1   0.45            0        0        9
#> 2   0.55            0        0        9
#> 3   0.65            0        0        9
```

------------------------------------------------------------------------

## Complete pipeline (one block)

``` r
library(LDxBlocks)

# 1. Open streaming backend (VCF auto-converts to GDS cache on first call)
be     <- read_geno("mydata.vcf.gz")

# 2. Block detection — chromosome-by-chromosome, gc() after each
blocks <- run_Big_LD_all_chr(be, method = "r2", CLQcut = 0.70, n_threads = 8L)
summarise_blocks(blocks)

# 3. Haplotypes — pass backend directly: streams one chromosome at a time,
#    full genome is NEVER loaded into RAM simultaneously.
haps <- extract_haplotypes(be, be$snp_info, blocks, min_snps = 5)
div  <- compute_haplotype_diversity(haps)
write_haplotype_diversity(div, "diversity.csv")

feat <- build_haplotype_feature_matrix(haps, top_n = 5,
                                        encoding = "additive_012",
                                        scale_features = TRUE)
write_haplotype_numeric(feat, "hap_matrix.csv")
write_haplotype_hapmap(feat, "hap_matrix.hmp.txt")

close_backend(be)
```

------------------------------------------------------------------------

## References

- Kim S-A et al. (2018) *Bioinformatics* **34**(4):588–596.
- VanRaden PM (2008) *J. Dairy Sci.* **91**(11):4414–4423.
- Calus MPL et al. (2008) *Genetics* **178**(1):553–561.
- de Roos APW et al. (2009) *Genetics* **183**(4):1545–1553.
- Nei M (1973) *PNAS* **70**(12):3321–3323.
- Tong J et al. (2024) *Theor Appl Genet* **137**:274.
  <https://doi.org/10.1007/s00122-024-04784-w>
