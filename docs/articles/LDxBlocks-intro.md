# Introduction to LDxBlocks

## 1. Overview

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

## 2. Why LDxBlocks? Relationship to the original Big-LD

LDxBlocks is built on the clique-based segmentation algorithm of Kim et
al. (2018), which introduced interval graph modelling of LD bins as a
principled alternative to sliding-window approaches. The mathematical
core – CLQD bin assignment, maximum-weight independent set block
construction, and Bron-Kerbosch clique enumeration – is preserved
exactly. LDxBlocks extends that foundation to address three limitations
of the original implementation that become critical for modern breeding
and genomics programmes.

### 2.1 Computational bottleneck

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

### 2.2 Memory wall

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

### 2.3 Pipeline gap

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

### 2.4 What is preserved from the original

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

### 2.5 Comparison with related tools

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
| Output | Block table | Block table + gene regions | Block table | Block table + 19 downstream functions |

------------------------------------------------------------------------

## 3. Example data

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

## 4. Step 1 — Reading genotype data

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

## 5. Step 2 — Phenotype input format

The `blues` argument accepts pre-adjusted phenotype means from a mixed
model. The format is flexible — any of four structures are accepted:

**Format 1: Named numeric vector** (single trait, simplest)

``` r
blues <- c(ind001 = 4.21, ind002 = 3.87, ind003 = 5.14)
# names = individual IDs matching rownames(geno_matrix)
```

**Format 2: Data frame, single trait**

This is the most common format when reading output from ASReml-R, lme4,
or SpATS. Any column names are accepted — just tell the function which
columns to use:

``` r
blues <- read.csv("blues.csv")  # e.g. columns: "Genotype", "YLD_BLUE"
res   <- run_haplotype_prediction(geno, snp_info, blocks,
                                   blues    = blues,
                                   id_col   = "Genotype",   # ID column name
                                   blue_col = "YLD_BLUE")   # BLUE column name
```

**Format 3: Data frame, multiple traits**

``` r
blues_mt <- read.csv("blues_mt.csv")  # columns: id, YLD, DIS, PHT
res_mt   <- run_haplotype_prediction(geno, snp_info, blocks,
                                      blues           = blues_mt,
                                      id_col          = "id",
                                      blue_cols       = c("YLD", "DIS", "PHT"),
                                      importance_rule = "any")
```

When `blue_cols = NULL` (default), all numeric non-ID columns are used
automatically.

**Format 4: Named list** (different individuals per trait)

``` r
blues <- list(
  YLD = c(ind001 = 4.21, ind002 = 3.87),
  DIS = c(ind001 = 0.32, ind003 = 0.28)  # ind002 missing, ind003 added
)
```

**Column requirements and ID matching**

| Element | Requirement |
|----|----|
| ID column | Any name — set via `id_col`. Must match `rownames(geno_matrix)` exactly (case-sensitive). |
| BLUE column(s) | Any name — set via `blue_col` (single) or `blue_cols` (multi). Must be numeric. |
| Missing values | `NA` BLUEs are silently dropped from GEBV estimation. |
| ID mismatch | Intersection of phenotyped and genotyped individuals is used; unmatched individuals trigger an informative message. |

The package ships with `ldx_blues` and `inst/extdata/example_blues.csv`
showing the expected format (columns: `id`, `YLD`, `RES`):

``` r
blues_file <- system.file("extdata", "example_blues.csv", package = "LDxBlocks")
blues      <- read.csv(blues_file)
head(blues, 3)
#>       id     YLD     RES
#> 1 ind001  0.8909 -0.8208
#> 2 ind002 -0.8711  0.9927
#> 3 ind003 -0.9154 -0.9956
```

The most common source of errors is an ID mismatch between the phenotype
file and the genotype matrix. Verify with:

``` r
data(ldx_geno)
all(blues$id %in% rownames(ldx_geno))  # should be TRUE
```

------------------------------------------------------------------------

## 6. Step 3 — LD block detection

### 6.1 Genome-wide (recommended)

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

### 6.2 Single chromosome

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

### 6.3 WGS panels: CLQmode = “Leiden” and max_bp_distance

For WGS-density panels (\> 500k SNPs), the default `CLQmode = "Density"`
(Bron-Kerbosch clique enumeration) can suffer exponential blowup — a
single 1500-SNP window with dense LD can stall for hours. Two parameters
fix this:

``` r
# Recommended for WGS panels (> 500 k SNPs)
blocks_wgs <- run_Big_LD_all_chr(
  be,
  CLQmode         = "Leiden",    # polynomial O(n log n) — guaranteed connected
  max_bp_distance = 500000L,     # skip r² for pairs > 500 kb (near-O(p))
  CLQcut          = 0.70,        # sparser LD graph
  subSegmSize     = 500L,        # smaller windows
  leng            = 50L,         # narrower boundary scan at WGS density
  checkLargest    = TRUE,        # belt-and-suspenders guard
  n_threads       = n_threads
)
```

------------------------------------------------------------------------

## 7. Step 4 — Summarising and visualising blocks

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

## 8. Step 5 — Haplotype extraction

### 8.1 Unphased mode (default)

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

### 8.2 Phased mode

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

## 9. Step 6 — Haplotype diversity

[`compute_haplotype_diversity()`](https://FAkohoue.github.io/LDxBlocks/reference/compute_haplotype_diversity.md)
works with both phased and unphased strings. For phased data each
individual contributes two gamete observations, doubling the effective
sample size for frequency estimates:

``` r
div <- compute_haplotype_diversity(haps)
head(div)
#>                block_id CHR start_bp end_bp n_snps n_ind n_haplotypes        He
#> 1    block_1_1000_25535   1     1000  25535     25   120           44 0.8488796
#> 2  block_1_81986_100878   1    81986 100878     20   120           43 0.8460784
#> 3 block_1_156776_181114   1   156776 181114     25   120           50 0.8826331
#> 4    block_2_1000_29445   2     1000  29445     30   120           59 0.8866947
#> 5  block_2_85463_104532   2    85463 104532     20   120           44 0.8535014
#> 6 block_2_160237_178996   2   160237 178996     20   120           39 0.8484594
#>    Shannon n_eff_alleles freq_dominant sweep_flag phased
#> 1 2.768290         6.321     0.3583333      FALSE  FALSE
#> 2 2.654193         6.212     0.3000000      FALSE  FALSE
#> 3 2.966925         8.018     0.2916667      FALSE  FALSE
#> 4 3.108972         8.285     0.2833333      FALSE  FALSE
#> 5 2.649231         6.510     0.2833333      FALSE  FALSE
#> 6 2.604423         6.305     0.3083333      FALSE  FALSE
```

| Metric             | Formula                | Interpretation              |
|--------------------|------------------------|-----------------------------|
| Richness ($`k`$)   | unique strings         | High = diverse              |
| $`H_e`$            | Nei (1973) corrected   | Expected heterozygosity     |
| Shannon ($`H'`$)   | $`-\sum p_i \log p_i`$ | Entropy                     |
| $`n_{\text{eff}}`$ | $`1/\sum p_i^2`$       | Effective number of alleles |
| $`f_{\max}`$       | $`\max_i p_i`$         | Near 1.0 = possible sweep   |
| `sweep_flag`       | $`f_{\max} \geq 0.90`$ | Logical sweep indicator     |

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

## 10. Step 7 — Post-GWAS QTL region definition

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
#>  [6] n_sig_markers lead_snp      lead_p        lead_beta     sig_snps     
#> [11] sig_betas     traits        n_traits      pleiotropic  
#> <0 rows> (or 0-length row.names)
```

------------------------------------------------------------------------

## 11. Step 8 — Feature matrix for genomic prediction

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

# Presence/absence 0/1 — for kernel methods or random forest
feat_pa <- build_haplotype_feature_matrix(
  haplotypes = haps,
  top_n      = 5L,
  encoding   = "presence_01"
)
dim(feat_pa)
#> [1] 120  45
```

Write the feature matrix:

``` r
# Numeric: rows=haplotypes, cols=individuals, values=0/1/2/NA dosage
write_haplotype_numeric(feat_add, "hap_matrix_numeric.csv",
                         haplotypes = haps, snp_info = ldx_snp_info)

# Character: nucleotide sequence per individual per haplotype allele
write_haplotype_character(haplotypes = haps, snp_info = ldx_snp_info,
                           out_file = "hap_matrix_character.txt")
```

### 11.1 Haplotype prediction pipeline (Tong et al. 2025)

When pre-adjusted phenotype values (BLUEs or adjusted entry means) are
available, the full Tong et al. (2025) haplotype stacking pipeline runs
in a single call. `blues` accepts any of the four formats described in
Step 2 (named numeric vector, single-trait data frame, multi-trait data
frame, or named list). See Step 2 for column name requirements:

``` r
blues <- read.csv("blues.csv")  # columns: id, YLD
pred  <- run_haplotype_prediction(
  geno_matrix = ldx_geno,
  snp_info    = ldx_snp_info,
  blocks      = blocks,
  blues       = blues,
  id_col      = "id",
  blue_col    = "YLD"
)
pred$block_importance[pred$block_importance$important, ]
sort(pred$gebv, decreasing = TRUE)  # ranked GEBVs

# Integrate GWAS + variance + diversity evidence
priority <- rank_haplotype_blocks(
  diversity   = pred$diversity,
  qtl_regions = qtl,          # from define_qtl_regions()
  pred_result = pred
)
priority$ranked_blocks[priority$ranked_blocks$rank_score >= 0.9, ]
```

### 11.2 Multi-trait prediction

Pass a data frame with multiple trait columns — or a named list of named
numeric vectors — to
[`run_haplotype_prediction()`](https://FAkohoue.github.io/LDxBlocks/reference/run_haplotype_prediction.md)
to analyse all traits simultaneously with a single shared GRM. All
traits are fitted via
[`rrBLUP::kin.blup()`](https://rdrr.io/pkg/rrBLUP/man/kin.blup.html) per
trait using the shared GRM, making cross-trait block importance values
directly comparable. Block importance is then aggregated across traits
so rankings are not specific to any single trait:

``` r
blues_df <- read.csv("blues.csv")  # columns: id, YLD, DIS, PHT
res_mt <- run_haplotype_prediction(
  geno_matrix     = ldx_geno,
  snp_info        = ldx_snp_info,
  blocks          = blocks,
  blues           = blues_df,
  id_col          = "id",
  blue_cols       = c("YLD", "DIS", "PHT"),
  importance_rule = "any"   # flag block if important for >= 1 trait
)
res_mt$solver_used   # 'rrBLUP'
# Cross-trait block importance
res_mt$block_importance[res_mt$block_importance$important_any, 
  c("block_id", "var_scaled_YLD", "var_scaled_DIS",
    "var_scaled_mean", "n_traits_important")]
```

### 11.3 Building a haplotype GRM for GBLUP

``` r
G_hap <- tcrossprod(feat_add) / ncol(feat_add)
dim(G_hap)
#> [1] 120 120
round(range(diag(G_hap)), 3)
#> [1] 0.300 2.969
```

------------------------------------------------------------------------

## 12. Step 9 — Parameter auto-tuning

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
#> $min_freq
#> [1] 0.01
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

## 13. Complete pipeline (one block)

``` r
library(LDxBlocks)

# Step 1. Open streaming backend (VCF auto-converts to GDS cache on first call)
be     <- read_geno("mydata.vcf.gz")

# Step 2. Prepare phenotype data (BLUEs from mixed model — see Step 2 of this vignette)
blues <- read.csv("blues.csv")  # columns: id (or any name), trait1, trait2 ...

# Step 3. Block detection
# Use Leiden for WGS panels: polynomial O(n log n) and guaranteed
# connected communities (Traag et al. 2019). Louvain can produce
# disconnected communities and is not recommended.
blocks <- run_Big_LD_all_chr(
  be,
  method          = "r2",
  CLQcut          = 0.70,
  CLQmode         = "Leiden",     # polynomial; guaranteed connected communities
  max_bp_distance = 500000L,      # skip pairs > 500 kb; near-O(p)
  subSegmSize     = 500L,
  leng            = 50L,
  n_threads       = 8L
)
summarise_blocks(blocks)

# 3. Haplotypes — pass backend directly: streams one chromosome at a time,
#    full genome is NEVER loaded into RAM simultaneously.
haps <- extract_haplotypes(be, be$snp_info, blocks, min_snps = 5)
div  <- compute_haplotype_diversity(haps)
write_haplotype_diversity(div, "diversity.csv")

feat <- build_haplotype_feature_matrix(haps,
                                        encoding = "additive_012",
                                        scale_features = TRUE)

# Numeric dosage matrix: rows=haplotypes, cols=individuals
# values: 0/1/2/NA for phased data; 0/1/NA for unphased data
write_haplotype_numeric(feat, "hap_matrix_numeric.csv",
                         haplotypes = haps, snp_info = be$snp_info)

# Nucleotide character matrix: shows actual allele sequence per individual
write_haplotype_character(haplotypes = haps, snp_info = be$snp_info,
                           out_file = "hap_matrix_character.txt")

# Decode haplotype strings to nucleotide sequences (for biological interpretation)
decoded <- decode_haplotype_strings(haps, snp_info = be$snp_info)
head(decoded[, c("block_id","hap_rank","dosage_string",
                  "nucleotide_sequence","frequency")])

close_backend(be)
```

------------------------------------------------------------------------

## 14. References

- Kim S-A, Cho C-S, Kim S-R, Bull SB, Yoo Y-J (2018). A new haplotype
  block detection method for dense genome sequencing data based on
  interval graph modeling and dynamic programming. *Bioinformatics*
  **34**(4):588-596. <https://doi.org/10.1093/bioinformatics/btx609>
- Difabachew YF, Frisch M, Langstroff AL, Stahl A, Wittkop B, Snowdon
  RJ, Koch M, Kirchhoff M, Csélényi L, Wolf M, Förster J, Weber S, Okoye
  UJ, Zenke-Philippi C (2023). Genomic prediction with haplotype blocks
  in wheat. *Frontiers in Plant Science* **14**:1168547.
  <https://doi.org/10.3389/fpls.2023.1168547>
- Weber SE, Frisch M, Snowdon RJ, Voss-Fels KP (2023). Haplotype blocks
  for genomic prediction: a comparative evaluation in multiple crop
  datasets. *Frontiers in Plant Science* **14**:1217589.
  <https://doi.org/10.3389/fpls.2023.1217589>
- Pook T, Schlather M, de los Campos G, Mayer M, Schoen CC, Simianer H
  (2019). HaploBlocker: Creation of subgroup-specific haplotype blocks
  and libraries. *Genetics* **212**(4):1045-1061.
  <https://doi.org/10.1534/genetics.119.302283>
- Tong J, Tarekegn ZT, Jambuthenne D, Alahmad S, Periyannan S, Hickey L,
  Dinglasan E, Hayes B (2024). Stacking beneficial haplotypes from the
  Vavilov wheat collection to accelerate breeding for multiple disease
  resistance. *Theoretical and Applied Genetics* **137**:274.
  <https://doi.org/10.1007/s00122-024-04784-w>
- Tong J et al. (2025). Haplotype stacking to improve stability of
  stripe rust resistance in wheat. *Theoretical and Applied Genetics*
  **138**:267. <https://doi.org/10.1007/s00122-025-05045-0>
- Mangin B, Siberchicot A, Nicolas S, Doligez A, This P, Cierco-Ayrolles
  C (2012). Novel measures of linkage disequilibrium that correct the
  bias due to population structure and relatedness. *Heredity*
  **108**(3):285-291. <https://doi.org/10.1038/hdy.2011.73>
- VanRaden PM (2008). Efficient methods to compute genomic predictions.
  *Journal of Dairy Science* **91**(11):4414-4423.
  <https://doi.org/10.3168/jds.2007-0980>
- Covarrubias-Pazaran G (2016). Genome-assisted prediction of
  quantitative traits using the R package sommer. *PLOS ONE*
  **11**:e0156744. <https://doi.org/10.1371/journal.pone.0156744>
- Calus MPL, Meuwissen THE, de Roos APW, Veerkamp RF (2008). Accuracy of
  genomic selection using different methods to define haplotypes.
  *Genetics* **178**(1):553-561.
  <https://doi.org/10.1534/genetics.107.080838>
- de Roos APW, Hayes BJ, Goddard ME (2009). Reliability of genomic
  predictions across multiple populations. *Genetics*
  **183**(4):1545-1553. <https://doi.org/10.1534/genetics.109.104935>
- Nei M (1973). Analysis of gene diversity in subdivided populations.
  *Proceedings of the National Academy of Sciences*
  **70**(12):3321-3323. <https://doi.org/10.1073/pnas.70.12.3321>
- Blondel VD, Guillaume J-L, Lambiotte R, Lefebvre E (2008). Fast
  unfolding of communities in large networks. *Journal of Statistical
  Mechanics: Theory and Experiment* **2008**:P10008.
  <https://doi.org/10.1088/1742-5468/2008/10/P10008>
- Traag VA, Waltman L, van Eck NJ (2019). From Louvain to Leiden:
  guaranteeing well-connected communities. *Scientific Reports*
  **9**:5233. <https://doi.org/10.1038/s41598-019-41695-z>
