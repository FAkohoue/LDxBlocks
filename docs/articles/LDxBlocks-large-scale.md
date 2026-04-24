# Large-Scale Analysis: GDS and PLINK BED Backends

All code chunks in this vignette have `eval = FALSE` because they
require either optional packages (SNPRelate, BEDMatrix) or large data
files not shipped with the package. Copy and adapt them to your own
data.

------------------------------------------------------------------------

## 1. Why scale matters: the original Big-LD limitation

The original `Big_LD()` function requires the complete genotype dataset
in RAM simultaneously. For modern WGS panels this is untenable:

| Dataset              | SNPs       | Individuals | RAM required (naive) |
|----------------------|------------|-------------|----------------------|
| Typical SNP chip     | 50,000     | 5,000       | ~2 GB                |
| Medium WGS panel     | 1,000,000  | 500         | ~4 GB                |
| Large WGS panel      | 3,000,000  | 204         | ~5 GB                |
| Rice mega-GWAS panel | 3,000,000  | 5,000       | ~120 GB              |
| Bovine WGS           | 10,000,000 | 2,000       | ~160 GB              |

## 2. How LDxBlocks solves it: never-full-genome memory model

| Format | Memory strategy | Peak RAM for 3M × 204 panel |
|----|----|----|
| VCF / HapMap | Auto-converted to GDS cache; streams per window | ~60 MB per window |
| Numeric CSV | Two-pass chunked reader; 50,000-row chunks | ~50 MB per chunk |
| GDS | [`snpgdsGetGeno()`](https://rdrr.io/pkg/SNPRelate/man/snpgdsGetGeno.html) per window | ~60 MB per window |
| PLINK BED | `BEDMatrix` OS memory-mapping | ~60 MB per window |
| bigmemory | `filebacked.big.matrix`; OS pages on demand | ~0.8 MB per window |

Peak RAM formula:

    RAM_peak ~ n_samples × subSegmSize × 8 bytes
             = 204 × 1500 × 8 = 2.4 MB  (3M-SNP rice panel)
             = 5000 × 1500 × 8 = 60 MB  (5000-sample cattle panel)

## 3. GDS backend (SNPRelate)

``` r
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("SNPRelate")
```

``` r
library(SNPRelate)
SNPRelate::snpgdsVCF2GDS(
  vcf.fn      = "mydata.vcf.gz",
  out.fn      = "mydata.gds",
  method      = "biallelic.only",
  snpfirstdim = FALSE
)
```

``` r
be_gds <- read_geno("mydata.gds")

blocks <- run_Big_LD_all_chr(
  be_gds,
  method       = "r2",
  CLQcut       = 0.70,
  leng         = 100L,
  subSegmSize  = 1500L,
  n_threads    = 16L,
  min_snps_chr = 100L,
  verbose      = TRUE
)
close_backend(be_gds)
```

## 4. WGS-scale acceleration

### 4.1 Leiden community detection (recommended)

``` r
blocks <- run_Big_LD_all_chr(
  be_gds,
  CLQmode         = "Leiden",      # polynomial O(n log n); guaranteed connected
  CLQcut          = 0.70,
  max_bp_distance = 500000L,       # skip pairs > 500 kb
  subSegmSize     = 500L,
  leng            = 50L,
  checkLargest    = TRUE,
  n_threads       = 16L,
  min_snps_chr    = 100L
)
```

### 4.2 Sparse LD computation (max_bp_distance)

At 500 kb, 70–90% of pairs are skipped, reducing O(p²) to near-O(p):

``` r
blocks <- run_Big_LD_all_chr(
  be_gds,
  CLQmode         = "Leiden",
  max_bp_distance = 500000L,
  CLQcut          = 0.70,
  n_threads       = 16L
)
```

### 4.3 bigmemory backend

``` r
be_gds <- read_geno("mydata.gds")
be_bm  <- read_geno_bigmemory(
  be_gds,
  backingfile = "mydata_bm",
  backingpath = "/data/ldxblocks",
  type        = "char",
  verbose     = TRUE
)
close_backend(be_gds)
blocks <- run_Big_LD_all_chr(be_bm, CLQmode = "Leiden", CLQcut = 0.70)
close_backend(be_bm)
```

## 5. LD decay analysis on large panels

``` r
decay <- compute_ld_decay(
  geno         = be_gds,
  sampling     = "random",
  r2_threshold = "both",
  n_pairs      = 50000L,
  max_dist     = 5000000L,
  fit_model    = "loess",
  n_threads    = 8L,
  verbose      = TRUE
)
decay$critical_r2_param
decay$decay_dist
qtl <- define_qtl_regions(gwas, blocks, snp_info, ld_decay = decay)
```

## 6. PLINK BED backend

``` r
install.packages("BEDMatrix")
```

``` r
be_bed <- read_geno("mydata.bed")
blocks <- run_Big_LD_all_chr(be_bed, method = "r2", CLQcut = 0.65,
                              n_threads = 8L, verbose = FALSE)
close_backend(be_bed)
```

## 7. Parameter selection

### 7.1 subSegmSize

| Individuals (n) | subSegmSize (w) | Window RAM |
|-----------------|-----------------|------------|
| 1,000           | 1,500           | 12 MB      |
| 5,000           | 1,500           | 60 MB      |
| 5,000           | 5,000           | 200 MB     |
| 20,000          | 1,500           | 240 MB     |

### 7.2 leng

For dense WGS panels, reduce `leng` to 50–100:

``` r
blocks <- run_Big_LD_all_chr(be_gds, method = "r2", CLQcut = 0.70,
                              leng = 50L, subSegmSize = 3000L, n_threads = 16L)
```

## 8. Haplotype analysis on large datasets

``` r
be     <- read_geno("wgs_panel.gds")
blocks <- run_Big_LD_all_chr(be, method = "r2", CLQcut = 0.70,
                              n_threads = 16L, verbose = FALSE)
haps   <- extract_haplotypes(be, be$snp_info, blocks, min_snps = 5)
div    <- compute_haplotype_diversity(haps)
feat   <- build_haplotype_feature_matrix(haps, encoding = "additive_012",
                                          scale_features = TRUE)$matrix
blues  <- read.csv("blues.csv")
pred   <- run_haplotype_prediction(
  geno_matrix = as.matrix(read_chunk(be, seq_len(be$n_snps))),
  snp_info    = be$snp_info,
  blocks      = blocks,
  blues       = blues,
  id_col      = "id",
  blue_col    = "YLD"
)
close_backend(be)
```

## 9. Troubleshooting

**GDS handle errors on Windows.** Avoid FORK-based parallelism; use
`n_threads > 1` (OpenMP within a single process) instead.

**“fewer than min_snps_chr SNPs” messages.** Set `min_snps_chr = 100` or
higher to skip scaffold chromosomes automatically.

**CLQD hangs or takes hours.** Switch to `CLQmode = "Leiden"` with
`CLQcut = 0.70` and `max_bp_distance = 500000L`.

**Processing stalls after “Segment N/N done”.** Update to v0.3.1+ which
routes overlap resolution through `resolve_overlap_cpp()` — a 15,000×
reduction in per-SNP LD computation cost.

## 10. References

- Lawrence M et al. (2013). *PLOS Computational Biology*
  **9**(8):e1003118.
- VanRaden PM (2008). *Journal of Dairy Science* **91**(11):4414-4423.
