# Large-Scale Analysis: GDS and PLINK BED Backends

All code chunks in this vignette have `eval = FALSE` because they
require either optional packages (SNPRelate, BEDMatrix) or large data
files that are not shipped with the package. Copy and adapt them to your
own data.

------------------------------------------------------------------------

## 1. Why scale matters: the original Big-LD limitation

The original `Big_LD()` function (Kim et al. 2018) accepts a plain R
numeric matrix – the entire genotype dataset must be in RAM
simultaneously. For a 500-individual panel with 50,000 SNPs this works
fine. For modern breeding programmes operating with whole-genome
sequencing data, the numbers become untenable:

| Dataset              | SNPs       | Individuals | RAM required (naive) |
|----------------------|------------|-------------|----------------------|
| Typical SNP chip     | 50,000     | 5,000       | ~2 GB                |
| Medium WGS panel     | 1,000,000  | 500         | ~4 GB                |
| Large WGS panel      | 3,000,000  | 204         | ~5 GB                |
| Rice mega-GWAS panel | 3,000,000  | 5,000       | ~120 GB              |
| Bovine WGS           | 10,000,000 | 2,000       | ~160 GB              |

The 120 GB rice example is from the motivating use case for LDxBlocks:
`second_variant_calling_merged.snpsOnly.maf05.miss20.vcf.gz`, a 3M-SNP,
204-sample WGS panel where the original Big-LD would require loading the
entire genotype matrix before a single block boundary is computed.

## 2. How LDxBlocks solves it: never-full-genome memory model

LDxBlocks enforces a strict memory contract: **the full genotype matrix
is never held in RAM at once, for any dataset size or format.**

| Format | Memory strategy | Peak RAM for 3M x 204 panel |
|----|----|----|
| VCF / HapMap | Auto-converted to GDS cache; streams per chromosome window | ~60 MB per window |
| Numeric CSV | Two-pass chunked pre-allocated reader; 50,000-row chunks | ~50 MB per chunk |
| GDS | [`snpgdsGetGeno()`](https://rdrr.io/pkg/SNPRelate/man/snpgdsGetGeno.html) per window | ~60 MB per window |
| PLINK BED | `BEDMatrix` OS memory-mapping; page faults load only requested bytes | ~60 MB per window |
| bigmemory | `filebacked.big.matrix`; OS pages on demand via [`read_geno_bigmemory()`](https://FAkohoue.github.io/LDxBlocks/reference/read_geno_bigmemory.md) | ~0.8 MB per window (204 × 500 × 8) |

The peak RAM formula for file-backed backends is:

    RAM_peak ~ n_samples x subSegmSize x 8 bytes
             = 204 x 1500 x 8 = 2.4 MB per window (3M-SNP rice panel)
             = 5000 x 1500 x 8 = 60 MB per window (5000-sample cattle panel)

Compare to the naive approach: `204 x 3,000,000 x 8 = 4.9 GB` just for
the matrix, before any LD computation allocates additional memory.

In all cases
[`run_Big_LD_all_chr()`](https://FAkohoue.github.io/LDxBlocks/reference/run_Big_LD_all_chr.md)
calls `read_chunk(backend, col_idx)` to obtain exactly the SNP window it
needs, frees it with [`rm()`](https://rdrr.io/r/base/rm.html) and
`gc(FALSE)` before moving to the next window, and never accumulates more
than one chromosome in RAM at any point.

Two file-backed backends handle the on-disk streaming:

| Backend | Format           | Package required                    |
|---------|------------------|-------------------------------------|
| `"gds"` | SNPRelate GDS    | `BiocManager::install("SNPRelate")` |
| `"bed"` | PLINK binary BED | `install.packages("BEDMatrix")`     |

Both backends are read-only memory-mapped: at any point only the columns
requested by the current
[`read_chunk()`](https://FAkohoue.github.io/LDxBlocks/reference/read_chunk.md)
call are in RAM.

## 3. Memory model

Peak RAM with a file-backed backend:

``` math
\text{RAM}_{\text{peak}} \approx n \times w \times 8 \text{ bytes}
```

where $`n`$ = individuals and $`w`$ = `subSegmSize`. For $`n = 5{,}000`$
and the default `subSegmSize = 1500`:

``` math
5{,}000 \times 1{,}500 \times 8 = 60 \text{ MB per window}
```

Compare to the plain-matrix path:
$`5{,}000 \times 10{,}000{,}000 \times 8
= 400 \text{ GB}`$.

------------------------------------------------------------------------

## 4. GDS backend (SNPRelate)

### 4.1 Installation

``` r
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("SNPRelate")
```

### 4.2 Converting VCF to GDS

Convert once; use for all subsequent analyses. The GDS format supports
random-access compressed storage, making repeated window reads fast:

``` r
library(SNPRelate)
# Convert VCF to GDS once; all subsequent analyses stream from the GDS file
SNPRelate::snpgdsVCF2GDS(
  vcf.fn      = "mydata.vcf.gz",
  out.fn      = "mydata.gds",
  method      = "biallelic.only",
  snpfirstdim = FALSE
)
```

### 4.3 Block detection from GDS

``` r
be_gds <- read_geno("mydata.gds")
be_gds
# LDxBlocks backend
#   Type    : gds
#   Samples : 5000
#   SNPs    : 10241680
#   Chr     : 1, 2, ..., 29

blocks <- run_Big_LD_all_chr(
  be_gds,
  method       = "r2",
  CLQcut       = 0.70,
  leng         = 100L,      # reduce for dense WGS panels
  subSegmSize  = 1500L,     # controls peak RAM
  n_threads    = 16L,
  min_snps_chr = 100L,      # skip scaffold chromosomes
  verbose      = TRUE
)

close_backend(be_gds)
```

### 4.4 What happens inside read_chunk() for GDS

Each call executes:

``` r
dos <- SNPRelate::snpgdsGetGeno(
  gds,
  snp.id      = snp_int_ids[col_idx],
  sample.id   = sample_ids,
  snpfirstdim = FALSE,   # returns samples x SNPs
  with.id     = FALSE
)
dos <- 2L - dos  # convert REF count to ALT dosage convention
```

For `subSegmSize = 1500` and a 30-chromosome bovine genome with 10 M
markers, this is approximately
$`\lceil 10{,}000{,}000 / 1{,}500 \rceil \approx 6{,}700`$ GDS reads
total.

------------------------------------------------------------------------

## 5. WGS-scale acceleration (v0.3.1)

Three features address the exponential-blowup and I/O bottlenecks
observed when running the default parameters on 3M+ SNP panels.

### 5.1 Louvain/Leiden community detection (CLQmode) – Leiden recommended

The default `CLQmode = "Density"` uses Bron-Kerbosch clique enumeration,
which has exponential worst-case complexity. On a 3M-SNP WGS panel with
`CLQcut = 0.5` and `subSegmSize = 1500`, a single window can find 4+
million maximal cliques and run for hours. `CLQmode = "Louvain"` and
`CLQmode = "Leiden"` replace this with O(n log n) community detection:

``` r
# Recommended for any WGS panel > 500 k SNPs per chromosome
blocks <- run_Big_LD_all_chr(
  be_gds,
  CLQmode         = "Leiden",    # polynomial — guaranteed connected communities
  CLQcut          = 0.70,        # sparser graph reduces community count
  max_bp_distance = 500000L,     # sparse LD: skip pairs > 500 kb
  subSegmSize     = 500L,        # smaller windows for safety
  leng            = 50L,         # narrower boundary scan at WGS density
  checkLargest    = TRUE,        # extra guard if Density mode is used
  n_threads       = 16L,
  min_snps_chr    = 100L
)
```

Edge weights in the community graph are inversely proportional to
base-pair distance, so local LD structure is respected. Communities with
a single member become singletons (NA), exactly matching the Density
path output format.

### 5.2 Sparse LD computation (max_bp_distance)

When `max_bp_distance > 0`, only SNP pairs within that physical distance
have their r² computed via `compute_r2_sparse_cpp()`. Pairs beyond the
threshold are set to 0 in the adjacency matrix — a valid assumption
because long-range r² decays to near-zero at WGS density. At 500 kb,
70–90% of pairs are skipped, reducing O(p²) LD computation to near-O(p):

``` r
# 500 kb cutoff: appropriate for most plant and animal WGS panels
blocks <- run_Big_LD_all_chr(
  be_gds,
  CLQmode         = "Leiden",    # Leiden preferred: guaranteed connected communities
  max_bp_distance = 500000L,
  CLQcut          = 0.70,
  n_threads       = 16L
)
```

### 5.3 bigmemory backend (read_geno_bigmemory())

For cases where even the filtered matrix exceeds available RAM,
[`read_geno_bigmemory()`](https://FAkohoue.github.io/LDxBlocks/reference/read_geno_bigmemory.md)
writes the genotype matrix to a binary file (`filebacked.big.matrix`)
and wraps it in the standard backend interface.
[`read_chunk()`](https://FAkohoue.github.io/LDxBlocks/reference/read_chunk.md)
retrieves columns via OS page faults — only the requested bytes are
loaded. Storage type `"char"` (1 byte per cell) saves 8× over double for
0/1/2 dosage data. Backing files persist across sessions:

``` r
# Requires: install.packages('bigmemory')
be_gds <- read_geno('mydata.gds')

# First run: convert to memory-mapped store (~5 min for 3M SNPs)
be_bm <- read_geno_bigmemory(
  be_gds,
  backingfile = 'mydata_bm',       # .bin + .desc created here
  backingpath = '/data/ldxblocks', # persistent directory
  type        = 'char',            # 1 byte/cell; 8x smaller than double
  verbose     = TRUE
)
close_backend(be_gds)

# Subsequent runs: reattach without reloading
be_bm <- read_geno_bigmemory(
  '/data/ldxblocks/mydata_bm.desc',
  snp_info = be_bm$snp_info  # supply separately on reattach
)

# All standard functions work identically
blocks <- run_Big_LD_all_chr(be_bm, CLQmode = 'Leiden', CLQcut = 0.70)
close_backend(be_bm)
```

------------------------------------------------------------------------

## 6. LD decay analysis on large panels

[`compute_ld_decay()`](https://FAkohoue.github.io/LDxBlocks/reference/compute_ld_decay.md)
uses the same backend interface as block detection. It never loads a
full chromosome – pair indices are sampled from SNP positions only, then
[`read_chunk()`](https://FAkohoue.github.io/LDxBlocks/reference/read_chunk.md)
is called once for the unique SNPs involved:

``` r
# For GDS or bigmemory backends, read_chunk() is called for sampled SNPs only.
# For 50k random pairs on a 500k-SNP chromosome: ~10k unique SNPs -> ~40 MB.
decay <- compute_ld_decay(
  geno         = be_gds,      # any backend type
  sampling     = "random",
  r2_threshold = "both",      # parametric threshold: unlinked-marker 95th pctile
  n_pairs      = 50000L,
  max_dist     = 5000000L,    # 5 Mb
  fit_model    = "loess",
  n_threads    = 8L,          # OpenMP for compute_r2_sparse_cpp
  verbose      = TRUE
)
decay$critical_r2_param   # background kinship LD level (use rV2 if > 0.05)
decay$decay_dist           # per-chromosome decay distances for candidate gene search

# Visualise the decay curve
if (requireNamespace("ggplot2", quietly = TRUE))
  print(plot_ld_decay(decay, plot_threshold = TRUE, facet = TRUE))

# Pass to define_qtl_regions for biologically justified GWAS windows
qtl <- define_qtl_regions(gwas, blocks, snp_info, ld_decay = decay)
# qtl$candidate_region_start / end ready for BioMart / Ensembl Plants
```

------------------------------------------------------------------------

## 7. PLINK BED backend (BEDMatrix)

### 4.1 Installation

``` r
install.packages("BEDMatrix")
```

### 6.2 Block detection from BED

The `.bim` and `.fam` files must exist at the same path stem:

``` r
be_bed <- read_geno("mydata.bed")   # reads .bim and .fam automatically
be_bed
# LDxBlocks backend
#   Type    : bed
#   Samples : 2847
#   SNPs    : 658423
#   Chr     : 1, 2, ..., 26

blocks <- run_Big_LD_all_chr(
  be_bed,
  method    = "r2",
  CLQcut    = 0.65,
  n_threads = 8L,
  verbose   = FALSE
)

close_backend(be_bed)
```

`BEDMatrix` memory-maps the `.bed` file at the OS level. Row and column
indexing triggers page faults that load only the requested bytes from
disk.

------------------------------------------------------------------------

## 7. Parameter selection for large panels

### 7.1 subSegmSize

The single most important memory parameter. The C++ r² kernel allocates
an $`n \times w`$ genotype window and a $`w \times w`$ LD matrix:

| Individuals (n) | subSegmSize (w) | Window RAM | LD matrix RAM |
|-----------------|-----------------|------------|---------------|
| 1,000           | 1,500           | 12 MB      | 18 MB         |
| 5,000           | 1,500           | 60 MB      | 18 MB         |
| 5,000           | 5,000           | 200 MB     | 200 MB        |
| 20,000          | 1,500           | 240 MB     | 18 MB         |

For very large $`n`$, keep `subSegmSize` at 1,500 or reduce it. For
moderate $`n`$, increasing to 3,000–5,000 reduces GDS reads at the cost
of more RAM.

### 7.2 leng

The boundary-scan half-window in SNPs. `boundary_scan_cpp()` scans
`2 * leng` columns per candidate cut position, with cost scaling as
$`O(p \times \text{leng}^2)`$ per chromosome. For a 500 k-SNP chromosome
with `leng = 200` this is ~20 billion element comparisons. Reduce `leng`
to 50–100 for very dense panels:

``` r
# Dense WGS panel: reduce leng, increase subSegmSize
blocks <- run_Big_LD_all_chr(
  be_gds,
  method      = "r2",
  CLQcut      = 0.70,
  leng        = 50L,      # narrower boundary scan
  subSegmSize = 3000L,    # fewer GDS reads per chromosome
  n_threads   = 16L
)
```

### 7.3 n_threads

OpenMP thread scaling is efficient up to approximately 8–16 threads for
typical window sizes. Benchmark on your hardware:

``` r
times <- sapply(c(1L, 2L, 4L, 8L, 16L), function(nt) {
  system.time(
    run_Big_LD_all_chr(be_gds, method = "r2", CLQcut = 0.7,
                        n_threads = nt, verbose = FALSE)
  )["elapsed"]
})
names(times) <- paste0(c(1,2,4,8,16), " threads")
round(times, 1)
```

------------------------------------------------------------------------

## 8. Haplotype analysis on large datasets

For large panels, extract haplotypes one chromosome at a time to avoid
loading the full genome matrix:

``` r
be     <- read_geno("wgs_panel.gds")
blocks <- run_Big_LD_all_chr(be, method = "r2", CLQcut = 0.70,
                              n_threads = 16L, verbose = FALSE)

# Pass the backend directly — extract_haplotypes() detects LDxBlocks_backend
# and streams one chromosome at a time internally.
# Each chromosome is loaded, processed, rm()'d and gc(FALSE) called before
# the next chromosome. The full genome is NEVER in RAM simultaneously.
haps <- extract_haplotypes(be, be$snp_info, blocks, min_snps = 5)
div  <- compute_haplotype_diversity(haps)

# Feature matrix and output (haps is a list of strings, always small in RAM)
feat <- build_haplotype_feature_matrix(haps,
                                        encoding = "additive_012",
                                        scale_features = TRUE)
write_haplotype_numeric(feat, "haplotype_matrix_dosage.csv",
                         haplotypes = haps, snp_info = be$snp_info)
write_haplotype_character(haplotypes = haps, snp_info = be$snp_info,
                           out_file = "haplotype_matrix_nucleotide.txt")
write_haplotype_diversity(div, "haplotype_diversity.csv")

# Genomic prediction — supply BLUEs in any of the four supported formats.
# For format details and ID matching rules, see the Introduction vignette
# (Step 2: Phenotype input format) or ?run_haplotype_prediction.
blues <- read.csv("blues.csv")  # columns: id, YLD (any column names)
pred  <- run_haplotype_prediction(
  geno_matrix = as.matrix(read_chunk(be, seq_len(be$n_snps))),
  snp_info    = be$snp_info,
  blocks      = blocks,
  blues       = blues,
  id_col      = "id",      # name of the ID column in blues
  blue_col    = "YLD"      # name of the BLUE column (single trait)
)

close_backend(be)
```

For chromosome-specific analyses (e.g. a single chromosome only):

``` r
# Filter blocks and extract by chromosome — backend still streams
haps_chr1 <- extract_haplotypes(be, be$snp_info, blocks,
                                 chr = "1", min_snps = 5)
```

------------------------------------------------------------------------

## 9. Troubleshooting

**GDS handle errors on Windows.** SNPRelate GDS file handles are not
fork-safe. Avoid FORK-based parallelism
([`parallel::mclapply()`](https://rdrr.io/r/parallel/mcdummies.html));
use `n_threads > 1` (OpenMP within a single process) instead.

**BEDMatrix and missing data.** BEDMatrix codes missing genotypes as
`NA`. These are correctly handled by `compute_r2_cpp()` via mean
imputation per column. No pre-imputation is required.

**“fewer than min_snps_chr SNPs” messages.** Many reference genomes have
hundreds of unplaced scaffolds with a handful of markers. Set
`min_snps_chr = 100` or higher to skip them automatically.

**Out-of-memory on boundary scan.** If `boundary_scan_cpp()` runs out of
RAM on very large chromosomes, reduce `leng` to 30–50. This rarely
affects block quality because large-block boundaries are evident even
with a narrow scan window.

**Very long chromosomes (\> 2 M SNPs).** For panels with \> 2 M SNPs per
chromosome, set `subSegmSize = 500` and `leng = 30`. The algorithm will
make more GDS reads but use far less RAM per step.

**CLQD hangs or takes hours on a single chromosome.** This is
Bron-Kerbosch clique enumeration blowup. The
`[CLQD] Found N maximal cliques` message will show N in the millions.
Switch to `CLQmode = "Leiden"` — it runs in polynomial time and finishes
the same window in seconds. Also set `CLQcut = 0.70` and
`max_bp_distance = 500000L` to further reduce graph density.
`checkLargest = TRUE` adds a safety guard for the Density/Maximal modes
if you need exact clique enumeration.

------------------------------------------------------------------------

## 10. References

- Lawrence M, Huber W, Pages H, et al. (2013). Software for computing
  and annotating genomic ranges. *PLOS Computational Biology*
  **9**(8):e1003118. <https://doi.org/10.1371/journal.pcbi.1003118>
- VanRaden PM (2008). Efficient methods to compute genomic predictions.
  *Journal of Dairy Science* **91**(11):4414-4423.
  <https://doi.org/10.3168/jds.2007-0980>
