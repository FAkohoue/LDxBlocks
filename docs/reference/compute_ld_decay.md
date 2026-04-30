# Compute LD Decay and Chromosome-Specific Decay Distances

Calculates the decay of linkage disequilibrium (r\\^2\\ or rV\\^2\\)
with physical distance for each chromosome, fits optional decay models
(Hill-Weir nonlinear or LOESS), and determines the distance at which LD
drops below a user-specified or data-driven critical threshold. The
output can be passed directly to
[`define_qtl_regions()`](https://FAkohoue.github.io/LDxBlocks/reference/define_qtl_regions.md)
to define chromosome-specific candidate gene windows around GWAS hits.

## Usage

``` r
compute_ld_decay(
  geno,
  snp_info = NULL,
  method = c("r2", "rV2"),
  kin_method = c("chol", "eigen"),
  sampling = c("random", "sliding_window", "both"),
  window_snps = 50L,
  n_pairs = 50000L,
  max_dist = 5000000L,
  r2_threshold = "both",
  n_unlinked = 100000L,
  pctile = 95,
  fit_model = c("loess", "nonlinear", "both", "none"),
  n_bins = 100L,
  by_chr = TRUE,
  chr = NULL,
  seed = 42L,
  n_threads = 1L,
  verbose = TRUE
)
```

## Arguments

- geno:

  `LDxBlocks_backend` (from
  [`read_geno`](https://FAkohoue.github.io/LDxBlocks/reference/read_geno.md)
  or
  [`read_geno_bigmemory`](https://FAkohoue.github.io/LDxBlocks/reference/read_geno_bigmemory.md)),
  a numeric dosage matrix (individuals x SNPs), or a file path accepted
  by
  [`read_geno`](https://FAkohoue.github.io/LDxBlocks/reference/read_geno.md).

- snp_info:

  Data frame with columns `SNP`, `CHR`, `POS`. Auto-extracted from the
  backend when `geno` is a backend object.

- method:

  LD metric: `"r2"` (default) or `"rV2"` (kinship-adjusted; requires
  AGHmatrix and ASRgenomics).

- kin_method:

  Whitening method for `rV2`: `"chol"` (default) or `"eigen"`.

- sampling:

  Pair-sampling strategy: `"random"` (default, WGS-safe),
  `"sliding_window"` (TASSEL approach), or `"both"`.

- window_snps:

  Integer. Window size in SNPs for `"sliding_window"`. Default `50L`.

- n_pairs:

  Integer. Maximum SNP pairs per chromosome for `"random"` sampling.
  Default `50000L`.

- max_dist:

  Integer. Maximum physical distance (bp) between a pair. Default
  `5000000L` (5 Mb).

- r2_threshold:

  Threshold for decay distance: a numeric value (e.g. `0.1`),
  `"parametric"` (95th percentile of unlinked r\\^2\\), `"both"`, or
  `NULL`. Default `"both"`.

- n_unlinked:

  Integer. Cross-chromosome pairs for parametric threshold. Default
  `100000L`.

- pctile:

  Numeric (0-100). Percentile of unlinked r\\^2\\ distribution for the
  parametric threshold. Default `95`.

- fit_model:

  Decay model: `"loess"` (default), `"nonlinear"` (Hill-Weir), `"both"`,
  or `"none"`.

- n_bins:

  Integer. Distance bins for the smoothed decay curve. Default `100L`.

- by_chr:

  Logical. Return per-chromosome decay distances. Default `TRUE`.

- chr:

  Character vector of chromosomes to process. `NULL` = all.

- seed:

  Integer. RNG seed. Default `42L`.

- n_threads:

  Integer. OpenMP threads for C++ r\\^2\\ kernel. Default `1L`.

- verbose:

  Logical. Print timestamped progress. Default `TRUE`.

## Value

A named list of class `LDxBlocks_decay`:

- `pairs`:

  Data frame: `CHR`, `dist_bp`, `r2`.

- `decay_curve`:

  Binned decay curve per chromosome: `CHR`, `dist_bp`, `r2_mean`,
  `r2_median`, `n_pairs`, and optionally `r2_loess`, `r2_hw`.

- `decay_dist`:

  Per-chromosome decay distance: `CHR`, `decay_dist_bp`,
  `decay_dist_kb`, `threshold_used`, `r2_col_used`, `censored` (`TRUE`
  when the smoothed curve never drops below the threshold within
  `max_dist` – the returned distance is a lower bound, not an estimate;
  extend `max_dist` or lower the threshold).

- `decay_dist_genome`:

  Genome-wide median decay distance (bp).

- `critical_r2`:

  Active threshold value.

- `critical_r2_fixed`:

  Always set to 0.1 as a reference value (standard threshold used in
  most GWAS papers). This field is populated regardless of which
  `r2_threshold` was requested. It is the *active* threshold only when
  `r2_threshold = 0.1` or `r2_threshold = "both"` with no parametric
  result available.

- `critical_r2_param`:

  Parametric threshold (`NULL` if not computed).

- `unlinked_r2`:

  Cross-chromosome r\\^2\\ values used for the parametric threshold
  (`NULL` if not computed).

- `model_params`:

  Hill-Weir C parameter per chromosome (`NULL` if not fitted).

- `n_pairs_used`:

  Named integer: pairs computed per chromosome.

- `method`:

  LD metric used.

- `sampling`:

  Sampling strategy used.

- `call`:

  Matched function call.

## Details

**C++ acceleration:** All pairwise r\\^2\\ computations use the
package's compiled C++ kernels (Armadillo BLAS via RcppArmadillo):

- Random sampling: `compute_r2_sparse_cpp()` computes all pairs within
  `max_dist` in a single C++ call on the compact submatrix, replacing an
  R pair-by-pair loop. `n_threads` is forwarded to enable OpenMP
  parallelism.

- Sliding window: column preparation is done once per window via
  [`apply()`](https://rdrr.io/r/base/apply.html); within-window
  [`stats::cor()`](https://rdrr.io/r/stats/cor.html) is used (windows
  are small, typically 50 SNPs).

- Parametric threshold: `col_r2_cpp()` computes cross-chromosome r\\^2\\
  in Armadillo BLAS.

**Memory-efficient backend access:** For GDS and bigmemory backends the
function never loads a full chromosome. Instead it samples pair indices
from SNP positions (no genotype I/O), collects only the unique SNP
indices involved, and calls
[`read_chunk()`](https://FAkohoue.github.io/LDxBlocks/reference/read_chunk.md)
once for that compact set. For 50,000 random pairs on a 500,000-SNP
chromosome with 500 individuals, peak RAM is approximately 40 MB rather
than 2 GB.

**Important: `method = "rV2"` is not memory-safe at WGS scale.** The
kinship whitening matrix requires a full-genome read of all SNPs once
before per-chromosome processing begins. For panels with millions of
markers this can be the dominant memory step. Use `method = "r2"` for LD
decay estimation; reserve `method = "rV2"` for block detection.

**Sampling strategies:**

- `"sliding_window"`:

  Computes r\\^2\\ within a moving window of `window_snps` SNPs, reading
  only those columns per step. Uses an **adaptive stride**: the window
  is placed at approximately `n_pairs / (window_snps*(window_snps-1)/2)`
  evenly-spaced positions across the chromosome, so total pairs stays
  near `n_pairs` regardless of marker density. For `window_snps=50` on a
  500k-SNP chromosome with `n_pairs=50000`: stride=12,195, only 41
  `read_chunk` calls instead of 10,000. A **density guard**
  automatically switches chromosomes with \> 200,000 SNPs to `"random"`
  sampling (consecutive SNPs at WGS density are always in near- perfect
  LD, making the sliding-window curve uninformative). Matches the TASSEL
  approach (Bradbury et al. 2007); recommended for chip-density panels
  (\< 200k SNPs per chromosome).

- `"random"`:

  Randomly samples up to `n_pairs` SNP pairs per chromosome within
  `max_dist` bp, then loads only the unique SNPs involved. Bounded RAM
  regardless of panel density; recommended for WGS panels (\> 500k SNPs
  per chromosome).

- `"both"`:

  Sliding window for the decay curve shape; random sampling for the
  parametric threshold estimation.

**Critical r\\^2\\ threshold:**

- Fixed numeric (e.g. `0.1`):

  Standard value used in most GWAS papers. The distance where the decay
  curve crosses this value defines the LD block extent for candidate
  gene searches.

- `"parametric"`:

  95th percentile of r\\^2\\ between markers on *different* chromosomes
  (unlinked markers). In a structured population, unlinked r\\^2\\ \> 0
  due to kinship. This threshold captures the background LD attributable
  to structure, above which LD is genuine. A high parametric threshold
  (\> 0.05) is a strong indicator that `method = "rV2"` should be used
  for block detection.

- `"both"`:

  Returns both; uses the parametric as the active threshold for decay
  distance estimation.

**Decay model:** The Hill-Weir (1988) expectation relates expected
r\\^2\\ to distance \\d\\ (in Mb) via \\C = 4 N_e r d\\. LDxBlocks fits
this per chromosome with C as the free parameter via nonlinear least
squares. The LOESS smooth is a non-parametric alternative robust to
non-standard decay shapes.

**Connection to
[`define_qtl_regions()`](https://FAkohoue.github.io/LDxBlocks/reference/define_qtl_regions.md)
and GWAS integration:** The output of `compute_ld_decay()` can be passed
directly to
[`define_qtl_regions`](https://FAkohoue.github.io/LDxBlocks/reference/define_qtl_regions.md)`(ld_decay = decay)`
to replace fixed block boundaries with chromosome-specific LD-aware
windows.

- Without `ld_decay`:

  [`define_qtl_regions()`](https://FAkohoue.github.io/LDxBlocks/reference/define_qtl_regions.md)
  asks *which blocks CONTAIN a significant GWAS SNP?* The search window
  equals the block boundary. GWAS hits in gaps between blocks or near
  block edges are missed.

- With `ld_decay`:

  Asks *which blocks are IN LD with a significant GWAS SNP?* The search
  window is extended by the chromosome-specific decay distance on both
  sides. Adds columns `candidate_region_start`, `candidate_region_end`,
  and `candidate_region_size_kb` ready for BioMart/Ensembl Plants
  candidate gene queries.

This directly feeds
[`integrate_gwas_haplotypes()`](https://FAkohoue.github.io/LDxBlocks/reference/integrate_gwas_haplotypes.md):
the GWAS biological-evidence layer (layer 1 of the three-evidence score)
awards a point to blocks that *contain* a GWAS hit. With `ld_decay`,
blocks *near* a hit (within the decay distance) also receive the
evidence point, making the three-evidence scoring in
[`rank_haplotype_blocks()`](https://FAkohoue.github.io/LDxBlocks/reference/rank_haplotype_blocks.md)
more sensitive and biologically accurate. The recommended workflow is
therefore:

1.  `decay <- compute_ld_decay(be, ...)`

2.  `qtl <- define_qtl_regions(..., ld_decay = decay)`

3.  `ranked <- integrate_gwas_haplotypes(qtl, pred, ...)`

## See also

[`define_qtl_regions`](https://FAkohoue.github.io/LDxBlocks/reference/define_qtl_regions.md),
[`plot_ld_decay`](https://FAkohoue.github.io/LDxBlocks/reference/plot_ld_decay.md),
[`run_Big_LD_all_chr`](https://FAkohoue.github.io/LDxBlocks/reference/run_Big_LD_all_chr.md)

## Examples

``` r
# \donttest{
data(ldx_geno,     package = "LDxBlocks")
data(ldx_snp_info, package = "LDxBlocks")

# Random sampling -- memory-safe for any panel density
decay <- compute_ld_decay(
  geno         = ldx_geno,
  snp_info     = ldx_snp_info,
  sampling     = "random",
  r2_threshold = "both",
  fit_model    = "loess",
  n_pairs      = 5000L,
  verbose      = FALSE
)
decay$critical_r2_param    # background kinship-induced LD
#>        95% 
#> 0.03718883 
decay$critical_r2_fixed    # standard 0.1 threshold
#> [1] 0.1
decay$decay_dist           # per-chromosome decay distances
#>   CHR decay_dist_bp decay_dist_kb threshold_used r2_col_used censored
#> 1   1         39181         39.18     0.03718883    r2_loess    FALSE
#> 2   2         50902         50.90     0.03718883    r2_loess    FALSE
#> 3   3         29560         29.56     0.03718883    r2_loess    FALSE

# Sliding window (TASSEL approach) + Hill-Weir model
decay_hw <- compute_ld_decay(
  geno         = ldx_geno,
  snp_info     = ldx_snp_info,
  sampling     = "sliding_window",
  window_snps  = 50L,
  r2_threshold = 0.1,
  fit_model    = "nonlinear",
  verbose      = FALSE
)

# WGS backend -- only sampled columns are loaded per chromosome
if (FALSE) { # \dontrun{
be <- read_geno("my_wgs.vcf.gz")
decay_wgs <- compute_ld_decay(
  geno         = be,
  sampling     = "random",
  n_pairs      = 50000L,
  max_dist     = 5000000L,
  r2_threshold = "both",
  fit_model    = "loess",
  n_threads    = 8L
)
close_backend(be)
} # }

# Pass to define_qtl_regions for chromosome-specific candidate windows
data(ldx_blocks, package = "LDxBlocks")
data(ldx_gwas,   package = "LDxBlocks")
qtl <- define_qtl_regions(ldx_gwas, ldx_blocks, ldx_snp_info,
                          ld_decay    = decay,
                          p_threshold = NULL,
                          trait_col   = "trait")
qtl[, c("block_id", "lead_marker", "candidate_region_start",
         "candidate_region_end", "candidate_region_size_kb")]
#>                block_id lead_marker candidate_region_start candidate_region_end
#> 1    block_1_1000_25027      rs1005                      0                44369
#> 2   block_1_81064_99022      rs1048                  58031               136393
#> 3 block_1_155368_179371      rs1070                 129829               208191
#> 4    block_2_1000_30023      rs2004                      0                55303
#> 5  block_2_86236_105290      rs2050                  49548               151352
#> 6    block_3_1000_19068      rs3004                      0                33383
#>   candidate_region_size_kb
#> 1                     78.4
#> 2                     78.4
#> 3                     78.4
#> 4                    101.8
#> 5                    101.8
#> 6                     59.1
# }
```
