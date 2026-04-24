# End-to-End Haplotype Block Pipeline

Single-call wrapper: one genotype file in, complete haplotype dataset
out. Handles format detection, MAF filtering, call-rate filtering,
imputation, LD block detection, optional Beagle phasing, haplotype
extraction, diversity analysis, feature matrix construction, and output
writing.

## Usage

``` r
run_ldx_pipeline(
  geno_source,
  out_dir = ".",
  out_blocks,
  out_diversity,
  out_hap_matrix,
  hap_format = c("numeric", "character"),
  phase = FALSE,
  beagle_jar = NULL,
  beagle_threads = 1L,
  java_path = "java",
  beagle_java_mem_gb = NULL,
  beagle_args = "",
  beagle_ref_panel = NULL,
  beagle_map_file = NULL,
  beagle_chrom = NULL,
  beagle_seed = NULL,
  beagle_burnin = NULL,
  beagle_iterations = NULL,
  beagle_window = NULL,
  beagle_overlap = NULL,
  maf_cut = 0.05,
  impute = c("mean_rounded", "mode", "none"),
  min_callrate = 0,
  CLQcut = 0.5,
  method = c("r2", "rV2"),
  kin_method = "chol",
  CLQmode = c("Density", "Maximal", "Louvain", "Leiden"),
  leng = 200L,
  subSegmSize = 1500L,
  clstgap = 40000L,
  split = FALSE,
  appendrare = FALSE,
  singleton_as_block = FALSE,
  checkLargest = FALSE,
  digits = -1L,
  n_threads = 1L,
  min_snps_chr = 10L,
  max_bp_distance = 0L,
  min_snps_block = 3L,
  top_n = NULL,
  min_freq = 0.01,
  scale_hap_matrix = FALSE,
  chr = NULL,
  clean_malformed = FALSE,
  use_bigmemory = FALSE,
  bigmemory_path = tempdir(),
  bigmemory_type = "char",
  verbose = TRUE
)
```

## Arguments

- geno_source:

  File path or `LDxBlocks_backend`. Supported formats: VCF
  (`.vcf`/`.vcf.gz`), HapMap (`.hmp.txt`), CSV, GDS, PLINK BED.
  `phase = TRUE` requires VCF/VCF.gz.

- out_dir:

  Output directory. Default `"."`. When `phase = TRUE`, `beagle.jar` is
  expected here and the phased VCF is written here.

- out_blocks:

  Path for LD block table CSV.

- out_diversity:

  Path for haplotype diversity table CSV.

- out_hap_matrix:

  Path for haplotype genotype matrix file.

- hap_format:

  `"numeric"` (default) or `"character"`.

- phase:

  If `TRUE`, phase with Beagle. Default `FALSE`.

- beagle_jar:

  Path to `beagle.jar`. Default: `file.path(out_dir, "beagle.jar")`.

- beagle_threads:

  Beagle threads. Default `1L`.

- java_path:

  Java executable. Default `"java"`.

- beagle_java_mem_gb:

  JVM heap in GB (`-Xmx`). Default `NULL`.

- beagle_args:

  Extra Beagle argument string. Default `""`.

- beagle_ref_panel:

  Phased reference VCF path. Default `NULL`.

- beagle_map_file:

  Genetic map file path. Default `NULL`.

- beagle_chrom:

  Restrict Beagle to one chromosome. Default `NULL` (inherits `chr` if
  set).

- beagle_seed:

  Integer seed for reproducibility. Default `NULL`.

- beagle_burnin:

  Burn-in iterations. Default `NULL`.

- beagle_iterations:

  Phasing iterations. Default `NULL`.

- beagle_window:

  Window size. Default `NULL`.

- beagle_overlap:

  Window overlap. Default `NULL`.

- maf_cut:

  Minimum MAF. Default `0.05`.

- impute:

  `"mean_rounded"` (default), `"mode"`, or `"none"`.

- min_callrate:

  Minimum per-SNP call rate. Default `0.0`.

- CLQcut:

  r\\^2\\ threshold for CLQD. Default `0.5`.

- method:

  `"r2"` (default) or `"rV2"`.

- kin_method:

  `"chol"` (default) or `"eigen"`.

- CLQmode:

  `"Density"` (default), `"Maximal"`, `"Louvain"`, or `"Leiden"`.

- leng:

  Boundary scan half-window (SNPs). Default `200L`.

- subSegmSize:

  Max SNPs per sub-segment. Default `1500L`.

- clstgap:

  Max bp gap within clique. Default `40000L`.

- split:

  Split cliques at largest gap. Default `FALSE`.

- appendrare:

  Append rare SNPs to nearest block. Default `FALSE`.

- singleton_as_block:

  Return singletons as blocks. Default `FALSE`.

- checkLargest:

  Dense-core pre-pass. Default `FALSE`.

- digits:

  Round r\\^2\\ (`-1L` = off). Default `-1L`.

- n_threads:

  OpenMP threads. Default `1L`.

- min_snps_chr:

  Skip chromosomes below this SNP count. Default `10L`.

- max_bp_distance:

  Max bp for r\\^2\\ (`0L` = all). Default `0L`.

- min_snps_block:

  Min SNPs per haplotype block. Default `3L`.

- top_n:

  Max alleles per block (`NULL` = all above `min_freq`). Default `NULL`.

- min_freq:

  Min haplotype allele frequency. Default `0.01`.

- scale_hap_matrix:

  Scale haplotype matrix columns. Default `FALSE`.

- chr:

  Chromosomes to process (`NULL` = all). Default `NULL`.

- clean_malformed:

  Remove malformed VCF lines. Default `FALSE`.

- use_bigmemory:

  File-backed bigmemory store. Default `FALSE`.

- bigmemory_path:

  Directory for backing files. Default
  [`tempdir()`](https://rdrr.io/r/base/tempfile.html).

- bigmemory_type:

  `"char"` (default), `"short"`, or `"double"`.

- verbose:

  Print progress. Default `TRUE`.

## Value

Named list (invisibly): `blocks`, `diversity`, `hap_matrix`,
`hap_matrix_info`, `haplotypes`, `geno_matrix`, `snp_info_filtered`,
`phased_vcf`, `phased_backend_desc`, `phase_method`, `n_blocks`,
`n_hap_columns`.

## Phasing modes

- `phase = FALSE` (default):

  Dosage-pattern haplotypes extracted directly from the imputed matrix.
  No external tools required. Fast. Suitable for genomic prediction.
  Each block entry is a multi-SNP dosage string - not a true gametic
  haplotype. Frequencies are individual-level pattern proportions.

- `phase = TRUE`:

  Beagle 5.x called on the original input VCF after LD block detection,
  producing true statistically-inferred gametic haplotypes using
  population-LD across all markers. Haplotype strings become `"g1|g2"`.
  Frequencies are computed over \\2N\\ gamete observations. Recommended
  for diversity analysis and biologically interpretable results.

  **Requirements:** `geno_source` must be VCF/VCF.gz. Place `beagle.jar`
  in `out_dir` or supply via `beagle_jar`. Download:
  <https://faculty.washington.edu/browning/beagle/beagle.html>

## Why Beagle and why a cleaned VCF

Beagle 5.x performs chromosome-level statistical phasing using
population LD across all markers simultaneously, producing true inferred
gametic haplotypes. This is the only supported phasing method in
LDxBlocks.

When `phase = TRUE`, the pipeline does **not** pass `geno_source`
directly to Beagle. Instead it first writes a cleaned VCF
(`<geno_source_stem>_cleaned.vcf.gz` in `out_dir`) containing exactly
the SNPs that survived MAF filtering, call-rate filtering, chromosome
subsetting, and malformed-record removal. This guarantees that the
phased VCF's SNP set is identical to the marker set used for LD block
detection, so `.read_and_cache_phased_vcf()` can align phased gametes to
blocks without SNP-count mismatches.

## See also

[`run_Big_LD_all_chr`](https://FAkohoue.github.io/LDxBlocks/reference/run_Big_LD_all_chr.md),
[`extract_haplotypes`](https://FAkohoue.github.io/LDxBlocks/reference/extract_haplotypes.md),
[`build_haplotype_feature_matrix`](https://FAkohoue.github.io/LDxBlocks/reference/build_haplotype_feature_matrix.md),
[`run_haplotype_prediction`](https://FAkohoue.github.io/LDxBlocks/reference/run_haplotype_prediction.md)

## Examples

``` r
# \donttest{
geno_file <- system.file("extdata", "example_genotypes_numeric.csv",
                          package = "LDxBlocks")
res <- run_ldx_pipeline(
  geno_source = geno_file, out_dir = tempdir(),
  out_blocks = tempfile(fileext=".csv"),
  out_diversity = tempfile(fileext=".csv"),
  out_hap_matrix = tempfile(fileext=".csv"),
  phase = FALSE, maf_cut = 0.05, CLQcut = 0.5,
  leng = 10L, subSegmSize = 80L, verbose = FALSE
)
#> [LDxBlocks] Hap QC: n_snps range [19, 30] | blocks=90 | NA_matrix=0
#> [LDxBlocks] Pipeline QC: all checks passed.
if (FALSE) { # \dontrun{
# With Beagle phasing (place beagle.jar in out_dir first):
res2 <- run_ldx_pipeline(
  geno_source = "data.vcf.gz", out_dir = "results/",
  out_blocks = "results/blocks.csv",
  out_diversity = "results/diversity.csv",
  out_hap_matrix = "results/hap_matrix.csv",
  phase = TRUE, beagle_threads = 4L,
  beagle_java_mem_gb = 8, beagle_seed = 42L
)
} # }# }
```
