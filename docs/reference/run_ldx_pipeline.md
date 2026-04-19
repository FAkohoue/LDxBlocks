# End-to-End Haplotype Block Pipeline

A single-call wrapper that takes one genotype file and produces a
complete haplotype-based dataset ready for genomic prediction.
Internally handles format detection, optional GDS conversion for large
files, MAF filtering, genome-wide LD block detection, haplotype
extraction, diversity analysis, and output writing.

The user provides only the input file path and desired output paths. All
intermediate steps run transparently.

## Usage

``` r
run_ldx_pipeline(
  geno_source,
  out_blocks,
  out_diversity,
  out_hap_matrix,
  hap_format = c("numeric", "character"),
  maf_cut = 0.05,
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
  min_snps_block = 3L,
  top_n = NULL,
  min_freq = 0.01,
  scale_hap_matrix = FALSE,
  chr = NULL,
  verbose = TRUE,
  max_bp_distance = 0L,
  clean_malformed = FALSE,
  use_bigmemory = FALSE,
  bigmemory_path = tempdir(),
  bigmemory_type = "char"
)
```

## Arguments

- geno_source:

  Path to a genotype file, OR an `LDxBlocks_backend` object already
  created by
  [`read_geno`](https://FAkohoue.github.io/LDxBlocks/reference/read_geno.md)
  or
  [`read_geno_bigmemory`](https://FAkohoue.github.io/LDxBlocks/reference/read_geno_bigmemory.md).
  Passing a backend skips the internal
  [`read_geno()`](https://FAkohoue.github.io/LDxBlocks/reference/read_geno.md)
  call, allowing you to use any pre-built backend including file-backed
  `bigmemory` stores. Supported file formats when a path is supplied:
  numeric dosage CSV (\`.csv\`), HapMap (\`.hmp.txt\`), VCF (\`.vcf\`,
  \`.vcf.gz\`), SNPRelate GDS (\`.gds\`), PLINK BED (\`.bed\`). Format
  is detected automatically from the file extension.

- out_blocks:

  Path for the LD block table CSV. Columns: \`CHR\`, \`start\`, \`end\`,
  \`start.rsID\`, \`end.rsID\`, \`start.bp\`, \`end.bp\`, \`length_bp\`,
  \`length_snps\`, \`block_name\`.

- out_diversity:

  Path for the haplotype diversity table CSV. Columns: `block_id`,
  `CHR`, `start_bp`, `end_bp`, `n_snps`, `n_ind`, `n_haplotypes`, `He`
  (Nei 1973 corrected), `Shannon`, `n_eff_alleles` (effective number of
  alleles = 1/\\\sum p_i^2\\), `freq_dominant`, `sweep_flag` (`TRUE`
  when freq_dominant \>= 0.90, indicating a possible selective sweep or
  strong founder effect), `phased`.

- out_hap_matrix:

  Path for the haplotype genotype matrix. Format is controlled by
  \`hap_format\`. See section \*\*Haplotype genotype matrix\*\*.

- hap_format:

  Output format for the haplotype matrix:

  - `"numeric"` (default) - Tab-delimited. Rows = haplotype alleles,
    columns = individuals. Metadata columns: `hap_id`, `CHR`,
    `start_bp`, `end_bp`, `n_snps`, `alleles` (nucleotide sequence of
    this allele), `frequency`. Individual columns contain 0/1/2/NA
    dosage.

  - `"character"` - Tab-delimited. Same row/column orientation.
    Individual cells contain the nucleotide sequence of the haplotype
    allele if the individual carries it, `"-"` if absent, `"."` if
    missing.

- maf_cut:

  Minimum minor allele frequency. SNPs below this threshold are removed
  before block detection. Default \`0.05\`.

- CLQcut:

  r^2 threshold for clique edges in CLQD. Higher values produce tighter,
  smaller blocks. Default \`0.5\`.

- method:

  LD metric: `"r2"` (default) or `"rV2"` (requires kinship matrix; see
  [`run_Big_LD_all_chr`](https://FAkohoue.github.io/LDxBlocks/reference/run_Big_LD_all_chr.md)).

- kin_method:

  Whitening method for `rV2`: `"chol"` (Cholesky, default, faster) or
  `"eigen"` (eigendecomposition, more stable for near-singular kinship
  matrices).

- CLQmode:

  Clique scoring mode: `"Density"` (default, prefers compact
  high-density cliques – recommended for most analyses) or `"Maximal"`
  (prefers the largest cliques regardless of span).

- leng:

  Boundary scan half-window in SNPs. Default \`200L\`. Reduce to 50-100
  for very dense WGS panels.

- subSegmSize:

  Maximum SNPs per CLQD sub-segment. Controls peak RAM: \`subSegmSize x
  n_individuals x 8\` bytes. Default \`1500L\`.

- clstgap:

  Maximum base-pair gap allowed within a clique before it is split.
  Default `40000L` (40 kb). Increase for populations with long-range LD
  (e.g. inbred lines); decrease for high-recombination panels.

- split:

  Logical. If `TRUE`, split cliques whose SNP span exceeds `clstgap` bp
  at the largest internal gap. Default `FALSE`.

- appendrare:

  Logical. If `TRUE`, SNPs that fail MAF filtering are appended to the
  nearest block after detection. Default `FALSE`.

- singleton_as_block:

  Logical. If `TRUE`, SNPs that receive no clique assignment (singletons
  at recombination hotspots) are returned as single-SNP blocks in the
  block table. Default `FALSE` (original Big-LD behaviour – singletons
  are silently dropped).

- checkLargest:

  Logical. If `TRUE`, apply a dense-core pre-pass before clique
  enumeration on sub-segments with \>=500 SNPs to prevent exponential
  blowup. Default `FALSE`.

- digits:

  Integer. Round r^2 values to this many decimal places before clique
  detection. `-1L` (default) disables rounding.

- n_threads:

  OpenMP threads for the C++ LD kernel. Default \`1L\`.

- min_snps_chr:

  Skip chromosomes with fewer post-filter SNPs than this. Default
  \`10L\`. Increase to skip unplaced scaffolds.

- min_snps_block:

  Minimum SNPs per haplotype block. Blocks smaller than this are
  excluded from haplotype analysis. Default \`3L\`.

- top_n:

  Integer or `NULL`. Maximum haplotype alleles per block in the output
  matrix. `NULL` (default) retains all alleles above `min_freq` – no
  cap. Set an integer (e.g. `5L`) only to limit column count for memory
  reasons.

- min_freq:

  Minimum haplotype allele frequency (0–1). Alleles observed at lower
  frequency than this threshold are dropped before building the feature
  matrix and output files. Default `0.01` (1%). Rare alleles below this
  threshold cannot be estimated reliably in typical training sets and
  add noise. Lower values retain more rare alleles; higher values (e.g.
  `0.05`) match the MAF filter applied to SNPs.

- scale_hap_matrix:

  Logical. If \`TRUE\`, scale the haplotype matrix columns to zero mean
  and unit variance before writing. Useful for GBLUP-style models.
  Default \`FALSE\`.

- chr:

  Character vector of chromosome names to process. \`NULL\` (default)
  processes all chromosomes.

- verbose:

  Logical. Print timestamped progress. Default \`TRUE\`.

- max_bp_distance:

  Integer. Maximum bp distance between a SNP pair for its r\\^2\\ to be
  computed in the LD graph. `0L` (default) computes all pairs.
  Recommended for WGS panels: `500000L` (500 kb). Reduces O(p\\^2\\) LD
  computation to near-O(p) when set.

- clean_malformed:

  Logical. Stream-clean the input file before reading by removing any
  lines whose column count does not match the header. Needed for files
  from NGSEP and some variant callers. Default `FALSE`.

- use_bigmemory:

  Logical. If `TRUE` and `geno_source` is a file path, the source file
  is first loaded into a file-backed
  [`bigmemory::big.matrix`](https://rdrr.io/pkg/bigmemory/man/big.matrix.html)
  before block detection. Only the OS pages needed for each sub-segment
  window are loaded into RAM, keeping peak memory proportional to
  `n_samples x subSegmSize` rather than the full genome. Requires the
  bigmemory package. Default `FALSE`. Ignored when `geno_source` is
  already an `LDxBlocks_backend` (the supplied backend is used as-is).

- bigmemory_path:

  Directory where the bigmemory backing files (`.bin` and `.desc`) and
  SNP info cache (`_snpinfo.rds`) are written or read from. Defaults to
  [`tempdir()`](https://rdrr.io/r/base/tempfile.html) so files are
  cleaned up at the end of the R session. Supply a persistent directory
  (e.g. next to the output CSVs) when you want to reattach the backing
  file on a restart without re-reading the VCF.

- bigmemory_type:

  Storage type for the `big.matrix`: `"char"` (1 byte per cell, 8x
  smaller than double, sufficient for 0/1/2 dosage), `"short"` (2
  bytes), or `"double"` (8 bytes). Default `"char"`.

## Value

A named list (invisibly) with elements:

- `blocks`:

  Data frame of LD blocks from `run_Big_LD_all_chr`.

- `diversity`:

  Data frame of per-block haplotype diversity metrics: `block_id`,
  `CHR`, `start_bp`, `end_bp`, `n_snps`, `n_ind`, `n_haplotypes`, `He`
  (sample-size corrected expected heterozygosity), `Shannon`,
  `n_eff_alleles`, `freq_dominant`, `sweep_flag`, `phased`.

- `hap_matrix`:

  Numeric matrix (individuals x haplotype allele columns) – the
  dimensionality-reduced genotype matrix for genomic prediction. Always
  returned as a numeric R matrix regardless of `hap_format` (which only
  controls the *file* written to `out_hap_matrix`). Dosage values:
  0/1/2/NA for phased data; 0/1/NA for unphased data. Compatible with
  rrBLUP, BGLR, sommer, ASReml-R.

- `haplotypes`:

  Named list of per-block haplotype dosage strings from
  [`extract_haplotypes()`](https://FAkohoue.github.io/LDxBlocks/reference/extract_haplotypes.md).
  Pass to
  [`decode_haplotype_strings()`](https://FAkohoue.github.io/LDxBlocks/reference/decode_haplotype_strings.md)
  for nucleotide sequences, or to
  [`rank_haplotype_blocks()`](https://FAkohoue.github.io/LDxBlocks/reference/rank_haplotype_blocks.md)
  for evidence-based ranking.

- `snp_info_filtered`:

  Data frame of SNP metadata after MAF filtering: `SNP`, `CHR`, `POS`,
  and any additional columns from the input file.

- `geno_matrix`:

  Numeric matrix (individuals x SNPs) of MAF-filtered genotypes
  (0/1/2/NA). Needed directly by
  [`tune_LD_params`](https://FAkohoue.github.io/LDxBlocks/reference/tune_LD_params.md)
  and
  [`run_haplotype_prediction`](https://FAkohoue.github.io/LDxBlocks/reference/run_haplotype_prediction.md)
  – avoids reloading the genotype file after the pipeline completes.

- `n_blocks`:

  Integer. Total LD blocks detected genome-wide.

- `n_hap_columns`:

  Integer. Total haplotype allele columns after `min_freq` filtering –
  the effective number of predictors for genomic prediction.

## File-backed memory-mapped genotype store (`use_bigmemory`)

When `use_bigmemory = TRUE` the pipeline converts the source file into a
[`bigmemory::big.matrix`](https://rdrr.io/pkg/bigmemory/man/big.matrix.html)
backed by two binary files on disk before running any analysis. Only the
OS pages needed for each sub-segment window are loaded into RAM; the
rest of the genome stays on disk. Peak RAM is proportional to
`n_samples x subSegmSize x bytes_per_cell` rather than the full genome
matrix.

**Backing files created in `bigmemory_path`:**

- `ldxblocks_bm.bin`:

  Raw binary genotype data. Size = n_samples x n_snps x bytes_per_cell.
  With `type = "char"` and 204 samples x 3M SNPs: ~0.6 GB.

- `ldxblocks_bm.desc`:

  Tiny text descriptor that bigmemory uses to memory-map the `.bin`
  file. Never edit this manually.

- `ldxblocks_bm_snpinfo.rds`:

  Cached SNP metadata (SNP, CHR, POS, REF, ALT). bigmemory does not
  store metadata, so this file is saved separately to enable restart
  without re-reading the source file.

**Restart behaviour:** if all three files already exist in
`bigmemory_path` when the pipeline is called, the `.bin` is reattached
instantly via
[`bigmemory::attach.big.matrix()`](https://rdrr.io/pkg/bigmemory/man/attach.big.matrix.html)
– the source VCF is not touched. To force a rebuild, delete the `.bin`
and `.desc` files.

**Choosing `bigmemory_type`:**

- `"char"` (default, recommended):

  1 signed byte per cell (range –128..127). Genotype dosage values 0, 1,
  2 fit without loss. 8x smaller than `"double"`.

- `"short"`:

  2 bytes per cell. Not needed for 0/1/2 dosage; use only if your
  pipeline stores values outside –128..127 in the same matrix.

- `"double"`:

  8 bytes per cell. Same as a standard R
  [`matrix()`](https://rdrr.io/r/base/matrix.html). Use only if
  downstream code requires `double` precision and the RAM saving is not
  needed.

**When to use `use_bigmemory = TRUE`:**

- Peak RAM from GDS streaming exceeds available node memory.

- You need fast restart: the `.bin` persists across R sessions (supply a
  persistent `bigmemory_path`, not
  [`tempdir()`](https://rdrr.io/r/base/tempfile.html)).

- Multiple pipeline runs share the same panel (one `.bin`, many readers
  – bigmemory memory-maps are safe for concurrent access).

## Scale behaviour

Files are processed via the \`LDxBlocks_backend\` streaming interface.
For VCF, HapMap, and numeric CSV files the backend reads one chromosome
window at a time - peak RAM equals one \`subSegmSize\`-SNP window
regardless of total marker count. For very large files (\> 2 M SNPs) the
GDS backend via \`SNPRelate\` is used automatically if the \`SNPRelate\`
package is installed.

## Haplotype genotype matrix

Both output formats have haplotype alleles as rows and individuals as
columns, preceded by metadata columns: `hap_id`, `CHR`, `start_bp`,
`end_bp`, `n_snps`, `alleles` (nucleotide sequence of this allele),
`frequency`.

- `"numeric"`: individual cells are haplotype dosage values:

  - **Phased data**: 0/1/2/NA – 0 = neither gamete carries this allele,
    1 = one gamete carries it (heterozygous), 2 = both gametes carry it
    (homozygous).

  - **Unphased data**: 0/1/NA – 0 = absent, 1 = present (individual's
    block-level string matches this allele exactly). The value 2 is not
    used for unphased data because the two chromosomes cannot be
    distinguished – an individual homozygous for this allele and one
    heterozygous for it produce different observable strings and are
    treated as different alleles.

  Compatible with rrBLUP, BGLR, sommer, ASReml-R.

- `"character"`: individual cells are the full nucleotide sequence of
  the allele (e.g. `"AGTTA"`) if the individual carries it, `"-"` if
  absent, `"."` if missing. Heterozygous positions use IUPAC ambiguity
  codes (R=A/G, Y=C/T, S=G/C, W=A/T, K=G/T, M=A/C).

## See also

[`run_Big_LD_all_chr`](https://FAkohoue.github.io/LDxBlocks/reference/run_Big_LD_all_chr.md),
[`extract_haplotypes`](https://FAkohoue.github.io/LDxBlocks/reference/extract_haplotypes.md),
[`compute_haplotype_diversity`](https://FAkohoue.github.io/LDxBlocks/reference/compute_haplotype_diversity.md),
[`build_haplotype_feature_matrix`](https://FAkohoue.github.io/LDxBlocks/reference/build_haplotype_feature_matrix.md),
[`rank_haplotype_blocks`](https://FAkohoue.github.io/LDxBlocks/reference/rank_haplotype_blocks.md),
[`run_haplotype_prediction`](https://FAkohoue.github.io/LDxBlocks/reference/run_haplotype_prediction.md),
[`tune_LD_params`](https://FAkohoue.github.io/LDxBlocks/reference/tune_LD_params.md)

## Examples

``` r
# \donttest{
# Path A: supply a file path (original interface, unchanged)
geno_file <- system.file("extdata", "example_genotypes_numeric.csv",
                         package = "LDxBlocks")

res <- run_ldx_pipeline(
  geno_source    = geno_file,
  out_blocks     = tempfile(fileext = ".csv"),
  out_diversity  = tempfile(fileext = ".csv"),
  out_hap_matrix = tempfile(fileext = ".csv"),
  hap_format     = "numeric",
  maf_cut        = 0.05,
  CLQcut         = 0.5,
  leng           = 10L,
  subSegmSize    = 80L,
  min_snps_block = 3L,
  # top_n       = 5L,   # optional integer cap; NULL = all above min_freq
  # hap_format  = "character",  # alternative to "numeric"
  verbose        = FALSE
)

head(res$blocks)
#>   start end start.rsID end.rsID start.bp end.bp CHR length_bp
#> 1     1  25     rs1001   rs1025     1000  25027   1     24028
#> 2    31  50     rs1031   rs1050    81064  99022   1     17959
#> 3    56  80     rs1056   rs1080   155368 179371   1     24004
#> 4     1  30     rs2001   rs2030     1000  30023   2     29024
#> 5    36  55     rs2036   rs2055    86236 105290   2     19055
#> 6    61  80     rs2061   rs2080   161515 180473   2     18959
head(res$diversity)
#>                block_id CHR start_bp end_bp n_snps n_ind n_haplotypes        He
#> 1    block_1_1000_25027   1     1000  25027     25   120           30 0.9163866
#> 2   block_1_81064_99022   1    81064  99022     20   120           21 0.8988796
#> 3 block_1_155368_179371   1   155368 179371     25   120           26 0.9159664
#> 4    block_2_1000_30023   2     1000  30023     30   120           31 0.9161064
#> 5  block_2_86236_105290   2    86236 105290     20   120           22 0.8955182
#> 6 block_2_161515_180473   2   161515 180473     20   120           25 0.9155462
#>    Shannon n_eff_alleles freq_dominant sweep_flag phased
#> 1 2.745726        10.959     0.1833333      FALSE  FALSE
#> 2 2.476482         9.207     0.1833333      FALSE  FALSE
#> 3 2.667890        10.909     0.1666667      FALSE  FALSE
#> 4 2.772723        10.926     0.2000000      FALSE  FALSE
#> 5 2.485949         8.933     0.2166667      FALSE  FALSE
#> 6 2.642908        10.860     0.1500000      FALSE  FALSE
dim(res$hap_matrix)
#> [1] 120  90

# Path B: supply a pre-built bigmemory backend
if (requireNamespace("bigmemory", quietly = TRUE)) {
  # read_geno_bigmemory() accepts a file path directly:
  # it calls read_geno() internally and wraps the result.
  be_bm <- read_geno_bigmemory(
    source      = geno_file,
    backingfile = tempfile("ldxbm"),
    type        = "char"
  )
  res2 <- run_ldx_pipeline(
    geno_source    = be_bm,
    out_blocks     = tempfile(fileext = ".csv"),
    out_diversity  = tempfile(fileext = ".csv"),
    out_hap_matrix = tempfile(fileext = ".csv"),
    CLQcut         = 0.5, leng = 10L, subSegmSize = 70L
  )
  close_backend(be_bm)
}
#> [bigmemory] Opening 'example_genotypes_numeric.csv' via read_geno() ...
#> [read_geno] Reading numeric dosage (chunked): C:/Users/fakohoue/AppData/Local/R/win-library/4.5/LDxBlocks/extdata/example_genotypes_numeric.csv 
#> [read_geno]   230 SNPs x 120 samples
#> [bigmemory] Allocating 120 x 230 big.matrix (type = 'char') ...
#> [bigmemory]   chr 1 loaded (80/230 SNPs)
#> [bigmemory]   chr 2 loaded (160/230 SNPs)
#> [bigmemory]   chr 3 loaded (230/230 SNPs)
#> [bigmemory] Done. Backing file: C:\Users\fakohoue\AppData\Local\Temp\RtmpIFdZjq/ldxbm74a863ad255d.bin
#> [23:24:29] Using pre-built backend: bigmemory | 120 ind | 230 SNPs
#> [23:24:29] MAF filtering (>= 0.05) ...
#> [MAF filter] Computing MAF for 230 SNPs ...
#> [MAF filter] 230 / 230 SNPs pass MAF >= 0.05
#> [23:24:29] Loading filtered genotype matrix ...
#> [23:24:29] Genotype matrix: 120 x 230
#> [23:24:29] Running genome-wide LD block detection ...
#> 
#> [run_Big_LD_all_chr] Processing 1 ...
#> [Big_LD] Subsegmenting via C++ boundary scan...
#> [Big_LD]   Cut at SNP 24
#> [Big_LD]   Cut at SNP 29
#> [Big_LD]   Cut at SNP 50
#> [Big_LD]   Cut at SNP 55
#> [CLQD] Building graph on 24 SNPs...
#> [CLQD] Found 6 maximal cliques.
#> [Big_LD] Segment 1/5 done.
#> [CLQD] Building graph on 5 SNPs...
#> [CLQD] Found 0 maximal cliques.
#> [Big_LD] Segment 2/5 done.
#> [CLQD] Building graph on 21 SNPs...
#> [CLQD] Found 6 maximal cliques.
#> [Big_LD] Segment 3/5 done.
#> [CLQD] Building graph on 5 SNPs...
#> [CLQD] Found 0 maximal cliques.
#> [Big_LD] Segment 4/5 done.
#> [CLQD] Building graph on 25 SNPs...
#> [CLQD] Found 6 maximal cliques.
#> [Big_LD] Segment 5/5 done.
#> 
#> [run_Big_LD_all_chr] Processing 2 ...
#> [Big_LD] Subsegmenting via C++ boundary scan...
#> [Big_LD]   Cut at SNP 30
#> [Big_LD]   Cut at SNP 35
#> [Big_LD]   Cut at SNP 55
#> [Big_LD]   Cut at SNP 60
#> [CLQD] Building graph on 30 SNPs...
#> [CLQD] Found 6 maximal cliques.
#> [Big_LD] Segment 1/5 done.
#> [CLQD] Building graph on 5 SNPs...
#> [CLQD] Found 0 maximal cliques.
#> [Big_LD] Segment 2/5 done.
#> [CLQD] Building graph on 20 SNPs...
#> [CLQD] Found 5 maximal cliques.
#> [Big_LD] Segment 3/5 done.
#> [CLQD] Building graph on 5 SNPs...
#> [CLQD] Found 0 maximal cliques.
#> [Big_LD] Segment 4/5 done.
#> [CLQD] Building graph on 20 SNPs...
#> [CLQD] Found 5 maximal cliques.
#> [Big_LD] Segment 5/5 done.
#> 
#> [run_Big_LD_all_chr] Processing 3 ...
#> [CLQD] Building graph on 70 SNPs...
#> [CLQD] Found 16 maximal cliques.
#> [Big_LD] Segment 1/1 done.
#> [23:24:30] Detected 9 LD blocks
#> [23:24:30] Block table written: C:\Users\fakohoue\AppData\Local\Temp\RtmpIFdZjq\file74a81bfc267b.csv
#> [23:24:30] Extracting haplotypes (min_snps = 3) ...
#> [23:24:30] Haplotypes extracted for 9 blocks
#> [23:24:30] Computing haplotype diversity ...
#> [23:24:30] Diversity table written: C:\Users\fakohoue\AppData\Local\Temp\RtmpIFdZjq\file74a84268e17.csv
#> [23:24:30] Building haplotype feature matrix (top_n = , min_freq = 0.01) ...
#> [23:24:30] Haplotype matrix: 120 individuals x 90 haplotype allele columns
#> [23:24:30] Writing haplotype matrix (format = numeric) ...
#> [write_haplotype_numeric] C:\Users\fakohoue\AppData\Local\Temp\RtmpIFdZjq\file74a817981499.csv (90 haplotypes x 120 individuals)
#> [23:24:30] Haplotype matrix written: C:\Users\fakohoue\AppData\Local\Temp\RtmpIFdZjq\file74a817981499.csv
#> [23:24:30] Pipeline complete.
#> [23:24:30]   Blocks:              9
#> [23:24:30]   Haplotype blocks:    9
#> [23:24:30]   Haplotype columns:   90
#> [23:24:30]   Individuals:         120
# }
```
