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
  geno_file,
  out_blocks,
  out_diversity,
  out_hap_matrix,
  hap_format = c("numeric", "character"),
  maf_cut = 0.05,
  CLQcut = 0.5,
  method = c("r2", "rV2"),
  kin_method = "chol",
  CLQmode = "Density",
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
  clean_malformed = FALSE
)
```

## Arguments

- geno_file:

  Path to genotype file. Supported formats: numeric dosage CSV
  (\`.csv\`), HapMap (\`.hmp.txt\`), VCF (\`.vcf\`, \`.vcf.gz\`),
  SNPRelate GDS (\`.gds\`), PLINK BED (\`.bed\`). Format is detected
  automatically from the file extension.

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

- clean_malformed:

  Logical. Stream-clean the input file before reading by removing any
  lines whose column count does not match the header. Needed for files
  from NGSEP and some variant callers. Default `FALSE`.

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

  Numeric matrix (individuals x haplotype allele columns) — the
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
  — avoids reloading the genotype file after the pipeline completes.

- `n_blocks`:

  Integer. Total LD blocks detected genome-wide.

- `n_hap_columns`:

  Integer. Total haplotype allele columns after `min_freq` filtering —
  the effective number of predictors for genomic prediction.

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

  - **Phased data**: 0/1/2/NA — 0 = neither gamete carries this allele,
    1 = one gamete carries it (heterozygous), 2 = both gametes carry it
    (homozygous).

  - **Unphased data**: 0/1/NA — 0 = absent, 1 = present (individual's
    block-level string matches this allele exactly). The value 2 is not
    used for unphased data because the two chromosomes cannot be
    distinguished — an individual homozygous for this allele and one
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
geno_file <- system.file("extdata", "example_genotypes_numeric.csv",
                         package = "LDxBlocks")

res <- run_ldx_pipeline(
  geno_file      = geno_file,
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
#> 1     1  25     rs1001   rs1025     1000  25535   1     24536
#> 2    31  50     rs1031   rs1050    81986 100878   1     18893
#> 3    56  80     rs1056   rs1080   156776 181114   1     24339
#> 4     1  30     rs2001   rs2030     1000  29445   2     28446
#> 5    36  55     rs2036   rs2055    85463 104532   2     19070
#> 6    61  80     rs2061   rs2080   160237 178996   2     18760
head(res$diversity)
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
dim(res$hap_matrix)
#> [1] 120  85
# }
```
