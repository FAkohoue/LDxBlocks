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
  hap_format = c("numeric", "hapmap"),
  maf_cut = 0.05,
  CLQcut = 0.5,
  method = c("r2", "rV2"),
  leng = 200L,
  subSegmSize = 1500L,
  n_threads = 1L,
  min_snps_chr = 10L,
  min_snps_block = 3L,
  top_n = 5L,
  scale_hap_matrix = FALSE,
  chr = NULL,
  verbose = TRUE
)
```

## Arguments

- geno_file:

  Path to genotype file. Supported formats: numeric dosage CSV
  (\`.csv\`), HapMap (\`.hmp.txt\`), VCF (\`.vcf\`, \`.vcf.gz\`),
  SeqArray GDS (\`.gds\`), PLINK BED (\`.bed\`). Format is detected
  automatically from the file extension.

- out_blocks:

  Path for the LD block table CSV. Columns: \`CHR\`, \`start\`, \`end\`,
  \`start.rsID\`, \`end.rsID\`, \`start.bp\`, \`end.bp\`, \`length_bp\`,
  \`length_snps\`, \`block_name\`.

- out_diversity:

  Path for the haplotype diversity table CSV. Columns: \`block_id\`,
  \`n_ind\`, \`n_haplotypes\`, \`He\`, \`Shannon\`, \`freq_dominant\`.

- out_hap_matrix:

  Path for the haplotype genotype matrix. Format is controlled by
  \`hap_format\`. See section \*\*Haplotype genotype matrix\*\*.

- hap_format:

  Output format for the haplotype matrix:

  - \`"numeric"\` (default) — CSV with rows = individuals, columns =
    haplotype alleles coded 0/2/NA.

  - \`"hapmap"\` — HapMap format with rows = haplotype alleles, columns
    = individuals, nucleotide encoding.

- maf_cut:

  Minimum minor allele frequency. SNPs below this threshold are removed
  before block detection. Default \`0.05\`.

- CLQcut:

  r² threshold for clique edges in CLQD. Higher values produce tighter,
  smaller blocks. Default \`0.5\`.

- method:

  LD metric: \`"r2"\` (default) or \`"rV2"\` (requires kinship matrix;
  see \`Big_LD()\`).

- leng:

  Boundary scan half-window in SNPs. Default \`200L\`. Reduce to 50–100
  for very dense WGS panels.

- subSegmSize:

  Maximum SNPs per CLQD sub-segment. Controls peak RAM: \`subSegmSize ×
  n_individuals × 8\` bytes. Default \`1500L\`.

- n_threads:

  OpenMP threads for the C++ LD kernel. Default \`1L\`.

- min_snps_chr:

  Skip chromosomes with fewer post-filter SNPs than this. Default
  \`10L\`. Increase to skip unplaced scaffolds.

- min_snps_block:

  Minimum SNPs per haplotype block. Blocks smaller than this are
  excluded from haplotype analysis. Default \`3L\`.

- top_n:

  Number of top haplotype alleles to retain per block in the output
  matrix. Default \`5L\`.

- scale_hap_matrix:

  Logical. If \`TRUE\`, scale the haplotype matrix columns to zero mean
  and unit variance before writing. Useful for GBLUP-style models.
  Default \`FALSE\`.

- chr:

  Character vector of chromosome names to process. \`NULL\` (default)
  processes all chromosomes.

- verbose:

  Logical. Print timestamped progress. Default \`TRUE\`.

## Value

A named list (invisibly) with elements:

- \`blocks\`:

  \`data.frame\` of LD blocks.

- \`diversity\`:

  \`data.frame\` of per-block haplotype diversity metrics.

- \`hap_matrix\`:

  Numeric matrix of haplotype dosages (individuals × haplotype alleles).
  \`NULL\` if written to file only.

- \`snp_info_filtered\`:

  \`data.frame\` of SNP metadata after MAF filter.

- \`n_blocks\`:

  Integer. Total blocks detected.

- \`n_hap_columns\`:

  Integer. Total haplotype allele columns.

## Scale behaviour

Files are processed via the \`LDxBlocks_backend\` streaming interface.
For VCF, HapMap, and numeric CSV files the backend reads one chromosome
window at a time — peak RAM equals one \`subSegmSize\`-SNP window
regardless of total marker count. For very large files (\> 2 M SNPs) the
GDS backend via \`SeqArray\` is used automatically if the \`SeqArray\`
package is installed.

## Haplotype genotype matrix

The haplotype matrix has one row per individual and one column per
haplotype allele (top-\`top_n\` haplotypes per block). Each cell
contains the dosage of that haplotype allele encoded as:

- \`0\` — individual does not carry this haplotype

- \`2\` — individual carries this haplotype (homozygous)

- \`NA\` — missing data in block

This encoding is directly compatible with genomic prediction software
(ASReml-R, rrBLUP, BGLR, GBLUP) without further transformation.

## See also

\[run_Big_LD_all_chr()\], \[extract_haplotypes()\],
\[compute_haplotype_diversity()\], \[build_haplotype_feature_matrix()\]

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
  top_n          = 3L,
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
dim(res$hap_matrix)
#> [1] 120  27
# }
```
