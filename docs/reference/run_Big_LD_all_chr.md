# Genome-Wide LD Block Detection by Chromosome

Applies
[`Big_LD`](https://FAkohoue.github.io/LDxBlocks/reference/Big_LD.md)
chromosome by chromosome, collects results and returns a single tidy
data frame annotated with chromosome and block length. This is the
recommended entry point for genome-wide analyses.

## Usage

``` r
run_Big_LD_all_chr(
  geno_matrix,
  snp_info = NULL,
  CLQcut = 0.5,
  method = c("r2", "rV2"),
  n_threads = 1L,
  clstgap = 40000,
  leng = 200,
  subSegmSize = 1500,
  MAFcut = 0.05,
  appendrare = FALSE,
  singleton_as_block = FALSE,
  checkLargest = FALSE,
  CLQmode = "Density",
  kin_method = "chol",
  split = FALSE,
  digits = -1L,
  seed = NULL,
  min_snps_chr = 10L,
  chr = NULL,
  clean_malformed = FALSE,
  verbose = FALSE
)
```

## Arguments

- geno_matrix:

  Numeric matrix (individuals x SNPs; values 0/1/2), spanning all
  chromosomes.

- snp_info:

  Data frame with columns `CHR`, `SNP`, `POS`. Column order is flexible;
  columns are matched by name.

- CLQcut, clstgap, leng, subSegmSize, MAFcut, appendrare, checkLargest,
  CLQmode, kin_method, split, digits, seed, verbose:

  Forwarded to
  [`Big_LD`](https://FAkohoue.github.io/LDxBlocks/reference/Big_LD.md).
  See that function's documentation for details.

- method:

  Character. LD metric: `"r2"` (default, standard squared Pearson
  correlation) or `"rV2"` (kinship-adjusted). See
  [`Big_LD`](https://FAkohoue.github.io/LDxBlocks/reference/Big_LD.md)
  for details.

- n_threads:

  Integer. Number of OpenMP threads for the C++ LD kernel. Default `1L`.
  Increase for multi-core systems.

- singleton_as_block:

  Logical. If `TRUE`, SNPs that pass MAF filtering but are not assigned
  to any clique are returned as single-SNP blocks (`start == end`,
  `length_bp == 1`). Default `FALSE`. See
  [`Big_LD`](https://FAkohoue.github.io/LDxBlocks/reference/Big_LD.md)
  for details.

- min_snps_chr:

  Integer. Chromosomes with fewer SNPs than this after MAF filtering are
  skipped. Default `10L`. Increase to skip small scaffolds in
  whole-genome datasets.

- chr:

  Character vector or `NULL`. If supplied, only the named chromosomes
  are processed. Labels must match the values in `snp_info$CHR` after
  normalisation (no `chr` prefix). E.g. `chr = c("1","3","5")` or
  `chr = "X"`. Default `NULL` processes all chromosomes.

- clean_malformed:

  Logical. If `TRUE`, stream-clean the input file before reading by
  removing lines whose column count does not match the header. Only
  relevant when `geno_matrix` is a file path wrapped into a backend.
  Default `FALSE`.

## Value

A `data.frame` with columns: `start`, `end`, `start.rsID`, `end.rsID`,
`start.bp`, `end.bp`, `CHR`, `length_bp`. Rows are sorted by `CHR` then
`start.bp`.

## See also

[`Big_LD`](https://FAkohoue.github.io/LDxBlocks/reference/Big_LD.md),
[`tune_LD_params`](https://FAkohoue.github.io/LDxBlocks/reference/tune_LD_params.md),
[`extract_haplotypes`](https://FAkohoue.github.io/LDxBlocks/reference/extract_haplotypes.md)

## Examples

``` r
# \donttest{
# Use the package example data — 120 individuals, 230 SNPs, 3 chromosomes,
# 9 simulated LD blocks (3 per chromosome).
data(ldx_geno,     package = "LDxBlocks")
data(ldx_snp_info, package = "LDxBlocks")
blocks <- run_Big_LD_all_chr(
  ldx_geno, ldx_snp_info,
  method = "r2", CLQcut = 0.55, leng = 15L,
  subSegmSize = 100L, verbose = FALSE
)
head(blocks)
#>   start end start.rsID end.rsID start.bp end.bp CHR length_bp
#> 1     1  25     rs1001   rs1025     1000  25535   1     24536
#> 2    31  50     rs1031   rs1050    81986 100878   1     18893
#> 3    56  80     rs1056   rs1080   156776 181114   1     24339
#> 4     1  30     rs2001   rs2030     1000  29445   2     28446
#> 5    36  55     rs2036   rs2055    85463 104532   2     19070
#> 6    61  80     rs2061   rs2080   160237 178996   2     18760
summarise_blocks(blocks)
#>      CHR n_blocks min_bp median_bp  mean_bp max_bp total_bp_covered
#> 1      1        3  18893     24339 22589.33  24536            67768
#> 2      2        3  18760     19070 22092.00  28446            66276
#> 3      3        3  18995     19653 19481.33  19796            58444
#> 4 GENOME        9  18760     19653 21387.56  28446           192488
# }
```
