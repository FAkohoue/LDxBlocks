# Genome-Wide LD Block Detection by Chromosome

Applies `Big_LD()` chromosome by chromosome, collects results and
returns a single tidy data frame annotated with chromosome and block
length. This is the recommended entry point for genome-wide analyses.

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
  max_bp_distance = 0L,
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

  Forwarded to `Big_LD()`. See that function's documentation for
  details.

- method:

  Character. LD metric: `"r2"` (default, standard squared Pearson
  correlation) or `"rV2"` (kinship-adjusted). See `Big_LD()` for
  details.

- n_threads:

  Integer. Number of OpenMP threads for the C++ LD kernel. Default `1L`.
  Increase for multi-core systems.

- singleton_as_block:

  Logical. If `TRUE`, SNPs that pass MAF filtering but are not assigned
  to any clique are returned as single-SNP blocks (`start == end`,
  `length_bp == 1`). Default `FALSE`. See `Big_LD()` for details.

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

- max_bp_distance:

  Integer. Maximum base-pair distance between a SNP pair for its r\\^2\\
  to be computed. Pairs beyond this distance are set to zero in the
  adjacency matrix (assumed to be in negligible LD). `0L` (default)
  disables this and computes all pairs (original behaviour). Recommended
  value for WGS panels: `500000L` (500 kb). Has no effect when `CLQmode`
  is `"Louvain"` or `"Leiden"` and the window spans less than
  `max_bp_distance`. Requires sorted SNP positions within each
  sub-segment (guaranteed by `run_Big_LD_all_chr`).

## Value

A `data.frame` with columns: `start`, `end`, `start.rsID`, `end.rsID`,
`start.bp`, `end.bp`, `CHR`, `length_bp`. Rows are sorted by `CHR` then
`start.bp`.

## See also

`Big_LD()`,
[`tune_LD_params`](https://FAkohoue.github.io/LDxBlocks/reference/tune_LD_params.md),
[`extract_haplotypes`](https://FAkohoue.github.io/LDxBlocks/reference/extract_haplotypes.md)

## Examples

``` r
# \donttest{
# Use the package example data -- 120 individuals, 230 SNPs, 3 chromosomes,
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
#> 1     1  25     rs1001   rs1025     1000  25027   1     24028
#> 2    31  50     rs1031   rs1050    81064  99022   1     17959
#> 3    56  80     rs1056   rs1080   155368 179371   1     24004
#> 4     1  30     rs2001   rs2030     1000  30023   2     29024
#> 5    36  55     rs2036   rs2055    86236 105290   2     19055
#> 6    61  80     rs2061   rs2080   161515 180473   2     18959
summarise_blocks(blocks)
#>      CHR n_blocks min_bp median_bp  mean_bp max_bp total_bp_covered
#> 1      1        3  17959     24004 21997.00  24028            65991
#> 2      2        3  18959     19055 22346.00  29024            67038
#> 3      3        3  18069     18356 18385.00  18730            55155
#> 4 GENOME        9  17959     18959 20909.33  29024           188184
# }
```
