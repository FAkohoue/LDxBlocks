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
set.seed(42)
geno <- matrix(sample(0:2, 80 * 300, replace = TRUE), 80, 300)
rownames(geno) <- paste0("ind", 1:80)
colnames(geno) <- paste0("rs", 1:300)
snp_info <- data.frame(
  CHR = rep(paste0("chr", 1:3), each = 100),
  SNP = colnames(geno),
  POS = unlist(lapply(1:3, function(i) sort(sample(1e6:200e6, 100))))
)
blocks <- run_Big_LD_all_chr(geno, snp_info, CLQcut = 0.6, verbose = FALSE)
head(blocks)
#> data frame with 0 columns and 0 rows
# }
```
