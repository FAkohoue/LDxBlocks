# Summarise LD Block Characteristics

Computes per-chromosome and genome-wide summary statistics for a block
table returned by
[`run_Big_LD_all_chr`](https://FAkohoue.github.io/LDxBlocks/reference/run_Big_LD_all_chr.md).

## Usage

``` r
summarise_blocks(blocks)
```

## Arguments

- blocks:

  Data frame of LD blocks. Must contain at least `start.bp`, `end.bp`.
  If a `CHR` column is present, per-chromosome summaries are produced.

## Value

A `data.frame` with one row per chromosome (plus one row for the
genome-wide summary) and columns:

- `CHR`:

  Chromosome (or `"GENOME"`).

- `n_blocks`:

  Number of blocks.

- `min_bp`, `median_bp`, `mean_bp`, `max_bp`:

  Block size distribution in base pairs.

- `total_bp_covered`:

  Sum of block sizes (bp).

## Examples

``` r
# \donttest{
blocks <- data.frame(
  CHR      = c("chr1","chr1","chr2"),
  start.bp = c(1000, 500000, 2000),
  end.bp   = c(80000, 650000, 200000)
)
summarise_blocks(blocks)
#>      CHR n_blocks min_bp median_bp  mean_bp max_bp total_bp_covered
#> 1   chr1        2  79001    114501 114501.0 150001           229002
#> 2   chr2        1 198001    198001 198001.0 198001           198001
#> 3 GENOME        3  79001    150001 142334.3 198001           427003
# }
```
