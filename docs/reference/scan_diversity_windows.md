# Sliding-Window Genome-Wide Diversity Scan

Computes haplotype diversity metrics (He, Shannon entropy,
n_eff_alleles, dominant frequency) in sliding windows across the genome,
independently of LD block boundaries. Useful for identifying diversity
valleys (bottlenecks, selective sweeps) and comparing wild/elite panels
without needing pre-defined blocks.

## Usage

``` r
scan_diversity_windows(
  geno_matrix,
  snp_info,
  window_bp = 1000000L,
  step_bp = 500000L,
  min_snps_win = 5L,
  missing_val = NA
)
```

## Arguments

- geno_matrix:

  Numeric matrix (individuals x SNPs), 0/1/2/NA.

- snp_info:

  Data frame with `SNP`, `CHR`, `POS`.

- window_bp:

  Integer. Window size in base pairs. Default `1e6L` (1 Mb).

- step_bp:

  Integer. Step size in base pairs. Default `5e5L` (500 kb, i.e. 50%
  overlap).

- min_snps_win:

  Integer. Minimum SNPs in a window to compute diversity (windows with
  fewer are skipped). Default `5L`.

- missing_val:

  Numeric. Value representing missing data in `geno_matrix`. Default
  `NA`.

## Value

Data frame with one row per sliding window, sorted by `CHR` then
`win_start`. Columns:

- `CHR`:

  Chromosome label.

- `win_start`, `win_end`:

  Window boundaries (bp).

- `win_mid`:

  Window midpoint (bp).

- `n_snps`:

  Number of SNPs in the window.

- `n_ind`:

  Number of individuals with non-missing data.

- `n_haplotypes`:

  Number of distinct haplotype strings.

- `He`:

  Nei (1973) expected heterozygosity, sample-size corrected.

- `Shannon`:

  Shannon entropy of haplotype frequencies.

- `n_eff_alleles`:

  Effective number of alleles (1/sum(p_i^2)).

- `freq_dominant`:

  Frequency of the most common haplotype.

- `sweep_flag`:

  Logical; TRUE when freq_dominant \>= 0.90.

## See also

[`compute_haplotype_diversity`](https://FAkohoue.github.io/LDxBlocks/reference/compute_haplotype_diversity.md),
[`compare_haplotype_populations`](https://FAkohoue.github.io/LDxBlocks/reference/compare_haplotype_populations.md)

## Examples

``` r
# \donttest{
data(ldx_geno, ldx_snp_info, package = "LDxBlocks")
scan <- scan_diversity_windows(
  geno_matrix  = ldx_geno,
  snp_info     = ldx_snp_info,
  window_bp    = 50000L,
  step_bp      = 25000L,
  min_snps_win = 3L
)
# Plot He across chromosome 1
chr1 <- scan[scan$CHR == "1", ]
plot(chr1$win_mid / 1e3, chr1$He, type = "l",
     xlab = "Position (kb)", ylab = "He",
     main = "Haplotype diversity scan - chr 1")

# }
```
