# LD Block Segmentation (r² or rV², C++ accelerated)

Core per-chromosome LD block detection. Two LD metrics are supported:

- `method = "r2"` (default):

  Standard squared Pearson correlation, computed by the C++ Armadillo
  kernel. No kinship matrix is required. Suitable for large datasets
  (hundreds of thousands to millions of markers) and unstructured or
  mildly structured populations.

- `method = "rV2"`:

  Kinship-adjusted squared correlation (Kim et al. 2018). Requires
  computing and inverting a GRM via AGHmatrix + ASRgenomics. Recommended
  for highly related populations (livestock, inbred lines) with moderate
  marker counts (\< 200 k per chr).

For genome-wide analyses use
[`run_Big_LD_all_chr`](https://FAkohoue.github.io/LDxBlocks/reference/run_Big_LD_all_chr.md).
For automatic parameter tuning use
[`tune_LD_params`](https://FAkohoue.github.io/LDxBlocks/reference/tune_LD_params.md).

## Usage

``` r
Big_LD(
  geno,
  SNPinfo,
  method = c("r2", "rV2"),
  CLQcut = 0.5,
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
  n_threads = 1L,
  seed = NULL,
  verbose = FALSE
)
```

## Arguments

- geno:

  Numeric matrix (individuals x SNPs), 0/1/2. Row names used as
  individual IDs; auto-generated if absent.

- SNPinfo:

  Data frame: col 1 = rsID, col 2 = bp position.

- method:

  Character. `"r2"` (default) or `"rV2"`.

- CLQcut:

  Numeric in (0,1\]. LD threshold for clique edges. Default 0.5.

- clstgap:

  Integer. Max bp gap within clique (split=TRUE). Default 40000.

- leng:

  Integer. Boundary-scan half-window (SNPs). Default 200.

- subSegmSize:

  Integer. Max SNPs per CLQD call. Default 1500.

- MAFcut:

  Numeric. Minor allele frequency minimum. Default 0.05.

- appendrare:

  Logical. Append rare SNPs after block detection. Default FALSE.

- singleton_as_block:

  Logical. If `TRUE`, every SNP that passes MAF filtering but has
  pairwise r² below `CLQcut` with all neighbours (i.e. is not assigned
  to any clique) is returned as a single-SNP block with `start == end`
  and `length_bp == 1`. Default `FALSE`. These blocks are excluded from
  haplotype analysis by the default `min_snps = 3` threshold in
  [`extract_haplotypes()`](https://FAkohoue.github.io/LDxBlocks/reference/extract_haplotypes.md),
  but are useful for auditing coverage and for single-SNP feature
  engineering.

- checkLargest:

  Logical. Dense-core pre-pass for large windows. Default FALSE.

- CLQmode:

  `"Density"` (default) or `"Maximal"`.

- kin_method:

  Character. GRM whitening: `"chol"` (default) or `"eigen"`. Ignored
  when `method = "r2"`.

- split:

  Logical. Split cliques at gaps \> clstgap. Default FALSE.

- digits:

  Integer. LD rounding: `-1` (default, no rounding) or a positive
  integer.

- n_threads:

  Integer. OpenMP threads for C++ LD kernel. Default 1.

- seed:

  Integer or NULL. Set for reproducibility.

- verbose:

  Logical. Progress messages.

## Value

data.frame with columns: start, end, start.rsID, end.rsID, start.bp,
end.bp.

## References

Kim S-A et al. (2018) GENETICS 209(3):855-868.  
VanRaden PM (2008) J. Dairy Sci. 91(11):4414-4423.

## See also

[`run_Big_LD_all_chr`](https://FAkohoue.github.io/LDxBlocks/reference/run_Big_LD_all_chr.md),
[`CLQD`](https://FAkohoue.github.io/LDxBlocks/reference/CLQD.md),
[`compute_r2`](https://FAkohoue.github.io/LDxBlocks/reference/compute_r2.md),
[`compute_rV2`](https://FAkohoue.github.io/LDxBlocks/reference/compute_rV2.md),
[`tune_LD_params`](https://FAkohoue.github.io/LDxBlocks/reference/tune_LD_params.md)

## Examples

``` r
# \donttest{
set.seed(1)
m <- 80; b1 <- 60; b2 <- 60
make_block <- function(m, size, p, flip = 0.03) {
  seed_snp <- rbinom(m, 2, p)
  M <- matrix(seed_snp, nrow = m, ncol = size)
  for (j in seq_len(size)) {
    idx <- sample.int(m, max(1L, floor(flip * m)))
    M[idx, j] <- pmin(2L, pmax(0L, M[idx, j] + sample(c(-1L, 1L), length(idx), TRUE)))
  }
  M
}
G <- cbind(make_block(m, b1, 0.30), make_block(m, b2, 0.60))
colnames(G) <- paste0("rs", seq_len(ncol(G)))
rownames(G) <- paste0("ind", seq_len(nrow(G)))
pos <- c(seq(1, by = 1000, length.out = b1), seq(5e6, by = 1000, length.out = b2))
SNPinfo <- data.frame(SNP = colnames(G), POS = pos)
blocks <- Big_LD(G, SNPinfo, method = "r2", CLQcut = 0.6,
                 leng = 30, subSegmSize = 120, verbose = FALSE)
head(blocks)
#>   start end start.rsID end.rsID start.bp  end.bp
#> 1     1  60        rs1     rs60    1e+00   59001
#> 2    61 120       rs61    rs120    5e+06 5059000
# }
```
