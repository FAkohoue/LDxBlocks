# Auto-Tune LD Block Detection Parameters

Performs a grid search over
[`Big_LD`](https://FAkohoue.github.io/LDxBlocks/reference/Big_LD.md)
parameters and selects the combination that minimises, in order of
priority:

1.  Unassigned GWAS markers (markers not falling in any block).

2.  Forced assignments (nearest-block fall-back).

3.  Number of blocks (parsimony).

4.  Deviation from `target_bp_band` (biological plausibility).

If `prefer_perfect = TRUE` (default), combinations achieving zero
unassigned and zero forced assignments are prioritised among the above.

After selecting the best parameter set, `tune_LD_params` runs
[`run_Big_LD_all_chr`](https://FAkohoue.github.io/LDxBlocks/reference/run_Big_LD_all_chr.md)
on all chromosomes and assigns every GWAS marker to a block, returning
the final blocks and assignments.

## Usage

``` r
tune_LD_params(
  geno_matrix,
  snp_info,
  gwas_df,
  grid = NULL,
  chromosomes = NULL,
  target_bp_band = c(50000, 5e+05),
  parallel = FALSE,
  seed = NULL,
  prefer_perfect = TRUE,
  return_all_perfect = TRUE
)
```

## Arguments

- geno_matrix:

  Numeric matrix (individuals x SNPs; 0/1/2), genome-wide.

- snp_info:

  Data frame with columns `SNP`, `CHR`, `POS`.

- gwas_df:

  Data frame with columns `Marker`, `CHR`, `POS`.

- grid:

  Optional data frame of parameter combinations. Each row is one
  combination; columns must match parameter names of `Big_LD`. If `NULL`
  (default), a sensible grid over `CLQcut` (4 values) and `min_freq` (2
  values) is used, giving 8 combinations. Both are treated as
  hyperparameters following Weber SE et al. (2023, *Front. Plant Sci.*
  **14**:1217589,
  [doi:10.3389/fpls.2023.1217589](https://doi.org/10.3389/fpls.2023.1217589)
  ), who show that no single threshold is universally optimal across
  datasets and traits.

- chromosomes:

  Optional character vector of chromosome names to include in tuning.
  `NULL` uses all chromosomes in `snp_info`.

- target_bp_band:

  Length-2 numeric vector: preferred median block size range in base
  pairs. Default `c(5e4, 5e5)` (50 kb - 500 kb).

- parallel:

  Logical. If `TRUE`, uses
  [`future.apply::future_lapply`](https://future.apply.futureverse.org/reference/future_lapply.html)
  for parallelism (user must set a `future` plan before calling).
  Default `FALSE`.

- seed:

  Integer seed for reproducibility. Default `NULL`.

- prefer_perfect:

  Logical. Give priority to parameter sets with zero unassigned / zero
  forced GWAS markers. Default `TRUE`.

- return_all_perfect:

  Logical. Include a table of all zero-zero combinations in the returned
  list. Default `TRUE`.

## Value

A named list:

- `best_params`:

  Named list of the selected parameters.

- `score_table`:

  Data frame of all grid combinations and their scores.

- `perfect_table`:

  Data frame of all zero-zero combinations (or `NULL` if none found /
  `return_all_perfect = FALSE`).

- `final_blocks`:

  Block table produced with `best_params` on all chromosomes (output of
  [`run_Big_LD_all_chr`](https://FAkohoue.github.io/LDxBlocks/reference/run_Big_LD_all_chr.md)).

- `gwas_assigned`:

  Input `gwas_df` with an added column `LD_block`. Entries ending in `*`
  denote forced assignments.

## References

Weber SE, Frisch M, Snowdon RJ, Voss-Fels KP (2023). Haplotype blocks
for genomic prediction: a comparative evaluation in multiple crop
datasets. *Frontiers in Plant Science* **14**:1217589.
[doi:10.3389/fpls.2023.1217589](https://doi.org/10.3389/fpls.2023.1217589)

Difabachew YF et al. (2023). Genomic prediction with haplotype blocks in
wheat. *Frontiers in Plant Science* **14**:1168547.
[doi:10.3389/fpls.2023.1168547](https://doi.org/10.3389/fpls.2023.1168547)

## See also

[`run_Big_LD_all_chr`](https://FAkohoue.github.io/LDxBlocks/reference/run_Big_LD_all_chr.md),
[`Big_LD`](https://FAkohoue.github.io/LDxBlocks/reference/Big_LD.md)
