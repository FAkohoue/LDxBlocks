# Within-Block Pairwise SNP Epistasis Scan

Tests pairwise SNP x SNP interactions within haplotype blocks that were
identified as significant by
[`test_block_haplotypes`](https://FAkohoue.github.io/LDxBlocks/reference/test_block_haplotypes.md).
For each significant block, all C(p, 2) SNP pairs are tested for
interaction using the model:

\$\$y = \mu + a_i x_i + a_j x_j + aa\_{ij}(x_i x_j) + \varepsilon\$\$

on GRM-corrected REML residuals (same null model as
[`test_block_haplotypes`](https://FAkohoue.github.io/LDxBlocks/reference/test_block_haplotypes.md)),
ensuring all tests are population-structure-corrected.

Restricting to significant blocks avoids the genome-wide explosion of
pairwise tests: for 15 significant blocks with ~200 SNPs each, the total
number of tests is ~300,000, compared to ~4.4 billion for an
unrestricted genome-wide scan.

## Usage

``` r
scan_block_epistasis(
  assoc,
  geno_matrix,
  snp_info,
  blocks,
  blues,
  haplotypes,
  trait = NULL,
  sig_blocks = NULL,
  min_freq = 0.05,
  max_snps_per_block = 300L,
  sig_threshold = 0.05,
  sig_metric = c("p_simplem_sidak", "p_simplem", "p_bonf", "p_fdr"),
  meff_percent_cut = 0.995,
  id_col = "id",
  blue_col = "blue",
  verbose = TRUE
)
```

## Arguments

- assoc:

  Output of
  [`test_block_haplotypes`](https://FAkohoue.github.io/LDxBlocks/reference/test_block_haplotypes.md).
  Used to identify significant blocks and obtain pre-computed GRM
  residuals.

- geno_matrix:

  Numeric matrix (individuals x SNPs) or `LDxBlocks_backend`. The
  imputed, MAF-filtered genotype matrix from Job 1 (`res$geno_matrix`).

- snp_info:

  Data frame with columns `SNP`, `CHR`, `POS`.

- blocks:

  LD block table from `run_Big_LD_all_chr`.

- blues:

  Pre-adjusted phenotype means (same format as `test_block_haplotypes`).
  Named numeric vector or named list.

- haplotypes:

  Named list from `extract_haplotypes`. Used to re-build the GRM for
  REML residual computation.

- trait:

  Character. Which trait to use for residual computation when `blues` is
  a named list. Default `NULL` uses the first trait.

- sig_blocks:

  Character vector. Block IDs to scan. `NULL` (default) uses all blocks
  with `significant_omnibus = TRUE` in `assoc$block_tests`.

- min_freq:

  Numeric. Minimum MAF for SNPs within the block. Default `0.05`.

- max_snps_per_block:

  Integer. Maximum SNPs per block before switching to random subsampling
  of pairs. Default `300L` (C(300,2)=44,850 pairs per block). Set `NULL`
  to always use all SNPs.

- sig_threshold:

  Numeric. Significance threshold. Default `0.05`.

- sig_metric:

  Character. Which correction drives the primary `significant` flag. One
  of:

  - `"p_simplem_sidak"` (default, recommended) – simpleM Sidak.

  - `"p_simplem"` – simpleM Bonferroni.

  - `"p_bonf"` – plain Bonferroni (p x n_pairs).

  All three p-value columns are always present regardless of this
  choice.

- meff_percent_cut:

  Numeric. Variance cutoff for simpleM Meff. Default `0.995`.

- id_col:

  Character. ID column when blues is a data frame.

- blue_col:

  Character. Phenotype column when blues is a data frame.

- verbose:

  Logical. Default `TRUE`.

## Value

A named list of class `LDxBlocks_epistasis`:

- `results`:

  Data frame. One row per tested SNP pair per block per trait. Always
  present columns: `block_id`, `CHR`, `start_bp`, `end_bp`, `trait`,
  `SNP_i`, `SNP_j`, `POS_i`, `POS_j`, `dist_bp`, `aa_effect`, `SE`,
  `t_stat`, `p_wald`, `Meff` (simpleM effective test count from
  interaction eigenspectrum), `p_bonf` (Bonferroni: p x n_pairs),
  `p_simplem` (simpleM Bonferroni: p x Meff), `p_simplem_sidak` (simpleM
  Sidak: 1-(1-p)^Meff), `significant` (driven by `sig_metric`),
  `significant_bonf`, `significant_simplem`,
  `significant_simplem_sidak`.

- `scan_summary`:

  Data frame. One row per block: number of pairs tested, number
  significant, minimum p-value.

- `n_blocks_scanned`:

  Integer.

- `n_pairs_total`:

  Integer. Total pairwise tests performed.

## References

Cordell HJ (2009). Detecting gene-gene interactions that underlie human
diseases. *Nature Reviews Genetics* **10**:392-404.

## See also

[`test_block_haplotypes`](https://FAkohoue.github.io/LDxBlocks/reference/test_block_haplotypes.md),
[`scan_block_by_block_epistasis`](https://FAkohoue.github.io/LDxBlocks/reference/scan_block_by_block_epistasis.md),
[`fine_map_epistasis_block`](https://FAkohoue.github.io/LDxBlocks/reference/fine_map_epistasis_block.md)
