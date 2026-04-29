# Between-Block Haplotype Allele Epistasis Scan

Tests interactions between haplotype allele dosage columns from
different LD blocks. For each significant allele (from
[`test_block_haplotypes`](https://FAkohoue.github.io/LDxBlocks/reference/test_block_haplotypes.md)),
tests its interaction with every haplotype allele at all other blocks:

\$\$y = \mu + \alpha_i x_i + \alpha_j x_j + \gamma\_{ij}(x_i x_j) +
\varepsilon\$\$

Restricting to significant alleles x all blocks gives O(n_sig x
n_total_alleles) tests – tractable at WGS scale. For your panel: ~25
significant alleles x 17,943 alleles = ~450,000 tests.

This is a **trans-haplotype epistasis scan**: it tests whether the
effect of a resistance haplotype at one block depends on the haplotype
background at another block – a form of genetic background dependence
that single-SNP and even single-block analyses cannot detect.

## Usage

``` r
scan_block_by_block_epistasis(
  assoc,
  haplotypes,
  blues,
  blocks,
  trait = NULL,
  sig_alleles = NULL,
  min_freq = 0.05,
  top_n = NULL,
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

- haplotypes:

  Named list from `extract_haplotypes`.

- blues:

  Pre-adjusted phenotype means.

- blocks:

  LD block table.

- trait:

  Character. Trait to use. Default `NULL` uses first trait.

- sig_alleles:

  Data frame with columns `block_id` and `allele` identifying the query
  alleles. `NULL` (default) uses all `significant = TRUE` rows from
  `assoc$allele_tests`.

- min_freq:

  Numeric. Minimum allele frequency. Default `0.05`.

- top_n:

  Integer or `NULL`. Alleles per block in feature matrix.

- sig_threshold:

  Numeric. Significance threshold. Default `0.05`.

- id_col:

  Character. ID column when blues is a data frame.

- blue_col:

  Character. Phenotype column when blues is a data frame.

- verbose:

  Logical. Default `TRUE`.

## Value

A named list of class `LDxBlocks_block_epistasis`:

- `results`:

  Data frame sorted by p_wald. Columns: `block_i`, `allele_i`,
  `block_j`, `allele_j`, `CHR_i`, `CHR_j`, `same_chr`, `aa_effect`,
  `SE`, `t_stat`, `p_wald`, `p_bonf`, `significant`.

- `n_tests`:

  Integer. Total interaction tests performed.

- `n_sig_alleles`:

  Integer. Number of query alleles.

## References

Mackay TFC (2014). Epistasis and quantitative traits: using model
organisms to study gene-gene interactions. *Nature Reviews Genetics*
**15**:22-33.

## See also

[`scan_block_epistasis`](https://FAkohoue.github.io/LDxBlocks/reference/scan_block_epistasis.md),
[`test_block_haplotypes`](https://FAkohoue.github.io/LDxBlocks/reference/test_block_haplotypes.md),
[`fine_map_epistasis_block`](https://FAkohoue.github.io/LDxBlocks/reference/fine_map_epistasis_block.md)
