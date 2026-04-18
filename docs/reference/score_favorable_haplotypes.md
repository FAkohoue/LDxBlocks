# Score Individual Haplotype Portfolios Against Known Allele Effects

Scores each individual's genome-wide haplotype composition against a
table of per-allele effects (from genomic prediction or association
analysis), producing a per-block score and a genome-wide stacking index.
This is the primary tool for translating haplotype genomic prediction
results into actionable breeding rankings.

**Scoring rule:** For each individual \\i\\ and each LD block \\b\\, the
block score is: \$\$s\_{ib} = \sum\_{j} \hat{\alpha}\_j \cdot
d\_{ij}\$\$ where \\\hat{\alpha}\_j\\ is the known effect of allele
\\j\\ and \\d\_{ij}\\ is the dosage (copies) of allele \\j\\ carried by
individual \\i\\ (0/1 for unphased data; 0/1/2 for phased data). The
genome-wide stacking index is the sum of block scores across all scored
blocks, optionally normalised to \[0, 1\].

## Usage

``` r
score_favorable_haplotypes(
  haplotypes,
  allele_effects,
  min_freq = 0.02,
  missing_string = ".",
  normalize = TRUE
)
```

## Arguments

- haplotypes:

  Named list produced by
  [`extract_haplotypes`](https://FAkohoue.github.io/LDxBlocks/reference/extract_haplotypes.md).
  Must carry a `block_info` attribute. Individual IDs are taken from
  `names(haplotypes[[1]])` (names of the first block's haplotype
  vector). All individuals present in `haplotypes` are scored, including
  those not in `allele_effects` (they receive a score of 0 for missing
  blocks).

- allele_effects:

  Data frame specifying known per-allele additive effects. **Required
  columns:**

  - `block_id` (character) — Block identifier matching
    `names(haplotypes)`.

  - `allele` (character) — Haplotype allele string matching the strings
    in `haplotypes[[block_id]]`.

  - `allele_effect` (numeric) — Effect size. Positive = allele increases
    trait value; negative = decreases it. Units are the same as the
    phenotype scale used to estimate effects.

  Accepts the output of
  [`decompose_block_effects`](https://FAkohoue.github.io/LDxBlocks/reference/decompose_block_effects.md)
  directly. When effects from
  [`test_block_haplotypes`](https://FAkohoue.github.io/LDxBlocks/reference/test_block_haplotypes.md)
  are used, filter to one trait first
  (`allele_effects[allele_effects$trait == "YLD", ]`) and rename
  `effect` to `allele_effect`. Blocks or alleles in `haplotypes` with no
  matching entry in `allele_effects` contribute a score of 0 for those
  individuals.

- min_freq:

  Numeric in (0, 1). Minimum allele frequency in the full panel. Alleles
  with population frequency below this threshold are excluded from
  scoring even if they appear in `allele_effects`. Default `0.02`. This
  prevents rare private alleles from dominating scores based on
  unreliable frequency estimates.

- missing_string:

  Character. The string used in haplotype vectors to denote missing
  genotype data (e.g. for individuals with insufficient marker coverage
  at a block). Individuals with this string at a block receive `NA` for
  that block's score (excluded from mean but still contribute to
  `n_blocks_scored` = 0 for that block). Default `"."`.

- normalize:

  Logical. If `TRUE` (default), the genome-wide stacking index is
  linearly scaled to \[0, 1\] across all individuals:
  \$\$\text{index}\_i = (S_i - S\_{\min}) / (S\_{\max} - S\_{\min})\$\$
  where \\S_i\\ is the raw sum of block scores. This makes the index
  interpretable as a percentile rank within the panel and comparable
  across panels with different numbers of scored blocks. Set `FALSE` to
  return raw summed effects (in phenotype units), which is preferable
  when comparing candidate sets across different evaluation sets.

## Value

Data frame with one row per individual, sorted ascending by `rank` (best
candidates first). Contains the following columns:

- `id`:

  Character. Individual identifier, taken from `names(haplotypes[[1]])`.

- `stacking_index`:

  Numeric. Genome-wide sum of per-block haplotype scores. When
  `normalize = TRUE`, scaled to \[0, 1\] where 1.0 = the individual with
  the most favourable genome-wide haplotype combination in the panel and
  0.0 = the least favourable. When `normalize = FALSE`, units are
  phenotype units × allele dosage.

- `n_blocks_scored`:

  Integer. Number of LD blocks for which at least one allele effect was
  available in `allele_effects` and at least one allele passed
  `min_freq`. Maximum possible value equals the number of unique
  `block_id` values in `allele_effects` that overlap with
  `names(haplotypes)`.

- `mean_block_score`:

  Numeric. Raw (unnormalised) mean per-block score:
  `stacking_index_raw / n_blocks_scored`. Useful for comparing
  candidates across panels with different numbers of scored blocks,
  since it is independent of `n_blocks_scored`.

- `rank`:

  Integer. Rank based on `stacking_index`, with 1 indicating the
  individual with the highest (most favourable) index. Ties are broken
  by `"min"` (tied individuals share the lower rank number).

- `score_<block_id>`:

  Numeric. One column per scored LD block, named `score_` followed by
  the block identifier (e.g. `score_block_1_1000_103000`). Contains the
  raw per-block score for each individual (effect × dosage sum across
  alleles). Zero when the individual carries no alleles with known
  effects at that block. Used to identify which genomic regions drive an
  individual's stacking index.

## See also

[`decompose_block_effects`](https://FAkohoue.github.io/LDxBlocks/reference/decompose_block_effects.md),
[`summarize_parent_haplotypes`](https://FAkohoue.github.io/LDxBlocks/reference/summarize_parent_haplotypes.md),
[`run_haplotype_prediction`](https://FAkohoue.github.io/LDxBlocks/reference/run_haplotype_prediction.md)

## Examples

``` r
# \donttest{
data(ldx_geno, ldx_snp_info, ldx_blocks, ldx_blues, package = "LDxBlocks")
haps    <- extract_haplotypes(ldx_geno, ldx_snp_info, ldx_blocks)
hap_mat <- build_haplotype_feature_matrix(haps)
G       <- compute_haplotype_grm(hap_mat)
pred    <- run_haplotype_prediction(ldx_geno, ldx_snp_info, ldx_blocks,
                                     blues    = setNames(ldx_blues$YLD,
                                                         ldx_blues$id),
                                     verbose  = FALSE)
ae <- decompose_block_effects(haps, ldx_snp_info, ldx_blocks,
                               snp_effects = pred$snp_effects)
scores <- score_favorable_haplotypes(haps, allele_effects = ae)
head(scores[order(scores$rank), c("id","stacking_index","rank")])
#>            id stacking_index rank
#> ind083 ind083       1.000000    1
#> ind025 ind025       0.944427    2
#> ind095 ind095       0.892077    3
#> ind108 ind108       0.862556    4
#> ind049 ind049       0.856506    5
#> ind068 ind068       0.847426    6
# }
```
