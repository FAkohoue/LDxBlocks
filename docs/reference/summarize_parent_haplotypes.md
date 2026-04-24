# Summarise Haplotype Allele Inventory Per Candidate Parent

Produces a tidy long-format allele inventory for each candidate parent
individual, reporting which haplotype alleles they carry at each LD
block and in what dosage. This is the primary decision support table for
haplotype stacking crosses: it reveals which parents carry complementary
rare alleles at important blocks, which candidates are genomically
redundant (same alleles at all blocks), and which blocks should be
targeted in crossing schemes.

The output includes one row per (individual x block x allele)
combination, including rows where the individual carries 0 copies of an
allele. This ensures that all candidates can be compared on the same
rows for any given block.

## Usage

``` r
summarize_parent_haplotypes(
  haplotypes,
  candidate_ids = NULL,
  allele_effects = NULL,
  blocks = NULL,
  min_freq = 0.02,
  missing_string = "."
)
```

## Arguments

- haplotypes:

  Named list produced by
  [`extract_haplotypes`](https://FAkohoue.github.io/LDxBlocks/reference/extract_haplotypes.md).
  Must carry a `block_info` attribute. All blocks in the list are
  included in the inventory unless filtered by `candidate_ids` or
  `min_freq`.

- candidate_ids:

  Character vector of individual IDs to include in the inventory. Must
  match names in the haplotype vectors (i.e. `names(haplotypes[[1]])`).
  `NULL` (default) includes all individuals present in `haplotypes`.
  Supply the `id` column from
  [`score_favorable_haplotypes`](https://FAkohoue.github.io/LDxBlocks/reference/score_favorable_haplotypes.md)
  (e.g. top 20 by rank) to focus the inventory on selection candidates.

- allele_effects:

  Optional data frame of per-allele effects to join to the inventory.
  When supplied, the `allele_effect` column is populated by matching on
  `block_id` and `allele`. `NULL` (default) leaves `allele_effect` as
  `NA` throughout. **Required columns when not `NULL`:**

  - `block_id` (character) - matching `names(haplotypes)`.

  - `allele` (character) - matching haplotype allele strings.

  - `allele_effect` (numeric) - effect value per allele.

  Accepts the output of
  [`decompose_block_effects`](https://FAkohoue.github.io/LDxBlocks/reference/decompose_block_effects.md)
  directly. For
  [`test_block_haplotypes`](https://FAkohoue.github.io/LDxBlocks/reference/test_block_haplotypes.md)
  output, filter to one trait and rename `effect` to `allele_effect`
  first.

- blocks:

  LD block table from
  [`run_Big_LD_all_chr`](https://FAkohoue.github.io/LDxBlocks/reference/run_Big_LD_all_chr.md).
  Used only as an alternative source of block coordinate metadata if the
  `block_info` attribute of `haplotypes` is missing or incomplete.
  Default `NULL` (uses `block_info` attribute).

- min_freq:

  Numeric in (0, 1). Minimum allele frequency in the full panel (all
  individuals in `haplotypes`, not just `candidate_ids`). Alleles below
  this threshold are excluded from the inventory entirely. Default
  `0.02`. This prevents the inventory from being dominated by private
  alleles observed in only one or two individuals.

- missing_string:

  Character. String used in haplotype vectors to indicate missing
  genotype data. Individuals with this string at a block are skipped for
  that block (no rows emitted). Default `"."`.

## Value

Data frame in long format. One row per individual x block x allele
combination, including rows where `dosage = 0` (allele absent). Sorted
ascending by `id`, `CHR`, `start_bp`, then descending by `dosage` (so
carried alleles appear before absent ones for each individual-block
combination). Contains 10 columns:

- `id`:

  Character. Individual identifier.

- `block_id`:

  Character. LD block identifier matching `names(haplotypes)`.

- `CHR`:

  Character. Chromosome label.

- `start_bp`:

  Integer. Block start coordinate (base pairs).

- `end_bp`:

  Integer. Block end coordinate (base pairs).

- `allele`:

  Character. Haplotype allele string (dosage-coded concatenation of SNP
  dosages within the block, e.g. `"012102"`).

- `dosage`:

  Integer. Number of copies of this allele carried by the individual: 0
  = absent; 1 = one copy (heterozygous for phased data, or present for
  unphased data); 2 = two copies (homozygous, phased data only). For
  unphased input (no `"|"` in strings), dosage is always 0 or 1 - the
  two chromosomes cannot be distinguished, so an individual matching
  this allele string exactly receives dosage 1 regardless of whether
  they are truly homozygous or heterozygous.

- `allele_freq`:

  Numeric in \[0, 1\]. Population frequency of this allele in the full
  panel (all individuals, not just candidates). Computed as
  `table(valid_haplotypes) / length(valid_haplotypes)` before any
  `candidate_ids` filtering is applied.

- `allele_effect`:

  Numeric or `NA`. Effect value from `allele_effects` if supplied and
  the allele-block combination is matched; `NA` otherwise.

- `is_rare`:

  Logical. `TRUE` when `allele_freq < 0.10`. Convenience flag for
  filtering the inventory to rare alleles that may be private to
  specific parents and therefore valuable for introgression.

## See also

[`score_favorable_haplotypes`](https://FAkohoue.github.io/LDxBlocks/reference/score_favorable_haplotypes.md),
[`decompose_block_effects`](https://FAkohoue.github.io/LDxBlocks/reference/decompose_block_effects.md),
[`compare_haplotype_populations`](https://FAkohoue.github.io/LDxBlocks/reference/compare_haplotype_populations.md)

## Examples

``` r
# \donttest{
data(ldx_geno, ldx_snp_info, ldx_blocks, ldx_blues, package = "LDxBlocks")
haps <- extract_haplotypes(ldx_geno, ldx_snp_info, ldx_blocks)
pred <- run_haplotype_prediction(ldx_geno, ldx_snp_info, ldx_blocks,
                                  blues   = setNames(ldx_blues$YLD,
                                                     ldx_blues$id),
                                  verbose = FALSE)
ae <- decompose_block_effects(haps, ldx_snp_info, ldx_blocks,
                               snp_effects = pred$snp_effects)
# Inventory for top 5 candidates
top5 <- rownames(ldx_geno)[1:5]
inv  <- summarize_parent_haplotypes(haps, candidate_ids = top5,
                                     allele_effects = ae)
# Show blocks where candidates carry different alleles
inv[inv$dosage > 0, c("id","block_id","allele","dosage","allele_effect")]
#>         id              block_id                         allele dosage
#> 10  ind001    block_1_1000_25027      2022222002222002220020220      1
#> 60  ind001   block_1_81064_99022           22000002220202000020      1
#> 102 ind001 block_1_155368_179371      0110110010000000001010100      1
#> 156 ind001    block_2_1000_30023 111122112222221222211122211112      1
#> 259 ind001 block_2_161515_180473           11211121201111111221      1
#> 305 ind001    block_3_1000_19068           02000022220222000202      1
#> 353 ind001   block_3_74532_93854           20220202222222222220      1
#> 399 ind001 block_3_149647_168376           11112111211111112011      1
#> 64  ind002   block_1_81064_99022           11122211120121121200      1
#> 165 ind002    block_2_1000_30023 111011101001101100012000110001      1
#> 215 ind002  block_2_86236_105290           11100211112111111102      1
#> 270 ind002 block_2_161515_180473           22222220202222002222      1
#> 317 ind002    block_3_1000_19068           12111122221222111112      1
#> 357 ind002   block_3_74532_93854           01010020101211111011      1
#> 412 ind002 block_3_149647_168376           22112010110011102121      1
#> 27  ind003    block_1_1000_25027      1011121012111011210020110      1
#> 76  ind003   block_1_81064_99022           12111102120112111110      1
#> 129 ind003 block_1_155368_179371      1221121121212121112101112      1
#> 171 ind003    block_2_1000_30023 000022002222220222200222220002      1
#> 230 ind003  block_2_86236_105290           22200222222222222002      1
#> 271 ind003 block_2_161515_180473           00020020020222200222      1
#> 321 ind003    block_3_1000_19068           01000011210111110101      1
#> 374 ind003   block_3_74532_93854           21211212222112212210      1
#> 422 ind003 block_3_149647_168376           22202000220022202020      1
#> 33  ind004    block_1_1000_25027      0020000000000000002202000      1
#> 85  ind004   block_1_81064_99022           11122221211221121111      1
#> 138 ind004 block_1_155368_179371      1210120110111011112101111      1
#> 184 ind004    block_2_1000_30023 100012001111111111101111120001      1
#> 235 ind004  block_2_86236_105290           11100211112111111102      1
#> 283 ind004 block_2_161515_180473           00110021110111210221      1
#> 328 ind004    block_3_1000_19068           00011100101000210000      1
#> 381 ind004   block_3_74532_93854           11120111212222222121      1
#> 425 ind004 block_3_149647_168376           11022121101100012112      1
#> 50  ind005    block_1_1000_25027      2022222002222002220020220      1
#> 92  ind005   block_1_81064_99022           02222212111122222111      1
#> 142 ind005 block_1_155368_179371      0110110010000000001010100      1
#> 195 ind005    block_2_1000_30023 111011101001101100012000110001      1
#> 250 ind005  block_2_86236_105290           22200222222222222002      1
#> 293 ind005 block_2_161515_180473           00110021110111210221      1
#> 339 ind005    block_3_1000_19068           01000011210111110101      1
#> 394 ind005   block_3_74532_93854           21211212222112212210      1
#> 439 ind005 block_3_149647_168376           22112010110011102121      1
#>     allele_effect
#> 10      -0.020711
#> 60       0.042963
#> 102      0.060325
#> 156      1.265583
#> 259     -0.237346
#> 305     -0.543943
#> 353     -0.195629
#> 399      0.000421
#> 64       0.108826
#> 165      0.424797
#> 215     -0.298645
#> 270     -0.302974
#> 317     -0.068503
#> 357     -0.057814
#> 412      0.021767
#> 27      -0.054256
#> 76       0.125705
#> 129      0.601534
#> 171      1.043313
#> 230     -0.564878
#> 271     -0.512270
#> 321     -0.159478
#> 374     -0.106051
#> 422      0.059172
#> 33      -0.127989
#> 85       0.279913
#> 138      0.386147
#> 184      0.365229
#> 235     -0.298645
#> 283     -0.341994
#> 328      0.333311
#> 381     -0.153709
#> 425     -0.036984
#> 50      -0.020711
#> 92       0.379535
#> 142      0.060325
#> 195      0.424797
#> 250     -0.564878
#> 293     -0.341994
#> 339     -0.159478
#> 394     -0.106051
#> 439      0.021767
# }
```
