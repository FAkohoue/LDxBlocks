# Harmonize Haplotype Allele Labels Across Panels or Analysis Runs

Ensures that haplotype allele labels are biologically comparable across
different datasets, analysis runs, or training/validation splits.
Without harmonization, the allele string `"010110"` in one panel is not
guaranteed to correspond to the same biological haplotype in another
panel if block boundaries, SNP ordering, or allele encoding differ
between runs.

The function anchors allele identity to a **reference dictionary** built
from a training/reference panel. New (target) haplotypes are then
matched against this dictionary:

1.  **Exact match**: the allele string exists verbatim in the reference
    dictionary -\> labelled with the reference allele label.

2.  **Nearest-Hamming match**: no exact match -\> labelled with the most
    similar reference allele (minimum Hamming distance). If the minimum
    Hamming distance exceeds `max_hamming`, the allele is labelled
    `"<novel>"`.

3.  **Novel**: distance \> `max_hamming` -\> `"<novel>"`.

## Usage

``` r
harmonize_haplotypes(
  haplotypes_target,
  haplotypes_ref,
  min_freq_ref = 0.02,
  max_hamming = NULL,
  missing_string = "."
)
```

## Arguments

- haplotypes_target:

  Named list from
  [`extract_haplotypes`](https://FAkohoue.github.io/LDxBlocks/reference/extract_haplotypes.md)
  (the panel to harmonize - validation set, new environment, etc.).

- haplotypes_ref:

  Named list from
  [`extract_haplotypes`](https://FAkohoue.github.io/LDxBlocks/reference/extract_haplotypes.md)
  (the reference panel - training set, base population, etc.). Must
  cover the same blocks as `haplotypes_target` (extra blocks in either
  panel are silently skipped).

- min_freq_ref:

  Numeric. Only alleles above this frequency in the reference panel form
  the dictionary. Default `0.02`.

- max_hamming:

  Integer. Maximum Hamming distance for a nearest-neighbour match;
  alleles beyond this distance are labelled `"<novel>"`. Default `NULL`
  (no limit - always assigns to nearest reference allele).

- missing_string:

  Character. Missing data placeholder. Default `"."`.

## Value

Named list of the same structure as `haplotypes_target`, with allele
strings replaced by their reference-anchored equivalents. The
`block_info` attribute from `haplotypes_target` is preserved. A
`harmonization_report` attribute is attached: a data frame with one row
per block reporting `n_exact`, `n_nearest`, `n_novel`, and
`mean_hamming_dist` for matched alleles.

## See also

[`extract_haplotypes`](https://FAkohoue.github.io/LDxBlocks/reference/extract_haplotypes.md),
[`collapse_haplotypes`](https://FAkohoue.github.io/LDxBlocks/reference/collapse_haplotypes.md),
[`build_haplotype_feature_matrix`](https://FAkohoue.github.io/LDxBlocks/reference/build_haplotype_feature_matrix.md)

## Examples

``` r
# \donttest{
data(ldx_geno, ldx_snp_info, ldx_blocks, package = "LDxBlocks")
# Split into training (70 pct) and validation (30 pct)
n    <- nrow(ldx_geno)
idx  <- sample(n)
ref_geno  <- ldx_geno[idx[1:round(n*0.7)], ]
tgt_geno  <- ldx_geno[idx[(round(n*0.7)+1):n], ]
haps_ref  <- extract_haplotypes(ref_geno, ldx_snp_info, ldx_blocks)
haps_tgt  <- extract_haplotypes(tgt_geno, ldx_snp_info, ldx_blocks)
haps_harm <- harmonize_haplotypes(haps_tgt, haps_ref)
attr(haps_harm, "harmonization_report")
#>                                    block_id n_exact n_nearest n_novel
#> block_1_1000_25027       block_1_1000_25027      24        12       0
#> block_1_81064_99022     block_1_81064_99022      32         4       0
#> block_1_155368_179371 block_1_155368_179371      31         5       0
#> block_2_1000_30023       block_2_1000_30023      28         8       0
#> block_2_86236_105290   block_2_86236_105290      31         5       0
#> block_2_161515_180473 block_2_161515_180473      33         3       0
#> block_3_1000_19068       block_3_1000_19068      30         6       0
#> block_3_74532_93854     block_3_74532_93854      33         3       0
#> block_3_149647_168376 block_3_149647_168376      31         5       0
#>                       mean_hamming_dist
#> block_1_1000_25027             4.000000
#> block_1_81064_99022            1.000000
#> block_1_155368_179371          1.400000
#> block_2_1000_30023             4.250000
#> block_2_86236_105290           1.000000
#> block_2_161515_180473          1.000000
#> block_3_1000_19068             1.166667
#> block_3_74532_93854            1.000000
#> block_3_149647_168376          1.000000
# }
```
