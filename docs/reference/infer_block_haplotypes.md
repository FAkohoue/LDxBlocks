# Infer Structured Block-Level Diplotypes Per Individual

Converts the raw haplotype strings produced by
[`extract_haplotypes`](https://FAkohoue.github.io/LDxBlocks/reference/extract_haplotypes.md)
into a structured per-individual, per-block diplotype table with
explicit `hap1` / `hap2` gamete alleles, a diplotype code, a
heterozygosity flag, and a phase-ambiguity flag.

When the input is a **phased list** (from
[`read_phased_vcf`](https://FAkohoue.github.io/LDxBlocks/reference/read_phased_vcf.md)
or
[`phase_with_pedigree`](https://FAkohoue.github.io/LDxBlocks/reference/phase_with_pedigree.md)),
`hap1` and `hap2` are the true gametic strings and `phase_ambiguous` is
always `FALSE`.

When the input is **unphased** (diploid allele strings), heterozygous
individuals have two possible phase assignments. The function flags
these with `phase_ambiguous = TRUE`; `hap1` / `hap2` are set to `NA` for
those individuals unless `resolve_unphased = TRUE`, in which case the
most-frequent gametic split consistent with observed allele frequencies
is imputed (maximum-parsimony heuristic — not statistically rigorous,
use
[`phase_with_beagle`](https://FAkohoue.github.io/LDxBlocks/reference/phase_with_beagle.md)
for rigorous phasing).

## Usage

``` r
infer_block_haplotypes(
  haplotypes,
  resolve_unphased = FALSE,
  missing_string = "."
)
```

## Arguments

- haplotypes:

  Named list from
  [`extract_haplotypes`](https://FAkohoue.github.io/LDxBlocks/reference/extract_haplotypes.md).
  If produced from a phased list input, strings contain `"|"` separators
  (e.g. `"011|100"`). If unphased, strings are plain diploid allele
  codes (e.g. `"012201"`).

- resolve_unphased:

  Logical. For unphased heterozygous genotypes, impute the most
  parsimonious hap1/hap2 split using population allele frequencies
  (parsimony heuristic). Default `FALSE` (leave as `NA`).

- missing_string:

  Character. Missing data placeholder. Default `"."`.

## Value

Data frame with one row per individual × block combination:

- `block_id`:

  Block identifier.

- `CHR`, `start_bp`, `end_bp`:

  Block coordinates.

- `id`:

  Individual ID.

- `hap1`:

  Gamete 1 allele string (`NA` if unresolvable).

- `hap2`:

  Gamete 2 allele string (`NA` if unresolvable).

- `diplotype`:

  Canonical diplotype code: the two gamete alleles sorted alphabetically
  and joined with `"/"`, e.g. `"010/110"`.

- `heterozygous`:

  Logical; `TRUE` when the individual is biologically heterozygous. For
  phased input: `hap1 ≠ hap2`. For unphased input: any dosage position
  equals 1, regardless of whether phase was resolved.

- `phase_ambiguous`:

  Logical; `TRUE` for unphased heterozygous individuals when
  `resolve_unphased = FALSE`. These individuals have
  `heterozygous = TRUE` but `hap1 = hap2 = NA`.

- `missing`:

  Logical; `TRUE` when the allele string contains the missing
  placeholder.

## See also

[`extract_haplotypes`](https://FAkohoue.github.io/LDxBlocks/reference/extract_haplotypes.md),
[`phase_with_beagle`](https://FAkohoue.github.io/LDxBlocks/reference/phase_with_beagle.md),
[`phase_with_pedigree`](https://FAkohoue.github.io/LDxBlocks/reference/phase_with_pedigree.md),
[`collapse_haplotypes`](https://FAkohoue.github.io/LDxBlocks/reference/collapse_haplotypes.md)

## Examples

``` r
# \donttest{
data(ldx_geno, ldx_snp_info, ldx_blocks, package = "LDxBlocks")
haps <- extract_haplotypes(ldx_geno, ldx_snp_info, ldx_blocks)
dip  <- infer_block_haplotypes(haps)
# Proportion of heterozygous diplotypes per block
tapply(dip$heterozygous, dip$block_id, mean, na.rm = TRUE)
#>    block_1_1000_25027 block_1_155368_179371   block_1_81064_99022 
#>             0.7500000             0.7833333             0.6916667 
#>    block_2_1000_30023 block_2_161515_180473  block_2_86236_105290 
#>             0.7750000             0.7916667             0.7166667 
#>    block_3_1000_19068 block_3_149647_168376   block_3_74532_93854 
#>             0.7583333             0.8000000             0.7000000 
# }
```
