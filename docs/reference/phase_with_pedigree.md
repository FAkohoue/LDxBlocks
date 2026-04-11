# Pedigree-Based Allele Transmission Phasing

Assigns gametic phase using Mendelian transmission within
parent-offspring trios. Exact when parents are homozygous; ambiguous
positions are NA or randomly resolved when resolve_het=TRUE.

## Usage

``` r
phase_with_pedigree(
  dosage_mat,
  pedigree,
  resolve_het = FALSE,
  seed = 1L,
  verbose = TRUE
)
```

## Arguments

- dosage_mat:

  Numeric matrix (individuals x SNPs, 0/1/2/NA). Rownames must match
  pedigree\$id.

- pedigree:

  Data frame with columns id, sire, dam. Use NA or "0" for unknown
  parents.

- resolve_het:

  Randomly resolve ambiguous het transmissions. Default FALSE.

- seed:

  RNG seed when resolve_het=TRUE. Default 1L.

- verbose:

  Logical. Default TRUE.

## Value

List: hap1 (sire gamete), hap2 (dam gamete), dosage, n_resolved,
n_ambiguous, phased=TRUE.
