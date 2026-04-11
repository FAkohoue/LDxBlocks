# Write Haplotype Diversity Table

Write Haplotype Diversity Table

## Usage

``` r
write_haplotype_diversity(
  diversity,
  out_file,
  append_summary = TRUE,
  verbose = TRUE
)
```

## Arguments

- diversity:

  Data frame from compute_haplotype_diversity().

- out_file:

  Output CSV path.

- append_summary:

  Append genome-wide mean row. Default TRUE.

- verbose:

  Logical. Default TRUE.

## Value

Invisibly returns out_file.
