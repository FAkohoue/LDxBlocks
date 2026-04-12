# Write Haplotype Character (Nucleotide) Matrix

Writes a matrix where each cell contains the nucleotide sequence of the
haplotype allele carried by each individual. Rows are haplotype alleles,
columns are individuals. This is the most interpretable format: you can
read directly which nucleotides define each haplotype allele and which
individuals carry it.

## Usage

``` r
write_haplotype_character(
  haplotypes,
  snp_info,
  out_file,
  min_freq = 0.01,
  top_n = NULL,
  missing_string = ".",
  verbose = TRUE
)
```

## Arguments

- haplotypes:

  List from
  [`extract_haplotypes`](https://FAkohoue.github.io/LDxBlocks/reference/extract_haplotypes.md).

- snp_info:

  Data frame with `SNP`, `CHR`, `POS`, `REF`, `ALT`.

- out_file:

  Output file path (tab-delimited).

- min_freq:

  Minimum haplotype frequency. Default `0.01`.

- top_n:

  Integer or `NULL`. Cap alleles per block. `NULL` (default) keeps all
  above `min_freq`.

- missing_string:

  Missing genotype marker. Default `"."`.

- verbose:

  Logical. Default `TRUE`.

## Value

Invisibly returns `out_file`.

## Details

The cell value for individual i at haplotype allele h is:

- The nucleotide sequence (e.g. `"AGTTA"`) if the individual carries
  that allele (dosage = 2 for unphased, or present in either gamete for
  phased).

- `"-"` if the individual does not carry that allele.

- `"."` if the individual has missing data in that block.
