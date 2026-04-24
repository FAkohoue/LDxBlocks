# Read a Pre-Phased VCF

Reads a phased VCF/VCF.gz file and extracts the two haplotype allele
matrices, allele dosage matrix, SNP metadata, and sample IDs.

This function expects phased genotype calls in the VCF sample columns,
using the pipe separator, for example `0|0`, `0|1`, `1|0`, `1|1`, or
`.|.`. FORMAT subfields are allowed; for example, `0|1:35:99` is parsed
as `0|1`.

The returned `hap1`, `hap2`, and `dosage` matrices are stored as SNPs x
individuals.

## Usage

``` r
read_phased_vcf(vcf_file, min_maf = 0, verbose = TRUE)
```

## Arguments

- vcf_file:

  Character. Path to a phased `.vcf` or `.vcf.gz` file.

- min_maf:

  Numeric. Minimum minor allele frequency used to filter variants after
  dosage calculation. Must be between 0 and 0.5. Default is `0.0`,
  meaning no MAF filtering.

- verbose:

  Logical. If `TRUE`, progress messages are printed. Default is `TRUE`.

## Value

A list with elements:

- `hap1`:

  Integer matrix of first haplotype alleles, SNPs x individuals.

- `hap2`:

  Integer matrix of second haplotype alleles, SNPs x individuals.

- `dosage`:

  Integer matrix of allele dosages, SNPs x individuals.

- `snp_info`:

  Data frame with columns `SNP`, `CHR`, `POS`, `REF`, and `ALT`.

- `sample_ids`:

  Character vector of sample IDs.

- `phased`:

  Logical. Always `TRUE`.

- `phase_method`:

  Character. Always `"vcf_phased"`.
