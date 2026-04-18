# Decode Haplotype Strings to Nucleotide Sequences

Converts the dosage-encoded haplotype strings produced by
[`extract_haplotypes`](https://FAkohoue.github.io/LDxBlocks/reference/extract_haplotypes.md)
(e.g. `"02110"`) into nucleotide sequences (e.g. `"AGTT?"`) using the
REF and ALT alleles of each SNP in the block.

## Usage

``` r
decode_haplotype_strings(
  haplotypes,
  snp_info,
  min_freq = 0.01,
  top_n = NULL,
  missing_string = "."
)
```

## Arguments

- haplotypes:

  List from
  [`extract_haplotypes`](https://FAkohoue.github.io/LDxBlocks/reference/extract_haplotypes.md).

- snp_info:

  Data frame with columns `SNP`, `CHR`, `POS`, `REF`, `ALT`. Must
  contain all SNPs in the blocks.

- min_freq:

  Minimum haplotype frequency to include. Default `0.01`.

- top_n:

  Integer or `NULL`. Maximum alleles per block. `NULL` (default) retains
  all above `min_freq`.

- missing_string:

  Missing genotype marker. Default `"."`.

## Value

A data frame with columns:

- block_id:

  Block identifier.

- CHR:

  Chromosome.

- start_bp, end_bp:

  Block boundaries.

- hap_rank:

  Rank by frequency (1 = most common).

- hap_id:

  Column name as it appears in the feature matrix.

- dosage_string:

  Raw dosage string e.g. `"02110"`.

- nucleotide_sequence:

  Decoded nucleotide string e.g. `"AGTT?"`.

- frequency:

  Observed frequency across non-missing individuals.

- n_carriers:

  Number of individuals carrying this haplotype (dosage \> 0).

- snp_positions:

  Semicolon-separated CHR:POS of each SNP in the block.

- snp_alleles:

  Semicolon-separated REF/ALT for each SNP.

## Details

Each character in a haplotype string is the dosage at one SNP in the
block:

- `"0"` = homozygous REF -\> REF nucleotide (e.g. `A`)

- `"1"` = heterozygous -\> IUPAC ambiguity code (e.g. `R` for A/G)

- `"2"` = homozygous ALT -\> ALT nucleotide (e.g. `G`)

- `"."` = missing -\> `N`

The result is a data frame with one row per unique haplotype allele per
block, showing its nucleotide sequence, frequency, and the REF/ALT at
each SNP position. This is the most interpretable representation of what
each haplotype allele actually encodes biologically.

## Examples

``` r
data(ldx_geno, ldx_snp_info, ldx_blocks, package = "LDxBlocks")
haps <- extract_haplotypes(ldx_geno, ldx_snp_info, ldx_blocks, min_snps = 3)
decoded <- decode_haplotype_strings(haps, ldx_snp_info)
head(decoded[, c("block_id","hap_rank","dosage_string",
                  "nucleotide_sequence","frequency")])
#>             block_id hap_rank             dosage_string
#> 1 block_1_1000_25027        1 1011121012111011210020110
#> 2 block_1_1000_25027        2 0111021121110121201121101
#> 3 block_1_1000_25027        3 0010010011000010101111000
#> 4 block_1_1000_25027        4 2022222002222002220020220
#> 5 block_1_1000_25027        5 1122122111221112211121211
#> 6 block_1_1000_25027        6 0000020022000020200020000
#>         nucleotide_sequence frequency
#> 1 SCWRRTYAKTSSWARYCMGACCWRG      0.22
#> 2 GSWRGTYMTWSSTRAYCAKWCSWAR      0.14
#> 3 GCWAGYTAKWGCTARCYAKWSSTAG      0.12
#> 4 CCAGATCAGTCGAAGTCCGACCAGG      0.12
#> 5 SSAGRTCMKWCGWRRTCMKWCSARR      0.11
#> 6 GCTAGTTATTGCTAACCAGACCTAG      0.10
```
