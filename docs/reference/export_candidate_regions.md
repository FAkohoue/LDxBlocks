# Export Candidate Gene Regions to BED, CSV, or biomaRt Format

Converts the output of
[`define_qtl_regions`](https://FAkohoue.github.io/LDxBlocks/reference/define_qtl_regions.md)
into formats ready for downstream annotation: standard BED (for
BEDtools, UCSC), a plain CSV, or a named list suitable for direct use
with `biomaRt::getBM()` filters.

## Usage

``` r
export_candidate_regions(
  qtl_regions,
  format = c("bed", "csv", "biomart"),
  out_file = NULL,
  chr_prefix = "",
  use_lead_snp = TRUE,
  padding_bp = 0L
)
```

## Arguments

- qtl_regions:

  Data frame from
  [`define_qtl_regions`](https://FAkohoue.github.io/LDxBlocks/reference/define_qtl_regions.md).

- format:

  Character. Output format: `"bed"`, `"csv"`, or `"biomart"`. Default
  `"bed"`.

- out_file:

  Character path. If `NULL` (default), returns the object invisibly
  without writing to disk.

- chr_prefix:

  Character. Prefix to add to chromosome names in BED output (e.g.
  `"chr"` for UCSC, `""` for Ensembl). Default `""`.

- use_lead_snp:

  Logical. If `TRUE` and `ld_decay` was supplied to
  `define_qtl_regions`, use the LD-extended
  `candidate_region_start`/`candidate_region_end` columns instead of
  block boundaries. Default `TRUE`.

- padding_bp:

  Integer. Additional bp to add on each side of each region (applied
  after the LD extension if any). Default `0L`.

## Value

Invisibly:

- `"bed"`:

  Data frame with columns `chrom`, `start` (0-based), `end`, `name`,
  `score`, `strand`.

- `"csv"`:

  The `qtl_regions` data frame, optionally written to `out_file`.

- `"biomart"`:

  Named list with elements `chromosome_name`, `start`, `end` suitable
  for
  `biomaRt::getBM(filters = c("chromosome_name","start","end"), ...)`.

## See also

[`define_qtl_regions`](https://FAkohoue.github.io/LDxBlocks/reference/define_qtl_regions.md),
[`compute_ld_decay`](https://FAkohoue.github.io/LDxBlocks/reference/compute_ld_decay.md)

## Examples

``` r
# \donttest{
data(ldx_geno, ldx_snp_info, ldx_blocks, ldx_gwas, package = "LDxBlocks")
qtl <- define_qtl_regions(ldx_gwas, ldx_blocks, ldx_snp_info,
                           p_threshold = NULL)
# BED file
bed <- export_candidate_regions(qtl, format = "bed", chr_prefix = "chr")
head(bed)
#>   chrom  start    end                  name score strand
#> 1  chr1    999  25027    block_1_1000_25027     3      .
#> 2  chr1  81063  99022   block_1_81064_99022     3      .
#> 3  chr1 155367 179371 block_1_155368_179371     3      .
#> 4  chr2    999  30023    block_2_1000_30023     4      .
#> 5  chr2  86235 105290  block_2_86236_105290     3      .
#> 6  chr3    999  19068    block_3_1000_19068     4      .
# biomaRt-ready list
bm  <- export_candidate_regions(qtl, format = "biomart")
# biomaRt::getBM(attributes = c("ensembl_gene_id","external_gene_name"),
#                filters = c("chromosome_name","start","end"),
#                values  = bm, mart = my_mart)
# }
```
