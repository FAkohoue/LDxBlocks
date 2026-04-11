# Read Genotype Data into an LDxBlocks Backend

Reads genotype data from disk (or accepts an in-memory matrix) and
returns a unified backend object that all LDxBlocks functions can
consume. The format is auto-detected from the file extension unless
overridden with `format`.

## Usage

``` r
read_geno(
  path,
  format = NULL,
  snp_info = NULL,
  sample_ids = NULL,
  sep = ",",
  na_strings = c("NA", "N", "NN", "./.", ".", ""),
  gds_cache = NULL,
  verbose = FALSE
)
```

## Arguments

- path:

  Character. Path to the genotype file, OR an R numeric matrix when
  `format = "matrix"`.

- format:

  Character or `NULL`. One of `"numeric"`, `"hapmap"`, `"vcf"`, `"gds"`,
  `"bed"`, `"matrix"`. `NULL` (default) auto-detects from extension.

- snp_info:

  Data frame with columns `SNP`, `CHR`, `POS` (and optionally `REF`,
  `ALT`). Required only when `format = "matrix"`.

- sample_ids:

  Character vector. Override sample IDs extracted from the file. Length
  must equal number of samples.

- sep:

  Character. Field separator for `"numeric"` format.

- na_strings:

  Character vector. Strings treated as NA. Default
  `c("NA", "N", "NN", "./.", ".", "")`.

- gds_cache:

  Character or `NULL`. Path where a GDS cache file should be written
  when `format = "vcf"` and SNPRelate is available. If `NULL` (default),
  the GDS file is placed next to the VCF with a `.gds` extension. Set to
  `FALSE` to disable auto-conversion and read the VCF fully into memory
  instead. When the cache file already exists it is reused without
  re-converting (fast subsequent calls). Ignored for all non-VCF
  formats. Default `","`.

- verbose:

  Logical. Print progress messages.

## Value

An object of class `"LDxBlocks_backend"` with elements: `type`,
`n_samples`, `n_snps`, `sample_ids`, `snp_info` (data.frame
SNP/CHR/POS/REF/ALT). Use
[`read_chunk`](https://FAkohoue.github.io/LDxBlocks/reference/read_chunk.md)
to extract genotype slices and
[`close_backend`](https://FAkohoue.github.io/LDxBlocks/reference/close_backend.md)
to release file handles.

## Supported formats

- `"numeric"`:

  CSV or TXT, one row per SNP. Required columns: `SNP`, `CHR`, `POS`,
  `REF`, `ALT`. Remaining columns are sample dosages in {0, 1, 2, NA}.
  Extension: `.csv`, `.txt`.

- `"hapmap"`:

  Standard 11-column HapMap header followed by sample columns with
  two-character nucleotide calls (`AA`, `AT`, `NN` for missing).
  Extension: `.hmp.txt`.

- `"vcf"`:

  VCF v4.2. Both phased (`0|1`) and unphased (`0/1`) GT fields are
  accepted. Multi-allelic sites use first ALT. Missing (`./.`) becomes
  `NA`. Extension: `.vcf`, `.vcf.gz`.

- `"gds"`:

  SNPRelate GDS file. Requires `BiocManager::install("SNPRelate")`.
  Extension: `.gds`.

- `"bed"`:

  PLINK binary BED. Companion `.bim` and `.fam` files must exist at the
  same path stem. Extension: `.bed`.

- `"matrix"`:

  In-memory R numeric matrix (individuals x SNPs, 0/1/2). Must supply
  `snp_info` separately.

## See also

[`read_chunk`](https://FAkohoue.github.io/LDxBlocks/reference/read_chunk.md),
[`close_backend`](https://FAkohoue.github.io/LDxBlocks/reference/close_backend.md),
[`run_Big_LD_all_chr`](https://FAkohoue.github.io/LDxBlocks/reference/run_Big_LD_all_chr.md)

## Examples

``` r
# \donttest{
# In-memory matrix (existing workflow, unchanged)
G   <- matrix(sample(0:2, 100*500, replace=TRUE), 100, 500)
rownames(G) <- paste0("ind", 1:100)
colnames(G) <- paste0("rs",  1:500)
info <- data.frame(SNP=colnames(G), CHR=rep(1:5, each=100),
                   POS=rep(seq(1e4,5e6,length.out=100),5))
be <- read_geno(G, format="matrix", snp_info=info)
be$n_snps   # 500
#> [1] 500
mat <- read_chunk(be, 1:20)
dim(mat)    # 100 x 20
#> [1] 100  20
close_backend(be)
# }
```
