# Open a bigmemory-backed Genotype Store

Creates a `big.matrix` on disk from any supported genotype source and
wraps it in the standard `LDxBlocks_backend` interface. Subsequent calls
to
[`read_chunk()`](https://FAkohoue.github.io/LDxBlocks/reference/read_chunk.md)
retrieve columns via memory-mapped I/O – the OS pages in only the
requested bytes on demand, so peak RAM is proportional to the columns
accessed, not the full matrix.

This is useful when:

- The filtered genotype matrix (individuals x SNPs) is too large to hold
  in RAM (\> 4-8 GB for typical server configurations).

- The pipeline needs to be restarted after a previous run – the `.bin`
  and `.desc` files persist on disk and can be reused without re-loading
  the source VCF/GDS.

- Multiple R sessions or future workers need simultaneous read access to
  the same matrix (bigmemory's file-backed store is safe for concurrent
  reads).

## Usage

``` r
read_geno_bigmemory(
  source,
  snp_info = NULL,
  backingfile = tempfile("ldxbm_"),
  backingpath = tempdir(),
  type = "char",
  verbose = TRUE
)
```

## Arguments

- source:

  Either an `LDxBlocks_backend` object (any format), a plain R matrix,
  or a path to a previously saved `.desc` file (reattach without
  reloading).

- snp_info:

  Data frame with `SNP`, `CHR`, `POS`. Required when `source` is a plain
  matrix or a `.desc` file (bigmemory does not store metadata). Optional
  and ignored when `source` is a file path or an `LDxBlocks_backend` –
  SNP info is obtained automatically from those sources.

- backingfile:

  Character. Stem for the `.bin` and `.desc` files. Default: a tempfile.
  Supply a persistent path to reuse across sessions.

- backingpath:

  Character. Directory for backing files. Default: tempdir().

- type:

  Storage type: `"char"` (1 byte per cell, values 0-2 fit, saves 8x vs
  double), `"short"` (2 bytes), or `"double"`. Default `"char"`.

- verbose:

  Logical. Default `TRUE`.

## Value

An `LDxBlocks_backend` object with `type = "bigmemory"`. Use
`read_chunk(be, col_idx)` and `close_backend(be)` as normal.

## Memory model

`big.matrix` stores the matrix as a raw binary file (`.bin`) with a
companion descriptor (`.desc`). The OS memory-maps the file:
`read_chunk(be, col_idx)` calls `bigmemory::as.matrix(bm[, col_idx])`
which triggers page faults that load only the requested column pages.
This is equivalent to the GDS streaming model but works for any input
format and avoids repeated `snpgdsGetGeno()` calls.

## Examples

``` r
if (FALSE) { # \dontrun{
# Convert a GDS backend to a persistent bigmemory store
be_gds <- read_geno("mydata.gds")
be_bm  <- read_geno_bigmemory(be_gds,
                              backingfile = "mydata_bm",
                              backingpath = "/data/ldxblocks")
close_backend(be_gds)

# All subsequent runs reattach without reloading
be_bm2 <- read_geno_bigmemory("/data/ldxblocks/mydata_bm.desc")
blocks  <- run_Big_LD_all_chr(be_bm2, CLQcut = 0.70)
close_backend(be_bm2)

# From a plain matrix
be_bm <- read_geno_bigmemory(ldx_geno, snp_info = ldx_snp_info)
} # }
```
