# Close an LDxBlocks Backend and Release File Handles

Closes any open file connections held by the backend. For in-memory
backends (`"matrix"`, `"numeric"`, `"hapmap"`, `"vcf"`) this is a no-op.
For `"gds"` backends it calls
[`SeqArray::seqClose()`](https://rdrr.io/pkg/SeqArray/man/seqClose.html).
For `"bed"` backends the memory-mapped file is released.

## Usage

``` r
close_backend(backend)
```

## Arguments

- backend:

  An `"LDxBlocks_backend"` object.

## Value

Invisibly `NULL`.
