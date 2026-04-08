# Clique-Based LD Block Detection

Partitions SNPs in a genomic window into LD cliques. The LD matrix is
computed by the C++ Armadillo kernel; adjacency is built by
`build_adj_matrix_cpp`; cliques are found by igraph.

## Usage

``` r
CLQD(
  subgeno,
  subSNPinfo,
  adj_subgeno,
  CLQcut = 0.5,
  clstgap = 40000,
  CLQmode = c("Density", "Maximal"),
  codechange = FALSE,
  checkLargest = FALSE,
  split = FALSE,
  digits = -1L,
  n_threads = 1L,
  verbose = FALSE
)
```

## Arguments

- subgeno:

  Numeric matrix (individuals x SNPs), 0/1/2 coded.

- subSNPinfo:

  Data frame: col 1 = rsID, col 2 = bp position.

- adj_subgeno:

  n x p matrix from
  [`prepare_geno()`](https://FAkohoue.github.io/LDxBlocks/reference/prepare_geno.md).

- CLQcut:

  Numeric in (0,1\]. LD threshold for edges. Default 0.5.

- clstgap:

  Integer. Max bp gap within clique (split=TRUE).

- CLQmode:

  "Density" (default) or "Maximal".

- codechange:

  Logical. Sign-flip heuristic. Default FALSE.

- checkLargest:

  Logical. Dense-core pre-pass for \>= 500 SNPs.

- split:

  Logical. Split cliques over gaps. Default FALSE.

- digits:

  Integer. Rounding for LD. -1 = none (default).

- n_threads:

  Integer. OpenMP threads. Default 1.

- verbose:

  Logical.

## Value

Integer vector length p: clique ID per SNP, NA = singleton.
