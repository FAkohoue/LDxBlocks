# Clique-Based LD Block Detection

Partitions SNPs in a genomic window into LD cliques. Three algorithms
are available via `CLQmode`:

- `"Density"` (default):

  Standard Big-LD clique scoring – edges built from the full r? matrix,
  cliques found by Bron-Kerbosch
  ([`igraph::max_cliques`](https://r.igraph.org/reference/cliques.html)),
  scored by density = size / span. Exact, but exponential worst-case.
  Use `checkLargest = TRUE` to guard against blowup on dense WGS panels.

- `"Maximal"`:

  Same algorithm as Density but scored by clique size only (original
  Big-LD maximality criterion).

- `"Louvain"`:

  Louvain community detection
  ([`igraph::cluster_louvain`](https://r.igraph.org/reference/cluster_louvain.html)).
  Polynomial time O(n log n) – never exponential. Recommended for WGS
  panels (\>500 SNPs per window) or any time `CLQmode = "Density"` is
  slow. Communities are returned as bin IDs matching the Density output
  format so all downstream code is unaffected.

- `"Leiden"`:

  Leiden community detection
  ([`igraph::cluster_leiden`](https://r.igraph.org/reference/cluster_leiden.html)).
  Stricter than Louvain – produces well-connected communities without
  disconnected nodes. Requires igraph \>= 1.3.0. Slightly slower than
  Louvain but higher quality communities. Recommended over Louvain when
  resolution matters.

When `max_bp_distance` is supplied (\> 0), only SNP pairs within that
physical distance have their r? computed (`compute_r2_sparse_cpp`). This
reduces the O(p?) LD computation to near-O(p) for WGS panels where
long-range LD is negligible (typically \> 500 kb). Applied before any
clique algorithm.

## Usage

``` r
CLQD(
  subgeno,
  subSNPinfo,
  adj_subgeno,
  CLQcut = 0.5,
  clstgap = 40000,
  CLQmode = c("Density", "Maximal", "Louvain", "Leiden"),
  codechange = FALSE,
  checkLargest = FALSE,
  split = FALSE,
  digits = -1L,
  n_threads = 1L,
  max_bp_distance = 0L,
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

  Algorithm: `"Density"` (default, exact, exponential worst-case),
  `"Maximal"` (exact, maximality criterion), `"Louvain"` (polynomial,
  recommended for WGS), `"Leiden"` (polynomial, stricter communities,
  igraph \>= 1.3.0).

- codechange:

  Logical. Sign-flip heuristic. Default FALSE.

- checkLargest:

  Logical. Dense-core pre-pass for \>= 500 SNPs when using Density or
  Maximal mode. Default FALSE. Has no effect when CLQmode is Louvain or
  Leiden (already polynomial).

- split:

  Logical. Split cliques over gaps. Default FALSE.

- digits:

  Integer. Rounding for LD. -1 = none (default).

- n_threads:

  Integer. OpenMP threads. Default 1.

- max_bp_distance:

  Integer. Maximum base-pair distance between a SNP pair for its r? to
  be computed. Pairs beyond this distance are assumed to be in
  negligible LD and set to 0. Default `0L` = disabled (compute all
  pairs, original behaviour). Recommended value for WGS panels:
  `500000L` (500 kb). Has no effect when the window spans less than
  `max_bp_distance`.

- verbose:

  Logical.

## Value

Integer vector length p: clique/community ID per SNP, NA = singleton.
