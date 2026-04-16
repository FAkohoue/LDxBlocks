# Plot LD Decay Curve

Plots r\\^2\\ vs physical distance from
[`compute_ld_decay`](https://FAkohoue.github.io/LDxBlocks/reference/compute_ld_decay.md).

## Usage

``` r
plot_ld_decay(
  x,
  plot_points = FALSE,
  plot_threshold = TRUE,
  plot_decay_dist = TRUE,
  max_dist_kb = NULL,
  facet = FALSE,
  ...
)
```

## Arguments

- x:

  `LDxBlocks_decay` object.

- plot_points:

  Logical. Plot raw pair r\\^2\\ as semi-transparent points. Default
  `FALSE`.

- plot_threshold:

  Logical. Horizontal line at critical r\\^2\\. Default `TRUE`.

- plot_decay_dist:

  Logical. Vertical lines at per-chromosome decay distances. Default
  `TRUE`.

- max_dist_kb:

  Numeric. X-axis limit in kb. Default `NULL` (data range).

- facet:

  Logical. One panel per chromosome. Default `FALSE`.

- ...:

  Ignored.

## Value

A `ggplot2` object.
