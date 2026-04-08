# ─────────────────────────────────────────────────────────────────────────────
# utils.R  –  Summary statistics and visualisation helpers
# ─────────────────────────────────────────────────────────────────────────────


#' Summarise LD Block Characteristics
#'
#' @description
#' Computes per-chromosome and genome-wide summary statistics for a block table
#' returned by \code{\link{run_Big_LD_all_chr}} or \code{\link{Big_LD}}.
#'
#' @param blocks Data frame of LD blocks. Must contain at least
#'   \code{start.bp}, \code{end.bp}. If a \code{CHR} column is present,
#'   per-chromosome summaries are produced.
#'
#' @return A \code{data.frame} with one row per chromosome (plus one row for
#'   the genome-wide summary) and columns:
#'   \describe{
#'     \item{\code{CHR}}{Chromosome (or \code{"GENOME"}).}
#'     \item{\code{n_blocks}}{Number of blocks.}
#'     \item{\code{min_bp}, \code{median_bp}, \code{mean_bp}, \code{max_bp}}{
#'       Block size distribution in base pairs.}
#'     \item{\code{total_bp_covered}}{Sum of block sizes (bp).}
#'   }
#'
#' @examples
#' \donttest{
#' blocks <- data.frame(
#'   CHR      = c("chr1","chr1","chr2"),
#'   start.bp = c(1000, 500000, 2000),
#'   end.bp   = c(80000, 650000, 200000)
#' )
#' summarise_blocks(blocks)
#' }
#'
#' @export
summarise_blocks <- function(blocks) {
  if (!all(c("start.bp", "end.bp") %in% names(blocks)))
    stop("blocks must contain 'start.bp' and 'end.bp' columns.")

  blocks$length_bp <- as.numeric(blocks$end.bp) - as.numeric(blocks$start.bp) + 1L

  summarise_one <- function(df, chr_label) {
    len <- df$length_bp
    data.frame(
      CHR             = chr_label,
      n_blocks        = nrow(df),
      min_bp          = min(len),
      median_bp       = stats::median(len),
      mean_bp         = mean(len),
      max_bp          = max(len),
      total_bp_covered = sum(len),
      stringsAsFactors = FALSE
    )
  }

  if ("CHR" %in% names(blocks)) {
    chrs    <- unique(as.character(blocks$CHR))
    per_chr <- lapply(chrs, function(chr) {
      summarise_one(blocks[blocks$CHR == chr, ], chr)
    })
    genome  <- summarise_one(blocks, "GENOME")
    do.call(rbind, c(per_chr, list(genome)))
  } else {
    summarise_one(blocks, "ALL")
  }
}


#' Plot LD Block Structure Across Chromosomes
#'
#' @description
#' Produces a \code{ggplot2}-based overview of LD block boundaries and sizes.
#' Requires \code{ggplot2} to be installed.
#'
#' @param blocks Data frame of LD blocks (output of
#'   \code{\link{run_Big_LD_all_chr}}). Must contain \code{start.bp},
#'   \code{end.bp}, and \code{CHR}.
#' @param colour_by Character. One of \code{"length_bp"} (default, colours by
#'   block size) or \code{"CHR"} (distinct colour per chromosome).
#' @param highlight_blocks Optional character vector of block IDs (from
#'   \code{block_name} column) to highlight with a border. \code{NULL}
#'   (default) highlights nothing.
#' @param mb_scale Logical. If \code{TRUE} (default), x-axis is in Megabases.
#'
#' @return A \code{ggplot} object.
#'
#' @examples
#' \donttest{
#' if (requireNamespace("ggplot2", quietly = TRUE)) {
#'   blocks <- data.frame(
#'     CHR = c("chr1","chr1","chr2","chr2"),
#'     start.bp = c(1e5, 6e5, 2e5, 8e5),
#'     end.bp   = c(5e5, 9e5, 7e5, 1.5e6)
#'   )
#'   p <- plot_ld_blocks(blocks)
#'   print(p)
#' }
#' }
#'
#' @export
plot_ld_blocks <- function(
    blocks,
    colour_by        = c("length_bp", "CHR"),
    highlight_blocks = NULL,
    mb_scale         = TRUE
) {
  if (!requireNamespace("ggplot2", quietly = TRUE))
    stop("Package 'ggplot2' is required for plot_ld_blocks(). Install it with: install.packages('ggplot2')")

  colour_by <- match.arg(colour_by)
  blocks    <- as.data.frame(blocks)

  if (!all(c("start.bp", "end.bp", "CHR") %in% names(blocks)))
    stop("blocks must contain 'start.bp', 'end.bp', and 'CHR'.")

  blocks$length_bp  <- as.numeric(blocks$end.bp) - as.numeric(blocks$start.bp) + 1L
  blocks$CHR        <- factor(as.character(blocks$CHR),
                               levels = mixsort(unique(as.character(blocks$CHR))))

  scale_div <- if (mb_scale) 1e6 else 1
  x_label   <- if (mb_scale) "Position (Mb)" else "Position (bp)"

  blocks$start_x <- as.numeric(blocks$start.bp) / scale_div
  blocks$end_x   <- as.numeric(blocks$end.bp)   / scale_div
  blocks$chr_int <- as.integer(blocks$CHR)

  gg <- ggplot2::ggplot(blocks) +
    ggplot2::geom_rect(
      ggplot2::aes(xmin = start_x, xmax = end_x,
                   ymin = chr_int - 0.4, ymax = chr_int + 0.4,
                   fill = if (colour_by == "length_bp") length_bp else CHR),
      colour = NA, alpha = 0.85
    ) +
    ggplot2::scale_y_continuous(
      breaks = seq_along(levels(blocks$CHR)),
      labels = levels(blocks$CHR)
    ) +
    ggplot2::labs(
      x    = x_label,
      y    = "Chromosome",
      fill = if (colour_by == "length_bp") "Block size (bp)" else "Chromosome",
      title = "LD Block Structure"
    ) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      panel.grid.major.y = ggplot2::element_blank(),
      panel.grid.minor   = ggplot2::element_blank(),
      legend.position    = "right"
    )

  if (colour_by == "length_bp") {
    gg <- gg + ggplot2::scale_fill_viridis_c(option = "plasma", trans = "log10")
  }
  gg
}


# ── Internal helper: natural sort of chromosome names ────────────────────────
mixsort <- function(x) {
  x[order(
    suppressWarnings(as.integer(gsub("[^0-9]", "", x))),
    x
  )]
}
