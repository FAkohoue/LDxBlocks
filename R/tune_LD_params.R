# Unexported helpers ------------------------------------------------------
#' @noRd
.build_blocks_chr <- function(geno_chr, info_chr, params) {
  Big_LD(
    geno        = geno_chr,
    SNPinfo     = as.data.frame(info_chr),
    CLQcut      = params$CLQcut,
    clstgap     = params$clstgap,
    leng        = params$leng,
    subSegmSize = params$subSegmSize,
    MAFcut      = params$MAFcut,
    appendrare  = params$appendrare,
    checkLargest= params$checkLargest,
    CLQmode     = params$CLQmode,
    rV2method   = params$rV2method,
    split       = params$split,
    digits      = params$digits,
    seed        = params$seed %||% NULL,
    verbose     = FALSE
  )
}
`%||%` <- function(a, b) if (!is.null(a)) a else b

#' @noRd
.consolidate_blocks <- function(ld_blocks_by_chr) {
  if (length(ld_blocks_by_chr) == 0) return(NULL)
  all_blocks <- data.table::rbindlist(ld_blocks_by_chr, use.names = TRUE, fill = TRUE)
  stopifnot(all(c("start","end","start.rsID","end.rsID","start.bp","end.bp","CHR") %in% names(all_blocks)))
  all_blocks <- all_blocks[order(CHR, start.bp)]
  all_blocks[, block_idx := sequence(.N), by = CHR]
  all_blocks[, block_name := paste0("qtl", CHR, ".", block_idx)]
  all_blocks[, length_snps := end - start + 1L]
  all_blocks[, length_bp   := end.bp - start.bp + 1L]
  all_blocks[]
}

#' @noRd
.assign_gwas <- function(all_blocks, gwas_df) {
  if (is.null(all_blocks) || nrow(all_blocks) == 0L) {
    gwas_df$LD_block <- NA_character_
    return(list(gwas = gwas_df, n_unassigned = nrow(gwas_df), n_forced = 0L))
  }
  gwas <- gwas_df
  gwas$LD_block <- NA_character_

  idx_list <- lapply(seq_len(nrow(gwas)), function(i) {
    chr <- as.character(gwas$CHR[i]); pos <- as.numeric(gwas$POS[i])
    which(all_blocks$CHR == chr & all_blocks$start.bp <= pos & all_blocks$end.bp >= pos)
  })
  gwas$LD_block <- vapply(idx_list, function(idx) {
    if (length(idx) == 0) NA_character_ else paste(unique(all_blocks$block_name[idx]), collapse = ",")
  }, character(1))

  unassigned <- which(is.na(gwas$LD_block))
  n_forced <- 0L
  if (length(unassigned) > 0) {
    for (i in unassigned) {
      chr <- as.character(gwas$CHR[i]); pos <- as.numeric(gwas$POS[i])
      chr_blocks <- all_blocks[all_blocks$CHR == chr, ]
      if (nrow(chr_blocks) == 0) {
        gwas$LD_block[i] <- NA_character_
      } else {
        dists <- pmax(chr_blocks$start.bp - pos, pos - chr_blocks$end.bp)
        dists[(pos >= chr_blocks$start.bp) & (pos <= chr_blocks$end.bp)] <- 0
        nearest <- which.min(dists)
        gwas$LD_block[i] <- paste0(chr_blocks$block_name[nearest], "*")
        n_forced <- n_forced + 1L
      }
    }
  }
  list(gwas = gwas, n_unassigned = sum(is.na(gwas$LD_block)), n_forced = n_forced)
}

#' @noRd
.score_combo <- function(assign_res, all_blocks, target_bp = c(5e4, 5e5)) {
  if (is.null(all_blocks) || nrow(all_blocks) == 0L) {
    return(list(n_unassigned = Inf, n_forced = Inf, n_blocks = Inf, med_bp = Inf, penalty_bp = Inf))
  }
  med_bp <- stats::median(all_blocks$length_bp, na.rm = TRUE)
  L <- target_bp[1]; U <- target_bp[2]
  penalty_bp <- if (is.na(med_bp)) Inf else if (med_bp < L) L - med_bp else if (med_bp > U) med_bp - U else 0
  list(n_unassigned = assign_res$n_unassigned, n_forced = assign_res$n_forced,
       n_blocks = nrow(all_blocks), med_bp = med_bp, penalty_bp = penalty_bp)
}

# Exported tuner -----------------------------------------------------------

#' Auto-tune LD Parameters to Minimize Unassigned GWAS Markers
#'
#' @description
#' Grid-search parameter combinations for \code{Big_LD()} and select the set
#' that minimizes unassigned GWAS markers, then forced assignments, then number
#' of blocks, then deviation from a target median block size band. Optionally
#' prioritize perfect (0/0) solutions.
#'
#' @param geno_matrix Numeric matrix (individuals x SNPs; 0/1/2).
#' @param snp_info Data frame with columns \code{SNP}, \code{CHR}, \code{POS}.
#' @param gwas_df Data frame with columns \code{Marker}, \code{CHR}, \code{POS}.
#' @param grid Optional data frame of parameter combos; if \code{NULL}, a default grid is used.
#' @param chromosomes Optional character vector of CHR to tune on; \code{NULL}=all.
#' @param target_bp_band Preferred median block size band (bp), length-2 numeric.
#' @param parallel If TRUE, uses \code{future.apply::future_lapply} (user sets plan()).
#' @param seed Integer seed (default \code{NULL}); if supplied, sets RNG seed.
#' @param prefer_perfect If TRUE, select among 0/0 combos first (tie-break by n_blocks, penalty).
#' @param return_all_perfect If TRUE, include a table of all perfect combos in the result.
#'
#' @return List with \code{best_params}, \code{score_table}, \code{perfect_table},
#'   \code{final_blocks}, \code{gwas_assigned}.
#' @export
tune_LD_params <- function(
    geno_matrix,
    snp_info,
    gwas_df,
    grid = NULL,
    chromosomes = NULL,
    target_bp_band = c(5e4, 5e5),
    parallel = FALSE,
    seed = NULL,
    prefer_perfect = TRUE,
    return_all_perfect = TRUE
) {
  if (!is.null(seed)) set.seed(seed)

  req1 <- c("SNP","CHR","POS"); req2 <- c("Marker","CHR","POS")
  if (!all(req1 %in% names(snp_info))) stop("snp_info must contain: ", paste(req1, collapse=", "))
  if (!all(req2 %in% names(gwas_df)))  stop("gwas_df must contain: ", paste(req2, collapse=", "))

  if (parallel && !requireNamespace("future.apply", quietly = TRUE)) {
    stop("parallel=TRUE requires the 'future.apply' package. Install it or set parallel=FALSE.")
  }

  if (is.null(grid)) {
    grid <- expand.grid(
      CLQcut       = c(0.65, 0.70, 0.75, 0.80),
      clstgap      = c(2e6),
      leng         = c(1000),
      subSegmSize  = c(10000),
      split        = c(TRUE),
      checkLargest = c(FALSE),
      MAFcut       = c(0.05),
      CLQmode      = c("Density"),
      rV2method    = c("chol"),
      digits       = c(6),
      appendrare   = c(FALSE),
      stringsAsFactors = FALSE
    )
  }

  snp_info$CHR <- as.character(snp_info$CHR)
  if (is.null(chromosomes)) chromosomes <- unique(snp_info$CHR)
  chromosomes <- as.character(chromosomes)

  chr_index <- lapply(chromosomes, function(chr) which(snp_info$CHR == chr))
  names(chr_index) <- chromosomes

  eval_one_combo <- function(row_idx) {
    params <- as.list(grid[row_idx, , drop = FALSE])
    ld_by_chr <- vector("list", length(chromosomes)); names(ld_by_chr) <- chromosomes
    for (k in seq_along(chromosomes)) {
      chr <- chromosomes[k]
      idx <- chr_index[[k]]
      if (length(idx) < 10) next
      geno_chr <- geno_matrix[, idx, drop = FALSE]
      info_chr <- snp_info[idx, c("SNP","POS")]
      if (!is.numeric(info_chr$POS)) info_chr$POS <- as.numeric(info_chr$POS)
      blk <- try(.build_blocks_chr(geno_chr, info_chr, params), silent = TRUE)
      if (!inherits(blk, "try-error") && !is.null(blk)) { blk$CHR <- chr; ld_by_chr[[k]] <- blk }
    }
    ld_blocks  <- .consolidate_blocks(ld_by_chr)
    assign_res <- .assign_gwas(ld_blocks, gwas_df)
    score      <- .score_combo(assign_res, ld_blocks, target_bp = target_bp_band)

    data.frame(
      row             = row_idx,
      CLQcut          = params$CLQcut,
      clstgap         = params$clstgap,
      leng            = params$leng,
      subSegmSize     = params$subSegmSize,
      split           = params$split,
      checkLargest    = params$checkLargest,
      MAFcut          = params$MAFcut,
      CLQmode         = params$CLQmode,
      rV2method       = params$rV2method,
      digits          = params$digits,
      appendrare      = params$appendrare,
      n_unassigned    = assign_res$n_unassigned,
      n_forced        = assign_res$n_forced,
      n_blocks        = if (is.null(ld_blocks)) NA_integer_ else nrow(ld_blocks),
      median_block_bp = score$med_bp,
      penalty_bp      = score$penalty_bp,
      stringsAsFactors = FALSE
    )
  }

  scores_list <- if (parallel) {
    future.apply::future_lapply(seq_len(nrow(grid)), eval_one_combo)
  } else {
    lapply(seq_len(nrow(grid)), eval_one_combo)
  }
  scores <- dplyr::bind_rows(scores_list)
  if (!nrow(scores)) stop("No parameter combinations could be evaluated successfully.")

  perfect <- if (isTRUE(prefer_perfect)) dplyr::filter(scores, n_unassigned == 0L, n_forced == 0L) else NULL
  if (!is.null(perfect) && nrow(perfect) > 0) {
    best <- perfect |>
      dplyr::arrange(n_blocks, penalty_bp) |>
      dplyr::slice(1)
  } else {
    best <- scores |>
      dplyr::arrange(n_unassigned, n_forced, n_blocks, penalty_bp) |>
      dplyr::slice(1)
  }

  best_params <- as.list(best[1, c(
    "CLQcut","clstgap","leng","subSegmSize","split",
    "checkLargest","MAFcut","CLQmode","rV2method",
    "digits","appendrare"
  )])

  all_chr <- unique(as.character(snp_info$CHR))
  final_ld_by_chr <- vector("list", length(all_chr)); names(final_ld_by_chr) <- all_chr
  for (chr in all_chr) {
    idx <- which(snp_info$CHR == chr)
    if (length(idx) < 10) next
    geno_chr <- geno_matrix[, idx, drop = FALSE]
    info_chr <- snp_info[idx, c("SNP","POS")]
    if (!is.numeric(info_chr$POS)) info_chr$POS <- as.numeric(info_chr$POS)
    blk <- try(.build_blocks_chr(geno_chr, info_chr, best_params), silent = TRUE)
    if (!inherits(blk, "try-error") && !is.null(blk)) { blk$CHR <- chr; final_ld_by_chr[[chr]] <- blk }
  }
  final_blocks <- .consolidate_blocks(final_ld_by_chr)
  final_assign <- .assign_gwas(final_blocks, gwas_df)

  list(
    best_params    = best_params,
    score_table    = scores,
    perfect_table  = if (isTRUE(return_all_perfect) && !is.null(perfect) && nrow(perfect) > 0) perfect else NULL,
    final_blocks   = final_blocks,
    gwas_assigned  = final_assign$gwas
  )
}
