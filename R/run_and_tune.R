# -----------------------------------------------------------------------------
# run_Big_LD_all_chr.R  -  Chromosome-wise wrapper
# -----------------------------------------------------------------------------

#' Genome-Wide LD Block Detection by Chromosome
#'
#' @description
#' Applies \code{Big_LD()} chromosome by chromosome, collects results and
#' returns a single tidy data frame annotated with chromosome and block length.
#' This is the recommended entry point for genome-wide analyses.
#'
#' @param geno_matrix Numeric matrix (individuals x SNPs; values 0/1/2),
#'   spanning all chromosomes.
#' @param snp_info Data frame with columns \code{CHR}, \code{SNP}, \code{POS}.
#'   Column order is flexible; columns are matched by name.
#' @param method Character. LD metric: \code{"r2"} (default, standard squared
#'   Pearson correlation) or \code{"rV2"} (kinship-adjusted). See
#'   \code{Big_LD()} for details.
#' @param n_threads Integer. Number of OpenMP threads for the C++ LD kernel.
#'   Default \code{1L}. Increase for multi-core systems.
#' @param min_snps_chr Integer. Chromosomes with fewer SNPs than this after
#'   MAF filtering are skipped. Default \code{10L}. Increase to skip small
#'   scaffolds in whole-genome datasets.
#' @param chr Character vector or \code{NULL}. If supplied, only the named
#'   chromosomes are processed. Labels must match the values in
#'   \code{snp_info$CHR} after normalisation (no \code{chr} prefix). E.g.
#'   \code{chr = c("1","3","5")} or \code{chr = "X"}. Default \code{NULL}
#'   processes all chromosomes.
#' @param clean_malformed Logical. If \code{TRUE}, stream-clean the input file
#'   before reading by removing lines whose column count does not match the
#'   header. Only relevant when \code{geno_matrix} is a file path wrapped into
#'   a backend. Default \code{FALSE}.
#' @param singleton_as_block Logical. If \code{TRUE}, SNPs that pass MAF
#'   filtering but are not assigned to any clique are returned as single-SNP
#'   blocks (\code{start == end}, \code{length_bp == 1}). Default \code{FALSE}.
#'   See \code{Big_LD()} for details.
#' @param max_bp_distance Integer. Maximum base-pair distance between a SNP
#'   pair for its r\eqn{^2} to be computed. Pairs beyond this distance are set
#'   to zero in the adjacency matrix (assumed to be in negligible LD).
#'   \code{0L} (default) disables this and computes all pairs (original
#'   behaviour). Recommended value for WGS panels: \code{500000L} (500 kb).
#'   Has no effect when \code{CLQmode} is \code{"Louvain"} or
#'   \code{"Leiden"} and the window spans less than
#'   \code{max_bp_distance}. Requires sorted SNP positions within each
#'   sub-segment (guaranteed by \code{run_Big_LD_all_chr}).
#' @param CLQcut,clstgap,leng,subSegmSize,MAFcut,appendrare,checkLargest,CLQmode,kin_method,split,digits,seed,verbose
#'   Forwarded to \code{Big_LD()}. See that function's documentation for
#'   details.
#'
#' @return A \code{data.frame} with columns:
#'   \code{start}, \code{end}, \code{start.rsID}, \code{end.rsID},
#'   \code{start.bp}, \code{end.bp}, \code{CHR}, \code{length_bp}.
#'   Rows are sorted by \code{CHR} then \code{start.bp}.
#'
#' @seealso \code{Big_LD()}, \code{\link{tune_LD_params}},
#'   \code{\link{extract_haplotypes}}
#'
#' @examples
#' \donttest{
#' # Use the package example data -- 120 individuals, 230 SNPs, 3 chromosomes,
#' # 9 simulated LD blocks (3 per chromosome).
#' data(ldx_geno,     package = "LDxBlocks")
#' data(ldx_snp_info, package = "LDxBlocks")
#' blocks <- run_Big_LD_all_chr(
#'   ldx_geno, ldx_snp_info,
#'   method = "r2", CLQcut = 0.55, leng = 15L,
#'   subSegmSize = 100L, verbose = FALSE
#' )
#' head(blocks)
#' summarise_blocks(blocks)
#' }
#'
#' @export
run_Big_LD_all_chr <- function(
    geno_matrix,
    snp_info        = NULL,
    CLQcut          = 0.5,
    method          = c("r2", "rV2"),
    n_threads       = 1L,
    clstgap         = 40000,
    leng            = 200,
    subSegmSize     = 1500,
    MAFcut          = 0.05,
    appendrare      = FALSE,
    singleton_as_block = FALSE,
    checkLargest    = FALSE,
    CLQmode         = "Density",
    kin_method      = "chol",
    split           = FALSE,
    digits          = -1L,
    seed            = NULL,
    min_snps_chr    = 10L,
    chr             = NULL,
    clean_malformed = FALSE,
    max_bp_distance = 0L,
    verbose         = FALSE
) {
  # -- Accept either a plain matrix+snp_info OR an LDxBlocks_backend ----------
  if (inherits(geno_matrix, "LDxBlocks_backend")) {
    backend  <- geno_matrix
    snp_info <- backend$snp_info
  } else {
    # Wrap plain matrix into backend for uniform downstream code
    if (is.null(snp_info))
      stop("snp_info must be supplied when geno_matrix is a plain matrix.")
    req_cols <- c("CHR", "SNP", "POS")
    if (!all(req_cols %in% names(snp_info)))
      stop("snp_info must contain columns: ", paste(req_cols, collapse = ", "))
    backend <- read_geno(geno_matrix, format = "matrix", snp_info = snp_info,
                         clean_malformed = clean_malformed)
  }

  chromosomes <- unique(as.character(snp_info$CHR))

  # -- Optional: restrict to user-specified chromosomes ----------------------
  if (!is.null(chr)) {
    chr_req <- as.character(chr)
    missing_chr <- setdiff(chr_req, chromosomes)
    if (length(missing_chr))
      warning("Chromosomes not found in data and will be skipped: ",
              paste(missing_chr, collapse = ", "), call. = FALSE)
    chromosomes <- intersect(chr_req, chromosomes)
    if (!length(chromosomes))
      stop("None of the requested chromosomes exist in the data.", call. = FALSE)
  }

  ld_blocks_all <- vector("list", length(chromosomes))
  names(ld_blocks_all) <- chromosomes

  for (chr_i in chromosomes) {
    if (isTRUE(verbose)) cat("\n[run_Big_LD_all_chr] Processing", chr_i, "...\n")
    idx      <- which(snp_info$CHR == chr_i)
    geno_chr <- read_chunk(backend, idx)           # works for ALL backend types
    info_chr <- as.data.frame(snp_info[idx, c("SNP", "POS")])
    if (!is.numeric(info_chr$POS)) info_chr$POS <- as.numeric(info_chr$POS)

    if (ncol(geno_chr) < as.integer(min_snps_chr)) {
      if (isTRUE(verbose)) cat("  Skipping", chr_i, "- fewer than", min_snps_chr, "SNPs.\n")
      next
    }

    blk <- tryCatch(
      Big_LD(geno = geno_chr, SNPinfo = info_chr,
             CLQcut = CLQcut, clstgap = clstgap, leng = leng,
             subSegmSize = subSegmSize, MAFcut = MAFcut,
             appendrare = appendrare, singleton_as_block = singleton_as_block,
             checkLargest = checkLargest,
             CLQmode = CLQmode, kin_method = kin_method,
             split = split, method = method[1], n_threads = n_threads,
             digits = digits, seed = seed,
             max_bp_distance = max_bp_distance, verbose = verbose),
      error = function(e) {
        if (isTRUE(verbose))
          cat("  Error in Big_LD for", chr_i, ":", conditionMessage(e), "\n")
        NULL
      }
    )
    if (!is.null(blk) && nrow(blk) > 0L) { blk$CHR <- chr_i; ld_blocks_all[[chr_i]] <- blk }

    # Free chromosome genotype matrix immediately after use.
    # gc(FALSE) is called to release allocator pressure without a full GC cycle.
    # For GDS/BED backends geno_chr is already a fresh extraction; for matrix
    # backends it is a column slice. Either way it can be freed here.
    rm(geno_chr, info_chr); gc(FALSE)
  }

  all_blocks <- data.table::rbindlist(ld_blocks_all, use.names = TRUE, fill = TRUE)
  if (nrow(all_blocks)) {
    all_blocks$start.rsID <- as.character(all_blocks$start.rsID)
    all_blocks$end.rsID   <- as.character(all_blocks$end.rsID)
    all_blocks$length_bp  <- all_blocks$end.bp - all_blocks$start.bp + 1L
    data.table::setorder(all_blocks, CHR, start.bp)
  }
  as.data.frame(all_blocks)
}


# -----------------------------------------------------------------------------
# tune_LD_params.R  -  Grid-search auto-tuner
# -----------------------------------------------------------------------------

# -- Unexported helpers --------------------------------------------------------
`%||%` <- function(a, b) if (!is.null(a)) a else b

.build_blocks_chr <- function(geno_chr, info_chr, params) {
  Big_LD(
    geno         = geno_chr,
    SNPinfo      = as.data.frame(info_chr),
    CLQcut       = params$CLQcut,
    clstgap      = params$clstgap,
    leng         = params$leng,
    subSegmSize  = params$subSegmSize,
    MAFcut       = params$MAFcut,
    appendrare   = params$appendrare,
    checkLargest = params$checkLargest,
    CLQmode      = params$CLQmode,
    kin_method = params$kin_method,
    split        = params$split,
    digits       = params$digits,
    seed         = params$seed %||% NULL,
    verbose      = FALSE
  )
}

.consolidate_blocks <- function(ld_blocks_by_chr) {
  if (length(ld_blocks_by_chr) == 0L) return(NULL)
  ld_blocks_by_chr <- Filter(function(x) !is.null(x) && nrow(x) > 0L, ld_blocks_by_chr)
  if (length(ld_blocks_by_chr) == 0L) return(NULL)
  all_blocks <- data.table::rbindlist(ld_blocks_by_chr, use.names = TRUE, fill = TRUE)
  if (nrow(all_blocks) == 0L) return(NULL)
  req_cols   <- c("start", "end", "start.rsID", "end.rsID", "start.bp", "end.bp", "CHR")
  if (!all(req_cols %in% names(all_blocks))) return(NULL)
  all_blocks <- all_blocks[order(CHR, start.bp)]
  all_blocks[, block_idx  := sequence(.N), by = CHR]
  all_blocks[, block_name := paste0("qtl", CHR, ".", block_idx)]
  all_blocks[, length_snps := end - start + 1L]
  all_blocks[, length_bp   := end.bp - start.bp + 1L]
  all_blocks[]
}

.assign_gwas <- function(all_blocks, gwas_df) {
  if (is.null(all_blocks) || nrow(all_blocks) == 0L) {
    gwas_df$LD_block <- NA_character_
    return(list(gwas = gwas_df, n_unassigned = nrow(gwas_df), n_forced = 0L))
  }
  gwas          <- gwas_df
  gwas$LD_block <- NA_character_

  idx_list <- lapply(seq_len(nrow(gwas)), function(i) {
    chr <- as.character(gwas$CHR[i])
    pos <- as.numeric(gwas$POS[i])
    which(all_blocks$CHR == chr &
            all_blocks$start.bp <= pos &
            all_blocks$end.bp   >= pos)
  })
  gwas$LD_block <- vapply(idx_list, function(idx) {
    if (length(idx) == 0L) NA_character_
    else paste(unique(all_blocks$block_name[idx]), collapse = ",")
  }, character(1L))

  unassigned <- which(is.na(gwas$LD_block))
  n_forced   <- 0L
  for (i in unassigned) {
    chr        <- as.character(gwas$CHR[i])
    pos        <- as.numeric(gwas$POS[i])
    chr_blocks <- all_blocks[all_blocks$CHR == chr, ]
    if (nrow(chr_blocks) == 0L) next
    dists      <- pmax(chr_blocks$start.bp - pos, pos - chr_blocks$end.bp)
    dists[(pos >= chr_blocks$start.bp) & (pos <= chr_blocks$end.bp)] <- 0L
    nearest    <- which.min(dists)
    gwas$LD_block[i] <- paste0(chr_blocks$block_name[nearest], "*")
    n_forced   <- n_forced + 1L
  }
  list(gwas         = gwas,
       n_unassigned = sum(is.na(gwas$LD_block)),
       n_forced     = n_forced)
}

.score_combo <- function(assign_res, all_blocks, target_bp = c(5e4, 5e5)) {
  if (is.null(all_blocks) || nrow(all_blocks) == 0L)
    return(list(n_unassigned = Inf, n_forced = Inf,
                n_blocks = Inf, med_bp = Inf, penalty_bp = Inf))
  med_bp     <- stats::median(all_blocks$length_bp, na.rm = TRUE)
  L <- target_bp[1L]; U <- target_bp[2L]
  penalty_bp <- if (is.na(med_bp)) Inf
  else if (med_bp < L) L - med_bp
  else if (med_bp > U) med_bp - U
  else 0
  list(n_unassigned = assign_res$n_unassigned,
       n_forced     = assign_res$n_forced,
       n_blocks     = nrow(all_blocks),
       med_bp       = med_bp,
       penalty_bp   = penalty_bp)
}


#' Auto-Tune LD Block Detection Parameters
#'
#' @description
#' Performs a grid search over \code{Big_LD()} parameters and selects the
#' combination that minimises, in order of priority:
#' \enumerate{
#'   \item Unassigned GWAS markers (markers not falling in any block).
#'   \item Forced assignments (nearest-block fall-back).
#'   \item Number of blocks (parsimony).
#'   \item Deviation from \code{target_bp_band} (biological plausibility).
#' }
#' If \code{prefer_perfect = TRUE} (default), combinations achieving zero
#' unassigned and zero forced assignments are prioritised among the above.
#'
#' After selecting the best parameter set, \code{tune_LD_params} runs
#' \code{\link{run_Big_LD_all_chr}} on all chromosomes and assigns every GWAS
#' marker to a block, returning the final blocks and assignments.
#'
#' @param geno_matrix Numeric matrix (individuals x SNPs; 0/1/2), genome-wide.
#' @param snp_info Data frame with columns \code{SNP}, \code{CHR}, \code{POS}.
#' @param gwas_df Data frame with columns \code{Marker}, \code{CHR}, \code{POS}.
#' @param grid Optional data frame of parameter combinations. Each row is one
#'   combination; columns must match parameter names of \code{Big_LD}. If
#'   \code{NULL} (default), a sensible grid over \code{CLQcut} (4 values)
#'   and \code{min_freq} (2 values) is used, giving 8 combinations. Both
#'   are treated as hyperparameters following Weber SE et al. (2023,
#'   \emph{Front. Plant Sci.} \strong{14}:1217589,
#'   \doi{10.3389/fpls.2023.1217589}), who show that no single threshold
#'   is universally optimal across datasets and traits.
#' @param chromosomes Optional character vector of chromosome names to include
#'   in tuning. \code{NULL} uses all chromosomes in \code{snp_info}.
#' @param target_bp_band Length-2 numeric vector: preferred median block size
#'   range in base pairs. Default \code{c(5e4, 5e5)} (50 kb - 500 kb).
#' @param parallel Logical. If \code{TRUE}, uses
#'   \code{future.apply::future_lapply} for parallelism (user must set a
#'   \code{future} plan before calling). Default \code{FALSE}.
#' @param seed Integer seed for reproducibility. Default \code{NULL}.
#' @param prefer_perfect Logical. Give priority to parameter sets with zero
#'   unassigned / zero forced GWAS markers. Default \code{TRUE}.
#' @param return_all_perfect Logical. Include a table of all zero-zero
#'   combinations in the returned list. Default \code{TRUE}.
#'
#' @return A named list:
#'   \describe{
#'     \item{\code{best_params}}{Named list of the selected parameters.}
#'     \item{\code{score_table}}{Data frame of all grid combinations and their
#'       scores.}
#'     \item{\code{perfect_table}}{Data frame of all zero-zero combinations
#'       (or \code{NULL} if none found / \code{return_all_perfect = FALSE}).}
#'     \item{\code{final_blocks}}{Block table produced with \code{best_params}
#'       on all chromosomes (output of \code{\link{run_Big_LD_all_chr}}).}
#'     \item{\code{gwas_assigned}}{Input \code{gwas_df} with an added column
#'       \code{LD_block}. Entries ending in \code{*} denote forced assignments.}
#'   }
#'
#' @references
#' Weber SE, Frisch M, Snowdon RJ, Voss-Fels KP (2023). Haplotype
#' blocks for genomic prediction: a comparative evaluation in multiple
#' crop datasets. \emph{Frontiers in Plant Science} \strong{14}:1217589.
#' \doi{10.3389/fpls.2023.1217589}
#'
#' Difabachew YF et al. (2023). Genomic prediction with haplotype
#' blocks in wheat. \emph{Frontiers in Plant Science} \strong{14}:1168547.
#' \doi{10.3389/fpls.2023.1168547}
#'
#' @seealso \code{\link{run_Big_LD_all_chr}}, \code{Big_LD()}
#'
#' @export
tune_LD_params <- function(
    geno_matrix,
    snp_info,
    gwas_df,
    grid               = NULL,
    chromosomes        = NULL,
    target_bp_band     = c(5e4, 5e5),
    parallel           = FALSE,
    seed               = NULL,
    prefer_perfect     = TRUE,
    return_all_perfect = TRUE
) {
  if (!is.null(seed)) set.seed(seed)

  req1 <- c("SNP", "CHR", "POS"); req2 <- c("Marker", "CHR", "POS")
  if (!all(req1 %in% names(snp_info)))
    stop("snp_info must contain: ", paste(req1, collapse = ", "))
  if (!all(req2 %in% names(gwas_df)))
    stop("gwas_df must contain: ",  paste(req2, collapse = ", "))

  if (parallel && !requireNamespace("future.apply", quietly = TRUE))
    stop("parallel = TRUE requires the 'future.apply' package.")

  if (is.null(grid)) {
    grid <- expand.grid(
      CLQcut       = c(0.65, 0.70, 0.75, 0.80),
      min_freq     = c(0.01, 0.05),
      clstgap      = 2e6,
      leng         = 1000,
      subSegmSize  = 10000,
      split        = TRUE,
      checkLargest = FALSE,
      MAFcut       = 0.05,
      CLQmode      = "Density",
      kin_method   = "chol",
      digits       = 6,
      appendrare   = FALSE,
      stringsAsFactors = FALSE
    )
  }

  snp_info$CHR   <- as.character(snp_info$CHR)
  if (is.null(chromosomes)) chromosomes <- unique(snp_info$CHR)
  chromosomes    <- as.character(chromosomes)
  chr_index      <- stats::setNames(
    lapply(chromosomes, function(chr) which(snp_info$CHR == chr)),
    chromosomes
  )

  eval_one_combo <- function(row_idx) {
    params     <- as.list(grid[row_idx, , drop = FALSE])
    if (is.null(params$min_freq)) params$min_freq <- 0.01  # backward compat
    ld_by_chr  <- vector("list", length(chromosomes))
    names(ld_by_chr) <- chromosomes

    for (k in seq_along(chromosomes)) {
      chr  <- chromosomes[k]
      idx  <- chr_index[[k]]
      if (length(idx) < 10L) next
      geno_chr <- geno_matrix[, idx, drop = FALSE]
      info_chr <- snp_info[idx, c("SNP", "POS")]
      if (!is.numeric(info_chr$POS)) info_chr$POS <- as.numeric(info_chr$POS)
      blk <- try(.build_blocks_chr(geno_chr, info_chr, params), silent = TRUE)
      if (!inherits(blk, "try-error") && !is.null(blk) && nrow(blk) > 0L) {
        blk$CHR <- chr
        ld_by_chr[[k]] <- blk
      }
    }
    ld_blocks  <- .consolidate_blocks(ld_by_chr)
    assign_res <- .assign_gwas(ld_blocks, gwas_df)
    score      <- .score_combo(assign_res, ld_blocks, target_bp = target_bp_band)

    data.frame(
      row             = row_idx,
      CLQcut          = params$CLQcut,
      min_freq        = params$min_freq,
      clstgap         = params$clstgap,
      leng            = params$leng,
      subSegmSize     = params$subSegmSize,
      split           = params$split,
      checkLargest    = params$checkLargest,
      MAFcut          = params$MAFcut,
      CLQmode         = params$CLQmode,
      kin_method       = params$kin_method,
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
  if (!nrow(scores)) stop("No parameter combinations could be evaluated.")

  perfect <- if (isTRUE(prefer_perfect))
    dplyr::filter(scores, n_unassigned == 0L, n_forced == 0L) else NULL

  if (!is.null(perfect) && nrow(perfect) > 0L) {
    best <- dplyr::slice(dplyr::arrange(perfect, n_blocks, penalty_bp), 1L)
  } else {
    best <- dplyr::slice(
      dplyr::arrange(scores, n_unassigned, n_forced, n_blocks, penalty_bp), 1L)
  }

  best_params <- as.list(best[1L, c(
    "CLQcut", "min_freq", "clstgap", "leng", "subSegmSize", "split",
    "checkLargest", "MAFcut", "CLQmode", "kin_method", "digits", "appendrare"
  )])

  # -- Final run with best params (all chromosomes) ---------------------------
  all_chr        <- unique(as.character(snp_info$CHR))
  final_ld_by_chr <- vector("list", length(all_chr))
  names(final_ld_by_chr) <- all_chr

  for (chr in all_chr) {
    idx <- which(snp_info$CHR == chr)
    if (length(idx) < 10L) next
    geno_chr <- geno_matrix[, idx, drop = FALSE]
    info_chr <- snp_info[idx, c("SNP", "POS")]
    if (!is.numeric(info_chr$POS)) info_chr$POS <- as.numeric(info_chr$POS)
    blk <- try(.build_blocks_chr(geno_chr, info_chr, best_params), silent = TRUE)
    if (!inherits(blk, "try-error") && !is.null(blk)) {
      blk$CHR <- chr; final_ld_by_chr[[chr]] <- blk
    }
  }

  final_blocks <- .consolidate_blocks(final_ld_by_chr)
  final_assign <- .assign_gwas(final_blocks, gwas_df)

  list(
    best_params   = best_params,
    score_table   = scores,
    perfect_table = if (isTRUE(return_all_perfect) && !is.null(perfect) &&
                        nrow(perfect) > 0L) perfect else NULL,
    final_blocks  = final_blocks,
    gwas_assigned = final_assign$gwas
  )
}
