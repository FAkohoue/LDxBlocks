# ==============================================================================
# breeding_decisions.R
# Haplotype-based breeding decision tools.
#
# 1. score_favorable_haplotypes()    Score each individual's haplotype
#                                    portfolio against known allele effects.
# 2. summarize_parent_haplotypes()   Per-candidate-parent allele inventory
#                                    for haplotype stacking decisions.
# ==============================================================================


# ==============================================================================
# 1. score_favorable_haplotypes
# ==============================================================================

#' Score Individual Haplotype Portfolios Against Known Allele Effects
#'
#' @description
#' Scores each individual's genome-wide haplotype composition against a table
#' of per-allele effects (from genomic prediction or association analysis),
#' producing a per-block score and a genome-wide stacking index. This is the
#' primary tool for translating haplotype genomic prediction results into
#' actionable breeding rankings.
#'
#' \strong{Scoring rule:} For each individual \eqn{i} and each LD block
#' \eqn{b}, the block score is:
#' \deqn{s_{ib} = \sum_{j} \hat{\alpha}_j \cdot d_{ij}}
#' where \eqn{\hat{\alpha}_j} is the known effect of allele \eqn{j} and
#' \eqn{d_{ij}} is the dosage (copies) of allele \eqn{j} carried by individual
#' \eqn{i} (0/1 for unphased data; 0/1/2 for phased data). The genome-wide
#' stacking index is the sum of block scores across all scored blocks,
#' optionally normalised to [0, 1].
#'
#' @param haplotypes Named list produced by \code{\link{extract_haplotypes}}.
#'   Must carry a \code{block_info} attribute. Individual IDs are taken from
#'   \code{names(haplotypes[[1]])} (names of the first block's haplotype vector).
#'   All individuals present in \code{haplotypes} are scored, including those
#'   not in \code{allele_effects} (they receive a score of 0 for missing blocks).
#'
#' @param allele_effects Data frame specifying known per-allele additive effects.
#'   \strong{Required columns:}
#'   \itemize{
#'     \item \code{block_id} (character) — Block identifier matching
#'       \code{names(haplotypes)}.
#'     \item \code{allele} (character) — Haplotype allele string matching the
#'       strings in \code{haplotypes[[block_id]]}.
#'     \item \code{allele_effect} (numeric) — Effect size. Positive = allele
#'       increases trait value; negative = decreases it. Units are the same
#'       as the phenotype scale used to estimate effects.
#'   }
#'   Accepts the output of \code{\link{decompose_block_effects}} directly.
#'   When effects from \code{\link{test_block_haplotypes}} are used, filter
#'   to one trait first (\code{allele_effects[allele_effects$trait == "YLD", ]})
#'   and rename \code{effect} to \code{allele_effect}. Blocks or alleles in
#'   \code{haplotypes} with no matching entry in \code{allele_effects}
#'   contribute a score of 0 for those individuals.
#'
#' @param min_freq Numeric in (0, 1). Minimum allele frequency in the full
#'   panel. Alleles with population frequency below this threshold are excluded
#'   from scoring even if they appear in \code{allele_effects}. Default
#'   \code{0.02}. This prevents rare private alleles from dominating scores
#'   based on unreliable frequency estimates.
#'
#' @param missing_string Character. The string used in haplotype vectors to
#'   denote missing genotype data (e.g. for individuals with insufficient
#'   marker coverage at a block). Individuals with this string at a block
#'   receive \code{NA} for that block's score (excluded from mean but still
#'   contribute to \code{n_blocks_scored} = 0 for that block). Default
#'   \code{"."}.
#'
#' @param normalize Logical. If \code{TRUE} (default), the genome-wide
#'   stacking index is linearly scaled to [0, 1] across all individuals:
#'   \deqn{\text{index}_i = (S_i - S_{\min}) / (S_{\max} - S_{\min})}
#'   where \eqn{S_i} is the raw sum of block scores. This makes the index
#'   interpretable as a percentile rank within the panel and comparable
#'   across panels with different numbers of scored blocks. Set \code{FALSE}
#'   to return raw summed effects (in phenotype units), which is preferable
#'   when comparing candidate sets across different evaluation sets.
#'
#' @return Data frame with one row per individual, sorted ascending by
#'   \code{rank} (best candidates first). Contains the following columns:
#' \describe{
#'   \item{\code{id}}{Character. Individual identifier, taken from
#'     \code{names(haplotypes[[1]])}.}
#'   \item{\code{stacking_index}}{Numeric. Genome-wide sum of per-block
#'     haplotype scores. When \code{normalize = TRUE}, scaled to [0, 1] where
#'     1.0 = the individual with the most favourable genome-wide haplotype
#'     combination in the panel and 0.0 = the least favourable. When
#'     \code{normalize = FALSE}, units are phenotype units × allele dosage.}
#'   \item{\code{n_blocks_scored}}{Integer. Number of LD blocks for which at
#'     least one allele effect was available in \code{allele_effects} and
#'     at least one allele passed \code{min_freq}. Maximum possible value
#'     equals the number of unique \code{block_id} values in
#'     \code{allele_effects} that overlap with \code{names(haplotypes)}.}
#'   \item{\code{mean_block_score}}{Numeric. Raw (unnormalised) mean per-block
#'     score: \code{stacking_index_raw / n_blocks_scored}. Useful for comparing
#'     candidates across panels with different numbers of scored blocks, since
#'     it is independent of \code{n_blocks_scored}.}
#'   \item{\code{rank}}{Integer. Rank based on \code{stacking_index}, with 1
#'     indicating the individual with the highest (most favourable) index.
#'     Ties are broken by \code{"min"} (tied individuals share the lower rank
#'     number).}
#'   \item{\code{score_<block_id>}}{Numeric. One column per scored LD block,
#'     named \code{score_} followed by the block identifier (e.g.
#'     \code{score_block_1_1000_103000}). Contains the raw per-block score for
#'     each individual (effect × dosage sum across alleles). Zero when the
#'     individual carries no alleles with known effects at that block. Used to
#'     identify which genomic regions drive an individual's stacking index.}
#' }
#'
#' @examples
#' \donttest{
#' data(ldx_geno, ldx_snp_info, ldx_blocks, ldx_blues, package = "LDxBlocks")
#' haps    <- extract_haplotypes(ldx_geno, ldx_snp_info, ldx_blocks)
#' hap_mat <- build_haplotype_feature_matrix(haps)
#' G       <- compute_haplotype_grm(hap_mat)
#' pred    <- run_haplotype_prediction(ldx_geno, ldx_snp_info, ldx_blocks,
#'                                      blues    = setNames(ldx_blues$YLD,
#'                                                          ldx_blues$id),
#'                                      verbose  = FALSE)
#' ae <- decompose_block_effects(haps, ldx_snp_info, ldx_blocks,
#'                                snp_effects = pred$snp_effects)
#' scores <- score_favorable_haplotypes(haps, allele_effects = ae)
#' head(scores[order(scores$rank), c("id","stacking_index","rank")])
#' }
#' @seealso \code{\link{decompose_block_effects}},
#'   \code{\link{summarize_parent_haplotypes}},
#'   \code{\link{run_haplotype_prediction}}
#' @export
score_favorable_haplotypes <- function(
    haplotypes,
    allele_effects,
    min_freq       = 0.02,
    missing_string = ".",
    normalize      = TRUE
) {
  req <- c("block_id", "allele", "allele_effect")
  miss <- setdiff(req, names(allele_effects))
  if (length(miss))
    stop("allele_effects missing columns: ", paste(miss, collapse=", "),
         call. = FALSE)

  bi <- attr(haplotypes, "block_info")

  # Get individual IDs from the first block present
  all_ids <- NULL
  for (bn in names(haplotypes)) {
    all_ids <- names(haplotypes[[bn]])
    if (!is.null(all_ids)) break
  }
  if (is.null(all_ids))
    stop("Could not determine individual IDs from haplotypes.", call. = FALSE)

  # Build per-block per-individual score matrix
  scored_blocks <- intersect(names(haplotypes),
                             unique(allele_effects$block_id))
  if (!length(scored_blocks)) {
    warning("No block IDs overlap between haplotypes and allele_effects.",
            call. = FALSE)
    return(data.frame(id = all_ids, stacking_index = 0,
                      n_blocks_scored = 0L, mean_block_score = 0,
                      rank = seq_along(all_ids)))
  }

  score_mat <- matrix(NA_real_, nrow = length(all_ids),
                      ncol = length(scored_blocks),
                      dimnames = list(all_ids, scored_blocks))

  for (bn in scored_blocks) {
    hap <- haplotypes[[bn]]
    ae_bn <- allele_effects[allele_effects$block_id == bn, ]
    if (!nrow(ae_bn)) next

    # Allele frequency filter
    valid <- hap[!grepl(missing_string, hap, fixed = TRUE)]
    if (!length(valid)) next
    freq_tbl <- table(valid) / length(valid)
    ae_bn <- ae_bn[ae_bn$allele %in%
                     names(freq_tbl)[freq_tbl >= min_freq], ]
    if (!nrow(ae_bn)) next

    # Build effect lookup: allele_string -> effect
    eff_map <- stats::setNames(ae_bn$allele_effect, ae_bn$allele)

    for (ind in all_ids) {
      s <- hap[ind]
      if (is.na(s) || grepl(missing_string, s, fixed = TRUE)) next

      # Check if this is phased (contains |)
      if (grepl("|", s, fixed = TRUE)) {
        parts <- strsplit(s, "|", fixed = TRUE)[[1]]
        score_i <- sum(vapply(parts, function(p) {
          if (p %in% names(eff_map)) eff_map[[p]] else 0
        }, numeric(1L)))
      } else {
        # Unphased: check if the string matches any known allele
        score_i <- if (s %in% names(eff_map)) eff_map[[s]] else 0
      }
      score_mat[ind, bn] <- score_i
    }
  }

  # Aggregate
  block_sums <- rowSums(score_mat, na.rm = TRUE)
  n_scored   <- rowSums(!is.na(score_mat))
  mean_score <- ifelse(n_scored > 0, block_sums / n_scored, 0)

  if (normalize && stats::var(block_sums) > 1e-10) {
    min_s <- min(block_sums); max_s <- max(block_sums)
    stack_idx <- (block_sums - min_s) / (max_s - min_s)
  } else {
    stack_idx <- block_sums
  }

  out <- data.frame(
    id              = all_ids,
    stacking_index  = round(stack_idx, 6),
    n_blocks_scored = n_scored,
    mean_block_score = round(mean_score, 6),
    stringsAsFactors = FALSE
  )
  out$rank <- rank(-out$stacking_index, ties.method = "first")

  # Append per-block columns
  score_df <- as.data.frame(round(score_mat, 6))
  names(score_df) <- paste0("score_", names(score_df))
  out <- cbind(out, score_df)

  out[order(out$rank), ]
}


# ==============================================================================
# 2. summarize_parent_haplotypes
# ==============================================================================

#' Summarise Haplotype Allele Inventory Per Candidate Parent
#'
#' @description
#' Produces a tidy long-format allele inventory for each candidate parent
#' individual, reporting which haplotype alleles they carry at each LD block
#' and in what dosage. This is the primary decision support table for
#' haplotype stacking crosses: it reveals which parents carry complementary
#' rare alleles at important blocks, which candidates are genomically
#' redundant (same alleles at all blocks), and which blocks should be
#' targeted in crossing schemes.
#'
#' The output includes one row per (individual × block × allele) combination,
#' including rows where the individual carries 0 copies of an allele. This
#' ensures that all candidates can be compared on the same rows for any
#' given block.
#'
#' @param haplotypes Named list produced by \code{\link{extract_haplotypes}}.
#'   Must carry a \code{block_info} attribute. All blocks in the list are
#'   included in the inventory unless filtered by \code{candidate_ids} or
#'   \code{min_freq}.
#'
#' @param candidate_ids Character vector of individual IDs to include in the
#'   inventory. Must match names in the haplotype vectors (i.e.
#'   \code{names(haplotypes[[1]])}). \code{NULL} (default) includes all
#'   individuals present in \code{haplotypes}. Supply the \code{id} column
#'   from \code{\link{score_favorable_haplotypes}} (e.g. top 20 by rank) to
#'   focus the inventory on selection candidates.
#'
#' @param allele_effects Optional data frame of per-allele effects to join to
#'   the inventory. When supplied, the \code{allele_effect} column is populated
#'   by matching on \code{block_id} and \code{allele}. \code{NULL} (default)
#'   leaves \code{allele_effect} as \code{NA} throughout.
#'   \strong{Required columns when not \code{NULL}:}
#'   \itemize{
#'     \item \code{block_id} (character) — matching \code{names(haplotypes)}.
#'     \item \code{allele} (character) — matching haplotype allele strings.
#'     \item \code{allele_effect} (numeric) — effect value per allele.
#'   }
#'   Accepts the output of \code{\link{decompose_block_effects}} directly.
#'   For \code{\link{test_block_haplotypes}} output, filter to one trait
#'   and rename \code{effect} to \code{allele_effect} first.
#'
#' @param blocks LD block table from \code{\link{run_Big_LD_all_chr}}.
#'   Used only as an alternative source of block coordinate metadata if the
#'   \code{block_info} attribute of \code{haplotypes} is missing or incomplete.
#'   Default \code{NULL} (uses \code{block_info} attribute).
#'
#' @param min_freq Numeric in (0, 1). Minimum allele frequency in the full
#'   panel (all individuals in \code{haplotypes}, not just \code{candidate_ids}).
#'   Alleles below this threshold are excluded from the inventory entirely.
#'   Default \code{0.02}. This prevents the inventory from being dominated by
#'   private alleles observed in only one or two individuals.
#'
#' @param missing_string Character. String used in haplotype vectors to
#'   indicate missing genotype data. Individuals with this string at a block
#'   are skipped for that block (no rows emitted). Default \code{"."}.
#'
#' @return Data frame in long format. One row per individual × block × allele
#'   combination, including rows where \code{dosage = 0} (allele absent).
#'   Sorted ascending by \code{id}, \code{CHR}, \code{start_bp}, then
#'   descending by \code{dosage} (so carried alleles appear before absent ones
#'   for each individual–block combination). Contains 10 columns:
#' \describe{
#'   \item{\code{id}}{Character. Individual identifier.}
#'   \item{\code{block_id}}{Character. LD block identifier matching
#'     \code{names(haplotypes)}.}
#'   \item{\code{CHR}}{Character. Chromosome label.}
#'   \item{\code{start_bp}}{Integer. Block start coordinate (base pairs).}
#'   \item{\code{end_bp}}{Integer. Block end coordinate (base pairs).}
#'   \item{\code{allele}}{Character. Haplotype allele string (dosage-coded
#'     concatenation of SNP dosages within the block, e.g. \code{"012102"}).}
#'   \item{\code{dosage}}{Integer. Number of copies of this allele carried by
#'     the individual: 0 = absent; 1 = one copy (heterozygous for phased data,
#'     or present for unphased data); 2 = two copies (homozygous, phased data
#'     only). For unphased input (no \code{"|"} in strings), dosage is always
#'     0 or 1 — the two chromosomes cannot be distinguished, so an individual
#'     matching this allele string exactly receives dosage 1 regardless of
#'     whether they are truly homozygous or heterozygous.}
#'   \item{\code{allele_freq}}{Numeric in [0, 1]. Population frequency of this
#'     allele in the full panel (all individuals, not just candidates). Computed
#'     as \code{table(valid_haplotypes) / length(valid_haplotypes)} before
#'     any \code{candidate_ids} filtering is applied.}
#'   \item{\code{allele_effect}}{Numeric or \code{NA}. Effect value from
#'     \code{allele_effects} if supplied and the allele–block combination
#'     is matched; \code{NA} otherwise.}
#'   \item{\code{is_rare}}{Logical. \code{TRUE} when
#'     \code{allele_freq < 0.10}. Convenience flag for filtering the inventory
#'     to rare alleles that may be private to specific parents and therefore
#'     valuable for introgression.}
#' }
#'
#' @examples
#' \donttest{
#' data(ldx_geno, ldx_snp_info, ldx_blocks, ldx_blues, package = "LDxBlocks")
#' haps <- extract_haplotypes(ldx_geno, ldx_snp_info, ldx_blocks)
#' pred <- run_haplotype_prediction(ldx_geno, ldx_snp_info, ldx_blocks,
#'                                   blues   = setNames(ldx_blues$YLD,
#'                                                      ldx_blues$id),
#'                                   verbose = FALSE)
#' ae <- decompose_block_effects(haps, ldx_snp_info, ldx_blocks,
#'                                snp_effects = pred$snp_effects)
#' # Inventory for top 5 candidates
#' top5 <- rownames(ldx_geno)[1:5]
#' inv  <- summarize_parent_haplotypes(haps, candidate_ids = top5,
#'                                      allele_effects = ae)
#' # Show blocks where candidates carry different alleles
#' inv[inv$dosage > 0, c("id","block_id","allele","dosage","allele_effect")]
#' }
#' @seealso \code{\link{score_favorable_haplotypes}},
#'   \code{\link{decompose_block_effects}},
#'   \code{\link{compare_haplotype_populations}}
#' @export
summarize_parent_haplotypes <- function(
    haplotypes,
    candidate_ids  = NULL,
    allele_effects = NULL,
    blocks         = NULL,
    min_freq       = 0.02,
    missing_string = "."
) {
  bi <- attr(haplotypes, "block_info")

  # Validate allele_effects if supplied
  if (!is.null(allele_effects)) {
    req <- c("block_id", "allele", "allele_effect")
    miss <- setdiff(req, names(allele_effects))
    if (length(miss))
      stop("allele_effects missing: ", paste(miss, collapse=", "), call.=FALSE)
  }

  rows <- list()

  for (bn in names(haplotypes)) {
    hap <- haplotypes[[bn]]
    blk_meta <- if (!is.null(bi)) bi[bi$block_id == bn, , drop=FALSE]
    else NULL

    # Determine individuals to summarise
    inds <- if (!is.null(candidate_ids)) {
      intersect(names(hap), candidate_ids)
    } else {
      names(hap)
    }
    if (!length(inds)) next

    # Population allele frequencies (full panel)
    valid_all <- hap[!grepl(missing_string, hap, fixed=TRUE)]
    if (!length(valid_all)) next
    freq_tbl <- table(valid_all) / length(valid_all)
    common_alleles <- names(freq_tbl)[freq_tbl >= min_freq]
    if (!length(common_alleles)) next

    # Per-individual allele dosage
    for (ind in inds) {
      s <- hap[ind]
      if (is.na(s) || grepl(missing_string, s, fixed=TRUE)) next

      is_phased <- grepl("|", s, fixed=TRUE)

      if (is_phased) {
        parts <- strsplit(s, "|", fixed=TRUE)[[1]]
        # Count copies of each common allele
        for (al in common_alleles) {
          dosage <- sum(parts == al)
          ae_val <- if (!is.null(allele_effects)) {
            ae_row <- allele_effects[allele_effects$block_id == bn &
                                       allele_effects$allele == al, ]
            if (nrow(ae_row) > 0) ae_row$allele_effect[1] else NA_real_
          } else NA_real_

          rows[[length(rows)+1L]] <- data.frame(
            id = ind,
            block_id = bn,
            CHR = if (!is.null(blk_meta) && nrow(blk_meta)>0) blk_meta$CHR[1] else NA_character_,
            start_bp = if (!is.null(blk_meta) && nrow(blk_meta)>0) blk_meta$start_bp[1] else NA_integer_,
            end_bp   = if (!is.null(blk_meta) && nrow(blk_meta)>0) blk_meta$end_bp[1] else NA_integer_,
            allele   = al,
            dosage   = dosage,
            allele_freq  = round(as.numeric(freq_tbl[al]), 4),
            allele_effect = ae_val,
            is_rare  = as.numeric(freq_tbl[al]) < 0.10,
            stringsAsFactors = FALSE
          )
        }
      } else {
        # Unphased: dosage = 1 if string matches, 0 otherwise
        for (al in common_alleles) {
          dosage <- if (s == al) 1L else 0L
          ae_val <- if (!is.null(allele_effects)) {
            ae_row <- allele_effects[allele_effects$block_id == bn &
                                       allele_effects$allele == al, ]
            if (nrow(ae_row) > 0) ae_row$allele_effect[1] else NA_real_
          } else NA_real_

          rows[[length(rows)+1L]] <- data.frame(
            id = ind,
            block_id = bn,
            CHR = if (!is.null(blk_meta) && nrow(blk_meta)>0) blk_meta$CHR[1] else NA_character_,
            start_bp = if (!is.null(blk_meta) && nrow(blk_meta)>0) blk_meta$start_bp[1] else NA_integer_,
            end_bp   = if (!is.null(blk_meta) && nrow(blk_meta)>0) blk_meta$end_bp[1] else NA_integer_,
            allele   = al,
            dosage   = dosage,
            allele_freq  = round(as.numeric(freq_tbl[al]), 4),
            allele_effect = ae_val,
            is_rare  = as.numeric(freq_tbl[al]) < 0.10,
            stringsAsFactors = FALSE
          )
        }
      }
    }
  }

  if (!length(rows)) return(data.frame())
  out <- do.call(rbind, rows)
  out[order(out$id, out$CHR, out$start_bp, -out$dosage), ]
}
