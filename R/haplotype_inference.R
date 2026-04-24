# ==============================================================================
# haplotype_inference.R
# True haplotype inference, collapsing, and harmonization layer.
#
# These functions address the gap between phase-free diploid allele strings
# (produced by extract_haplotypes() from unphased dosage data) and true
# gametic haplotype states needed for diplotype-aware analyses.
#
# 1. infer_block_haplotypes()   Structured diplotype table (hap1/hap2/diplotype
#                               code/ambiguity) from phased or unphased input.
# 2. collapse_haplotypes()      Merge rare alleles into biologically meaningful
#                               groups rather than dropping them.
# 3. harmonize_haplotypes()     Make allele labels transferable across panels,
#                               environments, or training/validation splits.
# ==============================================================================


# ==============================================================================
# 1. infer_block_haplotypes
# ==============================================================================

#' Infer Structured Block-Level Diplotypes Per Individual
#'
#' @description
#' Converts the raw haplotype strings produced by \code{\link{extract_haplotypes}}
#' into a structured per-individual, per-block diplotype table with explicit
#' \code{hap1} / \code{hap2} gamete alleles, a diplotype code, a heterozygosity
#' flag, and a phase-ambiguity flag.
#'
#' When the input is a \strong{phased list} (from \code{\link{read_phased_vcf}}
#' or a pre-phased VCF via \code{\link{read_phased_vcf}}), \code{hap1} and \code{hap2} are the
#' true gametic strings and \code{phase_ambiguous} is always \code{FALSE}.
#'
#' When the input is \strong{unphased} (diploid allele strings), heterozygous
#' individuals have two possible phase assignments. The function flags these
#' with \code{phase_ambiguous = TRUE}; \code{hap1} / \code{hap2} are set to
#' \code{NA} for those individuals unless \code{resolve_unphased = TRUE}, in
#' which case the most-frequent gametic split consistent with observed allele
#' frequencies is imputed (maximum-parsimony heuristic - not statistically
#' rigorous, use \code{\link{phase_with_beagle}} for rigorous phasing).
#'
#' @param haplotypes    Named list from \code{\link{extract_haplotypes}}.
#'   If produced from a phased list input, strings contain \code{"|"}
#'   separators (e.g. \code{"011|100"}). If unphased, strings are plain
#'   diploid allele codes (e.g. \code{"012201"}).
#' @param resolve_unphased Logical. For unphased heterozygous genotypes,
#'   impute the most parsimonious hap1/hap2 split using population allele
#'   frequencies (parsimony heuristic). Default \code{FALSE} (leave as
#'   \code{NA}).
#' @param missing_string Character. Missing data placeholder. Default \code{"."}.
#'
#' @return Data frame with one row per individual x block combination:
#' \describe{
#'   \item{\code{block_id}}{Block identifier.}
#'   \item{\code{CHR}, \code{start_bp}, \code{end_bp}}{Block coordinates.}
#'   \item{\code{id}}{Individual ID.}
#'   \item{\code{hap1}}{Gamete 1 allele string (\code{NA} if unresolvable).}
#'   \item{\code{hap2}}{Gamete 2 allele string (\code{NA} if unresolvable).}
#'   \item{\code{diplotype}}{Canonical diplotype code: the two gamete alleles
#'     sorted alphabetically and joined with \code{"/"}, e.g.
#'     \code{"010/110"}.}
#'   \item{\code{heterozygous}}{Logical; \code{TRUE} when the individual is
#'     biologically heterozygous. For phased input: \code{hap1 != hap2}.
#'     For unphased input: any dosage position equals 1, regardless of
#'     whether phase was resolved.}
#'   \item{\code{phase_ambiguous}}{Logical; \code{TRUE} for unphased
#'     heterozygous individuals when \code{resolve_unphased = FALSE}.
#'     These individuals have \code{heterozygous = TRUE} but
#'     \code{hap1 = hap2 = NA}.}
#'   \item{\code{missing}}{Logical; \code{TRUE} when the allele string
#'     contains the missing placeholder.}
#' }
#'
#' @examples
#' \donttest{
#' data(ldx_geno, ldx_snp_info, ldx_blocks, package = "LDxBlocks")
#' haps <- extract_haplotypes(ldx_geno, ldx_snp_info, ldx_blocks)
#' dip  <- infer_block_haplotypes(haps)
#' # Proportion of heterozygous diplotypes per block
#' tapply(dip$heterozygous, dip$block_id, mean, na.rm = TRUE)
#' }
#' @seealso \code{\link{extract_haplotypes}}, \code{\link{phase_with_beagle}},
#'   \code{\link{collapse_haplotypes}}
#' @export
infer_block_haplotypes <- function(
    haplotypes,
    resolve_unphased = FALSE,
    missing_string   = "."
) {
  bi <- attr(haplotypes, "block_info")
  if (is.null(bi))
    stop("haplotypes must carry a block_info attribute from extract_haplotypes().",
         call. = FALSE)

  rows <- lapply(seq_along(haplotypes), function(i) {
    bn  <- names(haplotypes)[i]
    hap <- haplotypes[[bn]]
    blk <- bi[bi$block_id == bn, , drop = FALSE]
    if (!nrow(blk)) return(NULL)

    is_phased <- isTRUE(blk$phased[1L])
    inds      <- names(hap)

    if (is_phased) {
      # Strings are "hap1_str|hap2_str"
      split    <- strsplit(hap, "|", fixed = TRUE)
      h1       <- vapply(split, `[[`, character(1L), 1L)
      h2       <- vapply(split, function(x) if (length(x) >= 2) x[[2]] else NA_character_,
                         character(1L))
      missing  <- grepl(missing_string, hap, fixed = TRUE)
      het      <- h1 != h2 & !missing
      phase_amb <- rep(FALSE, length(hap))

    } else {
      # Unphased: diploid allele string - phase is unknown for hets
      missing <- grepl(missing_string, hap, fixed = TRUE)

      # For homozygous individuals: hap1 = hap2 = the allele string
      # (since all positions are 0 or 2, each gamete is 0 or 1 = half dosage)
      # For heterozygous: at least one position is "1" in the dosage string
      is_het <- vapply(hap, function(s) {
        if (grepl(missing_string, s, fixed = TRUE)) return(FALSE)
        any(strsplit(s, "")[[1]] == "1")
      }, logical(1L))

      h1        <- rep(NA_character_, length(hap))
      h2        <- rep(NA_character_, length(hap))
      phase_amb <- rep(FALSE, length(hap))

      # Homozygous: gametes = half dosage (2->1, 0->0)
      hom_idx <- which(!is_het & !missing)
      if (length(hom_idx)) {
        for (ii in hom_idx) {
          chars      <- strsplit(hap[[ii]], "")[[1]]
          gamete_str <- paste(as.integer(chars) %/% 2L, collapse = "")
          h1[ii] <- h2[ii] <- gamete_str
        }
      }

      # Heterozygous
      het_idx <- which(is_het & !missing)
      if (length(het_idx)) {
        if (resolve_unphased) {
          # Maximum-parsimony heuristic: use population allele frequency
          # to pick the most likely gamete assignment.
          # Collect all valid strings for frequency table
          valid_haps <- hap[!is_het & !missing]
          freq_tbl   <- if (length(valid_haps) > 0) table(valid_haps) else NULL

          for (ii in het_idx) {
            chars <- strsplit(hap[[ii]], "")[[1]]
            # Positions with dosage 1 are heterozygous sites
            het_pos <- which(chars == "1")
            # All possible gamete pairs consistent with dosage
            n_het  <- length(het_pos)
            # Enumerate all 2^(n_het-1) phase configurations (cap at 10 sites)
            if (n_het > 10L) {
              # Too many combinations: leave ambiguous
              phase_amb[ii] <- TRUE
              next
            }
            configs <- as.matrix(expand.grid(rep(list(0:1), n_het)))
            best_h1  <- NA_character_; best_h2 <- NA_character_
            best_sc  <- -Inf

            for (ci in seq_len(nrow(configs))) {
              g1 <- chars; g2 <- chars
              for (pi in seq_len(n_het)) {
                g1[het_pos[pi]] <- as.character(configs[ci, pi])
                g2[het_pos[pi]] <- as.character(1L - configs[ci, pi])
              }
              # Replace dosage=2 positions with 1 in each gamete
              g1[chars == "2"] <- "1"; g2[chars == "2"] <- "1"
              g1[chars == "0"] <- "0"; g2[chars == "0"] <- "0"
              s1 <- paste(g1, collapse = "")
              s2 <- paste(g2, collapse = "")
              score <- if (!is.null(freq_tbl)) {
                (if (s1 %in% names(freq_tbl)) freq_tbl[[s1]] else 0L) +
                  (if (s2 %in% names(freq_tbl)) freq_tbl[[s2]] else 0L)
              } else 0L
              if (score > best_sc) {
                best_sc <- score; best_h1 <- s1; best_h2 <- s2
              }
            }
            h1[ii] <- best_h1; h2[ii] <- best_h2
            phase_amb[ii] <- (best_sc == 0L)  # no frequency support -> still ambiguous
          }
        } else {
          phase_amb[het_idx] <- TRUE
        }
      }

      het <- is_het
    }

    # Canonical diplotype: sort hap1/hap2 alphabetically and join with "/"
    diplotype <- ifelse(
      !is.na(h1) & !is.na(h2),
      paste(pmin(h1, h2), pmax(h1, h2), sep = "/"),
      NA_character_
    )
    # For phased input, heterozygous = resolved hap1 != hap2.
    # For unphased input, preserve biological heterozygosity from the dosage
    # string (is_het) so that phase_ambiguous individuals are still flagged
    # as heterozygous = TRUE even when hap1/hap2 are NA.
    het <- if (is_phased) {
      !is.na(h1) & !is.na(h2) & (h1 != h2)
    } else {
      is_het & !missing
    }

    data.frame(
      block_id      = bn,
      CHR           = blk$CHR[1L],
      start_bp      = blk$start_bp[1L],
      end_bp        = blk$end_bp[1L],
      id            = inds,
      hap1          = h1,
      hap2          = h2,
      diplotype     = diplotype,
      heterozygous  = het,
      phase_ambiguous = phase_amb,
      missing       = missing,
      stringsAsFactors = FALSE
    )
  })

  rows <- Filter(Negate(is.null), rows)
  if (!length(rows)) return(data.frame())
  do.call(rbind, rows)
}


# ==============================================================================
# 2. collapse_haplotypes
# ==============================================================================

#' Collapse Rare Haplotype Alleles Into Biologically Meaningful Groups
#'
#' @description
#' Rather than dropping rare haplotype alleles below a frequency threshold
#' (as \code{\link{build_haplotype_feature_matrix}} does), this function merges
#' them into the most appropriate existing allele. Three strategies are
#' supported:
#'
#' \describe{
#'   \item{\code{"rare_to_other"}}{All alleles below \code{min_freq} are
#'     pooled into a single catch-all \code{"<other>"} category. Simple and
#'     lossless for total frequency; best when rare alleles are truly
#'     heterogeneous.}
#'   \item{\code{"nearest"}}{Each rare allele is merged with the most similar
#'     common allele (minimum Hamming distance). Preserves biological
#'     similarity; best when rare alleles are likely sequencing errors or
#'     very recent recombinants of existing haplotypes.}
#'   \item{\code{"tree_based"}}{Builds a UPGMA dendrogram from pairwise
#'     Hamming distances and merges rare alleles by cutting the tree at the
#'     level that produces the fewest groups while keeping all groups above
#'     \code{min_freq}. Computationally heavier but most biologically
#'     principled.}
#' }
#'
#' @param haplotypes    Named list from \code{\link{extract_haplotypes}}.
#' @param min_freq      Numeric. Alleles at or below this frequency are
#'   considered rare. Default \code{0.05}.
#' @param collapse      Character. Collapsing strategy: \code{"rare_to_other"},
#'   \code{"nearest"}, or \code{"tree_based"}. Default \code{"nearest"}.
#' @param missing_string Character. Missing data placeholder. Default \code{"."}.
#' @param keep_labels   Logical. If \code{TRUE} (default), a \code{"label_map"}
#'   attribute is attached to the output list: a named list per block giving
#'   the original -> collapsed allele mapping.
#'
#' @return Named list of the same structure as the input \code{haplotypes},
#'   with rare allele strings replaced by their collapsed equivalents.
#'   The \code{block_info} attribute is preserved. If \code{keep_labels = TRUE},
#'   a \code{label_map} attribute is also attached (used by
#'   \code{\link{harmonize_haplotypes}}).
#'
#' @examples
#' \donttest{
#' data(ldx_geno, ldx_snp_info, ldx_blocks, package = "LDxBlocks")
#' haps      <- extract_haplotypes(ldx_geno, ldx_snp_info, ldx_blocks)
#' haps_col  <- collapse_haplotypes(haps, min_freq = 0.05, collapse = "nearest")
#' # Compare diversity before/after
#' div_before <- compute_haplotype_diversity(haps)
#' div_after  <- compute_haplotype_diversity(haps_col)
#' summary(div_after$He - div_before$He)  # small reduction expected
#' }
#' @seealso \code{\link{extract_haplotypes}},
#'   \code{\link{harmonize_haplotypes}},
#'   \code{\link{build_haplotype_feature_matrix}}
#' @export
collapse_haplotypes <- function(
    haplotypes,
    min_freq       = 0.05,
    collapse       = c("nearest", "rare_to_other", "tree_based"),
    missing_string = ".",
    keep_labels    = TRUE
) {
  collapse <- match.arg(collapse)
  bi       <- attr(haplotypes, "block_info")
  label_maps <- list()

  result <- lapply(seq_along(haplotypes), function(i) {
    bn  <- names(haplotypes)[i]
    hap <- haplotypes[[bn]]

    valid   <- hap[!grepl(missing_string, hap, fixed = TRUE)]
    if (!length(valid)) { label_maps[[bn]] <<- list(); return(hap) }

    tbl    <- table(valid)
    freq   <- as.numeric(tbl) / sum(tbl)
    names(freq) <- names(tbl)

    common <- names(freq)[freq > min_freq]
    rare   <- names(freq)[freq <= min_freq]

    if (!length(rare)) { label_maps[[bn]] <<- list(); return(hap) }
    if (!length(common)) {
      # All alleles are rare: keep as-is (no common anchor to merge into)
      label_maps[[bn]] <<- list()
      return(hap)
    }

    # Build mapping: rare allele -> replacement
    mapping <- switch(collapse,

                      rare_to_other = {
                        stats::setNames(rep("<other>", length(rare)), rare)
                      },

                      nearest = {
                        # Hamming distance of each rare allele to each common allele
                        .hamming <- function(a, b) {
                          ca <- strsplit(a, "")[[1]]; cb <- strsplit(b, "")[[1]]
                          if (length(ca) != length(cb)) return(Inf)
                          sum(ca != cb)
                        }
                        vapply(rare, function(r) {
                          dists <- vapply(common, function(c) .hamming(r, c), numeric(1L))
                          common[which.min(dists)]
                        }, character(1L))
                      },

                      tree_based = {
                        all_al  <- c(common, rare)
                        n_al    <- length(all_al)
                        ham_mat <- matrix(0L, n_al, n_al, dimnames = list(all_al, all_al))
                        for (a in seq_len(n_al - 1L)) {
                          for (b in seq(a + 1L, n_al)) {
                            ca <- strsplit(all_al[a], "")[[1]]
                            cb <- strsplit(all_al[b], "")[[1]]
                            d  <- if (length(ca) == length(cb)) sum(ca != cb) else nchar(all_al[a])
                            ham_mat[a, b] <- ham_mat[b, a] <- d
                          }
                        }
                        hc      <- stats::hclust(stats::as.dist(ham_mat), method = "average")
                        # Cut tree to find the coarsest grouping where all common alleles
                        # remain as singletons (don't merge common with common)
                        # Simple approach: cut at height = min distance between any two common alleles
                        if (length(common) > 1L) {
                          com_idx <- match(common, all_al)
                          com_dists <- ham_mat[com_idx, com_idx]
                          diag(com_dists) <- Inf
                          cut_h <- min(com_dists) - 0.5
                        } else {
                          cut_h <- max(ham_mat) * 0.5
                        }
                        groups <- stats::cutree(hc, h = max(cut_h, 0.5))
                        # For each rare allele, find its group and the most frequent common
                        # allele in that group
                        vapply(rare, function(r) {
                          g <- groups[r]
                          group_members <- names(groups)[groups == g]
                          group_common  <- intersect(group_members, common)
                          if (!length(group_common)) {
                            # No common allele in group; find nearest common overall
                            dists <- ham_mat[r, common]
                            common[which.min(dists)]
                          } else {
                            # Most frequent common allele in the group
                            group_common[which.max(freq[group_common])]
                          }
                        }, character(1L))
                      }
    )  # end switch

    label_maps[[bn]] <<- as.list(mapping)

    # Apply mapping to haplotype strings
    new_hap <- hap
    for (r in names(mapping)) {
      new_hap[new_hap == r] <- mapping[[r]]
    }
    new_hap
  })

  names(result) <- names(haplotypes)
  attr(result, "block_info") <- bi
  if (keep_labels) attr(result, "label_map") <- label_maps
  result
}


# ==============================================================================
# 3. harmonize_haplotypes
# ==============================================================================

#' Harmonize Haplotype Allele Labels Across Panels or Analysis Runs
#'
#' @description
#' Ensures that haplotype allele labels are biologically comparable across
#' different datasets, analysis runs, or training/validation splits. Without
#' harmonization, the allele string \code{"010110"} in one panel is not
#' guaranteed to correspond to the same biological haplotype in another panel
#' if block boundaries, SNP ordering, or allele encoding differ between runs.
#'
#' The function anchors allele identity to a \strong{reference dictionary}
#' built from a training/reference panel. New (target) haplotypes are then
#' matched against this dictionary:
#' \enumerate{
#'   \item \strong{Exact match}: the allele string exists verbatim in the
#'     reference dictionary -> labelled with the reference allele label.
#'   \item \strong{Nearest-Hamming match}: no exact match -> labelled with
#'     the most similar reference allele (minimum Hamming distance). If the
#'     minimum Hamming distance exceeds \code{max_hamming}, the allele is
#'     labelled \code{"<novel>"}.
#'   \item \strong{Novel}: distance > \code{max_hamming} -> \code{"<novel>"}.
#' }
#'
#' @param haplotypes_target Named list from \code{\link{extract_haplotypes}}
#'   (the panel to harmonize - validation set, new environment, etc.).
#' @param haplotypes_ref    Named list from \code{\link{extract_haplotypes}}
#'   (the reference panel - training set, base population, etc.).
#'   Must cover the same blocks as \code{haplotypes_target} (extra blocks
#'   in either panel are silently skipped).
#' @param min_freq_ref  Numeric. Only alleles above this frequency in the
#'   reference panel form the dictionary. Default \code{0.02}.
#' @param max_hamming   Integer. Maximum Hamming distance for a
#'   nearest-neighbour match; alleles beyond this distance are labelled
#'   \code{"<novel>"}. Default \code{NULL} (no limit - always assigns to
#'   nearest reference allele).
#' @param missing_string Character. Missing data placeholder. Default \code{"."}.
#'
#' @return Named list of the same structure as \code{haplotypes_target}, with
#'   allele strings replaced by their reference-anchored equivalents. The
#'   \code{block_info} attribute from \code{haplotypes_target} is preserved.
#'   A \code{harmonization_report} attribute is attached: a data frame with
#'   one row per block reporting \code{n_exact}, \code{n_nearest},
#'   \code{n_novel}, and \code{mean_hamming_dist} for matched alleles.
#'
#' @examples
#' \donttest{
#' data(ldx_geno, ldx_snp_info, ldx_blocks, package = "LDxBlocks")
#' # Split into training (70 pct) and validation (30 pct)
#' n    <- nrow(ldx_geno)
#' idx  <- sample(n)
#' ref_geno  <- ldx_geno[idx[1:round(n*0.7)], ]
#' tgt_geno  <- ldx_geno[idx[(round(n*0.7)+1):n], ]
#' haps_ref  <- extract_haplotypes(ref_geno, ldx_snp_info, ldx_blocks)
#' haps_tgt  <- extract_haplotypes(tgt_geno, ldx_snp_info, ldx_blocks)
#' haps_harm <- harmonize_haplotypes(haps_tgt, haps_ref)
#' attr(haps_harm, "harmonization_report")
#' }
#' @seealso \code{\link{extract_haplotypes}}, \code{\link{collapse_haplotypes}},
#'   \code{\link{build_haplotype_feature_matrix}}
#' @export
harmonize_haplotypes <- function(
    haplotypes_target,
    haplotypes_ref,
    min_freq_ref   = 0.02,
    max_hamming    = NULL,
    missing_string = "."
) {
  bi_tgt <- attr(haplotypes_target, "block_info")
  bi_ref <- attr(haplotypes_ref,    "block_info")

  if (is.null(bi_tgt) || is.null(bi_ref))
    stop("Both haplotype lists must carry a block_info attribute from extract_haplotypes().",
         call. = FALSE)

  # Blocks present in both panels
  common_blocks <- intersect(names(haplotypes_target), names(haplotypes_ref))
  if (!length(common_blocks))
    stop("No common block IDs between target and reference haplotypes.\n",
         "Ensure both were extracted with the same blocks table.", call. = FALSE)

  report_rows <- list()

  result <- lapply(names(haplotypes_target), function(bn) {
    hap_tgt <- haplotypes_target[[bn]]

    if (!bn %in% common_blocks) {
      # Block not in reference: return as-is
      report_rows[[bn]] <<- data.frame(block_id = bn, n_exact = NA_integer_,
                                       n_nearest = NA_integer_, n_novel = NA_integer_,
                                       mean_hamming_dist = NA_real_, stringsAsFactors = FALSE)
      return(hap_tgt)
    }

    hap_ref <- haplotypes_ref[[bn]]
    valid_ref <- hap_ref[!grepl(missing_string, hap_ref, fixed = TRUE)]
    if (!length(valid_ref)) {
      report_rows[[bn]] <<- data.frame(block_id = bn, n_exact = 0L,
                                       n_nearest = 0L, n_novel = 0L, mean_hamming_dist = NA_real_,
                                       stringsAsFactors = FALSE)
      return(hap_tgt)
    }

    # Build reference dictionary: alleles above min_freq_ref
    tbl_ref  <- table(valid_ref)
    freq_ref <- as.numeric(tbl_ref) / sum(tbl_ref)
    names(freq_ref) <- names(tbl_ref)
    dict     <- names(freq_ref)[freq_ref >= min_freq_ref]

    if (!length(dict)) dict <- names(freq_ref)  # fallback: use all

    # Hamming distance helper (vectorised over dict)
    .ham_to_dict <- function(s) {
      chars_s <- strsplit(s, "")[[1]]
      vapply(dict, function(d) {
        chars_d <- strsplit(d, "")[[1]]
        if (length(chars_s) != length(chars_d)) return(Inf)
        sum(chars_s != chars_d)
      }, numeric(1L))
    }

    n_exact <- 0L; n_nearest <- 0L; n_novel <- 0L
    hamming_dists <- numeric(0)

    new_hap <- vapply(hap_tgt, function(s) {

      if (grepl(missing_string, s, fixed = TRUE)) return(s)

      # 1. Exact match
      if (s %in% dict) {
        n_exact <<- n_exact + 1L
        return(s)
      }

      # 2. Nearest-Hamming match
      dists  <- .ham_to_dict(s)
      min_d  <- min(dists)

      if (!is.null(max_hamming) && min_d > max_hamming) {
        n_novel <<- n_novel + 1L
        return("<novel>")
      }

      nearest <- dict[which.min(dists)]
      hamming_dists <<- c(hamming_dists, min_d)
      n_nearest <<- n_nearest + 1L
      nearest

    }, character(1L))

    report_rows[[bn]] <<- data.frame(
      block_id          = bn,
      n_exact           = n_exact,
      n_nearest         = n_nearest,
      n_novel           = n_novel,
      mean_hamming_dist = if (length(hamming_dists)) mean(hamming_dists) else NA_real_,
      stringsAsFactors  = FALSE
    )

    new_hap
  })

  names(result) <- names(haplotypes_target)
  attr(result, "block_info")           <- bi_tgt
  attr(result, "harmonization_report") <- do.call(rbind, report_rows)
  result
}
