# ==============================================================================
# ld_decay.R
# LD decay analysis: compute r2 vs physical distance, fit decay curves,
# determine critical r2 thresholds, and export chromosome-specific decay
# distances for use in define_qtl_regions() candidate gene window definition.
#
# MEMORY MODEL
# ============
# The naive approach (read full chromosome, compute all pairs) loads
# n_snps_chr x n_ind x 8 bytes per chromosome -- up to 2 GB for a 500k-SNP
# chromosome with 500 individuals. This is impractical for WGS panels.
#
# This implementation uses a position-first, index-sample, minimal-read strategy:
#
#   Random sampling:
#     1. Load SNP positions only (snp_info -- no genotypes needed).
#     2. Sample pair (i, j) indices such that |pos[i] - pos[j]| <= max_dist.
#     3. Collect the UNIQUE SNP indices involved in sampled pairs.
#     4. read_chunk(be, unique_indices) -- loads only those columns.
#     5. Compute r2 for each pair from the compact submatrix.
#     Peak RAM: n_unique x n_ind x 8 bytes (typically << full chromosome).
#     For n_pairs=50k on a 500k-SNP chr: ~10k unique SNPs -> ~40 MB.
#
#   Sliding window:
#     Reads exactly window_snps columns at a time via read_chunk().
#     Peak RAM: window_snps x n_ind x 8 bytes per step (~200 KB at w=50).
#
#   Parametric (cross-chromosome):
#     Samples a small fixed number of SNP indices per chromosome and calls
#     read_chunk() only for those indices. Never loads a full chromosome.
#
# References:
#   Hill WG, Weir BS (1988). Variances and covariances of squared linkage
#     disequilibria in finite populations. Theoretical Population Biology
#     33:54-78.
#   Remington DL et al. (2001). Structure of linkage disequilibrium and
#     phenotypic associations in the maize genome. PNAS 98:11479-11484.
# ==============================================================================

#' Compute LD Decay and Chromosome-Specific Decay Distances
#'
#' @description
#' Calculates the decay of linkage disequilibrium (r\eqn{^2} or rV\eqn{^2})
#' with physical distance for each chromosome, fits optional decay models
#' (Hill-Weir nonlinear or LOESS), and determines the distance at which LD
#' drops below a user-specified or data-driven critical threshold. The output
#' can be passed directly to \code{\link{define_qtl_regions}()} to define
#' chromosome-specific candidate gene windows around GWAS hits.
#'
#' @details
#' \strong{C++ acceleration:}
#' All pairwise r\eqn{^2} computations use the package's compiled C++
#' kernels (Armadillo BLAS via RcppArmadillo):
#' \itemize{
#'   \item{Random sampling: \code{compute_r2_sparse_cpp()} computes all
#'     pairs within \code{max_dist} in a single C++ call on the compact
#'     submatrix, replacing an R pair-by-pair loop. \code{n_threads}
#'     is forwarded to enable OpenMP parallelism.}
#'   \item{Sliding window: column preparation is done once per window
#'     via \code{apply()}; within-window \code{stats::cor()} is used
#'     (windows are small, typically 50 SNPs).}
#'   \item{Parametric threshold: \code{col_r2_cpp()} computes
#'     cross-chromosome r\eqn{^2} in Armadillo BLAS.}
#' }
#'
#' \strong{Memory-efficient backend access:}
#' For GDS and bigmemory backends the function never loads a full chromosome.
#' Instead it samples pair indices from SNP positions (no genotype I/O),
#' collects only the unique SNP indices involved, and calls
#' \code{read_chunk()} once for that compact set. For 50,000 random pairs
#' on a 500,000-SNP chromosome with 500 individuals, peak RAM is
#' approximately 40 MB rather than 2 GB.
#'
#' \strong{Important: \code{method = "rV2"} is not memory-safe at WGS scale.}
#' The kinship whitening matrix requires a full-genome read of all SNPs once
#' before per-chromosome processing begins. For panels with millions of markers
#' this can be the dominant memory step. Use \code{method = "r2"} for LD
#' decay estimation; reserve \code{method = "rV2"} for block detection.
#'
#' \strong{Sampling strategies:}
#' \describe{
#'   \item{\code{"sliding_window"}}{Computes r\eqn{^2} within a moving
#'     window of \code{window_snps} SNPs, reading only those columns per
#'     step. Uses an \strong{adaptive stride}: the window is placed at
#'     approximately \code{n_pairs / (window_snps*(window_snps-1)/2)}
#'     evenly-spaced positions across the chromosome, so total pairs
#'     stays near \code{n_pairs} regardless of marker density. For
#'     \code{window_snps=50} on a 500k-SNP chromosome with
#'     \code{n_pairs=50000}: stride=12,195, only 41 \code{read_chunk}
#'     calls instead of 10,000. A \strong{density guard} automatically
#'     switches chromosomes with > 200,000 SNPs to \code{"random"}
#'     sampling (consecutive SNPs at WGS density are always in near-
#'     perfect LD, making the sliding-window curve uninformative).
#'     Matches the TASSEL approach (Bradbury et al. 2007); recommended
#'     for chip-density panels (< 200k SNPs per chromosome).}
#'   \item{\code{"random"}}{Randomly samples up to \code{n_pairs} SNP pairs
#'     per chromosome within \code{max_dist} bp, then loads only the unique
#'     SNPs involved. Bounded RAM regardless of panel density; recommended
#'     for WGS panels (> 500k SNPs per chromosome).}
#'   \item{\code{"both"}}{Sliding window for the decay curve shape;
#'     random sampling for the parametric threshold estimation.}
#' }
#'
#' \strong{Critical r\eqn{^2} threshold:}
#' \describe{
#'   \item{Fixed numeric (e.g. \code{0.1})}{Standard value used in most GWAS
#'     papers. The distance where the decay curve crosses this value defines
#'     the LD block extent for candidate gene searches.}
#'   \item{\code{"parametric"}}{95th percentile of r\eqn{^2} between markers
#'     on \emph{different} chromosomes (unlinked markers). In a structured
#'     population, unlinked r\eqn{^2} > 0 due to kinship. This threshold
#'     captures the background LD attributable to structure, above which LD
#'     is genuine. A high parametric threshold (> 0.05) is a strong indicator
#'     that \code{method = "rV2"} should be used for block detection.}
#'   \item{\code{"both"}}{Returns both; uses the parametric as the active
#'     threshold for decay distance estimation.}
#' }
#'
#' \strong{Decay model:}
#' The Hill-Weir (1988) expectation relates expected r\eqn{^2} to distance
#' \eqn{d} (in Mb) via \eqn{C = 4 N_e r d}. LDxBlocks fits this per chromosome with C as the
#' free parameter via nonlinear least squares. The LOESS smooth is a
#' non-parametric alternative robust to non-standard decay shapes.
#'
#' \strong{Connection to \code{define_qtl_regions()} and GWAS integration:}
#' The output of \code{compute_ld_decay()} can be passed directly to
#' \code{\link{define_qtl_regions}(ld_decay = decay)} to replace
#' fixed block boundaries with chromosome-specific LD-aware windows.
#' \describe{
#'   \item{Without \code{ld_decay}}{\code{define_qtl_regions()} asks
#'     \emph{which blocks CONTAIN a significant GWAS SNP?} The search
#'     window equals the block boundary. GWAS hits in gaps between blocks
#'     or near block edges are missed.}
#'   \item{With \code{ld_decay}}{Asks \emph{which blocks are IN LD with
#'     a significant GWAS SNP?} The search window is extended by the
#'     chromosome-specific decay distance on both sides. Adds columns
#'     \code{candidate_region_start}, \code{candidate_region_end},
#'     and \code{candidate_region_size_kb} ready for BioMart/Ensembl
#'     Plants candidate gene queries.}
#' }
#' This directly feeds \code{\link{integrate_gwas_haplotypes}()}: the
#' GWAS biological-evidence layer (layer 1 of the three-evidence score)
#' awards a point to blocks that \emph{contain} a GWAS hit. With
#' \code{ld_decay}, blocks \emph{near} a hit (within the decay distance)
#' also receive the evidence point, making the three-evidence scoring
#' in \code{rank_haplotype_blocks()} more sensitive and biologically
#' accurate. The recommended workflow is therefore:
#' \enumerate{
#'   \item \code{decay <- compute_ld_decay(be, ...)}
#'   \item \code{qtl <- define_qtl_regions(..., ld_decay = decay)}
#'   \item \code{ranked <- integrate_gwas_haplotypes(qtl, pred, ...)}
#' }
#'
#' @param geno \code{LDxBlocks_backend} (from \code{\link{read_geno}} or
#'   \code{\link{read_geno_bigmemory}}), a numeric dosage matrix
#'   (individuals x SNPs), or a file path accepted by \code{\link{read_geno}}.
#' @param snp_info Data frame with columns \code{SNP}, \code{CHR}, \code{POS}.
#'   Auto-extracted from the backend when \code{geno} is a backend object.
#' @param method LD metric: \code{"r2"} (default) or \code{"rV2"}
#'   (kinship-adjusted; requires AGHmatrix and ASRgenomics).
#' @param kin_method Whitening method for \code{rV2}: \code{"chol"} (default)
#'   or \code{"eigen"}.
#' @param sampling Pair-sampling strategy: \code{"random"} (default, WGS-safe),
#'   \code{"sliding_window"} (TASSEL approach), or \code{"both"}.
#' @param window_snps Integer. Window size in SNPs for \code{"sliding_window"}.
#'   Default \code{50L}.
#' @param n_pairs Integer. Maximum SNP pairs per chromosome for
#'   \code{"random"} sampling. Default \code{50000L}.
#' @param max_dist Integer. Maximum physical distance (bp) between a pair.
#'   Default \code{5000000L} (5 Mb).
#' @param r2_threshold Threshold for decay distance:
#'   a numeric value (e.g. \code{0.1}), \code{"parametric"} (95th percentile
#'   of unlinked r\eqn{^2}), \code{"both"}, or \code{NULL}.
#'   Default \code{"both"}.
#' @param n_unlinked Integer. Cross-chromosome pairs for parametric threshold.
#'   Default \code{100000L}.
#' @param pctile Numeric (0-100). Percentile of unlinked r\eqn{^2} distribution
#'   for the parametric threshold. Default \code{95}.
#' @param fit_model Decay model: \code{"loess"} (default), \code{"nonlinear"}
#'   (Hill-Weir), \code{"both"}, or \code{"none"}.
#' @param n_bins Integer. Distance bins for the smoothed decay curve.
#'   Default \code{100L}.
#' @param by_chr Logical. Return per-chromosome decay distances.
#'   Default \code{TRUE}.
#' @param chr Character vector of chromosomes to process. \code{NULL} = all.
#' @param seed Integer. RNG seed. Default \code{42L}.
#' @param n_threads Integer. OpenMP threads for C++ r\eqn{^2} kernel.
#'   Default \code{1L}.
#' @param verbose Logical. Print timestamped progress. Default \code{TRUE}.
#'
#' @return A named list of class \code{LDxBlocks_decay}:
#' \describe{
#'   \item{\code{pairs}}{Data frame: \code{CHR}, \code{dist_bp}, \code{r2}.}
#'   \item{\code{decay_curve}}{Binned decay curve per chromosome:
#'     \code{CHR}, \code{dist_bp}, \code{r2_mean}, \code{r2_median},
#'     \code{n_pairs}, and optionally \code{r2_loess}, \code{r2_hw}.}
#'   \item{\code{decay_dist}}{Per-chromosome decay distance:
#'     \code{CHR}, \code{decay_dist_bp}, \code{decay_dist_kb},
#'     \code{threshold_used}, \code{r2_col_used}, \code{censored}
#'     (\code{TRUE} when the smoothed curve never drops below the threshold
#'     within \code{max_dist} -- the returned distance is a lower bound,
#'     not an estimate; extend \code{max_dist} or lower the threshold).}
#'   \item{\code{decay_dist_genome}}{Genome-wide median decay distance (bp).}
#'   \item{\code{critical_r2}}{Active threshold value.}
#'   \item{\code{critical_r2_fixed}}{Always set to 0.1 as a reference
#'     value (standard threshold used in most GWAS papers). This field is
#'     populated regardless of which \code{r2_threshold} was requested.
#'     It is the \emph{active} threshold only when \code{r2_threshold = 0.1}
#'     or \code{r2_threshold = "both"} with no parametric result available.}
#'   \item{\code{critical_r2_param}}{Parametric threshold (\code{NULL} if
#'     not computed).}
#'   \item{\code{unlinked_r2}}{Cross-chromosome r\eqn{^2} values used for
#'     the parametric threshold (\code{NULL} if not computed).}
#'   \item{\code{model_params}}{Hill-Weir C parameter per chromosome
#'     (\code{NULL} if not fitted).}
#'   \item{\code{n_pairs_used}}{Named integer: pairs computed per chromosome.}
#'   \item{\code{method}}{LD metric used.}
#'   \item{\code{sampling}}{Sampling strategy used.}
#'   \item{\code{call}}{Matched function call.}
#' }
#'
#' @examples
#' \donttest{
#' data(ldx_geno,     package = "LDxBlocks")
#' data(ldx_snp_info, package = "LDxBlocks")
#'
#' # Random sampling -- memory-safe for any panel density
#' decay <- compute_ld_decay(
#'   geno         = ldx_geno,
#'   snp_info     = ldx_snp_info,
#'   sampling     = "random",
#'   r2_threshold = "both",
#'   fit_model    = "loess",
#'   n_pairs      = 5000L,
#'   verbose      = FALSE
#' )
#' decay$critical_r2_param    # background kinship-induced LD
#' decay$critical_r2_fixed    # standard 0.1 threshold
#' decay$decay_dist           # per-chromosome decay distances
#'
#' # Sliding window (TASSEL approach) + Hill-Weir model
#' decay_hw <- compute_ld_decay(
#'   geno         = ldx_geno,
#'   snp_info     = ldx_snp_info,
#'   sampling     = "sliding_window",
#'   window_snps  = 50L,
#'   r2_threshold = 0.1,
#'   fit_model    = "nonlinear",
#'   verbose      = FALSE
#' )
#'
#' # WGS backend -- only sampled columns are loaded per chromosome
#' \dontrun{
#' be <- read_geno("my_wgs.vcf.gz")
#' decay_wgs <- compute_ld_decay(
#'   geno         = be,
#'   sampling     = "random",
#'   n_pairs      = 50000L,
#'   max_dist     = 5000000L,
#'   r2_threshold = "both",
#'   fit_model    = "loess",
#'   n_threads    = 8L
#' )
#' close_backend(be)
#' }
#'
#' # Pass to define_qtl_regions for chromosome-specific candidate windows
#' data(ldx_blocks, package = "LDxBlocks")
#' data(ldx_gwas,   package = "LDxBlocks")
#' qtl <- define_qtl_regions(ldx_gwas, ldx_blocks, ldx_snp_info,
#'                           ld_decay    = decay,
#'                           p_threshold = NULL,
#'                           trait_col   = "trait")
#' qtl[, c("block_id", "lead_marker", "candidate_region_start",
#'          "candidate_region_end", "candidate_region_size_kb")]
#' }
#'
#' @seealso \code{\link{define_qtl_regions}}, \code{\link{plot_ld_decay}},
#'   \code{\link{run_Big_LD_all_chr}}
#' @export
compute_ld_decay <- function(
    geno,
    snp_info        = NULL,
    method          = c("r2", "rV2"),
    kin_method      = c("chol", "eigen"),
    sampling        = c("random", "sliding_window", "both"),
    window_snps     = 50L,
    n_pairs         = 50000L,
    max_dist        = 5000000L,
    r2_threshold    = "both",
    n_unlinked      = 100000L,
    pctile          = 95,
    fit_model       = c("loess", "nonlinear", "both", "none"),
    n_bins          = 100L,
    by_chr          = TRUE,
    chr             = NULL,
    seed            = 42L,
    n_threads       = 1L,
    verbose         = TRUE
) {
  cl         <- match.call()
  method     <- match.arg(method)
  kin_method <- match.arg(kin_method)
  sampling   <- match.arg(sampling)
  fit_model  <- match.arg(fit_model)

  # -- Validate r2_threshold ------------------------------------------------
  use_param <- FALSE; use_fixed <- FALSE; fixed_val <- 0.1
  if (is.null(r2_threshold)) {
    # no threshold -- raw pairs only
  } else if (is.numeric(r2_threshold)) {
    if (length(r2_threshold) != 1L || is.na(r2_threshold) ||
        r2_threshold < 0 || r2_threshold > 1)
      stop("r2_threshold: when numeric, must be a single value in [0, 1].",
           call. = FALSE)
    use_fixed <- TRUE; fixed_val <- r2_threshold
  } else if (identical(r2_threshold, "parametric")) {
    use_param <- TRUE
  } else if (identical(r2_threshold, "both")) {
    use_param <- TRUE; use_fixed <- TRUE
  } else {
    stop("r2_threshold must be a numeric value in [0,1], ",
         "'parametric', 'both', or NULL.", call. = FALSE)
  }

  .log <- function(...) if (verbose)
    message(sprintf("[%s] [ld_decay] %s", format(Sys.time(), "%H:%M:%S"),
                    paste0(...)))

  # -- Open backend ---------------------------------------------------------
  own_backend <- FALSE
  if (inherits(geno, "LDxBlocks_backend")) {
    be <- geno
    if (is.null(snp_info)) snp_info <- be$snp_info
  } else if (is.matrix(geno) || is.data.frame(geno)) {
    if (is.null(snp_info))
      stop("snp_info must be supplied when geno is a matrix.", call. = FALSE)
    be <- read_geno(as.matrix(geno), format = "matrix", snp_info = snp_info)
    own_backend <- TRUE
  } else if (is.character(geno) && length(geno) == 1L) {
    be <- read_geno(geno)
    own_backend <- TRUE
    if (is.null(snp_info)) snp_info <- be$snp_info
  } else {
    stop("geno must be an LDxBlocks_backend, matrix, or file path.", call. = FALSE)
  }
  if (own_backend) on.exit(close_backend(be), add = TRUE)

  snp_info$CHR <- .norm_chr(snp_info$CHR)
  chrs <- if (!is.null(chr)) .norm_chr(as.character(chr)) else
    unique(snp_info$CHR)

  n_ind <- be$n_samples

  # -- rV2 whitening (genome-wide, computed once) ---------------------------
  # For rV2 we must whiten the genotype vectors before computing r2.
  # The whitening matrix V_inv_sqrt is n_ind x n_ind and applies to every
  # column read from any chromosome -- it is independent of which SNPs are loaded.
  V_inv_sqrt <- NULL
  if (method == "rV2") {
    .log("Computing VanRaden GRM for rV2 whitening (loads all SNPs once) ...")
    geno_all <- read_chunk(be, seq_len(be$n_snps))
    prep     <- prepare_geno(geno_all, method = "rV2",
                             kin_method = kin_method, verbose = verbose)
    V_inv_sqrt <- prep$V_inv_sqrt
    rm(geno_all, prep); gc(FALSE)
    .log("Whitening matrix ready (", n_ind, " x ", n_ind, ")")
  }

  # -- CORE HELPERS ----------------------------------------------------------

  # Apply optional whitening to a genotype column vector
  .prep_col <- function(x) {
    x[is.na(x)] <- mean(x, na.rm = TRUE)
    if (!is.null(V_inv_sqrt)) as.numeric(V_inv_sqrt %*% x) else x
  }

  # r2 between two ALREADY-PREPARED vectors.
  # Uses stats::cor() -- col_r2_cpp() operates on a full matrix column
  # against all others, which is used in the parametric section below.
  .col_r2 <- function(x, y) {
    r <- stats::cor(x, y)
    if (is.na(r)) 0 else r * r
  }

  # -- RANDOM SAMPLING: compute_r2_sparse_cpp + minimal read -----------------
  # STRATEGY:
  #   1. Sample pair (i,j) indices from SNP positions ONLY (no genotype I/O).
  #      For chromosomes with >200k SNPs, only a random subset of anchors is
  #      used to keep pair candidate enumeration O(n_pairs), not O(p^2).
  #   2. Collect UNIQUE SNP indices involved in sampled pairs.
  #   3. read_chunk() for ONLY those columns (n_unique << n_chr).
  #   4. Mean-impute + whiten the compact submatrix ONCE.
  #   5. compute_r2_sparse_cpp() on the compact submatrix with
  #      max_bp_dist=max_d to get all pairs within distance in one C++ call.
  #      This replaces the R pair-by-pair loop with Armadillo BLAS.
  #   6. Subsample to n_max pairs if the C++ result exceeds the budget.
  # Peak RAM: n_unique x n_ind x 8 bytes (typ. ~40 MB for 500k-SNP chr).
  .random_pairs <- function(chr_idx, pos_chr, n_max, max_d, seed_i) {
    p <- length(chr_idx)
    if (p < 2L) return(data.frame(dist_bp = integer(0), r2 = numeric(0)))
    set.seed(seed_i)

    # Step 1: sample anchor SNPs from positions to identify candidate pairs
    # without loading genotypes. For dense panels, sample a manageable
    # subset; for sparse panels use all SNPs.
    n_anchors  <- min(p, max(2L * ceiling(sqrt(n_max)), 500L))
    anchor_idx <- sort(sample(p, n_anchors))

    # Step 2: for each anchor, find partners within max_d bp
    pairs_ij <- lapply(anchor_idx, function(i) {
      j_cands <- which(pos_chr > pos_chr[i] & pos_chr <= pos_chr[i] + max_d)
      if (!length(j_cands)) return(NULL)
      cbind(i, sample(j_cands, min(length(j_cands), 5L)))
    })
    pairs_ij <- pairs_ij[!vapply(pairs_ij, is.null, logical(1L))]
    if (!length(pairs_ij))
      return(data.frame(dist_bp = integer(0), r2 = numeric(0)))
    mat_ij <- unique(do.call(rbind, pairs_ij))
    if (nrow(mat_ij) > n_max)
      mat_ij <- mat_ij[sample(nrow(mat_ij), n_max), , drop = FALSE]

    # Step 3: read ONLY the unique SNP columns needed
    uniq_local <- sort(unique(c(mat_ij[, 1L], mat_ij[, 2L])))
    global_idx <- chr_idx[uniq_local]
    .log("  random: reading ", length(global_idx), " / ", p,
         " unique SNPs for ~", nrow(mat_ij), " pairs")
    G_sub <- read_chunk(be, global_idx)   # n_ind x n_unique

    # Step 4: prepare columns ONCE (impute NAs, whiten if rV2)
    G_prep   <- apply(G_sub, 2L, .prep_col)
    col_sd   <- apply(G_prep, 2L, stats::sd, na.rm = TRUE)
    valid    <- col_sd > 1e-6
    rm(G_sub); gc(FALSE)

    # Step 5: compute_r2_sparse_cpp on the compact matrix in one C++ call.
    # This uses Armadillo BLAS internally -- far faster than an R pair loop.
    # bp positions are taken from pos_chr for the unique SNPs.
    bp_local  <- as.integer(pos_chr[uniq_local])
    # Keep only valid (non-monomorphic) columns
    valid_col <- which(valid)
    if (length(valid_col) < 2L)
      return(data.frame(dist_bp = integer(0), r2 = numeric(0)))
    G_valid   <- G_prep[, valid_col, drop = FALSE]
    bp_valid  <- bp_local[valid_col]

    sp <- compute_r2_sparse_cpp(
      G_valid,
      bp          = bp_valid,
      max_bp_dist = as.integer(max_d),
      threshold   = 0.0,
      n_threads   = n_threads
    )
    rm(G_prep, G_valid); gc(FALSE)

    if (length(sp$row) == 0L)
      return(data.frame(dist_bp = integer(0), r2 = numeric(0)))

    # Step 6: subsample if C++ returned more pairs than n_max
    dist_bp <- abs(bp_valid[sp$col] - bp_valid[sp$row])
    if (length(sp$r2) > n_max) {
      keep    <- sample(length(sp$r2), n_max)
      sp$r2   <- sp$r2[keep]
      dist_bp <- dist_bp[keep]
    }
    data.frame(dist_bp = as.integer(dist_bp),
               r2      = pmax(0, pmin(1, sp$r2)))
  }

  # -- SLIDING WINDOW: adaptive stride, bounded I/O calls ---------------------
  # Problem: window_snps=50 on a 500k-SNP chromosome = 10,000 read_chunk calls
  # and 12M pairs -- endless even with per-window reads.
  #
  # Solution: compute a stride that places ~n_pairs total pairs regardless of
  # chromosome length:
  #   pairs_per_window = window_snps * (window_snps - 1) / 2
  #   n_windows_needed = ceiling(n_pairs / pairs_per_window)
  #   stride = max(window_snps, floor(p / n_windows_needed))
  #
  # For p=500k, n_pairs=50k, window=50:
  #   pairs_per_window = 1225
  #   n_windows = ceil(50000/1225) = 41
  #   stride = floor(500000/41) = 12,195  -- only 41 read_chunk calls
  #
  # Density guard: for chromosomes > max_slide_snps, auto-switch to random.
  # Sliding window at WGS density is physically meaningless anyway --
  # consecutive SNPs 200 bp apart are always in perfect LD.
  .max_slide_snps <- 200000L   # ~200k SNPs: chip/low-density boundary

  .sliding_pairs <- function(chr_idx, pos_chr, win, max_d, n_max) {
    p <- length(chr_idx)
    if (p < 2L) return(data.frame(dist_bp = integer(0), r2 = numeric(0)))

    # Density guard
    if (p > .max_slide_snps) {
      .log("  chr has ", p, " SNPs (> ", .max_slide_snps,
           ") -- auto-switching sliding_window to random sampling")
      return(.random_pairs(chr_idx, pos_chr, n_max, max_d, seed + p))
    }

    # Compute adaptive stride
    pairs_per_win   <- win * (win - 1L) / 2L
    n_wins_needed   <- max(1L, ceiling(n_max / pairs_per_win))
    stride          <- max(win, floor(p / n_wins_needed))
    starts          <- seq(1L, p - 1L, by = stride)
    .log("  sliding_window: ", length(starts), " windows (stride=",
         stride, ", window=", win, ")")

    out_dist <- integer(0); out_r2 <- numeric(0)
    for (i in starts) {
      j_end   <- min(p, i + win - 1L)
      win_idx <- i:j_end
      G_win   <- read_chunk(be, chr_idx[win_idx])
      # Precompute all window columns once; apply variance guard
      G_wprep  <- apply(G_win, 2L, .prep_col)
      col_sd_w <- apply(G_wprep, 2L, stats::sd, na.rm = TRUE)
      rm(G_win); gc(FALSE)
      nw <- length(win_idx)
      for (a in seq_len(nw - 1L)) {
        if (col_sd_w[a] <= 1e-6) next
        for (b in seq(a + 1L, nw)) {
          if (col_sd_w[b] <= 1e-6) next
          d <- abs(pos_chr[win_idx[b]] - pos_chr[win_idx[a]])
          if (d > max_d) next
          out_dist <- c(out_dist, as.integer(d))
          out_r2   <- c(out_r2,   .col_r2(G_wprep[, a], G_wprep[, b]))
        }
      }
      rm(G_wprep)
    }
    data.frame(dist_bp = out_dist, r2 = pmax(0, pmin(1, out_r2)))
  }

  # -- MAIN LOOP: per chromosome ---------------------------------------------
  all_pairs    <- list()
  rn_pairs_list <- list()  # extra random pairs when sampling='both'
  n_pairs_used <- integer(length(chrs))
  names(n_pairs_used) <- chrs

  for (ci in seq_along(chrs)) {
    cb      <- chrs[ci]
    idx_chr <- which(snp_info$CHR == cb)
    p_chr   <- length(idx_chr)
    if (p_chr < 2L) { .log("chr ", cb, ": < 2 SNPs -- skipping"); next }
    # Sort by physical position -- REQUIRED for pair distance logic.
    # Both random and sliding_window assume pos_chr[i] < pos_chr[j]
    # implies index i comes before j. Unsorted input silently produces
    # wrong pair distances.
    pos_order <- order(snp_info$POS[idx_chr])
    idx_chr   <- idx_chr[pos_order]
    pos_chr   <- snp_info$POS[idx_chr]
    .log("chr ", cb, ": ", p_chr, " SNPs")

    pairs_chr <- switch(sampling,
                        random         = .random_pairs(idx_chr, pos_chr, n_pairs, max_dist, seed + ci),
                        sliding_window = .sliding_pairs(idx_chr, pos_chr, window_snps,
                                                        max_dist, n_pairs),
                        both = {
                          # 'both': sliding window provides the decay curve shape;
                          # random pairs are stored separately for parametric threshold.
                          # They are NOT merged into one pool -- that would mix two sampling
                          # designs and bias the curve shape estimate.
                          # Random pairs for this chromosome are stored in rn_pairs_list
                          # and used only for parametric threshold estimation below.
                          rn <- .random_pairs(idx_chr, pos_chr, n_pairs, max_dist, seed + ci + 1000L)
                          if (nrow(rn) > 0L) { rn$CHR <- cb; rn_pairs_list[[cb]] <- rn }
                          .sliding_pairs(idx_chr, pos_chr, window_snps, max_dist, n_pairs)
                        }
    )

    if (nrow(pairs_chr) > 0L) {
      pairs_chr$CHR    <- cb
      all_pairs[[cb]]  <- pairs_chr
      n_pairs_used[ci] <- nrow(pairs_chr)
      .log("  ", nrow(pairs_chr), " pairs computed")
    }
    # No gc() needed here -- .random_pairs and .sliding_pairs call gc() internally
  }

  pairs_df <- if (length(all_pairs) > 0L)
    do.call(rbind, all_pairs) else
      data.frame(CHR = character(), dist_bp = integer(), r2 = numeric())
  rownames(pairs_df) <- NULL

  # -- PARAMETRIC THRESHOLD: cross-chromosome pairs -------------------------
  # Never loads a full chromosome -- samples n_per_chr_side indices and
  # calls read_chunk() only for those SNPs.
  unlinked_r2 <- NULL; critical_r2_param <- NULL
  if (use_param && length(chrs) >= 2L) {
    .log("Parametric threshold: sampling ", n_unlinked,
         " cross-chromosome pairs ...")
    set.seed(seed + 9999L)
    # How many SNPs to sample per chromosome side per pair
    n_side    <- max(10L, ceiling(sqrt(n_unlinked / choose(length(chrs), 2L))))
    samp_idx  <- lapply(chrs, function(cb) {
      idx <- which(snp_info$CHR == cb)
      if (!length(idx)) return(NULL)
      sample(idx, min(n_side * 2L, length(idx)))
    })
    names(samp_idx) <- chrs

    unlinked_vals <- numeric(0)
    chr_pairs     <- utils::combn(chrs, 2L, simplify = FALSE)

    for (cp in chr_pairs) {
      idx1 <- samp_idx[[cp[1]]]; idx2 <- samp_idx[[cp[2]]]
      if (is.null(idx1) || is.null(idx2)) next
      # Minimal reads: only the sampled SNPs
      G1 <- read_chunk(be, idx1)   # n_ind x n_side
      G2 <- read_chunk(be, idx2)
      # Prepare both chunks once (impute + whiten)
      G1p <- apply(G1, 2L, .prep_col)
      G2p <- apply(G2, 2L, .prep_col)
      rm(G1, G2)
      # Use col_r2_cpp: for each sampled column of G1p, compute r2
      # against all columns of a combined [G1p | G2p] matrix, then
      # extract the cross-chromosome entries. This leverages Armadillo
      # BLAS instead of an R loop over individual cor() calls.
      k   <- min(ncol(G1p), ncol(G2p), 50L)
      s1  <- sample(ncol(G1p), k)
      n1  <- ncol(G1p)
      G12 <- cbind(G1p, G2p)   # combined matrix for col_r2_cpp
      for (i in seq_len(k)) {
        r2v <- as.numeric(col_r2_cpp(G12, s1[i]))
        # r2v[n1+1 .. ncol(G12)] are cross-chromosome r2 values
        unlinked_vals <- c(unlinked_vals, r2v[(n1 + 1L):ncol(G12)])
      }
      rm(G1p, G2p, G12); gc(FALSE)
    }

    unlinked_r2       <- unlinked_vals
    critical_r2_param <- stats::quantile(unlinked_vals, pctile / 100,
                                         na.rm = TRUE)
    .log("Parametric threshold (", pctile, "th pctile): ",
         round(critical_r2_param, 5))
  }

  critical_r2 <- if (use_param && !is.null(critical_r2_param))
    critical_r2_param else if (use_fixed) fixed_val else NULL

  # -- BIN INTO DECAY CURVE -------------------------------------------------
  .bin_curve <- function(df, n_b, max_d) {
    if (!nrow(df)) return(NULL)
    # Bin over the actual observed distance range so that all bins are
    # populated even when the data span << max_dist (e.g. small test panels).
    actual_max <- max(df$dist_bp, na.rm = TRUE)
    range_d <- min(max_d, if (actual_max > 0) actual_max else max_d)
    breaks <- seq(0, range_d, length.out = n_b + 1L)
    mids   <- (breaks[-1] + breaks[-(n_b + 1)]) / 2
    df$bin <- pmin(findInterval(df$dist_bp, breaks, rightmost.closed = TRUE), n_b)
    do.call(rbind, lapply(seq_len(n_b), function(b) {
      sub <- df[df$bin == b, , drop = FALSE]
      if (!nrow(sub)) return(NULL)
      data.frame(dist_bp   = mids[b],
                 r2_mean   = mean(sub$r2,             na.rm = TRUE),
                 r2_median = stats::median(sub$r2,    na.rm = TRUE),
                 n_pairs   = nrow(sub))
    }))
  }

  # -- FIT MODELS & COMPUTE DECAY DISTANCES ---------------------------------
  .hw_expect <- function(C, n, d_mb)
    (10 + C * d_mb) / ((2 + C * d_mb) * (11 + C * d_mb)) + 1 / n

  decay_curve_list <- list()
  decay_dist_list  <- list()
  model_params     <- list()

  for (cb in chrs) {
    df_chr <- pairs_df[pairs_df$CHR == cb, , drop = FALSE]
    if (!nrow(df_chr)) next
    curve <- .bin_curve(df_chr, n_bins, max_dist)
    if (is.null(curve)) next
    curve$CHR <- cb

    if (fit_model %in% c("loess", "both") && nrow(curve) >= 5L) {
      tryCatch({
        lo <- stats::loess(r2_mean ~ dist_bp, data = curve,
                           span = 0.3, degree = 2)
        curve$r2_loess <- pmax(0, stats::predict(lo, curve$dist_bp))
      }, error = function(e) {
        warning("[ld_decay] LOESS fit failed for chr ", cb,
                ": ", conditionMessage(e), call. = FALSE)
      })
    }

    if (fit_model %in% c("nonlinear", "both") && nrow(df_chr) >= 10L) {
      tryCatch({
        fit <- stats::nls(
          r2 ~ .hw_expect(C, n_ind, dist_bp / 1e6),
          data    = df_chr[df_chr$dist_bp > 0L, ],
          start   = list(C = 1),
          control = stats::nls.control(maxiter = 100, warnOnly = TRUE)
        )
        C_hat <- stats::coef(fit)[["C"]]
        model_params[[cb]] <- list(C = C_hat, n_ind = n_ind)
        curve$r2_hw <- pmax(0, .hw_expect(C_hat, n_ind, curve$dist_bp / 1e6))
      }, error = function(e) {
        warning("[ld_decay] Hill-Weir nonlinear fit failed for chr ", cb,
                ": ", conditionMessage(e),
                "\n  -> LOESS or binned mean will be used as fallback.",
                call. = FALSE)
      })
    }

    decay_curve_list[[cb]] <- curve

    if (!is.null(critical_r2)) {
      r2_col <- if ("r2_loess" %in% names(curve)) "r2_loess" else
        if ("r2_hw"    %in% names(curve)) "r2_hw"    else "r2_mean"
      smooth <- curve[[r2_col]]
      above  <- which(!is.na(smooth) & smooth >= critical_r2)
      below  <- which(!is.na(smooth) & smooth <  critical_r2)

      # Censored = threshold never meaningfully crossed:
      # either the smooth curve never rises above threshold (no decay)
      # or it never drops below (threshold too low for this data).
      censored <- length(above) == 0L || length(below) == 0L

      if (censored) {
        # Decay distance not reached: either r2 never exceeds threshold
        # (threshold too high) or never drops below (threshold too low).
        # Return max observed distance as a lower/upper bound.
        d_bp <- max(curve$dist_bp, na.rm = TRUE)
        warning("[ld_decay] chr ", cb, ": LD does not decay below threshold ",
                round(critical_r2, 4), " within max_dist=",
                format(max_dist, big.mark=","),
                " bp. decay_dist is censored (lower bound).", call. = FALSE)
      } else {
        # Linear interpolation between last-above and first-below bins
        # for sub-bin precision (avoids coarse bin-midpoint rounding).
        last_above <- max(above)
        first_below <- min(below)
        if (last_above < first_below) {
          # Normal case: monotone crossing between two consecutive bins
          x1 <- curve$dist_bp[last_above];  y1 <- smooth[last_above]
          x2 <- curve$dist_bp[first_below]; y2 <- smooth[first_below]
          # Linear interpolation: solve y = y1 + (y2-y1)/(x2-x1)*(x-x1) for x
          d_bp <- x1 + (critical_r2 - y1) * (x2 - x1) / (y2 - y1)
        } else {
          # Non-monotone curve (e.g. centromeric bump); use first below
          d_bp <- curve$dist_bp[first_below]
        }
      }

      decay_dist_list[[cb]] <- data.frame(
        CHR              = cb,
        decay_dist_bp    = as.integer(round(d_bp)),
        decay_dist_kb    = round(d_bp / 1000, 2),
        threshold_used   = critical_r2,
        r2_col_used      = r2_col,
        censored         = censored,
        stringsAsFactors = FALSE
      )
    }
  }

  decay_curve   <- if (length(decay_curve_list))
    do.call(rbind, decay_curve_list) else NULL
  rownames(decay_curve) <- NULL

  decay_dist    <- if (length(decay_dist_list))
    do.call(rbind, decay_dist_list) else NULL
  rownames(decay_dist) <- NULL

  decay_dist_genome <- if (!is.null(decay_dist))
    stats::median(decay_dist$decay_dist_bp, na.rm = TRUE) else NULL

  if (verbose && !is.null(decay_dist)) {
    .log("Decay distances:")
    for (i in seq_len(nrow(decay_dist)))
      .log("  chr ", decay_dist$CHR[i], ": ",
           format(decay_dist$decay_dist_bp[i], big.mark = ","), " bp  (",
           decay_dist$decay_dist_kb[i], " kb)")
    .log("Genome-wide median: ",
         format(round(decay_dist_genome), big.mark = ","), " bp")
  }

  structure(list(
    pairs              = pairs_df,
    decay_curve        = decay_curve,
    decay_dist         = decay_dist,
    decay_dist_genome  = decay_dist_genome,
    critical_r2        = critical_r2,
    critical_r2_fixed  = if (use_fixed) fixed_val else 0.1,
    critical_r2_param  = critical_r2_param,
    unlinked_r2        = unlinked_r2,
    model_params       = if (length(model_params)) model_params else NULL,
    n_pairs_used       = n_pairs_used,
    method             = method,
    sampling           = sampling,
    pctile             = pctile,
    call               = cl
  ), class = c("LDxBlocks_decay", "list"))
}


#' @export
print.LDxBlocks_decay <- function(x, ...) {
  cat("LDxBlocks LD Decay Analysis\n")
  cat("  Method      :", x$method, "\n")
  cat("  Sampling    :", x$sampling, "\n")
  cat("  Pairs total :", nrow(x$pairs), "\n")
  if (!is.null(x$critical_r2_fixed))
    cat("  r2 threshold (fixed)     :", round(x$critical_r2_fixed, 4), "\n")
  if (!is.null(x$critical_r2_param))
    cat("  r2 threshold (parametric):", round(x$critical_r2_param, 4),
        paste0("[", x$pctile, "th pctile unlinked]\n"))
  if (!is.null(x$decay_dist)) {
    cat("  Decay distances:\n")
    for (i in seq_len(nrow(x$decay_dist)))
      cat("    chr", x$decay_dist$CHR[i], ":",
          x$decay_dist$decay_dist_kb[i], "kb\n")
    cat("  Genome-wide median :",
        round(x$decay_dist_genome / 1000, 1), "kb\n")
  }
  invisible(x)
}


#' Plot LD Decay Curve
#'
#' @description
#' Plots r\eqn{^2} vs physical distance from \code{\link{compute_ld_decay}}.
#'
#' @param x \code{LDxBlocks_decay} object.
#' @param plot_points Logical. Plot raw pair r\eqn{^2} as semi-transparent
#'   points. Default \code{FALSE}.
#' @param plot_threshold Logical. Horizontal line at critical r\eqn{^2}.
#'   Default \code{TRUE}.
#' @param plot_decay_dist Logical. Vertical lines at per-chromosome decay
#'   distances. Default \code{TRUE}.
#' @param max_dist_kb Numeric. X-axis limit in kb. Default \code{NULL} (data
#'   range).
#' @param facet Logical. One panel per chromosome. Default \code{FALSE}.
#' @param ... Ignored.
#' @return A \code{ggplot2} object.
#' @export
plot_ld_decay <- function(x, plot_points = FALSE, plot_threshold = TRUE,
                          plot_decay_dist = TRUE, max_dist_kb = NULL,
                          facet = FALSE, ...) {
  if (!inherits(x, "LDxBlocks_decay"))
    stop("x must be an LDxBlocks_decay object.", call. = FALSE)
  if (!requireNamespace("ggplot2", quietly = TRUE))
    stop("ggplot2 required: install.packages('ggplot2')", call. = FALSE)

  curve <- x$decay_curve
  if (is.null(curve) || !nrow(curve))
    stop("No decay curve -- rerun with fit_model != 'none'.", call. = FALSE)

  has_fitted <- "r2_loess" %in% names(curve) || "r2_hw" %in% names(curve)
  if (!has_fitted)
    stop("No fitted decay curve: fit_model was 'none'. ",
         "Rerun compute_ld_decay() with fit_model = 'loess' or 'nonlinear'.",
         call. = FALSE)

  curve$dist_kb <- curve$dist_bp / 1000
  r2_col        <- if ("r2_loess" %in% names(curve)) "r2_loess" else
    if ("r2_hw"    %in% names(curve)) "r2_hw"    else "r2_mean"
  curve$r2_plot <- curve[[r2_col]]

  # -- Numeric chromosome ordering ----------------------------------------------
  # Convert CHR to an ordered factor so ggplot2 uses numeric order
  # (1, 2, 3, ..., 12) rather than lexicographic order (1, 10, 11, ..., 9).
  all_chrs  <- unique(as.character(curve$CHR))
  chr_num   <- suppressWarnings(as.integer(all_chrs))
  chr_order <- if (!any(is.na(chr_num))) {
    all_chrs[order(chr_num)]          # pure numeric chromosomes: 1,2,...,12
  } else {
    # Mixed (e.g. 1A, 2B): sort numerically on leading digits, alpha on suffix
    chr_lead <- suppressWarnings(as.integer(sub("[^0-9].*", "", all_chrs)))
    if (!any(is.na(chr_lead)))
      all_chrs[order(chr_lead, all_chrs)]
    else
      sort(all_chrs)
  }
  curve$CHR <- factor(curve$CHR, levels = chr_order)
  if (!is.null(x$pairs) && nrow(x$pairs) > 0L)
    x$pairs$CHR <- factor(x$pairs$CHR, levels = chr_order)
  if (!is.null(x$decay_dist))
    x$decay_dist$CHR <- factor(x$decay_dist$CHR, levels = chr_order)

  xlim <- c(0, if (!is.null(max_dist_kb)) max_dist_kb else
    max(curve$dist_kb, na.rm = TRUE))

  # -- Shape cycling: genome overlay uses colour + shape so chromosomes are
  # distinguishable even in greyscale or for colourblind readers.
  # ggplot2 has 6 default solid shapes (15-20); cycle them for > 6 chromosomes.
  n_chr      <- length(chr_order)
  shape_vals <- rep(c(15L, 16L, 17L, 18L, 8L, 7L), length.out = n_chr)
  names(shape_vals) <- chr_order

  # Genome overlay (facet=FALSE): colour + shape on lines + points for each chr.
  # Per-chromosome facet (facet=TRUE): all lines same colour (no distinction
  # needed per panel), legend suppressed to save space.
  if (!facet) {
    p <- ggplot2::ggplot(curve,
                         ggplot2::aes(x = dist_kb, y = r2_plot,
                                      colour = CHR, shape = CHR, group = CHR))
  } else {
    p <- ggplot2::ggplot(curve,
                         ggplot2::aes(x = dist_kb, y = r2_plot, group = CHR))
  }

  if (plot_points && nrow(x$pairs) > 0L) {
    pts <- x$pairs; pts$dist_kb <- pts$dist_bp / 1000
    pts <- pts[pts$dist_kb <= xlim[2], ]
    if (nrow(pts) > 10000L) pts <- pts[sample(nrow(pts), 10000L), ]
    if (!facet) {
      p <- p + ggplot2::geom_point(
        data = pts,
        ggplot2::aes(x = dist_kb, y = r2, colour = CHR),
        alpha = 0.1, size = 0.4, inherit.aes = FALSE)
    } else {
      p <- p + ggplot2::geom_point(
        data = pts,
        ggplot2::aes(x = dist_kb, y = r2),
        colour = "steelblue", alpha = 0.1, size = 0.4, inherit.aes = FALSE)
    }
  }

  if (!facet) {
    # Genome overlay: line + shape markers spaced along curve so each chr
    # is identifiable without colour alone. geom_line for the curve,
    # geom_point at ~10 evenly-spaced positions per chr for the shape marker.
    n_mark <- min(10L, max(3L, floor(nrow(curve) / n_chr)))
    mark_idx <- unique(round(seq(1, nrow(curve), length.out = n_mark * n_chr)))
    curve_marks <- curve[mark_idx, , drop = FALSE]

    p <- p +
      ggplot2::geom_line(linewidth = 0.8, na.rm = TRUE) +
      ggplot2::geom_point(
        data = curve_marks,
        ggplot2::aes(x = dist_kb, y = r2_plot, colour = CHR, shape = CHR),
        size = 2.2, na.rm = TRUE) +
      ggplot2::scale_shape_manual(values = shape_vals, name = "Chromosome") +
      ggplot2::scale_colour_discrete(name = "Chromosome") +
      ggplot2::guides(
        colour = ggplot2::guide_legend(ncol = 2L, override.aes = list(size = 2.5)),
        shape  = ggplot2::guide_legend(ncol = 2L)
      )
  } else {
    # Per-chromosome facet: single colour line, no legend needed
    p <- p +
      ggplot2::geom_line(colour = "steelblue", linewidth = 0.8, na.rm = TRUE)
  }

  p <- p +
    ggplot2::coord_cartesian(xlim = xlim, ylim = c(0, NA)) +
    ggplot2::labs(x = "Physical distance (kb)", y = expression(r^2),
                  title = paste("LD Decay -", x$method)) +
    ggplot2::theme_minimal(base_size = 11)

  # Suppress legend entirely for per-chromosome faceted plot
  if (facet) p <- p + ggplot2::theme(legend.position = "none")

  if (plot_threshold && !is.null(x$critical_r2)) {
    p <- p +
      ggplot2::geom_hline(yintercept = x$critical_r2,
                          linetype = "dashed", colour = "grey40",
                          linewidth = 0.5) +
      ggplot2::annotate("text", x = xlim[2] * 0.02, y = x$critical_r2,
                        label = paste0("threshold = ",
                                       round(x$critical_r2, 3)),
                        hjust = 0, vjust = -0.4, size = 3,
                        colour = "grey40")
  }

  if (plot_decay_dist && !is.null(x$decay_dist)) {
    dd <- x$decay_dist; dd$dist_kb <- dd$decay_dist_bp / 1000
    dd <- dd[dd$dist_kb <= xlim[2], ]
    if (nrow(dd)) {
      if (!facet) {
        p <- p + ggplot2::geom_vline(
          data = dd,
          ggplot2::aes(xintercept = dist_kb, colour = CHR),
          linetype = "dotted", linewidth = 0.5, alpha = 0.7,
          inherit.aes = FALSE)
      } else {
        p <- p + ggplot2::geom_vline(
          data = dd,
          ggplot2::aes(xintercept = dist_kb),
          colour = "firebrick", linetype = "dotted",
          linewidth = 0.5, alpha = 0.8, inherit.aes = FALSE)
      }
    }
  }

  if (facet) p <- p + ggplot2::facet_wrap(~CHR, scales = "free_x")
  p
}

# Internal: strip chr/Chr/CHR prefix (shared with other modules)
.norm_chr <- function(x) sub("^[Cc][Hh][Rr]", "", as.character(x))
