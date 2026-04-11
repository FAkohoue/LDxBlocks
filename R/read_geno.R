# -----------------------------------------------------------------------------
# read_geno.R  -  Unified genotype reader for LDxBlocks
# -----------------------------------------------------------------------------
#
# Supported formats (auto-detected from extension or set via format =):
#
#   "numeric"  .csv / .txt   One row per SNP; cols: SNP, CHR, POS, REF, ALT,
#                            then one col per sample with values in {0,1,2,NA}
#
#   "hapmap"   .hmp.txt      Standard 11-col HapMap header then sample cols
#                            with two-character nucleotide calls (AA/AT/NN...)
#
#   "vcf"      .vcf          Standard VCF v4.2; GT field parsed to 0/1/2/NA;
#              .vcf.gz       phased (0|1) and unphased (0/1) both accepted;
#                            multi-allelic -> first ALT; missing ./. -> NA
#
#   "gds"      .gds          SeqArray GDS file (Bioconductor)
#
#   "bed"      .bed + .bim   PLINK binary BED (requires .bim and .fam)
#              + .fam
#
#   "matrix"   (in-memory)   Plain R numeric matrix already in memory
#
# Returns an object of class "LDxBlocks_backend" with a unified interface:
#   backend$n_samples  integer
#   backend$n_snps     integer
#   backend$sample_ids character vector
#   backend$snp_info   data.frame(SNP, CHR, POS, REF, ALT)
#   read_chunk(backend, col_idx)  ->  n x length(col_idx) numeric matrix
#   close_backend(backend)        ->  releases file handles
#
# Chromosome names are normalised: "chr1" -> "1" consistently.
# -----------------------------------------------------------------------------


# -- Internal: normalise chromosome names --------------------------------------
.norm_chr <- function(x) {
  x <- as.character(x)
  sub("^[Cc][Hh][Rr]", "", x)   # strip leading "chr" / "Chr" / "CHR"
}


# -- Internal: detect format from file path ------------------------------------
.detect_format <- function(path) {
  p <- tolower(path)
  if (grepl("\\.hmp\\.txt$", p))                        return("hapmap")
  if (grepl("\\.vcf\\.gz$", p) || grepl("\\.vcf$", p)) return("vcf")
  if (grepl("\\.gds$", p))                              return("gds")
  if (grepl("\\.bed$", p))                              return("bed")
  if (grepl("\\.csv$", p) || grepl("\\.txt$", p))       return("numeric")
  stop("Cannot detect genotype format from extension of '", path,
       "'. Set format = explicitly.")
}


# -----------------------------------------------------------------------------
# Main entry point
# -----------------------------------------------------------------------------

#' Read Genotype Data into an LDxBlocks Backend
#'
#' @description
#' Reads genotype data from disk (or accepts an in-memory matrix) and returns
#' a unified backend object that all \pkg{LDxBlocks} functions can consume.
#' The format is auto-detected from the file extension unless overridden with
#' \code{format}.
#'
#' @section Supported formats:
#' \describe{
#'   \item{\code{"numeric"}}{CSV or TXT, one row per SNP. Required columns:
#'     \code{SNP}, \code{CHR}, \code{POS}, \code{REF}, \code{ALT}. Remaining
#'     columns are sample dosages in \{0, 1, 2, NA\}. Extension: \code{.csv},
#'     \code{.txt}.}
#'   \item{\code{"hapmap"}}{Standard 11-column HapMap header followed by sample
#'     columns with two-character nucleotide calls (\code{AA}, \code{AT},
#'     \code{NN} for missing). Extension: \code{.hmp.txt}.}
#'   \item{\code{"vcf"}}{VCF v4.2. Both phased (\code{0|1}) and unphased
#'     (\code{0/1}) GT fields are accepted. Multi-allelic sites use first ALT.
#'     Missing (\code{./.}) becomes \code{NA}. Extension: \code{.vcf},
#'     \code{.vcf.gz}.}
#'   \item{\code{"gds"}}{SeqArray GDS file. Requires
#'     \code{BiocManager::install("SeqArray")}. Extension: \code{.gds}.}
#'   \item{\code{"bed"}}{PLINK binary BED. Companion \code{.bim} and
#'     \code{.fam} files must exist at the same path stem. Extension:
#'     \code{.bed}.}
#'   \item{\code{"matrix"}}{In-memory R numeric matrix (individuals x SNPs,
#'     0/1/2). Must supply \code{snp_info} separately.}
#' }
#'
#' @param path Character. Path to the genotype file, OR an R numeric matrix
#'   when \code{format = "matrix"}.
#' @param format Character or \code{NULL}. One of \code{"numeric"},
#'   \code{"hapmap"}, \code{"vcf"}, \code{"gds"}, \code{"bed"},
#'   \code{"matrix"}. \code{NULL} (default) auto-detects from extension.
#' @param snp_info Data frame with columns \code{SNP}, \code{CHR}, \code{POS}
#'   (and optionally \code{REF}, \code{ALT}). Required only when
#'   \code{format = "matrix"}.
#' @param sample_ids Character vector. Override sample IDs extracted from the
#'   file. Length must equal number of samples.
#' @param sep Character. Field separator for \code{"numeric"} format.
#' @param gds_cache Character or \code{NULL}. Path where a GDS cache file
#'   should be written when \code{format = "vcf"} and SeqArray is available.
#'   If \code{NULL} (default), the GDS file is placed next to the VCF with a
#'   \code{.gds} extension. Set to \code{FALSE} to disable auto-conversion
#'   and read the VCF fully into memory instead. When the cache file already
#'   exists it is reused without re-converting (fast subsequent calls).
#'   Ignored for all non-VCF formats.
#'   Default \code{","}.
#' @param na_strings Character vector. Strings treated as NA. Default
#'   \code{c("NA", "N", "NN", "./.", ".", "")}.
#' @param verbose Logical. Print progress messages.
#'
#' @return An object of class \code{"LDxBlocks_backend"} with elements:
#'   \code{type}, \code{n_samples}, \code{n_snps}, \code{sample_ids},
#'   \code{snp_info} (data.frame SNP/CHR/POS/REF/ALT). Use
#'   \code{\link{read_chunk}} to extract genotype slices and
#'   \code{\link{close_backend}} to release file handles.
#'
#' @seealso \code{\link{read_chunk}}, \code{\link{close_backend}},
#'   \code{\link{run_Big_LD_all_chr}}
#'
#' @examples
#' \donttest{
#' # In-memory matrix (existing workflow, unchanged)
#' G   <- matrix(sample(0:2, 100*500, replace=TRUE), 100, 500)
#' rownames(G) <- paste0("ind", 1:100)
#' colnames(G) <- paste0("rs",  1:500)
#' info <- data.frame(SNP=colnames(G), CHR=rep(1:5, each=100),
#'                    POS=rep(seq(1e4,5e6,length.out=100),5))
#' be <- read_geno(G, format="matrix", snp_info=info)
#' be$n_snps   # 500
#' mat <- read_chunk(be, 1:20)
#' dim(mat)    # 100 x 20
#' close_backend(be)
#' }
#'
#' @export
read_geno <- function(
    path,
    format      = NULL,
    snp_info    = NULL,
    sample_ids  = NULL,
    sep         = ",",
    na_strings  = c("NA", "N", "NN", "./.", ".", ""),
    gds_cache   = NULL,
    verbose     = FALSE
) {
  # -- Dispatch on format ------------------------------------------------------
  if (is.matrix(path) || is.data.frame(path)) {
    fmt <- "matrix"
  } else {
    fmt <- if (!is.null(format)) tolower(format) else .detect_format(path)
  }

  # -- Auto-convert VCF or HapMap to GDS for streaming access ---------------
  # Triggered when SeqArray is available and gds_cache != FALSE.
  # GDS is placed next to the source file (same directory, .gds extension)
  # unless the user supplies a custom path via gds_cache. Subsequent calls
  # reuse the cached GDS without re-converting.
  if (fmt %in% c("vcf", "hapmap") && !isFALSE(gds_cache) &&
      requireNamespace("SeqArray", quietly = TRUE)) {

    cache_path <- if (is.null(gds_cache)) {
      sub("\\.(vcf(\\.gz)?|hmp\\.txt)$", ".gds", path, ignore.case = TRUE)
    } else {
      as.character(gds_cache)
    }

    if (!file.exists(cache_path)) {
      message("[read_geno] Converting ", toupper(fmt),
              " to GDS for streaming access.")
      message("[read_geno] GDS cache: ", cache_path)
      message("[read_geno] This runs once; subsequent calls reuse the cache.")

      if (fmt == "vcf") {
        SeqArray::seqVCF2GDS(
          vcf.fn         = path,
          out.fn         = cache_path,
          storage.option = "ZIP_RA",
          verbose        = isTRUE(verbose)
        )
      } else {
        # HapMap: SeqArray has no direct reader.
        # Read via data.table, decode nucleotide calls to dosage, write GDS.
        .hapmap_to_gds(path, cache_path, na_strings, verbose)
      }

    } else if (isTRUE(verbose)) {
      message("[read_geno] Reusing GDS cache: ", cache_path)
    }

    fmt  <- "gds"
    path <- cache_path
  }


  be <- switch(fmt,
               numeric = .read_numeric(path, sep, na_strings, verbose),
               hapmap  = .read_hapmap(path, na_strings, verbose),
               vcf     = .read_vcf(path, na_strings, verbose),
               gds     = .read_gds(path, verbose),
               bed     = .read_bed(path, verbose),
               matrix  = .wrap_matrix(path, snp_info, verbose),
               stop("Unknown format '", fmt, "'. Choose: numeric, hapmap, vcf, gds, bed, matrix.")
  )

  # -- Override sample IDs if supplied ----------------------------------------
  if (!is.null(sample_ids)) {
    if (length(sample_ids) != be$n_samples)
      stop("length(sample_ids) [", length(sample_ids), "] != n_samples [",
           be$n_samples, "]")
    be$sample_ids <- as.character(sample_ids)
  }

  # -- Normalise chromosome names ----------------------------------------------
  be$snp_info$CHR <- .norm_chr(be$snp_info$CHR)

  structure(be, class = "LDxBlocks_backend")
}


# -----------------------------------------------------------------------------
# Format readers
# -----------------------------------------------------------------------------

# -- 1. Numeric dosage (CSV/TXT) -----------------------------------------------
.read_numeric <- function(path, sep, na_strings, verbose,
                          chunk_rows = 50000L) {
  if (isTRUE(verbose)) cat("[read_geno] Reading numeric dosage (chunked):", path, "\n")
  .require_pkg("data.table", "numeric dosage reader")

  # -- Pass 1: header-only scan (zero data rows) ----------------------------
  # Determines column names and exact row count without loading any data.
  # Reading only column 1 for the row count costs essentially nothing.
  meta_cols <- c("SNP","CHR","POS","REF","ALT")
  hdr       <- data.table::fread(path, sep = sep, nrows = 0L,
                                 check.names = FALSE, data.table = TRUE)
  all_cols  <- names(hdr)
  missing   <- setdiff(c("SNP","CHR","POS"), all_cols)
  if (length(missing))
    stop("numeric dosage file missing: ", paste(missing, collapse=", "))
  if (!"REF" %in% all_cols) all_cols <- c(all_cols, "REF")
  if (!"ALT" %in% all_cols) all_cols <- c(all_cols, "ALT")

  samp_cols <- setdiff(names(hdr), meta_cols)
  if (!length(samp_cols)) stop("No sample columns found in numeric dosage file.")
  n_samp <- length(samp_cols)

  n_snp <- nrow(data.table::fread(path, sep = sep, select = 1L,
                                  header = TRUE, data.table = TRUE))
  if (isTRUE(verbose))
    cat("[read_geno]  ", n_snp, "SNPs x", n_samp, "samples\n")

  # -- Pre-allocate output: SNPs as columns, samples as rows ----------------
  # This is the ONLY large allocation. Peak RAM = pre-alloc + one chunk,
  # not 2x the full file (which naive fread -> as.matrix produces).
  geno_mat           <- matrix(NA_real_, nrow = n_samp, ncol = n_snp)
  rownames(geno_mat) <- samp_cols
  snp_info_list      <- vector("list", ceiling(n_snp / chunk_rows))
  snps_filled        <- 0L; chunk_idx <- 0L

  # -- Pass 2: fill chunks --------------------------------------------------
  repeat {
    if (snps_filled >= n_snp) break
    chunk_idx <- chunk_idx + 1L

    chunk <- data.table::fread(
      path, sep = sep,
      nrows      = chunk_rows,
      skip       = snps_filled + 1L,  # +1 skips the header line
      header     = FALSE,
      col.names  = names(hdr),
      check.names = FALSE,
      na.strings = na_strings,
      data.table = TRUE
    )
    if (!nrow(chunk)) break

    c_start <- snps_filled + 1L
    c_end   <- snps_filled + nrow(chunk)

    # Numeric block: transpose directly into pre-allocated slice
    blk <- as.matrix(chunk[, samp_cols, with = FALSE])  # n_SNPs x n_samples
    storage.mode(blk) <- "numeric"
    geno_mat[, c_start:c_end] <- t(blk)  # transpose: fill as n_samples x n_SNPs
    rm(blk)

    # REF/ALT fill if columns absent in source file
    if (!"REF" %in% names(chunk)) chunk[, REF := NA_character_]
    if (!"ALT" %in% names(chunk)) chunk[, ALT := NA_character_]

    snp_info_list[[chunk_idx]] <- data.frame(
      SNP = as.character(chunk$SNP),
      CHR = .norm_chr(chunk$CHR),
      POS = as.integer(chunk$POS),
      REF = as.character(chunk$REF),
      ALT = as.character(chunk$ALT),
      stringsAsFactors = FALSE
    )
    snps_filled <- c_end
    rm(chunk); gc(FALSE)  # release allocator pressure between chunks
  }

  snp_info           <- do.call(rbind, snp_info_list[seq_len(chunk_idx)])
  colnames(geno_mat) <- snp_info$SNP

  .make_matrix_backend(geno_mat, snp_info, samp_cols, "numeric")
}


# -- 2. HapMap ----------------------------------------------------------------
.read_hapmap <- function(path, na_strings, verbose) {
  if (isTRUE(verbose)) cat("[read_geno] Reading HapMap:", path, "\n")
  .require_pkg("data.table", "HapMap reader")

  dt <- data.table::fread(path, na.strings = na_strings,
                          showProgress = verbose, data.table = TRUE)

  # Standard HapMap header columns (case-insensitive match)
  hdr_map <- c("rs#" = "SNP", "alleles" = "alleles", "chrom" = "CHR",
               "pos" = "POS", "strand" = "strand",
               "assembly#" = "assembly", "center" = "center",
               "protlsid" = "protLSID", "assaylsid" = "assayLSID",
               "panellsid" = "panelLSID", "qccode" = "QCcode")
  old_names <- tolower(names(dt))
  for (i in seq_along(hdr_map)) {
    idx <- which(old_names == names(hdr_map)[i])
    if (length(idx)) data.table::setnames(dt, names(dt)[idx], hdr_map[i])
  }

  required <- c("SNP", "CHR", "POS")
  if (!all(required %in% names(dt)))
    stop("HapMap file missing columns: ", paste(setdiff(required, names(dt)), collapse = ", "))

  # Parse alleles column "A/T" -> REF=A, ALT=T
  if ("alleles" %in% names(dt)) {
    parts <- strsplit(as.character(dt$alleles), "/", fixed = TRUE)
    dt[, REF := vapply(parts, function(x) x[1L], character(1L))]
    dt[, ALT := vapply(parts, function(x) if (length(x) > 1L) x[2L] else NA_character_, character(1L))]
  } else {
    dt[, REF := NA_character_]; dt[, ALT := NA_character_]
  }

  meta_cols <- c("SNP","CHR","POS","alleles","strand","assembly","center",
                 "protLSID","assayLSID","panelLSID","QCcode","REF","ALT")
  samp_cols <- setdiff(names(dt), meta_cols)
  if (length(samp_cols) == 0) stop("No sample columns found in HapMap file.")

  snp_info <- as.data.frame(dt[, .(SNP, CHR, POS, REF, ALT)])

  # Decode two-character nucleotide calls -> dosage
  # Strategy: heterozygous = 1, homozygous ref = 0, homozygous alt = 2, NN = NA
  ref_alleles <- snp_info$REF
  alt_alleles <- snp_info$ALT

  geno_raw <- as.matrix(dt[, ..samp_cols])   # SNPs x samples (character)
  n_snp    <- nrow(geno_raw)
  n_samp   <- ncol(geno_raw)
  geno_dos <- matrix(NA_real_, n_snp, n_samp,
                     dimnames = list(snp_info$SNP, samp_cols))

  for (j in seq_len(n_snp)) {
    ref <- ref_alleles[j]; alt <- alt_alleles[j]
    for (i in seq_len(n_samp)) {
      call <- geno_raw[j, i]
      if (is.na(call) || call %in% c("NN", "N", "--", "NA")) {
        geno_dos[j, i] <- NA_real_
      } else {
        a1 <- substr(call, 1L, 1L); a2 <- substr(call, 2L, 2L)
        dos <- sum(c(a1, a2) == alt, na.rm = TRUE)
        geno_dos[j, i] <- dos
      }
    }
  }

  # Transpose to samples x SNPs
  geno_mat <- t(geno_dos)
  .make_matrix_backend(geno_mat, snp_info, samp_cols, "hapmap")
}


# -- 3. VCF / VCF.GZ ---------------------------------------------------------
.read_vcf <- function(path, na_strings, verbose) {
  if (isTRUE(verbose)) cat("[read_geno] Reading VCF:", path, "\n")
  .require_pkg("data.table", "VCF reader")

  # Read raw lines, skip meta (##) lines
  if (grepl("\\.gz$", path)) {
    con  <- gzcon(file(path, "rb"))
    lines <- readLines(con); close(con)
  } else {
    lines <- readLines(path)
  }

  meta_lines <- grepl("^##", lines)
  hdr_line   <- which(grepl("^#CHROM", lines))[1L]
  if (is.na(hdr_line)) stop("VCF file has no #CHROM header line.")

  data_lines <- lines[(hdr_line + 1L):length(lines)]
  data_lines <- data_lines[nchar(data_lines) > 0]

  # Parse header
  hdr      <- strsplit(sub("^#", "", lines[hdr_line]), "\t")[[1L]]
  fixed    <- c("CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT")
  samp_hdr <- hdr[(length(fixed) + 1L):length(hdr)]
  if (length(samp_hdr) == 0) stop("VCF has no sample columns.")

  # Parse body with data.table for speed
  dt <- data.table::fread(text = c(paste(hdr, collapse="\t"), data_lines),
                          sep = "\t", header = TRUE, na.strings = na_strings,
                          showProgress = verbose, data.table = TRUE)
  data.table::setnames(dt, names(dt)[1L], "CHROM")   # handle leading #

  snp_info <- data.frame(
    SNP = ifelse(is.na(dt$ID) | dt$ID == ".", paste0(dt$CHROM,":",dt$POS), dt$ID),
    CHR = as.character(dt$CHROM),
    POS = as.integer(dt$POS),
    REF = as.character(dt$REF),
    ALT = as.character(dt$ALT),
    stringsAsFactors = FALSE
  )

  # Parse GT field -> dosage
  # GT is always first sub-field in FORMAT column
  n_snp  <- nrow(dt)
  n_samp <- length(samp_hdr)
  geno_dos <- matrix(NA_real_, n_snp, n_samp,
                     dimnames = list(snp_info$SNP, samp_hdr))

  for (i in seq_len(n_samp)) {
    gt_raw <- as.character(dt[[samp_hdr[i]]])
    gt_raw <- sub(":.*", "", gt_raw)   # keep only GT sub-field
    geno_dos[, i] <- .parse_gt(gt_raw)
  }

  geno_mat <- t(geno_dos)   # samples x SNPs
  .make_matrix_backend(geno_mat, snp_info, samp_hdr, "vcf")
}

# Parse GT strings ("0/0", "0/1", "1|1", "./.", "0|1") -> dosage 0/1/2/NA
.parse_gt <- function(gt) {
  gt  <- as.character(gt)
  sep <- ifelse(grepl("|", gt, fixed = TRUE), "|", "/")
  vapply(gt, function(g) {
    if (is.na(g) || g %in% c(".", "./.", ".|.")) return(NA_real_)
    s <- strsplit(g, "[|/]")[[1L]]
    a <- suppressWarnings(as.integer(s))
    if (any(is.na(a))) return(NA_real_)
    # dosage = number of ALT alleles (any non-zero integer treated as ALT)
    as.numeric(sum(a > 0L))
  }, numeric(1L), USE.NAMES = FALSE)
}


# -- 4. GDS (SeqArray) --------------------------------------------------------
.read_gds <- function(path, verbose) {
  if (isTRUE(verbose)) cat("[read_geno] Opening GDS:", path, "\n")
  .require_pkg("SeqArray", "GDS reader",
               install = "BiocManager::install('SeqArray')")

  gds <- SeqArray::seqOpen(path, readonly = TRUE)

  sample_ids <- SeqArray::seqGetData(gds, "sample.id")
  var_ids    <- SeqArray::seqGetData(gds, "variant.id")
  chrom      <- SeqArray::seqGetData(gds, "chromosome")
  pos        <- SeqArray::seqGetData(gds, "position")

  # rsID if present, else CHR:POS
  rsid <- tryCatch(SeqArray::seqGetData(gds, "annotation/id"),
                   error = function(e) paste0(chrom, ":", pos))
  rsid[is.na(rsid) | rsid == "." | rsid == ""] <- paste0(chrom, ":", pos)[
    is.na(rsid) | rsid == "." | rsid == ""]

  # REF / ALT
  ref <- tryCatch(SeqArray::seqGetData(gds, "$ref"),
                  error = function(e) rep(NA_character_, length(var_ids)))
  alt <- tryCatch({
    a <- SeqArray::seqGetData(gds, "$alt")
    # alt may be a list for multi-allelic; take first
    if (is.list(a)) vapply(a, function(x) x[1L], character(1L)) else a
  }, error = function(e) rep(NA_character_, length(var_ids)))

  snp_info <- data.frame(SNP = rsid, CHR = chrom, POS = pos,
                         REF = ref, ALT = alt, stringsAsFactors = FALSE)

  list(
    type       = "gds",
    n_samples  = length(sample_ids),
    n_snps     = length(var_ids),
    sample_ids = as.character(sample_ids),
    snp_info   = snp_info,
    .gds       = gds,
    .var_ids   = var_ids
  )
}


# -- 5. PLINK BED -------------------------------------------------------------
.read_bed <- function(path, verbose) {
  if (isTRUE(verbose)) cat("[read_geno] Opening PLINK BED:", path, "\n")
  .require_pkg("BEDMatrix", "PLINK BED reader",
               install = "install.packages('BEDMatrix')")

  stem <- sub("\\.bed$", "", path, ignore.case = TRUE)
  bim  <- paste0(stem, ".bim")
  fam  <- paste0(stem, ".fam")
  if (!file.exists(bim)) stop("Cannot find .bim file: ", bim)
  if (!file.exists(fam)) stop("Cannot find .fam file: ", fam)

  bed <- BEDMatrix::BEDMatrix(path)

  bim_dt <- utils::read.table(bim, header = FALSE, stringsAsFactors = FALSE,
                              col.names = c("CHR","SNP","cM","POS","ALT","REF"))
  fam_dt <- utils::read.table(fam, header = FALSE, stringsAsFactors = FALSE,
                              col.names = c("FID","IID","PID","MID","SEX","PHENO"))

  snp_info   <- data.frame(SNP = bim_dt$SNP, CHR = bim_dt$CHR,
                           POS = bim_dt$POS, REF = bim_dt$REF,
                           ALT = bim_dt$ALT, stringsAsFactors = FALSE)
  sample_ids <- paste0(fam_dt$FID, "_", fam_dt$IID)

  list(
    type       = "bed",
    n_samples  = nrow(fam_dt),
    n_snps     = nrow(bim_dt),
    sample_ids = sample_ids,
    snp_info   = snp_info,
    .bed       = bed
  )
}


# -- 6. Plain R matrix (existing / backward-compatible) -----------------------
.wrap_matrix <- function(mat, snp_info, verbose) {
  if (!is.matrix(mat)) mat <- as.matrix(mat)
  if (is.null(snp_info))
    stop("snp_info must be supplied when format = 'matrix'.")
  req <- c("SNP", "CHR", "POS")
  if (!all(req %in% names(snp_info)))
    stop("snp_info must contain: ", paste(req, collapse = ", "))
  if (ncol(mat) != nrow(snp_info))
    stop("ncol(mat) [", ncol(mat), "] != nrow(snp_info) [", nrow(snp_info), "]")

  if (!"REF" %in% names(snp_info)) snp_info$REF <- NA_character_
  if (!"ALT" %in% names(snp_info)) snp_info$ALT <- NA_character_

  ids <- rownames(mat)
  if (is.null(ids)) ids <- paste0("ind", seq_len(nrow(mat)))

  .make_matrix_backend(mat, snp_info[, c("SNP","CHR","POS","REF","ALT")],
                       ids, "matrix")
}


# -- Internal: build in-memory backend list ------------------------------------
.make_matrix_backend <- function(mat, snp_info, sample_ids, type) {
  list(
    type       = type,
    n_samples  = nrow(mat),
    n_snps     = ncol(mat),
    sample_ids = as.character(sample_ids),
    snp_info   = snp_info,
    .mat       = mat
  )
}



# -- Internal: check / suggest package installation ----------------------------
.require_pkg <- function(pkg, context, install = NULL) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    msg <- sprintf("Package '%s' is required for %s.", pkg, context)
    if (!is.null(install)) msg <- paste0(msg, "\nInstall with: ", install)
    stop(msg)
  }
}


# -----------------------------------------------------------------------------
# Backend interface: read_chunk() and close_backend()
# -----------------------------------------------------------------------------

#' Extract a Genotype Slice from an LDxBlocks Backend
#'
#' @description
#' Returns an \code{n_samples x length(col_idx)} numeric matrix of dosage
#' values (0/1/2/NA) for the selected SNP columns. Works identically for all
#' supported backend types.
#'
#' @param backend An object of class \code{"LDxBlocks_backend"} as returned
#'   by \code{\link{read_geno}}.
#' @param col_idx Integer vector of column indices (1-based SNP positions).
#'
#' @return Numeric matrix (n_samples x length(col_idx)).
#'
#' @export
read_chunk <- function(backend, col_idx) {
  if (!inherits(backend, "LDxBlocks_backend"))
    stop("backend must be an LDxBlocks_backend object from read_geno().")

  switch(backend$type,
         numeric = ,
         hapmap  = ,
         vcf     = ,
         matrix  = backend$.mat[, col_idx, drop = FALSE],

         gds = {
           .require_pkg("SeqArray", "GDS chunk reader")
           var_ids <- backend$.var_ids[col_idx]
           if (!is.null(backend$.sample_filter)) {
             SeqArray::seqSetFilter(backend$.gds,
                                    sample.id  = backend$.sample_filter,
                                    variant.id = var_ids, verbose = FALSE)
           } else {
             SeqArray::seqSetFilter(backend$.gds,
                                    variant.id = var_ids, verbose = FALSE)
           }
           # "$dosage" convention varies across SeqArray versions (REF vs ALT count).
           # Compute ALT dosage directly from raw genotype [ploidy x n_samp x n_var]:
           #   allele encoding  0 = REF,  1 = ALT,  -1 = missing
           # apply over dims (2,3) = sample x variant, counting 1s = ALT alleles.
           gt  <- SeqArray::seqGetData(backend$.gds, "genotype")
           dos <- apply(gt, c(2L, 3L), function(a) {
             # SeqArray encodes missing as NA_integer_ (not -1).
             # anyNA() correctly handles both; sum(a==1L) counts ALT alleles.
             if (anyNA(a)) NA_real_ else sum(a == 1L)
           })
           storage.mode(dos) <- "numeric"
           SeqArray::seqResetFilter(backend$.gds, verbose = FALSE)
           if (!is.matrix(dos)) dos <- matrix(dos, nrow = backend$n_samples)
           rownames(dos) <- backend$sample_ids
           colnames(dos) <- backend$snp_info$SNP[col_idx]
           dos
         },

         bed = {
           .require_pkg("BEDMatrix", "PLINK BED chunk reader")
           idx_rows <- backend$.sample_filter %||% seq_len(backend$n_samples)
           m <- backend$.bed[idx_rows, col_idx, drop = FALSE]
           storage.mode(m) <- "numeric"
           m
         },

         stop("Unknown backend type: ", backend$type)
  )
}


#' Close an LDxBlocks Backend and Release File Handles
#'
#' @description
#' Closes any open file connections held by the backend. For in-memory
#' backends (\code{"matrix"}, \code{"numeric"}, \code{"hapmap"}, \code{"vcf"})
#' this is a no-op. For \code{"gds"} backends it calls
#' \code{SeqArray::seqClose()}. For \code{"bed"} backends the memory-mapped
#' file is released.
#'
#' @param backend An \code{"LDxBlocks_backend"} object.
#'
#' @return Invisibly \code{NULL}.
#'
#' @export
close_backend <- function(backend) {
  if (!inherits(backend, "LDxBlocks_backend")) return(invisible(NULL))
  if (backend$type == "gds" && !is.null(backend$.gds)) {
    try(SeqArray::seqClose(backend$.gds), silent = TRUE)
  }
  invisible(NULL)
}


#' Print Method for LDxBlocks Backend
#' @param x An \code{LDxBlocks_backend} object.
#' @param ... Further arguments passed to or from other methods (unused).
#' @export
print.LDxBlocks_backend <- function(x, ...) {
  cat("LDxBlocks backend\n")
  cat("  Type      :", x$type, "\n")
  cat("  Samples   :", x$n_samples, "\n")
  cat("  SNPs      :", x$n_snps, "\n")
  cat("  Chr       :", paste(unique(x$snp_info$CHR), collapse = ", "), "\n")
  invisible(x)
}


#' Summary Method for LDxBlocks Backend
#' @param object An \code{LDxBlocks_backend} object.
#' @param ... Further arguments passed to or from other methods (unused).
#' @export
summary.LDxBlocks_backend <- function(object, ...) {
  chrs <- table(.norm_chr(object$snp_info$CHR))
  cat("LDxBlocks backend summary\n")
  cat("  Format  :", object$type, "\n")
  cat("  Samples :", object$n_samples, "\n")
  cat("  SNPs    :", object$n_snps, "\n")
  cat("  Chromosomes:\n")
  for (nm in names(chrs))
    cat(sprintf("    chr%-4s  %d SNPs\n", nm, chrs[[nm]]))
  invisible(object)
}


# -- Internal: HapMap -> GDS converter -----------------------------------------
# SeqArray has no native HapMap reader. This function streams through a
# HapMap file row-by-row via data.table, decodes two-character nucleotide
# calls (AA, AT, TT, NN) to 0/1/2/NA dosage, and writes a SeqArray GDS
# using gdsfmt's low-level node API.
#
# Format assumption (standard HapMap):
#   Column 1  : rs#        (SNP ID)
#   Column 2  : alleles    "REF/ALT"
#   Column 3  : chrom      (chromosome)
#   Column 4  : pos        (position)
#   Columns 12+: sample columns with two-character nucleotide calls
#
.hapmap_to_gds <- function(hmp_path, gds_path, na_strings, verbose) {
  # Strategy: read HapMap in full with data.table, convert nucleotide calls to
  # 0/1/2, write a minimal VCF to a temp file, then use seqVCF2GDS() which
  # IS a real SeqArray export. This avoids SeqArray::seqNewGDS() which does
  # not exist in the package.
  .require_pkg("data.table", "HapMap-to-GDS converter")
  .require_pkg("SeqArray",   "HapMap-to-GDS converter")
  .require_pkg("gdsfmt",     "HapMap-to-GDS converter")

  if (isTRUE(verbose)) message("[.hapmap_to_gds] Reading HapMap: ", basename(hmp_path))
  dt <- data.table::fread(hmp_path, na.strings = na_strings,
                          showProgress = isTRUE(verbose), data.table = TRUE)

  # Normalise column names
  hdr_map <- c("rs#"="SNP","chrom"="CHR","pos"="POS","alleles"="alleles")
  nms_low <- tolower(names(dt))
  for (k in names(hdr_map)) {
    idx <- which(nms_low == k)
    if (length(idx)) data.table::setnames(dt, names(dt)[idx], hdr_map[k])
  }
  if (!all(c("SNP","CHR","POS","alleles") %in% names(dt)))
    stop("HapMap file missing rs#, chrom, pos, or alleles columns.", call.=FALSE)

  parts <- strsplit(as.character(dt$alleles), "/", fixed=TRUE)
  ref_v  <- vapply(parts, `[`, character(1L), 1L)
  alt_v  <- vapply(parts, function(x) if (length(x)>1L) x[2L] else NA_character_,
                   character(1L))

  meta_cols <- c("SNP","CHR","POS","alleles","strand","assembly#","center",
                 "protLSID","assayLSID","panelLSID","QCcode")
  meta_idx  <- which(tolower(names(dt)) %in% tolower(meta_cols))
  samp_cols <- names(dt)[-meta_idx]
  if (!length(samp_cols)) stop("No sample columns found in HapMap file.", call.=FALSE)

  n_snp  <- nrow(dt)
  n_samp <- length(samp_cols)
  if (isTRUE(verbose))
    message("[.hapmap_to_gds] ", n_snp, " SNPs x ", n_samp, " samples - decoding calls")

  # Decode nucleotide calls to GT strings (0/0, 0/1, 1/1, ./.)
  gt_mat <- matrix(NA_character_, nrow=n_snp, ncol=n_samp)
  for (i in seq_len(n_snp)) {
    r <- ref_v[i]; a <- alt_v[i]
    if (is.na(r) || is.na(a)) { gt_mat[i,] <- "./."; next }
    calls <- as.character(dt[i, ..samp_cols, drop=TRUE])
    rr <- paste0(r,r); ra1 <- paste0(r,a); ra2 <- paste0(a,r); aa <- paste0(a,a)
    gt <- rep("./.", n_samp)
    gt[calls==rr]  <- "0/0"
    gt[calls==ra1 | calls==ra2] <- "0/1"
    gt[calls==aa]  <- "1/1"
    gt_mat[i,] <- gt
  }

  # Write minimal VCF to temp file
  vcf_tmp <- tempfile(fileext=".vcf")
  on.exit(unlink(vcf_tmp), add=TRUE)

  if (isTRUE(verbose)) message("[.hapmap_to_gds] Writing temp VCF: ", basename(vcf_tmp))
  header <- c(
    "##fileformat=VCFv4.2",
    "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">",
    paste(c("#CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT",
            samp_cols), collapse="	")
  )
  data_rows <- paste(
    as.character(dt$CHR), as.integer(dt$POS),
    as.character(dt$SNP), ref_v, alt_v,
    ".", "PASS", ".", "GT",
    apply(gt_mat, 1L, paste, collapse="	"),
    sep="	"
  )
  writeLines(c(header, data_rows), vcf_tmp)

  if (isTRUE(verbose)) message("[.hapmap_to_gds] Converting VCF -> GDS: ", basename(gds_path))
  SeqArray::seqVCF2GDS(
    vcf.fn         = vcf_tmp,
    out.fn         = gds_path,
    storage.option = "ZIP_RA",
    verbose        = isTRUE(verbose)
  )
  if (isTRUE(verbose)) message("[.hapmap_to_gds] Done.")
  invisible(gds_path)
}
