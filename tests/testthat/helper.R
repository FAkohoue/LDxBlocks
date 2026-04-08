## tests/testthat/helper.R
## ─────────────────────────────────────────────────────────────────────────────
## Shared helpers loaded automatically by testthat before every test file.
## Functions here are available in all test-*.R files without explicit source().
## ─────────────────────────────────────────────────────────────────────────────

## Generate a small random genotype matrix with rownames and colnames.
## Used in test-cpp.R and test-io.R to create fixtures without relying on
## inst/extdata (which may not exist in source-only test runs).
make_geno <- function(n = 60, p = 30, seed = 1L) {
  set.seed(seed)
  G <- matrix(sample(0:2, n * p, replace = TRUE), nrow = n, ncol = p)
  rownames(G) <- paste0("ind", seq_len(n))
  colnames(G) <- paste0("rs",  seq_len(p))
  G
}

## Minimal SNP info data.frame matching a matrix produced by make_geno().
make_snpinfo <- function(p = 30, chr = "1") {
  data.frame(
    SNP = paste0("rs",  seq_len(p)),
    CHR = chr,
    POS = seq(1000L, by = 5000L, length.out = p),
    stringsAsFactors = FALSE
  )
}

## Write a minimal numeric dosage CSV to a tempfile and return the path.
## Caller is responsible for unlink(path) after use.
write_numeric_csv <- function(G, info, sep = ",") {
  if (!"REF" %in% names(info)) info$REF <- "A"
  if (!"ALT" %in% names(info)) info$ALT <- "T"
  meta <- info[, c("SNP","CHR","POS","REF","ALT")]
  df   <- cbind(meta, as.data.frame(t(G)))
  path <- tempfile(fileext = ".csv")
  write.table(df, path, sep = sep, row.names = FALSE, quote = FALSE)
  path
}

## Write a minimal VCF to a tempfile and return the path.
write_vcf <- function(G, info) {
  gt_enc <- function(g) vapply(g, function(x) {
    if (is.na(x)) return("./.")
    c("0" = "0/0", "1" = "0/1", "2" = "1/1")[as.character(x)]
  }, character(1))

  path  <- tempfile(fileext = ".vcf")
  lines <- c(
    "##fileformat=VCFv4.2",
    paste(c("#CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT",
            rownames(G)), collapse = "\t")
  )
  ref <- if ("REF" %in% names(info)) info$REF else rep("A", nrow(info))
  alt <- if ("ALT" %in% names(info)) info$ALT else rep("T", nrow(info))

  for (i in seq_len(nrow(info))) {
    lines <- c(lines, paste(c(
      info$CHR[i], info$POS[i], info$SNP[i],
      ref[i], alt[i], ".", "PASS", ".", "GT",
      gt_enc(G[, i])
    ), collapse = "\t"))
  }
  writeLines(lines, path)
  path
}

## Write a minimal HapMap file to a tempfile and return the path.
write_hapmap <- function(G, info) {
  ref <- if ("REF" %in% names(info)) info$REF else rep("A", nrow(info))
  alt <- if ("ALT" %in% names(info)) info$ALT else rep("T", nrow(info))

  hmp_decode <- function(g, r, a) {
    vapply(g, function(x) {
      if (is.na(x)) return("NN")
      c("0" = paste0(r, r), "1" = paste0(r, a), "2" = paste0(a, a))[as.character(x)]
    }, character(1))
  }

  sample_calls <- mapply(
    function(r, a, col) hmp_decode(col, r, a),
    ref, alt, as.list(as.data.frame(t(G))),
    SIMPLIFY = FALSE
  )
  sample_mat <- t(do.call(cbind, sample_calls))   # individuals x SNPs

  hdr <- data.frame(
    "rs#"       = info$SNP,
    alleles     = paste0(ref, "/", alt),
    chrom       = info$CHR,
    pos         = info$POS,
    strand      = "+",
    "assembly#" = "NA",
    center      = "NA",
    protLSID    = "NA",
    assayLSID   = "NA",
    panelLSID   = "NA",
    QCcode      = "NA",
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  out <- cbind(hdr, sample_mat)
  colnames(out)[12:ncol(out)] <- rownames(G)

  path <- tempfile(fileext = ".hmp.txt")
  write.table(out, path, sep = "\t", row.names = FALSE, quote = FALSE)
  path
}
