## tests/testthat/test-io.R
## ─────────────────────────────────────────────────────────────────────────────
## Tests for the read_geno() / read_chunk() / close_backend() I/O layer.
## Each format is tested by writing a minimal file to a tempfile, reading it
## with read_geno(), and verifying the backend's properties and read_chunk()
## output — independent of inst/extdata so they pass even without built data.
## ─────────────────────────────────────────────────────────────────────────────

library(testthat)
library(LDxBlocks)

# ── Local helper: tiny synthetic genotype matrix ──────────────────────────────
# Uses smaller defaults than make_geno() from helper.R (12 x 8 vs 60 x 30)
# to keep I/O tests fast. Local because test-io.R should be self-contained.
make_tiny <- function(n = 12, p = 8, seed = 7L) {
  set.seed(seed)
  G <- matrix(sample(0:2, n * p, replace = TRUE), n, p)
  rownames(G) <- paste0("s", seq_len(n))
  colnames(G) <- paste0("rs", seq_len(p))
  G
}

# ── Local flat-file writers ───────────────────────────────────────────────────
# These are intentionally local so test-io.R is self-contained and does not
# depend on helper.R write_* functions, which use write.table with sep="\t".

write_numeric_csv <- function(G, info, path) {
  if (!"REF" %in% names(info)) info$REF <- "A"
  if (!"ALT" %in% names(info)) info$ALT <- "T"
  meta <- info[, c("SNP","CHR","POS","REF","ALT")]
  df   <- cbind(meta, as.data.frame(t(G)))
  write.csv(df, path, row.names = FALSE, quote = FALSE)
}

write_hapmap <- function(G, info, path) {
  hmp_decode_snp <- function(snp_dos, ref, alt) {
    vapply(snp_dos, function(x) {
      if (is.na(x)) return("NN")
      switch(as.character(as.integer(x)),
             "0" = paste0(ref, ref),
             "1" = paste0(ref, alt),
             "2" = paste0(alt, alt),
             "NN")
    }, character(1))
  }
  n_snp <- nrow(info); n_ind <- nrow(G)
  calls <- matrix("NN", nrow = n_snp, ncol = n_ind)
  for (si in seq_len(n_snp))
    calls[si, ] <- hmp_decode_snp(G[, si], info$REF[si], info$ALT[si])
  hdr <- data.frame(
    "rs#"=info$SNP, alleles=paste0(info$REF,"/",info$ALT),
    chrom=info$CHR, pos=info$POS,
    strand="+", "assembly#"="NA", center="NA", protLSID="NA",
    assayLSID="NA", panelLSID="NA", QCcode="NA",
    check.names=FALSE, stringsAsFactors=FALSE
  )
  out <- cbind(hdr, as.data.frame(calls, stringsAsFactors=FALSE))
  colnames(out)[12:ncol(out)] <- rownames(G)
  write.table(out, path, sep="\t", row.names=FALSE, quote=FALSE)
}

write_vcf <- function(G, info, path) {
  gt_enc <- function(g) vapply(g, function(x) {
    if (is.na(x)) return("./.")
    c("0"="0/0","1"="0/1","2"="1/1")[as.character(x)]
  }, character(1))
  lines <- c(
    "##fileformat=VCFv4.2",
    "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">",  # required by SNPRelate
    paste(c("#CHROM","POS","ID","REF","ALT","QUAL","FILTER",
            "INFO","FORMAT",rownames(G)), collapse="\t")
  )
  if (!"REF" %in% names(info)) info$REF <- "A"
  if (!"ALT" %in% names(info)) info$ALT <- "T"
  for (i in seq_len(nrow(info))) {
    lines <- c(lines, paste(c(
      info$CHR[i], info$POS[i], info$SNP[i],
      info$REF[i], info$ALT[i], ".", "PASS", ".", "GT",
      gt_enc(G[, i])
    ), collapse="\t"))
  }
  writeLines(lines, path)
}

# ── Format: matrix ────────────────────────────────────────────────────────────

test_that("matrix backend: n_samples, n_snps, type are correct", {
  G    <- make_tiny()
  info <- data.frame(SNP=colnames(G), CHR="1",
                     POS=seq(1000L, by=2000L, length.out=8L))
  be   <- read_geno(G, format="matrix", snp_info=info)
  expect_equal(be$type,       "matrix")
  expect_equal(be$n_samples,  12L)
  expect_equal(be$n_snps,     8L)
  expect_equal(be$sample_ids, paste0("s", 1:12))
  close_backend(be)
})

test_that("matrix backend: read_chunk returns correct values", {
  G    <- make_tiny()
  info <- data.frame(SNP=colnames(G), CHR="1",
                     POS=seq(1000L, by=2000L, length.out=8L))
  be    <- read_geno(G, format="matrix", snp_info=info)
  chunk <- read_chunk(be, 2:5)
  expect_equal(dim(chunk), c(12L, 4L))
  expect_equal(chunk, G[, 2:5])
  close_backend(be)
})

test_that("matrix backend: snp_info without REF/ALT gets NA columns added", {
  G    <- make_tiny()
  info <- data.frame(SNP=colnames(G), CHR="2",
                     POS=seq(1000L, by=1000L, length.out=8L))
  be   <- read_geno(G, format="matrix", snp_info=info)
  expect_true("REF" %in% names(be$snp_info))
  expect_true("ALT" %in% names(be$snp_info))
  expect_true(all(is.na(be$snp_info$REF)))
  close_backend(be)
})

test_that("matrix backend: chromosome normalisation strips 'chr' prefix", {
  G    <- make_tiny()
  info <- data.frame(SNP=colnames(G), CHR="chr3",
                     POS=seq(1000L, by=1000L, length.out=8L))
  be   <- read_geno(G, format="matrix", snp_info=info)
  expect_equal(unique(be$snp_info$CHR), "3")
  close_backend(be)
})

test_that("matrix backend: missing snp_info throws informative error", {
  expect_error(read_geno(make_tiny(), format="matrix"), "snp_info must be supplied")
})

test_that("matrix backend: dimension mismatch throws error", {
  G    <- make_tiny()
  info <- data.frame(SNP=paste0("rs",1:5), CHR="1", POS=1:5)  # 5 rows vs 8 cols
  expect_error(read_geno(G, format="matrix", snp_info=info), "ncol\\(mat\\)")
})

# ── Format: numeric dosage CSV ────────────────────────────────────────────────

test_that("numeric CSV: basic read properties", {
  G    <- make_tiny()
  info <- data.frame(SNP=colnames(G), CHR="1",
                     POS=seq(1000L,by=2000L,length.out=8L),
                     REF="A", ALT="T", stringsAsFactors=FALSE)
  tmp  <- tempfile(fileext=".csv")
  write_numeric_csv(G, info, tmp)
  be <- read_geno(tmp)
  expect_equal(be$type,      "numeric")
  expect_equal(be$n_snps,    8L)
  expect_equal(be$n_samples, 12L)
  expect_equal(sort(be$sample_ids), sort(paste0("s", 1:12)))
  close_backend(be); unlink(tmp)
})

test_that("numeric CSV: read_chunk values match original matrix", {
  G    <- make_tiny()
  info <- data.frame(SNP=colnames(G), CHR="1",
                     POS=seq(1000L,by=2000L,length.out=8L),
                     REF="A", ALT="T", stringsAsFactors=FALSE)
  tmp  <- tempfile(fileext=".csv")
  write_numeric_csv(G, info, tmp)
  be      <- read_geno(tmp)
  chunk   <- read_chunk(be, 1:4)
  row_ord <- match(be$sample_ids, rownames(G))
  expect_equal(chunk, G[row_ord, 1:4], ignore_attr=TRUE)
  close_backend(be); unlink(tmp)
})

test_that("numeric CSV: NA values survive round-trip", {
  G       <- make_tiny()
  G[2, 3] <- NA
  info <- data.frame(SNP=colnames(G), CHR="1",
                     POS=seq(1000L,by=2000L,length.out=8L),
                     REF="A", ALT="T", stringsAsFactors=FALSE)
  tmp  <- tempfile(fileext=".csv")
  write_numeric_csv(G, info, tmp)
  be    <- read_geno(tmp)
  chunk <- read_chunk(be, 3L)
  expect_true(any(is.na(chunk)))
  close_backend(be); unlink(tmp)
})

# ── Format: HapMap ────────────────────────────────────────────────────────────

test_that("HapMap: dosage decoding is correct for all three calls and NA", {
  G    <- matrix(c(0L,1L,2L,NA_integer_), nrow=1, ncol=4)
  rownames(G) <- "s1"; colnames(G) <- paste0("rs",1:4)
  info <- data.frame(SNP=colnames(G), CHR="1", POS=1:4*1000L,
                     REF=c("A","G","C","T"), ALT=c("T","C","G","A"),
                     stringsAsFactors=FALSE)
  tmp <- tempfile(fileext=".hmp.txt")
  write_hapmap(G, info, tmp)
  be    <- read_geno(tmp)
  chunk <- read_chunk(be, 1:4)
  expect_equal(as.integer(chunk["s1", 1]), 0L)
  expect_equal(as.integer(chunk["s1", 2]), 1L)
  expect_equal(as.integer(chunk["s1", 3]), 2L)
  expect_true(is.na(chunk["s1", 4]))
  close_backend(be); unlink(tmp)
})

test_that("HapMap: backend properties match input", {
  G    <- make_tiny()
  info <- data.frame(SNP=colnames(G), CHR="2",
                     POS=seq(1000L,by=3000L,length.out=8L),
                     REF="G", ALT="C", stringsAsFactors=FALSE)
  tmp <- tempfile(fileext=".hmp.txt")
  write_hapmap(G, info, tmp)
  be <- read_geno(tmp)
  expect_true(be$type %in% c("hapmap","gds"))
  expect_equal(be$n_snps,    8L)
  expect_equal(be$n_samples, 12L)
  close_backend(be); unlink(tmp)
})

# ── Format: VCF ──────────────────────────────────────────────────────────────

test_that("VCF: unphased GT decoded correctly for all three dosages", {
  G    <- matrix(c(0L,1L,2L), nrow=1, ncol=3)
  rownames(G) <- "s1"; colnames(G) <- paste0("rs",1:3)
  info <- data.frame(SNP=colnames(G), CHR="1", POS=c(1000L,2000L,3000L),
                     REF="A", ALT="T", stringsAsFactors=FALSE)
  tmp <- tempfile(fileext=".vcf")
  write_vcf(G, info, tmp)
  be    <- read_geno(tmp)
  chunk <- read_chunk(be, 1:3)
  expect_equal(as.integer(chunk["s1", 1]), 0L)
  expect_equal(as.integer(chunk["s1", 2]), 1L)
  expect_equal(as.integer(chunk["s1", 3]), 2L)
  close_backend(be); unlink(tmp)
})

test_that("VCF: missing ./. becomes NA", {
  G    <- matrix(NA_integer_, nrow=1, ncol=2)
  rownames(G) <- "s1"; colnames(G) <- paste0("rs",1:2)
  info <- data.frame(SNP=colnames(G), CHR="1", POS=c(1000L,2000L),
                     REF="A", ALT="T", stringsAsFactors=FALSE)
  tmp <- tempfile(fileext=".vcf")
  write_vcf(G, info, tmp)
  be    <- read_geno(tmp)
  chunk <- read_chunk(be, 1:2)
  expect_true(all(is.na(chunk["s1", ])))
  close_backend(be); unlink(tmp)
})

test_that("VCF: chr prefix stripped from chromosome labels", {
  G    <- make_tiny(4L, 4L)
  info <- data.frame(SNP=colnames(G), CHR="chr5",
                     POS=seq(1000L,by=1000L,length.out=4L),
                     REF="A", ALT="T", stringsAsFactors=FALSE)
  tmp <- tempfile(fileext=".vcf")
  write_vcf(G, info, tmp)
  be  <- read_geno(tmp)
  expect_equal(unique(be$snp_info$CHR), "5")
  close_backend(be); unlink(tmp)
})

test_that("VCF: backend properties correct", {
  G    <- make_tiny()
  info <- data.frame(SNP=colnames(G), CHR="1",
                     POS=seq(1000L,by=2000L,length.out=8L),
                     REF="A", ALT="G", stringsAsFactors=FALSE)
  tmp <- tempfile(fileext=".vcf")
  write_vcf(G, info, tmp)
  be <- read_geno(tmp)
  expect_true(be$type %in% c("vcf","gds"))
  expect_equal(be$n_snps,    8L)
  expect_equal(be$n_samples, 12L)
  close_backend(be); unlink(tmp)
})

test_that("VCF: clean_malformed=TRUE removes lines with wrong field count", {
  G    <- make_tiny(4L, 3L)
  info <- data.frame(SNP=colnames(G), CHR="1",
                     POS=c(1000L,2000L,3000L),
                     REF="A", ALT="T", stringsAsFactors=FALSE)
  tmp <- tempfile(fileext=".vcf")
  write_vcf(G, info, tmp)
  # Append a malformed line (too few fields)
  cat("\n1\t99999\tBAD_SNP\tA\tT\t.\tPASS\n", file=tmp, append=TRUE)
  # clean_malformed=TRUE should silently drop the bad line and still parse
  expect_no_error({
    be <- read_geno(tmp, clean_malformed=TRUE)
    close_backend(be)
  })
  unlink(tmp)
})

# ── GDS backend ───────────────────────────────────────────────────────────────

test_that("GDS backend: read from SNPRelate GDS file", {
  skip_if_not_installed("SNPRelate")
  G    <- make_tiny()
  info <- data.frame(SNP=colnames(G), CHR="1",
                     POS=seq(1000L,by=2000L,length.out=8L),
                     REF="A", ALT="T", stringsAsFactors=FALSE)
  # Write VCF, convert to GDS, then read directly as GDS
  tmp_vcf <- tempfile(fileext=".vcf")
  tmp_gds <- tempfile(fileext=".gds")
  write_vcf(G, info, tmp_vcf)
  SNPRelate::snpgdsVCF2GDS(tmp_vcf, tmp_gds,
                           method="biallelic.only",
                           snpfirstdim=FALSE, verbose=FALSE)
  be <- read_geno(tmp_gds)
  expect_equal(be$type,      "gds")
  expect_equal(be$n_snps,    8L)
  expect_equal(be$n_samples, 12L)
  chunk <- read_chunk(be, 1:4)
  expect_equal(dim(chunk), c(12L, 4L))
  close_backend(be)
  unlink(c(tmp_vcf, tmp_gds))
})

# ── PLINK BED backend ─────────────────────────────────────────────────────────

test_that("BED backend: read from PLINK BED/BIM/FAM files", {
  skip_if_not_installed("BEDMatrix")
  G    <- make_tiny()
  info <- data.frame(SNP=colnames(G), CHR="1",
                     POS=seq(1000L,by=2000L,length.out=8L),
                     REF="A", ALT="T", stringsAsFactors=FALSE)
  # Write VCF then convert to PLINK BED via SNPRelate (requires SNPRelate)
  skip_if_not_installed("SNPRelate")
  tmp_vcf <- tempfile(fileext=".vcf")
  tmp_gds <- tempfile(fileext=".gds")
  tmp_bed <- tempfile()  # prefix only; .bed/.bim/.fam appended by PLINK tools
  write_vcf(G, info, tmp_vcf)
  SNPRelate::snpgdsVCF2GDS(tmp_vcf, tmp_gds,
                           method="biallelic.only",
                           snpfirstdim=FALSE, verbose=FALSE)
  SNPRelate::snpgdsGDS2BED(tmp_gds, tmp_bed, verbose=FALSE)
  be <- read_geno(paste0(tmp_bed, ".bed"))
  expect_equal(be$type,      "bed")
  expect_equal(be$n_snps,    8L)
  expect_equal(be$n_samples, 12L)
  chunk <- read_chunk(be, 1:3)
  expect_equal(dim(chunk), c(12L, 3L))
  close_backend(be)
  unlink(c(tmp_vcf, tmp_gds,
           paste0(tmp_bed, c(".bed",".bim",".fam"))))
})

# ── Auto-detection ────────────────────────────────────────────────────────────

test_that("format auto-detection works for CSV, HapMap, and VCF extensions", {
  G    <- make_tiny(6L, 4L)
  info <- data.frame(SNP=colnames(G), CHR="1",
                     POS=seq(1000L,by=1000L,length.out=4L),
                     REF="A", ALT="T", stringsAsFactors=FALSE)

  tmp_csv <- tempfile(fileext=".csv")
  write_numeric_csv(G, info, tmp_csv)
  be_csv <- read_geno(tmp_csv)
  expect_equal(be_csv$type, "numeric")
  close_backend(be_csv); unlink(tmp_csv)

  tmp_hmp <- tempfile(fileext=".hmp.txt")
  write_hapmap(G, info, tmp_hmp)
  be_hmp <- read_geno(tmp_hmp)
  expect_true(be_hmp$type %in% c("hapmap","gds"))
  close_backend(be_hmp); unlink(tmp_hmp)

  tmp_vcf <- tempfile(fileext=".vcf")
  write_vcf(G, info, tmp_vcf)
  be_vcf <- read_geno(tmp_vcf)
  expect_true(be_vcf$type %in% c("vcf","gds"))
  close_backend(be_vcf); unlink(tmp_vcf)
})

test_that("unknown extension without format= throws informative error", {
  tmp <- tempfile(fileext=".xyz")
  writeLines("dummy", tmp)
  expect_error(read_geno(tmp), "Cannot detect genotype format")
  unlink(tmp)
})

# ── print / summary S3 methods ────────────────────────────────────────────────

test_that("print.LDxBlocks_backend outputs type and dimension", {
  G    <- make_tiny()
  info <- data.frame(SNP=colnames(G), CHR="1",
                     POS=seq(1000L,by=2000L,length.out=8L))
  be   <- read_geno(G, format="matrix", snp_info=info)
  expect_output(print(be), "LDxBlocks backend")
  expect_output(print(be), "matrix")
  close_backend(be)
})

test_that("summary.LDxBlocks_backend shows per-chromosome SNP counts", {
  G    <- make_tiny()
  info <- data.frame(SNP=colnames(G), CHR=c(rep("1",4L),rep("2",4L)),
                     POS=seq(1000L,by=2000L,length.out=8L))
  be   <- read_geno(G, format="matrix", snp_info=info)
  expect_output(summary(be), "Chromosomes")
  close_backend(be)
})

# ── bigmemory backend ─────────────────────────────────────────────────────────

test_that("read_geno_bigmemory: accepts file path directly (no snp_info needed)", {
  skip_if_not_installed("bigmemory")
  G    <- make_tiny()
  info <- data.frame(SNP=colnames(G), CHR="1",
                     POS=seq(1000L,by=2000L,length.out=8L),
                     REF="A", ALT="T", stringsAsFactors=FALSE)
  tmp_csv <- tempfile(fileext=".csv")
  write_numeric_csv(G, info, tmp_csv)
  on.exit(unlink(tmp_csv))

  # File path route: read_geno_bigmemory() opens the CSV internally
  be_bm <- read_geno_bigmemory(
    source      = tmp_csv,
    backingfile = tempfile("ldxbm_io_"),
    backingpath = tempdir(),
    type        = "char",
    verbose     = FALSE
  )
  on.exit(close_backend(be_bm), add = TRUE)

  expect_s3_class(be_bm, "LDxBlocks_backend")
  expect_equal(be_bm$type, "bigmemory")
  expect_equal(be_bm$n_snps,    8L)
  expect_equal(be_bm$n_samples, 12L)
})

test_that("read_geno_bigmemory: read_chunk values match original (type='char')", {
  skip_if_not_installed("bigmemory")
  G    <- make_tiny()
  info <- data.frame(SNP=colnames(G), CHR="1",
                     POS=seq(1000L,by=2000L,length.out=8L),
                     REF="A", ALT="T", stringsAsFactors=FALSE)
  be_mat <- read_geno(G, format="matrix", snp_info=info)
  be_bm  <- read_geno_bigmemory(
    source      = be_mat,
    backingfile = tempfile("ldxbm_vals_"),
    backingpath = tempdir(),
    type        = "char",
    verbose     = FALSE
  )
  on.exit({ close_backend(be_bm); close_backend(be_mat) })

  chunk_mat <- read_chunk(be_mat, 1:8)
  chunk_bm  <- read_chunk(be_bm,  1:8)
  # Values should match (integer matrix from char big.matrix)
  expect_equal(chunk_bm[order(rownames(chunk_bm)), ],
               chunk_mat[order(rownames(chunk_mat)), ],
               ignore_attr = TRUE)
})

test_that("read_geno_bigmemory: NA genotypes correctly restored from char -128", {
  skip_if_not_installed("bigmemory")
  G       <- make_tiny()
  G[3, 5] <- NA_integer_   # inject NA at individual 3, SNP 5
  info <- data.frame(SNP=colnames(G), CHR="1",
                     POS=seq(1000L,by=2000L,length.out=8L),
                     REF="A", ALT="T", stringsAsFactors=FALSE)
  be_mat <- read_geno(G, format="matrix", snp_info=info)
  be_bm  <- read_geno_bigmemory(
    source      = be_mat,
    backingfile = tempfile("ldxbm_na_"),
    backingpath = tempdir(),
    type        = "char",
    verbose     = FALSE
  )
  on.exit({ close_backend(be_bm); close_backend(be_mat) })

  chunk <- read_chunk(be_bm, 5L)
  # Row for individual 3 should be NA
  ind3_row <- which(rownames(chunk) == "s3")
  expect_true(is.na(chunk[ind3_row, 1]))
  # Other individuals at SNP 5 should be non-NA (0, 1, or 2)
  expect_false(any(is.na(chunk[-ind3_row, ])))
})

test_that("read_geno_bigmemory: type='double' produces correct numeric values", {
  skip_if_not_installed("bigmemory")
  G    <- make_tiny(6L, 4L)
  info <- data.frame(SNP=colnames(G), CHR="1",
                     POS=seq(1000L,by=1000L,length.out=4L),
                     REF="A", ALT="T", stringsAsFactors=FALSE)
  be_mat <- read_geno(G, format="matrix", snp_info=info)
  be_bm  <- read_geno_bigmemory(
    source      = be_mat,
    backingfile = tempfile("ldxbm_dbl_"),
    backingpath = tempdir(),
    type        = "double",
    verbose     = FALSE
  )
  on.exit({ close_backend(be_bm); close_backend(be_mat) })

  chunk <- read_chunk(be_bm, 1:4)
  expect_equal(dim(chunk), c(6L, 4L))
  # All values should be 0, 1, or 2 (no spurious -128 or other artefacts)
  vals <- as.vector(chunk[!is.na(chunk)])
  expect_true(all(vals %in% c(0, 1, 2)))
})

test_that("read_geno_bigmemory: reattach from .desc file works", {
  skip_if_not_installed("bigmemory")
  G    <- make_tiny(6L, 4L)
  info <- data.frame(SNP=colnames(G), CHR="1",
                     POS=seq(1000L,by=1000L,length.out=4L),
                     REF="A", ALT="T", stringsAsFactors=FALSE)
  be_mat  <- read_geno(G, format="matrix", snp_info=info)
  bm_stem <- tempfile("ldxbm_reattach_")
  be_bm1  <- read_geno_bigmemory(
    source      = be_mat,
    backingfile = basename(bm_stem),
    backingpath = dirname(bm_stem),
    type        = "char",
    verbose     = FALSE
  )
  close_backend(be_bm1)
  close_backend(be_mat)

  # Reattach using the .desc file path
  desc_path <- paste0(bm_stem, ".desc")
  be_bm2 <- read_geno_bigmemory(
    source      = desc_path,
    snp_info    = info,            # required on reattach
    backingfile = basename(bm_stem),
    backingpath = dirname(bm_stem),
    verbose     = FALSE
  )
  on.exit(close_backend(be_bm2))

  expect_s3_class(be_bm2, "LDxBlocks_backend")
  expect_equal(be_bm2$n_snps, 4L)
  chunk <- read_chunk(be_bm2, 1:4)
  expect_equal(dim(chunk), c(6L, 4L))
})
