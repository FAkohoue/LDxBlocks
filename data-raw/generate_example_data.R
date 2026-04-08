## data-raw/generate_example_data.R
## ─────────────────────────────────────────────────────────────────────────────
## Generates all example datasets shipped with LDxBlocks.
## Run once with: source("data-raw/generate_example_data.R")
## Outputs land in data/ (as .rda) and inst/extdata/ (as flat files).
##
## Datasets produced
## ─────────────────
##   ldx_geno      numeric matrix 120 ind × 230 SNPs  (three-block structure per chromosome)
##   ldx_snp_info  data.frame SNP / CHR / POS / REF / ALT
##   ldx_blocks    data.frame  — pre-computed block table for examples/tests
##   ldx_gwas      data.frame  — 20 toy GWAS markers for tune_LD_params demos
##
## Flat-file copies
## ────────────────
##   inst/extdata/example_genotypes_numeric.csv   (numeric dosage)
##   inst/extdata/example_genotypes.hmp.txt       (HapMap)
##   inst/extdata/example_genotypes.vcf           (VCF)
## ─────────────────────────────────────────────────────────────────────────────

set.seed(20250407)

# ── Parameters ────────────────────────────────────────────────────────────────
N_IND   <- 120      # individuals
N_CHR   <- 3        # chromosomes
BLOCK_SIZES <- list(
  chr1 = c(25, 20, 25),   # three blocks on chr1 (sums to 70)
  chr2 = c(30, 20, 20),   # three blocks on chr2
  chr3 = c(20, 20, 20)    # three blocks on chr3
)
INTER_BLOCK_SNP <- 5   # low-LD SNPs between blocks (singleton stretches)

ALLELES <- matrix(c(
  "A","T",  "G","C",  "C","G",  "T","A",  "A","G",
  "G","A",  "C","T",  "T","C",  "A","C",  "G","T"
), ncol = 2, byrow = TRUE)

# ── Helper: simulate one LD block ─────────────────────────────────────────────
make_block <- function(n_ind, n_snp, p, flip = 0.04) {
  founder <- rbinom(n_ind, 2, p)
  M <- matrix(founder, nrow = n_ind, ncol = n_snp)
  for (j in seq_len(n_snp)) {
    idx <- sample.int(n_ind, max(1L, floor(flip * n_ind)))
    delta <- sample(c(-1L, 1L), length(idx), replace = TRUE)
    M[idx, j] <- pmin(2L, pmax(0L, M[idx, j] + delta))
  }
  M
}

# ── Helper: simulate low-LD singleton SNPs ────────────────────────────────────
make_singleton <- function(n_ind, n_snp) {
  sapply(seq_len(n_snp), function(i) {
    p <- runif(1, 0.1, 0.4)
    rbinom(n_ind, 2, p)
  })
}

# ── Build genome ──────────────────────────────────────────────────────────────
geno_list    <- list()
snpinfo_list <- list()
global_pos   <- 1L

for (ci in seq_len(N_CHR)) {
  chr_label <- as.character(ci)
  bs        <- BLOCK_SIZES[[ci]]
  cur_pos   <- 1000L
  chr_geno  <- NULL
  snp_ids   <- character(0)
  snp_pos   <- integer(0)
  snp_ref   <- character(0)
  snp_alt   <- character(0)
  snp_count <- 0L

  for (bi in seq_along(bs)) {
    n_b <- bs[bi]
    p_b <- runif(1, 0.20, 0.45)
    blk <- make_block(N_IND, n_b, p_b, flip = 0.04)
    chr_geno <- cbind(chr_geno, blk)

    for (j in seq_len(n_b)) {
      snp_count <- snp_count + 1L
      al <- ALLELES[(snp_count %% nrow(ALLELES)) + 1L, ]
      snp_ids <- c(snp_ids, sprintf("rs%s%03d", chr_label, snp_count))
      snp_pos <- c(snp_pos, cur_pos)
      snp_ref <- c(snp_ref, al[1])
      snp_alt <- c(snp_alt, al[2])
      cur_pos <- cur_pos + sample(800L:1200L, 1L)
    }

    # Inter-block singletons (except after last block)
    if (bi < length(bs)) {
      sing <- make_singleton(N_IND, INTER_BLOCK_SNP)
      chr_geno <- cbind(chr_geno, sing)
      cur_pos <- cur_pos + 50000L   # genomic gap
      for (j in seq_len(INTER_BLOCK_SNP)) {
        snp_count <- snp_count + 1L
        al <- ALLELES[(snp_count %% nrow(ALLELES)) + 1L, ]
        snp_ids <- c(snp_ids, sprintf("rs%s%03d", chr_label, snp_count))
        snp_pos <- c(snp_pos, cur_pos)
        snp_ref <- c(snp_ref, al[1])
        snp_alt <- c(snp_alt, al[2])
        cur_pos <- cur_pos + sample(800L:1200L, 1L)
      }
    }
  }

  geno_list[[ci]] <- chr_geno
  snpinfo_list[[ci]] <- data.frame(
    SNP = snp_ids, CHR = chr_label,
    POS = snp_pos, REF = snp_ref, ALT = snp_alt,
    stringsAsFactors = FALSE
  )
}

ldx_geno     <- do.call(cbind, geno_list)
ldx_snp_info <- do.call(rbind, snpinfo_list)
rownames(ldx_geno) <- paste0("ind", sprintf("%03d", seq_len(N_IND)))
colnames(ldx_geno) <- ldx_snp_info$SNP

stopifnot(ncol(ldx_geno) == nrow(ldx_snp_info))
message("Simulated genotype matrix: ", nrow(ldx_geno), " ind x ", ncol(ldx_geno), " SNPs")

# ── Pre-compute block table (stored as example output) ────────────────────────
# Use in-memory path (no external deps at data-raw stage)
# We derive a simple block table by chromosome for illustration.
# Full run: blocks <- run_Big_LD_all_chr(ldx_geno, ldx_snp_info, CLQcut=0.6)
# We build a hand-crafted reference table matching the simulated structure.

make_block_row <- function(snp_info, snp_ids) {
  idx <- match(snp_ids, snp_info$SNP)
  data.frame(
    start      = min(idx), end = max(idx),
    start.rsID = snp_info$SNP[min(idx)],
    end.rsID   = snp_info$SNP[max(idx)],
    start.bp   = snp_info$POS[min(idx)],
    end.bp     = snp_info$POS[max(idx)],
    CHR        = snp_info$CHR[min(idx)],
    length_bp  = snp_info$POS[max(idx)] - snp_info$POS[min(idx)] + 1L,
    stringsAsFactors = FALSE
  )
}

# Derive expected block SNP ranges from simulation structure
chr1_snps <- ldx_snp_info$SNP[ldx_snp_info$CHR == "1"]
chr2_snps <- ldx_snp_info$SNP[ldx_snp_info$CHR == "2"]
chr3_snps <- ldx_snp_info$SNP[ldx_snp_info$CHR == "3"]

ldx_blocks <- rbind(
  # chr1: 3 blocks of sizes 25, 20, 25 with 5-SNP gaps between
  make_block_row(ldx_snp_info, chr1_snps[1:25]),
  make_block_row(ldx_snp_info, chr1_snps[31:50]),
  make_block_row(ldx_snp_info, chr1_snps[56:80]),
  # chr2: 3 blocks of sizes 30, 20, 20
  make_block_row(ldx_snp_info, chr2_snps[1:30]),
  make_block_row(ldx_snp_info, chr2_snps[36:55]),
  make_block_row(ldx_snp_info, chr2_snps[61:80]),
  # chr3: 3 blocks of sizes 20, 20, 20
  make_block_row(ldx_snp_info, chr3_snps[1:20]),
  make_block_row(ldx_snp_info, chr3_snps[26:45]),
  make_block_row(ldx_snp_info, chr3_snps[51:70])
)
rownames(ldx_blocks) <- NULL

message("Pre-computed block table: ", nrow(ldx_blocks), " blocks")

# ── Toy GWAS marker table for tune_LD_params demos ───────────────────────────
set.seed(42)
gwas_snp_idx <- c(
  sample(1:25, 3),           # from chr1 block1
  sample(31:50, 3),          # from chr1 block2
  sample(56:80, 3),          # from chr1 block3
  sample(which(ldx_snp_info$CHR=="2")[1:30], 4),
  sample(which(ldx_snp_info$CHR=="2")[36:55], 3),
  sample(which(ldx_snp_info$CHR=="3")[1:20], 4)
)
gwas_snp_idx <- sort(unique(gwas_snp_idx))[1:20]

ldx_gwas <- data.frame(
  Marker = ldx_snp_info$SNP[gwas_snp_idx],
  CHR    = ldx_snp_info$CHR[gwas_snp_idx],
  POS    = ldx_snp_info$POS[gwas_snp_idx],
  P      = runif(20, 1e-8, 0.001),
  stringsAsFactors = FALSE
)

message("GWAS marker table: ", nrow(ldx_gwas), " markers")

# ── Save .rda to data/ ────────────────────────────────────────────────────────
usethis::use_data(ldx_geno,     overwrite = TRUE, compress = "xz")
usethis::use_data(ldx_snp_info, overwrite = TRUE, compress = "xz")
usethis::use_data(ldx_blocks,   overwrite = TRUE, compress = "xz")
usethis::use_data(ldx_gwas,     overwrite = TRUE, compress = "xz")

message("Saved all .rda files to data/")

# ── Write flat-file copies to inst/extdata/ ──────────────────────────────────

## 1. Numeric dosage CSV (SNPs x samples: SNP CHR POS REF ALT + sample cols)
numeric_df <- cbind(
  ldx_snp_info,
  as.data.frame(t(ldx_geno))
)
write.csv(numeric_df,
          file = "inst/extdata/example_genotypes_numeric.csv",
          row.names = FALSE, quote = FALSE
)

## 2. HapMap format
# Iterate over SNPs (rows of ldx_geno after transposing back to SNPs x samples)
# ldx_geno is individuals x SNPs, so columns = SNPs
hmp_decode_snp <- function(snp_dosages, ref, alt) {
  # snp_dosages: numeric vector length n_ind
  vapply(snp_dosages, function(x) {
    if (is.na(x)) return("NN")
    switch(as.character(as.integer(x)),
           "0" = paste0(ref, ref),
           "1" = paste0(ref, alt),
           "2" = paste0(alt, alt),
           "NN"
    )
  }, character(1))
}

# Build sample call matrix: n_snp rows x n_ind cols
n_snp <- nrow(ldx_snp_info)
n_ind <- nrow(ldx_geno)
hmp_calls <- matrix("NN", nrow = n_snp, ncol = n_ind)
for (si in seq_len(n_snp)) {
  hmp_calls[si, ] <- hmp_decode_snp(
    ldx_geno[, si],
    ldx_snp_info$REF[si],
    ldx_snp_info$ALT[si]
  )
}

hmp_header <- data.frame(
  "rs#"       = ldx_snp_info$SNP,
  alleles     = paste0(ldx_snp_info$REF, "/", ldx_snp_info$ALT),
  chrom       = ldx_snp_info$CHR,
  pos         = ldx_snp_info$POS,
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
hmp_out <- cbind(hmp_header, as.data.frame(hmp_calls, stringsAsFactors = FALSE))
colnames(hmp_out)[12:ncol(hmp_out)] <- rownames(ldx_geno)
write.table(hmp_out,
            file = "inst/extdata/example_genotypes.hmp.txt",
            sep = "\t", row.names = FALSE, quote = FALSE
)

## 3. VCF
vcf_lines <- c(
  "##fileformat=VCFv4.2",
  paste0("##source=LDxBlocks example data, generated ", Sys.Date()),
  "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">",
  paste(c("#CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT",
          rownames(ldx_geno)), collapse = "\t")
)
gt_encode <- function(g) {
  vapply(g, function(x) {
    if (is.na(x)) return("./.")
    switch(as.character(x), "0"="0/0", "1"="0/1", "2"="1/1", "./.")
  }, character(1))
}
for (si in seq_len(nrow(ldx_snp_info))) {
  gt_vec <- gt_encode(ldx_geno[, si])
  vcf_lines <- c(vcf_lines, paste(c(
    ldx_snp_info$CHR[si], ldx_snp_info$POS[si], ldx_snp_info$SNP[si],
    ldx_snp_info$REF[si], ldx_snp_info$ALT[si], ".", "PASS", ".", "GT",
    gt_vec
  ), collapse = "\t"))
}
writeLines(vcf_lines, "inst/extdata/example_genotypes.vcf")

message("Written flat-file examples to inst/extdata/")
message("Data generation complete.")
