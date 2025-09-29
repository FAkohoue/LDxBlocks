#' Chromosome-Wise LD Block Detection Using Kinship-Adjusted \eqn{rV^2}
#'
#' @description
#' Executes \code{Big_LD()} chromosome by chromosome across the genome.
#'
#' @param geno_matrix Numeric genotype matrix (individuals x SNPs), genome-wide.
#' @param snp_info Data frame with columns \code{CHR}, \code{SNP}, \code{POS}.
#' @param CLQcut,clstgap,leng,subSegmSize,MAFcut,appendrare,checkLargest,CLQmode,rV2method,split,digits,seed,verbose
#'   Passed to \code{Big_LD()}.
#'
#' @return A data frame of LD blocks with columns:
#' \code{start,end,start.rsID,end.rsID,start.bp,end.bp,CHR,length_bp}
#'
#' @examples
#' \donttest{
#' set.seed(1)
#' geno <- matrix(sample(0:2, 80*300, replace=TRUE), 80, 300)
#' snp_info <- data.frame(
#'   CHR = rep(paste0("chr", 1:3), each = 100),
#'   SNP = paste0("rs", 1:300),
#'   POS = ave(sample(1e6, 300), rep(1:3, each=100), FUN = sort)
#' )
#' blocks <- run_Big_LD_all_chr(geno, snp_info, CLQcut = 0.6, verbose = FALSE)
#' head(blocks)
#' }
#' @export
run_Big_LD_all_chr <- function(
    geno_matrix, snp_info,
    CLQcut = 0.5, clstgap = 40000, leng = 200, subSegmSize = 1500, MAFcut = 0.05,
    appendrare = FALSE, checkLargest = FALSE, CLQmode = "Density",
    rV2method = "chol", split = FALSE, digits = 6, seed = NULL, verbose = FALSE
) {
  chromosomes <- unique(as.character(snp_info$CHR))
  ld_blocks_all <- vector("list", length(chromosomes)); names(ld_blocks_all) <- chromosomes

  for (chr in chromosomes) {
    if (isTRUE(verbose)) cat("\n[run_Big_LD_all_chr] Running Big_LD for", chr, "...\n")
    idx <- which(snp_info$CHR == chr)
    geno_chr <- geno_matrix[, idx, drop = FALSE]
    info_chr <- as.data.frame(snp_info[idx, c("SNP","POS")])

    if (ncol(geno_chr) < 10) {
      if (isTRUE(verbose)) cat("  Skipped", chr, "- too few SNPs.\n")
      next
    }
    block_chr <- tryCatch({
      Big_LD(geno = geno_chr, SNPinfo = info_chr, CLQcut = CLQcut, clstgap = clstgap,
             leng = leng, subSegmSize = subSegmSize, MAFcut = MAFcut, appendrare = appendrare,
             checkLargest = checkLargest, CLQmode = CLQmode, rV2method = rV2method, split = split,
             digits = digits, seed = seed, verbose = verbose)
    }, error = function(e) {
      if (isTRUE(verbose)) cat("  Error in Big_LD for", chr, ":", conditionMessage(e), "\n")
      NULL
    })
    if (!is.null(block_chr)) {
      block_chr$CHR <- chr
      ld_blocks_all[[chr]] <- block_chr
    }
  }

  all_blocks <- data.table::rbindlist(ld_blocks_all, use.names = TRUE, fill = TRUE)
  if (nrow(all_blocks)) {
    all_blocks$start.rsID <- as.character(all_blocks$start.rsID)
    all_blocks$end.rsID   <- as.character(all_blocks$end.rsID)
    all_blocks$length_bp  <- all_blocks$end.bp - all_blocks$start.bp + 1
  }
  as.data.frame(all_blocks)
}
