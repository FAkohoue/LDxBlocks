#' @title Chromosome-Wise LD Block Detection Using Kinship-Adjusted rV²
#'
#' @description
#' Executes the \code{Big_LD()} function chromosome by chromosome across the genome.
#' This wrapper is optimized to detect LD blocks using kinship-adjusted squared correlation (\eqn{rV^2}),
#' based on the Big-LD algorithm (Kim et al., 2017), extended for structured populations.
#'
#' @param geno_matrix A numeric genotype matrix (individuals × SNPs), including all SNPs genome-wide.
#' @param snp_info A \code{data.frame} or \code{data.table} with columns \code{CHR}, \code{SNP}, and \code{POS}, corresponding to chromosome, SNP ID, and physical position.
#'                 The number of rows must match the number of columns in \code{geno_matrix}.
#' @param CLQcut Threshold for rV² correlation (0–1) to define cliques. Default is 0.5.
#' @param clstgap Maximum physical gap (in bp) allowed within cliques. Default is 40000.
#' @param leng SNP window size to detect weak-LD segments. Default is 200.
#' @param subSegmSize Maximum SNPs allowed in a subregion before forced segmentation. Default is 1500.
#' @param MAFcut Minimum minor allele frequency for SNP inclusion. Default is 0.05.
#' @param appendrare Logical. If TRUE, rare SNPs (MAF < MAFcut) are appended after LD block construction. Default is FALSE.
#' @param checkLargest Logical. If TRUE, uses a dense-core decomposition to accelerate clique extraction for large SNP sets. Default is FALSE.
#' @param CLQmode Clique scoring method: either "Density" (default) or "Maximal".
#' @param rV2method Method to compute V^{-1/2}: either "chol" (Cholesky) or "eigen". Default is "chol".
#' @param split Logical. If TRUE, apply gap-based splitting to divide long cliques across distant SNPs.
#'
#' @details
#' This function iterates through each chromosome defined in \code{snp_info} and runs the internal \code{\link{Big_LD}} algorithm using kinship-adjusted squared correlations (\eqn{rV^2}).
#' The approach is robust to population structure and improves LD block detection in structured or related populations.
#'
#' Chromosomes with fewer than 10 SNPs are skipped with a message.
#'
#' @return A data.frame with one row per detected LD block and the following columns:
#' \describe{
#'   \item{start}{Start index of the block (in genome-wide SNP index)}
#'   \item{end}{End index of the block}
#'   \item{start.rsID}{SNP ID of the first SNP in the block}
#'   \item{end.rsID}{SNP ID of the last SNP in the block}
#'   \item{start.bp}{Physical position of the first SNP}
#'   \item{end.bp}{Physical position of the last SNP}
#'   \item{CHR}{Chromosome identifier}
#'   \item{length_bp}{Physical length of the block in base pairs}
#' }
#'
#' @references
#' Kim S, Scheet P, Kang EY. (2017). A fast and simple algorithm for detecting identity-by-descent segments. \emph{Bioinformatics}, 33(16):2526–2534. \doi{10.1093/bioinformatics/btx609}
#'
#' @seealso \code{\link{Big_LD}}, \code{\link{CLQD}}, \code{\link{compute_rV2}}, \code{\link{get_V_inv_sqrt}}
#'
#' @examples
#' # Assuming geno_matrix and snp_info are properly loaded
#' blocks <- run_Big_LD_all_chr(geno_matrix, snp_info, CLQcut = 0.5, split = TRUE)
#' head(blocks)
#'
#' @author Félicien Akohoue
#' @export
#' 
run_Big_LD_all_chr <- function(geno_matrix, snp_info, 
                               CLQcut = 0.5, 
                               clstgap = 40000, 
                               leng = 200, 
                               subSegmSize = 1500, 
                               MAFcut = 0.05, 
                               appendrare = FALSE, 
                               checkLargest = FALSE, 
                               CLQmode = "Density", 
                               rV2method = "chol", 
                               split = FALSE) {
  chromosomes <- unique(snp_info$CHR)
  ld_blocks_all <- list()
  
  for (chr in chromosomes) {
    cat("\nRunning Big_LD for", chr, "...\n")
    chr_snps <- which(snp_info$CHR == chr)
    geno_chr <- geno_matrix[, chr_snps, drop = FALSE]
    SNPinfo_chr <- snp_info[chr_snps, c("SNP", "POS")]
    SNPinfo_chr <- as.data.frame(SNPinfo_chr)
    
    if (ncol(geno_chr) < 10) {
      cat("  Skipped", chr, "- too few SNPs.\n")
      next
    }
    
    block_chr <- tryCatch({
      Big_LD(geno = geno_chr,
             SNPinfo = SNPinfo_chr,
             CLQcut = CLQcut, 
             clstgap = clstgap, 
             leng = leng, 
             subSegmSize = subSegmSize, 
             MAFcut = MAFcut, 
             appendrare = FALSE, 
             checkLargest = FALSE, 
             CLQmode =CLQmode, rV2method = rV2method, split = FALSE)
    }, error = function(e) {
      cat("  Error in Big_LD for", chr, ":", conditionMessage(e), "\n")
      return(NULL)
    })
    
    if (!is.null(block_chr)) {
      block_chr$CHR <- chr
      ld_blocks_all[[chr]] <- block_chr
    }
  }
  
  all_blocks <- do.call(rbind, ld_blocks_all)
  all_blocks$start.rsID <- as.character(all_blocks$start.rsID)
  all_blocks$end.rsID   <- as.character(all_blocks$end.rsID)
  all_blocks$length_bp  <- all_blocks$end.bp - all_blocks$start.bp + 1
  
  return(all_blocks)
}
