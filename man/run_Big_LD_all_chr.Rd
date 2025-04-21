\name{run_Big_LD_all_chr}
\alias{run_Big_LD_all_chr}
\title{
Chromosome-Wise LD Block Detection Using Kinship-Adjusted rV²
}
\description{
Executes the Big_LD() function chromosome by chromosome across the genome. Designed for large-scale LD block detection with kinship-adjusted squared correlations (\eqn{rV^2}).
}
\usage{
run_Big_LD_all_chr(geno_matrix, snp_info,
  CLQcut = 0.5,
  clstgap = 40000,
  leng = 200,
  subSegmSize = 1500,
  MAFcut = 0.05,
  appendrare = FALSE,
  checkLargest = FALSE,
  CLQmode = "Density",
  rV2method = "chol",
  split = FALSE)
}
\arguments{
  \item{geno_matrix}{Genotype matrix of size n × p (individuals × SNPs).}
  \item{snp_info}{A data frame or data table with columns: CHR (chromosome), SNP (ID), POS (base-pair position). Must align with geno_matrix columns.}
  \item{CLQcut}{Threshold for rV² correlation (0–1). Default is 0.5.}
  \item{clstgap}{Maximum gap (bp) within cliques. Default: 40000.}
  \item{leng}{SNP window size to scan for weak LD segments.}
  \item{subSegmSize}{Max SNPs allowed in a subregion before segmentation.}
  \item{MAFcut}{Minor allele frequency threshold. SNPs below this are excluded unless appendrare = TRUE.}
  \item{appendrare}{Logical. Append rare SNPs after block construction?}
  \item{checkLargest}{Enable heuristic for dense-core decomposition (for large SNP sets).}
  \item{CLQmode}{Clique scoring: "Density" (default) or "Maximal".}
  \item{rV2method}{Method for computing \eqn{V^{-1/2}}: "chol" or "eigen".}
  \item{split}{Split cliques across large physical gaps?}
}
\details{
This function wraps Big_LD for per-chromosome LD block detection, making the algorithm scalable across genome-wide data.
Chromosomes with fewer than 10 SNPs are skipped.
}
\value{
A data.frame with columns:
\itemize{
  \item \code{start}, \code{end}: SNP indices (global)
  \item \code{start.rsID}, \code{end.rsID}: SNP IDs at block boundaries
  \item \code{start.bp}, \code{end.bp}: Base-pair positions
  \item \code{CHR}: Chromosome
  \item \code{length_bp}: Physical length of the block (in bp)
}
}
\references{
Kim S, Scheet P, Kang EY. (2017). \emph{Bioinformatics}, 33(16):2526–2534. \doi{10.1093/bioinformatics/btx609}
}

\keyword{LD}
\keyword{genomics}


\references{Félicien Akohoue
%% ~put references to the literature/web site here ~
}
\author{Félicien Akohoue
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{\code{\link{Big_LD}}, \code{\link{compute_rV2}}, \code{\link{get_V_inv_sqrt}}, \code{\link{CLQD}}}

\examples{
blocks <- run_Big_LD_all_chr(geno_matrix, snp_info, CLQcut = 0.5, split = TRUE)
head(blocks)

##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or standard data sets, see data().

## The function is currently defined as

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

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
