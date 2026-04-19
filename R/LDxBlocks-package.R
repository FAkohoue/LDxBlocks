#' LDxBlocks: Genome-Wide LD Block Detection and Haplotype Analysis
#'
#' @description
#' Fast, scalable linkage disequilibrium (LD) block detection for genome-wide
#' SNP data using a C++/Armadillo computational core and two LD metrics:
#' standard r^2 and kinship-adjusted rV^2. Includes haplotype reconstruction,
#' diversity metrics, prediction feature matrices, and multi-format genotype I/O.
#'
#' @section Choosing between r^2 and rV^2:
#' \describe{
#'   \item{r^2 (\code{method = "r2"})}{Standard squared Pearson correlation.
#'     No kinship matrix required. C++ kernel via RcppArmadillo + OpenMP.
#'     Suitable for unstructured populations and large datasets (10 M+ markers).}
#'   \item{rV^2 (\code{method = "rV2"})}{Kinship-whitened squared correlation
#'     (Mangin et al. 2012, Heredity 108:285-291). Corrects LD inflation from relatedness. For livestock,
#'     inbred lines, or family-based cohorts with < 200 k markers per chromosome.}
#' }
#'
#' @section Pipeline overview:
#' \enumerate{
#'   \item Read genotype data: \code{\link{read_geno}} -- auto-detects CSV,
#'     HapMap, VCF, GDS, BED, or plain matrix. For WGS panels where peak RAM
#'     is a concern, use \code{\link{read_geno_bigmemory}} to build a
#'     memory-mapped store (requires \pkg{bigmemory}).
#'   \item Detect LD blocks chromosome-wise: \code{\link{run_Big_LD_all_chr}}.
#'   \item Optionally auto-tune parameters: \code{\link{tune_LD_params}}.
#'   \item Reconstruct haplotypes: \code{\link{extract_haplotypes}}.
#'   \item Decode haplotype strings to nucleotides:
#'     \code{\link{decode_haplotype_strings}}.
#'   \item Compute diversity metrics: \code{\link{compute_haplotype_diversity}}.
#'   \item Build prediction feature matrix:
#'     \code{\link{build_haplotype_feature_matrix}}.
#'   \item Compute LD decay and chromosome-specific decay distances:
#'     \code{\link{compute_ld_decay}}. Provides the critical r\eqn{^2}
#'     threshold (parametric: 95th percentile of unlinked-marker r\eqn{^2})
#'     and per-chromosome decay distances for candidate gene windows.
#'   \item Map GWAS hits to QTL regions (LD-aware windows when
#'     \code{ld_decay} is supplied): \code{\link{define_qtl_regions}}.
#'   \item Genomic prediction (Tong et al. 2024/2025):
#'     \code{\link{run_haplotype_prediction}} -- full pipeline from BLUEs to
#'     block importance; or step-by-step via
#'     \code{\link{compute_haplotype_grm}},
#'     \code{\link{backsolve_snp_effects}},
#'     \code{\link{compute_local_gebv}},
#'     \code{\link{rank_haplotype_blocks}}.
#'   \item Write outputs: \code{\link{write_haplotype_numeric}},
#'     \code{\link{write_haplotype_character}},
#'     \code{\link{write_haplotype_diversity}}.
#'   \item Or run everything at once: \code{\link{run_ldx_pipeline}}.
#' }
#'
#' @section Example data:
#' \describe{
#'   \item{\code{\link{ldx_geno}}}{120 x 230 genotype matrix, 3 chromosomes,
#'     9 simulated LD blocks.}
#'   \item{\code{\link{ldx_snp_info}}}{SNP metadata (SNP, CHR, POS, REF, ALT).}
#'   \item{\code{\link{ldx_blocks}}}{Reference block table for \code{ldx_geno}.}
#'   \item{\code{\link{ldx_gwas}}}{20 toy GWAS markers for tuning demos.}
#'   \item{\code{\link{ldx_blues}}}{Pre-adjusted BLUEs (id, YLD, RES) for
#'     \code{\link{run_haplotype_prediction}} demos. Available as an .rda
#'     dataset and \code{inst/extdata/example_blues.csv}.}
#' }
#'
#' @references
#' Kim S-A et al. (2018). A new haplotype block detection method for dense
#' genome sequencing data based on interval graph modeling and dynamic
#' programming. \emph{Bioinformatics} \strong{34}(4):588-596.
#' \doi{10.1093/bioinformatics/btx609}
#'
#' Mangin B et al. (2012). Novel measures of linkage disequilibrium that
#' correct the bias due to population structure and relatedness.
#' \emph{Heredity} \strong{108}(3):285-291. \doi{10.1038/hdy.2011.73}
#'
#' VanRaden PM (2008). Efficient methods to compute genomic predictions.
#' \emph{Journal of Dairy Science} \strong{91}(11):4414-4423.
#' \doi{10.3168/jds.2007-0980}
#'
#' Difabachew YF et al. (2023). Genomic prediction with haplotype blocks
#' in wheat. \emph{Frontiers in Plant Science} \strong{14}:1168547.
#' \doi{10.3389/fpls.2023.1168547}
#'
#' Weber SE, Frisch M, Snowdon RJ, Voss-Fels KP (2023). Haplotype blocks
#' for genomic prediction: a comparative evaluation in multiple crop datasets.
#' \emph{Frontiers in Plant Science} \strong{14}:1217589.
#' \doi{10.3389/fpls.2023.1217589}
#'
#' Pook T et al. (2019). HaploBlocker: Creation of subgroup-specific
#' haplotype blocks and libraries. \emph{Genetics} \strong{212}(4):1045-1061.
#' \doi{10.1534/genetics.119.302283}
#'
#' Tong J et al. (2024). Stacking beneficial haplotypes from the Vavilov
#' wheat collection to accelerate breeding for multiple disease resistance.
#' \emph{Theoretical and Applied Genetics} \strong{137}:274.
#' \doi{10.1007/s00122-024-04784-w}
#'
#' Tong J et al. (2025). Haplotype stacking to improve stability of stripe
#' rust resistance in wheat. \emph{Theoretical and Applied Genetics}
#' \strong{138}:267. \doi{10.1007/s00122-025-05045-0}
#'
#'
#' Blondel VD, Guillaume J-L, Lambiotte R, Lefebvre E (2008). Fast unfolding
#' of communities in large networks. \emph{Journal of Statistical Mechanics:
#' Theory and Experiment} \strong{2008}:P10008.
#' \doi{10.1088/1742-5468/2008/10/P10008}
#'
#' Traag VA, Waltman L, van Eck NJ (2019). From Louvain to Leiden: guaranteeing
#' well-connected communities. \emph{Scientific Reports} \strong{9}:5233.
#' \doi{10.1038/s41598-019-41695-z}
#'
#' Hill WG, Weir BS (1988). Variances and covariances of squared linkage
#' disequilibria in finite populations. \emph{Theoretical Population Biology}
#' \strong{33}(1):54-78. \doi{10.1016/0040-5809(88)90004-4}
#'
#' Remington DL et al. (2001). Structure of linkage disequilibrium and
#' phenotypic associations in the maize genome. \emph{PNAS}
#' \strong{98}(20):11479-11484. \doi{10.1073/pnas.201394398}
#' @docType package
#' @importFrom graphics legend
#' @importFrom utils write.table
#' @name LDxBlocks-package
#' @aliases LDxBlocks
#' @keywords package
#'
#' @importFrom rrBLUP kin.blup mixed.solve
#' @import Rcpp
#' @importFrom data.table fread fwrite rbindlist setnames setorder :=
#' @importFrom stats cor median quantile na.omit setNames var aggregate sd lm coef
#'   residuals pt chisq.test as.dist hclust cutree anova pf prcomp
#'   loess nls nls.control predict
#' @importFrom igraph graph_from_adjacency_matrix coreness max_cliques cliques components
#' @importFrom igraph cluster_louvain cluster_leiden membership as_edgelist
#' @importFrom igraph induced_subgraph is_connected mst E V vcount layout_nicely
#' @importFrom grDevices hcl.colors
#' @importFrom utils globalVariables head read.table tail
"_PACKAGE"

utils::globalVariables(c(
  "CHR", "start.bp", "end.bp", "start", "end",
  "block_idx", ".N", "block_name", "length_snps", "length_bp",
  "n_unassigned", "n_forced", "n_blocks", "penalty_bp",
  "REF", "ALT", "SNP", "POS", "..samp_cols", ".",
  "start_x", "end_x", "chr_int",
  "dist_kb", "r2_plot", "r2",
  "ci"  # phased-path loop variable in extract_haplotypes()
))
