#' @import data.table
#' @importFrom stats cov median quantile na.omit
#' @importFrom igraph graph_from_adjacency_matrix coreness max_cliques cliques
#' @importFrom utils globalVariables
NULL

utils::globalVariables(c(
  "CHR","start.bp","block_idx",".N","block_name","length_snps",
  "end","start","length_bp","end.bp",
  "n_unassigned","n_forced","n_blocks","penalty_bp"
))
