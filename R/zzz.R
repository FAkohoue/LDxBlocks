#' @useDynLib LDxBlocks, .registration = TRUE
#' @importFrom Rcpp evalCpp
#' @importFrom digest digest
NULL

.onUnload <- function(libpath) {
  library.dynam.unload("LDxBlocks", libpath)
}
