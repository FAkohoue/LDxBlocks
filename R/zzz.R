#' @useDynLib LDxBlocks, .registration = TRUE
#' @importFrom Rcpp evalCpp
NULL

.onUnload <- function(libpath) {
  library.dynam.unload("LDxBlocks", libpath)
}
