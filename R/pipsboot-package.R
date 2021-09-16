## usethis namespace: start
#' @useDynLib pipsboot, .registration = TRUE
#' @importFrom Rcpp sourceCpp
## usethis namespace: end
NULL

# cleans up dll on package unload
.onUnload <- function(libpath) {
  library.dynam.unload("pipsboot", libpath)
}
