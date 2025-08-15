#' Square root of a matrix, by spectral decomposition
#'
#' @description: Accepts a matrix as input and returns the square root using spectral decomposition
#' 
#' @param A: a matrix
#'
#' @return A matrix
#' @export
#' 
#' @examples
#' vcov <- matrix(c(2, 0.25, 0.25, 0.5), 2, 2)
#' mat_spec_sqrt(vcov)
#' 
mat_spec_sqrt <- function(A) {
  eDecomp <- eigen(A) 
  return(eDecomp$vectors %*% diag(sqrt(eDecomp$values)))
}
