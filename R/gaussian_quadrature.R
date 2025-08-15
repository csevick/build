

#' Basic Gaussian Quadrature
#'
#' @description: This function implements integration of a function by Gaussian quadrature, for supplied points and weights 
#'
#' @param FUNCT: name of a function to integrate
#' @param PW: a list of points and weights 
#' @param ... : further arguments to the user supplied function in FUNCT
#'
#' @return A scalar result of the integration
#'
#' @export
#' 
#' @examples
#' 
#' # gauss-hermite quadrature: integrate (1+x^2)*exp((-x^2)/2)/sqrt(2*pi)
#' f <- function(x) (1+x^2) 

#' # generate points and weights
#' ghpw <- gauss_herm_consts(3)

#' # integrate
#' gquad(f, ghpw)

#' # check
#' f2 <- function(x) (1+x^2)*exp((-x^2)/2)/sqrt(2*pi)

#' integrate(f2, -Inf, Inf)

gquad <- function(FUNCT, PW, ...) {
  
  as.numeric(PW$wts %*% apply(PW$pnts, MARGIN = 2, FUNCT, ... ))
  
}

#' coefficients of the Hermite polynomial (Salzer 1952)
#'
#' @description: This function generates coefficients of the Hermite polynomial (Salzer 1952)
#' Salzer H. E., Zucker R., Capuano R. "Table of the Zero and Weight Factors of the First Twenty Hermite Polynomials" Journal of Research of the National Bureau of Standards, Vol 48, No. 2, February 1952
#' @param n: number of points to generate
#'
#' @return a numeric vector
#'
#' @export
#' 
#' @examples
#' herm_coef(5)

herm_coef <- function(n) {
  
  H <- rep(0, n + 1 )
  for(k in 0:(n/2)) {
    H[n-2*k + 1] <- ((-1)^k)*factorial(n)*(2^(n-2*k))/( factorial(k)*factorial(n-2*k) ) 
  }
  H
}

#' Weights for the Gauss-Hermite quadrature points (Salzer 1952)
#'
#' @description: This function generates weights for the Gauss-Hermite quadrature points (Salzer 1952).
#' Salzer H. E., Zucker R., Capuano R. "Table of the Zero and Weight Factors of the First Twenty Hermite Polynomials" Journal of Research of the National Bureau of Standards, Vol 48, No. 2, February 1952
#' @param x: coefficients of the Hermite polynomial
#' @param n: Number of coefficients of the Hermite polynomial
#'
#' @return a numeric vector
#'
#' @export
#' 
#' @examples
#' hx <- Re(polyroot(herm_coef(4)))
#' sapply(hx, herm_wt, n = 4)

herm_wt <- function(x, n){
  sqrt(pi)*(2^(n-1))*factorial(n-1) / (n*(t(herm_coef(n-1)) %*% sapply(0:(n-1), function(p){x^p}))^2)
}

#' Points and weight for Gauss-Hermite quadrature
#'
#' @description: This function generates a list of points an weights for Gauss-Hermite quadrature. Since it is anticipated that general use will be to integrate with respect to a normal density, the points are pre-multiplied by the square root of 2 and the weights divided by the square root of pi. 
#' 
#' @param n: number of points to generate
#' @param ndim: number of points to generate the dimension of the integration
#'
#' @return A list with numeric matrix components:
#' 1) **pnts**: quadrature points (ndim by n)
#' 2) **wts**: quadrature weights (1 by n)
#' 
#' @export
#' 
#' @examples
#' gauss_herm_consts(5)

gauss_herm_consts <- function(n, ndim = 1) {
  
  pnts <- Re(polyroot(herm_coef(n)))
  pnts <- matrix(sort(pnts), nrow = 1, ncol = n)
  
  wts <-sapply(pnts, herm_wt, n = n)
  
  if(ndim>1) {
    # npoints by ndim matrix of all combinations of points
    pnts <- t(as.matrix(expand.grid(rep(list(as.vector(pnts)), ndim))))
    # creates a matrix of weights in the right order to reflect x, then calculates the row product
    wts <- apply(expand.grid(rep(list(wts), ndim)), 1, prod) 
  }
  
  list(
    pnts = pnts*sqrt(2)
  ,wts =  matrix(wts, nrow = 1, ncol = length(wts))/sqrt(pi)
  )
  
}

#' Center and scale quadrature points
#'
#' @description: This function accepts an object from  gauss_herm_consts() and scales the points for a supplied mean vector and variance matrix
#' 
#' @param pw: a point and weight object from gauss_herm_consts()
#' @param mu: a mean vector used to center the points
#' @param vcov: a variance matrix used to scale the points
#'
#' @return A centered and scaled version of the input from gauss_herm_consts()
#' 
#' @export
#' 
#' @examples
#' ghpw <- gauss_herm_consts(3, 2)
#' mu <- c(1, 2)
#' vcov <- matrix(c(2, 0.25, 0.25, 0.5), 2, 2)
#' ghpw_new <- gquad_cntr_scale(ghpw, mu, vcov)

gquad_cntr_scale <- function(pw, mu, vcov) {
  
  if(missing(mu)) {mu <- rep(0, ncol(vcov))}

  pw[[1]] <- 
  if(nrow(vcov)>1) {

    mu + mat_spec_sqrt(vcov) %*% pw[[1]]
    
  } else {
    mu + pw[[1]]*sqrt(as.numeric(vcov))
  }
  
  pw
}


 
