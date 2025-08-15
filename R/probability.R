
#' Probability density function for MV normal
#'
#' @description: Probability density function for MV normal for arbitrary mean vector and variance matrix.
#' If mu and sigma are not specified then a standard MV normal distribution, of the same dimension as x, will be assumed.
#' 
#' @param x: numeric - a vector of observed points for which to compute the pdf value (real numbers)
#' @param mu: numeric - a vector of distribution mean values (real numbers)
#' @param sigma: numeric - distribution variance matrix (\eqn{\gt 0})
#'
#' @return A scalar representing the height of the normal surface at x
#' @export
#' 
#' @examples
#' 
#' dmvnorm(c(0.25, 0.05), mu = c(0.1, 0.2), sigma = matrix(c(1, 0.25, 0.25, 0.8), 2, 2))

dmvnorm <- function(x # numeric vector
                    , mu = NULL # numeric vector of means
                    , sigma = NULL # variance matrix
                    ) {

  if(is.null(mu)) {mu <- rep(0, length(x))}
  if(is.null(sigma)) {sigma <- diag(length(x))}
  as.numeric(((2*pi)^(-length(x)/2))*(det(sigma)^-0.5) * exp(-0.5*t(x-mu) %*% solve(sigma) %*% (x-mu)))

}



#' Probability density function for a beta-binomial distribution
#'
#' @description: Probability density function for a beta-binomial distribution given by:
#' 
#' \eqn{{n \choose k} \frac{\Gamma(a+b)\Gamma(x+a)\Gamma(n-x+b)}{\Gamma(a)\Gamma(b)\Gamma(a+b+n)}}
#' 
#' Where \eqn{\Gamma()} is the gamma function.
#' 
#' Reference: Casella G., Berger R. L. "Statistical Inference, Second Edition" Duxbury Press, 2002
#' 
#' @param n: numeric scalar - the number of trials (\eqn{\ge 0})
#' @param x: numeric scalar - the number of successes (\eqn{\ge 0})
#' @param a: numeric scalar - shape parameter (\eqn{\gt 0})
#' @param b: numeric scalar - shape parameter (\eqn{\gt 0})
#'
#' @return A scalar representing the PDF value at chosen parameters
#' @export
#' 
#' @examples
#'
#' sapply(0:4, function(x){pbetabinom(4,x,1.5,1)})
#' 

pbetabinom <- function(n, x, a, b) {
  choose(n, x)*(gamma(a+b)/(gamma(a)*gamma(b)))*gamma(x+a)*gamma(n-x+b)/gamma(a+b+n)
}


