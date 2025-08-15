#' Generate a simulated data set from a GLMM model with non-random dropout

#' @description 
#' Generate a simulated data set from a GLMM model with non-random dropout. For N subjects (indexed by \eqn{i}),
#' each with \eqn{M} measurements (indexed by \eqn{j}), The model has the following form:
#' \deqn{f(y_{ij}|b_i)=\beta_0+\beta_1t_{ij}+\beta_2x_i+\beta_3t_{if}x_i+\beta_4c_i+
#' \beta_d(d_i)+\beta_t(d_i)t_{ij}+\beta_x(d_i)x_i+\beta_{tx}(d_i)t_{ij}x_i+\beta_c(d_i)c_i+
#' Z_ib_i
#' }
#' 
#' Where:
#' * \eqn{f()} s the logit transform
#' * \eqn{y_{ij}} is a 0/1 indicator of an outcome
#' * \eqn{b_i} is a vector of random coefficients, assumed to be normally distributed
#' * \eqn{Z_i} is the design matrix of random effects
#' * \eqn{\beta_k} are fixed model parameters
#' * \eqn{\beta_l(d_i)} are functions of dropout time (\eqn{d_i})
#' * \eqn{t_{ij}} is time since baseline for participant \eqn{i} at measure \eqn{j}
#' * \eqn{x_i} is the exposure status of participant \eqn{i}, (0/1)
#' * \eqn{c_i} is the status of the \eqn{i^{th}} participant for a confounding effect for exposure 
#'
#' Note that when an effect is independent of time the \eqn{j} subscript is suppressed, for clarity.
#'
#' When designing the dropout pattern, any fixed effect may be chosen to be associated with dropout time. Further, choice of association 
#' need not be the same for all effects.
#' 
#' Random coefficients can be defined for all fixed effects (in the random argument). Take care that the dimension of G
#' is appropriate for your selection. 
#' 
#' Last follow up visit is determined starting at the first follow up using a beta binomial distribution with \eqn{p\sim Beta(\alpha_1,\alpha_2)}
#' and \eqn{d_i\sim Binomial(M-1, p)+1}. Both study time and dropout time are divided by \eqn{M} so they are between zero and one.
#' 

#' @param B: numeric vector of fixed effect betas (intercept, time, exposure (0/1), time*exposure, confounder)
#' @param G: random variance matrix
#' @param random: character vector of random variables, if NULL just an intercept
#' @param p_exposed: scalar proportion with exposure
#' @param p_confounder: proportion of subjects with the confounding effect
#' @param conf: scalar confounder effect on exposure as proportion increase percent with exposure
#' @param N: number of subjects
#' @param M: number of measurement points after baseline
#' @param D_B: numeric vector of dropout betas (intercept, time, exposure (0/1), time*exposure, confounder)
#' @param D_F: numeric vector indicating dropout functions applied to each effect (intercept, time, exposure (0/1), time*exposure, confounder). 1 = linear, 2 = -exp()
#' @param D_parms: (shape, scale) for the beta distribution for dropout. To specify a different dropout pattern for exposed specify a 2 element list 
#' of parameters: list(c(1,1), c(3, 5)). Element 1 will be used for the unexposed, element 2 for the exposed.
#' 
 

#' @return Output:
#' * data frame of \eqn{N \cdot M} rows with the following columns:
#' * id: numeric subject ID
#' * exposure: (0/1) indicator of exposure status
#' * confounder: (0/1) indicator of confounding effect
#' * time: time since baseline
#' * time_exp: time by exposure interaction
#' * dtime: last time point with follow up
#' * y: outcome
#' 
#' By default, the complete data are included.To obtain the MNAR afflicted data: 
#' 
#' fullDat <- simuDat( your specifications here )
#' 
#' mnarDat <- subset(fullDat, time <= dtime)
#' 

 
#' @examples
#' set.seed(058475869)
#' fullData <- simuDat(
#'                 B = c(-2, -0.05, 1, 0.1, 1.2), 
#'                 G = matrix(c(1, 0.1, 0.1, 0.5), 2, 2),
#'                 random= 'time',
#'                 p_exposed = 0.5,
#'                 p_confounder = 0.25,
#'                 conf = 0.2,
#'                 D_B = c(0, -4, 0, 0, 0),
#'                 D_F = c(1, 1, 1, 1, 1),
#'                 D_parms = c(1, 1),
#'                 N = 1000,
#'                 M = 5
#'                 )
#' mnarData <- subset(fullData, time <= dtime)
#' 

#' @export
#' @md

simuDat <- function(
    
    B , # fixed effect betas (inercept, time, exposure (0/1), time*exposure, confounder)
    G , # random variance matrix
    random = NULL, # random variables, if NULL just an intercept
    p_exposed = 0.5, # proportion with exposure
    p_confounder = 0.25 , # proportion of subjects with the confounding effect
    conf = 0.2, # confounder effect on exposure as proportion increase percent with exposure
    N = 100 , # number of subjects
    M = 10,  # measurement points past baseline
    
    # dropout parameters
    D_B , # dropout betas (inercept, time, exposure (0/1), time*exposure, confounder)
    D_F = c(1, 1, 1, 1, 1), # dropout functions (intercept, time, exposure (0/1), time*exposure, confounder)
    D_parms = c(1, 1) # dropout parameters for the beta distribution
) {
  
  if(length(random)+1 != nrow(G)) {
    print('ERROR: missmatch between G and random arguments.')
  }
  # fixed data
  confounder <- rbinom(N, 1, p_confounder)
  
  exposed <- vapply(confounder, function(x){rbinom(1, 1, p_exposed + x*conf)}, integer(1))
  
  personLevel <- data.frame(id = 1:N,
                            exposure = exposed,
                            confounder = confounder)
  time <- data.frame (time = (0:M)/M)
  
  longData <- within(merge(personLevel, time), {
    time_exp <- time * exposure
  }) 
  
  longData <- longData[order(longData$id, longData$time), ]
  
  # fixed component (non-dropout)
  XB <- as.matrix(cbind(1, longData[, c('time','exposure','time_exp','confounder')])) %*% B
  
  # dropout coefficients
  if(!missing(D_B)) {
    DFuncs <- list(function(dtime, alpha){ alpha*dtime},
                   function(dtime, alpha){-exp(alpha*dtime)}
    )  
    
    if(class(D_parms) == 'numeric' & length(D_parms) == 2) {D_parms <- list(D_parms, D_parms)}
    else if(class(D_parms) == 'list' & length(D_parms) == 1) {D_parms[[2]] <- D_parms[[1]]}
    else if(class(D_parms) != 'list' | !(length(D_parms) %in% c(1,2))) {
      print('ERROR: D_parms is miss-specified')
      stop()
    } 
    
    
    
    u <- vapply(personLevel$exposure, function(x){(rbinom(1,M-1,rbeta(1, D_parms[[x+1]][1], D_parms[[x+1]][2]))+1)/M}, numeric(1))
    idx <- rep(1:length(u), times = M + 1)
    idx <- sort(idx)

    
    DB <- vapply(u, 
                 function(x) { 
                   c(
                     DFuncs[[D_F[1]]](x, D_B[1]),
                     DFuncs[[D_F[2]]](x, D_B[2]),
                     DFuncs[[D_F[3]]](x, D_B[3]),
                     DFuncs[[D_F[4]]](x, D_B[4]),
                     DFuncs[[D_F[5]]](x, D_B[5])
                   )
                   }, 
                 numeric(length(D_B)))
    DB <- t(DB)
    DB <- (DB[idx,] * as.matrix(cbind(1,longData[ , c('time','exposure','time_exp','confounder')]))) %*% rep(1, length(D_B))
    
    #u <- u[idx]
    
    longData$dtime <- u[idx]
    
  } else DB <- 0

    #  random coefficients
  if(!missing(G)) {
    rcoef <- MASS::mvrnorm(N, mu=rep(0,dim(G)[1]), Sigma = G)
    idx <- rep(1:nrow(rcoef), times = M + 1)
    idx <- sort(idx)
    rcoef <- rcoef[idx,]
    
    Zb <- if(is.null(random)) rcoef else (as.matrix(cbind(1, longData[, random])) * rcoef) %*% rep(1, nrow(G))
    
  } else Zb <- 0
  
  
  # calculate mean and generate random outcome
  eta <- XB + DB + Zb
  
  p <- 1/(1+exp(-eta))
  
  longData$y <- vapply(p, function(x){rbinom(1,1,x)}, integer(1))
  
  row.names(longData) <- 1:nrow(longData)

  return(longData)
}

