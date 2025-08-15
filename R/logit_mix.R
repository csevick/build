
#' Accepts a covariance vector and outputs covariance matrix
#'
#' @description: This function accepts a covariance vector and outputs covariance matrix.
#' The assumed structure can be diagonal or unstructured and can compute the matrix assuming that the 
#' parameter vector represents the Cholesky root of the final product.  
#' 
#' @param cvec: a vector of covariance parameters
#' @param cstruct: character - desired matrix structure ('diag'/'un')
#' @param chol: boolean. If true then cvec is considered to be the Cholesky root of the final product matrix.
#'
#' @return A matrix
#' @export
#' 
#' @examples
#' a <- c(0.9,1.2)
#' gmat(a, 'diag', chol = F)
#' gmat(a, 'diag', chol = T)
#' 
#' b <- c(2, 0.25, 0.5)
#' gmat(b, 'diag', chol = F)
#' gmat(b, 'un', chol = F)
#' gmat(b, 'un', chol = T)

gmat <- function(cvec, cstruct = 'un', chol=F) {
  if(tolower(cstruct) == 'un') {
    
    nrowG <- 1
    while (nrowG^2<length(cvec)) {
      nrowG <- nrowG+1
    }
    L <- matrix(0, nrowG, nrowG)
    L[lower.tri(L, diag = TRUE)] <- cvec
    G <- if(chol==T) {
      L%*%t(L)
    } else {
      M <- L
      M[upper.tri(L)] <- t(L[lower.tri(L)])
      M
    }
    
  } else if(tolower(cstruct) == 'diag') {
    G <- diag(length(cvec))*cvec
    if(chol==T) {G <- G^2}
    
  } else {
    stop(print('Unknown covariance structure selected, use un or diag'))
  }
  G
}

#' Compute the log-likelihood for a binary logistic model for a single cluster.
#'
#' @description: Compute the log-likelihood for a binary logistic model for a single cluster. Option to add ascertainment correction.
#' 
#' @param parms: a vector of parameters, fixed effects followed by covariance parameters
#' @param Xi: a matrix of observed fixed data for the cluster
#' @param Yi: a vector outcomes for the cluster
#' @param Zi: a matrix of random effects for the cluster
#' @param cstruct: covariance structure for the random effects
#' @param ghqpw: a list of points and weights to integrate the likelihood
#' @param acml: boolean - if T then use ascertainment correction
#' @param SP: matrix of sampling probabilities 3XD - D is the number of dropout categories
#' @param OS: Outcome strata of the cluster (1-3)
#' @param DS: integer dropout strata of the cluster (1-D).
#'
#' @return scalar
#' @export
#' 
#' @examples
#' ll_i_logitGlmm(parms, Xi, Yi, Zi, cstruct, ghqpw, acml = F, SP, OS, DS)

ll_i_logitGlmm <- function(parms, Xi, Yi, Zi, cstruct, ghqpw, acml = F
                           , SP = NULL # matrix of sampling probabilities 3XD
                           , OS = NULL # Outcome strata
                           , DS = NULL # dropout strata
){
  # get number of parameters
  nx <- ncol(Xi)
  nz <- ncol(Zi)
  
  # unpack parameters
  XB <- as.vector(Xi %*% parms[1:nx])
  G <- gmat(parms[(nx+1):length(parms)], cstruct = cstruct, chol=T)
  
  # adjust quadrature points for the covariance matrix
  ghqT <- gquad_cntr_scale(pw = ghqpw, vcov = G)
  
  #XB <- Xi %*% B
  #L <- 0
  #for(k in 1:length(ghqT$wts)) {
  #eta <- XB + Zi %*% ghqT$pnts[k,]
  #    p   <- 1/(1+exp(-(XB + Zi %*% ghqT$pnts[k,])))
  #    L <- L + ghqT$wts[k]*prod( (p^Yi)*((1-p)^(1-Yi)) )
  #}
  
  p   <- 1/(1+exp(-(XB + Zi %*% ghqT$pnts)))
  
  if (acml == F){
    L <- sum(ghqT$wts*exp(colSums(log( (p^Yi)*((1-p)^(1-Yi))))))
  }
  # Ascertainment correction, if applicable
  else if(acml == T) {
    # the never strata
    L1 <- L_i_logit_glmm(Yi = rep(0, length(Yi))
                         ,Pi = p
                         ,ghqpw = ghqT)
    # the always strata
    L3 <- L_i_logit_glmm(Yi = rep(1, length(Yi))
                         ,Pi = p
                         ,ghqpw = ghqT)
    # the sometimes strata
    L2 <- 1 - L1 - L3
    
    L <- SP[OS, DS]/(SP[1, DS]*L1 + SP[2, DS]*L2 + SP[3, DS]*L3)
  }
  else (
    stop('The value of the acml arg must be T/F')
  )
  
  log(L)
}

#' Compute the likelihood for a binary logistic model for a single cluster.
#'
#' @description: Compute the likelihood for a binary logistic model for a single cluster.
#' 
#' @param Yi: a vector outcomes for the cluster
#' @param Pi: a vector probabilities for each measurement
#' @param ghqpw: a list of points and weights to integrate the likelihood
#'
#' @return scalar
#' @export
#' 
#' @examples
#' L_i_logit_glmm(Yi, Pi, ghqpw)

L_i_logit_glmm <- function(Yi, Pi, ghqpw){
  
  # compute likelihood
  
  L <- ghqpw$wts %*% exp(colSums(log( (Pi^Yi)*((1-Pi)^(1-Yi)))))
  
  return(L)
}

#' Cluster -2*log-likelihood contribution for specified random parameters.
#'
#' @description: Cluster -2*log-likelihood contribution for specified random parameters.
#' 
#' @param re_i: random parameters for cluster i
#' @param G: covariance parameters of the random effects
#' @param XB: vector of fixed effect variables multiplied by corresponding vector of parameters
#' @param Yi: a vector outcomes for the cluster
#' @param Zi: a matrix of random effects for the cluster
#' 
#' @return scalar
#' @export
#' 
#' @examples
#' n2ll_rparm_i_logitGlmm(re_i, G, XB, Yi, Zi)

n2ll_rparm_i_logitGlmm <- function(re_i, G, XB, Yi, Zi){

  # compute log-likelihood
  Pi   <- 1/(1+exp(-(XB + Zi %*% re_i)))
  
  ll <- Yi %*% log(Pi) + (1-Yi) %*% log(1-Pi) - 0.5*(t(re_i) %*% solve(G) %*% re_i) - 0.5*log(det(G)) - (length(re_i))/2*log(2*pi)
  #L <- (exp(colSums(log( (Pi^Yi)*((1-Pi)^(1-Yi))))))*dmvnorm(re_i, sigma = G)
  as.numeric(-2*ll)
  #-2*log(as.numeric(L))
}


#' Compute the -2log-likelihood for a binary logistic model for a full dataset.
#'
#' @description: Compute the -2log-likelihood for a binary logistic model for a full dataset.
#' 
#' @param parms: a vector of parameters, fixed effects followed by covariance parameters
#' @param X: list of matrices of observed fixed data for the cluster
#' @param Y: a list of a vectors of outcomes
#' @param Z: a list of matrices of random effects 
#' @param cstruct: covariance structure for the random effects
#' @param ghqpw: a list of points and weights to integrate the likelihood
#' @param acml: boolean - if T then use ascertainment correction
#' @param SP: matrix of sampling probabilities 3XD - D is the number of dropout categories
#' @param OS: a vector of outcome strata of the cluster (1-3)
#' @param DS: a vector of integer dropout strata of the cluster (1-D).
#'
#' @return scalar
#' @export
#' 
#' @examples
#' n2ll_logit_glmm(parms, X, Y, Z, SW, cstruct, ghqpw, acml = F
#' , SP = NULL # matrix of sampling probabilities 3XD
#' , OS = NULL # Outcome strata
#' , DS = NULL # dropout strata
#' , ll_only = T
#' , adapt = F
#' )

n2ll_logit_glmm <- function(parms, X, Y, Z, SW, cstruct, ghqpw, acml = F
                           , SP = NULL # matrix of sampling probabilities 3XD
                           , OS = NULL # Outcome strata
                           , DS = NULL # dropout strata
                           , ll_only = T
                           , adapt = F
) {
  
  # get number of parameters
  nx <- ncol(X[[1]])
  nz <- ncol(Z[[1]])
  
  # unpack parameters
  B <- parms[1:nx]
  G <- gmat(parms[(nx+1):length(parms)], cstruct = cstruct, chol=T)
  
  # adjust quadrature points for the covariance matrix
  if(adapt == F) {ghqT <- gquad_cntr_scale(pw = ghqpw, vcov = G)}

  L <- rep(0, length(Z))
  if(acml == T) {AC <- rep(0, length(Z))} # Ascertainment correction, if applicable
  for(i in 1:length(Z)) {
    
    XB <- as.vector(X[[i]] %*% B)

    if(adapt == T) {
      # estimate random effects;
      re_opt <- optim(rep(0,nz)
                      ,n2ll_rparm_i_logitGlmm
                      ,G = G
                      ,XB = XB
                      ,Yi = Y[[i]]
                      ,Zi = Z[[i]]
                      ,hessian = T
                      ,method = "BFGS" #"Nelder-Mead"
      )
      ghqT <- list()
      
      V <- solve(re_opt$hessian/2)
      #ghqT <- gquad_cntr_scale(pw = ghqpw, mu = re_opt$par, vcov = V)
      sqrtV <- sqrt(V) #mat_spec_sqrt(V)
      ghqT$pnts <- re_opt$par + sqrtV %*% ghqpw$pnts
      ghqT$wts <- ((2*pi)^(nz/2))*sqrt(det(V)) * 
        apply(exp(( ghqpw$pnts^2)/2), 2, prod) * 
        apply(ghqT$pnts, 2, function(x){dmvnorm(x, sigma = G)}) *
        ghqpw$wts
        
    }
    
    Pi <- 1/(1+exp(-(XB + Z[[i]] %*% ghqT$pnts)))
    
    L[i] <- L_i_logit_glmm(Yi = Y[[i]]
                          ,Pi = Pi
                          ,ghqpw = ghqT)
    
    # Ascertainment correction, if applicable
    if(acml == T) {
      # the never strata
      L1 <- L_i_logit_glmm(Yi = rep(0, length(Y[[i]]))
                           ,Pi = Pi
                           ,ghqpw = ghqT)
      # the always strata
      L3 <- L_i_logit_glmm(Yi = rep(1, length(Y[[i]]))
                           ,Pi = Pi
                           ,ghqpw = ghqT)
      # the sometimes strata
      L2 <- 1 - L1 - L3
      
      AC[i] <- SP[OS[i], DS[i]]/(SP[1, DS[i]]*L1 + SP[2, DS[i]]*L2 + SP[3, DS[i]]*L3)
    }
  }
  
  ll <- ifelse(acml == F
               ,-2*(SW %*% log(L))
               ,-2*sum(log(AC * L))
               )
  
  if(ll_only == T) {
    return(ll)
  } else {
    list(n2ll = ll, likelihood_raw = L, acml_adj = AC)
  }
  
}



#' GLMM for binary outcomes, using cluster sample weighting and ascertainment correction 
#'
#' @description: GLMM for binary outcomes, using cluster sample weighting and ascertainment correction. Robust variance estimates of the fixed effects are available. Derivatives of the log-likelihood are computed using the numDeriv package.
#'  
#' 
#' @param df: a data frame
#' @param fixed: formula for the fixed effect portion of the model
#' @param random: formula for the fixed effect portion of the model 
#' @param cluster: character - name of the variable in df that defines a cluster
#' @param sweights: character - name of the variable in df that holds sampling weights (cluster level)
#' @param acml: boolean - if true than use ascertainment corrected likelihood
#' @param sprobs: a matrix of sampling probabilities dimension 3XD - D being the number of dropout categories
#' @param outcome_strata: character - name of the variable in df that contains the cluster level outcome strata
#' @param drop_strata: character - name of the variable in df that contains the cluster level dropout strata
#' @param cstruct: text (diag, un) - covariance structure for the random effects
#' @param robust: boolean - if true then compute robust variances for the fixed effects
#' @param rpred: boolean - if true add a prediction of the random cluster random effects to the output object
#' @param qpoints: numeric scalar - number of quadrature points to use when integrating the likelihood
#' @param starting: numeric vector - starting values fixed followed by random. choice of cstruct and random arguments determine form
#' @param opt_method: character - optimization method (see optim function)
#' @param opt_control: list - arguments to fine tune the optimization (see optim function)
#' @param adapt: boolean - if true then use adaptive quadrature - this is in development, slow and does not converge well 
#'  
#' @return A list of regression products
#' func_call: the function call that created the model object
#' data: a copy of the dataframe used in the analysis, sorted the cluster ID
#' opt_fit: the output from the optim function for the model optimization
#' vcov: the variance covariance matrix of the fixed effects and parameters of the variance components (Cholesky decomposition)
#' G: the estimated covariance matrix of the random effects
#' rvcov: (if robust = T) the robust covariance matrix of the fixed effects
#' fixed: a data frame containing the estimates (estimate) standard error (std_err) and test statistic (wald_stat) of the fixed effect parameter estimates
#" repred: (if repred = T) a matrix of predicted random effects. The number of rows are equal to the number of unique clusters and the columns to the columns in G 

#' @references
#' Gilbert, P. and Varadhan, R. (2025). numDeriv: Accurate Numerical Derivatives. R package version 2022.9-1.

#' @export
#' 
#' @examples
#' set.seed(957385)
#' dat = simuDat(
#'   B = c(-3, 2, 6, 2, 0),
#'   G = matrix(c(4), ncol = 1),
#'   random = NULL,
#'   p_exposed = 0.5,
#'   p_confounder = 0.25,
#'   conf = 0,
#'   N = 1000,
#'   M = 5,
#'   D_B = c(0, -2, -1, -0.5, 0),
#'   D_F = c(1, 1, 1, 1, 1),
#'   D_parms = c(1, 1)
#' ) 
#' mdat <- subset(dat, time <= dtime)
#' 
#' l_mix_fit <- logit_mix(df=mdat, fixed=y~time, random = ~ 1, cluster='id')

logit_mix <- function(df, fixed, random = ~ 1, cluster
                     , sweights = NULL
                     , acml= F, sprobs = NULL, outcome_strata = NULL, drop_strata = NULL
                     , cstruct='un', robust = F, repred = F, qpoints = 10, starting = NULL
                     , opt_method = "BFGS", opt_control = list()
                     , adapt = F) {
  
  # create a list to hold output
  output <- list()
  
  # store the function call meta data
  output$func_call <- match.call()
  
  # sort data by cluster
  df <- df[order(df[,cluster, drop=T]), ]
  
  # add data to output object 
  output$data <- df		
  
  # create data components
  CID <- as.factor(df[ , cluster, drop = T])
  
  X <- model.matrix(fixed, data = df)
  nx <- ncol(X)
  Z <- model.matrix(random, data = df)
  nz <- ncol(Z)
  Y <- df[ , as.character(fixed[[2]]), drop = T]
  
  # test for missingness
  if(sum(is.na(CID)) > 0) {
    print('There are missing cluster IDs. Please fix your missing issues before continuing.')
    stop()
  }
  if(sum(is.na(X)) > 0) {
    print('There are missing data values among fixed predictors. Please fix your missing issues before continuing.')
    stop()
  }
  if(sum(is.na(Z)) > 0) {
    print('There are missing data values among random predictors. Please fix your missing issues before continuing.')
    stop()
  }
  if(sum(is.na(Y)) > 0) {
    print('There are missing outcome values. Please fix your missing issues before continuing.')
    stop()
  }

  # break data elements in to lists by cluster
  X <- split.data.frame(X, CID) #by(X, CID, function(x){x})
  nclus <- length(X)
  Z <- split.data.frame(Z, CID)
  Y <- split(Y, CID)
  SW <-  by(df, df[,cluster], function(x){
    if(is.null(sweights)) {
      1
    } else {
      tmp <- unique(x[,c(cluster,sweights)])[,sweights, drop=T]
      if(length(tmp) != 1){
        stop("There must be one selection weight per cluster")
      } else {
        tmp
      }
    }
  })
  if(acml == T) {
  OS <- by(df, df[,cluster], function(x){
    if(is.null(outcome_strata)) {
      stop("ERROR: you must specify an outcome strata for ACML")
    } else {
      tmp <- unique(x[,c(cluster,outcome_strata)])[,outcome_strata, drop=T]
      if(length(tmp) != 1){
        stop("There must be one outcome stratum per cluster")
      } else {
        tmp
      }
    }
  })
  DS <-  by(df, df[,cluster], function(x){
    if(is.null(drop_strata)) {
      stop("ERROR: you must specify a dropout strata for ACML")
    } else {
      tmp <- unique(x[,c(cluster,drop_strata)])[,drop_strata, drop=T]
      if(length(tmp) != 1){
        stop("There must be one dropout stratum per cluster")
      } else {
        tmp
      }
    }
  })
  }

  # points and weights for Gaussian quadrature
  ghqpw <- gauss_herm_consts(qpoints, nz)
  
  # Initial parameter estimates
  if(is.null(starting)) {
    # X starting values
    #my <- mean(df[,Y,drop = T])
    xstart <- stats::glm(fixed, data = df, family = binomial())$coefficients#c(log(my/(1-my)), rep(0, length(X)))
    
    zstart <- if(cstruct == 'un' & nz>1) {
      tmp <- diag(nz)#matrix(0, nz, nz)
      as.vector(tmp[lower.tri(tmp, diag = TRUE)])
    } else {
      rep(1, nz)
    }
    starting <- c(xstart, zstart)
  } 
  nparms <- length(starting)
  
  output$opt_fit <- optim(starting ,
                   n2ll_logit_glmm,
                  X=X, Z=Z, Y=Y,SW = SW
                  , acml = acml
                  , SP = sprobs # matrix of sampling probabilities 3XD
                  , OS = OS # Outcome strata
                  , DS = DS # dropout strata
                  ,cstruct = cstruct,
                  ghqpw = ghqpw,
                  adapt = adapt,
                  hessian=T, method=opt_method,
                  control = opt_control)
  
  output$vcov  <- solve(output$opt_fit$hessian/2) 
  
  output$G  <- gmat(output$opt_fit$par[(nx+1):nparms], cstruct=cstruct, chol=T)
  row.names(output$G) <- c(colnames(Z[[1]]))
  colnames(output$G) <- c(colnames(Z[[1]]))
  
  # robust variance estimate
  #if(robust == T & acml == T) {print('WARNING: Robust variances are not supported, or computed, with ACML.')}
  if(robust == T #& acml == F
     ) {
    
    Ji <- matrix(0, nparms,nparms)
    for(i in 1:nclus) {
      
      
      SBhat <- SW[i]*numDeriv::grad(ll_i_logitGlmm
                                    , x=output$opt_fit$par 
                                    , Xi=X[[i]]
                                    , Yi=Y[[i]]
                                    , Zi=Z[[i]]
                                    , SP = sprobs # matrix of sampling probabilities 3XD
                                    , OS = OS[i] # Outcome strata
                                    , DS = DS[i] # dropout strata
                                    
                                    , cstruct = cstruct
                                    , ghqpw = ghqpw
                                    , acml = acml)#$D[1,1:nParms]
      
      
      Ji[1:nparms,1:nparms] <- Ji + SBhat%*%t(SBhat) # Rabe-Hesketh (2006)
    }
    output$rvcov <-  output$vcov%*%((nclus/(nclus-1))*Ji)%*%output$vcov
  } 
  
  ## summary tables
  fixed <- data.frame(estimate=output$opt_fit$par[1:nx],
                      std_err = sqrt(diag(
                        if(!is.null(output$rvcov)) {output$rvcov[1:nx, 1:nx]} else {output$vcov[1:nx, 1:nx]}
                      ))
  )
  fixed$wald_stat <- fixed$estimate/fixed$std_err
  row.names(fixed) <- c(colnames(X[[1]]))
  output$fixed <- fixed
  
  # random effect predictions
  #(re_i, G, XB, Yi, Zi)
  if(repred == T) {
    repred. <- list()
    for(i in 1:nclus) {
      tmp <- optim(rep(0,nz)
                   ,n2ll_rparm_i_logitGlmm
                   ,method = 'BFGS' #ifelse(nz == 1 , "Brent", "Nelder-Mead")
                   #,lower = ifelse(nz == 1, -10, -Inf)
                   #,upper = ifelse(nz == 1,  10,  Inf)
                   ,G = output$G
                   ,XB=X[[i]] %*% output$fixed$estimate
                   ,Yi=Y[[i]]
                   ,Zi=Z[[i]]
      )
      repred.[[i]] <- if(tmp$convergence == 0) {tmp$par} else {NULL}
    }
    output$repred <- do.call('rbind', repred.)
  }
  

  return(output)
}
