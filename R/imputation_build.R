# likelihood for a logit model

#' likelihood for a logit model, Independence 
#'
#' @description: Returns the likelihood for a single observation.
#'
#' @param Yi: binary outcome (1/0) for the ith observation 
#' @param Pi: probability of the outcome
#'
#' @return A scalar 
#'
#' @export
#' 
#' @examples
#' L_i_logit(1, 0.25)
#' L_i_logit(0, 0.25)

L_i_logit <- function(Yi, Pi){
  
  # compute likelihood
  
  L <-  (Pi^Yi)*((1-Pi)^(1-Yi))
  
  return(L)
}


#' Compute exposure probabilities conditional on outcomes 
#'
#' @description: A function to compute exposure probabilities conditional on outcomes.
#'
#' @param y_formula: a model formula for the fixed effect portion of the outcomes model (see buildr::logit_mix)
#' @param z_formula: a model formula for the random effect portion of the outcomes model (see buildr::logit_mix)
#' @param y_fit: output object from the outcomes model (see buildr::logit_mix) 
#' @param e_formula: a model formula for the marginal exposure  (see the stat::gml function)
#' @param e_fit: an output object from the marginal exposure model (see the stat::gml function)
#' @param full_data: full population data (with undetermined exposure)
#' @param exposure: character - name of the exposure variable in the sampled data
#' @param study_id: character - name of the cluster id variable in the population and sampled data
#' @param outcome: character - name of the study outcome variable in the population and sampled data
#' @param time: character - name of the time variable in the population and sampled data
#' @param qpoints: integer - number of quadrature points needed to integrate the likelihood
#'
#' @return A vector of probabilities of exposure, in participant order
#'
#' @export
#' 

cond_p_exposure <- function(y_formula
                            , z_formula = ~1
                            , y_fit
                            , e_formula
                            , e_fit
                            , full_data  
                            , exposure = 'exposure'
                            , study_id = 'id'
                            , outcome = 'y'
                            ,time = 'time'
                            , qpoints
) {
  
    

  # process sample data 
  CID <- as.factor(full_data[ , study_id, drop = T])
  
  X_e <- model.matrix(y_formula  , data = within(full_data, {exposure = 1}))
  X_u <- model.matrix(y_formula  , data = within(full_data, {exposure = 0}))
  Z <- model.matrix(z_formula, data = full_data)
  Y <- full_data[ , outcome, drop = T]
 X_p <- model.matrix(e_formula  , data = within(full_data[full_data[, time, drop = T]== 0, ] , {exposure = 0}))
  X_e <- split.data.frame(X_e, CID) #by(X, CID, function(x){x})
  X_u <- split.data.frame(X_u, CID) #by(X, CID, function(x){x})
  X_p <- split.data.frame(X_p, as.factor(unique(full_data$id))) #by(X, CID, function(x){x})
  Z <- split.data.frame(Z, CID)
  Y <- split(Y, CID)
  
  # unpack parameters
  B <- y_fit$fixed$estimate #  c(parms[1:nx], Be)
  G <- y_fit$G  #gmat(parms[(nx+1):length(parms)], cstruct = cstruct, chol=T)
  #print(G)
  # adjust quadrature points for the covariance matrix
  
  # points and weights for Gaussian quadrature
  ghqpw <- gauss_herm_consts(qpoints, ncol(Z[[1]]))
  ghqT <- gquad_cntr_scale(pw = ghqpw, vcov = y_fit$G)
  
  Le <- rep(0, length(Z))
  Lu <- rep(0, length(Z))
  Lpe <- rep(0, length(Z))
  Lpu <- rep(0, length(Z))
  cpe <- rep(0, length(Z))
  for(i in 1:length(Z)) {
    
    # exposure P 
    Pi_p <- as.numeric(1/(1+exp(-( X_p[[i]] %*% e_fit$coefficients ))))
    # exposure likelihood
    Lpe[i] <- L_i_logit(1, Pi_p )
    # non-exposure likelihood
    Lpu[i] <- L_i_logit(0, Pi_p )
    
    
    # Likelihood as exposed
    XB <- as.vector(X_e[[i]] %*% B)
    
    Pi <- 1/(1+exp(-(XB + Z[[i]] %*% ghqT$pnts)))
    
    Le[i] <- L_i_logit_glmm(Yi = Y[[i]]
                            ,Pi = Pi
                            ,ghqpw = ghqT)
    
    #likelihood as unexposed
    XB <- as.vector(X_u[[i]] %*% B)
    
    Pi <- 1/(1+exp(-(XB + Z[[i]] %*% ghqT$pnts)))
    
    Lu[i] <- L_i_logit_glmm(Yi = Y[[i]]
                            ,Pi = Pi
                            ,ghqpw = ghqT)
    
    #conditional probability of exposure
    cpe[i] <- Lpe[i]*Le[i] / (Lpe[i]*Le[i] + Lpu[i]*Lu[i])
  }
  
  return(cpe)
}


#' Imputation of a baseline binary exposure 
#' 
#' @description Imputation of a baseline binary exposure. Produces a list of imputed dataframes on the dataset named in the 'full_data' argument.
#'
#' @param y_formula: a model formula for the fixed effect portion of the outcomes model (see buildr::logit_mix)
#' @param e_formula: a model formula for the marginal exposure  (see the stat::gml function)
#' @param sample_subj: subject level dataset, from the sampled data (includes assessed exposure)
#' @param full_data: population data that was sampled from 
#' @param strata: character vector - names of the outcome and dropout variables (in that order)
#' @param qpoints: number of quadrature points to integrate the likelihood
#' @param study_id: character scalar - cluster identifier
#' @param z_formula: a model formula for the random effect portion of the outcomes model (see buildr::logit_mix)
#' @param exposure: character scalar - name of the binary exposure variable
#' @param outcome: character scalar - name of the outcome variable
#' @param time: character - name of the time effect
#' @param sweights character scalar - name of the variable holding sample weight data
#' @param n_imp: number of imputed datasets to produce
#' 
#' @return A list of imputed dataframes
#' @export


  
binary_boot_impute <- function(y_formula, e_formula, sample_subj, full_data, strata=c('oc_strata', 'dtime'), qpoints,
                               study_id = 'id', z_formula, exposure = 'exposure', outcome, time = 'time', sweights='sampWT', n_imp = 10) {
  
  # adapted from: "Flexible Imputation of Missing Data" by Stef Van Buuren

  if(exposure %in% names(full_data)){
    stop(paste('ERROR: ', exposure, ' cannot be a variable in the full data'))
  }
  
  add_strata <- setdiff(strata, names(full_data))
  data_C <- merge(sample_subj[, c(study_id, exposure, sweights, add_strata)], full_data, by = study_id, all.x = F, all.y = F)
  
  ires <- list()
  for(i in 1:n_imp) {
    # draw boot sample
    B <- strata_clus_boot(data_C, strata, study_id)
 
    # estimate outcome model in sample 
    y_fit <- logit_mix(B, fixed =  y_formula, random = ~ 1, cluster='.ClusID.'
                       , sweights = sweights,  qpoints = qpoints)
    evalit <<- y_fit
    # survey design object
    D <- survey::svydesign(ids=~1, strata = as.formula(paste('~', paste(strata, collapse = '+'), collapse = ' ')) 
                           ,weights = as.formula(paste0('~',sweights)), data = subset(B, time == 0 ))
 
    # estimate exposure model in complete data
    e_fit <- survey::svyglm(e_formula, design = D, family = quasibinomial())
    #fit_C <- stats::glm(formula, data = B, family = binomial())
    
   # conditional probability
    .P. <- cond_p_exposure(y_formula = y_formula, y_fit = y_fit, e_formula = e_formula, e_fit = e_fit
      , z_formula = z_formula
      , full_data = full_data # Full data, with undetermined exposure
      , exposure = exposure
      , study_id = study_id
      , outcome = outcome
      , time = time
      , qpoints = qpoints
    )  

    # imputation values
    #.U. <- runif(nrow(data_F))
    data_F <- merge(sample_subj[, c(study_id, exposure)], full_data, by = study_id, all.x = F, all.y = T)
    .imp1. <- sapply(.P., function(x){rbinom(1,1,x)}) #ifelse(.U. <= .P., 1, 0)

    # impute 
    tmp <- data_F[data_F[,time, drop = T] == 0, exposure] 
    .imp. <- rep(0, length(tmp))
    for(j in 1:length(tmp)) {
      .imp.[j] <- ifelse(is.na(tmp[j]) == T, .imp1.[j], tmp[j])
    }
    ires[[i]] <- cbind(subset(data_F, time == 0, c(id)), .imp. = .imp.)
    
  }
  
  return(ires)
}


#' Function to combine results from models on imputed data

#' @description 
#' A function to combine results from models on imputed data.
#' Returns a list, the first element is the mean vector, while element 2 
#' is the variance matrix

#' @param mu: a mXp matrix of parameter estimates (p is the number of parameters and m is the number of imputations)
#' @param vcov: a pXpXm tensor of covariance matrices

#' @return A two element list:
#' 1) a vector of the imputation adjusted estimates
#' 2) a variance covariance matrix for the estimates 
#' 
#' @export

combine_imp <- function(mu # a mXp matrix of parameter estimates (p is the number of parameters and m is the number of imputations)
                        , vcov # a pXpXm tensor of covariance matrices
){
  
  # Rubin's rules: "Flexible Imputation of Missing Data" by Stef Van Buuren
  
  mu_bar <- colMeans(mu)
  
  V_bar <- apply(vcov, 1:2, mean) + (1 + 1/nrow(mu)) * var(mu)
  
  list(mu_bar, V_bar)
  
}


#' A convenience function to run the imputation, analyze the data and combine results 
#' 
#' @description A convenience function to run the imputation, analyze the data and combine results
#'
#' @param y_formula: a model formula for the fixed effect portion of the outcomes model (see buildr::logit_mix)
#' @param e_formula: a model formula for the marginal exposure  (see the stat::gml function)
#' @param sample_subj: subject level dataset, from the sampled data (includes assessed exposure)
#' @param full_data: population data that was sampled from 
#' @param strata: character vector - names of the outcome and dropout variables (in that order)
#' @param qpoints: number of quadrature points to integrate the likelihood
#' @param study_id: character scalar - cluster identifier
#' @param z_formula: a model formula for the random effect portion of the outcomes model (see buildr::logit_mix)
#' @param exposure: character scalar - name of the binary exposure variable
#' @param outcome: character scalar - name of the outcome variable
#' @param time: character - name of the time effect (scaled such that 0 )
#' @param sweights: character scalar - name of the variable holding sample weight data
#' @param n_imp: number of imputed datasets to produce
#' 
#' @return A list. The first element contains a vector of dropout corrected regression coefficients. Element 2 is a matrix of regression coefficients, each row from a separate imputation. 
#' The third element is the corrected variance matrix for the coefficients in element 1. Element 4 is an array with three dimensions containing the variance matrix from each imputation. The third dimension indexes the imputation.
#' The final element is a list of the fit object from each imputation. 
#' 
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
#' s_frame <- build_sample_frame(data = mdat, study_id = 'id', outcome = 'y', drop_strata = 'dtime')
#' 
#' s_plan <- build_sample_plan(n_samp = c(50, 200, 50), N_bystrata = s_frame$strata_counts, drop_grad = 1)
#' 
#' s_sample <- strat_samp(s_frame$sample_frame,oc_strata = 'oc_strata', drop_strata = 'dtime', sample_plan=s_plan$sample_plan,sample_by = 'n')
#' 
#' # get exposure status
#'  sample <- merge(mdat[mdat$time ==0, ]
#'               , s_sample[, c('id'
#'                          , 'oc_strata'
#'                          , 'drop_strata'
#'                          , 'sampWT'
#'                         )]
#'               , by = 'id')
#' # MI analysis
#'  im_fit <- imputeit(y_formula = y ~ time + exposure + time:exposure 
#'                   , e_formula = exposure ~ 1
#'                   ,  strata='oc_strata'
#'                   , sample_subj  = sample
#'                   , full_data = mdat 
#'                   , study_id = 'id', sweights='sampWT', z_formula = ~1, exposure = 'exposure'
#'                   , outcome='y', time = 'time' 
#'                   , n_imp = 5
#'                   , qpoints = 20)



imputeit <- function(y_formula, e_formula
                     , sample_subj, full_data, strata 
                     , study_id , sweights, z_formula, exposure
                     , outcome, time 
                     , n_imp = 30
                     , qpoints = 10
                     , aparms =c('(Intercept)','time','exposure','time:exposure') # parameters with hypothesized dropout effects, needing adjustment by dparms
                     , dparms =c('dtime','time:dtime','exposure:dtime','time:exposure:dtime')# dropout parameters to adjust aparms 
                     , dtime = 'dtime' #name of the variable containing subject dropout time
                     , where = c(0,0,1,1) # numeric vector of dimension aparms, indicating which subgroup of exposure the dropout time proportion should be calculated
) {
  
  imp_data <- binary_boot_impute(y_formula=y_formula, e_formula=e_formula
                                 , sample_subj=sample_subj, full_data=full_data, strata = strata
                                 , n_imp = n_imp, study_id = study_id, sweights=sweights, z_formula=z_formula, exposure = exposure
                                 , outcome, time = time, qpoints = qpoints
                                 

                                 )

    # models run on imputed data 
    mres <- list()
    for(i in 1:length(imp_data)) {
      
      tmp_imp_dat <- imp_data[[i]][,c(study_id,'.imp.')]
      names(tmp_imp_dat) <- c(study_id, exposure)
      
      idata <- merge(full_data, tmp_imp_dat, by = study_id)

      tmp_fit <- try(logit_mix(idata, fixed =  y_formula, random = z_formula, cluster=study_id, sweights = NULL, qpoints = qpoints)
      )
      if(!inherits(tmp_fit, "try-error")) {
        tmpdcorr <- drop_adj_all(tmp_fit
                                                , aparms = aparms
                                                , dparms = dparms  
                                                , time = time 
                                                , exposure = exposure 
                                                , dtime = dtime
                                                , where = where )
        tmp_fit$dropout_correct <- tmpdcorr
        
        ## summary table
        fixed_dropcorrect <- data.frame(estimate=tmpdcorr$estimate,
                            std_err = sqrt(diag(tmpdcorr$vcov))
        
        )
        fixed_dropcorrect$wald_stat <- fixed_dropcorrect$estimate/fixed_dropcorrect$std_err
        row.names(fixed_dropcorrect) <- names(tmpdcorr$estimate)
        tmp_fit$fixed_dropcorrect <- fixed_dropcorrect
        
      }
      mres[[i]] <- tmp_fit  
    }

  mdim <- length(mres[[1]]$fixed_dropcorrect$estimate)
  mu <- matrix(NA, length(mres), mdim)
  colnames(mu) <- rownames(mres[[1]]$fixed_dropcorrect)
  vcov <- array(NA, c(mdim,mdim,length(mres)))
  for(i in 1:length(mres)) {
    mu[i,] <- mres[[i]]$fixed_dropcorrect$estimate
    vcov[,,i] <- mres[[i]]$dropout_correct$vcov[1:mdim,1:mdim]
  }
  
  iparms <- combine_imp(mu, vcov)
  
  return(list(mu_bar = iparms[[1]], MUs = mu, V_bar = iparms[[2]], VCOVs = vcov, fits = mres))
}

