
#' Function to apply the dropout correction to multiple parameter estimates

#' @description 
#' This function is only validated for simple conditional linear models

#' @param model_fit: a fit object from logit_mix
#' @param aparms: character vector - parameters with hypothesized dropout effects, needing adjustment by dparms
#' @param dparms: character vector - dropout parameters to adjust aparms
#' @param time: character - study time, with zero indicating baseline
#' @param exposure: character - name of the exposure variable, binary (0/1)
#' @param dtime: character - name of the variable containing subject dropout time
#' @param where: numeric vector of dimension aparms, indicating which subgroup of exposure the dropout time proportion should be calculated

#' @return A list with the following elements: dataframe of model fixed effect estimates, adjusted for dropout effects
#' 1) **estimate**: a vector of the model fixed effect estimates, adjusted for dropout
#' 2) **vcov**: dropout adjusted Variance-covariance matrix for **estimate**, computed using the delta method
#' 3) **lvec**: the vector used to compute **estimate**
#' 4) **bvec**: derivative vector of **lvec**
#' 5) **bvcov**: block diagonal matrix of the model estimated variance matrix fixed effects and the variance matrix of the dropout time proportions
#' @export


drop_adj_all <- function(model_fit
                         , aparms =c('(Intercept)','time','exposure','time:exposure') # parameters with hypothesized dropout effects, needing adjustment by dparms
                         , dparms =c('dtime','time:dtime','exposure:dtime','time:exposure:dtime')# dropout parameters to adjust aparms 
                         , time = 'time' # study time, with zero indicating baseline
                         , exposure = 'exposure' # name of the exposure variable, binary (0/1)
                         , dtime = 'dtime' #name of the variable containing subject dropout time
                         , where = c(0,0,1,1) # numeric vector of dimension aparms, indicating which subgroup of exposure the dropout time proportion should be calculated
) {
  
  tmpdat <- model_fit$data[which(model_fit$data[, time, drop = T] == 0),]
  
  effect_names <- row.names(model_fit$fixed)

  # proportion of drop times: exposed
  tmpdat1 <- tmpdat[which(tmpdat[, exposure, drop = T] == 1),]
  p_time1 <- if(is.null(model_fit$func_call$sweights)) {
    table(tmpdat1[,dtime,drop=T]) / nrow(tmpdat1)
  } else {tapply(tmpdat1[, model_fit$func_call$sweights, drop = T],
                 tmpdat1[,dtime,drop=T], sum) / sum(tmpdat1[, model_fit$func_call$sweights, drop = T])}
  times1 <- as.numeric(names(p_time1))

  # variance - covariance matrix for proportion of dropout times

  pcov1 <- (diag(p_time1) - p_time1 %*% t(p_time1)) / nrow(tmpdat1)

  # proportion of drop times: unexposed
  tmpdat0 <- tmpdat[which(tmpdat[, exposure, drop = T] == 0),]
  p_time0 <- if(is.null(model_fit$func_call$sweights)) {
    table(tmpdat0[,dtime,drop=T]) / nrow(tmpdat0)
  } else {tapply(tmpdat0[, model_fit$func_call$sweights, drop = T],
                 tmpdat0[,dtime,drop=T], sum) / sum(tmpdat0[, model_fit$func_call$sweights, drop = T])}
  times0 <- as.numeric(names(p_time0))

  # variance - covariance matrix for proportion of dropout times

  pcov0 <- (diag(p_time0) - p_time0 %*% t(p_time0)) / nrow(tmpdat0)

  vcov <- if(is.null(model_fit$rvcov)) {model_fit$vcov} else {model_fit$rvcov}
  vcov <- vcov[1:nrow(model_fit$fixed),  1:nrow(model_fit$fixed)]
  bvcov <- Matrix::bdiag(vcov, pcov1, pcov0)

  # compute vectors to compute corrected estimates and their derivatives for the delta method approximation
  lvec <- list()
  bvec <- list()
  for(i in 1:length(aparms)) {

    tmpvec <- rep(0, length(effect_names))
    tmpvec[match(aparms[i], effect_names)] <- 1
    if(where[i] == 0 ) {

      tmpvec[match(dparms[i], effect_names)] <- times0 %*% p_time0
      bvec[[i]] <- c(tmpvec, rep(0, ncol(pcov1)), times0 * sum(model_fit$fixed$estimate[match(aparms[i], effect_names)]))

    } else if(where[i] == 1 ) {

      tmpvec[match(dparms[i], effect_names)] <- times1 %*% p_time1
      bvec[[i]] <- c(tmpvec, times1 * sum(model_fit$fixed$estimate[match(aparms[i], effect_names)]), rep(0, ncol(pcov0)))

    } else {
      print('Values of where parameter should be 0 or 1')
      stop()
    }

    lvec[[i]] <- tmpvec

  }

  lvec <- do.call('cbind', lvec)
  bvec <- do.call('cbind', bvec)

  estimate <- as.numeric(t(as.numeric(model_fit$fixed$estimate)) %*% lvec)
  
  names(estimate) <- aparms
  
  tmp <- list(
    estimate = estimate
    ,vcov    = as.matrix(t(bvec) %*% bvcov %*% bvec)
    ,lvec    = lvec
    ,bvec    = bvec
    ,bvcov   = bvcov
  )


  tmp

}
