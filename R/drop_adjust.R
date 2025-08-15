#' Function to apply the dropout correction to individual parameter estimates

#' @description 
#' This function is only validated for simple conditional linear models

#' @param model_fit: a fit object from logit_mix
#' @param data: the data used to estimate the model in model_fit
#' @param var: character - variable to be adjusted
#' @param dvar: character - the dropout specific variable used to adjust
#' @param where: Boolean - subset specification for the correct group to calculate 
#' dropout time proportions
#' @param effect_name: character - if var is a combination then specify this for the new name
#' @param weights: character - name of the sampling weight column in 'data'

#' @return A one row dataframe containing the dropout adjusted fixed effect estimate. Columns include:
#'  1) **estimate**: estimate of the specified parameter, adjusted for dropout
#'  2) **std_err**: dropout adjusted standard error of **estimate**, computed using the delta method
#'  
#' @export

drop_adjust <- function(model_fit, data, var, dvar, where, effect_name, weights = NULL) {

  if(missing(data)) {data <- model_fit$data}
  
  tmpdat <- subset(data, time == 0)
  if(!missing(where)) {
    where <- deparse(substitute(where))
    tmpdat <- subset(tmpdat, subset = eval(parse(text = where)))
  }
  
  effect_names <- row.names(model_fit$fixed)
  
  # proportion of drop times
  p_time <- if(is.null(weights)) {
    table(tmpdat$dtime) / nrow(tmpdat)
  } else {tapply(tmpdat[, weights, drop = T], tmpdat$dtime, sum) / sum(tmpdat[, weights, drop = T])}
  times <- as.numeric(names(p_time))
  
  # variance - covariance matrix for proportion of dropout times 
  pcov <- (diag(p_time) - p_time %*% t(p_time)) / nrow(tmpdat)
  
  vcov <- if(is.null(model_fit$rvcov)) {model_fit$vcov} else {model_fit$rvcov}
  vcov <- vcov[1:nrow(model_fit$fixed), 1: 1:nrow(model_fit$fixed)]
  big_vcov <- Matrix::bdiag(vcov, pcov)
  
  vec <- rep(0, length(effect_names))
  vec[match(var, effect_names)] <- 1
  vec[match(dvar, effect_names)] <- times %*% p_time
  big_vec <- c(vec, times * sum(model_fit$fixed$estimate[match(var, effect_names)]))
  
  tmp <- data.frame(
    estimate = as.numeric(model_fit$fixed$estimate %*% vec),
    std_err   = as.numeric(sqrt(t(big_vec) %*% big_vcov %*% big_vec))
  )
  if(missing(effect_name)) row.names(tmp) <- var else row.names(tmp) <- effect_name
  tmp
  
}

