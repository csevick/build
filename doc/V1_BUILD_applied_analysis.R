## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)


## ----eval = F-----------------------------------------------------------------
# install.packages("devtools")

## ----eval = F-----------------------------------------------------------------
# 
# devtools::install_github("csevick/build")
# 

## -----------------------------------------------------------------------------

library(build)

set.seed(209587)

dat <- simuDat(
    B            = c(-3,2,5,2.5,0),
    G            = matrix(c(4), ncol = 1),
    random       = NULL,
    p_exposed    = 0.5,
    p_confounder = 0.25,
    conf         = 0,
    N            = 5000,
    M            = 5,
    D_B          = c(0, -2, -1, -0.5, 0),
    D_F          = c(1, 1, 1, 1, 1),
    D_parms      = list(c(1,1), c(1,1.5))
  ) 

# create a subset with non-random dropout
mdat <- subset(dat, time <= dtime)

## ----eval = F, include = F----------------------------------------------------
# \noindent
# Of course, the main data source is not supposed to have the exposure data. This is why we are using an outcome dependent sampling design to isolate a small subset of participants to submit to further analysis to determine their exposure status. The following, final, transformation simulates this applied situation.

## -----------------------------------------------------------------------------
# create look up table for exposure status
## this will represent the assessment results for the selected random sample 
exposure_status <- unique(mdat[c('id', 'exposure')])

# delete exposure status from the data source (we are not supposed to know!)
mdat <- subset(mdat, select = -exposure)


## -----------------------------------------------------------------------------
# create the sample frame
s_frame <- build_sample_frame(data        = mdat
                            , study_id    = 'id'
                            , outcome     = 'y'
                            , drop_strata = 'dtime')

# the strata cell counts are:
s_frame$strata_counts

# the marginal outcome strata counts are:
rowSums(s_frame$strata_counts)

# the marginal dropout time strata counts are:
colSums(s_frame$strata_counts)

# the first 6 rows of the sample frame
head(s_frame$sample_frame)

## -----------------------------------------------------------------------------

# create the sampling plan: marginal outcome samples at [75, 450, 75]
s_plan <- build_sample_plan(n_samp     = c(75, 450, 75)
                          , N_bystrata = s_frame$strata_counts
                          , drop_grad  = 1
                            )

# the [75, 450, 75] participants are distributed as:
s_plan$N_sample

## -----------------------------------------------------------------------------
# draw the stratified sample
set.seed(4580978)
sample_subj <- strat_samp(df = s_frame$sample_frame
                          , oc_strata = 'oc_strata'
                          , drop_strata = 'dtime'
                          , sample_plan = s_plan$sample_plan
                          , sample_by = 'n')

# The first 6 rows from the data frame holding the selected participants
head(sample_subj)

## -----------------------------------------------------------------------------
# get exposure status
sample_subj <- merge(sample_subj, exposure_status, by = 'id')


## ----eval = F, include = F----------------------------------------------------
# 
# \noindent
# In the final step, prior to analysis the selected participants are reunited with their longitudinal data, along with the newly derived exposure variable and sampling weights, and statistical analyses can begin.

## -----------------------------------------------------------------------------
sample <- merge(mdat
              , sample_subj[, c('id'
                              , 'oc_strata'
                              , 'drop_strata'
                              , 'sampWT'
                              , 'exposure')]
              , by = 'id')

# sort for good form
sample <- sample[order(sample$id, sample$time),]



## -----------------------------------------------------------------------------

# fit the model
WL_fit <- logit_mix(df = sample 
                  , fixed = y ~ time + exposure + time:exposure +
                                dtime + time:dtime + exposure:dtime + 
                                time:exposure:dtime
                  , random = ~ 1
                  , cluster = 'id'
                  , sweights = 'sampWT'
                  , robust = T
                  , qpoints = 20
                  )

# fixed effect estimates
knitr::kable(WL_fit$fixed, digits = 2)

## -----------------------------------------------------------------------------

edist0 <- with(subset(sample_subj, exposure == 0), table(dtime)/length(dtime))
edist0

## -----------------------------------------------------------------------------

edist1 <- with(subset(sample_subj, exposure == 1), table(dtime)/length(dtime))
edist1

## -----------------------------------------------------------------------------

slope_n_diff <- function(fit_obj) {
  
  # compute empirical dropout time distributions
  edist0 <- with(subset(fit_obj$data, exposure == 0 & time == 0)
               , table(dtime)/length(dtime))
  edist1 <- with(subset(fit_obj$data, exposure == 1 & time == 0)
               , table(dtime)/length(dtime))

  # extract fixed effects
  feff <- fit_obj$fixed$estimate
  names(feff) <- row.names(fit_obj$fixed)
  
  # unique ordered dropout times
  dtimes <- unique(fit_obj$data$dtime)
  dtimes <- sort(dtimes)
  
  # compute marginal effects
  results <- data.frame(
    slope_unexposed = feff['time'] + edist0 %*% dtimes * feff['time:dtime']
  , slope_exposed = feff['time'] + feff['time:exposure'] + 
                    edist1 %*% dtimes * feff['time:dtime'] + 
                    edist1 %*% dtimes * feff['time:exposure:dtime']
  )
  
  results$slope_diff <- with(results, slope_exposed - slope_unexposed)
  
  return(results)

}

# compute the observed marginal slopes and slope difference
obs_splope_n_diff <- slope_n_diff(WL_fit)
obs_splope_n_diff

## -----------------------------------------------------------------------------

# function to carryout a single analysis
boot_fun <- function(data) {
  
  #draw a bootstrap sample
  B <- strata_clus_boot(df = data
                      , strata = c('oc_strata','dtime')
                      , cluster = 'id')
  # fit the model
  ## note: robust standard errors are intensive and not needed in a bootstrap analysis
  B_WL_fit <- logit_mix(df = B
                      , fixed = y ~ time + exposure + time:exposure +
                                    dtime + time:dtime + exposure:dtime + 
                                    time:exposure:dtime
                      , random = ~ 1
                      , cluster = '.ClusID.'
                      , sweights = 'sampWT'
                      , robust = F
                      , qpoints = 20
                      )

  return(slope_n_diff(B_WL_fit))
}

boot_fun(sample)

## ----warning=F, message=F-----------------------------------------------------
library(doParallel)
library(doRNG)

# set up a cluster
n_cores <- detectCores() - 1
registerDoParallel(cores = n_cores)
  
# set the number of bootstrap iterations
n_boots <- 100
system.time({
B_result <- foreach(1:n_boots
                 , .packages = c('build')
                 , .options.RNG = 4879549
                 , .combine = 'rbind'
) %dorng% {
  boot_fun(sample)
}
})

#stop the cluster
stopImplicitCluster()

# compute bootstrap SE's
bse <- sqrt(diag(var(B_result)))

# compute test statistics
Z <- obs_splope_n_diff / bse

# compute p-values
p_val <- sapply(Z, function(Z) {pchisq(q = Z^2, df = 1, lower.tail = F)})

# compute 95% confidence limits
lcl <- obs_splope_n_diff - qnorm(0.975)*bse
ucl <- obs_splope_n_diff + qnorm(0.975)*bse

# summary table
boot_tab <- data.frame(
  estimate = t(obs_splope_n_diff)
 ,lcl_95 = t(lcl)
 ,ucl_95 = t(ucl)
 ,pvalue = format.pval(p_val, eps=0.001, digits = 2)
)
row.names(boot_tab) <- names(obs_splope_n_diff)
knitr::kable(boot_tab
             , digits = 2
             , caption = 'Results of bootstap analysis of marginal parameters')

