#' Sample selection

#' @description 
#' Sampling from a data frame without replacement

#' @param df: the name of a dataframe to select rows from
#' @param p: numeric (0-1) selection probability
#' @param n: integer number to be randomly selected

#' @return A dataframe consisting of a random sample of rows from the original.
#' 
#' @export
#' 
#' @examples
#' 
#' bwdat <- MASS::birthwt
#' 
#' sample_dat <- samp_fun(bwdat, n = 20)
#' 

samp_fun <- function(df, p=NULL, n=NULL) {
  
  if(!is.null(p) & !is.null(n)) {
    print('ERROR: p and n cannot both be specified')
    stop()
  }
  
  if(!is.null(n)) {
    return(df[sample(1:nrow(df), n),])
  }
  else if(!is.null(p)) {
    return(df[rbinom(nrow(df), 1, p)==1,])
  }
  else {
    print('ERROR: p or n must be given a value')
    stop()
  }
}

#' Stratified sample selection

#' @description 
#' Stratified sampling from a data frame without replacement.
#' Within each stratum a fixed proportion of rows will be drawn

#' @param df: the name of a dataframe to select rows from
#' @param oc_strata: character, variable in the dataframe that specifies the outcome strata membership
#' @param drop_strata: character, variable in the dataframe that specifies the dropout strata membership
#' @param sample_plan: dataframe - contains the sampling plan from 'build_sample_plan' 
#' @param sample_by: (p/n) specify whether to sample by counts or probability (default is 'n')
#' 
#' @return A dataframe consisting of a stratified random sample of rows from the original.
#' 
#' @export
#'
#' @examples
#' 
#' bwdat <- MASS::birthwt
#' 
#' table(bwdat$race)
#' 
#' s_plan <- data.frame(race=c(1,2,3), n = c(20, 25, 15), n_population = as.numeric(table(bwdat$race)))
#' s_plan <- within(s_plan, {p <- n/n_population})
#' 
#' sample_dat <- strat_samp(bwdat, oc_strata = 'race', sample_plan = s_plan,sample_by='n')
#' 

strat_samp <- function(df,oc_strata = NULL, drop_strata = NULL, sample_plan,sample_by = 'n'){
    
    strata = c(oc_strata, drop_strata)
    
    sample_plan <- subset(sample_plan, n>0 & !is.na(n))
    
    sample_by <- tolower(sample_by)
    if(sample_by == 'n'){
      p <- NULL
      n <- sample_plan[, 'n', drop = T]
    } else if (sample_by == 'p') {
      p <- sample_plan[, 'p', drop = T]
      n <- NULL
    } else {
      print('ERROR: sampling must be by n or p')
      stop()
    }
    lst = list()
    for(i in 1:nrow(sample_plan)){
        dfsub <- merge(df, sample_plan[i,c(strata, 'drop_strata')]) #df[strVec == grps[i], ,drop = F]
        dfSample <- samp_fun(dfsub, p = p[i], n = n[i])
        dfSample$sampWT <- 1/sample_plan[i , 'p', drop = T] #if(is.null(p)) {gsize[i]/n[i]} else {1/p[i]}
        lst[[i]] <- dfSample
        rm(dfsub)
    }
    return(do.call('rbind',lst))
}


#' BUILD sample frame 

#' @description 
#' Accepts a longitudinal dataset, computes the outcome strata and outputs a dataframe with the sampling frame for the analysis.
#' 
#' @param data: the name of a dataframe to select rows from
#' @param study_id: character, variable in the dataframe that identifies the subject
#' @param outcome: character, variable in the dataframe that specifies the study outcome
#' @param drop_strata: character - variable in the dataframe that identifies the dropout strata
#' 
#'  @return A list of two elements:
#'   1) sample_frame: a dataframe with one row per sampling unit. columns include the study ID from the original data, the dropout strata membership (drop_strata) and the derived outcome strata (oc_strata).
#'   2) strata_counts: a matrix of unit counts within each outcome (rows) and dropout (columns) strata.
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

build_sample_frame <- function(data, study_id, outcome, drop_strata){

  sample_frame <- by(data
                    , data[, study_id, drop=T]
                    , function(x){
                      n <- nrow(x)
                      y <- sum(x[, outcome, drop = T])
                      x$oc_strata <- if(is.na(y)) {NA} else if(y == 0) {1} else if(0 < y & y < n) {2} else {3}
                      x
                    }) |>
   do.call('rbind', args = _)
  sample_frame <- sample_frame[, c(study_id, drop_strata, 'oc_strata' )] |> unique()

    list(sample_frame = sample_frame
        , strata_counts = table(sample_frame[, c('oc_strata', drop_strata)])
    )
  
}

or_change <- function(p, or) {
  # update a probability by an odds ratio (Fliess 1981)
  
  return(or*p/(or*p+(1-p)))
  
}

gslope <- function(v){
  a <- 1:length(v)
  b <- a - mean(a)
  #s <- sign(b)
  #b2 <- ceiling(abs(b))*s
  
  return(b)
}

#' Determine strata oversampling odds ratios to use in build_sample_plan()

#' @description 
#' Accepts desired sample size across marginal outcome strata and a matrix of outcome and dropout stratum counts
#' as well as desired adjustment to the sampling gradient across the dropout strata. Returns the 
#' appropriate odds ratio oversampling for the outcome strata to be used in build_sample_plan(). 
#' The return value is a result object from the optim function.
#' 
#' @param n_samp: a three element vector (if specifying marginal outcome strata sampling goal)
#' @param N_bystrata: matrix of outcome by dropout counts. dimensions must be named
#' @param drop_grad: numeric scalar. a value of zero will indicate proportional sampling across the dropout strata. a non-zero value is a slope (expressed as an odds ratio to adjust the sampling gradient across dropout strata).
#' @param drop_unif: (T/F) if true then include the drop_grad as a parameter to optimize seeking a uniform sample (marginally) over the dropout strata
#' @param drop_weight: a vector of weights equivalent to the expected number of longitudinal measures per participant in the equivalent ordered dropout strata. this is only used when drop_unit = T and it is desired to balance the number of person-measures across dropout time
#' 
#' @return A object from the R optim function with par representing the estimated vector of log-odds parameters to transform the sampling probabilities for the ODS sample. Intended as input to the argument of the build_sample_plan() function.
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
#' s_or_plan <- build_strata_ors(n_samp = c(50, 200, 50), N_bystrata = s_frame$strata_counts, drop_grad = 1)
#'
#' s_plan <- build_sample_plan(n_samp = 300, N_bystrata = s_frame$strata_counts, oc_sample_adj = s_or_plan$par, drop_grad = 1)
#' 
#' s_sample <- strat_samp(s_frame$sample_frame,oc_strata = 'oc_strata', drop_strata = 'dtime', sample_plan=s_plan$sample_plan,sample_by = 'n')

build_strata_ors <- function(n_samp, N_bystrata, drop_unif = F, drop_weight = NULL, drop_grad = 1, control = list(maxit = 1000)) {
  
  if(is.null(drop_weight)) {drop_weight <- rep(1, ncol(N_bystrata))}
  
  par <-  if(drop_unif == T) {c(1,1,1,1)} else {c(1,1,1)}
  optim(par #c(1,1,1,1)
        ,function(par){
          
          if(drop_unif == T) {drop_grad <- par[4]}
          
          p_obs <- N_bystrata / sum(N_bystrata)
          
          n_samp_base <- p_obs * sum(n_samp)
          
          p_samp_base <- n_samp_base / N_bystrata
          for(i in 1:nrow(p_samp_base)){
            for(j in 1:ncol(p_samp_base)){
              p_samp_base[i,j] <- or_change(p_samp_base[i,j], par[i])
            }
          }
          
          # adjust sampling probabilities for change in the dropout gradient
          if(drop_grad != 1) {
            d_p_adj <- exp(gslope(1:ncol(p_samp_base))*log(drop_grad))
            
            for(i in 1:nrow(p_samp_base)){
              for(j in 1:ncol(p_samp_base)){
                p_samp_base[i,j] <- or_change(p_samp_base[i,j],  d_p_adj[j])
              }
            }
          }
          
          
          
          s1 <- sum((n_samp - rowSums(p_samp_base*N_bystrata))^2) 
          
          if(drop_unif == T) {
            unif_base <- rep(1/ncol(N_bystrata), ncol(N_bystrata))
            s1 <- s1 + sum((sum(n_samp)*unif_base*mean(drop_weight) - colSums(p_samp_base*N_bystrata)*drop_weight)^2)
          }
          
          s1
        }
        , control = control)
}

#' BUILD sample plan 

#' @description 
#' Accepts desired sample size and a matrix of outcome and dropout stratum counts
#' along with specification for stratum oversampling and outputs a dataframe with the 
#' needed stratum sampling numbers for the sampling function
#' 
#' @param n_samp: a scalar (if only specifying total desired sample size), or three element vector (if specifying marginal outcome strata sampling goal)
#' @param N_bystrata: a matrix of outcome by dropout counts. dimensions must be named
#' @param oc_sample_adj: a numeric vector of dim 3. Specified only if a scalar sample size is specified. Amount to adjust the observed outcome strata proportions expressed as an odds increase in the sampling proportion of the stratum. All three may be specified and the result is standardised to the sample size requested in 'n_samp'.
#' @param drop_grad: numeric scalar. a value of zero will indicate proportional sampling across the dropout strata. a non-zero value is a slope (expressed as an odds ratio to adjust the sampling gradient across dropout strata).
#' 
#' @return A list of three elements:
#' 1) sample_plan: a dataframe of all possible outcome and dropout strata (oc_strata, drop_strata), the number to sample from each combination (n), the total number in the population (n_population) and the resulting sampling probability (p)
#' 2) P_matrix: a matrix of sampling probabilities for that combination of outcome strata (row) and dropout strata (column).
#' N) N_sample: a matrix of sample sizes for that combination of outcome strata (row) and dropout strata (column).
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

build_sample_plan <- function(n_samp, N_bystrata, oc_sample_adj=c(NA, 1, 1), drop_grad = 1) {
  
  stopifnot(length(oc_sample_adj) == 3)
  
  stopifnot(length(n_samp) %in% c(1,3))
  
  p_obs <- N_bystrata / sum(N_bystrata)
  
  n_samp_base <- p_obs * sum(n_samp)
  
  p_samp_base <- n_samp_base / N_bystrata
  
  # if marginal sample numbers for the outcome strata are given then determine adjustment
  if(length(n_samp) == 3) {
    tmp_s <- n_samp / sum(n_samp)
    tmp_o <- rowSums(N_bystrata) / sum(N_bystrata)
    
    oc_sample_adj <- tmp_s*(1-tmp_o)/(tmp_o*(1-tmp_s))
  }
  
  # adjust sampling probabilities for desired marginal outcome numbers
  for(i in 1:nrow(p_samp_base)){
    for(j in 1:ncol(p_samp_base)){
      p_samp_base[i,j] <- or_change(p_samp_base[i,j], oc_sample_adj[i])
    }
  }
  
  # adjust sampling probabilities for change in the dropout gradient
  if(drop_grad != 1) {
    d_p_adj <- exp(gslope(1:ncol(p_samp_base))*log(drop_grad))
    
    for(i in 1:nrow(p_samp_base)){
      for(j in 1:ncol(p_samp_base)){
        p_samp_base[i,j] <- or_change(p_samp_base[i,j],  d_p_adj[j])
      }
    }
  }
  
  N_sample <- N_bystrata*p_samp_base
  N_sample[is.na(N_sample)] <- 0
  
  # standardize N sampled to reflect n_samp
  if(length(n_samp) == 3) {
    N_sample <- N_sample/rowSums(N_sample)*n_samp
  } else {
    N_sample <- N_sample/sum(N_sample)*n_samp
  }
  
  N_sample[N_sample > 0 & N_sample < 1] <- 1
  N_sample <- round(N_sample)
  P_sample <- N_sample / N_bystrata
  P_sample[is.na(P_sample)] <- 0
  if(!all(P_sample <= 1 & P_sample >= 0)) {
    print("ERROR: at least one sampling probability is not within [0, 1]")
    return(list(bad_P_sample = P_sample))
  }
  
  long_N <- as.data.frame(N_sample)
  long_NS <- as.data.frame(N_bystrata)
  long_dat <- merge(long_N, long_NS, by = names(long_N)[1:2])
  names(long_dat) <- c(names(long_N)[1:2], 'n', 'n_population')
  long_dat$p <- long_dat$n / long_dat$n_population
  long_dat$drop_strata <- do.call('c',by(long_dat
                                         ,long_dat[,'oc_strata', drop = T]
                                         ,\(x){1:nrow(x)}
                                         , simplify = F))
  return(list(sample_plan = long_dat, P_matrix = P_sample, N_sample = N_sample))
  
}
