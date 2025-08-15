#' Single non-parametric bootstrap

#' @description 
#' Returns a single iteration non-parametric bootstrap data frame. The new data frame will include an integer record key: .rsUnit.

#' @param df: a dataframe to select a bootstrap sample from 

#' @return A bootstrap resampled version of the submitted dataframe with an added psuedo-identifier (.rsUnit.)

#' @export
#' 
#' @examples
#' 
#' bwdat <- MASS::birthwt
#' 
#' bdat <- np_boot(bwdat)

np_boot <- function(df, ...) {
    toReturn <- df[sample(nrow(df), replace=TRUE), , drop=FALSE]
    toReturn$.rsUnit. <- 1:nrow(toReturn)
    return(toReturn)
}

#' Single stratified non-parametric bootstrap

#' @description 
#' Returns a stratified non-parametric bootstrap data frame.

#' @param df: a dataframe to select a bootstrap sample from 
#' @param strata: character, name of a variable that defines stratum  membership 

#' @return A bootstrap resampled (within specified strata) version of the submitted dataframe with an added psuedo-identifier (.rsUnit.)

#' @export
#' 
#' @examples
#' 
#' bwdat <- MASS::birthwt
#' 
#' sbdat <- strata_boot(bwdat, 'smoke')
#' 
#' table(bwdat$smoke)
#' table(sbdat$smoke)

strata_boot <- function(df,strata, ...){
  if (is.null(strata)) return(NULL)
  else {
    samp  <- by(df,df[,strata],np_boot)
    outDF <- do.call("rbind",samp) 
    return(outDF)
  }
}


#' Single non-parametric cluster bootstrap

#' @description 
#' returns a single iteration non-parametric cluster bootstrap data frame

#' @param df: a dataframe to select a bootstrap sample from 
#' @param cluster: character, name of a variable that defines cluster membership 

#' @return Output:
#' A cluster bootstrap version of the input data, with an additional variable meant as
#' a synthetic cluster ID that can be used validly as input to a function that uses 
#' clustering as part of the algorithm: .ClusID.

#' @export
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
#' cdat <- np_clus_boot(mdat, 'id')

np_clus_boot <- function(df,cluster, ...){
    
    #sample unique cluster id's with replacement
    ids <- data.frame(unique(df[cluster]) )
    # ids
    samp <- np_boot(ids)[,names(ids), drop = F]
    samp$.ClusID. <- 1:nrow(samp)
    # samp
    # merge back to source data to retrieve full cluster info
    merge(samp,df,by=cluster)
}

#' Single stratified non-parametric cluster bootstrap

#' @description 
#' returns a single iteration non-parametric stratified cluster bootstrap data frame. Cluster IDs in different strata are assumed to be separate entities.

#' @param df: a dataframe to select a bootstrap sample from 
#' @param strata: character, name of a variable that defines stratum  membership 
#' @param cluster: character, name of a variable that defines cluster membership 

#' @return Output:
#' A cluster bootstrap (within strata) version of the input data, with an additional variable meant as
#' a synthetic cluster ID that can be used validly as input to a function that uses 
#' clustering as part of the algorithm: .ClusID.

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
#' scdat <- strata_clus_boot(mdat, 'exposure', 'id')


strata_clus_boot <- function(df,strata,cluster, ...) {
    #sample unique strata/cluster id's with replacement
    ids <- data.frame(unique(df[c(cluster,strata)]) )
    # ids
    # select a stratified random sample of clusters
    samp  <- by(ids,ids[,strata],np_boot)
    # stack each dataframe in to one large
    stackDF <- do.call("rbind",samp)[, names(ids)]
    # assign each selected cluster an ID (to differentiate those selected >= 1 time)
    stackDF$.ClusID. <- 1:nrow(stackDF)
    # stackDF
    # re-merge with original data, then output
    merge(stackDF,df,by=c(cluster,strata))
}


