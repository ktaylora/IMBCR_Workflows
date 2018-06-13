# Title     : TODO
# Objective : TODO
# Created by: ktaylora
# Created on: 5/23/18

# load a previously built gdistsamp workflow for us to use for re-fitting our top models
# and then do a power analysis using bootstrapping to assess the potential
# influence of downsampling

load(commandArgs(trailingOnly=T))

# define our global constants
DOWNSAMPLING_THRESHOLD    <- 0.3
N_BS_REPLICATES           <- 999

#
# Local Functions
#

fit_gdistsamp <- function(lambdas=NULL, umdf=NULL, mixture="P"){
  if(length(lambdas)>1){
    cl <- parallel::makeCluster(parallel::detectCores()-1)
    parallel::clusterExport(cl, varlist=c("umdf","mixture"), envir=environment())
    unmarked_models <- parallel::parLapply(
        cl=cl,
        X=lambdas,
        fun=function(m){
          ret <- try(unmarked::gdistsamp(
            pformula=as.formula("~1"),
            lambdaformula=as.formula(paste("~", m, sep="")),
            phiformula=as.formula("~1"),
            data=umdf,
            se=T,
            K=max(rowSums(umdf@y)),
            keyfun="halfnorm",
            unitsOut="kmsq",
            mixture=mixture,
            output="abund"
          ))
         if(class(ret) == "try-error"){
           return(NA)
         } else {
           return(ret)
         }
      })
    parallel::stopCluster(cl);
    return(unmarked_models);
  } else {
     ret <- try(unmarked::gdistsamp(
      pformula=as.formula("~1"),
      lambdaformula=as.formula(paste("~", unlist(lambdas), sep="")),
      phiformula=as.formula("~1"),
      data=umdf,
      se=T,
      K=max(rowSums(umdf@y)),
      keyfun="halfnorm",
      unitsOut="kmsq",
      mixture=mixture,
      output="abund"
    ))
    if(class(ret) == "try-error"){
      return(NA)
    } else {
      return(ret)
    }
  }
}
#' a re-worked local implementation of the sample() function that can
#' accomodate AB's vision for minimum sampling density per IMBCR transect
downsample_stratified <- function(x=NULL, size=NULL, strata=NULL, return_rows=T){
  MAX_SEARCH_ITER = 100
  MIN_STRATUM_TRANSECTS = 2
  if(is.null(x)) stop("x= argument not specified")
  if(is.null(size)) stop("size= argument not specified")
  # if the user didn't specify strata, try and figure them out from
  # the transect strings specified by x=
  if(is.null(strata)){
    state_bcr <- sapply(
      strsplit(as.character(x), split="-"),
      FUN=function(x) paste(x[1],x[2],sep="-")
    )
    statum_code_minus_letter <- gsub(
      sapply(
        strsplit(as.character(x), split="-"),
        FUN=function(x) x[3]
      ),
      pattern="[0-9]",
      replacement=""
    )
    # merge the state-bcr and stratum strings
    strata <- paste(state_bcr, statum_code_minus_letter,sep="-")
  }
  # make a pivot table of our strata and figure out how many we have to part with from
  sample_counts <- sample(table(strata))
  # sanity-check -- do we have enough transects in the original dataset
  # to cover the size requested by the user?
  if(sum(sample_counts) < size) stop("size= argument is too large to attempt a stratification")
  full_imbcr_dataset <- data.frame(transect=x, strata=strata)
  rows_to_keep <- vector() # will contain rows from the full dataset we are going to retain
  i <- 1
  while(length(rows_to_keep) < size && i < MAX_SEARCH_ITER) {
    for(j in 1:length(unique(strata))){
      # do we have at least two transects in this stratum?
      if(sample_counts[j]>MIN_STRATUM_TRANSECTS){
         strata_count <- sample(MIN_STRATUM_TRANSECTS:sample_counts[j], size=1) # how many transects are we going to keep from this stratum?
         rows <- which(
           grepl(
             full_imbcr_dataset$transect,
             pattern=names(sample_counts[j])
           )
         )
         rows_to_keep <- unique(append(rows_to_keep, sample(rows, size=strata_count)))
      }
      # if not, then keep these transects in our final sample
      else {
         rows <- which(
           grepl(
             full_imbcr_dataset$transect,
             pattern=names(sample_counts[j])
           )
         )
         rows_to_keep <- unique(append(rows_to_keep, rows))
      }
      # sanity-check : do we satisfy our n rows requirement?
      if(length(rows_to_keep)>=size) break
    }
    # sanity-check : break out of an endless-loop
    i <- i+1;
    if(i == MAX_SEARCH_ITER) stop("failed to find a sample that satisfies size= argument -- this shouldn't happen")
  }
  if(return_rows){
    return(rows_to_keep)
  } else {
    full_imbcr_dataset$keep <- FALSE;
    full_imbcr_dataset$keep[rows_to_keep] <- TRUE;
    return(full_imbcr_dataset)
  }
}
#' randomly downsample a dataset N=replicates times and calculate the density
#' and standard error for each run. Will return the full table of all replicate
#' runs (with statistics)
bs_calc_power <- function(
  replicates=N_BS_REPLICATES,
  s="/global_workspace/imbcr_number_crunching/results/RawData_PLJV_IMBCR_20171017.csv",
  unmarked_models=NULL,
  original_formulas=NULL,
  top_model=NULL,
  downsample_fun=sample
){
  # read-in our IMBCR transect data
  s <- OpenIMBCR:::scrub_imbcr_df(
    OpenIMBCR:::imbcrTableToShapefile(s),
    four_letter_code = toupper(argv[2])
  )
  # calculate our detections and sampling effort
  detections <- OpenIMBCR:::calc_dist_bins(s)
  effort     <- as.vector(OpenIMBCR:::calc_transect_effort(s))
  # which formula are we going to test?
  formula <- original_formulas[top_model]
  # re-build our training data frame
  s <- OpenIMBCR:::calc_route_centroids(
      s=s,
      four_letter_code=toupper(argv[2])
  )
  # join with our habitat covariates
  s <- OpenIMBCR:::spatial_join(s, units)
  s$effort <- effort

  # set-up our cluster
  cl <- parallel::makeCluster(parallel::detectCores()-1)
  parallel::clusterExport(
    cl,
    varlist=c(
      "s","detections","vars","unmarked_models",
      "top_model","fit_gdistsamp", "downsample_fun",
      "formula"
    ),
    envir=environment()
  )
  parallel::clusterExport(
    cl,
    varlist=c("DOWNSAMPLING_THRESHOLD"),
    envir=globalenv()
  )
  # parallelize an IMBCR model re-fit operation across N=replicates bootstrap
  # replications -- return a table for each replicate that we will then
  # merge and return to the user
  replicates <- do.call(
    rbind,
    parallel::parLapply(
      cl=cl,
      X=1:replicates,
      fun=function(x){
        require(unmarked)
        # pre-allocate a table so that we can return at-least an NA value
        # if a run fails
        return_table <- data.frame(
          density=NA,
          density_downsampled=NA,
          se=NA,
          se_downsampled=NA,
          aic=NA,
          aic_downsampled=NA
        )
        # determine rows to keep that satisfy our DOWNSAMPLING_THRESHOLD
        if(identical(downsample_fun, sample)){
          sample <- downsample_fun(
              1:nrow(s@data),
              size=nrow(s@data)*(1-DOWNSAMPLING_THRESHOLD)
            )
        # assume we are using the min stratification strategy
        } else {
          sample <- downsample_fun(
            x=s$transect,
            size=nrow(s@data)*(1-DOWNSAMPLING_THRESHOLD)
          )
        }
        # randomly downsample our unmarked data.frame to the specified density
        umdf <- unmarked::unmarkedFrameGDS(
          y=as.matrix(detections$y),
          siteCovs=s@data[,c(vars,'effort')],
          dist.breaks=detections$breaks,
          numPrimary=1,
          survey="point",
          unitsIn="m"
        )
        umdf_downsampled <- unmarked::unmarkedFrameGDS(
          y=as.matrix(detections$y)[sample,],
          siteCovs=s@data[sample,c(vars,'effort')],
          dist.breaks=detections$breaks,
          numPrimary=1,
          survey="point",
          unitsIn="m"
        )
        # re-fit our top model using our regular and our downsampled dataset
        # use some fairly-robust exception handling here to capture when
        # we completely fail to re-fit a model -- this information is useful
        # to us for our power analysis
        try(m <- fit_gdistsamp(
          lambda=formula,
          umdf=umdf,
          mixture=unmarked_models[[top_model]]@mixture
        ))
        if(class(m) %in% c("try-error","logical")){
          return(return_table)
        }
        try(m_downsampled <- fit_gdistsamp(
            lambda=formula,
            umdf=umdf_downsampled,
            mixture=unmarked_models[[top_model]]@mixture
        ))
        if(class(m_downsampled) %in% c("try-error","logical")) {
          return(return_table)
        }
        # calculate density and standard error for our regular dataset
        density <- se <-
          unmarked::predict(
            m,
            type="lambda"
        )[,1:2]
        return_table$density <- mean(density[,1])
        return_table$se <-  mean(se[,2])
        return_table$aic <- m@AIC
        # now for our downsampled dataset
        density_downsampled <- se_downsampled <-
          unmarked::predict(
            m_downsampled,
            type="lambda"
        )[,1:2]
        return_table$density_downsampled <- mean(density_downsampled[,1])
        return_table$se_downsampled <-  mean(se_downsampled[,2])
        return_table$aic_downsampled <- m_downsampled@AIC
        # return to user
        return(return_table)
      }
    )
  )
  parallel::stopCluster(cl)
  return(replicates)
}
#' bootstrap a simple model of a known covariate to use as a basis 
#' for representing the power of a model given a fixed number of transects.
#' By default, this calculates the density + SE for 'total area of shrubland'
#' covariate
calc_power_of_cov <- function(
  replicates=N_BS_REPLICATES,
  var="shrub_ar",
  s="/global_workspace/imbcr_number_crunching/results/RawData_PLJV_IMBCR_20171017.csv",
  formula="poly(shrub_ar, 1, raw = T) + poly(grass_ar, 1, raw = T) + offset(log(effort))",
  four_letter_code=NULL,
  units=NULL,
  downsample_fun=sample,
  sample_size=NULL
){
  # read-in our IMBCR transect data
  s <- OpenIMBCR:::scrub_imbcr_df(
    OpenIMBCR:::imbcrTableToShapefile(s),
    four_letter_code = toupper(four_letter_code)
  )
  # calculate our detections and sampling effort
  detections <- OpenIMBCR:::calc_dist_bins(s)
  effort     <- as.vector(OpenIMBCR:::calc_transect_effort(s))
  # re-build our training data frame
  s <- OpenIMBCR:::calc_route_centroids(
      s=s,
      four_letter_code=toupper(four_letter_code)
    )
  # join with our habitat covariates
  s <- OpenIMBCR:::spatial_join(s, units)
  s$effort <- effort
  # set-up our cluster
  cl <- parallel::makeCluster(parallel::detectCores()-1)
  parallel::clusterExport(
    cl,
    varlist=c(
      "s","detections","fit_gdistsamp","formula","var"
    ),
    envir=environment()
  )
  replicates <- do.call(
    rbind,
    parallel::parLapply(
      cl=cl,
      X=1:replicates,
      fun=function(x){
        require(unmarked)  
        # which sample densities should we keep?
        low    <- which(s@data[,var] <= as.vector(quantile(s@data[,var], probs=0.33)))
        low <- try(sample(
            low, 
            size=sample_size/3
        ))
        # bin our low densities
        if(class(low) == "try-error"){
            low    <- which(s@data[,var] <= as.vector(quantile(s@data[,var], probs=0.33)))
            low <- sample(
            low,
            size=sample_size/3,
            replace=T
            )
        }
        # bin our medium densities
        medium <- which( s@data[,var] > as.vector(quantile(s@data[,var], probs=0.33)) & s@data[,var] <= as.vector(quantile(s@data[,var], probs=0.66)) )
        medium <- try(sample(
            medium, 
            size=sample_size/3
        ))
        if(class(medium) == "try-error"){
            medium <- which( s@data[,var] > as.vector(quantile(s@data[,var], probs=0.33)) & s@data[,var] <= as.vector(quantile(s@data[,var], probs=0.66)) )
            medium <- sample(
            medium,
            size=sample_size/3,
            replace=T
            )
        }
        # bin our high densities
        high <- which( s@data[,var] > as.vector(quantile(s@data[,var], probs=0.66)) )
        high <- try(sample(
            high, 
            size=sample_size/3
        ))
        if(class(high) == "try-error"){
            high <- which( s@data[,var] > as.vector(quantile(s@data[,var], probs=0.66)) )
            high <- sample(
            high,
            size=sample_size/3,
            replace=T
            )
        }
        # build our unmarked data.frame
        keep <- c(low, medium, high)
        umdf <- unmarked::unmarkedFrameGDS(
                y=as.matrix(detections$y)[keep,],
                siteCovs=s@data[keep,c(var,'grass_ar','effort')],
                dist.breaks=detections$breaks,
                numPrimary=1,
                survey="point",
                unitsIn="m"
                )
        # pre-allocate a table so that we can return at-least an NA value
        # if a run fails
        return_table <- data.frame(
            density=NA,
            se=NA,
            aic=NA
        )
        # fit a model
        try(m <- fit_gdistsamp(
            lambda=formula,
            umdf=umdf,
            mixture="NB"
        ))
        if(class(m) %in% c("try-error","logical")){
            return(return_table)
        }
        # record fit statistics in our return table
        density <- se <-
            unmarked::predict(
            m,
            type="lambda"
        )[,1:2]
        return_table$density <- mean(density[,1])
        return_table$se <-  mean(se[,2])
        return_table$aic <- m@AIC
        # return to user
        return(return_table)
    })
  )
  parallel::stopCluster(cl)
  return(replicates)
}
# for 2016, determine the loss of power from current stratification if we
# downsample by certain thresholds

DOWNSAMPLING_THRESHOLD <- 0.3
p_30_perc_reduction_2016 <- bs_calc_power(
  s="/global_workspace/imbcr_number_crunching/results/RawData_PLJV_IMBCR_20161201.csv",
  unmarked_model=unmarked_models,
  original_formulas=original_formulas,
  downsample_fun=downsample_stratified,
  top_model=as.numeric(row.names(model_selection_table@Full[1,]))
)

DOWNSAMPLING_THRESHOLD <- 0.2
p_20_perc_reduction_2016 <- bs_calc_power(
  s="/global_workspace/imbcr_number_crunching/results/RawData_PLJV_IMBCR_20161201.csv",
  unmarked_model=unmarked_models,
  original_formulas=original_formulas,
  downsample_fun=downsample_stratified,
  top_model=as.numeric(row.names(model_selection_table@Full[1,]))
)

DOWNSAMPLING_THRESHOLD <- 0.1
p_10_perc_reduction_2016 <- bs_calc_power(
  s="/global_workspace/imbcr_number_crunching/results/RawData_PLJV_IMBCR_20161201.csv",
  unmarked_model=unmarked_models,
  original_formulas=original_formulas,
  downsample_fun=downsample_stratified,
  top_model=as.numeric(row.names(model_selection_table@Full[1,]))
)

# for 2017 (loss of power from current dataset)

DOWNSAMPLING_THRESHOLD <- 0.3
p_30_perc_reduction_2017 <- bs_calc_power(
  s="/global_workspace/imbcr_number_crunching/results/RawData_PLJV_IMBCR_20171017.csv",
  unmarked_model=unmarked_models,
  original_formulas=original_formulas,
  downsample_fun=downsample_stratified,
  top_model=as.numeric(row.names(model_selection_table@Full[1,]))
)

DOWNSAMPLING_THRESHOLD <- 0.2
p_20_perc_reduction_2017 <- bs_calc_power(
  s="/global_workspace/imbcr_number_crunching/results/RawData_PLJV_IMBCR_20171017.csv",
  unmarked_model=unmarked_models,
  original_formulas=original_formulas,
  downsample_fun=downsample_stratified,
  top_model=as.numeric(row.names(model_selection_table@Full[1,]))
)

DOWNSAMPLING_THRESHOLD <- 0.1
p_10_perc_reduction_2017 <- bs_calc_power(
  s="/global_workspace/imbcr_number_crunching/results/RawData_PLJV_IMBCR_20171017.csv",
  unmarked_model=unmarked_models,
  original_formulas=original_formulas,
  downsample_fun=downsample_stratified,
  top_model=as.numeric(row.names(model_selection_table@Full[1,]))
)

# Alternative : for 2016 and 2017, what's the sample-size needed to capture a significant effect
# for an important variable -- in this case "shrub area"?

birds <- c("BRSP","CASP","GTTO","SATH") # Brewer's Sparrow, Cassin's Sparrow, Green-tailed Tohee, Sage Thrasher 

brsp_2016 <- calc_power_of_cov(
  four_letter_code=birds[1], 
  s="/global_workspace/imbcr_number_crunching/results/RawData_PLJV_IMBCR_20161201.csv", 
  sample_size=99, 
  units=units
)

brsp_2017 <- calc_power_of_cov(
  four_letter_code=birds[1], 
  s="/global_workspace/imbcr_number_crunching/results/RawData_PLJV_IMBCR_20171017.csv", 
  sample_size=99, 
  units=units
)

casp_2016 <- calc_power_of_cov(
  four_letter_code=birds[2], 
  s="/global_workspace/imbcr_number_crunching/results/RawData_PLJV_IMBCR_20161201.csv", 
  sample_size=99, 
  units=units
)

casp_2017 <- calc_power_of_cov(
  four_letter_code=birds[2], 
  s="/global_workspace/imbcr_number_crunching/results/RawData_PLJV_IMBCR_20171017.csv", 
  sample_size=99, 
  units=units
)

gtto_2016 <- calc_power_of_cov(
  four_letter_code=birds[3], 
  s="/global_workspace/imbcr_number_crunching/results/RawData_PLJV_IMBCR_20161201.csv", 
  sample_size=99, 
  units=units
)

gtto_2017 <- calc_power_of_cov(
  four_letter_code=birds[3], 
  s="/global_workspace/imbcr_number_crunching/results/RawData_PLJV_IMBCR_20171017.csv", 
  sample_size=99, 
  units=units
)

sath_2016 <- calc_power_of_cov(
  four_letter_code=birds[4], 
  s="/global_workspace/imbcr_number_crunching/results/RawData_PLJV_IMBCR_20161201.csv", 
  sample_size=99, 
  units=units
)

sath_2017 <- calc_power_of_cov(
  four_letter_code=birds[4], 
  s="/global_workspace/imbcr_number_crunching/results/RawData_PLJV_IMBCR_20171017.csv", 
  sample_size=99, 
  units=units
)

alternative_model_power_results <- rbind(
    data.frame(
      species="casp",
      density=median( mean(na.rm=T, casp_2016$density), mean(na.rm=T, casp_2017$density) ),
      se=median( mean(na.rm=T, casp_2016$se), mean(na.rm=T, casp_2017$se) ),
      z_score=median( ( mean(na.rm=T, casp_2016$density) / mean(na.rm=T, casp_2016$se) ) - (1.96/2), ( mean(na.rm=T, casp_2017$density) / mean(na.rm=T, casp_2017$se) ) - (1.96/2) ),
      power=1 - dnorm( median( ( mean(na.rm=T, casp_2016$density) / mean(na.rm=T, casp_2016$se) ) - (1.96/2), ( mean(na.rm=T, casp_2017$density) /   mean(na.rm=T, casp_2017$se) ) - (1.96/2) ) )
    ),
    data.frame(
      species="gtto",
      density=median( mean(na.rm=T, gtto_2016$density), mean(na.rm=T, gtto_2017$density) ),
      se=median( mean(na.rm=T, gtto_2016$se), mean(na.rm=T, gtto_2017$se) ),
      z_score=median( ( mean(na.rm=T, gtto_2016$density) / mean(na.rm=T, gtto_2016$se) ) - (1.96/2), ( mean(na.rm=T, gtto_2017$density) / mean(na.rm=T, gtto_2017$se) ) - (1.96/2) ),
      power=1 - dnorm( median( ( mean(na.rm=T, gtto_2016$density) / mean(na.rm=T, gtto_2016$se) ) - (1.96/2), ( mean(na.rm=T, gtto_2017$density) / mean(na.rm=T, gtto_2017$se) ) - (1.96/2) ) )
    ),
    data.frame(
      species="brsp",
      density=median( mean(na.rm=T, brsp_2016$density), mean(na.rm=T, brsp_2017$density) ),
      se=median( mean(na.rm=T, brsp_2016$se), mean(na.rm=T, brsp_2017$se) ),
      z_score=median( ( mean(na.rm=T, brsp_2016$density) / mean(na.rm=T, brsp_2016$se) ) - (1.96/2), ( mean(na.rm=T, brsp_2017$density) / mean(na.rm=T, brsp_2017$se) ) - (1.96/2) ),
      power=1 - dnorm( median( ( mean(na.rm=T, brsp_2016$density) / mean(na.rm=T, brsp_2016$se) ) - (1.96/2), ( mean(na.rm=T, brsp_2017$density) / mean(na.rm=T, brsp_2017$se) ) - (1.96/2) ) )
    ),
    data.frame(
      species="sath",
      density=median( mean(na.rm=T, sath_2016$density), mean(na.rm=T, sath_2017$density) ),
      se=median( mean(na.rm=T, sath_2016$se), mean(na.rm=T, sath_2017$se) ),
      z_score=median( ( mean(na.rm=T, sath_2016$density) / mean(na.rm=T, sath_2016$se) ) - (1.96/2), ( mean(na.rm=T, sath_2017$density) / mean(na.rm=T, sath_2017$se) ) - (1.96/2) ),
      power=1 - dnorm( median( ( mean(na.rm=T, sath_2016$density) / mean(na.rm=T, sath_2016$se) ) - (1.96/2), ( mean(na.rm=T, sath_2017$density) / mean(na.rm=T, sath_2017$se) ) - (1.96/2) ) )
    )
)

save(
  list=ls(pattern="^p_"),
  file=paste(
    "/home/ktaylora/power_analysis_run_results.rdata",
    sep=""
  ),
  compress=T
)
