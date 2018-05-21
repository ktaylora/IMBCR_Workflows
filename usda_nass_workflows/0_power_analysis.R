# load a previously built gdistsamp workflow for us to use for re-fitting our top models
# and then do a power analysis using bootstrapping to assess the potential
# influence of downsampling

load(commandArgs(trailingOnly=T))

# define our global constants
AIC_RESCALE_CONST         <- 100000
AIC_SUBSTANTIAL_THRESHOLD <- 8
CORRELATION_THRESHOLD     <- 0.65
AREA_OF_PLJV_KM           <- 646895.6
DOWNSAMPLING_THRESHOLD    <- 0.3
N_BS_REPLICATES           <- 999

#
# Local Functions
#

backscale_var <- function(var=NULL, df=NULL, m_scale=NULL){
  return(
      df[, var] *
      attr(m_scale, 'scaled:scale')[var] +
      attr(m_scale, 'scaled:center')[var]
    )
}

fit_gdistsamp <- function(lambdas=NULL, umdf=NULL, mixture="P"){
  if(length(lambdas)>1){
    cl <- parallel::makeCluster(parallel::detectCores()-1)
    parallel::clusterExport(cl, varlist=c("umdf","mixture"))
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
#' randomly downsample a dataset N=replicates times and calculate the density
#' and standard error for each run. Will return the full table of all replicate
#' runs (with statistics)
bs_calc_power <- function(
  replicates=N_BS_REPLICATES,
  s="/global_workspace/imbcr_number_crunching/results/RawData_PLJV_IMBCR_20171017.csv"
){
  # read-in our IMBCR transect data
  s <- OpenIMBCR:::scrub_imbcr_df(
    OpenIMBCR:::imbcrTableToShapefile(s),
    four_letter_code = toupper(argv[2])
  )
  # calculate our detections and sampling effort
  detections <- OpenIMBCR:::calc_dist_bins(s)
  effort     <- as.vector(OpenIMBCR:::calc_transect_effort(s))

  s <- calc_transect_summary_detections(
      s=s,
      name=toupper(argv[2]),
      field='est_abund'
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
      "top_model","fit_gdistsamp","DOWNSAMPLING_THRESHOLD",
      "formula"
    ),
    envir=environment()
  )
  replicates <- do.call(rbind, parallel::parLapply(
    cl=cl,
    X=1:replicates,
    fun=function(x){
      require(unmarked)
      # determine rows to keep that satisfy our DOWNSAMPLING_THRESHOLD
      sample <- sample(1:nrow(s@data), size=nrow(s@data)*(1-DOWNSAMPLING_THRESHOLD))
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
      # re-fit our top model using our downsampled dataset
      m <- fit_gdistsamp(
        lambda=formula,
        umdf=umdf,
        mixture=unmarked_models[[top_model]]@mixture
      )
      m_downsampled <- fit_gdistsamp(
          formula,
          umdf=umdf_downsampled,
          mixture=unmarked_models[[top_model]]@mixture
      )
      # calculate density and standard error
      density <- se <-
        unmarked::predict(
          m,
          type="lambda"
      )[,1:2]
      density <- mean(density[,1])
      se <-  mean(se[,2])
      aic <- m@AIC
      # now for our downsampled dataset
      density_downsampled <- se_downsampled <-
        unmarked::predict(
          m_downsampled,
          type="lambda"
      )[,1:2]
      density_downsampled <- mean(density_downsampled[,1])
      se_downsampled <-  mean(se_downsampled[,2])
      aic_downsampled <- m_downsampled@AIC

      return(data.frame(
        density=density,
        density_downsampled=density_downsampled,
        se=se,
        se_downsampled=se_downsampled,
        aic=aic,
        aic_downsampled=aic_downsampled
      ))
    }
  ));
  return(replicates)
}

# for 2016
n_detections_in_alternative_sample <- round(mean(sapply(
  X=1:N_BS_REPLICATES,
  FUN=function(x){
    sum(
      sample(
        rowSums(detections$y), size=DOWNSAMPLING_THRESHOLD*nrow(detections$y))
    )
  }
)))

# for 2017
