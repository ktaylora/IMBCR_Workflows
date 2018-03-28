AIC_RESCALE_CONST         <- 100000
AIC_SUBSTANTIAL_THRESHOLD <- 8
CORRELATION_THRESHOLD     <- 0.65
AREA_OF_PLJV_KM           <- 646895.6

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

# for 2016
top_model <- as.numeric(row.names(model_selection_table@Full[1,]))
density_2016 <- mean(
    unmarked::predict(unmarked_models[[top_model]], type="lambda")[,1]
  )
se_2016 <- mean(
    unmarked::predict(unmarked_models[[top_model]], type="lambda")[,2]
  )
aic_2016 <- unmarked_models[[top_model]]@AIC

# for 2017

formula <- original_formulas[top_model]

s <- OpenIMBCR:::scrub_imbcr_df(
    OpenIMBCR:::imbcrTableToShapefile(
        "/global_workspace/imbcr_number_crunching/results/RawData_PLJV_IMBCR_20171017.csv"
      ),
    four_letter_code = toupper(argv[2])
  )

detections <- OpenIMBCR:::calc_dist_bins(s)
effort     <- as.vector(OpenIMBCR:::calc_transect_effort(s))

s <- calc_transect_summary_detections(
    s=s,
    name=toupper(argv[2]),
    field='est_abund'
  )

s <- OpenIMBCR:::spatial_join(s, units)

if(unmarked_models[[top_model]]@mixture == "NB"){
  umdf <- unmarked::unmarkedFrameGDS(
    y=as.matrix(detections$y),
    siteCovs=s@data[,c(vars,'effort')],
    dist.breaks=detections$breaks,
    numPrimary=1,
    survey="point",
    unitsIn="m"
  )
  m <- fit_gdistsamp(formula, umdf=umdf, mixture="NB")
  density_2017 <- se_2017 <- 
    unmarked::predict(
      m, 
      type="lambda"
  )[,1:2]     
  density_2017 <- mean(density_2017[,1])
  se_2017 <-  mean(se_2017[,2])
  aic_2017 <- m@AIC
} else {
  umdf <- unmarked::unmarkedFrameGDS(
    y=as.matrix(detections$y),
    siteCovs=s@data[,c(vars,'effort')],
    dist.breaks=detections$breaks,
    numPrimary=1,
    survey="point",
    unitsIn="m"
  )
  m <- fit_gdistsamp(formula, umdf=umdf, mixture="P") 
  density_2017 <- se_2017 <- 
    unmarked::predict(
      m, 
      type="lambda"
  )[,1:2]     
  density_2017 <- mean(density_2017[,1])
  se_2017 <-  mean(se_2017[,2])
  aic_2017 <- m@AIC
}

data.frame(
    species=argv[2],
    dens_2016=density_2016,
    se_2016=se_2016,
    n_2016=round(AREA_OF_PLJV_KM*density_2016),
    dens_2017=density_2017,
    se_2017=se_2017, 
    n_se_2017=round(AREA_OF_PLJV_KM*se_2017)
  )
