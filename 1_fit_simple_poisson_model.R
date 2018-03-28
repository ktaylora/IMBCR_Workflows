options(warn = -1, error=traceback)

require(OpenIMBCR)
require(glmulti)
require(unmarked)

AIC_RESCALE_CONST         <- 100000
AIC_SUBSTANTIAL_THRESHOLD <- 8
CORRELATION_THRESHOLD     <- 0.65

calc_transect_summary_detections <- function(s=NULL, name=NULL, field='est_abund'){
  return(OpenIMBCR:::calc_route_centroids(
      s = s,
      four_letter_code = toupper(name),
      use='est_abund'
    ))
}

quadratics_to_keep <- function(m){
  quadratic_terms <- grepl(names(coefficients(m)), pattern = "[)]2")
  # test : are we negative and are we a quadratic term?
  keep <- (coefficients(m) < 0) * quadratic_terms
  quadratic_terms <- names(coefficients(m))[keep == 1]
    quadratic_terms <- gsub(quadratic_terms, pattern = "[)]2", replacement = ")")
  # test: are we a positive linear term and a negative quadratic
  direction_of_coeffs <- coefficients(m)/abs(coefficients(m))
  steps <- seq(2, length(direction_of_coeffs), by = 2)
  keep <- names(which(direction_of_coeffs[steps] + direction_of_coeffs[steps+1]  == 0))
    quadratic_terms <- gsub(keep, pattern = "[)]1", replacement = ")")
  # test are both our linear and quadratic terms negative? drop if so
  if (length(quadratic_terms) > 0){
    return(quadratic_terms)
  } else {
    return(NULL)
  }
}

pred_hn_det_from_distance <- function(x=NULL, dist=NULL){
  param <- exp(coef(x, type = "det"))
  return(as.vector(unmarked:::gxhn(x=dist, param)))
}

calc_emp_dispersion_statistic <- function(x = NULL, bs=999999){
  observed <- sd(x)/round(mean(x))
  predicted <- median(sapply(
    X=bs,
    FUN=function(i){
      predicted <- rpois(n = length(x), round(mean(x)))
      return( sd(predicted)/round(mean(predicted)) )
    }))
  return(round(observed/predicted, 2))
}

calc_intercept_statistics <- function(x=NULL){
  load(x); # from an rdata file
  PLJV_AREA = 646191.287039 # in square kilometers
  return(data.frame(
    spp=strsplit(r_data_file, split="_")[[1]][1],
    mean_p_det=mean(per_obs_det_probabilities, na.rm=T),
    upper_pred=max(unmarked::predict(intercept_m, type="state")[,1]),
    lower_pred=min(unmarked::predict(intercept_m, type="state")[,1]),
    mean_pred=mean(unmarked::predict(intercept_m, type="state")[,1]),
    median_pred=median(unmarked::predict(intercept_m, type="state")[,1]),
    max_detections=max(rowSums(detections$y)),
    n_hat=mean(unmarked::predict(intercept_m, type="state")[,1])*PLJV_AREA
  ))
}
#' wrapper function for calc_intercept_statistics that will accept a
#' vector of rdata filenames and build a summary table of results for
#' all birds
calc_descriptive_statistics_rdata_files <- function(x=NULL){
  cat(" -- processing:")
  descriptive_statistics <- do.call(rbind, lapply(
    X=r_data_files,
    FUN=function(x){
      cat(paste("[",which(r_data_files %in% x),"]",sep=""))
      return(
        calc_intercept_descriptive_statistics(x)
      )
    }
  ))
  cat("\n")

  write.csv(
      descriptive_statistics,
      "intercept_model_descriptive_statistics.csv",
      row.names=F
    )
}

refit_glm <- function(m, s=NULL, effort=NULL){
  glm(formula=m$formula,
      offset=log(effort),
      family=poisson,
      data=s@data)
}

#
# MAIN
#

# Read-in our IMBCR transect data

argv <- commandArgs(trailingOnly = T)

cat(" -- fitting a model for :", argv[2], "\n")

r_data_file <- tolower(paste(
      tolower(argv[2]),
      "_imbcr_pois_glm_workflow_",
      gsub(format(Sys.time(), "%b %d %Y"), pattern = " ", replacement = "_"),
      ".rdata",
      sep = ""
    ))

# Do we have lurking output in the CWD?
if (file.exists(r_data_file)) {
   cat(" -- found existing rdata file:", r_data_file, "(skipping)\n");
   stop();
}

s <- OpenIMBCR:::scrub_imbcr_df(
    OpenIMBCR:::imbcrTableToShapefile(
        "/global_workspace/imbcr_number_crunching/results/RawData_PLJV_IMBCR_20161201.csv"
      ),
    four_letter_code = toupper(argv[2])
  )

detections <- OpenIMBCR:::calc_dist_bins(s)
effort     <- as.vector(OpenIMBCR:::calc_transect_effort(s))

if(sum(s$radialdistance>0, na.rm=T) < 180){
  cat(" -- not enough detection to fit a model (quitting)\n")
  quit("no")
}

# fit an intercept-only detection function in unmarked

umdf <- unmarked::unmarkedFrameDS(
    y=as.matrix(detections$y),
    siteCovs=data.frame(effort=effort),
    dist.breaks=detections$breaks,
    survey="point",
    unitsIn="m"
  )

intercept_m <- unmarked::distsamp(
    formula = ~1 ~1+offset(log(effort)),
    data = umdf,
    se = T,
    keyfun = "halfnorm",
    unitsOut = "kmsq",
    output = "abund"
  )

per_obs_det_probabilities <- round(sapply(
    s$radialdistance,
    function(x) pred_hn_det_from_distance(intercept_m, dist=x)),
    2
  )


# calculate an estimate of abundance (accounting for p-det)
s$est_abund <- round(s$cl_count > 0 / per_obs_det_probabilities)

s <- calc_transect_summary_detections(
    s=s,
    name=toupper(argv[2]),
    field='est_abund'
  )

# keep a local copy around so we don't lose it while re-scaling
est_abund <- s$est_abund

# read-in habitat covariates
units <- OpenIMBCR:::readOGRfromPath(argv[1])

# spatial join of transect detection data with habitat covariates
# summarized by-unit
s <- OpenIMBCR:::spatial_join(s, units)

# add latitude and longitude
cat(" -- calculating spatial covariates\n")

coords <- as.data.frame(sp::spTransform(s,"+init=epsg:4326")@coords)
  colnames(coords) <- c("lon","lat")
    s@data <- cbind(s@data, coords)

coords <- as.data.frame(rgeos::gCentroid(sp::spTransform(units,"+init=epsg:4326"), byid=T)@coords)
  colnames(coords) <- c("lon","lat")
    units@data <- cbind(units@data, coords)

# define the covariates we are going to use in our analysis

vars <- c("grass_ar","shrub_ar","crp_ar","wetland_ar","pat_ct", "map", "mat")
#vars <- c("grass_ar","shrub_ar","crp_ar","wetland_ar","pat_ct", "lat", "lon")
#vars <- c("grass_ar","shrub_ar","crp_ar","wetland_ar","pat_ct")

# ensure a consistent scale for our input data (we will use this a lot)

s@data <- s@data[, vars]
  s@data <- s@data[, sapply(s@data[1,], FUN=is.numeric)]

m_scale <- scale(s@data)
s@data <- as.data.frame(scale(s@data))

s$effort <- effort

# make sure the scale of the units we are predicting into is consistent
# with the training data

units@data <- as.data.frame(
    scale(units@data[,vars], attr(m_scale, "scaled:center"), attr(m_scale, "scaled:scale"))
  )

# tack-on a fake effort variable for predict()
units$effort <- median(effort)

# drop strongly-correlated variables
x_cor_matrix <- cor(s@data[,vars])

cor_threshold <- abs(x_cor_matrix)>CORRELATION_THRESHOLD

passing <- vector()
failing <- vector()

for(i in 1:ncol(cor_threshold)){
  not_correlated <- sum(which(cor_threshold[,i])!=i) == 0
  if(not_correlated){
    passing <- append(passing, colnames(cor_threshold)[i])
  } else {
    failing <- append(failing, colnames(cor_threshold)[i])
  }
}

rm(cor_threshold);

if(length(failing)>0){
  cat("-- these variables dropped due to colinearity problems:", failing,"\n")
  dev.new();
  corrplot::corrplot(x_cor_matrix)
vars <- unlist(strsplit(readline("enter a comma-separated list of vars you want to keep:"), split=","))
} else {
  vars <- passing
}

# fit our full model
m <- glm(
    formula=paste(
      "est_abund~",
      paste(paste("poly(",vars,",2,raw=T)",sep=""), collapse="+",sep=""),
      "+offset(log(effort))",
      sep=""
    ),
    family=poisson,
    data=s@data
  )

# drop un-intuitive quadratics
quadratics <- quadratics_to_keep(m)

if(length(quadratics)>0){
  vars <- vars[!as.vector(sapply(
    vars,
    FUN=function(p){ sum(grepl(x=quadratics, pattern=p))>0
  }))]
  # use AIC to justify our proposed quadratic terms
  for(q in quadratics){
    lin_var <- gsub(
        q,
        pattern=", 2,",
        replacement=", 1,"
      )

    m_lin_var <- AIC(glm(
        formula=paste(
          "est_abund~",
          paste(
            c(paste("poly(", paste(vars, ", 1, raw=T)", sep=""), sep=""),
            lin_var,quadratics[!(quadratics %in% q)]), collapse="+"),
            "+offset(log(effort))",
            sep=""
        ),
        family=poisson,
        data=s@data
      ))+AIC_RESCALE_CONST
    m_quad_var <- AIC(glm(
        formula=paste(
          "est_abund~",
          paste(
            c(paste("poly(", paste(vars, ", 1, raw=T)", sep=""), sep=""), quadratics),
            collapse="+"
          ),
          "+offset(log(effort))",
          sep=""
        ),
        family=poisson,
        data=s@data
      ))+AIC_RESCALE_CONST
    # if we don't improve our AIC with the quadratic by at-least 8 aic units
    # (pretty substatial support), keep the linear version
    if( m_lin_var-m_quad_var < AIC_SUBSTANTIAL_THRESHOLD){
      quadratics <- quadratics[!(quadratics %in% q)]
      vars <- c(
        vars,
        gsub(
          gsub(lin_var, pattern="poly[(]", replacement=""),
          pattern=", [0-9], raw.*=*.T[)]",
          replacement="")
        )
    }
  }

  vars <- c(paste("poly(", paste(vars, ", 1, raw=T)", sep=""), sep=""), quadratics)

  # refit our full model minus un-intuitive quadratics
  m <- glm(
      formula=as.formula(paste(
        "est_abund~",
        paste(vars, collapse="+"),
        sep=""
      )),
      family=poisson,
      offset=log(effort),
      data=s@data
    )
}

# use model selection with interactions across our candidate variables
tests <- glmulti::glmulti(
  m,
  intercept=T,
  family=poisson,
  offset=log(effort),
  level=1,
  plotty=F)

# predict across our full run dataset
predicted <- units

# are there multiple top models to work from?
if(sum((tests@crits)-min(tests@crits) < AIC_SUBSTANTIAL_THRESHOLD) > 1){
  models <- which((tests@crits)-min(tests@crits) < AIC_SUBSTANTIAL_THRESHOLD)
  weights <- tests@crits[models]
    weights <- OpenIMBCR:::akaike_weights(weights)
  models <- lapply(
    X=models, 
    FUN=function(x){
      m <- tests@objects[[x]]
      return(refit_glm(m, s, effort))
  })
  models <- do.call(cbind, lapply(
    X=models, 
    FUN=function(x){
      m <- predict(
          x,
          offset=log(effort),
          family=poisson,
          type="response",
          newdata=units@data
      )
    }))
  predicted@data <- data.frame(
    pred=apply(models, MARGIN=1, FUN=weighted.mean, weights=weights)
  )
  rm(models)
# if there was only one top model, averaging won't work
} else {
  # re-fit our standard model
  tests <- 
  glm(formula=tests@objects[[1]]$formula,
      offset=log(effort),
      family=poisson,
      data=s@data)
  predicted@data <- data.frame(
    pred=as.vector(floor(suppressWarnings(predict(
      tests,
      newdata=units@data,
      offset=log(effort),
      type="response")))
    )
  )
}

if(inherits(tests, 'glm')){
  # censor any predictions greater than K (max)
  cat(
      " -- number of sites with prediction greater than predicted max(K):",
      sum(predicted$pred > max(round(predict(tests, type="response")))),
      "\n"
    )
  predicted$pred[( predicted$pred > max(round(predict(tests, type="response"))) )] <-
    max(round(predict(tests, type="response")))
} else {
  # censor any predictions greater than K (max)
  cat(
      " -- number of sites with prediction greater than predicted max(K):",
      sum(predicted$pred > max(round(predict(tests@objects[[1]], type="response")))),
      "\n"
    )
  predicted$pred[( predicted$pred > max(round(predict(tests@objects[[1]], type="response"))) )] <-
    max(round(predict(tests@objects[[1]], type="response")))
}


predicted$pred[predicted$pred<1] <- 0

rgdal::writeOGR(
  predicted,
  ".",
  tolower(paste(argv[2],
      "_imbcr_pois_glm_prediction_",
      gsub(format(Sys.time(), "%b %d %Y"), pattern=" ", replacement="_"),
      sep="")
  ),
  driver="ESRI Shapefile",
  overwrite=T
)

save(
    compress=T,
    list=ls(),
    file=r_data_file
  )
