#
# Author: Kyle Taylor (2017) PLJV
#
# Accepts a single argument at runtime specifying the full path to
# an optionally unattributed USNG file and then does some inferential
# statistics with the Sohl (2017) LULC dataset
#

require(rgdal)
require(ggplot2)
require(OpenIMBCR)

#
# Runtime Parameters
#

AIC_RESCALE_CONST         <- 100000
AIC_SUBSTANTIAL_THRESHOLD <- 8
CORRELATION_THRESHOLD     <- 0.65
NORMALIZE                 <- T

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

ggplot2_univar_density <- function(s=NULL, var=NULL, correction=1, ylab=NULL, xlab=NULL){
  xlim <- range( backscale_var(var, s@data, m_scale)/correction )

  xlim[1] <- xlim[1] - xlim[1]*0.1
  xlim[2] <- xlim[2] + xlim[2]*0.1

  rug <- data.frame(x=backscale_var(var, s@data, m_scale) / correction)
  median <- median(rug$x)

  ggplot(
    data=data.frame(
        mat=backscale_var(var, s@data, m_scale)
      )
      , aes(mat)
    ) +
    scale_x_continuous(limits=xlim) +
    geom_density(fill='#edf8b1', alpha=0.6, color="#97a837") +
    geom_rug(aes(x), sides="b", color="#97a837", alpha=0.8, size=0.25, data=rug) +
    xlab(xlab) +
    ylab(ylab) +
    geom_vline(xintercept = median, color="#97a837", size=0.5) +
    scale_color_brewer(type = 'seq', palette="YlGnBu") +
    theme(
        axis.text = element_text(face="bold",size=10),
        text = element_text(face="bold",size=12),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        #panel.background = element_blank(),
        #panel.grid.minor = element_blank(),
        #panel.grid.major = element_line("white"),
        axis.line=element_line("grey50")
      )
}

pts_to_landcover_metrics <- function(
    pts=NULL,
    r=NULL,
    grid_units=NULL,
    composition_statistics=NULL,
    force_valid_habitat_val=NULL)
{
  # automatically calculate 'patch count' automatically using the
  # habitat covariates passed by the user as composition_statistics
  configuration_statistics <- c(
    'pat_ct'
  )
  # whare are the valid raster cell values associated for all of our
  # habitat composition metrics?
    if(is.null(force_valid_habitat_val)) {
    valid_habitat_values <- eval(parse(
      text=paste("c(",paste(composition_statistics$src_raster_value[
        !grepl(composition_statistics$field_name, pattern="rd_ar")
      ], collapse = ","), ")", sep="")
    ))
  } else {
    valid_habitat_values <- force_valid_habitat_val
  }

  # subset our units shapefile by overlapping points in s
  grid_units_over <- !is.na(sp::over(
      units, spTransform(s, sp::CRS(raster::projection(units))))[,1]
    )
  units <- units[grid_units_over,]
  # zero-out our attribute table
  units@data <- data.frame(id=1:nrow(units))
  # split for list comprehension and parallelization
  e_units <- split(units, 1:nrow(units))
  # do our extractions
  e_units <- OpenIMBCR:::extract_by(e_units, r)
  # sanity-check : do we have enough units to cover our features in s?
  if(nrow(s) != length(e_units)){
    stop("we are missing some overlapping USNG units for features in s")
  }
  # calculate our landscape metric (total area) across our input dataset
  metrics <- lapply(
    X=1:nrow(composition_statistics),
    FUN=function(i){
      OpenIMBCR:::par_calc_stat(
        # using our 1 km unit raster extractions
        X=e_units,
        fun = OpenIMBCR:::calc_total_area,
        # using these PLJV landcover cell values in the reclassification
        from = eval(parse(text=as.character(composition_statistics[i, 2])))
      )
    }
  )

  metrics <- as.data.frame(do.call(cbind, metrics))

  # calculate our landscape configuration metric (patch count
  metrics[, ncol(metrics)+1] <-
  OpenIMBCR::par_calc_stat(
      # using our using our un-buffered unit raster extractions
      e_units,
      # parse the focal landscape configuration metric
      fun = OpenIMBCR::calc_patch_count,
      # using these PLJV landcover cell values in the supplemental
      # reclassification
      from = valid_habitat_values
    )
  # settle on our column names
  colnames(metrics) <-
    c(as.vector(composition_statistics$field_name), configuration_statistics)
  # return table of units to user
  s@data <- metrics
  return(s)
}

plot_hn_det <- function(x=NULL, breaks=NULL, add_bins=T){
  param <- exp(unmarked::coef(x, type = "det"))
  # first plot an empty canvas to the right dimensions
  plot(
      function(x) gxhn(x, param), 0, max(breaks),
  	  xlab = "Distance (m)", ylab = "Detection probability",
  	  lwd=1.5,
  	  col="white"
    )
  # add bins if requested
  if(add_bins){
    intensities <- colSums(x@data@y)
      intensities <- intensities/max(intensities)
    rect(
      x@data@dist.breaks[-length(x@data@dist.breaks)], 
      0, 
      x@data@dist.breaks[-1], 
      intensities,
      col="DarkGrey",
      border="white"
    )
  }
  # now add a red line for our detection function
  plot(
      function(x) gxhn(x, param), 0, max(breaks),
  	  xlab = "Distance (m)", ylab = "Detection probability",
  	  lwd=1.5,
  	  col="red",
  	  add=T
    )
  grid(); grid();
}

quadratics_to_keep <- function(m){
  direction_of_coeffs <- unmarked::coef(m)/abs(unmarked::coef(m))
  quadratic_terms <- grepl(names(unmarked::coef(m)), pattern = "[)]2")
  # bug-fix: drop any alpha or p parameters, they mess things up
  is_lambda <- grepl(tolower(names(unmarked::coef(m))), pattern="lam")
  direction_of_coeffs <- na.omit(direction_of_coeffs[is_lambda])
  quadratic_terms <- quadratic_terms[is_lambda]
  # test : are we negative and are we a quadratic term?
  keep <- (direction_of_coeffs < 0) * quadratic_terms
  quadratic_terms <- names(direction_of_coeffs)[keep == 1]
  # no negative quadratics? then leave
  if(length(quadratic_terms)==0){
    return(NULL)
  # negative quadratic? let's make sure the linear term is positive
  } else {
    quadratic_terms <- gsub(quadratic_terms, pattern = "[)]2|[)]2.0", replacement = ")")
    quadratic_terms <- gsub(quadratic_terms, pattern = "lambda[(]|lam[(]", replacement="")
    quadratic_terms <- gsub(quadratic_terms, pattern = "[)][)]", replacement=")")
    # test: are we a positive linear term and a negative quadratic
    steps <- seq(2, length(direction_of_coeffs), by = 2) # always skip the intercept
    keep <- names(which(direction_of_coeffs[steps] + direction_of_coeffs[steps+1]  == 0))
    keep <- gsub(keep, pattern="[)].[)]", replacement=")")
    keep <- gsub(keep, pattern = "lambda[(]|lam[(]", replacement="")
    if(length(keep)==0){
      return(NULL)
    }
    # are our negative quadratic(s) in the "keep" array?
    quadratic_terms <-
      quadratic_terms[gsub(quadratic_terms, pattern=" ", replacement="") %in% gsub(keep, pattern=" ", replacement="")]
    # test are both our linear and quadratic terms negative? drop if so
    if (length(quadratic_terms) > 0){
      return(quadratic_terms)
    } else {
      return(NULL)
    }
  }
}

calc_all_distsamp_combinations <- function(vars = NULL, poly_order=T){
  formulas <- gsub(OpenIMBCR:::mCombinations(
          siteCovs = paste("poly(",vars,",2,raw=T)",sep=""),
          availCovs = NULL,
          detCovs = NULL,
          offset = "offset(log(effort))")[,1],
        pattern = " ~1 ~1",
      replacement = ""
    )
  formulas <- gsub(
      formulas,
      pattern = " ",
      replacement = ""
    )
  formulas <- gsub(
      formulas,
      pattern = "[+]offset[(]log[(]effort[)][)]",
      replacement = ""
    )
  formulas <- gsub(formulas, pattern = "~", replacement = "")
  return(formulas)
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
#' this function has a ridiculous amount of complexity built into it and needs
#' to be refactored. But building these tests into simpler function(s) will
#' take some thinking.
aic_test_quadratic_terms_gdistsamp <- function(unmarked_models=NULL, original_formulas=NULL, umdf=NULL, mixture="P"){
  # here's out built-in aic threshold method that needs to be refactored
  aic_threshold_test <- function(i=NULL, quadratics=NULL, vars=NULL){
      # drop the lam() prefix
      quads <- gsub(gsub(quadratics[[i]], pattern="lambda[(]|lam[(]", replacement=""), pattern="[)][)]", replacement=")")
      v <- vars[i]
      # drop poly() notation from the list of all covariates using in this model
      v <- gsub(gsub(v, pattern="poly[(]", replacement=""), pattern="[,][0-2][,]raw=T[)]", replacement="")
      v <- unlist(strsplit(v, split="[+]"))
      if(length(quads)>0){
        # drops our linear variable from consideration in the quadratics list
        v <- v[!as.vector(sapply(v, FUN=function(p=NULL){ sum(grepl(x=quads, pattern=p))>0  }))]
        # use AIC to justify our proposed quadratic terms
        for(q in quads){
          lin_var <- gsub(
            q,
            pattern=",2",
            replacement=",1"
          )
          # lambda formula, with all variables (INCLUDING the focal variable)
          # specified as poly(var,1)
          lambda_formula <- ifelse(
            length(v)==0,
            # empty v?
            paste(
              "~",
              paste(
                c(
                  lin_var,
                  gsub(quads[!(quads %in% q)], pattern="2", replacement="1")
                ),
                collapse="+"
              ),
              "+offset(log(effort))",
              sep=""
            ),
            # valid v?
            paste(
              "~",
              paste(
                c(
                  paste("poly(", paste(v, ", 1, raw=T)", sep=""), sep=""),
                  lin_var,
                  gsub(quads[!(quads %in% q)], pattern="2", replacement="1")
                ),
                collapse="+"
              ),
              "+offset(log(effort))",
              sep=""
            )
          )
          m_lin_var <- try(OpenIMBCR:::AIC(unmarked::gdistsamp(
            pformula=as.formula("~1"),
            lambdaformula=as.formula(lambda_formula),
            phiformula=as.formula("~1"),
            data=umdf,
            se=T,
            keyfun="halfnorm",
            mixture=mixture,
            unitsOut="kmsq",
            output="abund"
          )) + AIC_RESCALE_CONST)
          if(class(m_lin_var) == "try-error"){
            warning(
              "we failed to converge on a solution when fitting the linear model:",
              lambda_formula
              )
            return(NA)
          }
          # lambda formula, with all variables (EXCEPT the focal variable)
          # specified as poly(var,1)
          lambda_formula <- ifelse(
            length(v)==0,
            # empty v?
            paste(
              "~",
              paste(
                c(
                   gsub(lin_var, pattern="1", replacement="2"),
                   gsub(quads[!(quads %in% q)], pattern="2", replacement="1")
                 ),
                 collapse="+"
             ),
              "+offset(log(effort))",
              sep=""
            ),
            # valid v?
            paste(
              "~",
              paste(
                c(
                  paste("poly(", paste(v, ", 1, raw=T)", sep=""), sep=""),
                  gsub(lin_var, pattern="1", replacement="2"),
                  gsub(quads[!(quads %in% q)], pattern="2", replacement="1")
                ),
                collapse="+"
              ),
              "+offset(log(effort))",
              sep=""
            )
          )
          m_quad_var <- try(OpenIMBCR:::AIC(unmarked::gdistsamp(
            pformula=as.formula("~1"),
            lambdaformula=as.formula(lambda_formula),
            phiformula=as.formula("~1"),
            data=umdf,
            se=T,
            keyfun="halfnorm",
            mixture=mixture,
            unitsOut="kmsq",
            output="abund"
          )) + AIC_RESCALE_CONST)
          if(class(m_quad_var) == "try-error"){
            warning(
              "we failed to converge on a solution when fitting the quadratic model:",
              lambda_formula
            )
            return(NA)
          }
          # if we don't improve our AIC with the quadratic by at-least 8 aic units
          # (pretty substatial support), keep the linear version
          if( m_lin_var-m_quad_var < AIC_SUBSTANTIAL_THRESHOLD ){
            quads <- quads[!(quads %in% q)]
            v <- c(
              v,
              gsub(
                x=gsub(lin_var, pattern="poly[(]", replacement=""),
                pattern=",[0-9],raw.*=*.T[)]",
                replacement=""
              )
            )
          }
        }
        v <- c(paste("poly(", paste(v, ",1,raw=T)", sep=""), sep=""), quads)
      } else {
        # if there were no valid quadratics to test, we will default to using
        # the linear form only
        v <- c(paste("poly(", paste(v, ",1,raw=T)", sep=""), sep=""), quads)
      }
      return(v)
  }
  # here is where we actually call our test
  cl <- parallel::makeCluster(parallel::detectCores()-1)
  # determine our run-time parameters
  quadratics <- lapply(unmarked_models, FUN=quadratics_to_keep)
        vars <- original_formulas
  # remove any padding around our vars (messes with string regular expressions)
  quadratics <- lapply(quadratics, FUN=function(x) gsub(x, pattern=" ", replacement=""))
  vars <- sapply(vars, FUN=function(x) gsub(vars, pattern=" ", replacement=""))
  # set-up our run and parallelize across our cores
  #parallel::clusterExport(cl, varlist=c("AIC_RESCALE_CONST", "AIC_SUBSTANTIAL_THRESHOLD"), envir=globalenv())
  #parallel::clusterExport(cl, varlist=c("umdf","aic_threshold_test"), envir=environment())
  #vars <- parallel::parLapply(
  #  cl=cl,
  #  X=1:length(original_formulas),
  #  fun=aic_threshold_test,
  #  quadratics=quadratics,
  #  vars=vars
  #)
  #parallel::stopCluster(cl);
  vars <- sapply(
    X=1:length(original_formulas),
    FUN=function(x){
      cat(".")
      aic_threshold_test(i=x, quadratics=quadratics, vars=vars)
    }
  )
  cat("\n")
  # sanity check -- never return a "poly(,[0-9],raw=T)" pattern
  # if this occurs, return the original formula with the linear term specified
  problems <- which(grepl(vars, pattern="poly[(],"))
  if( sum(problems) > 0 ){
    warning("we encountered problems finding polynomical terms for vars array:",vars)
    vars[problems] <- gsub(original_formulas[problems], pattern="2", replacement="1")
  }
  return(vars);
}

par_unmarked_predict <- function(unmarked_models=NULL, predict_df=NULL, type="lambda", weights=NULL){
  # set-up our run and parallelize across our cores
  cl <- parallel::makeCluster(parallel::detectCores()-1)
  steps <- round(seq(1,nrow(predict_df), length.out=parallel::detectCores()-1))
  # add 1 to the last step to accomodate our lapply splitting
  steps[length(steps)] <- steps[length(steps)]+1
  parallel::clusterExport(cl, varlist=c("predict_df","steps","type"), envir=environment())
  predicted <- lapply(
   X=unmarked_models,
   FUN=function(model){
     parallel::clusterExport(cl, varlist=c("model"), envir=environment())
     return(parallel::parLapply(
       cl=cl,
       X=1:(length(steps)-1),
       fun=function(i){
         return(unmarked::predict(
             model,
             type=type,
             se=T,
             newdata=predict_df[seq(steps[i], (steps[i+1]-1)),]
           ))
       }))
     }
  )
  rm(predict_df);
  parallel::stopCluster(cl); rm(cl);
  # join the individual rows from within each model run
  # and return a single data.frame for each model
  predicted <- lapply(
     predicted,
     FUN=function(x) do.call(rbind, x)
   )
   # each model run will have a predicted column and
   # a confidence interval -- we are only interested in
   # the predicted column right now
   predicted <- lapply(
       X=1:length(unmarked_models),
       FUN=function(x){ predicted[[x]]$Predicted }
   )
   if(is.null(weights)){
     return(predicted)
   } else {
     # if we have more than one model in the top models
     # table, let's average the results across our models
     # using AIC weighting parameter taken from 'unmarked'.
     if (length(predicted) > 1){
       # join our predictions across models into a single
       # matrix that we can apply a weight across
       predicted <- do.call(cbind, predicted)
       predicted <- sapply(
           1:nrow(predicted),
           FUN=function(i){
             weighted.mean(x=predicted[i, ], w = weights)
           }
         )
     } else {
       predicted <- unlist(predicted)
     }
     return(predicted)
   }
}
#' testing: standard PCA reconstruction that, for each variable, will find
#' the principal component that captures the greatest variance and then 
#' partition-out that component and use it's variance to reconstruct the
#' original covariate
pca_partial_reconstruction <- function(df=NULL, vars=NULL){
  # by default, accept scaled covariates
  m_pca <- prcomp(df[,vars])
  partialed_covs <- df[,vars]
  remaining <- 1:length(vars) # make sure we never use the same component for more than one variable
  for(var in vars){
    if(length(remaining)>1){
      col <- which.max(abs(m_pca$rotation[var,]))
      remaining <- remaining[ remaining != col ]
    } else {
      # if we only have one remaining rotation to use, use it (even 
      # if it's not the best fit)
      col <- remaining
    }
    x_hat <- m_pca$x[,col] %*% t(m_pca$rotation[,col])
      x_hat <- x_hat[,col] # retain only our partial mean for THIS component
    # make sure the sign matches our original cov
    if ( cor(x_hat, df[,var]) < 0 ){
     x_hat <- -1 * x_hat
    }
    # re-scale to the max of our original input dataset  
    x_hat <- ( x_hat - min(x_hat) ) / ( max(x_hat) - min(x_hat) )  * max(df[,var])
    # store our partialed covariate
    partialed_covs[,var] <- x_hat
  }
  return(partialed_covs)
}
#' testing : drop covariates of middling importance and only retain the best
#' and worst axes.
pca_partial_reconstruction_without_middle <- function(df=NULL, vars=NULL){
  # by default, accept un-scaled covariates -- we will do the scaling
  # and then return components in the raw scale of the input data
  m_pca <- prcomp(df[,vars])
  partialed_covs <- df[,vars]
  for(var in vars){
    col <- as.vector(c( 
        which.max(abs(m_pca$rotation[var,])), 
        which.min(abs(m_pca$rotation[var,])) 
      ))
    x_hat <- m_pca$x[,col] %*% t(m_pca$rotation[,col])
      x_hat <- x_hat[,var] # retain only our variance for our focal variable
    # re-scale to the max of our original input dataset  
    x_hat <- ( x_hat - min(x_hat) ) / ( max(x_hat) - min(x_hat) )  * max(df[,var])
    #x_hat <- scale(x_hat, center = mean(df[,var]), scale = F)
    # make sure the sign matches our original cov
    #x_hat <- x_hat * as.vector( cor(x_hat, df[,var])/abs(cor(x_hat, df[,var])) )
    # store our partialed covariate
    partialed_covs[,var] <- x_hat
  }
  return(partialed_covs)
}

#
# Main
#

argv <- commandArgs(trailingOnly = T)

# let's use some sane defaults just in-case the user didn't specify
# anything at runtime
if(length(argv) != 2){
  argv <- argv[1] # always assume we at-least offered up a bird code
  argv[2] <- "/global_workspace/iplan_imbcr_workspace/vector/units_attributed_nass_crp_2016_training_1km.shp"
  argv[3] <- "/global_workspace/terry_sohl_sres_btu_30m_landcover_predictions/extracted/gplcc_gcam_45_2014.tif"
}

# define the covariates we are going to use for model fitting
#vars <- c("grass_ar","shrub_ar","wetland_ar","pat_ct", "mat", "map")
#vars <- c("grass_ar","shrub_ar","mat", "map")
vars <- c("grass_ar","shrub_ar","pat_ct")

cat(" -- fitting a HDS model for :", argv[1], "\n")

r_data_file <- tolower(paste(
      tolower(argv[1]),
      "_imbcr_hds_iale_workflow_",
      gsub(format(Sys.time(), "%b %d %Y"), pattern = " ", replacement = "_"),
      ".rdata",
      sep = ""
    ))

s <- OpenIMBCR:::scrub_imbcr_df(
    OpenIMBCR:::imbcrTableToShapefile(
        "/global_workspace/imbcr_number_crunching/results/RawData_PLJV_IMBCR_20161201.csv"
        #"/global_workspace/imbcr_number_crunching/results/RawData_PLJV_IMBCR_20171017.csv"
      ),
    four_letter_code = toupper(argv[1])
  )

detections <- OpenIMBCR:::calc_dist_bins(s)
effort     <- as.vector(OpenIMBCR:::calc_transect_effort(s))
transects  <- as.character(unique(s$transectnum))

# collapse our thousands of station-level points into transect-level centroids
s <- OpenIMBCR:::calc_route_centroids(
        s = s,
        four_letter_code = toupper(argv[1])
  )

cat(
    " -- reading-in our USNG units shapefile to use for attributing our land",
    "cover rasters\n"
  )

units <- OpenIMBCR:::readOGRfromPath(argv[2])
# assume we have pre-calculated long-term climate covs, drop anything else
units@data <- units@data[, c('mat','map')]

# This is a standard USDA NASS (2016) dataset
#r <- raster::raster(paste("/gis_data/Landcover/NASS/Raster/",
#    "2016_nass_crp_test_merge.tif", sep=""
#  ))
cat(" -- reading in source raster data\n")
# This is the baseline scenario for GCAM v.4.5 (RCP 4.5) that we are going
# to peg our IMBCR models to
r_gcam_45_2014 <- raster::raster(
    argv[3]
  )

# Build a specification for the land cover raster cell values and categories
# we are going to use (configuration [patch count] is handled explicitly
# by our local function def (above)
composition_statistics <-
  data.frame(
      field_name=c(
        'grass_ar',
        'shrub_ar'
        #'wetland_ar'
      ),
      src_raster_value=c(
        'c(19,26)',
        '27'
        #'c(28,29)'
      )
  )

cat(
    " -- calculating landscape metrics for USNG units that overlap",
    "IMBCR transects\n"
  )

# Attribute our points using landscape metrics calculated from our transect
# points and whatever raster data occurs under a point's corresponding USNG
# unit
s <- pts_to_landcover_metrics(
       pts=s,
       grid_units=units,
       r=r_gcam_45_2014,
       composition_statistics=composition_statistics,
       force_valid_habitat_val=c(19,26) # only use 'grass' for our config. metric
     )

# there are pre-calculated variables in the units file that we want to keep
# (for climate conditions)
s <- s_original <- OpenIMBCR:::spatial_join(s, units)

# do some pca reconstruction to uncorrelate our variables
s@data <- s@data[, vars]
s@data[,vars] <- pca_partial_reconstruction(
    s_original@data, 
    c(vars, "map", "mat")
  )[, vars]

# ensure a consistent scale for our input data (we will use this a lot)
s@data <- s@data[, sapply(s@data[1,], FUN=is.numeric)]
m_scale <- scale(s@data)
s@data <- as.data.frame(scale(s@data))

# now tack on our transect sampling effort (don't scale this)
s$effort <- effort


#
# Build some bird models
#

cat(
    " -- fitting a first round of negbin models and testing quadratic terms",
    "(this may take 40-80 minutes)\n")

# build a unmarked data.frame from our running
# input data (s@data is attributed IMBCR grid centroids)

umdf <- unmarked::unmarkedFrameGDS(
  y=as.matrix(detections$y),
  siteCovs=s@data,
  dist.breaks=detections$breaks,
  numPrimary=1,
  survey="point",
  unitsIn="m"
)


# calculate exhaustive (all possible) variable mCombinations
# for model selection

original_formulas <- unmarked_models <- calc_all_distsamp_combinations(vars)
unmarked_models <- paste(unmarked_models, "+offset(log(effort))", sep="")

#
# Let's get a few intercept-only density estimates to compare against our
# predictions with covariates on lambda
#

cat(" -- fitting intercept models for evaluation:\n")

m_negbin_intercept <- fit_gdistsamp(
  "1+offset(log(effort))",
  umdf=umdf,
  mixture="NB"
)

m_pois_intercept <- fit_gdistsamp(
  "1+offset(log(effort))",
  umdf=umdf,
  mixture="P"
)

m_negbin_intercept_n <- unmarked::backTransform(
  m_negbin_intercept, 
  type="lambda")@estimate

m_pois_intercept_n <- unmarked::backTransform(
  m_pois_intercept, 
  type="lambda")@estimate

cat(" -- determining whether we MUST use the negative binomial mixture\n")


# make an over-fit model of all variables to stare at
# and wonder
m_negbin_full_model <-fit_gdistsamp(
  paste(
      paste(paste("poly(",vars,",2,raw=T)",sep=""), collapse="+"),
      "+offset(log(effort))",
      sep=""
    ),
  mixture="NB",
  umdf=umdf
)

m_pois_full_model <- fit_gdistsamp(
  paste(
      paste(paste("poly(",vars,",2,raw=T)",sep=""), collapse="+"),
      "+offset(log(effort))",
      sep=""
    ),
  mixture="P",
  umdf=umdf
)

# favor the poisson mixture -- but if our data have high variance,
# it may just fail to converge. In that case, use the negative binomial
if(class(m_pois_full_model) == "logical"){
  mixture <- "NB"
} else {
  mixture <- "P"
}

cat(" -- testing linear / quadratic terms using AIC:")

# takes 45-90 minutes -- note that we have offsets here
models <- fit_gdistsamp(unmarked_models, umdf, mixture=mixture)

# sometimes models fail to converge -- they are coded as NA's

original_formulas <- original_formulas[
    suppressWarnings( !sapply(models, FUN=function(x) sum(is.na(x)) >0 ) )
  ]

models <- suppressWarnings(
    models[ !sapply(models, FUN=function(x) sum(is.na(x)) >0 ) ]
  )

# aic_test will add an offset term onto our model formulas here
unmarked_models <- aic_test_quadratic_terms_gdistsamp(
  unmarked_models=models,
  original_formulas=original_formulas,
  umdf=umdf,
  mixture=mixture
)

# refit our models using the specification justified from testing AIC
# across linear vs quadratic terms -- also, add the offset term back onto
# the model formula
unmarked_models <- fit_gdistsamp(
  lambdas=lapply(
    X=unmarked_models,
    FUN=function(x){
      if (!grepl(paste(x, collapse="+"), pattern="offset")){
        return(paste(paste(x, collapse="+"), "+offset(log(effort))", sep=""))
      } else {
        return(paste(x, collapse="+"))
      }
    }),
  umdf=umdf,
  mixture=mixture
)

cat(" -- performing model averaging tasks\n")

# make a fitList
model_selection_table <- suppressWarnings(unmarked::modSel(unmarked::fitList(
  fits=unmarked_models
)))

MOD_SEL_THRESHOLD <- max(which(model_selection_table@Full$delta < AIC_SUBSTANTIAL_THRESHOLD))
# select the top models (by array position) that satisfy our threshold
MOD_SEL_THRESHOLD <- as.numeric(model_selection_table@Full$model)[1:MOD_SEL_THRESHOLD]

# copy our input training data over to a predict data.frame for our
# n_hat calculation

predict_df <- s@data

aic_weights <- model_selection_table@Full$AICwt[
    MOD_SEL_THRESHOLD
]

predicted <- par_unmarked_predict(
  unmarked_models=unmarked_models[MOD_SEL_THRESHOLD],
  predict_df=predict_df,
  type="lambda",
  weights=aic_weights
)

if(NORMALIZE){
  predicted <- (predicted/mean(predicted)) * 
    min(c(m_pois_intercept_n, m_negbin_intercept_n))
}

n_hat <- mean(predicted)
n_hat_sd <- sd(predicted)

rm(predicted)

cat(" -- writing to disk\n")


# make some plots
#ggplot2_univar_density(var="shrub_ar", xlab="Total Area Shrubland (km2)")

save(
    compress=T,
    list=(ls(pattern="[a-z]")),
    file=r_data_file
  )
