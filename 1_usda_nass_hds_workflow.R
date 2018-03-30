#
# This is not pretty -- I'm sorry; Kyle (2018)
#
# Accepts a single argument at run-time: the full path to an
# RData file produced from a glm workflow. This has all of the
# transect-level covariate data and detections pre-defined and
# cut's a few hundred lines of code.
#
# This workflow will produce a negative-binomial HDS model from
# unmarked (gdistsamp), do model averaging, and make a shapefile
# of predictions.
#

load(commandArgs(trailingOnly = T)[1])
# consider using options(error=traceback)
# options(warn = -1, error=traceback)

if (file.exists(gsub(
    r_data_file,
    pattern = "pois_glm",
    replacement = "hds_negbin"
  ))){
  cat("previous workspace for this bird found in CWD (skipping)\n")
  q("no")
}

AIC_RESCALE_CONST         <- 100000
AIC_SUBSTANTIAL_THRESHOLD <- 8
CORRELATION_THRESHOLD     <- 0.65

#
# Local accessory functions, some of which may overload what's
# saved in the loaded Rdata file
#

backscale_var <- function(var=NULL, df=NULL, m_scale=NULL){
  return(
      df[, var] *
      attr(m_scale, 'scaled:scale')[var] +
      attr(m_scale, 'scaled:center')[var]
    )
}

plot_hn_det <- function(x=NULL, breaks=NULL){
  param <- exp(unmarked::coef(x, type = "det"))
  plot(
      function(x) gxhn(x, param), 0, max(breaks),
  	  xlab = "Distance (m)", ylab = "Detection probability"
    )
  grid(); grid();
}

plot_model_pi <- function(tests = NULL, unmarked_m = NULL){
  hinge_m_pred <- density(as.vector(floor(
      predict(tests@objects[[1]], type = "response")))
    )
  hds_m_pred <- density(as.vector(floor(
      unmarked::predict(unmarked_m, type = "state")[, 1]))
    )

  xlim <- range(c(hinge_m_pred$x, hds_m_pred$x))
  ylim <- range(c(hinge_m_pred$y, hds_m_pred$y))

  dev.new();

  plot(
      hinge_m_pred,
      xlim = c(xlim[1] - (diff(xlim) * 0.1), xlim[2] + (diff(xlim) * 0.1)),
      ylim = c(ylim[1], ylim[2] + 0.01),
      col = "red",
      main = ""
    )

  lines(
      hds_m_pred,
      col = "blue"
    )

  grid(); grid();
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
  # here's our built-in aic threshold method that needs to be refactored
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
  for(var in vars){n
    col <- as.vector(c( 
        which.max(abs(m_pca$rotation[var,])), 
        which.min(abs(m_pca$rotation[var,])) 
      ))
    x_hat <- m_pca$x[,col] %*% t(m_pca$rotation[,col])
      x_hat <- x_hat[,var] # retain only our variance for our focal variable
    # make sure the sign matches our original cov
    if ( cor(x_hat, df[,var]) < 0 ){
     x_hat <- -1 * x_hat
    }
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

glm_to_distsamp <- function(m=NULL, umdf=NULL){
  umdf@siteCovs <- m$data

  to_use <- names(coefficients(m))
    to_use <- to_use[2:length(to_use)]
  # drop any )1 or )2 poly() postfixes
  to_use <- unique(gsub(to_use, pattern="[)][0-9]", replacement=")"))

  to_use <- paste(
      "~1~",
      paste(to_use, collapse="+"),
      "+offset(log(effort))",
      sep=""
    )

  return(unmarked::distsamp(
      formula=as.formula(to_use),
      data=umdf,
      se=T,
      keyfun="halfnorm",
      unitsOut="kmsq",
      output="abund"
    ))
}

fit_distsamp <- function(lambdas=NULL, umdf=NULL){
  unmarked_models <- lapply(
      X=unmarked_models,
      FUN=function(m){
        unmarked::distsamp(
          formula=as.formula(paste("~1~",
            lambdas,
            sep=""
          )),
          data=umdf,
          se=T,
          keyfun="halfnorm",
          unitsOut="kmsq",
          output="abund"
        )
    })
}


aic_test_quadratic_terms_distsamp <- function(unmarked_model=NULL, original_formulas=NULL, umdf=NULL){
  quadratic_terms <- lapply(
      unmarked_models,
      FUN=quadratics_to_keep
    )
  # poisson version
  unmarked_models <- mapply(
    FUN=function(...){
      # drop the lam() prefix
      quadratics <- gsub(gsub(quadratics, pattern="lam[(]", replacement=""), pattern="[)][)]", replacement=")")
      # drop poly() notation from the list of all covariates using in this model
      vars <- gsub(gsub(vars, pattern="poly[(]", replacement=""), pattern="*.[0-2].*..*", replacement="")
      if(length(quadratics)>0){
        vars <- vars[!as.vector(sapply(vars, FUN=function(p){ sum(grepl(x=quadratics, pattern=p))>0  }))]
        # use AIC to justify our proposed quadratic terms
        for(q in quadratics){
          lin_var <- gsub(
            q,
            pattern=", 2,",
            replacement=", 1,"
          )

          m_lin_var <- OpenIMBCR:::AIC(unmarked::distsamp(
            formula=as.formula(paste("~1~",
                          paste(
                            c(paste("poly(", paste(vars, ", 1, raw=T)", sep=""), sep=""),lin_var,quadratics[!(quadratics %in% q)]), collapse="+"),
                          "+offset(log(effort))",
                          sep=""
            )),
            data=umdf,
            se=T,
            keyfun="halfnorm",
            unitsOut="kmsq",
            output="abund"
          ))+AIC_RESCALE_CONST
          m_quad_var <- OpenIMBCR:::AIC(unmarked::distsamp(
            formula=as.formula(paste("~1~",
                          paste(c(paste("poly(", paste(vars, ", 1, raw=T)", sep=""), sep=""), quadratics), collapse="+"),
                          "+offset(log(effort))",
                          sep=""
            )),
            data=umdf,
            se=T,
            keyfun="halfnorm",
            unitsOut="kmsq",
            output="abund"
          )) + AIC_RESCALE_CONST
          # if we don't improve our AIC with the quadratic by at-least 8 aic units
          # (pretty substatial support), keep the linear version
          if( m_lin_var-m_quad_var < AIC_SUBSTANTIAL_THRESHOLD){
            quadratics <- quadratics[!(quadratics %in% q)]
            vars <- c(
              vars,
              gsub(gsub(lin_var, pattern="poly[(]", replacement=""),
              pattern=", [0-9], raw.*=*.T[)]",
              replacement="")
            )
          }
        }

        vars <- c(paste("poly(", paste(vars, ", 1, raw=T)", sep=""), sep=""), quadratics)
      }
      return(vars)
    },
    quadratics=quadratic_terms,
    vars=original_formulas
  )
  return(unmarked_models)
}


#
# MAIN
#

# define the vars we are going to use
#vars <- c("grass_ar","shrub_ar","crp_ar","wetland_ar","pat_ct", "lat", "lon")
#vars <- c("grass_ar","shrub_ar","crp_ar","wetland_ar","pat_ct")
vars <- c("grass_ar","shrub_ar","crp_ar","wetland_ar","pat_ct", "mat", "map")

# calculate exhaustive (all possible) variable mCombinations
# for model selection

original_formulas <- unmarked_models <- calc_all_distsamp_combinations(vars)
  unmarked_models <- paste(unmarked_models, "+offset(log(effort))", sep="")

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

cat(
    " -- fitting a first round of negbin models and testing quadratic terms",
    "(this may take 40-80 minutes)\n")

# make an over-fit model of all variables and an underfit model to stare at
# and wonder
m_pois_full_model <- try(fit_gdistsamp(
  paste(paste("poly(",vars,",2,raw=T)",sep=""), collapse="+"),
  umdf=umdf,
  mixture="P"
))

m_pois_intercept_model <- try(fit_gdistsamp(
  "1",
  umdf=umdf,
  mixture="P"
))

m_negbin_full_model <- try(fit_gdistsamp(
  paste(paste("poly(",vars,",2,raw=T)",sep=""), collapse="+"),
  umdf=umdf,
  mixture="NB"
))

m_negbin_intercept_model <- try(fit_gdistsamp(
  "1",
  umdf=umdf,
  mixture="NB"
))

# favor the poisson mixture -- but if our data have high variance,
# it may just fail to converge. In that case, use the negative binomial
if(class(m_pois_full_model) != "try-error"){
  if( 
      ( (m_negbin_full_model@AIC + AIC_RESCALE_CONST) - (m_pois_full_model@AIC + AIC_RESCALE_CONST) ) 
      < (2*AIC_SUBSTANTIAL_THRESHOLD) 
  ){
    mixture <- "NB"
  } else {
    mixture <- "P"
  }
} else {
  stop("we couldn't fit a standard poisson model to this bird -- this shouldn't happen")
}

cat(" -- testing linear / quadratic terms using AIC:")

# takes 45-90 minutes -- note that we have offsets here baked-in to the unmarked
# model formulas
models <- fit_gdistsamp(unmarked_models, umdf, mixture=mixture)

# sometimes models fail to converge -- they are coded as NA's

original_formulas <- original_formulas[
    suppressWarnings( !sapply(models, FUN=function(x) sum(is.na(x)) >0 ) )
  ]

models <- suppressWarnings(
    models[ !sapply(models, FUN=function(x) sum(is.na(x)) >0 ) ]
  )

# takes about 45 minutes
unmarked_models <- aic_test_quadratic_terms_gdistsamp(
  unmarked_models=models,
  original_formulas=original_formulas,
  umdf=umdf,
  mixture=mixture
)

# refit our models using the specification justified from testing AIC
# across linear vs quadratic terms
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

# how does our dispersion look?
dispersion <- sapply(
  X=unmarked_models,
  FUN=function(m){
    # calculate k from number of parameters in our lambda
    # formula -- could include intercept parameters too
    k <- sum(
      grepl(
        unlist(strsplit(as.character(m@formula)[3], split="[+]")),
        pattern="poly")
    )
    # our detection bins across all IMBCR transects
    observed <- unmarked::getY(m@data)
    df <- (length(observed)-sum(is.na(observed)))-k
    return(round(AICcmodavg::Nmix.chisq(m)$chi.square / df, 2))
  }
)


# make a fitList
model_selection_table <- suppressWarnings(unmarked::modSel(unmarked::fitList(
  fits=unmarked_models
)))

MOD_SEL_THRESHOLD <- max(which(model_selection_table@Full$delta < AIC_SUBSTANTIAL_THRESHOLD))
# select the top models (by array position) that satisfy our threshold
MOD_SEL_THRESHOLD <- as.numeric(model_selection_table@Full$model)[1:MOD_SEL_THRESHOLD]
# read from units attribute table and drop anything that isn't numeric
predict_df <- units@data

# standard model averaging from unmarked -- sllloowww
# predicted <- unmarked::predict(
#   unmarked::fitList(unmarked_models[MOD_SEL_THRESHOLD]),
#   type="lambda",
#   newdata=predict_df
# )

aic_weights <- model_selection_table@Full$AICwt[
    MOD_SEL_THRESHOLD
]

predicted <- par_unmarked_predict(
  unmarked_models=unmarked_models[MOD_SEL_THRESHOLD],
  predict_df=predict_df,
  type="lambda",
  weights=aic_weights
)

# fancy model averaging from AICcmodavg (useful for getting averaged betas)
# predict_df <- AICcmodavg:::modavgPred.AICunmarkedFitDS(
#     cand.set = unmarked_models[MOD_SEL_THRESHOLD],
#     parm.type = "grass_ar",
#     newdata = predict_df
#   )

pred <- units
pred@data <- data.frame(pred = predicted);

rm(predicted)

cat(" -- writing to disk\n")

rgdal::writeOGR(
  pred,
  ".",
  tolower(paste(argv[2],
      "_imbcr_hds_negbin_prediction_",
      gsub(format(Sys.time(), "%b %d %Y"), pattern=" ", replacement="_"),
      sep="")
  ),
  driver="ESRI Shapefile",
  overwrite=T
)

r_data_file <- gsub(
    r_data_file,
    pattern="pois_glm",
    replacement="hds_negbin"
  )

save(
    compress=T,
    list=(ls(pattern="[a-z]")),
    file=r_data_file
  )
