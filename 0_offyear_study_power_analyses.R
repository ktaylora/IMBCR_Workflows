require(OpenIMBCR)
require(rgeos)
require(raster)
require(rgdal)
require(unmarked)

N_BS_REPLICATES = 200 # number of replicates to use for our bootstrapped operations

# load a previous model fitting workspace from Kyle that contains useful
# data for what we are doing here (units shapefile, m_scale, etc...)
load("/global_workspace/iplan_models/stable/weme_imbcr_hds_negbin_workflow_feb_15_2018.rdata")

#
# LOCAL FUNCTIONS FOR POWER ANALYSES AND RELATED TASKS
#

#' shorthand function will reverse a scale() operation on a scaled data.frame,
#' using a previously-fit m_scale scale() object
backscale_var <- function(var=NULL, df=NULL, m_scale=NULL) {
  return(
      df[, var] *
      attr(m_scale, 'scaled:scale')[var] +
      attr(m_scale, 'scaled:center')[var]
    )
}
#' shorthand function that will call gdistsamp with a user-specified
#' unmarked data.frame
fit_gdistsamp <- function(lambdas=NULL, umdf=NULL, mixture="NB") {
  if (length(lambdas) > 1) {
    cl <- parallel::makeCluster(parallel::detectCores() - 1)
    parallel::clusterExport(cl, varlist=c("umdf","mixture"))
    unmarked_models <- parallel::parLapply(
        cl=cl,
        X=lambdas,
        fun=function(m) {
          ret <- try(unmarked::gdistsamp(
            pformula = as.formula("~1"),
            lambdaformula = as.formula(paste("~", m, sep="")),
            phiformula = as.formula("~1"),
            data = umdf,
            se = T,
            K = max(rowSums(umdf@y)),
            keyfun = "halfnorm",
            unitsOut = "kmsq",
            mixture = mixture,
            output = "abund"
          ))
         if (class(ret) == "try-error") {
           return(NA)
         } else {
           return(ret)
         }
      })
    parallel::stopCluster(cl);
    return(unmarked_models);
  } else {
     ret <- try(unmarked::gdistsamp(
      pformula = as.formula("~1"),
      lambdaformula = as.formula(paste("~", unlist(lambdas), sep="")),
      phiformula = as.formula("~1"),
      data = umdf,
      se = T,
      K = max(rowSums(umdf@y)),
      keyfun = "halfnorm",
      unitsOut = "kmsq",
      mixture = mixture,
      output = "abund"
    ))
    if (class(ret) == "try-error") {
      return(NA)
    } else {
      return(ret)
    }
  }
}
#' build an unmarked data.frame (umdf) for fitting a distsamp model from
#' a matrix of 1-km2 site covariates
build_umdf_for_spp <- function(
  path="/global_workspace/imbcr_number_crunching/results/RawData_PLJV_IMBCR_20171017.csv",
  four_letter_code=NULL,
  site_covs=NULL) {
  colnames(site_covs) <- tolower(colnames(site_covs))
  # this returns an IMBCR SpatialPointsDataFrame
  imbcr_detections <-
    OpenIMBCR:::scrub_imbcr_df(
      OpenIMBCR::imbcrTableToShapefile(path),
      four_letter_code=four_letter_code,
      drop_na="none"  # keep all NA values
    )
  # only consider THESE transects, if specified by user
  if (!is.null(site_covs)) {
    transect_ids <- unique(tolower(site_covs$transectnum))
    matching_transects <- tolower(imbcr_detections$transectnum) %in%
      tolower(site_covs$transectnum)
    if ( sum(tolower(unique(imbcr_detections$transectnum)) %in%
             tolower(site_covs$transectnum) )  != nrow(site_covs) ) {
      warning("not all transect_ids were found in the IMBCR data.frame -- only found:",
              sum(matching_transects), "/", length(transect_ids), sep = "")
    }
    imbcr_detections <- imbcr_detections[matching_transects,]
    transect_ids <- transect_ids[ tolower(transect_ids) %in%
                                    tolower(imbcr_detections$transectnum) ]
  }

  # define an arbitrary 10 breaks that will be used
  # to construct distance bins
  breaks <- append(0,as.numeric(quantile(as.numeric(
    imbcr_detections$radialdistance),
    na.rm=T,
    probs=seq(0.05,0.90,length.out=9))
  ))

  detections <- OpenIMBCR:::calc_dist_bins(
    imbcr_detections,
    breaks=breaks
  )

  if (is.null(site_covs)) {
    return(detections)
  } else {
    return( unmarked::unmarkedFrameGDS(
      y=as.matrix(detections$y),
      siteCovs = site_covs[tolower(site_covs$transectnum) %in%
                           tolower(transect_ids) , ],
      dist.breaks = detections$breaks,
      numPrimary = 1,
      survey = "point",
      unitsIn = "m"
    )
    )
  }

}
#' downsample MEsquite and JUniper transects using their joint percent-cover
#' estimate -- this will balance our zero shrub transects with our non-zero
#' shrub transects
downsample_transects_by <- function(bird_data=NULL, by="mean", m_scale=NULL,
                                    inflation_threshold=1.1) {
  # back-scale if needed
  if (min(bird_data[ , grepl(colnames(bird_data), pattern = by) ]) < 0) {
    bird_data <- bird_data[ , grepl(colnames(bird_data), pattern = by) ]
    for (var in colnames(bird_data)[ which(grepl(colnames(bird_data),
                                                 pattern = by)) ]) {
      bird_data[ , var] <- backscale_var(var, df = bird_data, m_scale)
    }
  }

  zero_transects <- as.vector(which(rowSums(
    bird_data[ , grepl(colnames(bird_data), pattern = by)]) == 0))
  non_zero_transects <- as.vector(which(rowSums(
    bird_data[ , grepl(colnames(bird_data), pattern = by)]) != 0))

  if ( length(zero_transects) > inflation_threshold *
       length(non_zero_transects) ) {
    zero_transects <- sample(zero_transects, size = length(non_zero_transects))
    return(c(non_zero_transects, zero_transects))
  } else {
    return(1:nrow(bird_data))
  }
}
#' estimate standard residual error from standard deviation using
#' an estimator that can account for the degrees of freedom used in
#' our model
est_residual_sd <- function(m=NULL, log=F) {
  return(
      sqrt(
        sum( ( m$y - predict(m, type=ifelse(log, NULL, "response")) )^2 ) /
        ( m$df.residual )
      )
  )
}
#' a robust estimator of residual standard error that can account for
#' degrees of freedom in a model
est_residual_se <- function(m=NULL) {
  # A robust SE estimator that takes into account df
  # (still assumes error and sd are proportional and normal)
  return( est_residual_sd(m, log = F) / sqrt( m$df.residual ) ) # SD -> SE
}
#' estimate the number of parameters used in a model
est_k_parameters <- function(m=NULL) {
  k <- 0
  intercept_terms <- 0
  if ( inherits(m, "unmarkedFitGDS") ) {
    intercept_terms <- 3 # lambda, p, and dispersion
    k <- length(unlist(strsplit(paste(as.character(m@formula),
                                      collapse = ""), split = "[+]")))
    if (k == 1) { # no '+' signs separating terms in the model?
      k <- k - 1
    }
    k <- k + intercept_terms
  } else if ( inherits(m, "glm") ) {
    intercept_terms <- 1
    k <- length(unlist(strsplit(paste(as.character(m$formula),
                                      collapse = ""), split="[+]")))
    if (k == 1) { # no '+' signs separating terms in the model?
      k <- k - 1
    }
    k <- k + intercept_terms
  }
  return(k)
}
#'
#'
est_deviance_from_likelihood <- function(m=NULL) {
  if ( inherits(m, "unmarkedFitGDS") ) {
    return(2*as.numeric(abs(m@negLogLik)))
  } else if ( inherits(m, "glm") ) {
    return(-2*as.numeric(logLik(m)))
  } else {
    return(NA)
  }
}
#' a robust estimator for residual sum-of-square error that can account for
#' degrees of freedom in a model
est_residual_mse <- function(m, log=T) {
  if ( inherits(m, "unmarkedFitGDS") ) {
    df <- length(m@data@y) - est_k_parameters(m)
    if (log) {
      residuals <- log(residuals(m)^2)
      residuals[is.infinite(residuals)] <- 0
    } else {
      residuals <- residuals(m)^2
    }
    sum_sq_resid_err <- sum(colSums(residuals), na.rm = T)
    # Mean Sum of Squares Error
    sum_sq_resid_mse <- sum_sq_resid_err / df
    return(sum_sq_resid_mse)
  } else {
    return(NA)
  }
}
#' a robust estimator for residual sum-of-square error that can account for
#' degrees of freedom in a model
est_residual_sse <- function(m, log=F) {
  if ( inherits(m, "unmarkedFitGDS") ) {
    # do we want to log-transform our residuals, e.g. for Poisson data?
    if (log) {
      residuals <- log(residuals(m)^2)
        residuals[is.infinite(residuals)] <- 0
      sum_sq_resid_err <- sum(colSums(residuals^2), na.rm = T)
    } else {
      sum_sq_resid_err <- sum(colSums(residuals(m)^2), na.rm = T)
    }
    return(sum_sq_resid_err)
  } else {
    return(NA)
  }
}
#' estimate power from the effect (mean - 0) and residual error of a model
#' using Cohen's D statistic
#' @export
est_cohens_d_power <- function(m=NULL, report=T, alpha=0.05, log=T) {
  m_power <- NA
  z_alpha <- round( (1 - alpha)*2.06, 2 )
  if ( inherits(m, "glm") ) {
    # R's GLM interface reports standard deviation of residuals by default,
    # this derivation of Cohen's power accomodates SD
    # note: using the probability distribution function for our test
    # this is: probability of obtaining a value < Z_mean(pred-0) - 1.96
    m_power <- pnorm( (abs(mean(predict(m)) - 0) / sigma(m) ) *
                       sqrt(m$df.residual) - z_alpha , lower.tail = T)
    m_power <- ifelse( round( m_power, 4) == 0, 1 / m$df.residual,
                       round( m_power, 4) )
  } else if ( inherits(m, "unmarkedFitGDS") ) {
    # unmarked's HDS interface reports standard error of residuals by default,
    # by way of the Hessian. This derivation of Cohen's power accomodates SE
    df <- length(m@data@y) - est_k_parameters(m)
    m_power <- colMeans(unmarked::predict(m, type = "lambda")) # lambda, se, ...
    # log-scale our effect and se?
    if (log) m_power <- log(m_power)
    # note: using the probability distribution function for our test
    # this is: probability of obtaining a value < Z_mean(pred-0) - 1.96
    m_power <- pnorm( ((m_power[1] - 0) / m_power[2]) - z_alpha , lower.tail = T)
    m_power <- ifelse( round( m_power, 4) == 0, 1 / df, round( m_power, 4) )
  }
  if (report) {
    cat(" ######################################################\n")
    cat("  Cohen's Power Analysis\n")
    cat(" ######################################################\n")
    cat("  -- 1-beta (power) :", round(m_power, 4) ,"\n")
    cat("  -- significantly different than zero? :",
        as.character( m_power > (1 - alpha) ), "\n")
  }
  return(list(power = as.vector(m_power)))
}
#' estimate McFadden's pseudo r-squared
#' Note that this function will work with model objects fit with glm() in 'R',
#' but is primarily intended to work with model objects fit using the 'unmarked'
#' R package. There is some fairly-involved exception handling that has gone
#' into working with AIC and negative log-likelihood values reported by
#' unmarked that should give you pause. I try to be as verbose as I can
#' with warnings when I fudge numbers reported by unmarked models.
#' @export
est_pseudo_rsquared <- function(m=NULL, method="likelihood") {
  if ( inherits(m, "unmarkedFitGDS") ) {
    intercept_m <- unmarked::update(
      m,
      "~1+offset(log(effort))",
      "~1",
      "~1",
      mixture = m@mixture
    )
    if (grepl(tolower(method), pattern = "likelihood")) {
      # this is a k-parameter "adjusted" mcfadden's pseudo r-squared
      # sometimes "negLogLik" for the intercept model returned by unmarked
      # are negative and can't be easily compared with the full model. I
      # think this is a bug.
      if (intercept_m@negLogLike < 0 && m@negLogLike > 0) {
        warnings(paste(c("The signs of the negative log-likelihoods for the",
                 "intercept and alternative models provided do not agree;",
                 "this is odd. Assuming the absolute value of the intercept",
                 "model is the true negative log-likelihood. Please ",
                 "review the models you are comparing."), collapse = " ")
        )
      } else if (intercept_m@negLogLike < 0 && m@negLogLike < 0) {
        ADJ_LOGLIKE_CONST <- 2*abs(intercept_m@negLogLike)
        intercept_m@negLogLike <- intercept_m@negLogLike + ADJ_LOGLIKE_CONST
        m@negLogLike <- m@negLogLike + ADJ_LOGLIKE_CONST
      }
      # sanity-check : is the AIC of our null model lower than our alternative?
      # this suggests no support for our alternative model
      if (OpenIMBCR:::AICc(intercept_m) < OpenIMBCR:::AICc(m)) {
       warning(paste(c("The intercept model has a lower AIC than the proposed",
                       "alternative model. This shouldn't happen. Assuming the",
                       "intercept model's negative log-likelihood is an",
                       "absolute value and that unmarked made a mistake",
                       "in it's reporting. Please review the model objects",
                       "to double check that these changes make sense."),
               collapse = " ")
       )
       if (intercept_m@negLogLike < 0) {
         intercept_m@negLogLike <- abs(intercept_m@negLogLike)
       } else {
         warning(paste(c("There's no evidence that the negative",
                         "log-likeihood value of our alternative model is",
                         "wrong, and the intercept AIC is lower than",
                         "the alternative model's AIC -- assuming no support",
                         "for explaining residual error and returning 0"),
                       sep = " ")
         )
         return(0.00)
       }
      }
      m_k_adj_loglik <-
        abs(m@negLogLike) - est_k_parameters(m)
      intercept_m_k_adj_loglik <-
        abs(intercept_m@negLogLike) - est_k_parameters(intercept_m)
      # warn user if the loglikelihood of our full model
      # is lower for our null model
      if (intercept_m_k_adj_loglik < m_k_adj_loglik) {
        warning(c("the intercept likelihood is lower than the alternative model;",
                "this shouldn't happen and it suggests that there is no support",
                "for adding covariates to your model"))
      }
      # this r-squared estimator is a little more robust than:
      # 1 - (AIC(alt)/AIC(intercept)); but captures the same thing
      r_squared <- (intercept_m_k_adj_loglik - m_k_adj_loglik) / m_k_adj_loglik
    } else if (grepl(tolower(method), pattern = "mse")) {
      intercept_m <- unmarked::update(
        m,
        "~1+offset(log(effort))",
        "~1",
        "~1",
        mixture = m@mixture
      )
      r_squared <- ( est_residual_mse(intercept_m) - est_residual_mse(m)  ) / est_residual_mse(intercept_m)
    }
  } else if ( inherits(m, "glm") ) {
    # a standard glm internally calculates deviance and null deviance --
    # use it by default
    r_squared <- 1 - ( m$deviance / m$null.deviance )
  }
  return( ifelse(r_squared < 0, 0, round(r_squared, 2)) )
}
#' Cohen's (1988) power for null and alternative models that leverages
#' residual variance explained for effects sizes for a take on the f-ratio
#' test
est_cohens_f_power <- function(m_0=NULL, m_1=NULL, alpha=0.05)
{
  r_1 <- ifelse( is.numeric(m_1), m_1, est_pseudo_rsquared(m_1) )
  r_0 <- ifelse( is.numeric(m_0), m_0, est_pseudo_rsquared(m_0) )
  # estimate an effect size (f statistic)
  f_effect_size <-  (r_1 - r_0) / (1 - r_1)
  u <- length(m_0@data@y) - est_k_parameters(m_0) # degrees of freedom for null model
  v <- length(m_1@data@y) - est_k_parameters(m_1) # degrees of freedom for alternative model
  lambda <- f_effect_size * (u + v + 1)
  # calling the f-distribution probability density function (1 - beta)
  power <- pf(
    qf(alpha, u, v, lower = FALSE),
    u,
    v,
    lambda,
    lower = FALSE
  )
  return(list(power = power))
}
#' bs_est_cohens_f_power
bs_est_cohens_f_power <- function(
  m_0_formula=NULL,
  m_1_formula=NULL,
  bird_data=NULL,
  n=154,
  replace=T,
  m_scale=NULL,
  type="gdistsamp")
{
  # is this a standard glm?
  if (grepl(tolower(type), pattern = "glm")) {
    return(NA)
    # are we fitting a hierarchical model?
  } else if (grepl(tolower(type), pattern = "gdistsamp")) {
    # set-up our workspace for a parallelized operation
    cl <- parallel::makeCluster(parallel::detectCores() - 1)
    parallel::clusterExport(
      cl,
      varlist = c("downsample_transects_by",
                  "est_pseudo_rsquared", "est_k_parameters", "est_residual_mse",
                  "est_residual_sse","N_BS_REPLICATES", "backscale_var",
                  "est_cohens_f_power"),
      envir = globalenv()
    )
    parallel::clusterExport(
      cl,
      varlist = c("bird_data","n","replace","m_0_formula","m_1_formula","m_scale"),
      envir = environment()
    )
    # parallelize our cohen's d bootstrap operation
    cohens_f_n <- unlist(parallel::parLapply(
      cl = cl,
      X = 1:N_BS_REPLICATES,
      fun = function(i) {
        # balance our transects and resample (with replacement) how ever
        # many rows are needed to satisfy our parameter n
        rows <-  sample(
          downsample_transects_by(
            bird_data = bird_data@siteCovs,
            m_scale = m_scale
          ),
          size = n,
          replace = replace
        )
        bird_data@siteCovs <- bird_data@siteCovs[rows,]
        bird_data@y <- bird_data@y[rows,]
        bird_data@tlength <- rep(1, size = n)
        m_0 <- try(unmarked::gdistsamp(
          pformula = as.formula("~1"),
          lambdaformula = as.formula(paste("~", unlist(m_0_formula), sep = "")),
          phiformula = as.formula("~1"),
          data = bird_data,
          se = T,
          K = max(rowSums(bird_data@y)),
          keyfun = "halfnorm",
          unitsOut = "kmsq",
          mixture = "NB",
          output = "abund"
        ))
        m_1 <- try(unmarked::gdistsamp(
          pformula = as.formula("~1"),
          lambdaformula = as.formula(paste("~", unlist(m_1_formula), sep = "")),
          phiformula = as.formula("~1"),
          data = bird_data,
          se = T,
          K = max(rowSums(bird_data@y)),
          keyfun = "halfnorm",
          unitsOut = "kmsq",
          mixture = "NB",
          output = "abund"
        ))
        if (class(m_0) == "try-error" || class(m_1) == "try-error") {
          return(NA)
        } else {
          return(est_cohens_f_power(m_0 = m_0, m_1 = m_1)$power)
        }
      }
    ))
    parallel::stopCluster(cl);
    rm(cl);
    # check for normality
    if ( round(abs(median(cohens_f_n, na.rm=T) - mean(cohens_f_n, na.rm=T)), 2) != 0 ) {
      warning("cohen's d statistic looks skewed")
    }
    return(round(mean(cohens_f_n, na.rm = T), 2))
  }
}
#' a bootstrapped implementation of the Cohen's D test
bs_est_cohens_d_power <- function(formula=NULL, bird_data=NULL, n=154,
                                  replace=T, m_scale=NULL, type="gdistsamp") {
  # is this a standard glm?
  if (grepl(tolower(type), pattern = "glm")) {
    return(NA)
    # are we fitting a hierarchical model?
  } else if (grepl(tolower(type), pattern = "gdistsamp")) {
    # set-up our workspace for a parallelized operation
    cl <- parallel::makeCluster(parallel::detectCores() - 1)
    parallel::clusterExport(
      cl,
      varlist = c("downsample_transects_by",
                "est_pseudo_rsquared", "est_k_parameters", "est_residual_mse",
                "est_residual_sse","N_BS_REPLICATES", "backscale_var",
                "est_cohens_d_power"),
      envir = globalenv()
    )
    parallel::clusterExport(
      cl,
      varlist = c("bird_data","n","replace","formula","m_scale"),
      envir = environment()
    )
    # parallelize our cohen's d bootstrap operation
    cohens_d_n <- unlist(parallel::parLapply(
      cl = cl,
      X = 1:N_BS_REPLICATES,
      fun = function(i) {
        # balance our transects and resample (with replacement) how ever
        # many rows are needed to satisfy our parameter n
        rows <-  sample(
          downsample_transects_by(
            bird_data = bird_data@siteCovs,
            m_scale = m_scale
          ),
          size = n,
          replace = replace
        )
        bird_data@siteCovs <- bird_data@siteCovs[rows,]
        bird_data@y <- bird_data@y[rows,]
        bird_data@tlength <- rep(1, size = n)
        m <- try(unmarked::gdistsamp(
          pformula = as.formula("~1"),
          lambdaformula = as.formula(paste("~", unlist(formula), sep = "")),
          phiformula = as.formula("~1"),
          data = bird_data,
          se = T,
          K = max(rowSums(bird_data@y)),
          keyfun = "halfnorm",
          unitsOut = "kmsq",
          mixture = "NB",
          output = "abund"
        ))
        if (class(m) == "try-error") {
          return(NA)
        } else {
          return(est_cohens_d_power(m, report = F)$power)
        }
      }
    ))
    parallel::stopCluster(cl);
    rm(cl);
    # check for normality
    if ( round(abs(median(cohens_d_n) - mean(cohens_d_n)), 2) != 0 ) {
      warning("cohen's d statistic looks skewed")
    }
    return(round(mean(cohens_d_n, na.rm = T), 2))
  }
}
#' bs_est_pseudo_rsquared . This is a bootstrapped implementation of our
#' mcfadden's pseudo r squared estimator. It is very much in testing.
#' @export
bs_est_pseudo_rsquared <- function(
  formula=NULL,
  type="glm",
  bird_data=NULL,
  n=NULL,
  m_scale=NULL,
  replace=T,
  method="likelihood"
  ) {
    # is this a standard glm?
    if (grepl(tolower(type), pattern = "glm")) {
      pseudo_r_squared_n <- sapply(
        X=1:N_BS_REPLICATES,
        FUN=function(i) {
            m <- glm(
                formula = formula,
                data = bird_data[sample(downsample_transects_by(
                  bird_data = bird_data,
                  m_scale = m_scale),
                  size = n,
                  replace = replace),],
                family = poisson()
            );
            return(est_pseudo_rsquared(m, method=method))
        }
      )
    # are we fitting a hierarchical model?
    } else if (grepl(tolower(type), pattern = "gdistsamp")) {
      cl <- parallel::makeCluster(parallel::detectCores()-1)
      parallel::clusterExport(
        cl,
        varlist = c("downsample_transects_by",
                  "est_pseudo_rsquared", "est_k_parameters", "est_residual_mse",
                  "est_residual_sse","N_BS_REPLICATES", "backscale_var"),
        envir = globalenv()
      )
      parallel::clusterExport(
        cl,
        varlist = c("bird_data", "n", "replace","formula","m_scale", "method"),
        envir = environment()
      )
      # parallelize our r-squared bootstrap operation
      pseudo_r_squared_n <- unlist(parallel::parLapply(
        cl = cl,
        X = 1:N_BS_REPLICATES,
        fun = function(i) {
            # balance our transects and resample (with replacement) how ever
            # many rows are needed to satisfy our parameter n
            rows <-  sample(
              downsample_transects_by(
                bird_data = bird_data@siteCovs,
                m_scale = m_scale
              ),
              size = n,
              replace = replace
            )
            bird_data@siteCovs <- bird_data@siteCovs[rows,]
            bird_data@y <- bird_data@y[rows,]
            bird_data@tlength <- rep(1, size=n)
            m <- try(unmarked::gdistsamp(
                  pformula=as.formula("~1"),
                  lambdaformula=as.formula(paste("~", unlist(formula), sep="")),
                  phiformula=as.formula("~1"),
                  data=bird_data,
                  se=T,
                  K=max(rowSums(bird_data@y)),
                  keyfun="halfnorm",
                  unitsOut="kmsq",
                  mixture="NB",
                  output="abund"
            ))
            if (class(m) == "try-error") {
              return(NA)
            } else {
              return(est_pseudo_rsquared(m, method=method))
            }
        }
      ))
      parallel::stopCluster(cl);
      rm(cl);
      return(pseudo_r_squared_n)
    } else {
      return(NA)
    }
}
#' A boot strapped implementation of pearson's chi square "Goodness of Fit" test
#' stolen in part from the AICcmodavg package, which has a broken implementation
# bs_est_chisq_gof_test <- function (formula=NULL, nsim=10) {
#   require(parallel)
#
#   m <-
#     #gof <- Nmix.gof.test(mod=m, nsim=600, plot.hist=F)
#     model.type <- AICcmodavg:::Nmix.chisq(m)$model.type
#
#   cl <- parallel::makeCluster(parallel::detectCores()-1)
#
#   parallel::clusterExport(cl, varlist=c("nsim","umdf","m"))
#
#   bs <- parallel::parLapply(
#     cl=cl,
#     X=1:nsim,
#     fun=function(i) {
#       ret <- try(unmarked::parboot(
#         m,
#         statistic = function(i) {
#           AICcmodavg:::Nmix.chisq(i)$chi.square, nsim = 1, parallel = F)
#         }
#       )
#       if (class(ret) == "try-error") {
#         return(NA)
#       } else {
#         return(list(t0=as.vector(ret@t0), t.star=as.vector(ret@t.star)))
#       }
#     })
#
#   t.star <- unlist(bs) [ which(names(unlist(bs)) == "t.star") ]
#   t0 <- unlist(bs) [ which(names(unlist(bs)) == "t0") ]
#
#   nsim <- length(t.star) # update our nsim to account for any failures
#
#   p.value <- sum(t.star >= t0) / nsim
#
#   if (p.value == 0) {
#     p.display <- paste("<", round(1/length(t.star), 5))
#   } else {
#     p.display = paste("=", round(p.value, digits = 4))
#   }
#
#   if (plot.hist) {
#     hist(out@t.star, main = as.expression(substitute("Bootstrapped "*chi^2*" fit statistic ("*nsim*" samples)",
#                                                      list(nsim = nsim))),
#          xlim = range(c(out@t.star, out@t0)), xlab = paste("Simulated statistic ", "(observed = ", round(out@t0, digits = 2), ")", sep = ""))
#     title(main = bquote(paste(italic(P), .(p.display))), line = 0.5)
#     abline(v = out@t0, lty = "dashed", col = "red")
#   }
#   # estimate our dispersion parameter
#   c.hat.est <- t0/mean(t.star)
#
#   gof.out <- list(model.type = model.type, chi.square = mean(t0), t.star = mean(t.star), p.value = p.display, c.hat.est = mean(c.hat.est), nsim = nsim)
#
#   return(gof.out)
# }

#
# MAIN WORKFLOW
#

# meghan's shrub covariates with transect id's
shrubs <- read.csv("/home/ktaylora/Downloads/grids_w_shrubs.csv")
# only keep the covariates that Meghan calculated... we are going to
# do a spatial query and merge in with the other supplemental data, below
shrubs <- shrubs[ , grepl(colnames(shrubs), pattern="TransectNum|mean")]
# a shapefile of the original 2017 IMBCR data
transect_data_2017 <- OpenIMBCR::imbcrTableToShapefile(
  paste(sep="",
  "/global_workspace/imbcr_number_crunching/results/RawData_PLJV_IMBCR_2017101",
  "7.csv")
)
# back-scale the units of the data.frame we pulled out of our load() workflow
# so that it uses the original scale of the variables
units$grass_ar <- round(backscale_var("grass_ar", df=units@data, m_scale), 3)
units$shrub_ar <- round(backscale_var("shrub_ar", df=units@data, m_scale), 3)
units$crp_ar <- round(backscale_var("crp_ar", df=units@data, m_scale), 3)
units$wetland_ar <- round(backscale_var("wetland_ar", df=units@data, m_scale), 3)
units$pat_ct <- round(backscale_var("pat_ct", df=units@data, m_scale))
units$map <- round(backscale_var("map", df=units@data, m_scale), 3)
units$mat <- round(backscale_var("mat", df=units@data, m_scale), 3)
# calculate transect centroids in the order of meghan's shrub transects
transects <- lapply(
  unique(as.character(shrubs$TransectNum)),
  FUN=function(transect) {
    rgeos::gCentroid(transect_data_2017[transect_data_2017$transectnum == transect,])
})
# bind into a single SpatialPoints object
transects <- spTransform(
  do.call(rbind, transects),
  CRS(projection(units))
)
# make it a data.frame with meghan's shrub transect id's
transects <- SpatialPointsDataFrame(
  transects,
  data=data.frame(transectnum=unique(as.character(shrubs$TransectNum)))
)
# which grid units overlap with meghan's transects?
which_is_overlapping_grid_units <- which(!is.na(sp::over(units, transects)[,1]))
# bind meghan's original data.frame with the overlapping grid units data.frame
bird_data <- as.data.frame(cbind(shrubs, units[which_is_overlapping_grid_units,]@data))
# re-calculate our sampling effort per-transect for use as an offset in the models
bird_data$effort <- as.numeric(sapply(
  X=unique(as.character(shrubs$TransectNum)),
  FUN=function(transect) {
    return(
      max(
        transect_data_2017@data[transect_data_2017$transectnum == transect , 'point'],
        na.rm=T)
    )
}))
# estimate transect-level bird counts by species (GRSP and WEME)
bird_data$grsp_count <- as.numeric(sapply(
  X=unique(as.character(shrubs$TransectNum)),
  FUN=function(transect) {
    # filter out the dataset so that we are only looking at the focal transect
    this_transect <- transect_data_2017@data[transect_data_2017$transectnum == transect,]
    # do we see our species of interest in this transect?
    is_species_present <- sum(this_transect$birdcode == "GRSP", na.rm=T) > 0
    # if so, sum the radial distance observations that are greater than 0, else return 0
    if (is_species_present) {
      count <- sum(this_transect[this_transect$birdcode == "GRSP",'radialdistance'] >= 0, na.rm=T)
    } else {
      count <- 0
    }
    return(count)
}))
# ditto
bird_data$weme_count <- as.numeric(sapply(
  X=unique(as.character(shrubs$TransectNum)),
  FUN=function(transect) {
    this_transect <- transect_data_2017@data[transect_data_2017$transectnum == transect,]
    is_species_present <- sum(this_transect$birdcode == "WEME", na.rm=T) > 0
    if (is_species_present) {
      count <- sum(this_transect[this_transect$birdcode == "WEME",'radialdistance'] >= 0, na.rm=T)
    } else {
      count <- 0
    }
    return(count)
}))

# scale our explanatory dataset appropriately
colnames(bird_data) <- tolower(colnames(bird_data))
  to_scale <- !( colnames(bird_data) %in%
              c("grsp_count", "weme_count", "effort", "transectnum") )
m_scale <- scale(bird_data[,to_scale])
bird_data[,to_scale] <- scale(bird_data[,to_scale])
# fit some Poisson GLMs
m_pois_grsp <- glm(
  grsp_count ~ poly(crp_ar, 1, raw=T) + poly(grass_ar, 1, raw=T) +
  poly(shrub_ar, 1, raw=T) + poly(pat_ct, 1, raw=T) + poly(map, 2, raw=T) +
  poly(mat, 2, raw=T) + offset(log(effort)),
  data=bird_data,
  family=poisson()
)
m_pois_grsp_w_shrub_covs <- glm(
   grsp_count ~ poly(crp_ar, 1, raw=T) + poly(grass_ar, 1, raw=T) +
   poly(shrub_ar, 1, raw=T) + poly(pat_ct, 1, raw=T) + poly(map, 2, raw=T) +
   poly(mat, 2, raw=T) + poly(mean_me_sh, 1, raw=T) + poly(mean_me_o, 1, raw=T) +
   poly(mean_ju_sh, 1, raw=T) + poly(mean_ju_o, 1, raw=T) + offset(log(effort)),
   data=bird_data,
   family=poisson()
)
#
# build an unmarked compatible data.frame from our above data.frame
#
distance_models <- list()
distance_models$grsp <- list(umdf=build_umdf_for_spp(
  four_letter_code="GRSP",
  site_covs=bird_data
))
distance_models$weme <- list(umdf=build_umdf_for_spp(
  four_letter_code="WEME",
  site_covs=bird_data
))
distance_models$nobo <- list(umdf=build_umdf_for_spp(
  four_letter_code="NOBO",
  site_covs=bird_data
))
distance_models$hola <- list(umdf=build_umdf_for_spp(
  four_letter_code="HOLA",
  site_covs=bird_data
))
distance_models$casp <- list(umdf=build_umdf_for_spp(
  four_letter_code="CASP",
  site_covs=bird_data
))
distance_models$grro <- list(umdf=build_umdf_for_spp(
  four_letter_code="GRRO",
  site_covs=bird_data
))
distance_models$stfl <- list(umdf=build_umdf_for_spp(
  four_letter_code="STFL",
  site_covs=bird_data
))
distance_models$scqu <- list(umdf=build_umdf_for_spp(
  four_letter_code="SCQU",
  site_covs=bird_data
))
distance_models$losh <- list(umdf=build_umdf_for_spp(
  four_letter_code="LOSH",
  site_covs=bird_data
))
#
# fit some HDS model in unmarked
#

# GRSP
full_model_formula <- paste(
  c("poly(crp_ar, 1, raw = T) + poly(grass_ar, 1, raw = T) +",
    "poly(shrub_ar, 1, raw = T) + poly(pat_ct, 1, raw = T) +",
    "poly(mean_me_sh, 1, raw=T) + poly(mean_me_o, 1, raw=T) +",
    "poly(mean_ju_sh, 1, raw=T) + poly(mean_ju_o, 1, raw=T) +",
    "offset(log(effort))"),
  collapse = ""
)
minus_mesq_juni_model_formula <- paste(
  c("poly(crp_ar, 1, raw = T) + poly(grass_ar, 1, raw = T) +",
    "poly(shrub_ar, 1, raw = T) + poly(pat_ct, 1, raw = T) +",
    "offset(log(effort))"),
  collapse = ""
)
distance_models$grsp$minus_mesq_juni_model <- fit_gdistsamp(
  lambdas = minus_mesq_juni_model_formula,
  umdf = distance_models$grsp$umdf,
  mixture = "NB"
)
distance_models$grsp$full_model <- fit_gdistsamp(
  lambdas = full_model_formula,
  umdf = distance_models$grsp$umdf,
  mixture = "NB"
)
distance_models$grsp$intercept_model <- unmarked::update(
  distance_models$grsp$full_model,
  "~1+offset(log(effort))",
  "~1",
  "~1",
  mixture = "NB"
)
distance_models$grsp$pseudo_r_squared_n_154 <- bs_est_pseudo_rsquared(
  formula = full_model_formula,
  bird_data = distance_models$grsp$umdf,
  n = 154,
  replace = F,
  m_scale = m_scale,
  type = "gdistsamp"
)
distance_models$grsp$cohens_d_n_154 <- bs_est_cohens_d_power(
  formula = full_model_formula,
  bird_data = distance_models$grsp$umdf,
  n = 154,
  replace = T,
  m_scale = m_scale
)
distance_models$grsp$cohens_d_n_360 <- bs_est_cohens_d_power(
  formula = full_model_formula,
  bird_data = distance_models$grsp$umdf,
  n = 360,
  replace = T,
  m_scale = m_scale
)
distance_models$grsp$cohens_d_n_540 <- bs_est_cohens_d_power(
  formula = full_model_formula,
  bird_data = distance_models$grsp$umdf,
  n = 540,
  replace = T,
  m_scale = m_scale
)
distance_models$grsp$cohens_d_n_720 <- bs_est_cohens_d_power(
  formula = full_model_formula,
  bird_data = distance_models$grsp$umdf,
  n = 720,
  replace = T,
  m_scale = m_scale
)
distance_models$grsp$cohens_f_n_154 <- bs_est_cohens_f_power(
  minus_mesq_juni_model_formula,
  full_model_formula,
  n = 154,
  bird_data = distance_models$grsp$umdf,
  m_scale = m_scale
)
distance_models$grsp$cohens_f_n_540 <- bs_est_cohens_f_power(
  minus_mesq_juni_model_formula,
  full_model_formula,
  n = 540,
  bird_data = distance_models$grsp$umdf,
  m_scale = m_scale
)
distance_models$grsp$cohens_f_n_720 <- bs_est_cohens_f_power(
  minus_mesq_juni_model_formula,
  full_model_formula,
  n = 720,
  bird_data = distance_models$grsp$umdf,
  m_scale = m_scale
)
# WEME
full_model_formula <- paste(
  c("poly(grass_ar, 2, raw = T) + ",
    "poly(shrub_ar, 1, raw = T) + ",
    "poly(pat_ct, 1, raw = T) + ",
    "poly(map, 2, raw=T) + ",
    "poly(mean_me_sh, 1, raw=T) + ",
    "poly(mean_me_o, 1, raw=T) + ",
    "poly(mean_ju_sh, 1, raw=T) + ",
    "poly(mean_ju_o, 1, raw=T) + ",
    "offset(log(effort))"),
  collapse = ""
)
minus_mesq_juni_model_formula <- paste(
  c("poly(grass_ar, 2, raw = T) + ",
    "poly(shrub_ar, 1, raw = T) + ",
    "poly(pat_ct, 1, raw = T) + ",
    "poly(map, 2, raw=T) + ",
    "offset(log(effort))"),
  collapse = ""
)
distance_models$weme$full_model <- fit_gdistsamp(
  lambdas=full_model_formula,
  umdf=distance_models$weme$umdf,
  mixture="NB"
)
distance_models$weme$intercept_model <- unmarked::update(
  distance_models$weme$full_model,
  "~1+offset(log(effort))",
  "~1",
  "~1",
  mixture="NB"
)
distance_models$weme$pseudo_r_squared_n_154 <- bs_est_pseudo_rsquared(
  formula = full_model_formula,
  bird_data = distance_models$weme$umdf,
  n = 154,
  replace = F,
  m_scale = m_scale,
  type = "gdistsamp",
  method = "mse"
)
distance_models$weme$cohens_d_n_154 <- bs_est_cohens_d_power(
  formula = full_model_formula,
  bird_data = distance_models$weme$umdf,
  n = 154,
  replace = T,
  m_scale = m_scale
)
distance_models$weme$cohens_d_n_360 <- bs_est_cohens_d_power(
  formula = full_model_formula,
  bird_data = distance_models$weme$umdf,
  n = 360,
  replace = T,
  m_scale = m_scale
)
distance_models$weme$cohens_d_n_540 <- bs_est_cohens_d_power(
  formula = full_model_formula,
  bird_data = distance_models$weme$umdf,
  n = 540,
  replace = T,
  m_scale = m_scale
)
distance_models$weme$cohens_d_n_720 <- bs_est_cohens_d_power(
  formula = full_model_formula,
  bird_data = distance_models$weme$umdf,
  n = 720,
  replace = T,
  m_scale = m_scale
)
distance_models$weme$cohens_f_n_154 <- bs_est_cohens_f_power(
  minus_mesq_juni_model_formula,
  full_model_formula,
  n = 154,
  bird_data = distance_models$weme$umdf,
  m_scale = m_scale
)
distance_models$weme$cohens_f_n_540 <- bs_est_cohens_f_power(
  minus_mesq_juni_model_formula,
  full_model_formula,
  n = 540,
  bird_data = distance_models$weme$umdf,
  m_scale = m_scale
)
distance_models$weme$cohens_f_n_720 <- bs_est_cohens_f_power(
  minus_mesq_juni_model_formula,
  full_model_formula,
  n = 720,
  bird_data = distance_models$weme$umdf,
  m_scale = m_scale
)
# NOBO
full_model_formula <- paste(
  c("poly(crp_ar, 1, raw = T) + ",
    "poly(grass_ar, 2, raw = T) + ",
    "poly(shrub_ar, 1, raw = T) + ",
    "poly(mat, 2, raw=T) + ",
    "poly(map, 2, raw=T) + ",
    "poly(mean_me_sh, 1, raw=T) + ",
    "poly(mean_me_o, 1, raw=T) + ",
    "poly(mean_ju_sh, 1, raw=T) + ",
    "poly(mean_ju_o, 1, raw=T) + ",
    "offset(log(effort))"),
  collapse = ""
)
minus_mesq_juni_model_formula <- paste(
  c("poly(crp_ar, 1, raw = T) + ",
    "poly(grass_ar, 2, raw = T) + ",
    "poly(shrub_ar, 1, raw = T) + ",
    "poly(mat, 2, raw=T) + ",
    "poly(map, 2, raw=T) + ",
    "offset(log(effort))"),
  collapse = ""
)
distance_models$nobo$full_model <- fit_gdistsamp(
  lambdas=full_model_formula,
  umdf=distance_models$nobo$umdf,
  mixture="NB"
)
distance_models$nobo$intercept_model <- unmarked::update(
  distance_models$nobo$full_model,
  "~1",
  mixture="NB"
)
distance_models$nobo$pseudo_r_squared_n_154 <- bs_est_pseudo_rsquared(
  formula = full_model_formula,
  bird_data = distance_models$nobo$umdf,
  n = 154,
  replace = F,
  m_scale = m_scale,
  type = "gdistsamp"
)
distance_models$nobo$cohens_d_n_154 <- bs_est_cohens_d_power(
  formula = full_model_formula,
  bird_data = distance_models$nobo$umdf,
  n = 154,
  replace = T,
  m_scale = m_scale
)
distance_models$nobo$cohens_d_n_360 <- bs_est_cohens_d_power(
  formula = full_model_formula,
  bird_data = distance_models$nobo$umdf,
  n = 360,
  replace = T,
  m_scale = m_scale
)
distance_models$nobo$cohens_d_n_540 <- bs_est_cohens_d_power(
  formula = full_model_formula,
  bird_data = distance_models$nobo$umdf,
  n = 540,
  replace = T,
  m_scale = m_scale
)
distance_models$nobo$cohens_d_n_720 <- bs_est_cohens_d_power(
  formula = full_model_formula,
  bird_data = distance_models$nobo$umdf,
  n = 720,
  replace = T,
  m_scale = m_scale
)
distance_models$nobo$cohens_f_n_154 <- bs_est_cohens_f_power(
  minus_mesq_juni_model_formula,
  full_model_formula,
  n = 154,
  bird_data = distance_models$nobo$umdf,
  m_scale = m_scale
)
distance_models$nobo$cohens_f_n_360 <- bs_est_cohens_f_power(
  minus_mesq_juni_model_formula,
  full_model_formula,
  n = 360,
  bird_data = distance_models$nobo$umdf,
  m_scale = m_scale
)
distance_models$nobo$cohens_f_n_540 <- bs_est_cohens_f_power(
  minus_mesq_juni_model_formula,
  full_model_formula,
  n = 540,
  bird_data = distance_models$nobo$umdf,
  m_scale = m_scale
)
distance_models$nobo$cohens_f_n_720 <- bs_est_cohens_f_power(
  minus_mesq_juni_model_formula,
  full_model_formula,
  n = 720,
  bird_data = distance_models$nobo$umdf,
  m_scale = m_scale
)
# HOLA
full_model_formula <- paste(
  c("poly(crp_ar, 1, raw = T) + ",
    "poly(grass_ar, 1, raw = T) + ",
    "poly(shrub_ar, 1, raw = T) + ",
    "poly(mean_me_sh, 1, raw=T) + ",
    "poly(mean_me_o, 1, raw=T) + ",
    "poly(mean_ju_sh, 1, raw=T) + ",
    "poly(mean_ju_o, 1, raw=T) + ",
    "offset(log(effort))"),
  collapse = ""
)
minus_mesq_juni_model_formula <- paste(
  c("poly(crp_ar, 1, raw = T) + ",
    "poly(grass_ar, 1, raw = T) + ",
    "poly(shrub_ar, 1, raw = T) + ",
    "offset(log(effort))"),
  collapse = ""
)
distance_models$hola$full_model <- fit_gdistsamp(
  lambdas=full_model_formula,
  umdf=distance_models$hola$umdf,
  mixture="NB"
)
distance_models$hola$intercept_model <- unmarked::update(
  distance_models$hola$full_model,
  "~1+offset(log(effort))",
  "~1",
  "~1",
  mixture="NB"
)
distance_models$hola$pseudo_r_squared_n_154 <- bs_est_pseudo_rsquared(
  formula = full_model_formula,
  bird_data = distance_models$hola$umdf,
  n = 154,
  replace = F,
  m_scale = m_scale,
  type = "gdistsamp"
)
distance_models$hola$cohens_d_n_154 <- bs_est_cohens_d_power(
  formula = full_model_formula,
  bird_data = distance_models$hola$umdf,
  n = 154,
  replace = T,
  m_scale = m_scale
)
distance_models$hola$cohens_d_n_360 <- bs_est_cohens_d_power(
  formula = full_model_formula,
  bird_data = distance_models$hola$umdf,
  n = 360,
  replace = T,
  m_scale = m_scale
)
distance_models$hola$cohens_d_n_540 <- bs_est_cohens_d_power(
  formula = full_model_formula,
  bird_data = distance_models$hola$umdf,
  n = 540,
  replace = T,
  m_scale = m_scale
)
distance_models$hola$cohens_d_n_720 <- bs_est_cohens_d_power(
  formula = full_model_formula,
  bird_data = distance_models$hola$umdf,
  n = 720,
  replace = T,
  m_scale = m_scale
)
distance_models$hola$cohens_f_n_154 <- bs_est_cohens_f_power(
  minus_mesq_juni_model_formula,
  full_model_formula,
  n = 154,
  bird_data = distance_models$hola$umdf,
  m_scale = m_scale
)
distance_models$hola$cohens_f_n_360 <- bs_est_cohens_f_power(
  minus_mesq_juni_model_formula,
  full_model_formula,
  n = 360,
  bird_data = distance_models$hola$umdf,
  m_scale = m_scale
)
distance_models$hola$cohens_f_n_540 <- bs_est_cohens_f_power(
  minus_mesq_juni_model_formula,
  full_model_formula,
  n = 540,
  bird_data = distance_models$hola$umdf,
  m_scale = m_scale
)
distance_models$hola$cohens_f_n_720 <- bs_est_cohens_f_power(
  minus_mesq_juni_model_formula,
  full_model_formula,
  n = 720,
  bird_data = distance_models$hola$umdf,
  m_scale = m_scale
)
# CASP
distance_models$casp$full_model <- fit_gdistsamp(
  lambdas=as.character(full_model_formula)[[3]],
  umdf=distance_models$casp$umdf,
  mixture="NB"
)
distance_models$casp$intercept_model <- unmarked::update(
  distance_models$casp$full_model,
  "~1+offset(log(effort))",
  "~1",
  "~1",
  mixture="NB"
)
distance_models$casp$pseudo_r_squared_n_154 <- bs_est_pseudo_rsquared(
  formula = full_model_formula,
  bird_data = distance_models$casp$umdf,
  replace = F,
  n = 154,
  m_scale = m_scale,
  type = "gdistsamp"
)
distance_models$casp$cohens_d_n_154 <- bs_est_cohens_d_power(
  formula = full_model_formula,
  bird_data = distance_models$casp$umdf,
  n = 154,
  replace = T,
  m_scale = m_scale
)
distance_models$casp$cohens_d_n_360 <- bs_est_cohens_d_power(
  formula = full_model_formula,
  bird_data = distance_models$casp$umdf,
  n = 360,
  replace = T,
  m_scale = m_scale
)
distance_models$casp$cohens_d_n_540 <- bs_est_cohens_d_power(
  formula = full_model_formula,
  bird_data = distance_models$casp$umdf,
  n = 540,
  replace = T,
  m_scale = m_scale
)
distance_models$casp$cohens_d_n_720 <- bs_est_cohens_d_power(
  formula = full_model_formula,
  bird_data = distance_models$casp$umdf,
  n = 720,
  replace = T,
  m_scale = m_scale
)
# GRRO
full_model_formula <- paste(
  c("poly(crp_ar, 1, raw = T) + ",
    "poly(grass_ar, 1, raw = T) + ",
    "poly(shrub_ar, 1, raw = T) + ",
    "poly(mat, 1, raw = T) + ",
    "poly(map, 1, raw = T) + ",
    "poly(mean_me_sh, 1, raw=T) + ",
    "poly(mean_me_o, 1, raw=T) + ",
    "poly(mean_ju_sh, 1, raw=T) + ",
    "poly(mean_ju_o, 1, raw=T) + ",
    "offset(log(effort))"),
  collapse = ""
)
minus_mesq_juni_model_formula <- paste(
  c("poly(crp_ar, 1, raw = T) + ",
    "poly(grass_ar, 1, raw = T) + ",
    "poly(shrub_ar, 1, raw = T) + ",
    "offset(log(effort))"),
  collapse = ""
)
distance_models$grro$full_model <- fit_gdistsamp(
  lambdas=as.character(full_model_formula)[[3]],
  umdf=distance_models$grro$umdf,
  mixture="NB"
)
distance_models$grro$intercept_model <- unmarked::update(
  distance_models$grro$full_model,
  "~1+offset(log(effort))",
  "~1",
  "~1",
  mixture="NB"
)
distance_models$grro$pseudo_r_squared_n_154 <- bs_est_pseudo_rsquared(
  formula = full_model_formula,
  bird_data = distance_models$grro$umdf,
  n = 154,
  replace = F,
  m_scale = m_scale,
  type = "gdistsamp"
)
distance_models$grro$cohens_d_n_154 <- bs_est_cohens_d_power(
  formula = full_model_formula,
  bird_data = distance_models$grro$umdf,
  n = 154,
  replace = T,
  m_scale = m_scale
)
distance_models$grro$cohens_d_n_360 <- bs_est_cohens_d_power(
  formula = full_model_formula,
  bird_data = distance_models$grro$umdf,
  n = 360,
  replace = T,
  m_scale = m_scale
)
distance_models$grro$cohens_d_n_540 <- bs_est_cohens_d_power(
  formula = full_model_formula,
  bird_data = distance_models$grro$umdf,
  n = 540,
  replace = T,
  m_scale = m_scale
)
distance_models$grro$cohens_d_n_720 <- bs_est_cohens_d_power(
  formula = full_model_formula,
  bird_data = distance_models$grro$umdf,
  n = 720,
  replace = T,
  m_scale = m_scale
)
# STFL
distance_models$stfl$full_model <- fit_gdistsamp(
  lambdas=as.character(full_model_formula)[[3]],
  umdf=distance_models$stfl$umdf,
  mixture="NB"
)
distance_models$stfl$intercept_model <- unmarked::update(
  distance_models$stfl$full_model,
  "~1+offset(log(effort))",
  "~1",
  "~1",
  mixture="NB"
)odels$stfl$pseudo_r_squared_n_154 <- bs_est_pseudo_rsquared(
  formula = full_model_formula,
  bird_data = distance_models$stfl$umdf@siteCovs,
  n = 154,
  replace = F,
  m_scale = m_scale,
  type = "gdistsamp"
)
distance_models$stfl$cohens_d_n_154 <- bs_est_cohens_d_power(
  formula = full_model_formula,
  bird_data = distance_models$stfl$umdf,
  n = 154,
  replace = T,
  m_scale = m_scale
)
distance_models$stfl$cohens_d_n_360 <- bs_est_cohens_d_power(
  formula = full_model_formula,
  bird_data = distance_models$stfl$umdf,
  n = 360,
  replace = T,
  m_scale = m_scale
)
distance_models$stfl$cohens_d_n_540 <- bs_est_cohens_d_power(
  formula = full_model_formula,
  bird_data = distance_models$stfl$umdf,
  n = 540,
  replace = T,
  m_scale = m_scale
)
distance_models$stfl$cohens_d_n_720 <- bs_est_cohens_d_power(
  formula = full_model_formula,
  bird_data = distance_models$stfl$umdf,
  n = 720,
  replace = T,
  m_scale = m_scale
)
# SCQU
distance_models$scqu$full_model <- fit_gdistsamp(
  lambdas=as.character(full_model_formula)[[3]],
  umdf=distance_models$scqu$umdf,
  mixture="NB"
)
distance_models$scqu$intercept_model <- unmarked::update(
  distance_models$scqu$full_model,
  "~1+offset(log(effort))",
  "~1",
  "~1",
  mixture="NB"
)
distance_models$scqu$pseudo_r_squared_n_154 <- bs_est_pseudo_rsquared(
  formula=full_model_formula,
  bird_data=distance_models$scqu$umdf,
  n = 154,
  replace = F,
  m_scale=m_scale,
  type="gdistsamp"
)
distance_models$scqu$cohens_d_n_154 <- bs_est_cohens_d_power(
  formula = full_model_formula,
  bird_data = distance_models$scqu$umdf,
  n = 154,
  replace = T,
  m_scale = m_scale
)
distance_models$scqu$cohens_d_n_360 <- bs_est_cohens_d_power(
  formula = full_model_formula,
  bird_data = distance_models$scqu$umdf,
  n = 360,
  replace = T,
  m_scale = m_scale
)
distance_models$scqu$cohens_d_n_540 <- bs_est_cohens_d_power(
  formula = full_model_formula,
  bird_data = distance_models$scqu$umdf,
  n = 540,
  replace = T,
  m_scale = m_scale
)
distance_models$scqu$cohens_d_n_720 <- bs_est_cohens_d_power(
  formula = full_model_formula,
  bird_data = distance_models$scqu$umdf,
  n = 720,
  replace = T,
  m_scale = m_scale
)
# LOSH
distance_models$losh$full_model <- fit_gdistsamp(
  lambdas=as.character(full_model_formula)[[3]],
  umdf=distance_models$losh$umdf,
  mixture="NB"
)
distance_models$losh$intercept_model <- unmarked::update(
  distance_models$losh$full_model,
  "~1+offset(log(effort))",
  "~1",
  "~1",
  mixture="NB"
)
distance_models$losh$pseudo_r_squared_n_154 <- bs_est_pseudo_rsquared(
  formula=full_model_formula,
  bird_data=distance_models$losh$umdf,
  n=154,
  replace = F,
  m_scale=m_scale,
  type="gdistsamp"
)
distance_models$losh$cohens_d_n_154 <- bs_est_cohens_d_power(
  formula = full_model_formula,
  bird_data = distance_models$losh$umdf,
  n = 154,
  replace = T,
  m_scale = m_scale
)
distance_models$losh$cohens_d_n_360 <- bs_est_cohens_d_power(
  formula = full_model_formula,
  bird_data = distance_models$losh$umdf,
  n = 360,
  replace = T,
  m_scale = m_scale
)
distance_models$losh$cohens_d_n_540 <- bs_est_cohens_d_power(
  formula = full_model_formula,
  bird_data = distance_models$losh$umdf,
  n = 540,
  replace = T,
  m_scale = m_scale
)
distance_models$losh$cohens_d_n_720 <- bs_est_cohens_d_power(
  formula = full_model_formula,
  bird_data = distance_models$losh$umdf,
  n = 720,
  replace = T,
  m_scale = m_scale
)