require(OpenIMBCR)
require(rgeos)
require(raster)
require(rgdal)
require(unmarked)

N_BS_REPLICATES = 200 # number of replicates to use for our bootstrapped operations

# load a previous model fitting workspace from Kyle that contains useful
# data for what we are doing here (units shapefile, m_scale, etc...)


#
# LOCAL FUNCTIONS FOR POWER ANALYSES AND RELATED TASKS
#

#' hidden function that will download explanatory data from google drive
import_shrub_data <- function(
  url=paste(c("https://docs.google.com/spreadsheets/d/e/2PACX-1vTVV7wxAN1ras6Nmq_FM3aq",
  "S8b4bJc0SUD4HAvLw86ud2wyWWBkFZDiOKA0vEOg64J8rRx2APllafNf/pub?gid=1603612468&",
  "single=true&output=csv"), collapse = ""),
  temp_filename="shrub_data_import.csv"
  )
{
  unlink(temp_filename, force = T); download.file(url, temp_filename, quiet = T)
  t <- read.csv(temp_filename)
    unlink(temp_filename, force = T)
  return(t)
}
#' hidden function remove transects that have NA values from model fitting
drop_na_transects <- function(umdf=NULL){
  ret <- umdf
  rows_to_keep <- rowSums(is.na(ret@siteCovs)) == 0
  ret@siteCovs <- ret@siteCovs[ rows_to_keep, ]
  ret@y <- ret@y[ rows_to_keep , ]
  ret@tlength <- rep(1, sum(rows_to_keep))
  return(ret)
}
#' hidden function that will shuffle and up-or-down sample an unmarked
#' data.frame
shuffle <- function(umdf=NULL, n=NULL, replace=T){
  ret <- umdf
  rows_to_keep <- 1:nrow(ret@siteCovs)
  if (!is.null(n)) {
    rows_to_keep <- sample(rows_to_keep, size = n, replace = replace)
  } else {
    rows_to_keep <- sample(rows_to_keep, replace=F)
  }
  ret@siteCovs <- ret@siteCovs[ rows_to_keep, ]
  ret@y <- ret@y[ rows_to_keep , ]
  ret@tlength <- rep(1, length(rows_to_keep))
  return(ret)
}
#' hidden shorthand function will reverse a scale() operation on a scaled data.frame,
#' using a previously-fit m_scale scale() object
backscale_var <- function(var=NULL, df=NULL, m_scale=NULL) {
  return(
      df[, var] *
      attr(m_scale, 'scaled:scale')[var] +
      attr(m_scale, 'scaled:center')[var]
    )
}
#' hidden shorthand function that will call gdistsamp with a user-specified
#' unmarked data.frame
fit_gdistsamp <- function(lambdas=NULL, umdf=NULL, mixture="NB", ...) {
  df <- umdf
  if (length(lambdas) > 1) {
    cl <- parallel::makeCluster(parallel::detectCores() - 1)
    parallel::clusterExport(cl, varlist=c("df","mixture"))
    unmarked_models <- parallel::parLapply(
        cl = cl,
        X = lambdas,
        fun = function(m, ...) {
          ret <- try(unmarked::gdistsamp(
            lambdaformula = paste("~", m, collapse = ""),
            pformula = "~1",
            phiformula = "~1",
            data = df,
            se = T,
            K = max(rowSums(df@y)),
            keyfun = "halfnorm",
            unitsOut = "kmsq",
            mixture = mixture,
            output = "abund",
            method = "Nelder-Mead",
            ...
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
      lambdaformula = paste("~", unlist(lambdas), collapse = ""),
      pformula = "~1",
      phiformula = "~1",
      data = df,
      se = T,
      K = max(rowSums(df@y)),
      keyfun = "halfnorm",
      unitsOut = "kmsq",
      mixture = mixture,
      output = "abund",
      ...
    ))
    if (class(ret) == "try-error") {
      return(NA)
    } else {
      return(ret)
    }
  }
}
fit_gdistsamp_intercept <- function(m){
  df <- m@data
  intercept_m <- unmarked::update(
    m,
    "~1",
    "~1",
    "~1",
    mixture = m@mixture,
    data = df,
    se = F
  )
  return(intercept_m)
}
#'
sum_detections <- function(imbcr_df=NULL, four_letter_code=NULL){
  four_letter_code <- toupper(four_letter_code)
  colnames(imbcr_df) <- tolower(colnames(imbcr_df))
  sum_detections <- as.numeric(sapply(
    X = unique(as.character(imbcr_df$transectnum)),
    FUN = function(transect) {
      # filter out the dataset so that we are only looking at the focal transect
      this_transect <- imbcr_df@data[imbcr_df$transectnum == transect,]
      # do we see our species of interest in this transect?
      is_species_present <- sum(this_transect$birdcode == four_letter_code, na.rm=T) > 0
      # if so, sum the radial distance observations that are greater than 0, else return 0
      if (is_species_present) {
        count <- sum(this_transect[this_transect$birdcode == four_letter_code,'radialdistance'] >= 0, na.rm=T)
      } else {
        count <- 0
      }
      return(count)
    }))
  return(sum_detections)
}
sum_all_bird_detections <- function(imbcr_df=NULL, transect_ids=NULL, byid=F){
  if (inherits(imbcr_df, "Spatial")) {
    imbcr_df <- imbcr_df@data
  }
  colnames(imbcr_df) <- tolower(colnames(imbcr_df))
  # drop unknowns
  imbcr_df <- imbcr_df[ imbcr_df$birdcode != "UNKN" , ]
  # is this a robust detection with a distance observation?
  imbcr_df <- imbcr_df[ imbcr_df$radialdistance >= 0 , ]
  if (!is.null(transect_ids)) {
    imbcr_df <- imbcr_df[ imbcr_df$transectnum %in% transect_ids , ]
  }
  if (byid) {
    ret <- lapply(
      X = unique(as.character(transect_ids)),
      FUN = function(transect){
        return(table(imbcr_df[ imbcr_df$transectnum == transect, 'birdcode' ]))
      }
    )
  } else {
    return(table(imbcr_df$birdcode))
  }
  return(ret)
}
#' estimate alpha species richness and report relative abundance 
est_alpha_species_richness_matrix <- function(imbcr_df=NULL, transect_ids=NULL, byid=T){
  if (inherits(imbcr_df, "Spatial")) {
    imbcr_df <- imbcr_df@data
  }
  if(is.null(transect_ids)){
    transect_ids <- imbcr_df$transecnum
  }
  richness <- do.call(
    rbind,
    sum_all_bird_detections(imbcr_df, transect_ids = transect_ids, byid = byid)
  )
  # account for sampling effort
  richness <- round(t(t(richness) / imbcr_df$effort), 3)
  alpha_richness <- rowSums(richness) # estimate of alpha species richness
  
  sort <- order(rowSums(richness, na.rm=T), decreasing = T)
  richness <- cbind(data.frame(richness), transectnum=as.character(transect_ids))
  
  richness <- richness[ sort ,  ] # sort transect richness by region-wide richness
  alpha_richness <- alpha_richness[ sort ]
  
  return(cbind(alpha_richness=alpha_richness, richness))
}
#' hidden shorthand function that will build an unmarked data.frame (umdf)
#' for fitting a distsamp model from a matrix of 1-km2 site covariates
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
#' @export
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
#' @export
est_residual_sd <- function(m=NULL, log=F) {
  return(
      sqrt(
        sum( ( m$y - predict(m, type=ifelse(log, NULL, "response")) )^2 ) /
        ( m$df.residual )
      )
  )
}
#' a robust estimator of residual standard error that can account for
#' degrees of freedom in a generalized linear model. Unmarked calculates
#' standard error directly (using the hessian matrix), so this is only
#' really applicable to a generalized linear model
#' @export
est_residual_se <- function(m=NULL) {
  if ( inherits(m, "glm") ) {
    # A robust SE estimator that takes into account df
    # (still assumes error and sd are proportional and normal)
    return( est_residual_sd(m, log = F) / sqrt( m$df.residual ) ) # SD -> SE
  } else {
    return(NA)
  }
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

#' a robust estimator for residual sum-of-square error that can account for
#' degrees of freedom in a model
est_residual_mse <- function(m, log=T) {
  if ( inherits(m, "unmarkedFitGDS") ) {
    if ( sum(is.na(m@data@siteCovs)) > 0 ) {
      warning(paste(c("NA values found in site covariates table -- degrees of",
                "freedom estimate may be off depending on what was used",
                "for modeling"), collapse=" "))
    }
    df <- length(m@data@y) - est_k_parameters(m)
    # do we want to log-transform our residuals, e.g. to normalize Poisson data?
    if (log) {
      residuals <-  na.omit((abs(log(unmarked::residuals(m)^2))))
      residuals[is.infinite(residuals)] <- 0
    } else {
      residuals <- na.omit(abs(unmarked::residuals(m)^2))
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
est_residual_sse <- function(m, log=T) {
  if ( inherits(m, "unmarkedFitGDS") ) {
    # do we want to log-transform our residuals, e.g. to normalize Poisson data?
    if (log) {
      residuals <- abs(log(residuals(m)^2))
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
#' using Cohen's (1988) D test statistic
#' @export
est_cohens_d_power <- function(m=NULL, report=T, alpha=0.05, log=T) {
  m_power <- NA
  z_alpha <- round( (1 - alpha)*2.06, 2 )
  if ( inherits(m, "glm") ) {
    # R's GLM interface reports standard deviation of residuals by default,
    # this derivation of Cohen's power accomodates SD
    # note: using the probability distribution function for our test
    # this is: probability of obtaining a value < Z_mean(pred-0) - 1.96
    m_power <- pnorm( (abs(mean(predict(m), na.rm=T) - 0) / sigma(m, na.rm=T) ) *
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
# calculate a fitted poisson regression deviance statistic,
# e.g., : https://goo.gl/KdEUUa
est_deviance <- function(m, method="residuals"){
  if (grepl(tolower(method), pattern = "resid")) {
    if ( inherits(m, "unmarkedFitGDS") ) {
      observed <- unmarked::getY(m@data)
      expected <- unmarked::fitted(m)
      # Deviance of full model : 2*sum(obs*log(obs/predicted)-(obs-predicted))
      dev.part <- ( observed * log(observed/expected) ) - (observed - expected)
      sum.dev.part <- sum(dev.part,na.rm=T)
      dev.sum <- 2*sum.dev.part
      return(dev.sum)
    } else {
      return(NA)
    }
} else if (grepl(tolower(method), pattern = "likelihood")) {
    if ( inherits(m, "unmarkedFitGDS") ) {
      return(2*as.numeric(abs(m@negLogLik)))
    } else if ( inherits(m, "glm") ) {
      return(-2*as.numeric(logLik(m)))
    } else {
      return(NA)
    }
  }
}
#' estimate McFadden's pseudo r-squared
#' Note that this function will work with model objects fit with glm() in 'R',
#' but is primarily intended to work with model objects fit using the 'unmarked'
#' R package. There is some fairly-involved exception handling that has gone
#' into working with AIC and negative log-likelihood values reported by
#' unmarked that should give you pause. I try to be as verbose as I can
#' with warnings when I fudge numbers reported by unmarked models.
#' @export
est_pseudo_rsquared <- function(m=NULL, method="deviance") {
  if ( inherits(m, "unmarkedFitGDS") ) {
    df <- m@data
    intercept_m <- try(update(
      m,
      "~1",
      "~1",
      "~1",
      mixture = m@mixture,
      data = df,
      se=F
    ))
    if (class(intercept_m) == "try-error") {
      warning("failed to get an intercept-only model to converge")
      return(NA)
    }
    # A Deviance of Residuals estimate
    if (grepl(tolower(method), pattern = "dev")) {
      r_squared <- (est_deviance(intercept_m) - est_deviance(m)) / est_deviance(intercept_m)
    } else if (grepl(tolower(method), pattern = "mse")) {
      r_squared <- (est_residual_mse(intercept_m) - est_residual_mse(m)) / est_residual_mse(intercept_m)
    } else if (grepl(tolower(method), pattern = "likelihood")) {
      m_k_adj_loglik <- m@negLogLike - est_k_parameters(m)
      intercept_m_k_adj_loglik <- intercept_m@negLogLike - est_k_parameters(intercept_m)
      # warn user if the loglikelihood of our full model
      # is lower for our null model
      if (intercept_m_k_adj_loglik < m_k_adj_loglik) {
        warning(c("the intercept likelihood is lower than the alternative model;",
                "this shouldn't happen and it suggests that there is no support",
                "for adding covariates to your model"))
      }
      r_squared <- (intercept_m_k_adj_loglik - m_k_adj_loglik) / intercept_m_k_adj_loglik
    } else {
      stop("unknown method")
    }
  } else if ( inherits(m, "glm") ) {
    # a standard glm internally calculates deviance and null deviance --
    # use it by default
    r_squared <- (m$null.deviance - m$deviance) / m$null.deviance
  }
  return( round(r_squared, 4) )
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
#' a bootstrapped implementation of the Cohens (1988) f-squared test
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
      varlist = c("downsample_transects_by", "est_deviance", "shuffle",
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
        bird_data <- shuffle(bird_data, n = n, replace = replace)
        m_0 <- try(unmarked::gdistsamp(
          pformula = as.formula("~1"),
          lambdaformula = as.formula(paste("~", unlist(m_0_formula), sep = "")),
          phiformula = as.formula("~1"),
          data = bird_data,
          se = F,
          K = max(rowSums(bird_data@y)),
          keyfun = "halfnorm",
          unitsOut = "kmsq",
          mixture = "NB",
          output = "abund",
          method = "Nelder-Mead"
        ))
        m_1 <- try(unmarked::gdistsamp(
          pformula = "~1",
          lambdaformula = paste("~", unlist(m_1_formula), collapse = ""),
          phiformula = "~1",
          data = bird_data,
          se = F,
          K = max(rowSums(bird_data@y)),
          keyfun = "halfnorm",
          unitsOut = "kmsq",
          mixture = "NB",
          output = "abund",
          method = "Nelder-Mead"
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
#' a bootstrapped implementation of the Cohen's (1988) D test
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
      varlist = c("downsample_transects_by", "shuffle",
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
        bird_data <- shuffle(bird_data, n = n, replace = replace)
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
          output = "abund",
          method = "Nelder-Mead"
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
    if ( round(abs(median(cohens_d_n, na.rm=T) - mean(cohens_d_n, na.rm=T)), 2) != 0 ) {
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
  type="gdistsamp",
  bird_data=NULL,
  n=NULL,
  m_scale=NULL,
  replace=T,
  method="deviance"
  ) {
    # is this a standard glm?
    if (grepl(tolower(type), pattern = "glm")) {
      bird_data <- shuffle(bird_data, n = n, replace = replace)
      pseudo_r_squared_n <- sapply(
        X=1:N_BS_REPLICATES,
        FUN=function(i) {
            m <- glm(
                formula = formula,
                data = bird_data,
                family = poisson()
            );
            return(est_pseudo_rsquared(m, method=method))
        }
      )
    # are we fitting a hierarchical model?
    } else if (grepl(tolower(type), pattern = "gdistsamp")) {
      cl <- parallel::makeCluster(parallel::detectCores() - 1)
      parallel::clusterExport(
        cl,
        varlist = c("downsample_transects_by", "est_deviance", "shuffle",
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
            bird_data <- shuffle(bird_data, n = n, replace = replace)
            m <- try(unmarked::gdistsamp(
                  pformula = "~1",
                  lambdaformula = paste("~", unlist(formula), collapse=""),
                  phiformula = "~1",
                  data = bird_data,
                  se = T,
                  K = max(rowSums(bird_data@y)),
                  keyfun = "halfnorm",
                  unitsOut = "kmsq",
                  mixture = "NB",
                  output = "abund",
                  method = "Nelder-Mead"
            ))
            if (class(m) == "try-error") {
              return(NA)
            } else {
              return(est_pseudo_rsquared(m, method = method))
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
#' bootstrap
#'
bootstrap <- function(
  formula=NULL,
  df=NULL,
  summary_function=NULL,
  n=NULL,
  replace=T,
  bs_replicates=N_BS_REPLICATES,
  fit_function=fit_gdistsamp,
  shuffle_function=shuffle)
{
  require(parallel)
  cl <- parallel::makeCluster(parallel::detectCores() - 1, outfile = '')
  parallel::clusterExport(
    cl,
    varlist = c("formula", "df", "n", "replace",
                "summary_function", "bs_replicates", "fit_function",
                "shuffle_function"),
    envir = environment()
  )
  test <- unlist(parallel::parLapply(
    cl = cl,
    X = 1:bs_replicates,
    fun = function(x){
      print("doing shuffle operation")
      df <- shuffle_function(df, replace = replace, n = n)
      print("doing fit operation")
      m <-  fit_function(formula, df)
      print("doing summary operation")
      ret <- summary_function(m, df);
      print("returning");
      return(ret);
    }
  ))
  parallel::stopCluster(cl);
  rm(cl);
  return(test)
}

#
# MAIN WORKFLOW
#

# the units dataset
units <- OpenIMBCR:::readOGRfromPath("/gis_data/Grids/units_attributed_v.3.0.shp")
# meghan's shrub covariates with transect id's
shrubs <- import_shrub_data()
# only keep the covariates that Meghan calculated... we are going to
# do a spatial query and merge in with the other supplemental data, below
shrubs <- shrubs[ , grepl(colnames(shrubs), pattern="transectnum|mean|morap")]
# a shapefile of the original 2017 IMBCR data
transect_data_2017 <- OpenIMBCR::imbcrTableToShapefile(
  paste(sep="",
  "/global_workspace/imbcr_number_crunching/results/RawData_PLJV_IMBCR_2017101",
  "7.csv")
)
# calculate transect centroids in the order of meghan's shrub transects
transects <- do.call(
  rbind,
  lapply(
    X=unique(as.character(shrubs$transectnum)),
    FUN=function(transect) {
      pt <- rgeos::gCentroid(transect_data_2017[transect_data_2017$transectnum == transect,])
      pt <- SpatialPointsDataFrame(pt, data=data.frame(transectnum=transect))
      return(pt)
  })
)
# bind into a single SpatialPoints object
transects <- OpenIMBCR:::spatial_join(transects, units)
transects@data <- merge(
  transects@data,
  shrubs,
  by="transectnum"
)

transects@data$effort <- as.numeric(sapply(
  X = unique(as.character(transects@data$transectnum)),
  FUN = function(transect) {
    return(
      length(unique(
        transect_data_2017@data[transect_data_2017$transectnum == transect , 'point'],
        na.rm = T))
    )
}))

# estimate transect-level bird counts by species (GRSP and WEME)
transects$grsp_count <- as.numeric(sapply(
  X=unique(as.character(shrubs$transectnum)),
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
transects$weme_count <- as.numeric(sapply(
  X=unique(as.character(shrubs$transectnum)),
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

# write to disk for validation
writeOGR(
  transects,
  dsn="/home/ktaylora/offyear_study_transects.json",
  layer="offyear_study_transects",
  driver="GeoJSON",
  overwrite=T
)

bird_data <- transects@data
bird_data_pre_scaling <- bird_data
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
m_pois_weme <- glm(
  weme_count ~ poly(crp_ar, 1, raw=T) + poly(grass_ar, 1, raw=T) +
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
m_pois_weme_w_shrub_covs <- glm(
  weme_count ~ poly(crp_ar, 1, raw=T) + poly(grass_ar, 1, raw=T) +
    poly(shrub_ar, 1, raw=T) + poly(pat_ct, 1, raw=T) + poly(map, 2, raw=T) +
    poly(mat, 2, raw=T) + poly(mean_me_sh, 1, raw=T) + poly(mean_me_o, 1, raw=T) +
    poly(mean_ju_sh, 1, raw=T) + poly(mean_ju_o, 1, raw=T) + offset(log(effort)),
  data=bird_data,
  family=poisson()
)
#
# Species Richness
#
richness <- do.call(
  rbind,
  sum_all_bird_detections(transect_data_2017, transect_ids = shrubs$transectnum, byid = T)
)
# account for sampling effort
richness <- round(t(t(richness) / bird_data$effort), 3)
alpha_richness <- rowSums(richness) # estimate of alpha species richness

sort <- order(rowSums(richness, na.rm=T), decreasing = T)
richness <- cbind(data.frame(richness), transectnum=as.character(shrubs$transectnum))

richness <- richness[ sort ,  ] # sort transect richness by region-wide richness
alpha_richness <- alpha_richness[ sort ]

richness_cor_matrix <- cbind(alpha_richness=alpha_richness, merge(richness, bird_data, by="transectnum"))

richness_cor_matrix <- richness_cor_matrix[ , !grepl(colnames(richness_cor_matrix), pattern="transect|count|effort") ]
richness_cor_matrix <- as.matrix(richness_cor_matrix)
# pretty up the labels
colnames(richness_cor_matrix) <- gsub(colnames(richness_cor_matrix), pattern="map", replacement="precipitation")
colnames(richness_cor_matrix) <- gsub(colnames(richness_cor_matrix), pattern="mat", replacement="temperature")
colnames(richness_cor_matrix) <- gsub(colnames(richness_cor_matrix), pattern="mean_me", replacement="mesquite")
colnames(richness_cor_matrix) <- gsub(colnames(richness_cor_matrix), pattern="mean_ju", replacement="ERC")
colnames(richness_cor_matrix) <- gsub(colnames(richness_cor_matrix), pattern="me_sh", replacement="mesquite")
colnames(richness_cor_matrix) <- gsub(colnames(richness_cor_matrix), pattern="ju_sh", replacement="ERC")
colnames(richness_cor_matrix) <- gsub(colnames(richness_cor_matrix), pattern="_o", replacement="_overstory")
colnames(richness_cor_matrix) <- gsub(colnames(richness_cor_matrix), pattern="_sh", replacement="_shrub_layer")
colnames(richness_cor_matrix) <- gsub(colnames(richness_cor_matrix), pattern="_ar", replacement="_area_NASS")
colnames(richness_cor_matrix) <- gsub(colnames(richness_cor_matrix), pattern="_ct", replacement="chiness")
colnames(richness_cor_matrix) <- gsub(colnames(richness_cor_matrix), pattern="morap", replacement="MoRAP")
colnames(richness_cor_matrix) <- gsub(colnames(richness_cor_matrix), pattern="_", replacement=" ")
colnames(richness_cor_matrix) <- gsub("(^|[[:space:]])([[:alpha:]])", "\\1\\U\\2", colnames(richness_cor_matrix), perl=TRUE)

corrplot::corrplot(cor(na.omit(richness_cor_matrix)),
                   method="color")

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
full_model_formula_imbcr_covs <- paste(
  c("poly(crp_ar, 1, raw = T) +",
    "poly(grass_ar, 1, raw = T) +",
    "poly(shrub_ar, 1, raw = T) +",
    "poly(pat_ct, 1, raw = T) +",
    "poly(map, 2, raw = T) +",
    "poly(mat, 2, raw=T) +",
    "poly(mean_me_sh, 1, raw=T) +",
    "poly(mean_me_o, 1, raw=T) +",
    "poly(mean_ju_sh, 1, raw=T) +",
    "poly(mean_ju_o, 1, raw=T) +",
    "offset(log(effort))"),
  collapse = ""
)
full_model_formula_morap_covs <- paste(
  c("poly(crp_ar, 1, raw = T) +",
    "poly(grass_ar, 1, raw = T) +",
    "poly(shrub_ar, 1, raw = T) +",
    "poly(pat_ct, 1, raw = T) +",
    "poly(map, 2, raw = T) +",
    "poly(mat, 2, raw=T) +",
    "poly(morap_me_sh, 1, raw=T) +",
    "poly(morap_ju_sh, 1, raw=T) +",
    "offset(log(effort))"),
  collapse = ""
)
null_formula <- paste(
  c("poly(crp_ar, 1, raw = T) +",
    "poly(grass_ar, 1, raw = T) +",
    "poly(shrub_ar, 1, raw = T) +",
    "poly(pat_ct, 1, raw = T) +",
    "poly(map, 2, raw = T) +",
    "poly(mat, 2, raw=T) +",
    "offset(log(effort))"),
  collapse = ""
)
distance_models$grsp$null_model <- fit_gdistsamp(
  lambdas = null_formula,
  umdf = distance_models$grsp$umdf,
  mixture = "NB"
)
distance_models$grsp$full_model_imbcr_covs <- fit_gdistsamp(
  lambdas = full_model_formula_imbcr_covs,
  umdf = distance_models$grsp$umdf,
  mixture = "NB"
)
distance_models$grsp$full_model_morap_covs <- fit_gdistsamp(
  lambdas = full_model_formula_morap_covs,
  umdf = drop_na_transects(distance_models$grsp$umdf),
  mixture = "NB"
)
distance_models$grsp$intercept_model_morap_covs <- fit_gdistsamp_intercept(
  m = distance_models$grsp$full_model_morap_covs
)
distance_models$grsp$intercept_model_imbcr_covs <- fit_gdistsamp_intercept(
  m = distance_models$grsp$full_model_imbcr_covs
)
distance_models$grsp$pseudo_r_squared_n_154 <- bs_est_pseudo_rsquared(
  formula = full_model_formula_imbcr_covs,
  bird_data = distance_models$grsp$umdf,
  n = 154,
  replace = T,
  m_scale = m_scale,
  type = "gdistsamp",
  method = "deviance"
)
distance_models$grsp$pseudo_r_squared_n_84 <- bs_est_pseudo_rsquared(
  formula = full_model_formula_morap_covs,
  bird_data = drop_na_transects(distance_models$grsp$umdf),
  n = 84,
  replace = T,
  m_scale = m_scale,
  type = "gdistsamp",
  method = "deviance"
)
distance_models$grsp$cohens_d_n_154 <- bs_est_cohens_d_power(
  formula = full_model_formula_imbcr_covs,
  bird_data = distance_models$grsp$umdf,
  n = 154,
  replace = T,
  m_scale = m_scale
)
distance_models$grsp$cohens_d_n_360 <- bs_est_cohens_d_power(
  formula = full_model_formula_imbcr_covs,
  bird_data = distance_models$grsp$umdf,
  n = 360,
  replace = T,
  m_scale = m_scale
)
distance_models$grsp$cohens_d_n_540 <- bs_est_cohens_d_power(
  formula = full_model_formula_imbcr_covs,
  bird_data = distance_models$grsp$umdf,
  n = 540,
  replace = T,
  m_scale = m_scale
)
distance_models$grsp$cohens_d_n_720 <- bs_est_cohens_d_power(
  formula = full_model_formula_imbcr_covs,
  bird_data = distance_models$grsp$umdf,
  n = 720,
  replace = T,
  m_scale = m_scale
)
distance_models$grsp$cohens_f_n_154 <- bs_est_cohens_f_power(
  null_formula,
  full_model_formula_imbcr_covs,
  n = 154,
  bird_data = distance_models$grsp$umdf,
  m_scale = m_scale
)
distance_models$grsp$cohens_f_n_360 <- bs_est_cohens_f_power(
  null_formula,
  full_model_formula_imbcr_covs,
  n = 360,
  bird_data = distance_models$grsp$umdf,
  m_scale = m_scale
)
distance_models$grsp$cohens_f_n_540 <- bs_est_cohens_f_power(
  null_formula,
  full_model_formula_imbcr_covs,
  n = 540,
  bird_data = distance_models$grsp$umdf,
  m_scale = m_scale
)
distance_models$grsp$cohens_f_n_720 <- bs_est_cohens_f_power(
  null_formula,
  full_model_formula_imbcr_covs,
  n = 720,
  bird_data = distance_models$grsp$umdf,
  m_scale = m_scale
)
# WEME
full_model_formula_imbcr_covs <- paste(
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

full_model_formula_morap_covs <- paste(
  c("poly(grass_ar, 2, raw = T) + ",
    "poly(shrub_ar, 1, raw = T) + ",
    "poly(pat_ct, 1, raw = T) + ",
    "poly(map, 2, raw=T) + ",
    "poly(morap_me_sh, 1, raw=T) + ",
    "poly(morap_ju_sh, 1, raw=T) + ",
    "offset(log(effort))"),
  collapse = ""
)
null_formula <- paste(
  c("poly(grass_ar, 2, raw = T) + ",
    "poly(shrub_ar, 1, raw = T) + ",
    "poly(pat_ct, 1, raw = T) + ",
    "poly(map, 2, raw=T) + ",
    "offset(log(effort))"),
  collapse = ""
)
distance_models$weme$null_model <- fit_gdistsamp(
  lambdas = null_formula,
  umdf = distance_models$weme$umdf,
  mixture = "NB"
)
distance_models$weme$full_model_imbcr_covs <- fit_gdistsamp(
  lambdas = full_model_formula_imbcr_covs,
  umdf = distance_models$weme$umdf,
  mixture = "NB"
)
distance_models$weme$full_model_morap_covs <- fit_gdistsamp(
  lambdas = full_model_formula_morap_covs,
  umdf = drop_na_transects(distance_models$weme$umdf),
  mixture = "NB"
)
distance_models$weme$intercept_model_imbcr_covs <- fit_gdistsamp_intercept(
  distance_models$weme$full_model_imbcr_covs
)
distance_models$weme$intercept_model_morap_covs <- fit_gdistsamp_intercept(
  distance_models$weme$full_model_morap_covs
)
distance_models$weme$pseudo_r_squared_n_154 <- bs_est_pseudo_rsquared(
  formula = full_model_formula_imbcr_covs,
  bird_data = distance_models$weme$umdf,
  n = 154,
  replace = T,
  m_scale = m_scale,
  type = "gdistsamp",
  method = "deviance"
)
distance_models$weme$pseudo_r_squared_n_84 <- bs_est_pseudo_rsquared(
  formula = full_model_formula_morap_covs,
  bird_data = drop_na_transects(distance_models$weme$umdf),
  n = 84,
  replace = T,
  m_scale = m_scale,
  type = "gdistsamp",
  method = "deviance"
)
distance_models$weme$cohens_d_n_154 <- bs_est_cohens_d_power(
  formula = full_model_formula_imbcr_covs,
  bird_data = distance_models$weme$umdf,
  n = 154,
  replace = T,
  m_scale = m_scale
)
distance_models$weme$cohens_d_n_360 <- bs_est_cohens_d_power(
  formula = full_model_formula_imbcr_covs,
  bird_data = distance_models$weme$umdf,
  n = 360,
  replace = T,
  m_scale = m_scale
)
distance_models$weme$cohens_d_n_540 <- bs_est_cohens_d_power(
  formula = full_model_formula_imbcr_covs,
  bird_data = distance_models$weme$umdf,
  n = 540,
  replace = T,
  m_scale = m_scale
)
distance_models$weme$cohens_d_n_720 <- bs_est_cohens_d_power(
  formula = full_model_formula_imbcr_covs,
  bird_data = distance_models$weme$umdf,
  n = 720,
  replace = T,
  m_scale = m_scale
)
distance_models$weme$cohens_f_n_154 <- bs_est_cohens_f_power(
  null_formula,
  full_model_formula_imbcr_covs,
  n = 154,
  bird_data = distance_models$weme$umdf,
  m_scale = m_scale
)
distance_models$weme$cohens_f_n_360 <- bs_est_cohens_f_power(
  null_formula,
  full_model_formula_imbcr_covs,
  n = 360,
  bird_data = distance_models$weme$umdf,
  m_scale = m_scale
)
distance_models$weme$cohens_f_n_540 <- bs_est_cohens_f_power(
  null_formula,
  full_model_formula_imbcr_covs,
  n = 540,
  bird_data = distance_models$weme$umdf,
  m_scale = m_scale
)
distance_models$weme$cohens_f_n_720 <- bs_est_cohens_f_power(
  null_formula,
  full_model_formula_imbcr_covs,
  n = 720,
  bird_data = distance_models$weme$umdf,
  m_scale = m_scale
)
# NOBO
full_model_formula_imbcr_covs <- paste(
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

full_model_formula_morap_covs <- paste(
  c("poly(crp_ar, 1, raw = T) + ",
    "poly(grass_ar, 2, raw = T) + ",
    "poly(shrub_ar, 1, raw = T) + ",
    "poly(mat, 2, raw=T) + ",
    "poly(map, 2, raw=T) + ",
    "poly(morap_me_sh, 1, raw=T) + ",
    "poly(morap_ju_sh, 1, raw=T) + ",
    "offset(log(effort))"),
  collapse = ""
)

null_formula <- paste(
  c("poly(crp_ar, 1, raw = T) + ",
    "poly(grass_ar, 2, raw = T) + ",
    "poly(shrub_ar, 1, raw = T) + ",
    "poly(mat, 2, raw=T) + ",
    "poly(map, 2, raw=T) + ",
    "offset(log(effort))"),
  collapse = ""
)

distance_models$nobo$null_model <- fit_gdistsamp(
  lambdas = null_formula,
  umdf = distance_models$nobo$umdf,
  mixture = "NB"
)
distance_models$nobo$full_model_imbcr_covs <- fit_gdistsamp(
  lambdas = full_model_formula_imbcr_covs,
  umdf = distance_models$nobo$umdf,
  mixture = "NB"
)
distance_models$nobo$full_model_morap_covs <- fit_gdistsamp(
  lambdas = full_model_formula_morap_covs,
  umdf = drop_na_transects(distance_models$nobo$umdf),
  mixture = "NB"
)
distance_models$nobo$intercept_model_imbcr_covs <- fit_gdistsamp_intercept(
  distance_models$nobo$full_model_imbcr_covs
)
distance_models$nobo$intercept_model_morap_covs <- fit_gdistsamp_intercept(
  distance_models$nobo$full_model_morap_covs
)
distance_models$nobo$pseudo_r_squared_n_154 <- bs_est_pseudo_rsquared(
  formula = full_model_formula_imbcr_covs,
  bird_data = distance_models$nobo$umdf,
  n = 154,
  replace = T,
  m_scale = m_scale,
  type = "gdistsamp",
  method = "deviance"
)
distance_models$nobo$pseudo_r_squared_n_84 <- bs_est_pseudo_rsquared(
  formula = full_model_formula_morap_covs,
  bird_data = drop_na_transects(distance_models$nobo$umdf),
  n = 84,
  replace = T,
  m_scale = m_scale,
  type = "gdistsamp",
  method = "deviance"
)
distance_models$nobo$cohens_d_n_154 <- bs_est_cohens_d_power(
  formula = full_model_formula_imbcr_covs,
  bird_data = distance_models$nobo$umdf,
  n = 154,
  replace = T,
  m_scale = m_scale
)
distance_models$nobo$cohens_d_n_360 <- bs_est_cohens_d_power(
  formula = full_model_formula_imbcr_covs,
  bird_data = distance_models$nobo$umdf,
  n = 360,
  replace = T,
  m_scale = m_scale
)
distance_models$nobo$cohens_d_n_540 <- bs_est_cohens_d_power(
  formula = full_model_formula_imbcr_covs,
  bird_data = distance_models$nobo$umdf,
  n = 540,
  replace = T,
  m_scale = m_scale
)
distance_models$nobo$cohens_d_n_720 <- bs_est_cohens_d_power(
  formula = full_model_formula_imbcr_covs,
  bird_data = distance_models$nobo$umdf,
  n = 720,
  replace = T,
  m_scale = m_scale
)
distance_models$nobo$cohens_f_n_154 <- bs_est_cohens_f_power(
  null_formula,
  full_model_formula_imbcr_covs,
  n = 154,
  bird_data = distance_models$nobo$umdf,
  m_scale = m_scale
)
distance_models$nobo$cohens_f_n_360 <- bs_est_cohens_f_power(
  null_formula,
  full_model_formula_imbcr_covs,
  n = 360,
  bird_data = distance_models$nobo$umdf,
  m_scale = m_scale
)
distance_models$nobo$cohens_f_n_540 <- bs_est_cohens_f_power(
  null_formula,
  full_model_formula_imbcr_covs,
  n = 540,
  bird_data = distance_models$nobo$umdf,
  m_scale = m_scale
)
distance_models$nobo$cohens_f_n_720 <- bs_est_cohens_f_power(
  null_formula,
  full_model_formula_imbcr_covs,
  n = 720,
  bird_data = distance_models$nobo$umdf,
  m_scale = m_scale
)
# HOLA
full_model_formula_imbcr_covs <- paste(
  c("poly(crp_ar, 1, raw = T) + ",
    "poly(grass_ar, 1, raw = T) + ",
    "poly(shrub_ar, 1, raw = T) + ",
    "poly(pat_ct, 1, raw = T) + ",
    "poly(mean_me_sh, 1, raw=T) + ",
    "poly(mean_me_o, 1, raw=T) + ",
    "poly(mean_ju_sh, 1, raw=T) + ",
    "poly(mean_ju_o, 1, raw=T) + ",
    "offset(log(effort))"),
  collapse = ""
)
full_model_formula_morap_covs <- paste(
  c("poly(crp_ar, 1, raw = T) + ",
    "poly(grass_ar, 1, raw = T) + ",
    "poly(shrub_ar, 1, raw = T) + ",
    "poly(pat_ct, 1, raw = T) + ",
    "poly(morap_me_sh, 1, raw=T) + ",
    "poly(morap_ju_sh, 1, raw=T) + ",
    "offset(log(effort))"),
  collapse = ""
)
null_formula <- paste(
  c("poly(crp_ar, 1, raw = T) + ",
    "poly(grass_ar, 1, raw = T) + ",
    "poly(shrub_ar, 1, raw = T) + ",
    "poly(pat_ct, 1, raw = T) + ",
    "offset(log(effort))"),
  collapse = ""
)
null_formula <- paste(
  c("poly(crp_ar, 1, raw = T) + ",
    "poly(grass_ar, 1, raw = T) + ",
    "poly(shrub_ar, 1, raw = T) + ",
    "offset(log(effort))"),
  collapse = ""
)
distance_models$hola$null_model <- fit_gdistsamp(
  lambdas = null_formula,
  umdf = distance_models$hola$umdf,
  mixture = "NB"
)
distance_models$hola$full_model_imbcr_covs <- fit_gdistsamp(
  lambdas = full_model_formula_imbcr_covs,
  umdf = distance_models$hola$umdf,
  mixture = "NB"
)
distance_models$hola$full_model_morap_covs <- fit_gdistsamp(
  lambdas = full_model_formula_morap_covs,
  umdf = drop_na_transects(distance_models$hola$umdf),
  mixture = "NB"
)
distance_models$hola$intercept_model_imbcr_covs <- fit_gdistsamp_intercept(
  distance_models$hola$full_model_imbcr_covs
)
distance_models$hola$intercept_model_morap_covs <- fit_gdistsamp_intercept(
  distance_models$hola$full_model_morap_covs
)
distance_models$hola$pseudo_r_squared_n_154 <- bs_est_pseudo_rsquared(
  formula = full_model_formula_imbcr_covs,
  bird_data = distance_models$hola$umdf,
  n = 154,
  replace = T,
  m_scale = m_scale,
  type = "gdistsamp",
  method = "deviance"
)
distance_models$hola$pseudo_r_squared_n_84 <- bs_est_pseudo_rsquared(
  formula = full_model_formula_morap_covs,
  bird_data = drop_na_transects(distance_models$hola$umdf),
  n = 154,
  replace = T,
  m_scale = m_scale,
  type = "gdistsamp",
  method = "deviance"
)
distance_models$hola$cohens_d_n_154 <- bs_est_cohens_d_power(
  formula = full_model_formula_imbcr_covs,
  bird_data = distance_models$hola$umdf,
  n = 154,
  replace = T,
  m_scale = m_scale
)
distance_models$hola$cohens_d_n_360 <- bs_est_cohens_d_power(
  formula = full_model_formula_imbcr_covs,
  bird_data = distance_models$hola$umdf,
  n = 360,
  replace = T,
  m_scale = m_scale
)
distance_models$hola$cohens_d_n_540 <- bs_est_cohens_d_power(
  formula = full_model_formula_imbcr_covs,
  bird_data = distance_models$hola$umdf,
  n = 540,
  replace = T,
  m_scale = m_scale
)
distance_models$hola$cohens_d_n_720 <- bs_est_cohens_d_power(
  formula = full_model_formula_imbcr_covs,
  bird_data = distance_models$hola$umdf,
  n = 720,
  replace = T,
  m_scale = m_scale
)
distance_models$hola$cohens_f_n_154 <- bs_est_cohens_f_power(
  null_formula,
  full_model_formula_imbcr_covs,
  n = 154,
  bird_data = distance_models$hola$umdf,
  m_scale = m_scale
)
distance_models$hola$cohens_f_n_360 <- bs_est_cohens_f_power(
  null_formula,
  full_model_formula_imbcr_covs,
  n = 360,
  bird_data = distance_models$hola$umdf,
  m_scale = m_scale
)
distance_models$hola$cohens_f_n_540 <- bs_est_cohens_f_power(
  null_formula,
  full_model_formula_imbcr_covs,
  n = 540,
  bird_data = distance_models$hola$umdf,
  m_scale = m_scale
)
distance_models$hola$cohens_f_n_720 <- bs_est_cohens_f_power(
  null_formula,
  full_model_formula_imbcr_covs,
  n = 720,
  bird_data = distance_models$hola$umdf,
  m_scale = m_scale
)
# CASP
full_model_formula_imbcr_covs <- paste(
  c("poly(crp_ar, 1, raw = T) + ",
    "poly(grass_ar, 2, raw = T) + ",
    "poly(shrub_ar, 1, raw = T) + ",
    "poly(pat_ct, 1, raw = T) + ",
    "poly(mat, 2, raw = T) + ",
    "poly(map, 1, raw = T) + ",
    "poly(mean_me_sh, 1, raw=T) + ",
    "poly(mean_me_o, 1, raw=T) + ",
    "poly(mean_ju_sh, 1, raw=T) + ",
    "poly(mean_ju_o, 1, raw=T) + ",
    "offset(log(effort))"),
  collapse = ""
)

full_model_formula_morap_covs <- paste(
  c("poly(crp_ar, 1, raw = T) + ",
    "poly(grass_ar, 2, raw = T) + ",
    "poly(shrub_ar, 1, raw = T) + ",
    "poly(pat_ct, 1, raw = T) + ",
    "poly(mat, 2, raw = T) + ",
    "poly(map, 1, raw = T) + ",
    "poly(morap_me_sh, 1, raw=T) + ",
    "poly(morap_ju_sh, 1, raw=T) + ",
    "offset(log(effort))"),
  collapse = ""
)
null_formula <- paste(
  c("poly(crp_ar, 1, raw = T) + ",
    "poly(grass_ar, 2, raw = T) + ",
    "poly(shrub_ar, 1, raw = T) + ",
    "poly(pat_ct, 1, raw = T) + ",
    "poly(mat, 2, raw = T) + ",
    "poly(map, 1, raw = T) + ",
    "offset(log(effort))"),
  collapse = ""
)
distance_models$casp$null_model <- fit_gdistsamp(
  lambdas = null_formula,
  umdf = distance_models$casp$umdf,
  mixture = "NB"
)
distance_models$casp$full_model_imbcr_covs <- fit_gdistsamp(
  lambdas = full_model_formula_imbcr_covs,
  umdf = distance_models$casp$umdf,
  mixture = "NB"
)
distance_models$casp$full_model_morap_covs <- fit_gdistsamp(
  lambdas = full_model_formula_morap_covs,
  umdf = drop_na_transects(distance_models$casp$umdf),
  mixture = "NB"
)
distance_models$casp$intercept_model_imbcr_covs <- fit_gdistsamp_intercept(
  distance_models$casp$full_model_imbcr_covs
)
distance_models$casp$intercept_model_morap_covs <- fit_gdistsamp_intercept(
  distance_models$casp$full_model_morap_covs
)
distance_models$casp$pseudo_r_squared_n_154 <- bs_est_pseudo_rsquared(
  formula = full_model_formula_imbcr_covs,
  bird_data = distance_models$casp$umdf,
  replace = T,
  n = 154,
  m_scale = m_scale,
  type = "gdistsamp"
)
distance_models$casp$pseudo_r_squared_morap_n_84 <- bs_est_pseudo_rsquared(
  formula = full_model_formula_morap_covs,
  bird_data = drop_na_transects(distance_models$casp$umdf),
  replace = T,
  n = 84,
  m_scale = m_scale,
  type = "gdistsamp"
)
distance_models$casp$cohens_d_n_154 <- bs_est_cohens_d_power(
  formula = full_model_formula_imbcr_covs,
  bird_data = distance_models$casp$umdf,
  n = 154,
  replace = T,
  m_scale = m_scale
)
distance_models$casp$cohens_d_n_360 <- bs_est_cohens_d_power(
  formula = full_model_formula_imbcr_covs,
  bird_data = distance_models$casp$umdf,
  n = 360,
  replace = T,
  m_scale = m_scale
)
distance_models$casp$cohens_d_n_540 <- bs_est_cohens_d_power(
  formula = full_model_formula_imbcr_covs,
  bird_data = distance_models$casp$umdf,
  n = 540,
  replace = T,
  m_scale = m_scale
)
distance_models$casp$cohens_d_n_720 <- bs_est_cohens_d_power(
  formula = full_model_formula_imbcr_covs,
  bird_data = distance_models$casp$umdf,
  n = 720,
  replace = T,
  m_scale = m_scale
)
distance_models$casp$cohens_f_n_154 <- bs_est_cohens_f_power(
  null_formula,
  full_model_formula_imbcr_covs,
  n = 154,
  bird_data = distance_models$casp$umdf,
  m_scale = m_scale
)
distance_models$casp$cohens_f_n_360 <- bs_est_cohens_f_power(
  null_formula,
  full_model_formula_imbcr_covs,
  n = 360,
  bird_data = distance_models$casp$umdf,
  m_scale = m_scale
)
distance_models$casp$cohens_f_n_540 <- bs_est_cohens_f_power(
  null_formula,
  full_model_formula_imbcr_covs,
  n = 540,
  bird_data = distance_models$casp$umdf,
  m_scale = m_scale
)
distance_models$casp$cohens_f_n_720 <- bs_est_cohens_f_power(
  null_formula,
  full_model_formula_imbcr_covs,
  n = 720,
  bird_data = distance_models$casp$umdf,
  m_scale = m_scale
)
# GRRO
distance_models$grro$cohens_d_n_154 <- bs_est_cohens_d_power(
  formula = full_model_formula_imbcr_covs,
  bird_data = distance_models$grro$umdf,
  n = 154,
  replace = T,
  m_scale = m_scale
)
distance_models$grro$cohens_d_n_360 <- bs_est_cohens_d_power(
  formula = full_model_formula_imbcr_covs,
  bird_data = distance_models$grro$umdf,
  n = 360,
  replace = T,
  m_scale = m_scale
)
distance_models$grro$cohens_d_n_540 <- bs_est_cohens_d_power(
  formula = full_model_formula_imbcr_covs,
  bird_data = distance_models$grro$umdf,
  n = 540,
  replace = T,
  m_scale = m_scale
)
distance_models$grro$cohens_d_n_720 <- bs_est_cohens_d_power(
  formula = full_model_formula_imbcr_covs,
  bird_data = distance_models$grro$umdf,
  n = 720,
  replace = T,
  m_scale = m_scale
)
distance_models$grro$cohens_f_n_154 <- bs_est_cohens_f_power(
  null_formula,
  full_model_formula_imbcr_covs,
  n = 154,
  bird_data = distance_models$grro$umdf,
  m_scale = m_scale
)
distance_models$grro$cohens_f_n_360 <- bs_est_cohens_f_power(
  null_formula,
  full_model_formula_imbcr_covs,
  n = 360,
  bird_data = distance_models$grro$umdf,
  m_scale = m_scale
)
distance_models$grro$cohens_f_n_540 <- bs_est_cohens_f_power(
  null_formula,
  full_model_formula_imbcr_covs,
  n = 540,
  bird_data = distance_models$grro$umdf,
  m_scale = m_scale
)
distance_models$grro$cohens_f_n_720 <- bs_est_cohens_f_power(
  null_formula,
  full_model_formula_imbcr_covs,
  n = 720,
  bird_data = distance_models$grro$umdf,
  m_scale = m_scale
)
# STFL
full_model_formula_imbcr_covs <- paste(
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
full_model_formula_morap_covs <- paste(
  c("poly(crp_ar, 1, raw = T) + ",
    "poly(grass_ar, 1, raw = T) + ",
    "poly(shrub_ar, 1, raw = T) + ",
    "poly(mat, 1, raw = T) + ",
    "poly(map, 1, raw = T) + ",
    "poly(morap_me_sh, 1, raw=T) + ",
    "poly(morap_ju_sh, 1, raw=T) + ",
    "offset(log(effort))"),
  collapse = ""
)
null_model_formula <- paste(
  c("poly(crp_ar, 1, raw = T) + ",
    "poly(grass_ar, 1, raw = T) + ",
    "poly(shrub_ar, 1, raw = T) + ",
    "poly(mat, 1, raw = T) + ",
    "poly(map, 1, raw = T) + ",
    "offset(log(effort))"),
  collapse = ""
)
distance_models$stfl$null_model <- fit_gdistsamp(
  lambdas = null_formula,
  umdf = distance_models$stfl$umdf,
  mixture = "NB"
)
distance_models$stfl$full_model_imbcr_covs <- fit_gdistsamp(
  lambdas = full_model_formula_imbcr_covs,
  umdf = distance_models$stfl$umdf,
  mixture = "NB"
)
distance_models$stfl$full_model_morap_covs <- fit_gdistsamp(
  lambdas = full_model_formula_morap_covs,
  umdf = drop_na_transects(distance_models$stfl$umdf),
  mixture = "NB"
)
distance_models$stfl$intercept_model_imbcr_covs <- fit_gdistsamp_intercept(
  distance_models$stfl$full_model_imbcr_covs
)
distance_models$stfl$intercept_model_morap_covs <- fit_gdistsamp_intercept(
  distance_models$stfl$full_model_morap_covs
)
distance_models$stfl$pseudo_r_squared_n_154 <- bs_est_pseudo_rsquared(
  formula = full_model_formula_imbcr_covs,
  bird_data = distance_models$stfl$umdf,
  replace = T,
  n = 154,
  m_scale = m_scale,
  type = "gdistsamp"
)
distance_models$stfl$pseudo_r_squared_morap_n_84 <- bs_est_pseudo_rsquared(
  formula = full_model_formula_morap_covs,
  bird_data = drop_na_transects(distance_models$stfl$umdf),
  replace = T,
  n = 84,
  m_scale = m_scale,
  type = "gdistsamp"
)
distance_models$stfl$cohens_d_n_154 <- bs_est_cohens_d_power(
  formula = full_model_formula_imbcr_covs,
  bird_data = distance_models$stfl$umdf,
  n = 154,
  replace = T,
  m_scale = m_scale
)
distance_models$stfl$cohens_d_n_360 <- bs_est_cohens_d_power(
  formula = full_model_formula_imbcr_covs,
  bird_data = distance_models$stfl$umdf,
  n = 360,
  replace = T,
  m_scale = m_scale
)
distance_models$stfl$cohens_d_n_540 <- bs_est_cohens_d_power(
  formula = full_model_formula_imbcr_covs,
  bird_data = distance_models$stfl$umdf,
  n = 540,
  replace = T,
  m_scale = m_scale
)
distance_models$stfl$cohens_d_n_720 <- bs_est_cohens_d_power(
  formula = full_model_formula_imbcr_covs,
  bird_data = distance_models$stfl$umdf,
  n = 720,
  replace = T,
  m_scale = m_scale
)
distance_models$stfl$cohens_f_n_154 <- bs_est_cohens_f_power(
  null_formula,
  full_model_formula_imbcr_covs,
  n = 154,
  bird_data = distance_models$stfl$umdf,
  m_scale = m_scale
)
distance_models$stfl$cohens_f_n_360 <- bs_est_cohens_f_power(
  null_formula,
  full_model_formula_imbcr_covs,
  n = 360,
  bird_data = distance_models$stfl$umdf,
  m_scale = m_scale
)
distance_models$stfl$cohens_f_n_540 <- bs_est_cohens_f_power(
  null_formula,
  full_model_formula_imbcr_covs,
  n = 540,
  bird_data = distance_models$stfl$umdf,
  m_scale = m_scale
)
distance_models$stfl$cohens_f_n_720 <- bs_est_cohens_f_power(
  null_formula,
  full_model_formula_imbcr_covs,
  n = 720,
  bird_data = distance_models$stfl$umdf,
  m_scale = m_scale
)
# SCQU
full_model_formula_imbcr_covs <- paste(
  c("poly(crp_ar, 1, raw = T) + ",
    "poly(grass_ar, 1, raw = T) + ",
    "poly(shrub_ar, 1, raw = T) + ",
    "poly(mat, 2, raw = T) + ",
    "poly(map, 1, raw = T) + ",
    "poly(mean_me_sh, 1, raw=T) + ",
    "poly(mean_me_o, 1, raw=T) + ",
    "poly(mean_ju_sh, 1, raw=T) + ",
    "poly(mean_ju_o, 1, raw=T) + ",
    "offset(log(effort))"),
  collapse = ""
)
full_model_formula_morap_covs <- paste(
  c("poly(crp_ar, 1, raw = T) + ",
    "poly(grass_ar, 1, raw = T) + ",
    "poly(shrub_ar, 1, raw = T) + ",
    "poly(mat, 2, raw = T) + ",
    "poly(map, 1, raw = T) + ",
    "poly(morap_me_sh, 1, raw=T) + ",
    "poly(morap_ju_sh, 1, raw=T) + ",
    "offset(log(effort))"),
  collapse = ""
)
null_formula <- paste(
  c("poly(crp_ar, 1, raw = T) + ",
    "poly(grass_ar, 1, raw = T) + ",
    "poly(shrub_ar, 1, raw = T) + ",
    "poly(mat, 2, raw = T) + ",
    "poly(map, 1, raw = T) + ",
    "offset(log(effort))"),
  collapse = ""
)
distance_models$scqu$null_model <- fit_gdistsamp(
  lambdas = null_formula,
  umdf = distance_models$scqu$umdf,
  mixture = "NB"
)
distance_models$scqu$full_model_imbcr_covs <- fit_gdistsamp(
  lambdas = full_model_formula_imbcr_covs,
  umdf = distance_models$scqu$umdf,
  mixture = "NB"
)
distance_models$scqu$full_model_morap_covs <- fit_gdistsamp(
  lambdas = full_model_formula_morap_covs,
  umdf = drop_na_transects(distance_models$scqu$umdf),
  mixture = "NB"
)
distance_models$scqu$intercept_model_imbcr_covs <- fit_gdistsamp_intercept(
  distance_models$scqu$full_model_imbcr_covs
)
distance_models$scqu$intercept_model_morap_covs <- fit_gdistsamp_intercept(
  distance_models$scqu$full_model_morap_covs
)
distance_models$scqu$pseudo_r_squared_n_154 <- bs_est_pseudo_rsquared(
  formula = full_model_formula_imbcr_covs,
  bird_data = distance_models$scqu$umdf,
  replace = T,
  n = 154,
  m_scale = m_scale,
  type = "gdistsamp"
)
distance_models$scqu$pseudo_r_squared_morap_n_84 <- bs_est_pseudo_rsquared(
  formula = full_model_formula_morap_covs,
  bird_data = drop_na_transects(distance_models$scqu$umdf),
  replace = T,
  n = 84,
  m_scale = m_scale,
  type = "gdistsamp"
)
distance_models$scqu$cohens_d_n_154 <- bs_est_cohens_d_power(
  formula = full_model_formula_imbcr_covs,
  bird_data = distance_models$scqu$umdf,
  n = 154,
  replace = T,
  m_scale = m_scale
)
distance_models$scqu$cohens_d_n_360 <- bs_est_cohens_d_power(
  formula = full_model_formula_imbcr_covs,
  bird_data = distance_models$scqu$umdf,
  n = 360,
  replace = T,
  m_scale = m_scale
)
distance_models$scqu$cohens_d_n_540 <- bs_est_cohens_d_power(
  formula = full_model_formula_imbcr_covs,
  bird_data = distance_models$scqu$umdf,
  n = 540,
  replace = T,
  m_scale = m_scale
)
distance_models$scqu$cohens_d_n_720 <- bs_est_cohens_d_power(
  formula = full_model_formula_imbcr_covs,
  bird_data = distance_models$scqu$umdf,
  n = 720,
  replace = T,
  m_scale = m_scale
)
distance_models$scqu$cohens_f_n_154 <- bs_est_cohens_f_power(
  null_formula,
  full_model_formula_imbcr_covs,
  n = 154,
  bird_data = distance_models$scqu$umdf,
  m_scale = m_scale
)
distance_models$scqu$cohens_f_n_360 <- bs_est_cohens_f_power(
  null_formula,
  full_model_formula_imbcr_covs,
  n = 360,
  bird_data = distance_models$scqu$umdf,
  m_scale = m_scale
)
distance_models$scqu$cohens_f_n_540 <- bs_est_cohens_f_power(
  null_formula,
  full_model_formula_imbcr_covs,
  n = 540,
  bird_data = distance_models$scqu$umdf,
  m_scale = m_scale
)
distance_models$scqu$cohens_f_n_720 <- bs_est_cohens_f_power(
  null_formula,
  full_model_formula_imbcr_covs,
  n = 720,
  bird_data = distance_models$scqu$umdf,
  m_scale = m_scale
)
# LOSH
full_model_formula_imbcr_covs <- paste(
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
full_model_formula_morap_covs <- paste(
  c("poly(crp_ar, 1, raw = T) + ",
    "poly(grass_ar, 1, raw = T) + ",
    "poly(shrub_ar, 1, raw = T) + ",
    "poly(mat, 1, raw = T) + ",
    "poly(map, 1, raw = T) + ",
    "poly(morap_me_sh, 1, raw=T) + ",
    "poly(morap_ju_sh, 1, raw=T) + ",
    "offset(log(effort))"),
  collapse = ""
)
null_formula <- paste(
  c("poly(crp_ar, 1, raw = T) + ",
    "poly(grass_ar, 1, raw = T) + ",
    "poly(shrub_ar, 1, raw = T) + ",
    "poly(mat, 1, raw = T) + ",
    "poly(map, 1, raw = T) + ",
    "offset(log(effort))"),
  collapse = ""
)
distance_models$losh$null_model <- fit_gdistsamp(
  lambdas = null_formula,
  umdf = distance_models$losh$umdf,
  mixture = "NB"
)
distance_models$losh$full_model_imbcr_covs <- fit_gdistsamp(
  lambdas = full_model_formula_imbcr_covs,
  umdf = distance_models$losh$umdf,
  mixture = "NB"
)
distance_models$losh$full_model_morap_covs <- fit_gdistsamp(
  lambdas = full_model_formula_morap_covs,
  umdf = drop_na_transects(distance_models$losh$umdf),
  mixture = "NB"
)
distance_models$losh$intercept_model_imbcr_covs <- fit_gdistsamp_intercept(
  distance_models$losh$full_model_imbcr_covs
)
distance_models$losh$intercept_model_morap_covs <- fit_gdistsamp_intercept(
  distance_models$losh$full_model_morap_covs
)
distance_models$losh$pseudo_r_squared_n_154 <- bs_est_pseudo_rsquared(
  formula = full_model_formula_imbcr_covs,
  bird_data = distance_models$losh$umdf,
  replace = T,
  n = 154,
  m_scale = m_scale,
  type = "gdistsamp"
)
distance_models$losh$pseudo_r_squared_morap_n_84 <- bs_est_pseudo_rsquared(
  formula = full_model_formula_morap_covs,
  bird_data = drop_na_transects(distance_models$losh$umdf),
  replace = T,
  n = 84,
  m_scale = m_scale,
  type = "gdistsamp"
)
distance_models$losh$cohens_d_n_154 <- bs_est_cohens_d_power(
  formula = full_model_formula_imbcr_covs,
  bird_data = distance_models$losh$umdf,
  n = 154,
  replace = T,
  m_scale = m_scale
)
distance_models$losh$cohens_d_n_360 <- bs_est_cohens_d_power(
  formula = full_model_formula_imbcr_covs,
  bird_data = distance_models$losh$umdf,
  n = 360,
  replace = T,
  m_scale = m_scale
)
distance_models$losh$cohens_d_n_540 <- bs_est_cohens_d_power(
  formula = full_model_formula_imbcr_covs,
  bird_data = distance_models$losh$umdf,
  n = 540,
  replace = T,
  m_scale = m_scale
)
distance_models$losh$cohens_d_n_720 <- bs_est_cohens_d_power(
  formula = full_model_formula_imbcr_covs,
  bird_data = distance_models$losh$umdf,
  n = 720,
  replace = T,
  m_scale = m_scale
)
distance_models$losh$cohens_f_n_154 <- bs_est_cohens_f_power(
  null_formula,
  full_model_formula_imbcr_covs,
  n = 154,
  bird_data = distance_models$losh$umdf,
  m_scale = m_scale
)
distance_models$losh$cohens_f_n_360 <- bs_est_cohens_f_power(
  null_formula,
  full_model_formula_imbcr_covs,
  n = 360,
  bird_data = distance_models$losh$umdf,
  m_scale = m_scale
)
distance_models$losh$cohens_f_n_540 <- bs_est_cohens_f_power(
  null_formula,
  full_model_formula_imbcr_covs,
  n = 540,
  bird_data = distance_models$losh$umdf,
  m_scale = m_scale
)
distance_models$losh$cohens_f_n_720 <- bs_est_cohens_f_power(
  null_formula,
  full_model_formula_imbcr_covs,
  n = 720,
  bird_data = distance_models$losh$umdf,
  m_scale = m_scale
)

save.image(paste("power_analysis_workflow_", as.character(cut(Sys.Date(), "month")), ".rdata", collapse="", sep=""), compress=T)