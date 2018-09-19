#
# RUNTIME ARGUMENTS
#

ARGV = commandArgs(trailingOnly = T)
MAXIMUM_DISTANCE_QUANTILE = 0.9 # censor observations out in the right-tail
BIRD_CODE = ifelse(is.na(ARGV[1]), "WEME", toupper(ARGV[1]))

#
# LOCAL FUNCTIONS
#

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
#' Cohen's (1988) power for null and alternative models that leverages
#' residual variance explained for effects sizes for a take on the f-ratio
#' test
#' @export
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
#' calculate a deviance statistic from a count (Poisson) model,
#' e.g., : https://goo.gl/KdEUUa
#' @export
est_deviance <- function(m, method="residuals"){
  # by default, use model residual error to estimate deviance
  if (grepl(tolower(method), pattern = "resid")) {
    if ( inherits(m, "unmarkedFit") ) {
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
  # alternatively, use the log-likelihood value returned from our optimization
  # from unmarked. This approach is advocated by B. Bolker, but the likelihood
  # values returned from hierarchical models can look strange. Use this
  # with caution.
  } else if (grepl(tolower(method), pattern = "likelihood")) {
    if ( inherits(m, "unmarkedFit") ) {
      return(2*as.numeric(abs(m@negLogLik)))
    } else if ( inherits(m, "glm") ) {
      return(-2*as.numeric(logLik(m)))
    } else {
      return(NA)
    }
  }
}
#' estimate mcfadden's pseudo r-squared. Note that this function will work with
#' model objects fit with glm() in 'R', but that it is primarily intended to
#' work with model objects fit using the 'unmarked' R package. There is some
#' exception handling that has gone into working with AIC and negative
#' log-likelihood values reported by unmarked that should give you pause.
#' I try to be as verbose as I can with warnings when I fudge numbers reported
#' by unmarked models.
#' @export
est_pseudo_rsquared <- function(m=NULL, method="deviance") {
  if ( inherits(m, "unmarkedFit") ) {
    df <- m@data
    intercept_m <- try(unmarked::update(
      m,
      as.formula("~1~1"),
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
#' calculate the sum of detections aggregated by minute-period across all stations
#' in a transect. There's a fair amount of complexity baked-in here, but it
#' does what it says.
#' @export
calc_pooled_cluster_count_by_transect <- function(
  imbcr_df=NULL,
  four_letter_code=NULL,
  use_cl_count_field=F
){
  # what is the four-letter bird code that we will parse an IMBCR data.frame with?
  four_letter_code <- toupper(four_letter_code)
  if (inherits(imbcr_df, "Spatial")) {
    imbcr_df <- imbcr_df@data
  }
  # default action: process our imbcr data.frame by-transect (and year)
  transects <- unique(as.character(imbcr_df$transectnum))
  ret <- lapply(
    X = transects,
    FUN = function(transect){
      this_transect <- imbcr_df[imbcr_df$transectnum == transect, ]
      # bug-fix : drop 88 codes here (fly-over before sampling began)
      this_transect <- this_transect[ this_transect$timeperiod != 88 , ]
      years <- unique(this_transect$year)
      # cast an empty array to use for our removal counts, one row
      # for each year in our dataset
      removal_matrix <- matrix(0, ncol = 6, nrow = length(years))
      offset <- NULL
      for (i in 1:length(years)) {
        # unique stations across our transect
        stations <- unique(
            this_transect[ this_transect$year == years[i], 'point']
          )
        # note number of stations sampled for our offset
        offset <- append(offset, length(stations))
        # did we observe our focal species at this transect, for this year?
        bird_was_seen <- four_letter_code %in%
          this_transect[ this_transect$year == years[i] , 'birdcode']
        if (bird_was_seen) {
          # subset our focal transect-year for our species of interest
          match <- this_transect$year == years[i] &
            this_transect$birdcode == four_letter_code
          # pool counts across minute periods for all stations sampled. Should
          # we use the cluster count field? If not, assume each detection
          # is a '1'
          if (!use_cl_count_field) {
            counts <- table(this_transect[ match , 'timeperiod'])
          # if we are using cluster counts, sum all counts by minute-period
          } else {
            counts <- sapply(
              X = unique(this_transect$timeperiod),
              FUN = function(minute_period){
                return(sum(this_transect[ match & this_transect$timeperiod == minute_period, 'cl_count'], na.rm = T))
              }
            )
            names(counts) <- as.numeric(unique(this_transect$timeperiod))
          }
          # store all minute-period counts > 0, if no observations made we
          # will retain our NA value from our empty array cast
          removal_matrix[ i, as.numeric(names(counts))] <- counts
        }
      }
      # return a named list for THIS that we can rbind later
      return(list(
        y = removal_matrix,
        data=data.frame(transectnum = transect, year = years, effort = offset)
      ))
    }
  )
  # clean-up our list-of-lists and return to user
  return(list(
    y = do.call(rbind, lapply(ret, FUN=function(x) x$y)),
    data = do.call(rbind, lapply(ret, FUN=function(x) x$data))
  ))
}
#' use a half-normal distance function fit in unmarked to adjust our count observations
#' @export
pred_hn_det_from_distance <- function(x=NULL, dist=NULL){
  param <- exp(unmarked::coef(x, type = "det"))
  return(as.vector(unmarked:::gxhn(x = dist, param)))
}
#' Fit an intercept-only distance model in unmarked (with adjustments for effort)
#' and return the model object to the user. This is useful for extracting and
#' predicting probability of detection values using pred_hn_det_from_distance
#' @export
fit_intercept_only_distance_model <- function(raw_transect_data=NULL, verify_det_curve=F){
  # scrub the imbcr data.frame for our focal species
  distance_detections <- OpenIMBCR:::scrub_imbcr_df(
    raw_transect_data,
    four_letter_code = BIRD_CODE
  )
  # estimate effort and calculate our distance bins and dummy covariates on
  # detection and abundance
  year <- distance_detections$year[!duplicated(distance_detections$transectnum)]
  effort <- as.vector(OpenIMBCR:::calc_transect_effort(distance_detections))
  y <- OpenIMBCR:::calc_dist_bins(distance_detections)
  # build an unmarked data.frame with a column for effort
  umdf <- unmarked::unmarkedFrameDS(
    y = as.matrix(y$y),
    siteCovs = data.frame(effort = effort, year=year),
    dist.breaks = y$breaks,
    survey = "point",
    unitsIn = "m"
  )
  # model specification
  intercept_distance_m <- unmarked::distsamp(
    formula = ~1 ~as.factor(year) + offset(log(effort)),
    data = umdf,
    se = T,
    keyfun = "halfnorm",
    unitsOut = "kmsq",
    output = "abund"
  )
  # verify our detection function visually?
  if (verify_det_curve) {
    OpenIMBCR:::plot_hn_det(
      intercept_distance_m
    )
  }
  # return to user
  return(intercept_distance_m)
}

#
# MAIN
#

require(unmarked)
require(OpenIMBCR)
require(raster)

raw_transect_data <- rgdal::readOGR(
  "all_grids.json",
  verbose=F
)

FOCAL_TRANSECTS <- as.character(raw_transect_data$transectnum)
YEAR_SAMPLED <- as.character(raw_transect_data$year)

# treat transect-years as separate site-level observations
raw_transect_data$transectnum <- paste(
  FOCAL_TRANSECTS,
  YEAR_SAMPLED,
  sep = "-"
)

# fit a simple distance model that we can extract a detection function from
intercept_distance_m <- fit_intercept_only_distance_model(raw_transect_data)
# censor observations that are out in the tails of our distribution
raw_transect_data <-
  raw_transect_data[
    raw_transect_data$radialdistance <=
      quantile(raw_transect_data$radialdistance, p = MAXIMUM_DISTANCE_QUANTILE) ,
  ]
# calculate probability of detection from radial distance observations for all
# birds; we will filter this down to just our focal species when we get to
# pooling our removal data (below)
per_obs_det_probabilities <- round(sapply(
  raw_transect_data$radialdistance,
  function(x) pred_hn_det_from_distance(intercept_distance_m, dist=x)),
  2
)
# bug-fix : don't divide by zero -- this shouldn't be needed, because of
# our censoring the right-tail of our distance observations; but it's here
# just in-case
per_obs_det_probabilities[ per_obs_det_probabilities < 0.01 ] <- 0.01

# estimate an adjusted cluster-count field value using the detection
# function fit above. These will be our counts aggregate by minute period,
# adjusted for imperfect (visual) detection
raw_transect_data$cl_count <- floor(1 / per_obs_det_probabilities)

# use the adjusted point counts for removal modeling
adj_removal_detections <- calc_pooled_cluster_count_by_transect(
  imbcr_df = raw_transect_data,
  four_letter_code = BIRD_CODE,
  use_cl_count_field = T
)

# merge our minute intervals into two-minute intervals
adj_removal_detections$y <- cbind(
  rowSums(adj_removal_detections$y[, c(1:2)]),
  rowSums(adj_removal_detections$y[, c(3:4)]),
  rowSums(adj_removal_detections$y[, c(5:6)])
)

# tack-on a categorical "ranch status" variable to use for the modeling
adj_removal_detections$data$ranch_status <-
  as.numeric(grepl(adj_removal_detections$data$transectnum, pattern="RANCH"))

umdf <- unmarked::unmarkedFrameMPois(
  y = adj_removal_detections$y,
  siteCovs = adj_removal_detections$data,
  type = "removal"
)

intercept_adj_removal_m <- unmarked::multinomPois(
  ~1 ~ as.factor(year) + offset(log(effort)),
  se = T,
  umdf
)

ranch_status_adj_removal_m <- unmarked::multinomPois(
  ~1 ~ as.factor(ranch_status) + as.factor(year) + offset(log(effort)),
  se = T,
  umdf
)
# propotion of variance explained by adding our ranch covariate?
cat(" -- null model r-squared:",
    est_pseudo_rsquared(intercept_adj_removal_m),
    "\n"
)
cat(" -- alternative model r-squared:",
    est_pseudo_rsquared(ranch_status_adj_removal_m),
    "\n"
)
# add-in our habitat covariates, calculated by-grid unit
FOCAL_TRANSECTS <- sapply(
  adj_removal_detections$data$transectnum,
  FUN=function(x) paste(unlist(strsplit(as.character(x), split="-"))[1:3], collapse="-")
)

YEAR_SAMPLED <- sapply(
  adj_removal_detections$data$transectnum,
  FUN=function(x) unlist(strsplit(as.character(x), split="-"))[4]
)

usng_units <- OpenIMBCR:::readOGRfromPath(
  paste("/home/ktaylora/Projects/nm_audubon_habitat_modeling/study_",
  "region_convex_hull_usng_units.shp", sep="")
)
all_imbcr_transects <- OpenIMBCR:::readOGRfromPath(
  paste("/home/ktaylora/Projects/nm_audubon_habitat_modeling/all_",
  "transects_for_mapping.shp", sep="")
)
# spatial join the 1-km2 USNG grid units with the IMBCR transects
# we are using for our analysis
transect_usng_units <- OpenIMBCR:::spatial_join(
  usng_units,
  all_imbcr_transects
)
# select for transects used in our original dataset -- some transects
# are sampled more than one year
match <- as.vector(unlist(sapply(
  FOCAL_TRANSECTS,
  FUN=function(transect) min(which(as.character(transect_usng_units$trnsctn) == transect))))
)
transect_usng_units <- transect_usng_units[ match , ]
# 2016 NASS-CDL is consistent with when sampling began on Audubon ranch
# transects
usda_nass <- raster::raster(
  "/gis_data/Landcover/NASS/Raster/2016_30m_cdls.tif"
)
# define the NASS-CDL values we are attributing as "habitat"
cat(" -- calculating habitat composition/configuration metrics\n")
area_statistics <-
  data.frame(
      field_name=c(
        'grass_ar',
        'shrub_ar',
        'wetland_ar'
      ),
      src_raster_value=c(
        '176',
        'c(64,152)',
        'c(195,190)'
      )
    )
configuration_statistics <- c(
    'pat_ct'
)
# process our NASS-CDL composition statistics iteratively
usda_nass_by_unit <- OpenIMBCR:::extract_by(
  lapply(1:nrow(transect_usng_units), FUN=function(i) transect_usng_units[i,]),
  usda_nass
)
for(i in 1:nrow(area_statistics)){
  focal <- OpenIMBCR:::binary_reclassify(
    usda_nass_by_unit,
    from=eval(parse(text=as.character(area_statistics[i, 2])))
  )
  transect_usng_units@data[, as.character(area_statistics[i, 1])] <-
    sapply(X=focal, FUN=OpenIMBCR:::calc_total_area)
}
# ditto for configuration statistics
cat(" -- building a habitat/not-habitat raster surface\n")
valid_habitat_values <- eval(parse(
    text=paste("c(",paste(area_statistics$src_raster_value[
      !grepl(area_statistics$field_name, pattern="rd_ar")
    ], collapse = ","), ")", sep="")
))
cat(" -- calculating patch configuration metrics\n")

focal <- OpenIMBCR:::binary_reclassify(
  usda_nass_by_unit,
  from=valid_habitat_values
)

transect_usng_units@data[, as.character(configuration_statistics[1])] <-
  unlist(sapply(X=focal, FUN=OpenIMBCR:::calc_patch_count))


# fit our habitat model

# these are the columns that we are going to merge-in for exposure to
# our models
cols <- c("grass_ar","shrub_ar","wetland_ar","pat_ct")
transect_usng_units@data <- transect_usng_units@data[,cols]
m_scale <- scale(transect_usng_units@data)
transect_usng_units@data <- data.frame(scale(transect_usng_units@data))

# merge-in our data.frame of site-covs
adj_removal_detections$data <- cbind(
  adj_removal_detections$data,
  transect_usng_units@data[,cols]
)

umdf <- unmarked::unmarkedFrameMPois(
  y = adj_removal_detections$y,
  siteCovs = adj_removal_detections$data,
  type = "removal"
)

full_model_formula <- paste(
  "~1 ~grass_ar+shrub_ar+pat_ct+as.factor(year)+as.factor(ranch_status)+offset(log(effort))"
)

full_model_minus_ranch_cov_formula <- paste(
  "~1 ~grass_ar+shrub_ar+pat_ct+as.factor(year)+offset(log(effort))"
)

full_model_ranch_status_adj_removal_m <- unmarked::multinomPois(
  as.formula(full_model_formula),
  se = T,
  umdf
)

full_model_adj_removal_m <- unmarked::multinomPois(
  as.formula(full_model_minus_ranch_cov_formula),
  se = T,
  umdf
)
# density
mean_density <- median(
  unmarked::predict(full_model_adj_removal_m, type="state")[,1])
# se
mean_density_se <- median(
  unmarked::predict(full_model_adj_removal_m, type="state")[,2])
# population size estimates
regional_pop_size_est <- round(218196 * mean_density)
regional_pop_size_est_se <- round(218196 * mean_density_se)
ranch_pop_size_est <- round(72.843416 * mean_density)
ranch_pop_size_est_se <- round(72.843416 * mean_density_se)

cat(
  " -- density:",
  mean_density, "(", mean_density_se,")\n",
  sep=""
)
cat(" -- regional pop:",
  round(regional_pop_size_est/1000000, 1),
  "(", round(regional_pop_size_est_se/1000000, 1),")\n",
  sep=""
)
cat(" -- ranch pop:",
  ranch_pop_size_est, "(", ranch_pop_size_est_se,")\n",
  sep=""
)

# flush our session to disk and exit
r_data_file <- tolower(paste(
  tolower(BIRD_CODE),
  "_imbcr_hinge_modeling_workflow_",
  gsub(format(Sys.time(), "%b %d %Y"), pattern = " ", replacement = "_"),
  ".rdata",
  sep = ""
))
save(
  compress = T,
  list = ls(),
  file = r_data_file
)
