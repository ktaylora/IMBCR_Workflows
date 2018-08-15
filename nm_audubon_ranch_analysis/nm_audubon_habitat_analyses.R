argv <- commandArgs(trailingOnly = T)
# calculate a fitted poisson regression deviance statistic,
# e.g., : https://goo.gl/KdEUUa
est_deviance <- function(m, method="residuals"){
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
#' estimate McFadden's pseudo r-squared
#' Note that this function will work with model objects fit with glm() in 'R',
#' but is primarily intended to work with model objects fit using the 'unmarked'
#' R package. There is some fairly-involved exception handling that has gone
#' into working with AIC and negative log-likelihood values reported by
#' unmarked that should give you pause. I try to be as verbose as I can
#' with warnings when I fudge numbers reported by unmarked models.
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
      # bug-fix : drop 88 codes here
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
          if (!use_cl_count) {
            counts <- table(this_transect[ match , 'timeperiod'])
          # if we are using cluster counts, sum all counts by minute-period
          } else {
            counts <- sapply(
              X = unique(this_transect$timeperiod),
              FUN = function(minute_period){
                return(sum(this_transect[ match , 'cl_count'], na.rm = T))
              }
            )
            names(counts) <- as.numeric(unique(this_transect$timeperiod))
          }
          # store all minute-period counts > 0, if no observations made we
          # will retain our NA value from our empty array cast
          removal_matrix[ i, as.numeric(names(counts))] <- counts
        } 
      }
      return(list(
        y = removal_matrix, 
        data=data.frame(transectnum = transect, year = years, effort = offset)
      ))
    }
  )
  return(list(
    y = do.call(rbind, lapply(ret, FUN=function(x) x$y)),
    data = do.call(rbind, lapply(ret, FUN=function(x) x$data))
  ))
}
#' use a half-normal distance function fit in unmarked to adjust our count observations
pred_hn_det_from_distance <- function(x=NULL, dist=NULL){
  param <- exp(unmarked::coef(x, type = "det"))
  return(as.vector(unmarked:::gxhn(x = dist, param)))
}
#
# MAIN
#

MAXIMUM_DISTANCE_QUANTILE = 0.9 # censor observations that are way out in the shoulder
BIRD_CODE = ifelse(is.null(argv[1]),"WEME", toupper(argv[1]))

raw_transect_data <- rgdal::readOGR(
  "all_grids.json"
)

# treat transect-years as separate observations 
raw_transect_data$transectnum <- paste(
  as.character(raw_transect_data$transectnum), 
  as.character(raw_transect_data$year), 
  sep="-"
)

# fit an intercept-only detection function in unmarked

distance_detections <- OpenIMBCR:::scrub_imbcr_df(
  raw_transect_data, 
  four_letter_code = "WEME"
)

effort <- as.vector(OpenIMBCR:::calc_transect_effort(distance_detections))
distance_detections <- OpenIMBCR:::calc_dist_bins(distance_detections)


umdf <- unmarked::unmarkedFrameDS(
  y = as.matrix(distance_detections$y),
  siteCovs = data.frame(effort = effort),
  dist.breaks = distance_detections$breaks,
  survey = "point",
  unitsIn = "m"
)

intercept_distance_m <- unmarked::distsamp(
  formula = ~1 ~1+offset(log(effort)),
  data = umdf,
  se = T,
  keyfun = "halfnorm",
  unitsOut = "kmsq",
  output = "abund"
)

# verify our detection function visually
OpenIMBCR:::plot_hn_det(
  intercept_distance_m, 
  breaks = distance_detections$breaks
)

# censor observations that are out in the tails of our distribution
raw_transect_data <- 
  raw_transect_data[ 
    raw_transect_data$radialdistance <= 
      quantile(raw_transect_data$radialdistance, p=MAXIMUM_DISTANCE_QUANTILE) , 
  ] 

# calculate probability of detection from radial distance observations for all 
# birds, we will filter this down to just our focal species when we get to 
# pooling our removal data (below)
per_obs_det_probabilities <- round(sapply(
  raw_transect_data$radialdistance,
  function(x) pred_hn_det_from_distance(intercept_m, dist=x)),
  2
)

# bug-fix : divide by zero
per_obs_det_probabilities[ per_obs_det_probabilities < 0.01 ] <- 0.01

# estimate an adjusted cluster-count field value using the detection 
# function fit above
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
  as.factor(
    letters[ as.numeric(grepl(adj_removal_detections$data$transectnum, pattern="RANCH"))+1 ]
  )

umdf <- unmarked::unmarkedFrameMPois(
  y = adj_removal_detections$y, 
  siteCovs = adj_removal_detections$data,
  type = "removal"
)

intercept_adj_removal_m <- unmarked::multinomPois(
  ~1 ~ offset(log(effort)), 
  se = T,
  umdf
)

ranch_status_adj_removal_m <- unmarked::multinomPois(
  ~1 ~ as.factor(ranch_status) + offset(log(effort)), 
  se = T,
  umdf
)

# propotion of variance explained by adding our ranch covariate?
est_pseudo_rsquared(ranch_status_adj_removal_m)