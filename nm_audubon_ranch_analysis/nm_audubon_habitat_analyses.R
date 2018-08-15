#' calculate the sum of detections aggregated by minute-period across all stations
#' in a transect
#' @export
calc_pooled_cluster_count_by_transect <- function(imbcr_df=NULL, four_letter_code=NULL, use_cl_count_field=F, calc_offset=T){
  four_letter_code <- toupper(four_letter_code)
  if (inherits(imbcr_df, "Spatial")) {
    imbcr_df <- imbcr_df@data
  }
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
        stations <- unique(this_transect[ this_transect$year == years[i], 'point'])
        if (calc_offset) {
          offset <- append(offset, length(stations))
        }
        # did we observe our focal species at this transect, for this year?
        if (four_letter_code %in% this_transect[ this_transect$year == years[i] , 'birdcode']) {
          # pool counts across minute periods for all stations sampled. Should
          # we use the cluster count field? If not, assume each detection
          # is a '1'
          if (!use_cl_count) {
            counts <- table(this_transect[ this_transect$year == years[i] & this_transect$birdcode == four_letter_code, 'timeperiod'])
          # if we are using cluster counts, sum all counts by minute-period
          } else {
            counts <- sapply(
              X = unique(this_transect$timeperiod),
              FUN = function(minute_period){
                return(sum(this_transect[this_transect$timeperiod == minute_period, 'cl_count'], na.rm = T))
              }
            )
            names(counts) <- as.numeric(unique(this_transect$timeperiod))
          }
          # store all minute-period counts > 0, if no observations made we
          # will retain our NA value from our empty array cast
          removal_matrix[ i, as.numeric(names(counts))] <- counts
        } 
      }
      return(list(y = removal_matrix, data=data.frame(transectnum = transect, year = years, effort = offset)))
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

distance_detections <- OpenIMBCR:::calc_dist_bins(raw_transect_data)
effort <- as.vector(OpenIMBCR:::calc_transect_effort(raw_transect_data))

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

# verify our detection function
OpenIMBCR:::plot_hn_det(
  intercept_m, 
  breaks = distance_detections$breaks
)

# censor observations that are out in the tails of our distribution
raw_transect_data <- raw_transect_data[ raw_transect_data$radialdistance <= quantile(raw_transect_data$radialdistance, p=MAXIMUM_DISTANCE_QUANTILE) , ] 

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

# estimate an adjusted cluster-count field value using the detection function fit above
raw_transect_data$cl_count <- floor(1 / per_obs_det_probabilities)

# use the adjusted point counts for removal modeling
test <- calc_pooled_cluster_count_by_transect(
  imbcr_df=raw_transect_data, 
  four_letter_code="WEME", 
  use_cl_count_field=T
)

# merge our minute intervals into two-minute intervals
test$y <- cbind(
  rowSums(test$y[, c(1:2)]), 
  rowSums(test$y[, c(3:4)]), 
  rowSums(test$y[, c(5:6)])
  )

umdf <- unmarked::unmarkedFrameMPois(
  y=test$y, 
  siteCovs=test$data,
  type="removal"
)

intercept_m <- unmarked::multinomPois(
  ~1 ~ offset(log(effort)), 
  se=T,
  umdf
)