pred_hn_det_from_distance <- function(x=NULL, dist=NULL){
  param <- exp(coef(x, type = "det"))
  return(as.vector(unmarked:::gxhn(x=dist, param)))
}

calc_transect_summary_detections <- function(s=NULL, name=NULL, field='est_abund'){
  return(OpenIMBCR:::calc_route_centroids(
      s = s,
      four_letter_code = toupper(name),
      use='est_abund'
    ))
}

#
# MAIN
#

require(unmarked)
require(rgdal)

argv <- commandArgs(trailingOnly = T)

s_2016 <- OpenIMBCR:::scrub_imbcr_df(
    OpenIMBCR:::imbcrTableToShapefile(
        "/global_workspace/imbcr_number_crunching/results/RawData_PLJV_IMBCR_20161201.csv"
      ),
    four_letter_code = toupper(argv[1])
  )

detections <- OpenIMBCR:::calc_dist_bins(s_2016)
effort     <- as.vector(OpenIMBCR:::calc_transect_effort(s_2016))

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
    s_2016$radialdistance,
    function(x) pred_hn_det_from_distance(intercept_m, dist=x)),
    2
  )

# calculate an estimate of abundance (accounting for p-det)
s_2016$est_abund <- round(s_2016$radialdistance > 0)
  
s_2016 <- calc_transect_summary_detections(
    s=s_2016,
    name=toupper(argv[1]),
    field='est_abund'
  )

#
# 2017
#

s_2017 <- OpenIMBCR:::scrub_imbcr_df(
    OpenIMBCR:::imbcrTableToShapefile(
        "/global_workspace/imbcr_number_crunching/results/RawData_PLJV_IMBCR_20171017.csv"
      ),
    four_letter_code = toupper(argv[1])
  )

detections <- OpenIMBCR:::calc_dist_bins(s_2017)
effort     <- as.vector(OpenIMBCR:::calc_transect_effort(s_2017))

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
    s_2017$radialdistance,
    function(x) pred_hn_det_from_distance(intercept_m, dist=x)),
    2
  )

# calculate an estimate of abundance (accounting for p-det)
s_2017$est_abund <- round(s_2017$radialdistance > 0)
  
s_2017 <- calc_transect_summary_detections(
    s=s_2017,
    name=toupper(argv[1]),
    field='est_abund'
  )


#
# merge
#

s_2016$year <- 2016
s_2017$year <- 2017

s <- rbind(s_2016, s_2017)

rgdal::writeOGR(
    s, 
    ".", 
    paste(tolower(argv[1]),"_transect_detections", sep=""), 
    driver="ESRI Shapefile", 
    overwrite=T
  )




