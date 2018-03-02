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

# consider using options(error=traceback)
options(warn = -1, error=traceback)

argv <- commandArgs(trailingOnly=T)

  
#
# Local accessory functions, some of which may overload what's
# saved in the loaded Rdata file
#

plot_hn_det <- function(x=NULL, breaks=NULL){
  param <- exp(unmarked::coef(x, type = "det"))
  plot(
      function(x) gxhn(x, param), 0, max(breaks),
  	  xlab = "Distance (m)", ylab = "Detection probability"
    )
  grid(); grid();
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

#
# MAIN
#

argv <- commandArgs(trailingOnly = T)

if(length(argv) < 2){
  warning("using default (2016) IMBCR CSV")
  imbcr_csv <- "/global_workspace/imbcr_number_crunching/results/RawData_PLJV_IMBCR_20161201.csv"
} else {
  imbcr_csv <- argv[2]
}

s <- OpenIMBCR:::scrub_imbcr_df(
    OpenIMBCR:::imbcrTableToShapefile(
        imbcr_csv
      ),
    four_letter_code = toupper(argv[1])
  )

detections <- OpenIMBCR:::calc_dist_bins(s)
effort     <- as.vector(OpenIMBCR:::calc_transect_effort(s))

umdf <- unmarked::unmarkedFrameDS(
    y=as.matrix(detections$y),
    siteCovs=data.frame(effort=effort),
    dist.breaks=detections$breaks,
    survey="point",
    unitsIn="m"
  )

m_pois_intercept_model <- unmarked::distsamp(
    formula = ~1 ~1+offset(log(effort)),
    data = umdf,
    se = T,
    keyfun = "halfnorm",
    unitsOut = "kmsq",
    output = "density"
  )

m_pois_predicted <- median(unmarked::predict(
      m_pois_intercept_model, type="lambda")
    )
    
m_pois_predicted <- data.frame(
    est=median(m_pois_predicted[,1]),
    se=median(m_pois_predicted[,2])
  )

umdf <- unmarked::unmarkedFrameGDS(
  y=as.matrix(detections$y),
  siteCovs=data.frame(effort=effort),
  dist.breaks=detections$breaks,
  numPrimary=1,
  survey="point",
  unitsIn="m"
)

m_negbin_intercept_model <- unmarked::gdistsamp(
    lambdaformula = ~1+offset(log(effort)),
    pformula = ~1,
    phiformula = ~1,
    data = umdf,
    se = T,
    keyfun = "halfnorm",
    unitsOut = "kmsq",
    output = "density"
  )

m_negbin_predicted <- median(unmarked::predict(
      m_negbin_intercept_model, type="lambda")
    )
    
m_negbin_predicted <- data.frame(
    est=median(m_negbin_predicted[,1]),
    se=median(m_negbin_predicted[,2])
  )

r_data_file <- tolower(paste(argv[1],
      "_imbcr_intecerpt_predictions_",
      gsub(format(Sys.time(), "%b %d %Y"), pattern=" ", replacement="_"),
      ".rdata",
      sep="")
  )

save(
    compress=T,
    list=(ls(pattern="[a-z]")),
    file=r_data_file
  )
