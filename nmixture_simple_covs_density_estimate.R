#
# This is not pretty -- I'm sorry; Kyle (2018)
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
mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}
ct <- function(x){
  # this is ridiculous -- we need a larger sample size
  return(mean(c(mode(x),mean(x), median(x))))
} 
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

umdf <- unmarked::unmarkedFrameGDS(
  y=as.matrix(detections$y),
  siteCovs=data.frame(
    lon=sapply(unique(s$transectnum), FUN=function(x){ mean(s@data[s$transectnum == x, 'lon']) }), 
    lat=sapply(unique(s$transectnum), FUN=function(x){ mean(s@data[s$transectnum == x, 'lat']) }), 
    effort=effort),
  dist.breaks=detections$breaks,
  numPrimary=1,
  survey="point",
  unitsIn="m"
)

#
# first, fit for the standard poisson mixture distribution
#
m_pois_intercept_model <- try(unmarked::gdistsamp(
    lambdaformula = ~1+offset(log(effort)),
    pformula = ~1,
    phiformula = ~1,
    data = umdf,
    se = T,
    keyfun = "halfnorm",
    unitsOut = "kmsq",
    output = "density",
    mixture="P"
  ))

m_pois_alternative_model_simple <- try(unmarked::gdistsamp(
    lambdaformula = ~poly(lat,1)+poly(lon,1)+offset(log(effort)),
    pformula = ~1,
    phiformula = ~1,
    data = umdf,
    se = T,
    keyfun = "halfnorm",
    unitsOut = "kmsq",
    output = "density",
    mixture="P"
  ))

m_pois_alternative_model_complex <- try(unmarked::gdistsamp(
    lambdaformula = ~poly(lat,4)+poly(lon,4)+offset(log(effort)),
    pformula = ~1,
    phiformula = ~1,
    data = umdf,
    se = T,
    keyfun = "halfnorm",
    unitsOut = "kmsq",
    output = "density",
    mixture="P"
  ))

if(class(m_pois_intercept_model) != "try-error"){
  m_pois_int_predicted <- unmarked::predict(
    m_pois_intercept_model, 
    type="lambda"
  )
  m_pois_int_predicted <- data.frame(
    est=ct(m_pois_int_predicted[,1]),
    se=ct(m_pois_int_predicted[,2])
  )
} else {
  m_pois_int_predicted <- data.frame(est=NA, se=NA)
}

if(class(m_pois_alternative_model_simple) != "try-error"){
  m_pois_alt_simple_predicted <- unmarked::predict(
    m_pois_alternative_model_simple, 
    type="lambda"
  )
  m_pois_alt_simple_predicted <- data.frame(
    est=ct(m_pois_alt_simple_predicted[,1]),
    se=ct(m_pois_alt_simple_predicted[,2])
  )
} else {
  m_pois_alt_simple_predicted <- data.frame(est=NA, se=NA)
}


if(class(m_pois_alternative_model_complex) != "try-error"){
  m_pois_alt_complex_predicted <- unmarked::predict(
    m_pois_alternative_model_complex, 
    type="lambda"
  )
  m_pois_alt_complex_predicted <- data.frame(
    est=ct(m_pois_alt_complex_predicted[,1]),
    se=ct(m_pois_alt_complex_predicted[,2])
  )
} else {
  m_pois_alt_complex_predicted <- data.frame(est=NA, se=NA)
}


#
# Now fit for the negative binomial mixture distribution
#

m_negbin_intercept_model <- unmarked::gdistsamp(
    lambdaformula = ~1+offset(log(effort)),
    pformula = ~1,
    phiformula = ~1,
    data = umdf,
    se = T,
    keyfun = "halfnorm",
    unitsOut = "kmsq",
    output = "density",
    mixture="NB"
  )

m_negbin_alternative_model_simple <- unmarked::gdistsamp(
    lambdaformula = ~poly(lat,1)+poly(lon,1)+offset(log(effort)),
    pformula = ~1,
    phiformula = ~1,
    data = umdf,
    se = T,
    keyfun = "halfnorm",
    unitsOut = "kmsq",
    output = "density",
    mixture="NB"
  )
  
m_negbin_alternative_model_complex <- unmarked::gdistsamp(
    lambdaformula = ~poly(lat,4)+poly(lon,4)+offset(log(effort)),
    pformula = ~1,
    phiformula = ~1,
    data = umdf,
    se = T,
    keyfun = "halfnorm",
    unitsOut = "kmsq",
    output = "density",
    mixture="NB"
  )


if(class(m_negbin_intercept_model) != "try-error"){
  m_negbin_int_predicted <- unmarked::predict(
    m_negbin_intercept_model, 
    type="lambda"
  )  
  m_negbin_int_predicted <- data.frame(
    est=ct(m_negbin_int_predicted[,1]),
    se=ct(m_negbin_int_predicted[,2])
  )
} else {
  m_negbin_int_predicted <- data.frame(est=NA, se=NA)
}

if(class(m_negbin_alternative_model_simple) != "try-error"){
  m_negbin_alt_simple_predicted <- unmarked::predict(
    m_negbin_alternative_model_simple, 
    type="lambda"
  )
  m_negbin_alt_simple_predicted <- data.frame(
    est=ct(m_negbin_alt_simple_predicted[,1]),
    se=ct(m_negbin_alt_simple_predicted[,2])
  )
} else {
  m_negbin_alt_simple_predicted <- data.frame(est=NA, se=NA)
}

if(class(m_negbin_alternative_model_simple) != "try-error"){
  m_negbin_alt_complex_predicted <- unmarked::predict(
    m_negbin_alternative_model_complex, 
    type="lambda"
  )
  m_negbin_alt_complex_predicted <- data.frame(
    est=ct(m_negbin_alt_complex_predicted[,1]),
    se=ct(m_negbin_alt_complex_predicted[,2])
  )
} else {
  m_negbin_alt_complex_predicted <- data.frame(est=NA, se=NA)
}  

# average our predictions across mixture and covariates using AIC

densities <- c(
    m_pois_int_predicted$est, 
    m_pois_alt_simple_predicted$est,
    m_pois_alt_complex_predicted$est,
    m_negbin_int_predicted$est,
    m_negbin_alt_simple_predicted$est,
    m_negbin_alt_complex_predicted$est
  )

std_errors <- c(
    m_pois_int_predicted$se, 
    m_pois_alt_simple_predicted$se,
    m_pois_alt_complex_predicted$se,
    m_negbin_int_predicted$se,
    m_negbin_alt_simple_predicted$se,
    m_negbin_alt_complex_predicted$se
  )


aic <- mean(c(
    ifelse(class(m_pois_intercept_model) == "try-error", NA, m_pois_intercept_model@AIC),
    ifelse(class(m_pois_alternative_model_simple) == "try-error", NA, m_pois_alternative_model_simple@AIC), 
    ifelse(class(m_pois_alternative_model_complex) == "try-error", NA, m_pois_alternative_model_complex@AIC), 
    ifelse(class(m_negbin_intercept_model) == "try-error", NA, m_negbin_intercept_model@AIC), 
    ifelse(class(m_negbin_alternative_model_simple) == "try-error", NA, m_negbin_alternative_model_simple@AIC), 
    ifelse(class(m_negbin_alternative_model_complex) == "try-error", NA, m_negbin_alternative_model_complex@AIC)
  ),
  na.rm=T
)

weights <- OpenIMBCR:::akaike_weights(
  aic_values=c(
    ifelse(class(m_pois_intercept_model) == "try-error", 0, m_pois_intercept_model@AIC),
    ifelse(class(m_pois_alternative_model_simple) == "try-error", 0, m_pois_alternative_model_simple@AIC), 
    ifelse(class(m_pois_alternative_model_complex) == "try-error", 0, m_pois_alternative_model_complex@AIC), 
    ifelse(class(m_negbin_intercept_model) == "try-error", 0, m_negbin_intercept_model@AIC), 
    ifelse(class(m_negbin_alternative_model_simple) == "try-error", 0, m_negbin_alternative_model_simple@AIC), 
    ifelse(class(m_negbin_alternative_model_complex) == "try-error", 0, m_negbin_alternative_model_complex@AIC)
  ), 
  precision=10
)

density_ensemble <- weighted.mean(densities, weights=weights, na.rm=T)
se_ensemble <- weighted.mean(std_errors, weights=weights, na.rm=T)

print(data.frame(spp=argv[1], density=density_ensemble, se=se_ensemble, aic=aic))
