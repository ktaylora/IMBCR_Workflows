est_betas_stderr <- function(model_selection_table=null, var=null){
  betas <- colnames(model_selection_table@Full)
    betas <- betas[grepl(betas, pattern="^lambda")]

  stderrs <- colnames(model_selection_table@Full)
    stderrs <- stderrs[grepl(stderrs, pattern="^SElambda")]  
  
  akaike_wt <- !is.na(
      model_selection_table@Full[, betas[ grepl(betas, pattern=var)] ]
    )
    
  akaike_wt <- mean(
      model_selection_table@Full[akaike_wt, 'AICwt' ],
      na.rm=T
    )
    
  cov_beta <- betas[grepl(betas, pattern=var)]
  cov_beta <- weighted.mean(
      model_selection_table@Full[,cov_beta], 
      weights=model_selection_table@Full$AICwt, 
      na.rm=T
    )
    
  cov_se <- stderrs[grepl(stderrs, pattern=var)]
  cov_se <- weighted.mean(
      model_selection_table@Full[,cov_se], 
      weights=model_selection_table@Full$AICwt, 
      na.rm=T
    )
  
  ret <- as.data.frame(matrix(c(cov_beta, cov_se), ncol=2))
    ret$wt <- akaike_wt
  names(ret) <- paste(var,c("_beta", "_se", "_wt"),sep="")
  return(ret)
}
#' hidden function that will accept a single unmarked model and a target 
#' covariate and generate a range of predictions. Returns the x/y values
#' to the user. This is usually called for generating plots from model
#' averaged data
predict_var_response <- function(m=NULL, var=NULL, length.out=100, type="lambda", add_exp_tail=F){
  if (! var %in% colnames(m@data@siteCovs) ) {
    stop("var=",var,"not found in siteCovariates table")
  }
  # figure out an acceptable range of x-values to predict across
  predicted_range <- range(m@data@siteCovs[,var])
  predicted_range <- seq(
      from = predicted_range[1], 
      to = predicted_range[2], 
      length.out = length.out
    )
  if(add_exp_tail){
    predicted_range <- c(
      predicted_range,
      exp(predicted_range[length(predicted_range)]*seq(0.3,5,by=0.2))
    )
  }
  # call unmarked and predict density of birds (by default)
  predicted <- sapply(
    X=1:length(predicted_range),
    FUN=function(i){
      t <- m@data@siteCovs[,!grepl(colnames(m@data@siteCovs), pattern=var)]
      t <- cbind(t, var=predicted_range[i])
      colnames(t) <- gsub(colnames(t), pattern="var", replacement=var)
      return(mean(unmarked::predict(m, newdata=t, type=type)[,1], na.rm=T))
    }
  )  
  # pretty-up a data.frame that we can work with
  ret <- data.frame(var=predicted_range, pred=predicted)
    colnames(ret) <- c(var, "pred")
  return(ret)
}
#' optionally generate a variable response plot for a user-defined covariate
#' and a list of unmarked model objects (and weights)
#' @export
plot_var_response <- function(m=NULL, var=NULL, plot=T, xlim=NULL, ylim=NULL, xlab=NULL, ylab="Density (birds/km2)", w=NULL, add_exp_tail=F, add_1_1_line=F, add_rug=F, grid=F, ...){
  # did the user provide a weight parameter?
  if(is.null(w)){
    w <- rep(1, length(m))
    warning(
        "no AIC model weighting vector was provided --",
        "assuming equal weights for all models"
      )
  }
  # build a cluster object to parallelize var prediction across all models
  cl <- parallel::makeCluster(parallel::detectCores()-1);
  parallel::clusterExport(cl, varlist=c("predict_var_response","plot_var_response"));
  # assume that the m object is a list of unmarked models
  predicted <- parallel::parLapply(
    cl=cl,
    X=m,
    fun=function(model){
      predict_var_response(model, var=var, add_exp_tail=add_exp_tail)
    }
  )
  # clean-up
  parallel::stopCluster(cl); rm(cl);
  # back-scale our
  x_var <- backscale_var(
      df=predicted[[1]], 
      var=var, 
      m_scale=m_scale
    )
  # use our weights parameter to average response across our models
  predicted <- apply(
    do.call(
        cbind, 
        lapply(predicted, FUN=function(x) matrix(ncol=1, x[,2]))
      ), 
      MARGIN=1, 
      FUN=weighted.mean, 
      w=w
  )
  # make some generic 'R' plots, if asked
  if(plot){
    dev.new()
    plot(
      predicted~x_var, 
      type="l", 
      col="white", 
      lwd=2.5, 
      xlab=xlab,
      xlim=xlim,
      ylim=ylim, 
      ylab=ylab,
      ...
    )
    if (add_1_1_line) {
      abline(lm(predicted~x_var), lty=5)
    }
    if(grid){
      grid()
      grid()
    }
    lines(predicted~x_var, type="l", col="red", lwd=2.5, xlab=xlab, ylab=ylab, ...)
    if(add_rug){
      rug(backscale_var(var=var, df=m[[1]]@data@siteCovs, m_scale))
    }
  # otherwise, boot it back to the user
  } else {
    return(data.frame(y=predicted,x=x_var))
  }
}

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



calc_density_mins <- function(){
  return(
    data.frame(
      d_2014=min(median(predicted_2014), mean(predicted_2014)),
      d_2050=min(median(predicted_2050), mean(predicted_2050)),
      d_2100=min(median(predicted_2100), mean(predicted_2100))
    ))
}

calc_sum_density_change <- function(){
  # let's find a less lossy way of capturing density changes than taking the difference of two mean's
  # or median's. Let's invert it, so the CT is taken of the difference
  mode <- function(x) {
    ux <- unique(x)
    ux[which.max(tabulate(match(x, ux)))]
  }
  ct <- function(x){
    # this is ridiculous -- we need a larger sample size
    return(mean(c(mode(x),mean(x), median(x))))
  } 
  d <- 
    data.frame(
      d_2014=0,
      d_2050=ct(predicted_2050-predicted_2014), # median decline (birds/transect)
      d_2100=ct(predicted_2100-predicted_2014) 
    )
  # now add the median difference to whatever the lowest measure of CT was
  return( d + ct(predicted_2014) )
}

y_range <- c(
    0, 
    floor(max(unmarked::predict(m_negbin_intercept, type="lambda")[,1]))
  )

#
# Grassland "Patch" Area Plot with theoretical extrapolation
#

plot_var_response(
    unmarked_models, 
    var="grass_ar", 
    xlab="Patch Area (square-km)", 
    ylab="Incidence",
    w=model_selection_table@Full$AICwt,
    log="x",
    xaxt = 'n',
    yaxt = 'n'
  )

ranges <- predict_var_response(unmarked_models[[1]], var="grass_ar")

axis(
    side=1, 
    at=log(exp(backscale_var(ranges, var="grass_ar", m_scale))), 
    labels=round(exp(backscale_var(ranges, var="grass_ar", m_scale)),2)
  )

axis(
   side=2,
   at=seq(min(ranges$pred), max(ranges$pred), length.out=10),
   labels=round(seq(min(ranges$pred), max(ranges$pred), length.out=10)/max(ranges$pred),2)
   )

#
# Grassland "Patchiness" count plot with theoretical extrapolation
#

plot_var_response(
    unmarked_models, 
    var="pat_ct", 
    xlab="Patch Count", 
    ylab="Incidence",
    w=model_selection_table@Full$AICwt,
    log="x",
    xaxt = 'n',
    yaxt = 'n'
  )

# p <- recordPlot()

ranges <- predict_var_response(unmarked_models[[1]], var="pat_ct")

axis(
    side=1, 
    at=log(exp(backscale_var(ranges, var="grass_ar", m_scale))), 
    labels=round(exp(backscale_var(ranges, var="grass_ar", m_scale)),2)
  )

axis(
   side=2,
   at=seq(min(ranges$pred), max(ranges$pred), length.out=10),
   labels=round(seq(min(ranges$pred), max(ranges$pred), length.out=10)/max(ranges$pred),2)
   )


