est_betas_stderr <- function(model_selection_table=null, var=null){
  betas <- colnames(model_selection_table@Full)
    betas <- betas[grepl(betas, pattern="^lambda")]

  stderrs <- colnames(model_selection_table@Full)
    stderrs <- stderrs[grepl(stderrs, pattern="^SElambda")]  
    
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
    names(ret) <- paste(var,c("_beta", "_se"),sep="")

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


