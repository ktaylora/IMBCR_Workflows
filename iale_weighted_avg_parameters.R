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
