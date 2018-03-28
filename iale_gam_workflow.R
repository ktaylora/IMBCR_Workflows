var_response_plot <- function(m=NULL, var=NULL, rug=T){
  dev.new();
  plot(
    m, 
    which=var, 
    type="l", 
    col="red", 
    lwd=1.5,
    rug=rug,
  );
  grid();
  grid();
}

mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}
ct <- function(x){
  # this is ridiculous -- we need a larger sample size
  return(mean(c(mode(x),mean(x), median(x))))
} 

#
# MAIN
#

require(mboost)

s_2014$response <- unmarked::predict(m_negbin_intercept, type="lambda")[,1]

m_gam <- mboost::gamboost(
    as.numeric(response)~
      bmono(shrub_ar, constraint="concave")+
      bmono(grass_ar, constraint="increasing") +
      bbs(pat_ct), 
    control = boost_control(center=F), 
    family=Gaussian(), 
    offset=log(effort),
    data=s_2014@data[,c('response',vars)]
  )

suppressWarnings(data.frame(
  mean_dens=ct(predict(m_gam, newdata=s_2014@data, type="response")),
  mean_error=ct(s_2014$response - predict(m_gam, newdata=s_2014@data, type="response")),
  cint=plotrix::std.error(s_2014$response - predict(m_gam, newdata=s_2014@data, type="response"))*1.96,
  pi_range=diff(range(s_2014$response - predict(m_gam, newdata=s_2014@data, type="response")))
))


var_response_plot(m=m_gam, var="grass_ar")
#dev.off()
#file.rename("Rplots.pdf", tolower(paste(argv[1],"_grass_ar.pdf",sep="")))


var_response_plot(m=m_gam, var="shrub_ar")
#dev.off()
#file.rename("Rplots.pdf", tolower(paste(argv[1],"_shrub_ar.pdf",sep="")))

var_response_plot(m=m_gam, var="pat_ct")
#dev.off()
#file.rename("Rplots.pdf", tolower(paste(argv[1],"_pat_ct.pdf",sep="")))

m_gam_importance_table <- as.data.frame(varimp(m_gam, percent=T))

p_gam_2014 <- suppressWarnings(predict(m_gam, newdata=s_2014@data, type="response"))
p_gam_2050 <- suppressWarnings(predict(m_gam, newdata=s_2050@data, type="response"))
p_gam_2100 <- suppressWarnings(predict(m_gam, newdata=s_2100@data, type="response"))

m_gam_prediction_table <- data.frame(
    dens_2014=ct(p_gam_2014),
    dens_2050=ct(p_gam_2014)+round(ct(p_gam_2014-p_gam_2050),2),
    dens_2100=ct(p_gam_2014)+round(ct(p_gam_2014-p_gam_2100),2)
  )
  
m_gam_prediction_table$perc_reduction <- round(
    1 - (m_gam_prediction_table$dens_2100/m_gam_prediction_table$dens_2014) , 3
  )
  
r_data_file <- tolower(paste(
      tolower(argv[1]),
      "_imbcr_gam_iale_workflow_",
      gsub(format(Sys.time(), "%b %d %Y"), pattern = " ", replacement = "_"),
      ".rdata",
      sep = ""
    ))
    
save(
    compress=T,
    list=(ls(pattern="[a-z]")),
    file=r_data_file
  )
