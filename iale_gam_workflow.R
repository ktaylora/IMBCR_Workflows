mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}

ct <- function(x){
  # this is ridiculous -- we need a larger sample size
  return(mean(c(mode(x),mean(x), median(x))))
} 

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

#
# MAIN
#

s_2014$response <- predict(m_negbin_intercept, type="lambda")[,1]

m_gam <- gamboost(
    as.numeric(response)~
      bbs(shrub_ar)+
      bbs(grass_ar)+
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
dev.off()
file.rename("Rplots.pdf", tolower(paste(argv[1],"_grass_ar.pdf",sep="")))


var_response_plot(m=m_gam, var="shrub_ar")
dev.off()
file.rename("Rplots.pdf", tolower(paste(argv[1],"_shrub_ar.pdf",sep="")))

var_response_plot(m=m_gam, var="pat_ct")
dev.off()
file.rename("Rplots.pdf", tolower(paste(argv[1],"_pat_ct.pdf",sep="")))

as.data.frame(varimp(m_gam))

p_gam_2014 <- suppressWarnings(predict(m_gam, newdata=s_2014@data, type="response"))
p_gam_2050 <- suppressWarnings(predict(m_gam, newdata=s_2050@data, type="response"))
p_gam_2100 <- suppressWarnings(predict(m_gam, newdata=s_2100@data, type="response"))

data.frame(
    dens_2014=ct(p_gam_2014),
    dens_2050=ct(p_gam_2014)+round(ct(p_gam_2014-p_gam_2050),2),
    dens_2100=ct(p_gam_2014)+round(ct(p_gam_2014-p_gam_2100),2)
  )
