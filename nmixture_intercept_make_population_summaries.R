# consider using options(error=traceback)
options(warn = -1, error=traceback)

AREA_OF_PLJV_KM <- 646895.6

r_data_files <- list.files(pattern="*.intecerpt_predictions.*rdata$")

species_pop_estimates <- do.call(rbind, lapply(
  X=r_data_files,
  FUN=function(f){
    load(f)
    return(data.frame(
        species=argv[1],
        negbin_dens=round(m_negbin_predicted$est,2),
        negbin_se=round(m_negbin_predicted$se, 2),
        pois_dens=round(m_pois_predicted$est,2),
        pois_se=round(m_pois_predicted$se,2),
        n_negbin=round(AREA_OF_PLJV_KM*m_negbin_predicted$est),
        n_se_negbin=round(AREA_OF_PLJV_KM*m_negbin_predicted$se),
        n_pois=round(AREA_OF_PLJV_KM*m_pois_predicted$est),
        n_se_pois=round(AREA_OF_PLJV_KM*m_pois_predicted$se)
      ))
  }
))

filename <- commandArgs(trailingOnly=T)

if(length(filename)>0){
  filename <- filename
} else {
  filename <- "population_estimates.csv"
}

write.csv(species_pop_estimates, filename, row.names=F)
