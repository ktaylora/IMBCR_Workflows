# consider using options(error=traceback)
options(warn = -1, error=traceback)

AREA_OF_PLJV_KM <- 646895.6

r_data_files <- list.files(pattern="*.intecerpt_predictions.*rdata$")

species_pop_estimates <- do.call(rbind, lapply(
  X=r_data_files,
  FUN=function(f){
    load(f)
    if(!exists("density_ensemble")|!exists("se_ensemble")){ return(NULL) }
    return(data.frame(
        species=argv[1],
        dens=round(density_ensemble, 2),
        se=round(se_ensemble,2),
        n=round(AREA_OF_PLJV_KM*round(density_ensemble, 2)),
        n_se=round(AREA_OF_PLJV_KM*se_ensemble)
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
