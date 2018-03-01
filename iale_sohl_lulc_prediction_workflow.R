require(raster)
require(rgdal)

#
# Runtime Parameters
#

AIC_RESCALE_CONST         <- 100000
AIC_SUBSTANTIAL_THRESHOLD <- 8
CORRELATION_THRESHOLD     <- 0.65
NORMALIZE                 <- T

ggplot2_multivariate_densities <- function(densities=NULL, var=NULL, correction=1, ylab=NULL, xlab=NULL){
  # color brewer colors
  xlim <- range(densities[,2]) / correction

  xlim[1] <- xlim[1] - abs(xlim[1])*0.1
  xlim[2] <- xlim[2] + abs(xlim[2])*0.1

  rugs <- densities
  rugs[,2] <- rugs[,2] / correction

  densities$year <- as.factor(densities$year)
  rugs$year <- as.factor(rugs$year)

  gg_plot <- ggplot(densities, aes_string(var, fill = 'year', col = 'year')) + geom_density(alpha=0.65)
  gg_plot <- gg_plot + scale_x_continuous(limits=xlim)
  gg_plot <- gg_plot + geom_rug(aes_string(var, color='year'), sides="b", size=0.25, alpha=0.8, data=rugs)
  gg_plot <- gg_plot + xlab(xlab) + ylab(ylab)
  gg_plot <- gg_plot + scale_fill_brewer(type = 'seq', palette="YlGnBu") + scale_color_brewer(type = 'seq', palette="YlGnBu")

  gg_plot +
    theme(
        axis.text = element_text(face="bold",size=10),
        text = element_text(face="bold",size=12),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        #panel.background = element_blank(),
        #panel.grid.minor = element_blank(),
        #panel.grid.major = element_line("white"),
        axis.line=element_line("grey50")
      )
}

mean_normalization <- function(x, alternative_means=NULL){
  x <- x/mean(x) *
    min(alternative_means)
}

#
# Main
#

load(commandArgs(trailingOnly=T))

r_data_file <- tolower(paste(
      tolower(argv[1]),
      "_imbcr_hds_iale_prediction_workflow_",
      gsub(format(Sys.time(), "%b %d %Y"), pattern = " ", replacement = "_"),
      ".rdata",
      sep = ""
    ))

r_gcam_45_2014 <- raster::raster(
    "/global_workspace/terry_sohl_sres_btu_30m_landcover_predictions/extracted/gplcc_gcam_45_2014.tif"
  )

r_gcam_45_2050 <- raster::raster(
    "/global_workspace/terry_sohl_sres_btu_30m_landcover_predictions/extracted/gplcc_gcam_45_2050.tif"
  )

r_gcam_45_2100 <- raster::raster(
    "/global_workspace/terry_sohl_sres_btu_30m_landcover_predictions/extracted/gplcc_gcam_45_2100.tif"
  )

s_2014 <- pts_to_landcover_metrics(
       pts=s,
       grid_units=units,
       r=r_gcam_45_2014,
       composition_statistics=composition_statistics,
       force_valid_habitat_val=c(19,26) # only use 'grass' for our config. metric
     )

s_2050 <- pts_to_landcover_metrics(
       pts=s,
       grid_units=units,
       r=r_gcam_45_2050,
       composition_statistics=composition_statistics,
       force_valid_habitat_val=c(19,26) # only use 'grass' for our config. metric
     )


s_2100 <- pts_to_landcover_metrics(
       pts=s,
       grid_units=units,
       r=r_gcam_45_2100,
       composition_statistics=composition_statistics,
       force_valid_habitat_val=c(19,26) # only use 'grass' for our config. metric
     )

grass_densities <- rbind(
    data.frame(year=2016, grass_ar=s_2014$grass_ar),
    data.frame(year=2050, grass_ar=s_2050$grass_ar),
    data.frame(year=2100, grass_ar=s_2100$grass_ar)
  )

# dev.new()
# ggplot2_multivariate_densities(grass_densities, var='grass_ar', xlab="Grassland [Total Area] (km2)")

shrub_densities <- rbind(
    data.frame(year=2016, shrub_ar=s_2014$shrub_ar),
    data.frame(year=2050, shrub_ar=s_2050$shrub_ar),
    data.frame(year=2100, shrub_ar=s_2100$shrub_ar)
  )

# dev.new()
# ggplot2_multivariate_densities(shrub_densities, var='shrub_ar', xlab="Shrubland [Total Area] (km2)")

# wetland_densities <- rbind(
#     data.frame(year=2016, wetland_ar=s_2014$wetland_ar),
#     data.frame(year=2050, wetland_ar=s_2050$wetland_ar),
#     data.frame(year=2100, wetland_ar=s_2100$wetland_ar)
#   )

# dev.new()
# ggplot2_multivariate_densities(wetland_densities, var='wetland_ar', xlab="Wetland [Total Area] (km2)")

patch_densities <- rbind(
    data.frame(year=2016, pat_ct=s_2014$pat_ct),
    data.frame(year=2050, pat_ct=s_2050$pat_ct),
    data.frame(year=2100, pat_ct=s_2100$pat_ct)
  )

# dev.new()
# ggplot2_multivariate_densities(patch_densities, var='pat_ct', xlab="Patchiness [Patch Count] (km2)")

# Re-build our model selection criterion
MOD_SEL_THRESHOLD <- max(which(model_selection_table@Full$delta < AIC_SUBSTANTIAL_THRESHOLD))
# select the top models (by array position) that satisfy our threshold
MOD_SEL_THRESHOLD <- as.numeric(model_selection_table@Full$model)[1:MOD_SEL_THRESHOLD]

#
# Scale our datasets so they are consistent for predict()
#

s_2014@data <- as.data.frame(
    scale(s_2014@data[,names(s_2014)], attr(m_scale, "scaled:center")[names(s_2014)], attr(m_scale, "scaled:scale")[names(s_2014)])
  )

#s_2014@data <- cbind(s_2014@data, s@data[,c('mat','map')])
s_2014$effort <- effort

predicted_2014 <- par_unmarked_predict(
  unmarked_models=unmarked_models[MOD_SEL_THRESHOLD],
  predict_df=s_2014@data,
  type="lambda",
  weights=aic_weights
)

if(NORMALIZE){
  predicted_2014 <- mean_normalization(
      predicted_2014, 
      c(m_pois_intercept_n, m_negbin_intercept_n)
    )
}

s_2050@data <- as.data.frame(
    scale(s_2050@data[,names(s_2050)], attr(m_scale, "scaled:center")[names(s_2050)], attr(m_scale, "scaled:scale")[names(s_2050)])
  )

#s_2050@data <- cbind(s_2050@data, s@data[,c('mat','map')])
s_2050$effort <- effort

predicted_2050 <- par_unmarked_predict(
  unmarked_models=unmarked_models[MOD_SEL_THRESHOLD],
  predict_df=s_2050@data,
  type="lambda",
  weights=aic_weights
)

if(NORMALIZE){
  predicted_2050 <- mean_normalization(
      predicted_2050,
      c(m_pois_intercept_n, m_negbin_intercept_n)
    )
}

s_2100@data <-  as.data.frame(
    scale(s_2100@data[,names(s_2100)], attr(m_scale, "scaled:center")[names(s_2100)], attr(m_scale, "scaled:scale")[names(s_2100)])
  )

#s_2100@data <- cbind(s_2100@data, s@data[,c('mat','map')])
s_2100$effort <- effort

predicted_2100 <- par_unmarked_predict(
  unmarked_models=unmarked_models[MOD_SEL_THRESHOLD],
  predict_df=s_2100@data,
  type="lambda",
  weights=aic_weights
)

if(NORMALIZE){
  predicted_2100 <- mean_normalization(
      predicted_2100, 
      c(m_pois_intercept_n, m_negbin_intercept_n)
    )
}

save(
    compress=T,
    list=(ls(pattern="[a-z]")),
    file=r_data_file
  )
