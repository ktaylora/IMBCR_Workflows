require(rgdal)
require(ggplot2)
require(OpenIMBCR)

backscale_var <- function(var=NULL, df=NULL, m_scale=NULL){
  return(df[, var] * attr(m_scale, 'scaled:scale')[var] + attr(m_scale, 'scaled:center')[var] )
}

ggplot2_univar_density <- function(var=NULL, correction=1, ylab=NULL, xlab=NULL){
  xlim <- range( backscale_var(var, s@data, m_scale)/correction )

  xlim[1] <- xlim[1] - xlim[1]*0.1
  xlim[2] <- xlim[2] + xlim[2]*0.1

  rug <- data.frame(x=backscale_var(var, s@data, m_scale) / correction)
  median <- median(rug$x)
  
  ggplot(
    data=data.frame(
        mat=backscale_var(var, s@data, m_scale)
      ) 
      , aes(mat)
    ) + 
    scale_x_continuous(limits=xlim) +
    geom_density(fill='#edf8b1', alpha=0.6, color="#97a837") +
    geom_rug(aes(x), sides="b", color="#97a837", alpha=0.8, size=0.25, data=rug) +
    xlab(xlab) +
    ylab(ylab) +
    geom_vline(xintercept = median, color="#97a837", size=0.5) +
    scale_color_brewer(type = 'seq', palette="YlGnBu") +
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

extract_landcover_by_point <- function(s=NULL, units=NULL, r=NULL){
  area_statistics <-
    data.frame(
        field_name=c(
          'grass_ar',
          'shrub_ar',
          'wetland_ar',
          'crp_ar'
        ),
        src_raster_value=c(
          '176',
          'c(64,152)',
          '195',
          '233'
        )
    )
  configuration_statistics <- c(
    'pat_ct'
  )
  valid_habitat_values <- eval(parse(
    text=paste("c(",paste(area_statistics$src_raster_value[
      !grepl(area_statistics$field_name, pattern="rd_ar")
    ], collapse = ","), ")", sep="")
  ))
  # subset our units shapefile by overlapping points in s
  grid_units_over <- !is.na(sp::over(
      units, spTransform(s, sp::CRS(raster::projection(units))))[,1]
    ) 
  units <- units[grid_units_over,]
  # zero-out our attribute table
  units@data <- data.frame(id=1:nrow(units))
  # split for list comprehension and parallelization
  e_units <- split(units, 1:nrow(units)) 
  # do our extractions
  e_units <- OpenIMBCR:::extract_by(e_units, r)
  # sanity-check : do we have enough units to cover our features in s?
  if(nrow(s) != length(e_units)){
    stop("we are missing some overlapping USNG units for features in s")
  }
  # calculate our landscape metric (total area) across our input dataset
  done <- lapply(
    X=1:nrow(area_statistics),
    FUN=function(i){
      OpenIMBCR:::par_calc_stat(
        # using our 1 km unit raster extractions
        X=e_units,
        fun = OpenIMBCR:::calc_total_area,
        # using these PLJV landcover cell values in the reclassification
        from = eval(parse(text=as.character(area_statistics[i, 2])))
      )
    }
  )
  
  done <- as.data.frame(do.call(cbind, done))
  units@data <- done
  
  # calculate our landscape configuration metric (patch count
  units@data[, 5] <-
  OpenIMBCR::par_calc_stat(
      # using our using our un-buffered unit raster extractions
      e_units,
      # parse the focal landscape configuration metric
      fun = OpenIMBCR::calc_patch_count,
      # using these PLJV landcover cell values in the supplemental
      # reclassification
      from = valid_habitat_values
    )
  # settle on our column names
  colnames(units@data) <- 
    c(as.vector(area_statistics$field_name), configuration_statistics)
  # return table of units to user
  return(units)
}

#
# Main
#

r <- raster::raster(paste("/gis_data/Landcover/NASS/Raster/",
    "2016_nass_crp_test_merge.tif", sep=""
  ))
  
r_gcam_45_rcp45_2020 <- raster::raster("/global_workspace/terry_sohl_sres_btu_30m_landcover_predictions/extracted/gplcc_gcam_45_rcp45_2020.tif")

test <- extract_landcover_by_point(s, units, r_gcam_45_rcp45_2020)

ggplot2_univar_density(var="shrub_ar", xlab="Total Area Shrubland (km2)")
