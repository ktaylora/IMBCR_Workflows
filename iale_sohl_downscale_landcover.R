load(commandArgs(trailingOnly=T))

transects <- OpenIMBCR:::scrub_imbcr_df(
    OpenIMBCR:::imbcrTableToShapefile(
        "/global_workspace/imbcr_number_crunching/results/RawData_PLJV_IMBCR_20161201.csv"
      ),
    four_letter_code = toupper(argv[1])
  )
  
s$transect  <- as.character(unique(transects$transectnum))

#
# load our vegetation data and summarize by transect
#

veg_data <- read.csv(
  list.files(
      "..", 
      pattern="raw_veg.*.csv", 
      recursive=T, 
      full.names=T
    ))

grass_col_names <- colnames(veg_data)[ grep(colnames(veg_data), pattern="gc_.*._grass$|gc_.*._herb$") ]
shrub_col_names <- colnames(veg_data)[ grep(colnames(veg_data), pattern="gc_woody|shrub_cover") ]
 bare_col_names <- colnames(veg_data)[ grep(colnames(veg_data), pattern="gc_bare_litter") ]

veg_data <- do.call(rbind, lapply(
  X=as.vector(unique(s$transect)),
  FUN=function(x){
    focal <- as.character(as.vector(veg_data$TransectNum)) %in% x
    if(sum(focal)>0){
      grass_cover <- as.vector(apply(veg_data[focal,grass_col_names], MARGIN=1, FUN=sum))
      shrub_cover <- as.vector(apply(veg_data[focal,shrub_col_names], MARGIN=1, FUN=sum))
      bare_ground <- as.vector(veg_data[focal,bare_col_names])
      # this will be coded with counts where grass[1], shrub[2], or bare[3] were
      # the dominant types at each station
      t <- table(apply(rbind(grass_cover, shrub_cover, bare_ground), MARGIN=2, which.max))
      names <- gsub(names(t), pattern="1", replacement="grass")
      names <- gsub(names, pattern="2", replacement="shrub")
      names <- gsub(names, pattern="3", replacement="bare")
      t <- data.frame(matrix(t, nrow=1))
      colnames(t) <- names
      ret <- data.frame(grass=0, shrub=0, bare=0)
      ret[ , which(colnames(t) %in% colnames(ret))] <- 
        t[,which(colnames(t) %in% colnames(ret))]
      return(ret)
    } else {
      return(data.frame(grass=0, shrub=0, bare=0))
    }
  }
))

# re-scale so we add-up to 100%
#veg_data <- round((veg_data / as.vector(apply(veg_data, MARGIN=1, FUN=sum))), 2) * 100   
#  veg_data[is.na(veg_data)] <- 0

# convert percent-cover to total area 
#veg_data <- ((veg_data/100)*(rgeos::gArea(units[1,])))*(10^-6)

climate_data <- OpenIMBCR:::spatial_join(s, units)
climate_data <- data.frame(scale(climate_data@data[,c('map','mat')]))

# correlation tests
cor(na.omit(cbind(s$grass_ar,scale(veg_data$grass))))
cor(na.omit(cbind(s$shrub_ar,scale(veg_data$shrub))))

# fit a summary model
fitting <- data.frame(
    grass=veg_data$grass, 
    shrub=veg_data$shrub,
    bare=veg_data$bare,
    grass_ar=s$grass_ar, 
    shrub_ar=s$shrub_ar,
    pat_ct=s$pat_ct,
    mat=climate_data$mat,
    map=climate_data$map
  )
  
m_grass_nb_all_covs <- (MASS::glm.nb(
    grass~poly(grass_ar,1,raw=T)+poly(mat,1,raw=T)+poly(map,1,raw=T)+offset(log(effort)), 
    data=fitting,
  ))
  
range(as.vector(round(predict(m_grass_nb, newdata=fitting, type="response")/16, 2)))

m_grass_nb <- (MASS::glm.nb(
    grass~poly(grass_ar,2,raw=T)+offset(log(effort)), 
    data=fitting,
  ))

grass_cover_upscaled <- round(as.vector(round(predict(
    m_grass_nb, 
    newdata=fitting, 
    type="response"))
  )/16, 2) 
  
grass_cover_upscaled <- grass_cover_upscaled * (rgeos::gArea(units[1,]))*(10^-6)
  
#range(as.vector(round(predict(m_grass_pois, newdata=fitting, type="response")/16, 2)))
    
m_shrub_nb_all_covs <- (MASS::glm.nb(
    shrub~poly(shrub_ar,2,raw=T)+poly(mat,1,raw=T)+poly(map,1,raw=T)+offset(log(effort)), 
    data=fitting,
  ))

m_shrub_nb <- (MASS::glm.nb(
    shrub~poly(shrub_ar,1,raw=T)+offset(log(effort)), 
    data=fitting,
  ))
  
shrub_cover_upscaled <- round(as.vector(round(predict(
    m_shrub_nb_all_covs, 
    newdata=fitting, 
    type="response"))
  )/16, 2)
  
shrub_cover_upscaled <- shrub_cover_upscaled * (rgeos::gArea(units[1,]))*(10^-6)
