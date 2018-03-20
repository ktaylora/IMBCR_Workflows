require(OpenIMBCR)

#
# Local Functions
#
write_route_detections <- function(name=NULL, s=NULL, write=T){
  if(write){
    archive_name <- paste(name,"_transect_detections", sep="")
    rgdal::writeOGR(
      s, 
      dsn=".", 
      layer=archive_name, 
      driver="ESRI Shapefile",
      overwrite=T
    )
    unlink(paste(archive_name,".zip",sep=""))
    system(
      paste(
        "7za a", 
        paste(archive_name,".zip",sep=""), 
        paste(archive_name,".*",sep=""),
        sep=" "
      )
    )
    rm <- list.files(".",pattern=archive_name)
      rm <- rm[!grepl(rm, pattern="zip$")]
    unlink(rm)
  } else {
    return(s)
  }
}

#
# MAIN
#

argv <- commandArgs(trailingOnly=T)

if(length(argv) < 2){
  warning("using default (2016) IMBCR CSV")
  imbcr_csv <- "/global_workspace/imbcr_number_crunching/results/RawData_PLJV_IMBCR_20161201.csv"
} else {
  imbcr_csv <- argv[2]
}

s <- OpenIMBCR:::scrub_imbcr_df(
    OpenIMBCR:::imbcrTableToShapefile(
        imbcr_csv
      ),
    four_letter_code = toupper(argv[1])
  )

detections <- OpenIMBCR:::calc_dist_bins(s)
  s <- OpenIMBCR:::calc_route_centroids(s)

s@data <- data.frame(detections=rowSums(detections$y))

write_route_detections(name=toupper(argv[1]), s=s, write=T)
