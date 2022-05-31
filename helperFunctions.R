#' Add Admin 1 and Admin 2 metadata
#'
#' @param data DHS survey data, processed via SUMMER::get_births
#' @param gps DHS gps data, containing GPS coordinates of sampling clusters
#' @param gadm1 GADM shapefile, admin1
#' @param gadm2 GADM shapefile, admin2

add_adm1_adm2_meta <- function (data, gps, gadm1 = NULL, gadm2 = NULL) {
  
  # detect points in the DHS GPS file with mis-specified coordinates and remove them if any
  wrong.points <- which(gps$LATNUM == 0.0 & gps$LONGNUM == 0.0)
  if (length(wrong.points) > 0) {
    message("Removing ", length(wrong.points), " wrong GPS points: (Longitude, Latitude) = (0, 0)")
  }
  
  # remove wrong points in the data if any
  data <- data[!(data$v001 %in% gps$DHSCLUST[wrong.points]),]
  
  # add latitude and longitude to data
  data <- merge(data, gps[, c("LATNUM", "LONGNUM", "DHSCLUST")], by.x = "v001", by.y = "DHSCLUST")
  data <- data[!is.na(data$LATNUM) & !is.na(data$LONGNUM), ]
  
  # retrieve GPS coordinates where data is sampled and convert to sp object
  points.frame <- as.data.frame(data[,c("LONGNUM", "LATNUM")]) 
  points.frame <- SpatialPoints(points.frame)
  
  # admin 2
  if (!is.null(gadm2)) {
    
    poly.over.adm2 <- SpatialPolygons(gadm2@polygons)
    proj4string(points.frame) <- proj4string(poly.over.adm2) <- 
      proj4string(gadm2)  <- 
      proj4string(gadm1)  
    admin2.key <- over(points.frame, poly.over.adm2)
    miss.frame.adm2 <- unique(points.frame@coords[which(is.na(admin2.key)),])
    
    if(dim(miss.frame.adm2)[1] != 0){
      miss.gadm2 <- dist2Line( miss.frame.adm2, poly.over.adm2)
      
      for(i in 1:dim(miss.gadm2)[1]){
        long.ids <- which(points.frame@coords[,c("LONGNUM")] %in% miss.frame.adm2[i,1])
        lat.ids <- which(points.frame@coords[,c("LATNUM")] %in% miss.frame.adm2[i,2])
        ids <- intersect(long.ids, lat.ids)
        admin2.key[ids] <- rep(miss.gadm2[i, 'ID'], length(ids))
      }
    }
    
    data$admin2 <- admin2.key
    data$admin2.char <- paste0("admin2_", admin2.key)
    data$admin2.name <- as.character(gadm2@data$NAME_2)[admin2.key]
    
  } else {
    data$admin2 <- data$admin2.name <- NA
    message("There is no Admin2 polygon to assign points to.")
  }
  
  # admin 1
  if (!is.null(gadm1)) {
    
    poly.over.adm1 <- SpatialPolygons(gadm1@polygons)
    proj4string(points.frame) <- proj4string(poly.over.adm1) <- 
      proj4string(gadm1) 
    admin1.key <- over(points.frame, poly.over.adm1)
    miss.frame.adm1 <- unique(points.frame@coords[which(is.na(admin1.key)),])
    
    if(dim(miss.frame.adm1)[1] != 0){
      miss.gadm1 <- dist2Line( miss.frame.adm1, poly.over.adm1)
      
      for(i in 1:dim(miss.gadm1)[1]){
        long.ids <- which(points.frame@coords[,c("LONGNUM")] %in% miss.frame.adm1[i,1])
        lat.ids <- which(points.frame@coords[,c("LATNUM")] %in% miss.frame.adm1[i,2])
        ids <- intersect(long.ids, lat.ids)
        admin1.key[ids] <- rep(miss.gadm1[i, 'ID'], length(ids))
      }
    }
    
    data$admin1 <- admin1.key
    data$admin1.char <- paste0("admin1_", admin1.key)
    data$admin1.name <- as.character(gadm1@data$NAME_1)[admin1.key]
    
  } else {
    data$admin2 <- data$admin2.name <- NA
    message("There is no Admin1 polygon to assign points to.")
  }  
  
  return(data)
}


