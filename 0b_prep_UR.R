
# compute urban fractions by year
# Based on code here: https://github.com/wu-thomas/SUMMER-DHS/tree/main/Rcode/UR%20Pipeline

library(SUMMER)
library(classInt)
library(RColorBrewer)
library(dplyr)
library(tidyr)
library(rgdal)
library(scales)
library(INLA)
library(survey)
library(ggplot2)
library(maptools)
library(gridExtra)
library(mgcv)
library(caret)
library(geosphere)
library(rgeos)
library(sqldf)
library(raster)
library(spdep)
library(haven)
library(labelled)
library(data.table)
library(sqldf)
library(sp)
library(gstat)
library(stringdist)
library(openxlsx)

# settings
country <- "Liberia"
country_code <- "LBR"
frame.year <- 2008
beg.year <- 2000
end.year <- 2020


# Shapefiles --------------------------------------------------------------

gadm.abbrev <- country_code

poly.path <- paste0("Data/", country, "/shapeFiles_gadm") # specify the folder of the country shape files
poly.layer.adm0 <- paste('gadm36', gadm.abbrev,
                         '0', sep = "_") # specify the name of the national shape file
poly.layer.adm1 <- paste('gadm36', gadm.abbrev,
                         '1', sep = "_") # specify the name of the admin1 shape file
poly.layer.adm2 <- paste('gadm36', gadm.abbrev,
                         '2', sep = "_") # specify the name of the admin2 shape file

poly.adm0 <- readOGR(dsn = poly.path,
                     layer = as.character(poly.layer.adm0)) # load the national shape file
poly.adm1 <- readOGR(dsn = poly.path,
                     layer = as.character(poly.layer.adm1)) # load the shape file of admin-1 regions
poly.adm2 <- readOGR(dsn = poly.path,
                     layer = as.character(poly.layer.adm2)) # load the shape file of admin-2 regions

proj4string(poly.adm0) <- proj4string(poly.adm1) <- proj4string(poly.adm2)

# create the adjacency matrix for admin1 regions.
admin1.mat <- poly2nb(SpatialPolygons(poly.adm1@polygons))
admin1.mat <- nb2mat(admin1.mat, zero.policy = TRUE)
colnames(admin1.mat) <- rownames(admin1.mat) <- paste0("admin1_", 1:dim(admin1.mat)[1])
admin1.names <- data.frame(GADM = poly.adm1@data$NAME_1,
                           Internal = rownames(admin1.mat))

# create the adjacency matrix for admin2 regions.
admin2.mat <- poly2nb(SpatialPolygons(poly.adm2@polygons))
admin2.mat <- nb2mat(admin2.mat, zero.policy = TRUE)
colnames(admin2.mat) <- rownames(admin2.mat) <- paste0("admin2_", 1:dim(admin2.mat)[1])
admin2.names <- data.frame(GADM = poly.adm2@data$NAME_2,
                           Internal = rownames(admin2.mat))


# Load worldpop -----------------------------------------------------------

### automated downloading, if not working, try manually download
pop.abbrev <- tolower(gadm.abbrev)
census.year <- frame.year

file <- paste0( pop.abbrev,'_ppp_',census.year,'_1km_Aggregated_UNadj.tif')
pop_dir <- paste0("Data/", country, "/Population")

if(!file.exists(paste0(pop_dir, "/", file))){
  
  url <- paste0("https://data.worldpop.org/GIS/Population/Global_2000_2020_1km_UNadj/", 
                census.year, "/", toupper(pop.abbrev),"/",      
                pop.abbrev,'_ppp_',census.year,'_1km_Aggregated_UNadj.tif')
  
  download.file(url, file, method = "libcurl",mode="wb")
}

# UNadjusted population counts
worldpop <- raster(paste0(pop_dir, "/",
                          tolower(gadm.abbrev),
                          '_ppp_',frame.year,
                          '_1km_Aggregated_UNadj.tif',sep=''))


# Load cluster data -------------------------------------------------------

# load DHS cluster information
mod.dat <- readRDS(paste0("Data/", country, "/nmr_data_prepped.rds"))

# correct urban rural numeric to character
if (is.numeric(mod.dat$urban)) {
  mod.dat$urban <- ifelse(mod.dat$urban == 1, "urban", "rural")
}

# clusters
mod.dat$survey <- mod.dat$survey_year
cluster_list<-mod.dat[!duplicated(mod.dat[c('cluster','survey',
                                            'LONGNUM','LATNUM')]),]

# keep only the clusters from DHS surveys using the recent census
# cluster_list<-cluster_list[cluster_list$survey %in% survey_year,]


# Correct urban cluster ---------------------------------------------------

# The function below corrects urban clusters that are misclassified to be rural due to jittering. Jittered location is not the exact location 
# but a randomly shifted location for the purpose of confidentiality. This may results in some urban clusters to be jittered to some rural areas, 
# which is unexpected to our classification algorithm and therefore correcting the misclassified clusters is of interest. 

# The codes below fulfills the process above by assigning the possibly misclassified urban clusters to the nearest most densely populated areas.
# It's generally true that the urban areas tend to have a higher population density so this process can alleviate the side effect
# of jittering.
constr_prior <- function(obs,jitter_r,prob_r,poly_admin,pop_ras){
  
  # for the cluster, find its coordinates and admin1 area
  admin1_index<-obs$admin2
  sp_xy<-SpatialPoints(as.data.frame(obs)[,c('LONGNUM','LATNUM')],
                       proj4string = CRS(proj4string(poly_admin)))
  #pt<-as.data.frame(obs)[,c('LONGNUM','LATNUM')]
  #colnames(pt)<-c('x','y')
  
  # generate gitter 
  #jitter_r<-2000
  cluster_buffer<-buffer(sp_xy, width=jitter_r)
  
  # extract pixels within the buffer
  temp_pop_cir<-mask(crop(pop_ras,cluster_buffer),
                     cluster_buffer)
  
  # put admin area restriction
  admin1_poly<-poly_admin[admin1_index,]
  temp_pop_admin<-mask(crop(temp_pop_cir,admin1_poly),
                       admin1_poly)
  
  
  # check whether need to adjust for constraint 
  cir_val<-values(temp_pop_cir)
  admin_val<-values(temp_pop_admin)
  
  admin_adj<-length(which(!is.na(cir_val)))!=length(which(!is.na(admin_val)))
  
  if(admin_adj){
    #normc<-admin1_normc(pt,jitter_r,admin1_poly,ntrial=1000)
    normc<-1
  }else{  normc<-1}
  
  
  
  ## prepare sample frame
  
  temp_val<-values(temp_pop_admin)
  pop_index<-which(!is.na(temp_val))
  
  
  temp_frame<-as.data.frame(coordinates(temp_pop_admin))
  pixel_candidate<-temp_frame[pop_index,]
  pixel_candidate$pop_den<-temp_val[pop_index]
  
  pixel_candidate$center_x<-obs$LONGNUM
  pixel_candidate$center_y<-obs$LATNUM
  pixel_candidate$dist<-diag(distm(pixel_candidate[,c('x','y')], 
                                   pixel_candidate[,c('center_x','center_y')]))
  
  pixel_candidate$unn_w<-pixel_candidate$pop_den*
    1/(2*pi * 2 * pixel_candidate$dist)*normc*prob_r
  
  pixel_candidate$normc<-normc
  return(pixel_candidate[,c("x","y",'normc','unn_w')])
  #return(pixel_candidate)
  
}

# only correct urban clusters 
# check if the strata are named 'urban' and 'rural'
urban_clus<-cluster_list[cluster_list$urban=='urban',]
rural_clus<-cluster_list[cluster_list$urban=='rural',]

urban_clus$x_adj<-NA
urban_clus$y_adj<-NA
rural_clus$x_adj<-rural_clus$LONGNUM
rural_clus$y_adj<-rural_clus$LATNUM

# points have to stay within the same admin2 region 
for( i in 1:dim(urban_clus)[1]){
  
  print(i)
  temp_frame<-constr_prior(urban_clus[i,],2000,1,poly.adm2,worldpop)
  p_mode = sqldf("SELECT * FROM temp_frame GROUP BY x,y ORDER BY SUM(unn_w) DESC LIMIT 1")
  urban_clus[i,]$x_adj<-p_mode$x
  urban_clus[i,]$y_adj<-p_mode$y
  
}


prep_dat<-rbind(urban_clus,rural_clus)
#xy <- as.matrix(prep_dat[c('x_adj','y_adj')])

# save(prep_dat, file = paste0("Data/", country, "/prep_dat.rda"))


# Adding covariates -------------------------------------------------------

crc_dat<-prep_dat

# set up corrected xy for clusters
xy_crc <- as.matrix(crc_dat[c('x_adj','y_adj')])
crc_dat$x<-crc_dat$x_adj # x_adj and y_adj: corrected coordinates
crc_dat$y<-crc_dat$y_adj

# extract covariates
crc_dat$pop_den<-raster::extract(worldpop,xy_crc)

# only retain part of the columns to reduce redundancy
col_select<-c('cluster','urban','admin1','admin2',
              'admin1.name','admin2.name',
              'admin1.char','admin2.char',
              'survey','pop_den','x','y')

crc_dat_final<-crc_dat[,col_select]
#save(crc_dat_final,file='prepared_dat/crc_dat.rda')


# Prepare data without urban correction -----------------------------------

uncrc_dat<-prep_dat

# set up uncorrected xy for clusters
xy_uncrc <- as.matrix(uncrc_dat[c('LONGNUM','LATNUM')])
uncrc_dat$x<-uncrc_dat$LONGNUM # x_adj and y_adj: corrected coordinates
uncrc_dat$y<-uncrc_dat$LATNUM

# extract covariates
uncrc_dat$pop_den<-raster::extract(worldpop,xy_uncrc)

# keep columns
col_select<-c('cluster','urban','admin1','admin2',
              'admin1.name','admin2.name',
              'admin1.char','admin2.char',
              'survey','pop_den','x','y')

uncrc_dat_final<-uncrc_dat[,col_select]
# save(uncrc_dat_final,file='prepared_dat/uncrc_dat.rda')


# Prepare U1 population ---------------------------------------------------

### download U1 population, if not working, try manually download
#! downloading might take a long time, especially for big countries

pop.year <- beg.year:end.year

options(timeout = 1000) # adjust this time, should be longer than each download
for(year in pop.year){ #pop.year <- beg.year:end.year   ## population surface year
  print(year)
  for(age in c(0)){
    for(sex in c("f", "m")){
      file <- paste0("Data/", country, "/Population/", pop.abbrev,'_', sex, '_', age, '_', year,'.tif')
      if(!file.exists(file)){
        url <- paste0("https://data.worldpop.org/GIS/AgeSex_structures/Global_2000_2020/", 
                      year, "/", toupper(pop.abbrev), "/", pop.abbrev, "_", 
                      sex, "_", age, "_", year, ".tif")
        download.file(url, file, method = "libcurl",mode="wb")
      }
    }
  }
}

years <- c(beg.year:end.year)

# aggregate the under-1 population spatial surface to resolution of 1km*1km
for ( t in 1:length(years)){
  
  year <- years[t]
  print(year)
  
  ## first sum up four rasters female 0-1, male 0-1
  f_0_name = paste0("Data/", country, "/Population/", pop.abbrev,'_f_0_',year,'.tif')
  m_0_name = paste0("Data/", country, "/Population/", pop.abbrev,'_m_0_',year,'.tif')
  
  pop_f_0<-raster(f_0_name)
  pop_m_0<-raster(m_0_name)
  
  proj4string(pop_f_0) <- proj4string(pop_m_0) <- proj4string(poly.adm1)
  
  pop_u1<- pop_f_0+pop_m_0
  
  writeRaster(pop_u1, overwrite=TRUE,
              paste0("Data/", country, "/Population/",
                     pop.abbrev,'_u1_',year,'_100m','.tif'))
  
  pop_surf<-raster(paste0("Data/", country, "/Population/",
                          pop.abbrev,'_u1_',year,'_100m','.tif'))
  
  pop_u1_aggregate <- aggregate(pop_surf, fact=10, sum)
  
  pop_grid<-as.data.frame(coordinates(worldpop))
  colnames(pop_grid)<-c('x','y')
  
  pop_grid$u1_pop <- raster::extract(pop_u1_aggregate, 
                                     pop_grid[c('x','y')])
  
  u1_pop<-worldpop
  values(u1_pop)<-pop_grid$u1_pop 
  
  writeRaster(u1_pop, overwrite=TRUE,
              paste0("Data/", country, "/Population/",
                     pop.abbrev,'_u1_',year,'_1km.tif'))
  
}


# Prepare national grid ---------------------------------------------------

## set up grid
urb_dat<-as.data.frame(coordinates(worldpop))
colnames(urb_dat)<-c('x','y')

## add population density
urb_dat$pop_den<-raster::extract(worldpop,urb_dat[c('x','y')])

## add admin1 region for pixel
points.frame <- as.data.frame(urb_dat[,c("x", "y")])
points.frame <- SpatialPoints(points.frame)
poly.over.adm1 <- SpatialPolygons(poly.adm1@polygons)
admin1.key <- over(points.frame, poly.over.adm1)
urb_dat$admin1<-admin1.key
urb_dat$admin1.char <- paste0("admin1_", admin1.key)

# save(urb_dat,file='prepared_dat/natl_grid.rda')


# Sampling frame urban proportion at admin1 -------------------------------

# read the xlsx file containing urban population fraction at admin1 level.
frame <- read.csv(paste0("Data/", country, "/frame_urb_prop.csv"))
names(frame) <- c("V1", "frac")

# greedy algorithm to match admin names 
adm1.ref <- expand.grid(tolower(frame[, 1]),
                        tolower(admin1.names$GADM)) # Distance matrix in long form
names(adm1.ref) <- c("frame_name","gadm_name")
### string distance,  jw=jaro winkler distance, try 'dl' if not working
adm1.ref$dist <- stringdist(adm1.ref$frame_name,
                            adm1.ref$gadm_name, method="jw") 

greedyAssign <- function(a,b,d){
  x <- numeric(length(a)) # assgn variable: 0 for unassigned but assignable, 
  # 1 for already assigned, -1 for unassigned and unassignable
  while(any(x==0)){
    min_d <- min(d[x==0]) # identify closest pair, arbitrarily selecting 1st if multiple pairs
    a_sel <- a[d==min_d & x==0][1] 
    b_sel <- b[d==min_d & a == a_sel & x==0][1] 
    x[a==a_sel & b == b_sel] <- 1
    x[x==0 & (a==a_sel|b==b_sel)] <- -1
  }
  cbind(a=a[x==1],b=b[x==1],d=d[x==1])
}

match_order<-data.frame(greedyAssign(adm1.ref$frame_name,
                                     adm1.ref$gadm_name,
                                     adm1.ref$dist))

# create reference table 
ref.tab <- admin1.names
ref.tab$matched_name <- frame$V1[match_order$a] ### check!!!
ref.tab$urb_frac <- frame$frac[match_order$a] 


# Admin 1 threshold -------------------------------------------------------

# index the grid
urb_dat$index <- c(1:nrow(urb_dat))
adm1_dat <- split( urb_dat , f = urb_dat$admin1 )

# This function computes the urban population threshold for a given admin1 area.
# This is done by keep counting the urban locations until the urban population
# fraction in the reference table is reached.
thresh_urb<-function(adm_grid,ref_tab){
  
  # sort grid population
  vals <- adm_grid$pop_den
  vals[is.na(vals)] <- 0
  sort.obj <- sort.int(vals, decreasing = TRUE, index.return = TRUE, method = 'shell')
  svals <- sort.obj$x
  svals.int <- sort.obj$ix
  
  # extract cutoff proportion based on admin1
  adm.idx <- adm_grid$admin1.char[1]
  cutoff <- ref_tab[ref_tab$Internal==adm.idx,]$urb_frac
  
  # determine population threshold and urban rural
  csvals <- cumsum(svals)/sum(svals)
  is.urb <- csvals <= cutoff
  org.isurb <- is.urb[invPerm(svals.int)]
  threshold <- min(vals[org.isurb == 1]) #cutoff
  
  # prepare return object (grid with urban/rural)
  adm_grid$threshold <- threshold
  adm_grid$urban <- as.numeric(org.isurb)
  #adm_grid[is.na(adm_grid$pop_den),]$urban<-NA
  
  return(adm_grid)
  
}

urb_list<-lapply(adm1_dat, FUN=thresh_urb,ref_tab=ref.tab)

urb_class <- do.call("rbind", urb_list)

urb_grid <- urb_dat
urb_grid$urb_ind <-NA
urb_grid[urb_class$index,]$urb_ind <- urb_class$urban

urb_surf<-worldpop
values(urb_surf)<-urb_grid$urb_ind

## save reference table along with calculated threshold 
thresh_ref <- urb_class[!duplicated(urb_class[,c('admin1')]),]
ref.tab$threshold <- thresh_ref$threshold # check whether the thresholds are sensible (shouldn't be NA or all 0)
write.xlsx(ref.tab,
           file=paste0('Data/', country, '/reference_table.xlsx'),
           row.names = FALSE)


# Check classification accuracy based on clusters -------------------------

### remove rows with missing covariates, could also build model with missing data
crc_dat<-crc_dat_final[complete.cases(crc_dat_final), ]
uncrc_dat<-uncrc_dat_final[complete.cases(uncrc_dat_final), ]

xy_crc <- as.matrix(crc_dat[c('x','y')])
xy_uncrc <- as.matrix(uncrc_dat[c('x','y')])

# extract the urban/rural prediction
crc_dat$urb_pred<-raster::extract(urb_surf,xy_crc)
uncrc_dat$urb_pred<-raster::extract(urb_surf,xy_uncrc)
pred_crc <- factor( ifelse(crc_dat$urb_pred ==1 ,"urban","rural" ))
pred_crc  <- relevel(pred_crc, "urban") # make sure levels are same 
pred_uncrc <- factor( ifelse(uncrc_dat$urb_pred ==1 ,"urban","rural" ))
pred_uncrc  <- relevel(pred_uncrc, "urban") # make sure levels are same 


# Load function for urban fraction ----------------------------------------

get_subnatl_frac<-function(adm.names,adm.idx,wp,poly_file,wp_adm=NULL,
                           urb_vec){
  
  poly_file <- spTransform(poly_file, wp@crs)
  
  if(is.null(wp_adm))
    wp_adm <- lapply(1:nrow(poly_file), function(x) {
      list(state_id = x, state_raster = mask(crop(wp,poly_file[x,]), poly_file[x,]))
    })
  
  pred_surf <- wp
  values(pred_surf)<-urb_vec
  
  urb_adm <- lapply(1:nrow(poly_file), function(x) {
    list(state_id = x, state_raster = mask(crop(pred_surf,poly_file[x,]), poly_file[x,]))
  })
  
  frac_vec<-vector()
  
  for(j in 1:length(adm.names)){
    urb_j<-urb_adm[[j]]
    wp_j<-wp_adm[[j]]
    
    val_urb_j<-values(urb_j$state_raster)
    val_wp_j<-values(wp_j$state_raster)
    
    frac_vec[j]<-sum(val_urb_j*val_wp_j,na.rm=TRUE)/
      sum(val_wp_j,na.rm=TRUE)
  }
  
  subnatl_frac<-data.frame(adm_name=adm.names,adm_idx=adm.idx,urb_frac=frac_vec)
  return(subnatl_frac)
}


# National U1 urban proportion --------------------------------------------

years <- c(beg.year:end.year)
natl.u1.urb <- vector()

for ( t in 1:length(years)){
  
  year <- years[t]
  print(year)
  
  # load U1 population at year t
  u1_pop<-raster(paste0("Data/", country, "/Population/",
                        country.abbrev,'_u1_',year,'_1km.tif'))
  
  # national urban fraction for U1 population at year t
  u1_natl <- sum(urb_grid$urb_ind*values(u1_pop),na.rm=TRUE)/
    sum(values(u1_pop),na.rm=TRUE)
  natl.u1.urb[t] <- u1_natl
  
}

natl.urb.weights <- data.frame(years= years, urban=natl.u1.urb)

if (!dir.exists(paste0('Results/', country))) {
  dir.create(paste0('Results/', country))
}
saveRDS(natl.urb.weights,
        paste0('Results/', country, '/natl_urban_weights.rds'))


# Subnational U1 urban proportion -----------------------------------------

years <- c(beg.year:end.year)
adm1.weight.frame <- data.frame()
adm2.weight.frame <- data.frame()

for ( t in 1:length(years)){
  
  year <- years[t]
  print(year)
  
  # load U1 population at year t
  u1_pop<-raster(paste0("Data/", country, "/Population/",
                        country.abbrev,'_u1_',year,'_1km.tif'))
  
  # admin1 urban fraction for U1 population at year t
  u1_urb_admin1<-get_subnatl_frac(adm.names = admin1.names$GADM,
                                  adm.idx = admin1.names$Internal,
                                  wp=u1_pop,
                                  poly_file = poly.adm1,
                                  wp_adm = NULL,
                                  urb_vec = urb_grid$urb_ind)
  
  u1_urb_admin1$years <- year
  adm1.weight.frame <- rbind(adm1.weight.frame,u1_urb_admin1)
  
  # admin2 urban fraction for U1 population at year t
  u1_urb_admin2<-get_subnatl_frac(adm.names = admin2.names$GADM,
                                  adm.idx = admin2.names$Internal,
                                  wp=u1_pop,
                                  poly_file = poly.adm2,
                                  wp_adm = NULL,
                                  urb_vec = urb_grid$urb_ind)
  
  u1_urb_admin2$years <- year
  adm2.weight.frame <- rbind(adm2.weight.frame,u1_urb_admin2)
  
  # save calculated urban fractions
  saveRDS(u1_urb_admin1,file=paste0('Results/', country, '/admin1_',
                                    year, '_urban_frac.rds'))
  saveRDS(u1_urb_admin2,file=paste0('Results/', country, '/admin2_',
                                    year, '_urban_frac.rds'))
  
}

# process admin 1 urban rural weights data frame
adm1.weight.frame <- adm1.weight.frame[,c('adm_idx','years','urb_frac')]
colnames(adm1.weight.frame) <- c('region','years','urban')
adm1.weight.frame$rural <- 1 - adm1.weight.frame$urban

# process admin 2 urban rural weights data frame
adm2.weight.frame <- adm2.weight.frame[,c('adm_idx','years','urb_frac')]
colnames(adm2.weight.frame) <- c('region','years','urban')
adm2.weight.frame$rural <- 1 - adm2.weight.frame$urban

# save weights frames
saveRDS(adm1.weight.frame,paste0('Results/', country, '/admin1_urban_weights.rds'))
saveRDS(adm2.weight.frame,paste0('Results/', country, '/admin2_urban_weights.rds'))



