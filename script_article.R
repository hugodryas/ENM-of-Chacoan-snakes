########################################################################################################################
# Script to replicate the article A not promising future for the snakes of the Gran Chaco, 
# impact of climate change in snakes' community of the Gran Chaco
# Most part of the script were gathered from a postgradute course in Biodiversity Modelling and Global Change by Profesor Mario R Moura

setwd("Set Working Directory")

# Install and load R packages needed to run the analysis:
needed_packages<-c("ade4", "agricolae", "ape", "BAT", "beepr", "betapart", "BuenColors", "cluster", "CommEcol", "data.table", "dismo",
                   "doParallel", "dplyr", "ENMTML", "FD", "flexclust", "foreach", "gdm", "geiger", "ggcorrplot", "ggplot2", "ggnewscale", 
                   "ggrepel", "ggsflabel", "ggspatial", "ggtree", "grDevices", "leafsync", "mapview", "metacom", "MSDM", "ncf",
                   "PerformanceAnalytics", "pgirmess", "phangorn", "plyr", "prettyGraphs", "psych", "raster", "RColorBrewer", 
                   "reshape2", "Rfast", "installr", "rgdal", "letsR", "rgeos", "rnaturalearth", "RStoolbox", "sf", "sp",
                   "SpatialPack", "spThin", "stars", "stringr", "svMisc", "tools", "usdm", "vegan", "viridis", "weights")
new.packages<-needed_packages[!(needed_packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages, repos = "http://cran.us.r-project.org")
lapply(needed_packages, require, character.only = TRUE)
raster::removeTmpFiles(0.00001) # remove temporary files
rm(needed_packages, new.packages)

# Install packages available outside the CRAN repository:
installr::install.Rtools(check = F)
devtools::install_github("yutannihilation/ggsflabel", force=T)
devtools::install_github("ropenscilabs/rnaturalearthdata", force=T)
devtools::install_github("andrefaa/ENMTML", force=T) 
devtools::install_github("sjevelazco/MSDM") 
devtools::install_github('taiyun/corrplot', force=T)

# Create a basic directory structure to receive the outputs of this script:
dir() # visualize the working directory
dir.create("Project Folder", showWarnings = F)
dir.create("Datasets", showWarnings = F)
dir.create("Shapefiles", showWarnings = F) 
dir.create("Project Folder/Predictors Folder", showWarnings = F)
dir.create("Project Folder/Projection Folder", showWarnings = F)

# Step 1: BUILD ECOLOGICAL NICHE MODELS (ENM) FOR THE TARGET SPECIES
########################################################################################################################
# Paso 1: CONSTRUIR MODELOS DE NICHO ECOLÓGICO (ENM) PARA LAS ESPECIES
rm(list=ls()); gc()

# Load species occurrence data:
CompiledDatabase<-data.table::fread("Project Folder2/Occurences Folder/CleannedOcc.txt", stringsAsFactors=TRUE, encoding="UTF-8")
# Run the complete dataset:
a<-Sys.time()
ENMTML(
  pred_dir="Project Folder2/Predictors Folder", # path to predictors
  proj_dir = "Project Folder2/Projection Folder",  # path to projections
  result_dir = "Project Folder2/Results",  # where results will be saved
  occ_file="Project Folder2/Occurences Folder/CleannedOcc.txt", # path to species occurrence data
  sp="Species", # name of the column holding species identity
  x="Longitude", # name of column containing longitude in decimal degrees
  y="Latitude", # name of column containing latitude in decimal degrees
  min_occ = 10, # minimum number of unique occurrences
  thin_occ = NULL, # spatial filtering to reduce autocorrelation in occurrence records
  eval_occ = NULL, # external database to perform model evaluation based on independent data (if available)
  colin_var = c(method="PCA"), # method to reduce collinearity
  imp_var = FALSE, # assess variable importance
  sp_accessible_area = c(method='BUFFER', type='1'), #  restrict for each species the accessible area (current time)
  pseudoabs_method=c(method="RND"), # pseudo-absence allocation method
  pres_abs_ratio = 0.5, #  presence-absence ratio
  part=c(method='BOOT', replicates='2', proportion='0.7'), # cross-validation partition method
  save_part = FALSE, # save results for individual partition folds
  save_final = TRUE,
  algorithm = c("BIO", "MXD", "GLM"), # selected algorithms
  thr=c('JACCARD'),
  msdm = c(method="OBR"), # apply spatial restrictions to model projection (future)
  ensemble = c(method="W_MEAN", metric='Jaccard'), # method used to ensemble different algorithms
  extrapolation = FALSE,
  cores = cl
)

# Identify which species were not modelled:
OccRecords<-data.table::fread("Project Folder/Results/Occurrences_Cleaned.txt", stringsAsFactors=TRUE, encoding="UTF-8")
OccRecords<-OccRecords[, .N, by=sp] # count number of rows per species

# Get the name of the species modelled:
TargetSpecies<-dir("Project Folder/Results/Ensemble/W_MEAN/MSDMPosterior/")
TargetSpecies<-gsub(".tif", "", TargetSpecies) # remove file extension from the name
TargetSpecies<-TargetSpecies[TargetSpecies!="BIN"]

# Create an operator to get the opposite of %in%:
'%ni%' <- Negate('%in%')
OccRecords[OccRecords$sp %ni% TargetSpecies,]
# 14 species were excluded from ENM modelling framework, all with < 10 records.

#####

# Step 2: EXTRACT THE PRESENCE-ABSENCE MATRIX BASED ON ENM OUTPUTS FOR THE CURRENT CLIMATE
########################################################################################################################
# Paso 2: EXTRAER LA MATRIZ DE PRESENCIA-AUSENCIA EN BASE A LAS SALIDAS ENM PARA EL CLIMA ACTUAL
rm(list=ls()); gc()

# Check the raster files produced to represent species distributions as binary maps (current time):
dir_current_ENMs<-("Project Folder/Results/Ensemble/W_MEAN/MSDMPosterior/BIN/")

# Create a raster stack with species distributions for the current time:
current_ENMs<-raster::stack(list.files(path = dir_current_ENMs, pattern='.tif', full.names=T))

# ENM output for current climate is not masked by the species accessible area (extent area mask).
ExtentMask<-raster::stack(list.files(path="Project Folder/Results/Extent_Masks/", pattern='.tif', full.names=T))

# Apply the extent mask of each species iteratively:
current_ENMs_masked<-list()
for(i in 1:length(ExtentMask@layers)){
  current_ENMs_masked[[i]]<-raster::mask(x = current_ENMs[[i]], mask = ExtentMask[[i]])
  print(i)
}
rm(current_ENMs)

# Convert the list of masked rasters to a raster stack:
current_ENMs_masked<-do.call(raster::stack, current_ENMs_masked)

# Crop the raster stack to consider only the study area:
study_area<-rgdal::readOGR(dsn="Shapefiles", layer='Ecoregions2017')
study_area<-study_area[study_area@data$ECO_NAME=="Dry Chaco",]
current_ENMs_masked<-raster::crop(current_ENMs_masked, study_area) # the output is a raster brick
current_ENMs_masked<-stack(current_ENMs_masked) # convert it back to raster stack

# Some species might be present in the background, but not in the study area. Remove them:
remove_species<-which(maxValue(current_ENMs_masked)!=1)
current_ENMs_masked<-raster::dropLayer(current_ENMs_masked, remove_species)

# Recall that the masked raster stack has latlong as coordinate system and WGS84 as datum:
crs(current_ENMs_masked)

# Let's reclassify raster cells to replace '0' by 'NA', while keeping the original crs.
current_ENMs_wgs84<-list()

# Reclassify each layer and replace it in the 'current_ENMs_wgs84' object:
for(i in 1:length(current_ENMs_masked@layers)){
  r1 <- reclassify(current_ENMs_masked[[i]], cbind(-Inf, 0, NA), right=TRUE)
  
  # Rebuild the reclassified raster to match extent of the original raster stack:
  current_ENMs_wgs84[[i]] <- raster(vals=values(r1),
                                    ext=extent(current_ENMs_masked),
                                    crs=crs(current_ENMs_masked),
                                    nrows=dim(current_ENMs_masked)[1],
                                    ncols=dim(current_ENMs_masked)[2])
  
  print(i); rm(r1)
  
}

# Rename layers:
names(current_ENMs_wgs84)<-names(current_ENMs_masked)

# The line below does the same reclassification, but it takes 3-5x more time:
# current_ENMs_masked[current_ENMs_masked <= 0] <- NA

# Visualize the differences between non-reclassified and reclassified rasters:
par(mfrow = c(1, 2))
plot(current_ENMs_masked[[1]])
plot(current_ENMs_wgs84[[1]])
par(mfrow = c(1, 1))
rm(current_ENMs_masked)

# Transform raster projection to an equal area coordinate reference system (crs):
equalareaproj<-"+proj=cea +lon_0=0 +lat_ts=30 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs" # cylindrical equal area projection
wgs84<-"+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0" 
current_ENMs_cea<-lapply(current_ENMs_wgs84, function(x) # for each list within the object 'current_ENMs_wgs84', apply the function below
  projectRaster(x, crs=equalareaproj, method="ngb"))
rm(current_ENMs_wgs84)

# Extract some basic information on CEA rasters representing the study area (for later use):
data_resolution<-res(current_ENMs_cea[[1]]) # get the raster resolution
data_extent<-extent(current_ENMs_cea[[1]]) # get the raster extent
data_dimension<-dim(current_ENMs_cea[[1]]) # get the raster dimension

# Extract the list of species presences per grid cell:
SppPresences_list <- lapply(current_ENMs_cea, function(x) rasterToPoints(x, function(x){x==1}))
# View(SppPresences_list[["Amerotyphlops_reticulatus"]])

# Unpack the list of data.frames and rename them according to the respective list name:
SppPresences_df <- plyr::ldply(SppPresences_list, .id= "Sp")
rm(SppPresences_list)

# Convert data.frame from long to wide format:
SppPresences_df <- reshape2::dcast(SppPresences_df, x + y ~ Sp) 
SppPresences_df[is.na(SppPresences_df)]<-0 # replace NA with 0's
PAM_Current<-data.frame(latlong=paste0(SppPresences_df$y, "_", SppPresences_df$x), # unique values for each cell
                        SppPresences_df) # get the PAM table
# The PAM table include pixels outside the study area.
# Let's create a grid cell system for the study area and select the target occurrences (cells).

# Create the grid cell system:
study_area_cea<-sp::spTransform(study_area, CRSobj=equalareaproj) # get the study area in the same crs

# For a raster including only grid cells within the study area, run the code line below:
r<-raster(vals=1, ext=data_extent, crs=equalareaproj, nrows=data_dimension[1], ncols=data_dimension[2])
study_area_raster<-rasterize(study_area_cea, r, background=NA)
grid_cells<-rasterToPolygons(study_area_raster, dissolve=F)
# Alternatively, for a squared raster around the study area run the code line below:
# grid_cells<-rasterToPolygons(r, dissolve=F)
plot(grid_cells)

# Add coordinates in the attribute table of the 'grid_cells' shapefile:
grid_as_points<-as.data.frame(coordinates(grid_cells))
grid_cells@data$CellId<-1:nrow(grid_cells@data)
grid_cells@data$x<-grid_as_points$V1
grid_cells@data$y<-grid_as_points$V2
grid_cells@data<-grid_cells@data[,-1]
View(grid_cells@data)

# Visualize the grid_cells using an interactive map:
study_area_sf<-st_transform(sf::st_as_sf(study_area), crs=equalareaproj) # certify data is the same projection
grid_cells_sf<-st_transform(sf::st_as_sf(grid_cells), crs=equalareaproj) # certify data is the same projection
mapview::mapview(study_area_sf) + 
  mapview(grid_cells_sf, alpha.regions = 0.5)

# Get the coordinates of grid cells within the study area:
grid_cells@data$latlong<-paste0(grid_cells@data$y, "_", grid_cells@data$x)
fwrite(grid_cells@data, "Datasets/PAM_GridCells.csv") # save to disk

# Merge the column 'CellId' of the PAM_Current in the PAM_future:
PAM_Current<-merge(grid_cells@data[,c("CellId", "latlong")], # x
                   PAM_Current, # y
                   by="latlong", # column to match
                   all.y=TRUE # keep all rows in the y dataframe
)

# Remove cells outside the study area:
PAM_Current<-as.data.table(PAM_Current[!is.na(PAM_Current$CellId),])
PAM_Current<-PAM_Current[order(PAM_Current$CellId, decreasing = FALSE), -c("latlong")]

# Save the Presence-Absence matrix for the current time in the folder "Datasets":
data.table::fwrite(PAM_Current, "Datasets/PAM_CurrentTime.csv")

# Save the grid cell system as a shapefile in the folder "Shapefiles":
writeOGR(grid_cells, "Shapefiles", "grid_cells", driver="ESRI Shapefile", overwrite_layer=TRUE)

#####

# Step 3: COMPUTE THE ENSEMBLE ENM OUTPUTS ACROSS GLOBAL CLIMATE PROJECTIONS
########################################################################################################################
# Paso 3: CALCULE LOS RESULTADOS DEL ENSEMBLE ENM A TRAVÉS DE LAS PROYECCIONES CLIMÁTICAS
rm(list=ls()); gc()

# Create a raster stack with species distributions for the future climate (ssp5 & RCP8.5):
ssp585_Proj1_ENMs<-raster::stack(list.files(path = dir_ssp585_Proj1, pattern='.tif', full.names=T))
ssp585_Proj2_ENMs<-raster::stack(list.files(path = dir_ssp585_Proj2, pattern='.tif', full.names=T))
ssp585_Proj3_ENMs<-raster::stack(list.files(path = dir_ssp585_Proj3, pattern='.tif', full.names=T))

# Create empty lists to receive the masked rasters:
ssp585_Proj1_Masked<-list()
ssp585_Proj2_Masked<-list()
ssp585_Proj3_Masked<-list()

# Apply the extent mask of each species iteratively:
ExtentMask<-raster::stack(list.files(path="Project Folder/Results/Extent_Masks/", pattern='.tif', full.names=T))
for(i in 1:length(ExtentMask@layers)){
  ssp585_Proj1_Masked[[i]]<-raster::mask(x = ssp585_Proj1_ENMs[[i]], mask = ExtentMask[[i]])
  ssp585_Proj2_Masked[[i]]<-raster::mask(x = ssp585_Proj2_ENMs[[i]], mask = ExtentMask[[i]])
  ssp585_Proj3_Masked[[i]]<-raster::mask(x = ssp585_Proj3_ENMs[[i]], mask = ExtentMask[[i]])
  svMisc::progress(i) # print loop progress
}

# Visualize a projection before and after the masking operation:
par(mfrow = c(1, 3))
plot(ssp585_Proj1_ENMs[[1]], main="Future Projection\n without mask", legend=FALSE)
plot(ExtentMask[[1]], main="Extent Area Mask\nSpecies accessible area", legend=FALSE)
plot(ssp585_Proj1_Masked[[1]], main="Future Projection\n Masked", legend=FALSE)
par(mfrow = c(1, 1))

# Remove unnecessary objects:
rm(ssp585_Proj1_ENMs, ssp585_Proj2_ENMs, ssp585_Proj3_ENMs, dir_ssp585_Proj1, dir_ssp585_Proj2, dir_ssp585_Proj3)

# Convert the list of rasters to raster stack:
ssp585_Proj1_Masked<-do.call(raster::stack, ssp585_Proj1_Masked)
ssp585_Proj2_Masked<-do.call(raster::stack, ssp585_Proj2_Masked)
ssp585_Proj3_Masked<-do.call(raster::stack, ssp585_Proj3_Masked)

# Get the average habitat suitability for each species in each future scenario:
ssp585_AvgGCM<-list() # empty list to receive averaging rasters
ssp585_AvgGCM_sd<-list() # empty list to receive averaging rasters
N_spp<-length(names(ExtentMask)) # get the number of species (=iterations)
for(i in 1:N_spp){
  
  # Get the average raster cell values across GCMs for each species:
  ssp585_AvgGCM[[i]]<-raster::calc(
    stack(ssp585_Proj1_Masked[[i]], 
          ssp585_Proj2_Masked[[i]],
          ssp585_Proj3_Masked[[i]]),
    fun=mean, na.rm=T) # it is important to set na.rm=T
  
  # Get the standard deviation of raster cell values across GCMs for each species:
  ssp585_AvgGCM_sd[[i]]<-raster::calc(
    stack(ssp585_Proj1_Masked[[i]], 
          ssp585_Proj2_Masked[[i]],
          ssp585_Proj3_Masked[[i]]),
    fun=sd, na.rm=T) # sd = standard deviation
  
  svMisc::progress(i) # print loop progress
}

# Use species names to rename lists:
names(ssp585_AvgGCM)<-names(ExtentMask)
names(ssp585_AvgGCM_sd)<-names(ExtentMask)
rm(ssp585_Proj1_Masked, ssp585_Proj2_Masked, ssp585_Proj3_Masked, ExtentMask) # clean workspace

# Visualize some maps of after averaging masked raster projections:
for(i in 1:10){
  par(mfrow = c(1, 2))
  plot(ssp585_AvgGCM[[i]], main="Average future Projection\n across GCMs", legend=FALSE)
  plot(ssp585_AvgGCM_sd[[i]], main="SD of future Projection\n across GCMs", legend=FALSE)
  par(mfrow = c(1, 1))
  Sys.sleep(3)
}

# Information about the thresholds used to create the presence-absence maps for ensembled models:
Thresholds_Ensemble<-data.table::fread("Project Folder/Results/Thresholds_Ensemble.txt", stringsAsFactors=TRUE, encoding="UTF-8")

# Select one threshold to binarize habitat suitability maps:
Thresholds_Ensemble<-Thresholds_Ensemble[THR=="JACCARD",]
head(Thresholds_Ensemble)

# Binarize the habitat suitability maps for each species: 
for (i in 1:N_spp) {
  
  # Get the average threshold value for the target species and binarize the habitat suitability:
  Thr_Alg <- as.numeric(Thresholds_Ensemble[i,]$THR_VALUE)
  ssp585_AvgGCM[[i]] <- ssp585_AvgGCM[[i]] >= Thr_Alg
  svMisc::progress(i) # print loop progress
}

# Visualize some maps after binarization:
for(i in 1:10){plot(ssp585_AvgGCM[[i]]); Sys.sleep(2)}

# Future projections were not subject to spatial restrictions within the species accessible area.
# Use the occurrence records to constrain suitable patches in the future.
# Load data on presence records only for species used in ENMTML function:
OccPresencesTrain<-data.table::fread("Project Folder/Results/Occurrences_Fitting.txt", 
                                     stringsAsFactors=TRUE, encoding="UTF-8")[PresAbse==1,c("sp" ,"x", "y")]
OccPresencesTest<-data.table::fread("Project Folder/Results/Occurrences_Evaluation.txt", 
                                    stringsAsFactors=TRUE, encoding="UTF-8")[PresAbse==1,c("sp" ,"x", "y")]
OccPresences<-unique(rbind(OccPresencesTrain, OccPresencesTest))
rm(OccPresencesTrain, OccPresencesTest)

# The MSDM::MSDM_Posteriori() function does not allow the use of user-defined binarization thresholds.
# Use an adapted version developed for this course.
source("AdaptedRFunctions.R")

# Correction of overpredicting of the models in the scenario RCP 4.5:
MSDM_output <- MSDM_adapted(records = OccPresences, # data.frame with presence records (sp, x, y)
                            SpNames = names(ssp585_AvgGCM), # vector with species names (= sp)
                            RasterList = ssp585_AvgGCM, # list of binary rasters
                            Dispersion_Dist = NULL) # distance threshold

# Get the value of dispersal distance threshold computed for each species:
MSDM_Dist.Threshold<-MSDM_output$Dispersion_t
MSDM_Dist.Threshold<-data.frame(sp=names(MSDM_Dist.Threshold), dist.t=unlist(MSDM_Dist.Threshold))

# Extract the list of rasters representing the MSDM and convert to raster stack:
MSDM_stack<-raster::stack(MSDM_output$MSDM[[1]])
for(i in 2:length(MSDM_output$MSDM)){
  MSDM_stack<-raster::stack(MSDM_stack, MSDM_output$MSDM[[i]])
  svMisc::progress(i, progress.bar=TRUE) # print loop progress
} # faster than using do.call
names(MSDM_stack)<-names(MSDM_output$MSDM)
rm(MSDM_output)

# Some species may have lost all suitable areas in the future:
remove_species<-which(maxValue(MSDM_stack)!=1)
MSDM_stack<-raster::dropLayer(MSDM_stack, remove_species) # if not removed, they will be represent by an empty raster next.
Nspp_FinalSet<-dim(MSDM_stack)[3] # get the number of species in the final set

# Save the species rasters produced with MSDM_Posterior approach:
new_dir<-"Project Folder/Results/Projection/Ensemble_MSDMPosterior/"
dir.create(new_dir, showWarnings = FALSE)
for(i in 1:Nspp_FinalSet){
  raster::writeRaster(MSDM_stack[[i]],
                      filename=paste0(new_dir, "/", names(MSDM_stack)[i], ".tif"),
                      format="GTiff", overwrite=TRUE)
  svMisc::progress(i, progress.bar=TRUE) # print loop progress
}

#####

# Step 4: EXTRACT THE PRESENCE-ABSENCE MATRIX BASED ON ENM OUTPUTS FOR THE FUTURE CLIMATE
########################################################################################################################
# Paso 4: EXTRAER LA MATRIZ DE PRESENCIA-AUSENCIA BASADA EN RESULTADOS DE ENM PARA EL CLIMA FUTURO
rm(list=ls()); gc()

# Check the raster files produced to represent species distributions as binary maps (current time):
dir_fut_ensemble<-("Project Folder/Results/Projection/Ensemble_MSDMPosterior/")

# Create a raster stack with species distributions for the current time:
Future_ENMs<-raster::stack(list.files(path = dir_fut_ensemble, pattern='.tif', full.names=T))

# Crop the raster stack to consider only the study area:
study_area<-rgdal::readOGR(dsn="Shapefiles", layer='Ecoregions2017')
study_area<-study_area[study_area@data$ECO_NAME=="Dry Chaco",]
Future_ENMs<-raster::crop(Future_ENMs, study_area) # the output is a raster brick
Future_ENMs<-stack(Future_ENMs) # convert it back to raster stack

# Some species might be present in the background, but not in the study area. Remove them:
remove_species<-which(maxValue(Future_ENMs)!=1)
Future_ENMs<-raster::dropLayer(Future_ENMs, remove_species)

# Let's reclassify raster cells to replace '0' by 'NA', while keeping the original crs.
Future_ENMs_wgs84<-list()

# Reclassify each layer and replace it in the 'Future_ENMs_wgs84' object:
for(i in 1:length(Future_ENMs@layers)){
  r1 <- reclassify(Future_ENMs[[i]], cbind(-Inf, 0, NA), right=TRUE)
  
  # Rebuild the reclassified raster to match extent of the original raster stack:
  Future_ENMs_wgs84[[i]] <- raster(vals=values(r1),
                                   ext=extent(Future_ENMs),
                                   crs=crs(Future_ENMs),
                                   nrows=dim(Future_ENMs)[1],
                                   ncols=dim(Future_ENMs)[2])
  svMisc::progress(i); rm(r1)
}

# Rename layers:
names(Future_ENMs_wgs84)<-names(Future_ENMs)

# The line below does the same reclassification, but it takes 3-5x more time:
# Future_ENMs[Future_ENMs <= 0] <- NA

# Visualize the differences between non-reclassified and reclassified rasters:
par(mfrow = c(1, 2))
plot(Future_ENMs[[1]])
plot(Future_ENMs_wgs84[[1]])
par(mfrow = c(1, 1))
rm(Future_ENMs)

# Transform raster projection to an equal area coordinate reference system (crs):
equalareaproj<-"+proj=cea +lon_0=0 +lat_ts=30 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs" # cylindrical equal area projection
wgs84<-"+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0" 
Future_ENMs_cea<-lapply(Future_ENMs_wgs84, function(x) # for each list within the object 'Future_ENMs_wgs84', apply the function below
  projectRaster(x, crs=equalareaproj, method="ngb"))
rm(Future_ENMs_wgs84)

# Extract the list of species presences per grid cell:
SppPresences_list <- lapply(Future_ENMs_cea, function(x) rasterToPoints(x, function(x){x==1}))

# Unpack the list of data.frames and rename them according to the respective list name:
SppPresences_df <- plyr::ldply(SppPresences_list, .id= "Sp")
rm(SppPresences_list)

# Convert data.frame from long to wide format:
SppPresences_df <- reshape2::dcast(SppPresences_df, x + y ~ Sp) 
SppPresences_df[is.na(SppPresences_df)]<-0 # replace NA with 0's
PAM_Future<-data.frame(latlong=paste0(SppPresences_df$y, "_", SppPresences_df$x), # unique values for each cell
                       SppPresences_df) # get the PAM table

# Load the presence-absence matrix for the study area:
grid_cells<-rgdal::readOGR(dsn="Shapefiles", layer='grid_cells')

# Merge the column 'CellId' of the PAM_Future in the PAM_future:
PAM_Future<-merge(grid_cells[,c("CellId", "latlong")], PAM_Future, by="latlong", all.y=TRUE)

PAM_Future<-merge(grid_cells@data[,c("CellId", "latlong")], # x
                  PAM_Future, # y
                  by="latlong", # column to match
                  all.y=TRUE # keep all rows in the y dataframe
)




# Remove cells outside the study area:
PAM_Future<-as.data.table(PAM_Future[!is.na(PAM_Future$CellId),])
PAM_Future<-PAM_Future[order(PAM_Future$CellId, decreasing = FALSE), -c("latlong")]

# Save the Presence-Absence matrix for the current time in the folder "Datasets":
data.table::fwrite(PAM_Future, "Datasets/PAM_FutureTime.csv")

#####

# Step 5: COMPUTE THE SPECIES GEOGRAPHICAL RANGE SHIFTS BETWEEN TEMPORALLY DIFFERENT SCENARIOS
########################################################################################################################
# Paso 5: CALCULAR LOS CAMBIOS EN EL ALCANCE GEOGRÁFICO DE LAS ESPECIES ENTRE ESCENARIOS TEMPORALMENTE DIFERENTES
rm(list=ls()); gc()

# Check the raster files produced to represent species distributions as binary maps (current time):
dir_current_ENMs<-("Project Folder/Results/Ensemble/W_MEAN/MSDMPosterior/BIN/")
dir_future_ENMs<-("Project Folder/Results/Projection/Ensemble_MSDMPosterior/")

# Create a raster stack with species distributions for the current time:
current_ENMs<-raster::stack(list.files(path = dir_current_ENMs, pattern='.tif', full.names=T))
future_ENMs<-raster::stack(list.files(path = dir_future_ENMs, pattern='.tif', full.names=T))

# Some species might be present before binarization, but not be included after binarization. Remove them:
remove_species_curr<-which(maxValue(current_ENMs)!=1)
current_ENMs<-raster::dropLayer(current_ENMs, remove_species_curr)
remove_species_fut<-which(maxValue(future_ENMs)!=1)
future_ENMs<-raster::dropLayer(future_ENMs, remove_species_fut)

# Reclassify raster cells to replace '0' by 'NA'.
current_ENMs_wgs84<-list()
future_ENMs_wgs84<-list()
for(i in 1:length(current_ENMs@layers)){
  r1 <- reclassify(current_ENMs[[i]], cbind(-Inf, 0, NA), right=TRUE)
  current_ENMs_wgs84[[i]] <- raster(vals=values(r1),
                                    ext=extent(current_ENMs),
                                    crs=crs(current_ENMs),
                                    nrows=dim(current_ENMs)[1],
                                    ncols=dim(current_ENMs)[2])
  rm(r1)
  svMisc::progress(i)
}
for(i in 1:length(future_ENMs@layers)){
  r2 <- reclassify(future_ENMs[[i]], cbind(-Inf, 0, NA), right=TRUE)
  future_ENMs_wgs84[[i]] <- raster(vals=values(r2),
                                   ext=extent(future_ENMs),
                                   crs=crs(future_ENMs),
                                   nrows=dim(future_ENMs)[1],
                                   ncols=dim(future_ENMs)[2])
  rm(r2)
  svMisc::progress(i)
}

# Rename layers:
names(current_ENMs_wgs84)<-names(current_ENMs)
names(future_ENMs_wgs84)<-names(future_ENMs)
rm(current_ENMs, future_ENMs)

# Convert the raster layers to polygons (it may take ~9 minutes):
current_ENMs_pol <- lapply(current_ENMs_wgs84, function(x) raster::rasterToPolygons(x, function(x){x == 1}, na.rm=TRUE, digits=3, dissolve=T))
future_ENMs_pol <- lapply(future_ENMs_wgs84, function(x) raster::rasterToPolygons(x, function(x){x == 1}, na.rm=TRUE, digits=3, dissolve=T))

# Check some summary information for the resulting polygons:
summary(current_ENMs_pol[[1]]); summary(current_ENMs_pol[[2]])
# Note that the extent differs across polygons
# The attribute table contains one column named 'layer' filled with 1's

# Refill the column 'layer' to contain species name:
for(i in 1:length(current_ENMs_pol)){current_ENMs_pol[[i]]$layer<-names(current_ENMs_pol)[[i]]}
for(i in 1:length(future_ENMs_pol)){future_ENMs_pol[[i]]$layer<-names(future_ENMs_pol)[[i]]}

# Transform all SpatialPolygons to an equal area projection:
equalareaproj<-"+proj=cea +lon_0=0 +lat_ts=30 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs" # cylindrical equal area projection
current_ENMs_cea <- lapply(current_ENMs_pol, function(x) sf::as_Spatial(st_transform(sf::st_as_sf(x), crs=equalareaproj)))
future_ENMs_cea <- lapply(future_ENMs_pol, function(x) sf::as_Spatial(st_transform(sf::st_as_sf(x), crs=equalareaproj)))
# current_ENMs_cea2 <- lapply(current_ENMs_pol, function(x) sp::spTransform(x, CRSobj=equalareaproj)) # also works
rm(current_ENMs_pol, future_ENMs_pol)

# For each species, calculate area of polygons (km2):
area_curr <- lapply(current_ENMs_cea , function(x) ceiling(raster::area(x)/1000000)) # divide by 1000000 to convert m2 to km2
area_curr <- plyr::ldply(area_curr, data.table,.id= "Sp") # convert list to dataframe
names(area_curr)<-c("Sp", "AreaKm2_Current")
area_fut <- lapply(future_ENMs_cea , function(x) ceiling(raster::area(x)/1000000)) # divide by 1000000 to convert m2 to km2
area_fut <- plyr::ldply(area_fut, data.table,.id= "Sp") # convert list to dataframe
names(area_fut)<-c("Sp", "AreaKm2_Future")

# Merge information for the current and future climates:
spp_area<-merge(area_curr, area_fut, by="Sp", all.x=T)
summary(spp_area)
rm(current_ENMs_cea, future_ENMs_cea, area_curr, area_fut)


### SAME AS ABOVE, BUT WITHIN THE STUDY AREA ONLY
# Repeat all steps above but now cropping the study area to your study area.
study_area<-rgdal::readOGR(dsn="Shapefiles", layer='Ecoregions2017')
study_area<-study_area[study_area@data$ECO_NAME=="Dry Chaco",]

# Convert the list of rasters into a raster stack object and crop to the study area extent:
current_stack<-raster::stack(current_ENMs_wgs84[[1]])
for(i in 2:length(current_ENMs_wgs84)){
  current_stack<-raster::stack(current_stack, current_ENMs_wgs84[[i]])
  svMisc::progress(i)}
names(current_stack)<-names(current_ENMs_wgs84)
future_stack<-raster::stack(future_ENMs_wgs84[[1]])
for(i in 2:length(future_ENMs_wgs84)){
  future_stack<-raster::stack(future_stack, future_ENMs_wgs84[[i]])
  svMisc::progress(i)}
names(future_stack)<-names(future_ENMs_wgs84)
rm(current_ENMs_wgs84, future_ENMs_wgs84)
current_stack<-raster::stack(raster::crop(current_stack, study_area))
future_stack<-raster::stack(raster::crop(future_stack, study_area))

# Remove pixels outside the study area by applying a mask to the raster layers:
# First, create the mask of the study area:
r<-raster(vals=1, ext=extent(current_stack), 
          crs=crs(current_stack), 
          nrows=dim(current_stack)[1], 
          ncols=dim(current_stack)[2])
study_area_raster<-rasterize(study_area, r, background=NA)
rm(study_area)

# Apply the study area mask to each species iteratively, for both time periods:
current_stack_masked<-list()
for(i in 1:length(current_stack@layers)){
  current_stack_masked[[i]]<-raster::mask(x = current_stack[[i]], mask = study_area_raster)
}
future_stack_masked<-list()
for(i in 1:length(future_stack@layers)){
  future_stack_masked[[i]]<-raster::mask(x = future_stack[[i]], mask = study_area_raster)
}

# Name the lists with their respective species name:
names(current_stack_masked)<-names(current_stack)
names(future_stack_masked)<-names(future_stack)

# Keep only those species with predicted presences within the study area:
keep_current_spp<-names(unlist(lapply(current_stack_masked, function(x){which(maxValue(x)!=0)})))
current_stack_masked<-current_stack_masked[names(current_stack_masked) %in% keep_current_spp]
keep_future_spp<-names(unlist(lapply(future_stack_masked, function(x){which(maxValue(x)!=0)})))
future_stack_masked<-future_stack_masked[names(future_stack_masked) %in% keep_future_spp]

# Convert the raster layers to polygons and refill the column 'layer' with the respective species name:
current_ENMs_pol <- lapply(current_stack_masked, function(x) raster::rasterToPolygons(x, function(x){x == 1}, na.rm=TRUE, digits=3, dissolve=T))
future_ENMs_pol <- lapply(future_stack_masked, function(x) raster::rasterToPolygons(x, function(x){x == 1}, na.rm=TRUE, digits=3, dissolve=T))
for(i in 1:length(current_ENMs_pol)){current_ENMs_pol[[i]]$layer<-names(current_ENMs_pol)[[i]]}
for(i in 1:length(future_ENMs_pol)){future_ENMs_pol[[i]]$layer<-names(future_ENMs_pol)[[i]]}

# Transform all SpatialPolygons to an equal area projection and compute the area (km2):
equalareaproj<-"+proj=cea +lon_0=0 +lat_ts=30 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs" # cylindrical equal area projection
current_ENMs_cea <- lapply(current_ENMs_pol, function(x) sf::as_Spatial(st_transform(sf::st_as_sf(x), crs=equalareaproj)))
future_ENMs_cea <- lapply(future_ENMs_pol, function(x) sf::as_Spatial(st_transform(sf::st_as_sf(x), crs=equalareaproj)))
rm(current_ENMs_pol, future_ENMs_pol)

# For each species, calculate area of polygons (km2) for the current scenario:
area_curr <- lapply(current_ENMs_cea , function(x) ceiling(raster::area(x)/1000000))
area_curr <- plyr::ldply(area_curr, data.table,.id= "Sp") # convert list to dataframe
names(area_curr)<-c("Sp", "WithinAreaKm2_Current")
area_fut <- lapply(future_ENMs_cea , function(x) ceiling(raster::area(x)/1000000))
area_fut <- plyr::ldply(area_fut, data.table,.id= "Sp") # convert list to dataframe
names(area_fut)<-c("Sp", "WithinAreaKm2_Future")

# Merge information of species range shifts and save it to disk:
spp_area<-merge(spp_area, area_curr, by="Sp", all.x=T)
spp_area<-merge(spp_area, area_fut, by="Sp", all.x=T)
spp_area[is.na(spp_area)]<-0 # replace NA with 0's
rm(current_ENMs_cea, future_ENMs_cea, area_curr, area_fut)
fwrite(spp_area, "Datasets/SpeciesRangeShifts.csv")
beepr::beep(8)

#####

# Step 6: COMPUTE BETA-DIVERSITY METRICS USING THE PRESENCE-ABSENCE MATRIX OF THE CURRENT SCENARIO
########################################################################################################################
# Paso 6: CALCULE LAS MÉTRICAS DE BETA-DIVERSIDAD UTILIZANDO LA MATRIZ DE PRESENCIA-AUSENCIA DEL ESCENARIO ACTUAL
rm(list=ls()); gc()

# Load the database:
PAM_Current<-data.table::fread("Datasets/PAM_CurrentTime.csv", stringsAsFactors=TRUE, encoding="UTF-8")

# Load spatial data for interactive visualization of the grid system in use:
study_area<-rgdal::readOGR(dsn="Shapefiles", layer='Ecoregions2017')
study_area<-study_area[study_area@data$ECO_NAME=="Dry Chaco",]
grid_cells<-rgdal::readOGR(dsn="Shapefiles", layer='grid_cells')

# Visualize the grid_cells using an interactive map:
equalareaproj<-"+proj=cea +lon_0=0 +lat_ts=30 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs" # cylindrical equal area projection
study_area_sf<-st_transform(sf::st_as_sf(study_area), crs=equalareaproj) # certify data is the same projection
grid_cells_sf<-st_transform(sf::st_as_sf(grid_cells), crs=equalareaproj) # certify data is the same projection
mapview::mapview(study_area_sf) + 
  mapview(grid_cells_sf, alpha.regions = 0.5)

# Remove cells outside the study area (if grid cells were defined using a broader extent/spatial configuration):
intersection <- st_intersection(study_area_sf, grid_cells_sf)
'%ni%' <- Negate('%in%')
PAM_Current[PAM_Current$CellId %ni% intersection$CellId,]$CellId # visualize cells to be discarded
PAM_Current<-PAM_Current[PAM_Current$CellId %in% intersection$CellId,]

# Convert the 'PAM_current' data.table to matrix format and remove column 'CellId':
PAM<-as.matrix(PAM_Current[,2:ncol(PAM_Current)])

# Load just one raster layer to extract information on data spatial resolution:
current_ENMs<-raster::raster("Project Folder/Results/Ensemble/W_MEAN/MSDMPosterior/BIN/Amerotyphlops_brongersmianus.tif")
current_ENMs<-raster::projectRaster(current_ENMs, crs=equalareaproj, method="ngb")

# Define the radius of the searching window around the focal grid cell (based on the raster resolution):
window_size<-ceiling(res(current_ENMs[[1]])[2]/1000)*1000 # meters

# Select one focal cell and its respective surrounding cells based on the window_size defined above:
summary(PAM_Current$CellId)
FocalCell <- CommEcol::select.window(xf=PAM[which(PAM_Current$CellId==100),1], 
                                     yf=PAM[which(PAM_Current$CellId==100),2], 
                                     radius=window_size, # use raster resolution
                                     xydata=PAM)# [,-c(1,2)]

# Extract the multiple-site beta-diversity metrics between the focal cell and its surroundings:  
BetaOutput<-betapart::beta.multi(x=FocalCell[,3:ncol(FocalCell)], # PAM matrix, without coords or CellID
                                 index.family="jac" # Dissimilarity Index
)
BetaOutput

# Grid cells located at the edge of the study area are surrounded by less neighbouring sites.
# The difference in the total number of sites might affect beta-diversity metrics.
# It is recommended to sample the same number of sites around all focal cells.
BetaOutput<-betapart::beta.sample(x=FocalCell[,3:ncol(FocalCell)], index.family="jac",
                                  sites = 5, # number of cells
                                  samples = 10 # how many times
)
BetaOutput$mean.values

# What happens if there are fewer cells that the number of 'sites' specified: 
BetaOutput<-betapart::beta.sample(x=FocalCell[,3:ncol(FocalCell)], index.family="jac",
                                  sites = 7, # number of cells
                                  samples = 10 # how many times
)
BetaOutput$mean.values

# Identify which grid cells have at least four neighbouring cells:
cl<-makePSOCKcluster(detectCores()-1, type="SOCK")
registerDoParallel(cl)
getDoParWorkers()
CellsToUse<-foreach(i = 1:nrow(PAM_Current), 
                    .export = 'c',
                    .packages = c("CommEcol", "data.table")) %dopar% {
                      
                      FocalCell <- CommEcol::select.window(xf=PAM[i,1], yf=PAM[i,2], 
                                                           radius=window_size, xydata=PAM)
                      
                      MyCells<-data.frame(CellId=PAM_Current[i,]$CellId, # CellId
                                          NCells=(nrow(FocalCell)-1)) # Number of neighbouring cells
                      
                      MyCells
                    }
CellsToUse<-rbindlist(CellsToUse)
CellsToUse<-CellsToUse[CellsToUse$NCells>=4,]

# Export objects that will be useful later:
save(CellsToUse, window_size, file = "RData/CellsToUse.RData")

# Extract beta-diversity metrics in parallel: 
BetaCurrent<-foreach(i = 1:nrow(PAM_Current), 
                     .export = 'c',
                     .packages = c("betapart", "CommEcol", "data.table")) %dopar% {
                       
                       # Compute beta-diversity only if the minimum number of cells is reached:
                       if(PAM_Current$CellId[i] %in% CellsToUse$CellId){
                         
                         # Select the focal cell and its surroundings:
                         FocalCell <- CommEcol::select.window(xf=PAM[i,1], yf=PAM[i,2], radius=window_size, xydata=PAM)
                         
                         # Certify that the regional cells have at least two species (GammaRichness >=2):
                         if(ncol(FocalCell)>=4){ # first two columns of FocalCell are lat-long.
                           
                           # Compute multiple-site beta-diversity:
                           BetaMulti<-betapart::beta.sample(x=FocalCell[,3:ncol(FocalCell)], index.family="jac", sites = 4,  samples = 10)
                           BetaMulti<-unlist(BetaMulti[[2]])
                           BetaMulti<-data.frame(CellId=PAM_Current$CellId[i],
                                                 AlphaRichness=rowSums(PAM_Current[i,4:ncol(PAM_Current)]),
                                                 GammaRichness=ncol(FocalCell),
                                                 Beta.JAC=BetaMulti[3],
                                                 Beta.JTU=BetaMulti[1],
                                                 Beta.JNE=BetaMulti[2])
                           
                           BetaMulti
                           
                         } # end of if condition (on minimum spp. richness)
                         
                       } # end of if condition (on CellsToUse)
                       
                     } # end of foreach

parallel::stopCluster(cl)

# Unbind lists and get the relative contribution of JTU component to JAC:
BetaCurrent<-rbindlist(BetaCurrent)
BetaCurrent$Beta.Ratio<-BetaCurrent$Beta.JTU/BetaCurrent$Beta.JAC

# Beta.ratio can assume the value of zero if Beta.JAC = 0. 
BetaCurrent[is.na(BetaCurrent$Beta.Ratio),]
BetaCurrent[BetaCurrent$Beta.JAC==0,]$Beta.Ratio<-0 # manually set values to zero

#####

# Step 7: COMPUTE BETA-DIVERSITY METRICS BETWEEN CURRENT AND FUTURE SCENARIOS
########################################################################################################################
# Paso 7: CALCULE LAS MÉTRICAS DE DIVERSIDAD BETA ENTRE ESCENARIOS ACTUALES Y FUTUROS
rm(list=setdiff(ls(), c("BetaCurrent", "CellsToUse", "intersection", "PAM_Current", "window_size")))

# Load the PAM for the future climate and remove cells outside study area:
PAM_Future<-data.table::fread("Datasets/PAM_FutureTime.csv", stringsAsFactors=TRUE, encoding="UTF-8")
#remove last two rows


PAM_Future<-PAM_Future[PAM_Future$CellId %in% intersection$CellId,]
PAM<-as.matrix(PAM_Future[,2:ncol(PAM_Future)])

# Extract beta-diversity metrics in parallel: 
cl<-makePSOCKcluster(detectCores()-1, type="SOCK")
registerDoParallel(cl)
getDoParWorkers()
BetaFuture<-foreach(i = 1:nrow(PAM_Future), 
                    .export = 'c',
                    .packages = c("betapart", "CommEcol", "data.table")) %dopar% {
                      
                      # Compute beta-diversity only if the minimum number of cells is reached:
                      if(PAM_Future$CellId[i] %in% CellsToUse$CellId){
                        
                        # Select the focal cell and its surroundings:
                        FocalCell <- CommEcol::select.window(xf=PAM[i,1], yf=PAM[i,2], radius=window_size, xydata=PAM)
                        
                        # Certify that the regional cells have at least two species (GammaRichness >=2):
                        if(ncol(FocalCell)>=4){ # first two columns of FocalCell are lat-long.
                          
                          # Compute multiple-site beta-diversity:
                          BetaMulti<-betapart::beta.sample(x=FocalCell[,3:ncol(FocalCell)], index.family="jac", sites = 4,  samples = 10)
                          BetaMulti<-unlist(BetaMulti[[2]])
                          BetaMulti<-data.frame(CellId=PAM_Future$CellId[i],
                                                AlphaRichness=rowSums(PAM_Future[i,4:ncol(PAM_Future)]),
                                                GammaRichness=ncol(FocalCell),
                                                Beta.JAC=BetaMulti[3],
                                                Beta.JTU=BetaMulti[1],
                                                Beta.JNE=BetaMulti[2])
                          
                          BetaMulti
                          
                        } # end of if condition (on minimum spp. richness)
                        
                      } # end of if condition (on CellsToUse)
                      
                    } # end of foreach

# Unbind lists and get the relative contribution of JTU component to JAC:
BetaFuture<-rbindlist(BetaFuture)
BetaFuture$Beta.Ratio<-BetaFuture$Beta.JTU/BetaFuture$Beta.JAC
BetaFuture[BetaFuture$Beta.JAC==0,]$Beta.Ratio<-0 # manually set values to zero

# Certify that the same set of cells is used to compare present and future PAM's:
BetaFuture<-BetaFuture[BetaFuture$CellId %in% BetaCurrent$CellId,] # remove cells present only in the future

# Compute spatio-temporal changes in species composition (biotic homogenization):
SpatioTemporalBeta<-data.frame(CellId=BetaCurrent$CellId,
                               DeltaAlphaRich = (BetaFuture$AlphaRichness-BetaCurrent$AlphaRichness),
                               DeltaGammaRich = (BetaFuture$GammaRichness-BetaCurrent$GammaRichness),
                               DeltaJAC = (BetaFuture$Beta.JAC-BetaCurrent$Beta.JAC),
                               DeltaBRatio = (BetaFuture$Beta.Ratio-BetaCurrent$Beta.Ratio)
)

# Temporal changes in species composition (beta-diversity):
TemporalBeta<-betapart::beta.temp(x=PAM_Current[,4:ncol(PAM_Current)],
                                  y=PAM_Future[,4:ncol(PAM_Future)],
                                  index.family = "jac") # ERROR ALERT

# It is necessary for both PAM to have exactly the same dimension.
# Identify which species are missing in the current and future climates:
'%ni%' <- Negate('%in%')
SppCurrent<-names(PAM_Current)[4:ncol(PAM_Current)]
SppFuture<-names(PAM_Future)[4:ncol(PAM_Future)]
SppMissingCur<-SppFuture[SppFuture %ni% SppCurrent] # no species occurs only in the future climate
SppMissingFut<-SppCurrent[SppCurrent %ni% SppFuture] # 17 species occur only in the current climate

# Create an empty data.frame with missing species:
SppFuture<-as.data.frame(matrix(nrow=nrow(PAM_Current), ncol=length(SppMissingFut)))
names(SppFuture)<-SppMissingFut
SppFuture[,]<-0 # insert zero values

# Bind the missing species to the data.frame in which they are missing (repeat the process if necessary):
PAM_Future2<-as.data.frame(cbind(PAM_Future, SppFuture))

# Put the columns in both PAM's in the same order:
PAM_Future2<-PAM_Future2[, colnames(PAM_Current)]
(colnames(PAM_Future2)==colnames(PAM_Current)) # double-check

# Temporal changes in species composition (pairwise beta-diversity):
TemporalBeta<-betapart::beta.temp(x=PAM_Current[,4:ncol(PAM_Current)],
                                  y=PAM_Future2[,4:ncol(PAM_Future2)],
                                  index.family = "jac")

# Add information about 'CellId' and species richness:
TemporalBeta<-data.frame(CellId=PAM_Current$CellId,
                         CurRichness=rowSums(PAM_Current[,4:ncol(PAM_Current)]),
                         FutRichness=rowSums(PAM_Future2[,4:ncol(PAM_Future2)]),
                         Beta.jac=TemporalBeta$beta.jac,
                         Beta.jtu=TemporalBeta$beta.jtu,
                         Beta.jne=TemporalBeta$beta.jne,
                         Beta.ratio=TemporalBeta$beta.jtu/TemporalBeta$beta.jac)

# Save the results:
fwrite(BetaCurrent, "Datasets/BetaCurrent.csv")
fwrite(BetaFuture, "Datasets/BetaFuture.csv")
fwrite(TemporalBeta, "Datasets/TemporalBeta.csv")
fwrite(SpatioTemporalBeta, "Datasets/SpatioTemporalBeta.csv")

######

# STEP 8: COMPUTE PHYLOGENETIC-BASED BIODIVERSITY METRICS AT THE ASSEMBLAGE-LEVEL
########################################################################################################################
# STEP 8: COMPUTE PHYLOGENETIC-BASED BIODIVERSITY METRICS AT THE ASSEMBLAGE-LEVEL
rm(list=ls()); gc()
source("AdaptedRFunctions.R")

# Load presence-absence data on each temporal scenario and also info on target cells (e.g. CellId field of TemporalBeta):
PAM_Current<-data.table::fread("Datasets/PAM_CurrentTime.csv", stringsAsFactors=TRUE, encoding="UTF-8")
PAM_Future<-data.table::fread("Datasets/PAM_FutureTime.csv", stringsAsFactors=TRUE, encoding="UTF-8")
PAM_Future<- PAM_Future[-c(3557, 3558),]
TemporalBeta<-data.table::fread("Datasets/TemporalBeta.csv", stringsAsFactors=TRUE, encoding="UTF-8")
PAM_Current<-PAM_Current[PAM_Current$CellId %in% TemporalBeta$CellId,]
PAM_Future<-PAM_Future[PAM_Future$CellId %in% TemporalBeta$CellId,]
rm(TemporalBeta)


# Certify that both PAMs have the same columns (species) in the same order:
'%ni%' <- Negate('%in%')
SppCurrent<-names(PAM_Current)[4:ncol(PAM_Current)] # species list in the current time
SppFuture<-names(PAM_Future)[4:ncol(PAM_Future)] # species list in the future time
SppMissingCur<-SppFuture[SppFuture %ni% SppCurrent] # species found only in the future climate (n = 0)
SppMissingFut<-SppCurrent[SppCurrent %ni% SppFuture] # species found only in the present climate(n = 15)
SppFuture<-as.data.frame(matrix(nrow=nrow(PAM_Current), ncol=length(SppMissingFut))) # empty dataframe
names(SppFuture)<-SppMissingFut # rename columns as species name
SppFuture[,]<-0 # insert zero values
PAM_Future<-as.data.frame(cbind(PAM_Future, SppFuture)) # bind data.frames
PAM_Future<-PAM_Future[, colnames(PAM_Current)] # order columns 
rm(SppCurrent, SppMissingCur, SppFuture, SppMissingFut)

# Load the phylogenetic tree:
phylo_tree<-ape::read.tree('Phylogenies/PrunedTree.tre')

# Compute the temporal pairwise phylogenetic beta-diversity between present and future scenarios:

TemporalPhyloBeta<-tempbetagrid(oc1 = PAM_Current[,4:ncol(PAM_Current)], # PAM data.frame current climate
                                oc2 = PAM_Future[,4:ncol(PAM_Future)], # PAM data.frame future climate
                                index = "jaccard", # also works with "sorensen"
                                phylotree = phylo_tree, # phylogenetic tree
                                phylobeta = TRUE # compute phylobeta
)

# Compute the phylobeta-diversity across species assemblages for the current time.
# To use a fast function for PD calculation, it is necessary to convert PAM to sparse-format.
pam<-as.data.frame(PAM_Current[,4:ncol(PAM_Current)])
rownames(pam)<-PAM_Current$CellId
pam_sparse<-phyloregion::dense2sparse(pam)

# Extract the PD and the respective Standardized Effect Size:
PD_Current<-phyloregion::PD_ses(x=pam_sparse, # PAM data (sparse format)
                                phy=phylo_tree, # phylogenetic tree
                                model="tipshuffle", # null model 
                                reps=1000)

# Same as above, but for the future scenario:
pam<-as.data.frame(PAM_Future[,4:ncol(PAM_Future)])
rownames(pam)<-PAM_Future$CellId
pam_sparse<-phyloregion::dense2sparse(pam)
PD_Future<-phyloregion::PD_ses(x=pam_sparse, # PAM data (sparse format)
                               phy=phylo_tree, # phylogenetic tree
                               model="tipshuffle", # null model 
                               reps=1000)
rm(pam, pam_sparse)

# Combine PD metrics and TemporalPhyloBeta:
TemporalPhyloBeta<-data.frame(CellId=PAM_Current$CellId,
                              
                              # Current scenario
                              SR_Curr=PD_Current$richness,
                              PD_Curr=PD_Current$PD_obs,
                              SESPD_Curr=PD_Current$zscore,
                              SESPD_p_Curr=PD_Current$pd_obs_p,
                              
                              # Future scenario
                              SR_Fut=PD_Future$richness,
                              PD_Fut=PD_Future$PD_obs,
                              SESPD_Fut=PD_Future$zscore,
                              SESPD_p_Fut=PD_Future$pd_obs_p,
                              
                              # Inter-scenarios
                              Beta.jac=TemporalPhyloBeta$beta,
                              Beta.jtu=TemporalPhyloBeta$turnover,
                              Beta.jne=TemporalPhyloBeta$nestedness,
                              Beta.ratio=TemporalPhyloBeta$turn.beta)

# Save the results:
fwrite(TemporalPhyloBeta, "Datasets/TemporalPhyloBeta.csv")

#####

# Step 9: CALCULATE TEMPORAL CHANGES IN SPATIAL GRADIENTS OF PHYLOBETA DIVERSITY
########################################################################################################################
# Paso 9: CALCULE LOS CAMBIOS TEMPORALES EN LOS GRADIENTES ESPACIALES DE LA DIVERSIDAD DE FILOBETA
rm(list=ls()); gc()

# The 'betapart' package has a multiple-site beta-diversity function: phylo.beta.multi().
# However, 'betapart' package does not include the phylogenetic version of beta.sample().
# Hence, it is necessary to apply phylo.beta.multi() iteratively.

# Load presence-absence data on each temporal scenario and also info on target cells (e.g. CellId field of TemporalBeta):
PAM_Current<-data.table::fread("Datasets/PAM_CurrentTime.csv", stringsAsFactors=TRUE, encoding="UTF-8")
PAM_Future<-data.table::fread("Datasets/PAM_FutureTime.csv", stringsAsFactors=TRUE, encoding="UTF-8")
TemporalBeta<-data.table::fread("Datasets/TemporalBeta.csv", stringsAsFactors=TRUE, encoding="UTF-8")
PAM_Current<-PAM_Current[PAM_Current$CellId %in% TemporalBeta$CellId,]
PAM_Future<-PAM_Future[PAM_Future$CellId %in% TemporalBeta$CellId,]
rm(TemporalBeta)

# Load the phylogenetic tree and R-objects specifying which CellId's to use:
phylo_tree<-ape::read.tree('Phylogenies/PrunedTree.tre')
load("RData/CellsToUse.RData")

# Prepare R for parallel computation:
cl<-makePSOCKcluster(detectCores()-1, type="SOCK")
registerDoParallel(cl)
getDoParWorkers()

# Extract phylobeta-diversity metrics for the current scenario in parallel: 
PAM<-as.matrix(PAM_Current[,2:ncol(PAM_Current)])
PhyloBetaCurrent<-foreach(i = 1:nrow(PAM_Current), 
                          .export = 'c',
                          .packages = c("betapart", "CommEcol", "data.table")) %dopar% {
                            
                            # Compute beta-diversity only if the minimum number of cells is reached:
                            if(PAM_Current$CellId[i] %in% CellsToUse$CellId){
                              
                              # Select the focal cell and its surroundings:
                              FocalCell <- CommEcol::select.window(xf=PAM[i,1], yf=PAM[i,2], radius=window_size, xydata=PAM)
                              
                              # Certify that the regional cells have at least two species (GammaRichness >=2):
                              if(ncol(FocalCell)>=4){ # first two columns of FocalCell are lat-long.
                                
                                # Subsample sites within the regional cells
                                BetaMulti<-list()
                                for(j in 1:10){ # resampling 10 times
                                  
                                  # Subsample sites within the regional cells:
                                  set.seed(j) 
                                  SampledCell<-FocalCell[sample(x = c(1:nrow(FocalCell)),
                                                                size = 4,
                                                                replace = FALSE),]
                                  
                                  # Compute the multiple-site phylobeta diversity: 
                                  BetaMulti[[j]]<-betapart::phylo.beta.multi(x=SampledCell[,c(3:ncol(SampledCell))], 
                                                                             tree=phylo_tree,
                                                                             index.family="jac")
                                } # end of for loop
                                
                                BetaMulti<-rbindlist(BetaMulti) # bind outputs for each iteration
                                BetaMulti<-colMeans(BetaMulti, na.rm=T) # average beta-diversity metrics across iterations
                                
                                # Store in a datamframe:
                                BetaMulti<-data.frame(CellId=PAM_Current$CellId[i],
                                                      AlphaRichness=rowSums(PAM_Current[i,4:ncol(PAM_Current)]),
                                                      GammaRichness=ncol(FocalCell),
                                                      PBeta.JAC=BetaMulti[3],
                                                      PBeta.JTU=BetaMulti[1],
                                                      PBeta.JNE=BetaMulti[2])
                                
                                # Return output:
                                BetaMulti
                                
                              } # end of if condition (on minimum spp. richness)
                              
                            } # end of if condition (on CellsToUse)
                            
                          } # end of foreach

# Extract phylobeta-diversity metrics for the future scenario in parallel: 
PAM<-as.matrix(PAM_Future[,2:ncol(PAM_Future)])
PhyloBetaFuture<-foreach(i = 1:nrow(PAM_Future), 
                         .export = 'c',
                         .packages = c("betapart", "CommEcol", "data.table")) %dopar% {
                           
                           # Compute beta-diversity only if the minimum number of cells is reached:
                           if(PAM_Future$CellId[i] %in% CellsToUse$CellId){
                             
                             # Select the focal cell and its surroundings:
                             FocalCell <- CommEcol::select.window(xf=PAM[i,1], yf=PAM[i,2], radius=window_size, xydata=PAM)
                             
                             # Certify that the regional cells have at least two species (GammaRichness >=2):
                             if(ncol(FocalCell)>=4){ # first two columns of FocalCell are lat-long.
                               
                               # Subsample sites within the regional cells
                               BetaMulti<-list()
                               for(j in 1:10){ # resampling 10 times
                                 
                                 # Subsample sites within the regional cells:
                                 set.seed(j) 
                                 SampledCell<-FocalCell[sample(x = c(1:nrow(FocalCell)),
                                                               size = 4,
                                                               replace = FALSE),]
                                 
                                 # Compute the multiple-site phylobeta diversity: 
                                 BetaMulti[[j]]<-betapart::phylo.beta.multi(x=SampledCell[,c(3:ncol(SampledCell))], 
                                                                            tree=phylo_tree,
                                                                            index.family="jac")
                               } # end of for loop
                               
                               BetaMulti<-rbindlist(BetaMulti) # bind outputs for each iteration
                               BetaMulti<-colMeans(BetaMulti, na.rm=T) # average beta-diversity metrics across iterations
                               
                               # Store in a datamframe:
                               BetaMulti<-data.frame(CellId=PAM_Future$CellId[i],
                                                     AlphaRichness=rowSums(PAM_Future[i,4:ncol(PAM_Future)]),
                                                     GammaRichness=ncol(FocalCell),
                                                     PBeta.JAC=BetaMulti[3],
                                                     PBeta.JTU=BetaMulti[1],
                                                     PBeta.JNE=BetaMulti[2])
                               
                               # Return output:
                               BetaMulti
                               
                             } # end of if condition (on minimum spp. richness)
                             
                           } # end of if condition (on CellsToUse)
                           
                         } # end of foreach

parallel::stopCluster(cl)

# Unbind lists and get the relative contribution of JTU component to JAC:
PhyloBetaCurrent<-rbindlist(PhyloBetaCurrent)
PhyloBetaCurrent$PBeta.Ratio<-PhyloBetaCurrent$PBeta.JTU/PhyloBetaCurrent$PBeta.JAC
PhyloBetaCurrent[PhyloBetaCurrent$PBeta.JAC==0,]$PBeta.Ratio<-0 # mannually set values to zero
PhyloBetaFuture<-rbindlist(PhyloBetaFuture)
PhyloBetaFuture$PBeta.Ratio<-PhyloBetaFuture$PBeta.JTU/PhyloBetaFuture$PBeta.JAC
PhyloBetaFuture[PhyloBetaFuture$PBeta.JAC==0,]$PBeta.Ratio<-0 # mannually set values to zero

#rows not present in one data.frame
a<- setdiff(PhyloBetaFuture$CellId, PhyloBetaCurrent$CellId)
+3510-3489
#remove extra cellds not present in current
PhyloBetaFuture<- PhyloBetaFuture[-c(a),]




# Compute spatio-temporal changes in species composition (biotic homogenization):
SpatioTemporalPhyloBeta<-data.frame(CellId=PhyloBetaCurrent$CellId,
                                    DeltaAlphaRich = (PhyloBetaFuture$AlphaRichness-PhyloBetaCurrent$AlphaRichness),
                                    DeltaGammaRich = (PhyloBetaFuture$GammaRichness-PhyloBetaCurrent$GammaRichness),
                                    DeltaPBeta = (PhyloBetaFuture$PBeta.JAC-PhyloBetaCurrent$PBeta.JAC),
                                    DeltaPBRatio = (PhyloBetaFuture$PBeta.Ratio-PhyloBetaCurrent$PBeta.Ratio)
)

# Save the results:
fwrite(PhyloBetaCurrent, "Datasets/PhyloBetaCurrent.csv")
fwrite(PhyloBetaFuture, "Datasets/PhyloBetaFuture.csv")
fwrite(SpatioTemporalPhyloBeta, "Datasets/SpatioTemporalPhyloBeta.csv")

#####

# Step 10 - BUILD BARPLOTS OF EXPECTED CHANGES IN SPECIES RANGE AND ASSEMBLAGE GENERALISM 
##########################################################################################################################
# Paso 10 - CONSTRUIR GRÁFICOS DE BARRAS DE CAMBIOS ESPERADOS EN LA GAMA DE ESPECIES Y EL GENERALISMO DE CONJUNTO 
rm(list=ls()); gc() # clean the workspace

### BARPLOT 1 - SPECIES RANGE SHIFTS
# Load results obtained at the species- and assemblage-level:
SppData<-fread("Datasets/SpeciesRangeShifts.csv", stringsAsFactors=TRUE, encoding="UTF-8")

# Get the relative change in habitat suitability between present and future scenarios:
SppData$Rel_Change_TA<-SppData$AreaKm2_Future/SppData$AreaKm2_Current
SppData$Rel_Change_WA<-SppData$WithinAreaKm2_Future/SppData$WithinAreaKm2_Current

# Categorize species according to their estimated relative gain or loss of range size in the future:
ecdf_function<-ecdf(SppData$Rel_Change_TA)
SppData<-SppData %>% mutate(RangeShiftCat = cut(Rel_Change_TA, 
                                                breaks=quantile(SppData$Rel_Change_TA,
                                                                c(0, # >50% range loss
                                                                  ecdf_function(0.5), # between 50-0% of range loss 
                                                                  ecdf_function(1), # between 50-100% of range gain 
                                                                  ecdf_function(1.5), 1)),
                                                include.lowest=TRUE))

# Change label of levels within the 'SppData$RangeShiftCat':
SppData$RangeShiftCat<-factor(SppData$RangeShiftCat, 
                              levels = c(levels(SppData$RangeShiftCat)[1], # Species that will loss >50%
                                         levels(SppData$RangeShiftCat)[2], # Species that will loss 0-50%
                                         levels(SppData$RangeShiftCat)[3], # Species that will gain 0-50%
                                         levels(SppData$RangeShiftCat)[4] # Species that will gain >50%
                              ),
                              labels = c("Loss >50%", "Loss <50%", "Gain <50%", "Gain >50%"))

# Same as above but considering just the species range portion within the study area:
summary(SppData$Rel_Change_WA) # check if there are species present in the background but absent in the study area
SppData[is.na(SppData$Rel_Change_WA),] # visualize the species identified in the line above
SppData_noNA<-SppData[!is.na(SppData$Rel_Change_WA),]
ecdf_function<-ecdf(SppData_noNA$Rel_Change_WA)
SppData_noNA<-SppData_noNA %>% mutate(RangeShiftCat_WA = cut(Rel_Change_WA,  breaks=quantile(SppData_noNA$Rel_Change_WA, c(0, ecdf_function(0.5), 
                                                                                                                           ecdf_function(1), ecdf_function(1.5), 1)), include.lowest=TRUE))
SppData_noNA$RangeShiftCat_WA<-factor(SppData_noNA$RangeShiftCat_WA, 
                                      levels = c(levels(SppData_noNA$RangeShiftCat_WA)[1], levels(SppData_noNA$RangeShiftCat_WA)[2],
                                                 levels(SppData_noNA$RangeShiftCat_WA)[3], levels(SppData_noNA$RangeShiftCat_WA)[4]),
                                      labels = c("Loss >50%", "Loss <50%", "Gain <50%", "Gain >50%"))

# Compute marginal totals per generalism level:
Subtotal_1<-as.data.frame(SppData_noNA %>% dplyr::group_by(RangeShiftCat) %>% dplyr::summarise(N_Spp=length(Sp)))
Subtotal_1$PropSpp<-Subtotal_1$N_Spp/sum(Subtotal_1$N_Spp, na.rm=T)
Subtotal_1$RefArea<-"BackgroundArea"
Subtotal_2<-as.data.frame(SppData_noNA %>% dplyr::group_by(RangeShiftCat_WA) %>% dplyr::summarise(N_Spp=length(Sp)))
Subtotal_2$PropSpp<-Subtotal_2$N_Spp/sum(Subtotal_2$N_Spp, na.rm=T)
Subtotal_2$RefArea<-"StudyArea"
names(Subtotal_2)[1]<-"RangeShiftCat"

# Bind the marginal totals in a single table:
Subtotals<-rbind(Subtotal_1, Subtotal_2)
rm(Subtotal_1, Subtotal_2)

# Prepare colours:
myColors<-c("#E15944", # dark blue
            "#9F1651", # light blue
            "#430154", # light red
            "#FFB340" # dark red
)


# Build the barplot:
MyPlot1<-ggplot(data=Subtotals, aes(x=RefArea, y=PropSpp, fill=RangeShiftCat)) +
  geom_bar(position = position_fill(reverse = TRUE), stat="identity") +
  xlab("") + ylab("Proportion of species") +
  coord_flip() +  # flip coordinates (puts labels on y axis)
  scale_fill_manual(values = myColors,
                    name=bquote(atop("Proportional", "range shift"))) +
  theme(panel.grid.minor = element_blank(), # remove minor gridlines
        panel.grid.major = element_blank(), # remove major gridlines
        panel.background = element_blank(), # white background
        axis.line = element_line(colour="black"), # axis lines aesthetitcs
        axis.text.y = element_text(hjust=1, vjust=0.5, angle=0, size=10),
        axis.text.x = element_text(hjust=0.5, vjust=0.5, angle=0, size=10),
        axis.ticks.y=element_blank(),
        axis.title.x=element_text(size=12, colour="black", face="bold", margin=margin(t=5, r=0, b=5, l=0)), # margin between axis.title and axis.values
        axis.title.y=element_blank(),
        legend.position="top",
        legend.title=element_text(size=10, colour="black", margin=margin(t=5, r=10, b=5, l=5)),
        plot.background=element_blank(),
  ) +
  geom_vline(xintercept=1.5, color= "gray") +
  scale_y_continuous(limits=c(0, 1), breaks=c(0, 0.25, 0.5, 0.75, 1), 
                     expand=expansion(mult=c(0, .05), add=c(0, 0)),
                     labels=scales::number_format(accuracy = 0.01, decimal.mark='.')) +
  scale_x_discrete(position="bottom", labels=c("Within the\nbackground\narea", "Within the\nstudy area")) +
  guides(fill=guide_legend(nrow=2, byrow=F)); MyPlot1


### BARPLOT 2 - ASSEMBLAGE-LEVEL BIOTIC CHANGES
# Create a variable to represent ecological generalism (wide vs narrow distribution):
SppData<-as.data.frame(SppData)
summary(SppData$AreaKm2_Current)
SppData$Generalism<-NA
SppData[SppData$AreaKm2_Current>50000,]$Generalism<-"Wide"
SppData[SppData$AreaKm2_Current<=50000,]$Generalism<-"Narrow" # ERROR ALERT
summary(as.factor(SppData$Generalism))
# No species has less than 50 thousand km2.

#  For illustrative purposes, we will use the median 'AreaKm2_Current' to define narrowly distributed species:
SppData<-SppData %>% mutate(Generalism = cut(WithinAreaKm2_Current, 
                                             breaks=c(0, summary(SppData$WithinAreaKm2_Current)[3], max(SppData$WithinAreaKm2_Current, na.rm=T)),
                                             include.lowest=TRUE))

# Relabel levels of 'Generalism':
SppData$Generalism<-factor(SppData$Generalism,  
                           levels = c(levels(SppData$Generalism)[1], levels(SppData$Generalism)[2]),
                           labels = c("Narrow", "Wide"))
summary(as.factor(SppData$Generalism))

# Compute the relative proportion of widely-distributed species in each assemblage:
PAM_Current<-fread("Datasets/PAM_CurrentTime.csv", stringsAsFactors=TRUE, encoding="UTF-8")
pam_long<-reshape2::melt(PAM_Current, id.vars= c("CellId", "x", "y"), variable.name="Sp", value.name="Presence") # convert the PAM to long format
pam_long<-pam_long[pam_long$Presence==1,]
pam_long<-merge(pam_long, SppData[,c("Sp", "Generalism")], by='Sp', all.x=TRUE, allow.cartesian=TRUE) # add the species-level data to the PAM
pam_long<-as.data.table(pam_long) # convert to data.table for subsequent fast calculations
pam_long<-pam_long[, .(.N), by = .(CellId, Generalism)] # compute the richness of widely distributed species per assemblage
pam_long<-pam_long[Generalism=="Wide", c("CellId", "N")] # remove duplicated assemblages (= richness of narrowly-distributed spp.)
names(pam_long)<-c('CellId', 'WideRichness')

# Load data on spatio-temporal changes in phylobeta diversity and merge information on richness of widely-distributed:
GridData<-fread("Datasets/SpatioTemporalPhyloBeta.csv", stringsAsFactors=TRUE, encoding="UTF-8")
Richness<-fread("Datasets/PhyloBetaCurrent.csv", stringsAsFactors=TRUE, encoding="UTF-8")
GridData<-merge(GridData, Richness[,c("AlphaRichness", "CellId")], by="CellId", all.x=T) # add current richness data
GridData<-merge(GridData, pam_long, by="CellId", all.x=T)
rm(pam_long, Richness)

# Get the proportion of widely distributed species per assemblage:
GridData$WideProp<-GridData$WideRichness/GridData$AlphaRichness
summary(GridData$WideProp)

# Identify assemblages with higher dominance of widely distributed species (using the first quartile of 'WideProp'):
GridData<-GridData %>% mutate(RangeSizeDominance = cut(WideProp, breaks=c(0, summary(GridData$WideProp)[1], 2), include.lowest=TRUE))


# Relabel levels of 'RangeSizeDominance':
GridData$RangeSizeDominance<-factor(GridData$RangeSizeDominance,  
                                    levels = c(levels(GridData$RangeSizeDominance)[1], levels(GridData$RangeSizeDominance)[2]),
                                    labels = c("Narrow", "Wide"))


# Create a categorical variable to represent different levels of biotic changes:
GridData<-as.data.frame(GridData)
GridData$BioticChangeCat<-NA
GridData[GridData$DeltaPBeta<=(quantile(GridData$DeltaPBeta, probs=0.1, na.rm=T)),]$BioticChangeCat<-"Lower10"
GridData[GridData$DeltaPBeta>=(quantile(GridData$DeltaPBeta, probs=0.9, na.rm=T)),]$BioticChangeCat<-"Upper10"
GridData[is.na(GridData$BioticChangeCat) & GridData$DeltaPBeta<0,]$BioticChangeCat<-"Homogenization"
GridData[is.na(GridData$BioticChangeCat) & GridData$DeltaPBeta>=0,]$BioticChangeCat<-"Heterogenization"

# Compute marginal totals per generalism level:
Subtotal_1<-as.data.frame(GridData %>% dplyr::group_by(RangeSizeDominance, BioticChangeCat) %>% dplyr::summarise(Subtotal=length(CellId)))
Subtotal_2<-as.data.frame(GridData %>% dplyr::group_by(RangeSizeDominance) %>% dplyr::summarise(Total=length(CellId)))
Subtotals<-merge(Subtotal_1, Subtotal_2, by="RangeSizeDominance", all.x=T)
Subtotals$PropCells<-Subtotals$Subtotal/Subtotals$Total
head(Subtotals)
rm(Subtotal_1, Subtotal_2)

fwrite(GridData, "Datasets/GridDataBioticChanges.csv")
# Order levels of 'BioticChangeCat'
Subtotals$BioticChangeCat<-factor(Subtotals$BioticChangeCat,
                                  levels=c("Lower10", "Homogenization", "Heterogenization", "Upper10"),
                                  labels=c("Lowest 10%", "Homogenization", "Heterogenization", "Uppest 10%")
)

# Build the barplot:
MyPlot2<-ggplot(data=Subtotals, aes(x=RangeSizeDominance, y=PropCells, fill=BioticChangeCat)) +
  geom_bar(position = position_fill(reverse = TRUE), stat="identity") +
  coord_flip() +  # flip coordinates (puts labels on y axis)
  xlab("") + ylab("Proportion of assemblages") +
  scale_fill_manual(values = myColors,
                    name=bquote(atop("Biotic change","Phylo"*beta*"-diversity"))) +
  theme(panel.grid.minor = element_blank(), # remove minor gridlines
        panel.grid.major = element_blank(), # remove major gridlines
        panel.background = element_blank(), # white background
        axis.line = element_line(colour="black"), # axis lines aesthetitcs
        axis.text.y = element_text(hjust=1, vjust=0.5, angle=0, size=10),
        axis.text.x = element_text(hjust=0.5, vjust=0.5, angle=0, size=10),
        axis.ticks.y=element_blank(),
        axis.title.x=element_text(size=12, colour="black", face="bold", margin=margin(t=5, r=0, b=5, l=0)), # margin between axis.title and axis.values
        axis.title.y=element_blank(),
        legend.position="top",
        legend.title=element_text(size=10, colour="black", margin=margin(t=5, r=10, b=5, l=5)),
        plot.background=element_blank(),
  ) +
  geom_vline(xintercept=1.5, color= "gray") +
  scale_y_continuous(limits=c(0, 1), breaks=c(0, 0.25, 0.5, 0.75, 1), 
                     expand=expansion(mult=c(0, .05), add=c(0, 0)),
                     labels=scales::number_format(accuracy = 0.01, decimal.mark='.')) +
  scale_x_discrete(position="bottom", labels=c("More widely\ndistributed\nspecies",
                                               "Less widely\ndistributed\nspecies")) +
  guides(fill=guide_legend(nrow=2, byrow=F)); MyPlot2

# Multipanel barplot:
Barplot<-ggpubr::ggarrange(MyPlot1, MyPlot2,
                           labels=c("A", "B"), align="hv",
                           font.label=list(size=14, colour = "black"), ncol=1, nrow=2); Barplot

# Save to disk:
ggsave("Figures/Barplots.png", plot=Barplot, width=5, height=7, units="in", bg = "transparent")
ggsave("Figures/Barplots.pdf", plot=Barplot, width=5, height=7, units="in", bg = "transparent")

#####
