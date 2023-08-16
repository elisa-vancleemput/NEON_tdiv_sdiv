##############################################################################
###                                                                        ###
###             Main script for downloading and gathering all data         ###
###                                                                        ###
##############################################################################
# Code accompanying
# Van Cleemput, E., Adler, P., Suding, K. in Global Ecology and Geography
# Making remote sense of biodiversity: What grassland characteristics make spectral diversity a good proxy for taxonomic diversity?

# In this script, we download and prepare all data for analyses, step by step.
# The output of this script is used as input for 00_Main_script_analyses
# All steps should be performed for each of the sites separately, by indicating one of the four-digit NEON site codes as follows
site_code <- "ABBY"         
site_code <- "CLBJ"         
site_code <- "CPER"         
site_code <- "KONZ"         
site_code <- "NIWO"         
site_code <- "NOGP"         
site_code <- "OAES"
site_code <- "SJER"
site_code <- "WOOD"

# Scripts supporting this code
# - 01_Download_preprocess_aop_imagery.R
# - 02_Download_veg_plots_data.R
# - 03_Extract_plot_level_raster_values.R
# - 04_Calculate_spectral_diversity.R

# This script contains 8 parts, each part can be run independently upon opening this script

# PART 1: Prepare site polygon(s) 
# PART 2: Download and preprocess aop images
# PART 3: Exploratory Visualization of tiles with polygons
# PART 4: Download NEON vegetation plots information
# PART 5: Mosaic and filter aop data
# PART 6: Extract plot level information from rasters
# PART 7: Calculate spectral diversity for base plots
# PART 8: Merge data from different data sources

# This code was created and run in the CyVerse environment for computational and data size reasons
# - Open the CyVerse Discovery Environment (personal account needed)
# - Load and save the required shapefiles in the Data store, 
#   from your local computer, e.g. using Cyberduck, 
# - Open the RStudio environment via the "Apps": 
#   use the "Rocker_RStudio_Geospatial_3.6.3_analysis1" docker of Tyson Swetnam 
#   While opening the docker, specify the input folder containing the required polygon shapefiles
# - Open the R project "rangeland-div-stab" and run this script
# Important note: Working in Cyverse made me create a 2-map system, see more under "General lines..."

# Part of this script is based on 2020-neon-aop-workshop provided in Cyverse
# and code provided by NEON, Batelle:
# https://www.neonscience.org/download-explore-neon-data 
# https://github.com/elisa-vancleemput/neon-aop-workshop.git

# This script is also inspired by code written by Victoria Scholl:
# https://github.com/earthlab/neon-veg.git

# by Elisa Van Cleemput, 2023
#######################################################
# General lines to be run before each part  -----------
#######################################################
rm(list=ls())

# Set the main directory, this depends on the local computer you are working on, or is ~ in the case of remotely working on cyverse
server = "local" # or Cyverse or local

# load general packages
if (!require('rgdal')) install.packages('rgdal'); library('rgdal')
if (!require('raster')) install.packages('raster'); library('raster')
if (!require('dplyr')) install.packages('dplyr'); library('dplyr')

# Set global option to NOT convert all character variables to factors
options(stringsAsFactors=F)

# ---- General function to create folders ---- 
check_create_dir <- function(new_dir){
  # check if directory exists. If it doesn't, create it. 
  if (!dir.exists(new_dir)){
    dir.create(new_dir)
  }
}

# ---- Create folders ---- 
# In this script I work with 2 "data folders"
# - Data (capital D): A map in the Cyverse Discovery Environment, where I store all data, and which can be accessed from the Cyverse Data Store
# - data (small d): A map created during downloading and analysis, of which the content is later copied to "Data"
# The reason for this 2-map system is that it is impossible to edit existing maps in Cyverse in an active R session
# In the remainder of the script, the code will always check first if the output already exists in the Cyverse Data Store, and 
# if not, it will be downloaded/processed/stored into the data folder.

# Attention: to be able to read base items (e.g., polygons) in CyVerse R, copy them to your personal CyVerse account, 
# e.g. via Cyberduck, and set this folder as input folder before starting the script

# Set the main directory, this depends on the local computer you are working on,or is defined by Cyvers
if (server == "Cyverse") {

  main_dir <- "~/work/data" 
  main_dir_input <- file.path(main_dir,  "input") # this is the Data Store
  main_dir_output <- file.path(main_dir,  "output") # this is where Cyverse will save output which can afterwards be moved to the Data Store

  dir_scripts <- "~/NEON_tdiv_sdiv"
  
} else {
  main_dir_local <- # FILL IN
  main_dir_external <- # FILL IN
  
  main_dir_input <- main_dir_external
  main_dir_output <- main_dir_external
  
  dir_scripts <- file.path(main_dir_local, "NEON_tdiv_sdiv")
  
}

# Create folders to store the raw and processed data
check_create_dir(paste0(main_dir_output,"/data"))
dir_data_raw <- file.path(main_dir_output, "data", "data_raw")
dir_data_out <- file.path(main_dir_output, "data", "data_output")
check_create_dir(dir_data_raw)
check_create_dir(dir_data_out)

# Permanent storage maps
Dir_data_raw <- file.path(main_dir_input,"Data","data_raw")
Dir_data_out <- file.path(main_dir_input,"Data","data_output")
Dir_data_poly <- file.path(main_dir_input,"Data","NEON_boundaries")
Dir_data_GEE <- file.path(main_dir_input,"Data","GEE_results")

setwd(dir_scripts)

#######################################################
## PART 1: Prepare site polygon(s) --------------------
#######################################################
# USER-DEFINED-INPUTS
# Define parameters for the AOP remote sensing data to download
site_code <- "ABBY"         # Four-digit NEON site code, character string type
aop_data_year <- "2017"     # Four-digit year in character string "YYYY" for AOP imagery to download example image

if(file.exists(paste0(file.path(Dir_data_out,paste(site_code,aop_data_year,sep="_")),
                      "/poly_sp_proj_", site_code, "_", aop_data_year, ".shp"))){
  print(paste0("Shapefile of field boundaries for ", site_code, " already created and saved"))
} else if(file.exists(paste0(file.path(dir_data_out,paste(site_code,aop_data_year,sep="_")),
                        "/poly_sp_proj_", site_code, "_", aop_data_year, ".shp"))){
  print(paste0("Shapefile of field boundaries for ", site_code, " already created and saved"))
} else {
  # Shapefile of sites: field boundaries provided by NEON
  poly_name <- "Field_Sampling_Boundaries_2020"
  poly_name_shp <- "terrestrialSamplingBoundaries"
  
  ## ---- Open the polygon files  --------------
  if (!require('utils')) install.packages('utils'); library('utils')
  list.files(Dir_data_poly)
  unzip(paste0(Dir_data_poly,"/",poly_name,".zip"), exdir=dir_data_raw)
  
  poly_sp <- readOGR(paste0(dir_data_raw, "/", poly_name_shp,".shp"))
  crs(poly_sp)
  
  # access an overview of the sites via 
  FieldSites <- poly_sp@data[,c("siteName","siteID","domainNumb")]
  FieldSites[order(FieldSites[,'siteID']), ]
  plot(poly_sp[poly_sp$siteID == site_code,], border="red")
  
  ## ---- Clean polygon dataset --------------
  poly_sp_clean <- poly_sp
  poly_sp_clean$AREA <- area(poly_sp_clean)
  poly_sp_clean$HECTARES <- poly_sp_clean$AREA * 0.0001
  
  # Not all polygons provided in the shapefile are meaningful for our study, and we will remove them first
  # The spatial polygon dataframe for some sites contain some very small polygons
  # We will remove all polygons that are smaller than 20 x 20 m² (= NEON base plot size)
  poly_sp_clean <- subset(poly_sp_clean, AREA > 20^2) 
  
  plot(poly_sp_clean, border="blue")
  
  ## ---- Download 1 corresponding AOP tiles to identify aop CRS system --------------
  # We want the polygon shapefile to be in the same coordinate system as the AOP are provided in.
  # Therefore we will download 1 image, to extract the crs from
  # A message in the console will display the total file size to be downloaded.
  # Proceed by typing "y" and pressing Enter. 
  dp_hs_refl <- "DP3.30006.001"
  
  source("01_Download_preprocess_aop_imagery.R") # contains the function download_neon_refl_crsexample and "download_neon_aop"
  f <- download_neon_refl_crsexample(Dir_data_raw,dir_data_raw,site_code,aop_data_year)
  
  h5ls(f, all = TRUE)
  Metadata <- h5read(f, paste0("/",site_code,"/Reflectance/Metadata/"))
  # View(Metadata)
  h5read(f, paste0("/",site_code,"/Reflectance/Metadata/Spectral_Data"))$Wavelength
  # View(h5ls(f,all=T))
  aop_crs <- h5read(f, paste0("/",site_code,"/Reflectance/Metadata/Coordinate_System"))$Proj4
  aop_epsg <- h5read(f, paste0("/",site_code,"/Reflectance/Metadata/Coordinate_System"))$'EPSG Code'
  
  ## ---- Allign crs of polygons with that of the raster image --------------
  poly_sp_proj <- spTransform(poly_sp_clean,CRS(paste0("+init=epsg:",aop_epsg))) # aop_epsg is a character string and needs to be in the form of class CRS
  
  ## ---- Store the cleaned and projected spatial polygon data frame --------------
  dir_data_out_site <- file.path(dir_data_out,paste(site_code,aop_data_year,sep="_"))
  check_create_dir(dir_data_out_site)
  # shapefiles have a field name character limit of 10 characters, truncate because otherwise writeOGR will do this in a less predictable way
  names(poly_sp_proj) <- strtrim(names(poly_sp_proj),10)
  writeOGR(obj=poly_sp_proj, dsn=dir_data_out_site, layer=paste0("poly_sp_proj_", site_code, "_", aop_data_year), driver="ESRI Shapefile")
}

### END of PART 1

#######################################################
## PART 2: Download and preprocess aop images ---------
#######################################################
# USER-DEFINED-INPUTS
# Define parameters for the AOP remote sensing data to download
site_code     <- "ABBY"     # Four-digit NEON site code, character string type
aop_data_year <- "2017"     # Four-digit year in character string "YYYY" for AOP imagery to download example image
buffer_val    <- 0 #[m]     # integer buffer size around coordinates to download example imag. W
# We set it to zero because we will explicitly indicate the names of the tiles to be downloaded
n_sg          <- 7          # Window for Savitzky-Golay filter

## ---- ADD MORE HERE AND ALSO TO THE FUNCTIONS AGAIN? Specify NEON products --------------
# List of Data Product IDs that can be downloaded
dp_chm <- "DP3.30015.001"               # LiDAR-derived Canopy Height Model
dp_veg_indices <- "DP3.30026.001"       # Spectrometer-derived vegetation indices
dp_lai <- "DP3.30012.001"               # Spectrometer-derived leaf area index
dp_hs_refl <- "DP3.30006.001"           # Spectometer-derived surface reflectance
dp_aspect_slope <- "DP3.30025.001"      # LiDAR-derived Aspect and Slope
dp_elevation <- "DP3.30024.001"         # LiDAR-derived DTM
dp_lai <- "DP3.30012.001"               # Spectrometer-derived leaf area index
dp_fpar <- "DP3.30014.001"              # Spectrometer-derived fraction of incident photosynthetically active radiation (400-700 nm) absorbed by the green elements of a vegetation canopy

## ---- Open cleaned and projected SpatialPolygonsDataFrame --------------
if(file.exists(paste0(Dir_data_out,"/",paste0(site_code,"_",aop_data_year,"/"),
                      paste0("poly_sp_proj_", site_code, "_", aop_data_year,".shp")))){
  filename <- paste0(Dir_data_out,"/",paste0(site_code,"_",aop_data_year,"/"),
                     paste0("poly_sp_proj_", site_code, "_", aop_data_year,".shp"))
  poly_sp_proj <- readOGR(filename)
} else {
  filename <- paste0(dir_data_out,"/",paste0(site_code,"_",aop_data_year,"/"),
                     paste0("poly_sp_proj_", site_code, "_", aop_data_year,".shp"))
  poly_sp_proj <- readOGR(filename)
}

## ---- Define tile coordinates of the area to be downloaded --------------

# 1) Define easting and northing: vectors containing the easting/northing UTM coordinates of the locations to download.
# UTM:  coordinates are measured as northings and eastings in meters
# NEON stores tiles of 1 km by 1 km, so
# easting and northing vectors should have steps of max. 1000 m in order to download all adjacent tiles
# We code 500 m steps

eastings <- c(seq(poly_sp_proj@bbox[1,"min"],poly_sp_proj@bbox[1,"max"],by=500),poly_sp_proj@bbox[1,"max"])
northings <- c(seq(poly_sp_proj@bbox[2,"min"],poly_sp_proj@bbox[2,"max"],by=500),poly_sp_proj@bbox[2,"max"])

# For testing purposes: Download tiles of 1 specific polygon
# i = 19 # 19 is good example for JORN
# eastings <- c(seq(poly_sp_proj[i,]@bbox[1,"min"],poly_sp_proj[i,]@bbox[1,"max"],by=500),poly_sp_proj[i,]@bbox[1,"max"])
# northings <- c(seq(poly_sp_proj[i,]@bbox[2,"min"],poly_sp_proj[i,]@bbox[2,"max"],by=500),poly_sp_proj[i,]@bbox[2,"max"])

# 2) Combine eastings and northings to coordinates
poly_coordinates <- tidyr::expand_grid(eastings, northings)

# 3) Summarize this list to a list of the 1km x 1km tile coordinates covering the site
setwd(dir_scripts)
source("01_Download_preprocess_aop_imagery.R")
tile_coordinates <- list_tiles_covering_poly_tibble(poly_coordinates)

# If interested (e.g., for sharing) save the tile coordinates
# saveRDS(tile_coordinates,
#         paste0(main_dir,"/data","/data_output/","tile_coordinates_",site_code,".rds"))
# tile_coordinates <- readRDS(paste0(main_dir,"/data","/data_output/","tile_coordinates_",site_code,".rds"))

## ---- Define folder to store downloaded data temporally, while stack is being made --------------
dir_data_out_temp <- file.path(dir_data_raw,"Temp_AOPtiles")
check_create_dir(dir_data_out_temp)
Dir_data_out_temp <- file.path(Dir_data_raw,"Temp_AOPtiles")

## ---- Download hyperspectral reflectance, chm and veg_indices together --------------
# We add chm and veg_indices to this download, because this information will later on be used for masking (PART 5)

# If the data is already stacked and stored in Cyverse we do not have to download the files again
# Attention: this if-loop only makes sense when entire sites are being downloaded and stacked at once
# Beware, the following lines will take a loooooong time!
if (length(list.files(paste0(Dir_data_out, "/", site_code, "_", aop_data_year, "/stacked_aop_data"),all.files = TRUE)) != 0){
  message("The data of ", paste0(site_code, "_", aop_data_year), " was already downloaded, stacked and stored in Cyverse")
} else {
  source("01_Download_preprocess_aop_imagery.R")
  download_neon_aop(dp_chm,dp_aspect_slope=NA,dp_elevation=NA, 
                    dp_veg_indices,dp_lai=NA,dp_fpar=NA,
                    dp_hs_refl,
                    site_code,aop_data_year,
                    veg_coordinates = tile_coordinates,
                    buffer_val,
                    dir_data_out_temp,
                    check_size=F,delete_originals=F)
  
  ## ---- Preprocess hyperspectral images and stack all tile products with the same coordinates
  # This script:
  #  - removes bad bands and applies Savitzky-Golay filter to hyperspectral signatures
  #  - combines all airborne remote sensing data layers per square-km tile into a single multi-layer raster
  #  - saves the resulting stacked tiles as an .rds file. 
  dir_data_out_site <- file.path(dir_data_out,paste(site_code,aop_data_year,sep="_"))
  check_create_dir(dir_data_out_site)
  
  stack_neon_aop(tile_coordinates, dir_data_out_temp, dir_data_out_site, n_sg)
  
  ## ---- Save ancillary data
  stack_neon_aop_ancillary(tile_coordinates,dir_data_out_temp, dir_data_out_site)
}

# Empty or remove the Temp_AOPtiles map before downloading and stacking data of a new site and year

# read the text file with wavelength data for subsequent steps
wavelengths <- as.numeric(unlist(read.table(file.path(dir_data_out_site, "wavelengths.txt"),
                                            skip = 1,
                                            col.names = 'wavelength')))

# !!!!!! ATTENTION !!!!!!
# Due to masking before smoothing the spectra, reflectance values near the atmospheric water absorption windows are
# artificially lowered/raised. Therefore it is recommended to omit some extra bands in this region.
# We used a smoothing window of 7, so it is best to remove 3 bands one each side of the boundaries.
# We incorporated this in the section where we calculated spectral diversity

## ---- Download aop derivative products together --------------
if (length(list.files(paste0(Dir_data_out, "/", site_code, "_", aop_data_year, "/stacked_aop_data_derivatives"),all.files = TRUE)) != 0){
  message("The derivative data of ", paste0(site_code, "_", aop_data_year), " was already downloaded, stacked and stored in Cyverse")
} else {
  source("01_Download_preprocess_aop_imagery.R") 
  download_neon_aop(dp_chm=NA,dp_aspect_slope,dp_elevation, 
                    dp_veg_indices=NA,dp_lai,dp_fpar,
                    dp_hs_refl=NA,
                    site_code,aop_data_year,
                    veg_coordinates = tile_coordinates,
                    buffer_val,
                    dir_data_out_temp,
                    check_size=F,delete_originals=F)
  
  dir_data_out_site <- file.path(dir_data_out,paste(site_code,aop_data_year,sep="_"))
  check_create_dir(dir_data_out_site)
  
  stack_neon_aop_derivatives(tile_coordinates, Dir_data_out_temp, dir_data_out_site)
}

# Empty or remove the Temp_AOPtiles map before downloading and stacking data of a new site and year

### END of PART 2

#########################################################################
## PART 3: Exploratory Visualization of tiles with polygons -------------
#########################################################################
# USER-DEFINED-INPUTS
# Define parameters for the AOP remote sensing data
site_code     <- "ABBY"     # Four-digit NEON site code, character string type
aop_data_year <- "2017"     # Four-digit year in character string "YYYY" for AOP imagery to download example image

## ---- Open cleaned and projected SpatialPolygonsDataFrame --------------
if(file.exists(paste0(Dir_data_out, "/",paste0(site_code,"_",aop_data_year,"/"),
                      paste0("poly_sp_proj_", site_code, "_", aop_data_year,".shp")))){
  filename <- paste0(Dir_data_out, "/",paste0(site_code,"_",aop_data_year,"/"),
                     paste0("poly_sp_proj_", site_code, "_", aop_data_year,".shp"))
  poly_sp_proj <- readOGR(filename)
} else {
  filename <- paste0(dir_data_out, "/",paste0(site_code,"_",aop_data_year,"/"),
                     paste0("poly_sp_proj_", site_code, "_", aop_data_year,".shp"))
  poly_sp_proj <- readOGR(filename)
}

## ---- Specify which polygon to visualize --------------
# -- Option 1: Visualize 1 specific polygon
i = 1
polyname <- paste0(site_code,"_",aop_data_year,"_",poly_sp_proj$NAME_1[i])

# -- Option 2: Visualize a set of polygons
i = c(1,2,3)
polyname <- paste0(site_code,"_",aop_data_year,"_",
                   paste(poly_sp_proj$NAME_1[i],collapse="_"))

# -- Option 3: Visualize all polygons of a site
i <- seq(1:length(poly_sp_proj))
polyname <- paste0(site_code,"_",aop_data_year,"_allpolygons")

## ---- Define coordinates of the polygon(s) to be visualized --------------
# Make a list of unique easting, northing coordinates (500 m steps)
eastings <- c(seq(poly_sp_proj[i,]@bbox[1,"min"],poly_sp_proj[i,]@bbox[1,"max"],by=500),poly_sp_proj@bbox[1,"max"])
northings <- c(seq(poly_sp_proj[i,]@bbox[2,"min"],poly_sp_proj[i,]@bbox[2,"max"],by=500),poly_sp_proj@bbox[2,"max"])
poly_coordinates <- tidyr::expand_grid(eastings, northings)

## ---- Generate a list of the 1km x 1km tile coordinates covering the polygon(s) --------------
source("01_Download_preprocess_aop_imagery.R") 
tile_coordinates <- list_tiles_covering_poly(poly_coordinates)

## ---- Specify the names of the tiles (.rds stacked AOP data file) to plot --------------
# If the tiles were already downloaded and stored previously, access them via the Data folder,
# if they were downloaded during the current R session, access them via the data folder
# Attention: this if-loop only makes sense when entire sites were downloaded and stacked at once (see above)
if (length(list.files(paste0(Dir_data_out, "/", site_code, "_", aop_data_year, "/stacked_aop_data"),all.files = TRUE)) != 0){
  Dir_data_out_site_year <- file.path(main_dir,"Data","data_output",paste(site_code,aop_data_year,sep="_"))
  stacked_aop_data_filename <- file.path(Dir_data_out_site_year,
                                         "stacked_aop_data",
                                         paste0("stacked_aop_data_",
                                                tile_coordinates, ".rds"))
  wavelengths <- scan(paste0(Dir_data_out_site_year,"/wavelengths.txt"), skip=1,quiet=T)
} else {
  dir_data_out_site_year <- file.path(dir_data_out,paste(site_code,aop_data_year,sep="_"))
  stacked_aop_data_filename <- file.path(dir_data_out_site_year,
                                         "stacked_aop_data",
                                         paste0("stacked_aop_data_",
                                                tile_coordinates, ".rds"))
  wavelengths <- scan(paste0(dir_data_out_site_year,"/wavelengths.txt"), skip=1,quiet=T)
}

## ---- Create output folders --------------
dir_data_out_site <- file.path(dir_data_out,paste(site_code,aop_data_year,sep="_"))
check_create_dir(dir_data_out_site)
dir_data_out_maps <- file.path(dir_data_out,paste(site_code,aop_data_year,sep="_"),"maps")
check_create_dir(dir_data_out_maps)
dir_data_out_sign <- file.path(dir_data_out,paste(site_code,aop_data_year,sep="_"),"sign")
check_create_dir(dir_data_out_sign)

## ---- Plot and optionally safe --------------
plot_aop_imagery(filename = stacked_aop_data_filename, polyname, wavelengths, 
                 save_plots = FALSE, dir_data_out_maps, 
                 save_sign = TRUE, dir_data_out_sign,
                 plot_rgb_hs = T, plot_chm = T, 
                 plot_slope = F, plot_aspect = F,plot_indices = T,
                 poly_sp_proj = poly_sp_proj[i,], add_poly=T)

# !!! The plotting and orientation looks correct, but the x-axis counts up.
# Remember that the tiles are in UTM projection

### END of PART 3

##############################################################
## PART 4: Download NEON vegetation plots information --------
##############################################################
# USER-DEFINED-INPUTS
# Define parameters for the AOP remote sensing data to download
site_code     <- "ABBY"     # Four-digit NEON site code, character string type
aop_data_year <- "2017"     # Four-digit year in character string "YYYY" for AOP imagery to download example image
radius_m_base <- 10         # Plot radius used for extracting NEON baseplots based on their centroid coordinates

## ---- Specify NEON products --------------
# List of Data Product IDs to be downloaded
dp_spcomp <- "DP1.10058.001"      # Plant presence and % cover, on subplot level

## ---- Download NEON 20m-plot and 1m-subplot level species composition data --------------
source("02_Download_veg_plots_data.R")

# Download plant species composition and cover
if (file.exists(paste0(Dir_data_out,"/",site_code,"_",aop_data_year,"/",paste0(site_code, "_", aop_data_year,"_veg_plots_spcomp.rda")))){
  message("The data of ", paste0(site_code, "_", aop_data_year), " was already downloaded and stored")
} else {
  veg_plots_spcomp_site_year <- Extract_plot_spcomp(site_code, aop_data_year, dp_spcomp)

  dir_data_out_site <- (file.path(dir_data_out,paste0(site_code,"_",aop_data_year)))
  check_create_dir(dir_data_out_site)
  setwd(dir_data_out_site)
  saveRDS(veg_plots_spcomp_site_year, file=paste0(site_code, "_", aop_data_year,"_veg_plots_spcomp.rda"))
}

## ---- Data exploration: get some insights by looking at summary information of sp comp data --------------
### some summary data:
# Total number of species records in the site. The same species might be present multiple times if encountered in multiple subplots
veg_plots_spcomp_site_year$table_observation_species_subplots %>% nrow()
veg_plots_spcomp_site_year$table_observation_species_baseplots %>% nrow()

# species richness per subplot of 1 m x 1 m
sampling_effort_summary_subplot <- veg_plots_spcomp_site_year$table_observation_species_subplots %>%
  group_by(plotID, subplotID, sampleID, endDate) %>%
  # count number of species within each subplot
  summarise(event_count = taxonID %>% unique() %>% length())
View(sampling_effort_summary_subplot)

# species richness per plot of 20 m x 20 m calculated as the total richness of all 1 m² subplots in a plot
sampling_effort_summary_plot <- veg_plots_spcomp_site_year$table_observation_species_subplots %>%
  group_by(plotID, month, endDate) %>%
  # count number of species within each plot
  summarise(event_count = taxonID %>% unique() %>% length())
View(sampling_effort_summary_plot)

# Notes: 
# - KONZ_017 was sampled on July 11 (n species = 48) and August 22 (n species = 33)
# - KONZ_027 was sampled in January and September


### END of PART 4

## ---- Open cleaned and projected SpatialPolygonsDataFrame --------------
if(file.exists(paste0(Dir_data_out,"/",paste0(site_code,"_",aop_data_year,"/"),
                      paste0("poly_sp_proj_", site_code, "_", aop_data_year,".shp")))){
  filename <- paste0(Dir_data_out,"/",paste0(site_code,"_",aop_data_year,"/"),
                     paste0("poly_sp_proj_", site_code, "_", aop_data_year,".shp"))
  poly_sp_proj <- readOGR(filename)
} else {
  filename <- paste0(dir_data_out,"/",paste0(site_code,"_",aop_data_year,"/"),
                     paste0("poly_sp_proj_", site_code, "_", aop_data_year,".shp"))
  poly_sp_proj <- readOGR(filename)
}

## ---- Locate NEON's baseplots --------------
# Generate sample of plots coinciding with NEON vegetation base and subplots
Dir_data_out_plots <- file.path(Dir_data_out,paste0(site_code,"_",aop_data_year),"plot_locations")

if(file.exists(file.path(Dir_data_out_plots, paste0(site_code, "_veg_plots_sub", ".rds")))){
  print("NEON base- and subplots already located ans saved")
} else {
  
  # Read raw data on species composition
  if (file.exists(paste0(Dir_data_out,"/",site_code,"_",aop_data_year,"/",paste0(site_code, "_", aop_data_year,"_veg_plots_spcomp.rda")))){
    filename1 <- paste0(Dir_data_out,"/",site_code,"_",aop_data_year,"/",paste0(site_code, "_", aop_data_year,"_veg_plots_spcomp.rda"))
    veg_plots_spcomp_site_year <- readRDS(filename1)
  } else {
    filename1 <- paste0(dir_data_out,"/",site_code,"_",aop_data_year,"/",paste0(site_code, "_", aop_data_year,"_veg_plots_spcomp.rda"))
    veg_plots_spcomp_site_year <- readRDS(filename1)
  }
  
  # Retrieve coordinates
  source("02_Download_veg_plots_data.R")
  veg_plots_base_coords <- get_plot_coords(veg_plots_spcomp_site_year, plot_type = "base", radius_m_base, poly_sp_proj)

  # Remove plots that fall out of our scope
  veg_plots_base_coords$veg_spdf_point$id <- veg_plots_base_coords$veg_spdf_point$plotID
  veg_plots_base_coords$veg_utm$id <- veg_plots_base_coords$veg_utm$plotID
  veg_plots_base_coords_filtered <- remove_unwanted_neon_plots(veg_plots_base_coords, poly_sp_proj, veg_plots_spcomp_site_year)
  
  # Just a quick overview:
  veg_plots_polyinfo <- cbind("plotID" = veg_plots_base_coords_filtered$veg_spdf_point$plotID, over(veg_plots_base_coords_filtered$veg_spdf_point, poly_sp_proj))
  
  # Save the locations
  dir_data_out_site <- file.path(dir_data_out,paste0(site_code,"_",aop_data_year))
  check_create_dir(dir_data_out_site)
  dir_data_out_plots <- file.path(dir_data_out,paste0(site_code,"_",aop_data_year),"plot_locations")
  check_create_dir(dir_data_out_plots)
  
  filename_base = file.path(dir_data_out_plots, paste0(site_code, "_veg_plots_base", ".rds"))
  saveRDS(veg_plots_base_coords_filtered, file = filename_base)  
}

### END of PART 4

############################################################
## PART 5: Mosaic and filter aop data ----------------------
############################################################
# USER-DEFINED-INPUTS
# Define parameters for the AOP remote sensing data to download
site_code     <- "ABBY"     # Four-digit NEON site code, character string type
aop_data_year <- "2017"     # Four-digit year in character string "YYYY" for AOP imagery to download example image
ndvi_tresh    <- 0.2           # NDVI treshold used for masking

## ---- Open cleaned and projected SpatialPolygonsDataFrame --------------
if(file.exists(paste0(Dir_data_out,"/",paste0(site_code,"_",aop_data_year,"/"),
                      paste0("poly_sp_proj_", site_code, "_", aop_data_year,".shp")))){
  filename <- paste0(Dir_data_out,"/",paste0(site_code,"_",aop_data_year,"/"),
                     paste0("poly_sp_proj_", site_code, "_", aop_data_year,".shp"))
  poly_sp_proj <- readOGR(filename)
} else {
  filename <- paste0(dir_data_out,"/",paste0(site_code,"_",aop_data_year,"/"),
                     paste0("poly_sp_proj_", site_code, "_", aop_data_year,".shp"))
  poly_sp_proj <- readOGR(filename)
}


## ---- Mosaic and save aop derivative data tiles and create site-level (mosaicked) mask --------------
Dir_data_out_proc <- file.path(Dir_data_out,paste(site_code,aop_data_year,sep="_"),
                               "stacked_aop_data_processed")
aop_data_deriv_filename = file.path(Dir_data_out_proc,
                                    paste0("aop_data_derivatives_mosaic_", 
                                           site_code, "_", aop_data_year, "_",
                                           ndvi_tresh, "ndviTresh.rds"))

if (file.exists(aop_data_deriv_filename)){
  aop_data_derivatives_mosaic <- readRDS(aop_data_deriv_moscaic_filename) 
} else {
  dir_data_out_proc <- file.path(dir_data_out,paste(site_code,aop_data_year,sep="_"),
                                 "stacked_aop_data_processed")
  aop_data_deriv_moscaic_filename = file.path(dir_data_out_proc,paste0("aop_data_derivatives_mosaic_", 
                                                                       site_code, "_", aop_data_year,"_",
                                                                       ndvi_tresh, "ndviTresh.rds"))
  aop_data_derivatives_mosaic <- readRDS(aop_data_deriv_moscaic_filename) 
} else {
  # -- Download all tiles covering all polygons: specify eastings and northings
  eastings <- c(seq(poly_sp_proj@bbox[1,"min"],poly_sp_proj@bbox[1,"max"],by=500),poly_sp_proj@bbox[1,"max"])
  northings <- c(seq(poly_sp_proj@bbox[2,"min"],poly_sp_proj@bbox[2,"max"],by=500),poly_sp_proj@bbox[2,"max"])
  
  # -- Actual download of the tiles
  poly_coordinates <- tidyr::expand_grid(eastings, northings)
  
  ## ---- Generate a list of the 1km x 1km tile coordinates covering this specific polygon or all polygons --------------
  source("01_Download_preprocess_aop_imagery.R")
  tile_coordinates <- list_tiles_covering_poly(poly_coordinates)
  
  ## ---- Specify which tiles (.rds stacked AOP data file) to mosaic --------------
  # We are only opening derivative products here
  # If the tiles were already downloaded and stored previously, access them via the Data folder,
  # if they were downloaded during the current R session, access them via the data folder
  # Attention: this if-loop only makes sense when entire sites were downloaded and stacked at once (see above)
  if (length(list.files(paste0(Dir_data_out, "/", site_code, "_", aop_data_year, "/stacked_aop_data_derivatives"),all.files = TRUE)) != 0){
    Dir_data_out_site_year <- file.path(Dir_data_out,paste(site_code,aop_data_year,sep="_"))
    stacked_aop_data_filename <- file.path(Dir_data_out_site_year,
                                           "stacked_aop_data_derivatives",
                                           paste0("stacked_aop_data_derivatives_",
                                                  tile_coordinates, ".rds"))
  } else {
    dir_data_out_site_year <- file.path(dir_data_out,paste(site_code,aop_data_year,sep="_"))
    stacked_aop_data_filename <- file.path(dir_data_out_site_year,
                                           "stacked_aop_data_derivatives",
                                           paste0("stacked_aop_data_derivatives_",
                                                  tile_coordinates, ".rds"))
  }
  
  ## ---- Load the raster tiles covering the polygon(s) --------------
  filename = stacked_aop_data_filename
  
  aop_data_list <- vector(mode = "list", length = length(filename))
  names(aop_data_list) <- filename
  system.time(
    for (f in 1:length(filename)){
      
      if (file.exists(filename[f])){
        print(paste0("Reading tile No. ", f, " of ", length(filename), " (", tile_coordinates[f], ")"))
        aop_data_list[[f]] <- readRDS(file = filename[f])
      } else {
        print(paste0("tile No. ", f, " of ", length(filename), " (", tile_coordinates[f], ") does not exist"))
        aop_data_list[[f]] <- NA
      }
    })
  
  # Remove list items that are NA
  aop_data_list <- Filter(Negate(anyNA),aop_data_list)
  
  # Specify whether a specific polygon or the entire site (= all polygons of a site) should be masked
  # i = 1
  # poly_sp_proj_set = poly_sp_proj[i,]
  poly_sp_proj_set = poly_sp_proj
  
  ## ---- Create mosaic and mask layer for the entire site --------------
  aop_data_derivatives_mosaic <- derivatives_preprocessing(aop_data_list, poly_sp_proj_set, ndvi_tresh = ndvi_tresh)
  
  # Check whether some layers appear in a temporary folder, and if so, read them in memory
  print("Reading layers in memory - ignore potential error messages")
  aop_data_derivatives_mosaic$aop_data_mosaicked <- readAll(aop_data_derivatives_mosaic$aop_data_mosaicked)
  aop_data_derivatives_mosaic$aop_data_masked <- readAll(aop_data_derivatives_mosaic$aop_data_masked)
  aop_data_derivatives_mosaic$combined_ndvi_chm_mask <- readAll(aop_data_derivatives_mosaic$combined_ndvi_chm_mask)
  
  ## ---- Save the result --------------
  dir_data_out_sitename <- file.path(dir_data_out,paste(site_code,aop_data_year,sep="_"))
  check_create_dir(dir_data_out_sitename)
  dir_data_out_proc <- file.path(dir_data_out_sitename,"stacked_aop_data_processed")
  check_create_dir(dir_data_out_proc)
  aop_data_deriv_filename = file.path(dir_data_out_proc,
                                      paste0("aop_data_derivatives_mosaic_", 
                                             site_code, "_", aop_data_year, "_",
                                             ndvi_tresh, "ndviTresh.rds"))
  saveRDS(aop_data_derivatives_mosaic, file = aop_data_deriv_filename)  
  
}


## ---- Mosaic and save aop ancillary data tiles and check pixel-level weather conditions  --------------
# During inspection, it seemed that all pixels had excellent weather conditions (<10% clouds), 
# so here we just mosaic the ancillary data tiles without running further code to mask pixels based on cloud coverage
Dir_data_out_proc <- file.path(Dir_data_out,paste(site_code,aop_data_year,sep="_"),
                               "stacked_aop_data_processed")
aop_data_deriv_filename = file.path(Dir_data_out_proc,
                                    paste0("aop_data_ancillary_mosaic_", 
                                           site_code, "_", aop_data_year, ".rds"))

if (file.exists(aop_data_deriv_filename)){
  
  print("Mosaic of aop derivative and ancillary tiles already created and saved")
  
} else {
  # -- Download all tiles covering all polygons: specify eastings and northings
  eastings <- c(seq(poly_sp_proj@bbox[1,"min"],poly_sp_proj@bbox[1,"max"],by=500),poly_sp_proj@bbox[1,"max"])
  northings <- c(seq(poly_sp_proj@bbox[2,"min"],poly_sp_proj@bbox[2,"max"],by=500),poly_sp_proj@bbox[2,"max"])
  
  # -- Actual download of the tiles
  poly_coordinates <- tidyr::expand_grid(eastings, northings)
  
  ## ---- Generate a list of the 1km x 1km tile coordinates covering this specific polygon or all polygons --------------
  source("01_Download_preprocess_aop_imagery.R")
  tile_coordinates <- list_tiles_covering_poly(poly_coordinates)
  
  ## ---- Specify which ancillary data tiles (.rds stacked AOP data file) to mosaic --------------
  # We are only opening ancillary products here
  # If the tiles were already downloaded and stored previously, access them via the Data folder,
  # if they were downloaded during the current R sesseion, access them via the data folder
  # Attention: this if-loop only makes sense when entire sites were downloaded and stacked at once (see above)
  if (length(list.files(paste0(Dir_data_out, "/", site_code, "_", aop_data_year, "/stacked_aop_data_ancillary"),all.files = TRUE)) != 0){
    Dir_data_out_site_year <- file.path(Dir_data_out,paste(site_code,aop_data_year,sep="_"))
    stacked_aop_data_filename <- file.path(Dir_data_out_site_year,
                                           "stacked_aop_data_ancillary",
                                           paste0("stacked_aop_data_ancillary_",
                                                  tile_coordinates, ".rds"))
  } else {
    dir_data_out_site_year <- file.path(dir_data_out,paste(site_code,aop_data_year,sep="_"))
    stacked_aop_data_filename <- file.path(dir_data_out_site_year,
                                           "stacked_aop_data_ancillary",
                                           paste0("stacked_aop_data_ancillary_",
                                                  tile_coordinates, ".rds"))
  }
  
  ## ---- Load the ancillary raster tiles covering the polygon(s) --------------
  filename = stacked_aop_data_filename
  
  aop_data_list <- vector(mode = "list", length = length(filename))
  names(aop_data_list) <- filename
  system.time(
    for (f in 1:length(filename)){
      print(paste0("Reading tile No. ", f, " of ", length(filename)))
      if (file.exists(filename[f])){
        aop_data_list[[f]] <- readRDS(file = filename[f])
      } else {
        aop_data_list[[f]] <- character(0)
      }
    })
  # Remove list items that are NA
  aop_data_list <- Filter(Negate(anyNA),aop_data_list)
  
  # Specify whether a specific polygon or the entire site (= all polygons of a site) should be masked
  # i = 1
  # poly_sp_proj_set = poly_sp_proj[i,]
  poly_sp_proj_set = poly_sp_proj

  ## ---- Create mosaic of weather and ancillary data --------------
  aop_data_ancil_mosaic <- ancillary_mosaicking(aop_data_list, poly_sp_proj_set)
  
  # Check whether some layers appear in a temporary folder, and if so, read them in memory
  print("Reading layers in memory - ignore potential error messages")
  aop_data_ancil_mosaic <- readAll(aop_data_ancil_mosaic)
  
  # plot(aop_data_ancil_mosaic$Weather_Quality_Indicator2)
  
  ## ---- Save the result --------------
  dir_data_out_sitename <- file.path(dir_data_out,paste(site_code,aop_data_year,sep="_"))
  check_create_dir(dir_data_out_sitename)
  dir_data_out_proc <- file.path(dir_data_out_sitename,"stacked_aop_data_processed")
  check_create_dir(dir_data_out_proc)
  aop_data_ancil_filename = file.path(dir_data_out_proc,
                                      paste0("aop_data_ancillary_mosaic_",
                                             site_code, "_", aop_data_year, ".rds"))
  saveRDS(aop_data_ancil_mosaic, file = aop_data_ancil_filename)
}

## ---- Spatial visualization for exploratory purposes --------------
# Open plot coordinates
if (file.exists(paste0(Dir_data_out,"/",site_code,"_",aop_data_year,"/",paste0(site_code, "_", aop_data_year,"_veg_plots_spcomp.rda")))){
  filename <- paste0(Dir_data_out,"/",site_code,"_",aop_data_year,"/",paste0(site_code, "_", aop_data_year,"_veg_plots_spcomp.rda"))
  veg_plots_spcomp_site_year <- readRDS(filename)
  folder_coords <- paste0(Dir_data_out,"/",site_code,"_",aop_data_year,"/plot_locations/")
  veg_plots_base_coords <- readRDS(paste0(folder_coords,paste0(site_code,"_veg_plots_base.rds")))
} else {
  filename <- paste0(dir_data_out,"/",site_code,"_",aop_data_year,"/",paste0(site_code, "_", aop_data_year,"_veg_plots_spcomp.rda"))
  veg_plots_spcomp_site_year <- readRDS(filename)
  folder_coords <- paste0(dir_data_out,"/",site_code,"_",aop_data_year,"/plot_locations/")
  veg_plots_base_coords <- readRDS(paste0(folder_coords,paste0(site_code,"_veg_plots_base.rds")))
}

veg_spdf_poly <- veg_plots_base_coords$veg_spdf_point


source("01_Download_preprocess_aop_imagery.R")
plot_deriv_products_plot_loc(aop_data_derivatives_mosaic$aop_data_masked$dtm,
                             "DTM", poly_sp_proj, veg_plots_sub_coords, brewer.pal(n=10,name="BrBG"))  
plot(veg_plots_base_coords$veg_spdf_point, add=T, pch=18, col="magenta")

plot_deriv_products_plot_loc(aop_data_derivatives_mosaic$aop_data_masked$aspect,
                             "Aspect", poly_sp_proj, veg_plots_coords, brewer.pal(n=10,name="RdYlGn"))  
plot_deriv_products_plot_loc(aop_data_derivatives_mosaic$combined_ndvi_chm_mask,
                             "Mask based on NDVI and chm tresholds", poly_sp_proj, veg_plots_coords, brewer.pal(n=10,name="RdYlGn"))  


# plot ancillary data mosaic
Dir_data_out_proc <- file.path(Dir_data_out,paste(site_code,aop_data_year,sep="_"),
                               "stacked_aop_data_processed")
aop_data_deriv_filename = file.path(Dir_data_out_proc,
                                    paste0("aop_data_ancillary_mosaic_",
                                           site_code, "_", aop_data_year, ".rds"))
if (file.exists(aop_data_deriv_filename)){
  
} else {
  dir_data_out_proc <- file.path(dir_data_out,paste(site_code,aop_data_year,sep="_"),
                                 "stacked_aop_data_processed")
  aop_data_deriv_filename = file.path(dir_data_out_proc,
                                      paste0("aop_data_ancillary_mosaic_",
                                             site_code, "_", aop_data_year, ".rds"))
  
}
aop_data_ancil_mosaic <- readRDS(aop_data_deriv_filename)

plot(aop_data_ancil_mosaic$Weather_Quality_Indicator2)

### END of PART 5

##################################################################
## PART 6: Extract plot level information from rasters  ------
##################################################################
# Extract plot level values from raster products
# - from Rangeland Analysis Platform (RAP): NPP
# - from NEON: information on NDVI et.

# USER-DEFINED-INPUTS
# Define parameters for the AOP remote sensing data to download
site_code     <- "ABBY"     # Four-digit NEON site code, character string type
aop_data_year <- "2017"     # Four-digit year in character string "YYYY" for AOP imagery to download example image
ndvi_tresh    <- 0.2        # NDVI treshold used for masking

## ---- Open plot coordinates --------------
if (file.exists(paste0(Dir_data_out,"/",site_code,"_",aop_data_year,"/",paste0(site_code, "_", aop_data_year,"_veg_plots_spcomp.rda")))){
  filename <- paste0(Dir_data_out,"/",site_code,"_",aop_data_year,"/",paste0(site_code, "_", aop_data_year,"_veg_plots_spcomp.rda"))
  veg_plots_spcomp_site_year <- readRDS(filename)
  folder_coords <- paste0(Dir_data_out,"/",site_code,"_",aop_data_year,"/plot_locations/")
  veg_plots_base_coords <- readRDS(paste0(folder_coords,paste0(site_code,"_veg_plots_base.rds")))
} else {
  filename <- paste0(dir_data_out,"/",site_code,"_",aop_data_year,"/",paste0(site_code, "_", aop_data_year,"_veg_plots_spcomp.rda"))
  veg_plots_spcomp_site_year <- readRDS(filename)
  folder_coords <- paste0(dir_data_out,"/",site_code,"_",aop_data_year,"/plot_locations/")
  veg_plots_base_coords <- readRDS(paste0(folder_coords,paste0(site_code,"_veg_plots_base.rds")))
}

## ---- Open cleaned and projected SpatialPolygonsDataFrame --------------
if(file.exists(paste0(Dir_data_out,"/",paste0(site_code,"_",aop_data_year,"/"),
                      paste0("poly_sp_proj_", site_code, "_", aop_data_year,".shp")))){
  filename <- paste0(Dir_data_out,"/",paste0(site_code,"_",aop_data_year,"/"),
                     paste0("poly_sp_proj_", site_code, "_", aop_data_year,".shp"))
  poly_sp_proj <- readOGR(filename)
} else {
  filename <- paste0(main_dir,"/data/data_output/",paste0(site_code,"_",aop_data_year,"/"),
                     paste0("poly_sp_proj_", site_code, "_", aop_data_year,".shp"))
  poly_sp_proj <- readOGR(filename)
}

## ---- Open RAP imagery --------------
setwd(Dir_data_GEE)
files_site <- list.files(pattern = glob2rx(paste0(site_code,"*.tif")), full.names = TRUE)

cover_status_site_year <- raster::brick(files_site[grep(paste0("cover_",aop_data_year), files_site)])
npp_status_site_year <- raster::brick(files_site[grep(paste0("npp_",aop_data_year), files_site)])
status_list <- list("cover_status_site_year" = cover_status_site_year,
                    "npp_status_site_year" = npp_status_site_year)

# Set the working directory back to where the scripts are located
setwd(dir_scripts)

## ---- Open site aop derivative mosaic for masking --------------
Dir_data_out_proc <- file.path(Dir_data_out,paste(site_code,aop_data_year,sep="_"),
                               "stacked_aop_data_processed")
aop_data_deriv_filename = file.path(Dir_data_out_proc,
                                    paste0("aop_data_derivatives_mosaic_", 
                                           site_code, "_", aop_data_year, "_",
                                           ndvi_tresh, "ndviTresh.rds"))

if (file.exists(aop_data_deriv_filename)){
  aop_data_derivatives_mosaic <- readRDS(aop_data_deriv_moscaic_filename) 
} else {
  dir_data_out_proc <- file.path(dir_data_out,paste(site_code,aop_data_year,sep="_"),
                                 "stacked_aop_data_processed")
  aop_data_deriv_moscaic_filename = file.path(dir_data_out_proc,paste0("aop_data_derivatives_mosaic_", 
                                                                       site_code, "_", aop_data_year,"_",
                                                                       ndvi_tresh, "ndviTresh.rds"))
  aop_data_derivatives_mosaic <- readRDS(aop_data_deriv_moscaic_filename) 
}


## ---- Extract and save plot level RAP product values --------------
source("03_Extract_plot_level_raster_values.R")

# Create output folder
Dir_data_out_RAP <- file.path(Dir_data_out,paste(site_code,aop_data_year,sep="_"),"plot_RAP_values")
dir_data_out_siteyear <- file.path(dir_data_out,paste(site_code,aop_data_year,sep="_"))
check_create_dir(dir_data_out_siteyear)
dir_data_out_RAP <- file.path(dir_data_out_siteyear,"plot_RAP_values")
check_create_dir(dir_data_out_RAP)


# 1 yr status
if (file.exists(paste0(Dir_data_out_RAP,"/",paste0(site_code,"_veg_plots_sub_status_", aop_data_year,"_",
                                                   ndvi_tresh, "ndviTresh.rds")))){
  message("The data of ", site_code, " was already downloaded and stored")
} else {
  veg_plots_base_stat <- get_RAP_product_values_plot("stat", status_list, aop_data_derivatives_mosaic,
                                                     veg_plots_base_coords$veg_spdf_poly, ID = "id")
  saveRDS(veg_plots_base_stat, file=paste0(dir_data_out_RAP,  "/",site_code,"_veg_plots_base_status_", aop_data_year, "_", ndvi_tresh, "ndviTresh.rds"))
}


## ---- Visualize RAP extraction--------------
if (!require('RColorBrewer')) install.packages('RColorBrewer'); library('RColorBrewer')
palette <- brewer.pal(n=10,name="RdYlGn")

# Remark: In these visualizations, pixels are not masked based on NDVI and chm treshold
# Masking DID happen prior to feature extraction though!
plot_RAP_products_plot_loc(status_list, veg_plots_base_coords$veg_spdf_poly, "Status 2017", poly_sp_proj, palette)

## ---- Extract and save plot level abiotic NEON product values --------------

Dir_data_out_proc <- file.path(Dir_data_out,paste(site_code,aop_data_year,sep="_"),"stacked_aop_data_processed")
if (file.exists(paste0(Dir_data_out_proc,"/",paste0("veg_plots_sub_neon_", site_code, "_", aop_data_year, "_", 
                                                    ndvi_tresh, "ndviTresh.rds")))){
  message("The data of ", site_code, " was already downloaded and stored")
} else {
  veg_plots_base_NEON <- get_NEON_product_values_plot(aop_data_derivatives_mosaic, veg_plots_base_coords$veg_spdf_poly, ID = "id")

  # Create output folder
  dir_data_out_siteyear <- file.path(dir_data_out,paste(site_code,aop_data_year,sep="_"))
  check_create_dir(dir_data_out_siteyear)
  dir_data_out_proc <- file.path(dir_data_out_siteyear,"stacked_aop_data_processed")
  check_create_dir(dir_data_out_proc)
  setwd(dir_data_out_proc)
  
  saveRDS(veg_plots_base_NEON, file=paste0(dir_data_out_proc, "/", "veg_plots_base_neon_", site_code, "_", aop_data_year,"_", ndvi_tresh, "ndviTresh_cropped.rds"))
}


### END of PART 6


#####################################################################################################
## PART 7: Calculate spectral diversity for base plots ------
#####################################################################################################
# This section involves the following steps
# 1. Clip plots from raster tiles containing hyperspectral stack of reflectance images
# 2. Remove unwanted pixels from clipped NEON tiles: based on height and NDVI treshold (use the same tresholds as before)
# 4. Spectral preprocessing steps per plot: brightness normalization, PCA

# USER-DEFINED-INPUTS
# Define parameters for the AOP remote sensing data
site_code     <- "ABBY"     # Four-digit NEON site code, character string type
aop_data_year <- "2017"     # Four-digit year in character string "YYYY" for AOP imagery to download example image
ndvi_tresh    <- 0.2        # NDVI treshold used for masking
parallel      <- "no"       # Indicate whether you want to process plots in parallel or not. If yes: Specicy cores a bit further down

## ---- Open plot coordinates --------------
if (file.exists(paste0(Dir_data_out,"/",site_code,"_",aop_data_year,"/plot_locations/",paste0(site_code,"_veg_plots_base.rds")))){
  folder_coords <- paste0(Dir_data_out,"/",site_code,"_",aop_data_year,"/plot_locations/")
  veg_plots_base_coords <- readRDS(paste0(folder_coords,paste0(site_code,"_veg_plots_base.rds")))
  veg_plots_base_spdf_poly <- veg_plots_base_coords$veg_spdf_poly
} else {
  folder_coords <- paste0(dir_data_out,"/",site_code,"_",aop_data_year,"/plot_locations/")
  veg_plots_base_coords <- readRDS(paste0(folder_coords,paste0(site_code,"_veg_plots_base.rds")))
  veg_plots_base_spdf_poly <- veg_plots_base_coords$veg_spdf_poly
}

## ---- Specify which plots to process --------------
# Option 1: Just 1 plot
# i = 1
# plotname <- veg_plots_base_spdf_poly$id[i]
# veg_plots_base_set <- veg_plots_base_spdf_poly[i,]

# Option 2: Multiple plots
# i <- c(1,2,3)
# veg_plots_base_set <- veg_plots_base_spdf_poly[i,]

# Option 3: Process all plots of a site
veg_plots_base_set = veg_plots_base_spdf_poly

## ---- Optional: Specify how many cores to use in parallel processing --------------
if(parallel == "yes") {
  if (!require('doParallel')) install.packages('doParallel'); library('doParallel')
  if (!require('doSNOW')) install.packages('doSNOW'); library('doSNOW')
  
  # detectCores()
  # # if (detectCores() <= length(filename)){
  # #   cores_nb <- detectCores() -1 #use all but one core
  # # } else {
  # #   cores_nb <- length(filename)
  # # }
  cores_nb <- 8
}


## ---- Pre-processing data per plot and save results --------------
baseplot_filename = file.path(Dir_data_out,paste(site_code,aop_data_year,sep="_"),
                              "stacked_aop_data_processed",
                              paste0("aop_data_hs_veg_plots_base_", site_code, "_", aop_data_year, "_",
                                     ndvi_tresh,"ndviTresh.rds"))
if (file.exists(baseplot_filename)){
  aop_veg_plots_base <- readRDS(baseplot_filename)
  
} else if (file.exists(file.path(dir_data_out,paste(site_code,aop_data_year,sep="_"),
                                 "stacked_aop_data_processed",
                                 paste0("aop_data_hs_veg_plots_base_", site_code, "_", aop_data_year, "_",
                                        ndvi_tresh,"ndviTresh.rds")))) {
  aop_veg_plots_base <- readRDS(file.path(dir_data_out,paste(site_code,aop_data_year,sep="_"),
                                          "stacked_aop_data_processed",
                                          paste0("aop_data_hs_veg_plots_base_", site_code, "_", aop_data_year, "_",
                                                 ndvi_tresh,"ndviTresh.rds")))
  
} else {
  source("04_Calculate_spectral_diversity.R")
  # Perform the following steps per plot
  #  - Crop the hyperspectral rasterstack (tile or mosaic of tiles) to the extent of the plot
  #  - mask pixels using NDVI and chm tresholds
  #  - perform PCA
  
  if(parallel == "yes"){
    aop_veg_plots_base <- get_preprocessed_plots(veg_plots_base_set, ndvi_tresh=ndvi_tresh, scaling_nb=scaling_nb, 
                                                 site_code, aop_data_year, Dir_data_out, dir_data_out)
  } else {
    aop_veg_plots_base <- get_preprocessed_plots_parallel(veg_plots_base_set, ndvi_tresh=ndvi_tresh, scaling_nb=scaling_nb, 
                                                          site_code, aop_data_year, Dir_data_out, dir_data_out, cores_nb)
  } 
  
  # Save the results
  dir_data_site_year <- file.path(dir_data_out,paste(site_code,aop_data_year,sep="_"))
  check_create_dir(dir_data_site_year)
  dir_data_out_proc <- file.path(dir_data_out,paste(site_code,aop_data_year,sep="_"),"stacked_aop_data_processed")
  check_create_dir(dir_data_out_proc)
  aop_veg_plots_base_filename = file.path(dir_data_out_proc,paste0("aop_data_hs_veg_plots_base_",site_code, "_", aop_data_year, "_", 
                                                                   ndvi_tresh,"ndviTresh.rds"))
  
  saveRDS(aop_veg_plots_base, file = aop_veg_plots_base_filename)  
  
}

## ---- Calculate and save spectral diversity per plot: CV, CVH and SAM --------------
source("04_Calculate_spectral_diversity.R")

### 1) Coefficient of variation = standard deviation/mean
aop_veg_plots_base_cv <- get_cv_plot(aop_veg_plots_base)

### 2) Convex hull volume
aop_veg_plots_base_chv <- get_chv_plot(aop_veg_plots_base)

### 3) Mean spectral Angle Measure 
aop_veg_plots_base_sam <- get_sam_plot(aop_veg_plots_base)

### Concatenate different diversity metrics and save
library(dplyr)
aop_veg_plots_base_sdiv <- aop_veg_plots_base_cv %>%
  inner_join(aop_veg_plots_base_chv, by=c("id")) %>%
  inner_join(aop_veg_plots_base_sam, by=c("id")) 

dir_data_site_year <- file.path(dir_data_out,paste(site_code,aop_data_year,sep="_"))
check_create_dir(dir_data_site_year)
dir_data_out_proc <- file.path(dir_data_out,paste(site_code,aop_data_year,sep="_"),"stacked_aop_data_processed")
check_create_dir(dir_data_out_proc)

aop_veg_plots_base_sdiv_filename <- file.path(dir_data_out_proc, 
                                              paste0("veg_plots_base_sdiv_", site_code,"_",aop_data_year,"_",
                                                     ndvi_tresh,"ndviTresh.rds"))
saveRDS(aop_veg_plots_base_sdiv, aop_veg_plots_base_sdiv_filename)

### END of PART 7

############################################################
## PART 8: Merge data from different data sources ------------
############################################################
# USER-DEFINED-INPUTS
site_code     <- "ABBY"     # Four-digit NEON site code, character string type
aop_data_year <- "2017"     # Four-digit year in character string "YYYY" for AOP imagery to download example image
ndvi_tresh    <- 0.2        # NDVI treshold used for masking  

if (!require('magrittr')) install.packages('magrittr'); library('magrittr')

## ---- Open plot coordinates for NEON baseplots --------------
if (file.exists(paste0(Dir_data_out,"/",site_code,"_",aop_data_year,"/plot_locations/",paste0(site_code,"_veg_plots_base.rds")))){
  folder_coords <- paste0(Dir_data_out,"/",site_code,"_",aop_data_year,"/plot_locations/")
  veg_plots_base_coords <- readRDS(paste0(folder_coords,paste0(site_code,"_veg_plots_base.rds")))
} else {
  folder_coords <- paste0(dir_data_out,"/",site_code,"_",aop_data_year,"/plot_locations/")
  veg_plots_base_coords <- readRDS(paste0(folder_coords,paste0(site_code,"_veg_plots_base.rds")))
}

## ---- Open RAP product values for NEON baseplots --------------
Dir_data_out_RAP <- file.path(Dir_data_out,paste(site_code,aop_data_year,sep="_"),"plot_RAP_values")
dir_data_out_siteyear <- file.path(dir_data_out,paste(site_code,aop_data_year,sep="_"))
dir_data_out_RAP <- file.path(dir_data_out_siteyear,"plot_RAP_values")

# 1 yr status
if (file.exists(paste0(Dir_data_out_RAP,"/",paste0(site_code,"_veg_plots_sub_status_", aop_data_year,"_", 
                                                   ndvi_tresh,"ndviTresh.rds")))){
  veg_plots_base_stat <- readRDS(paste0(Dir_data_out_RAP,"/",site_code,"_veg_plots_base_status_", aop_data_year,"_", ndvi_tresh,"ndviTresh.rds"))
} else {
  veg_plots_base_stat <- readRDS(paste0(dir_data_out_RAP,"/",site_code,"_veg_plots_base_status_", aop_data_year,"_", ndvi_tresh,"ndviTresh.rds"))
}

## ---- Open abiotic NEON product values for NEON baseplots --------------
Dir_data_out_proc <- file.path(Dir_data_out,paste(site_code,aop_data_year,sep="_"),"stacked_aop_data_processed")

if (file.exists(paste0(Dir_data_out_proc,"/",paste0("veg_plots_base_neon_", site_code, "_", aop_data_year, "_",
                                                    ndvi_tresh,"ndviTresh.rds")))){
  veg_plots_base_NEON <- readRDS(paste0(Dir_data_out_proc,"/","veg_plots_base_neon_", site_code, "_", aop_data_year,"_", ndvi_tresh,"ndviTresh.rds"))
} else {
  dir_data_out_siteyear <- file.path(dir_data_out,paste(site_code,aop_data_year,sep="_"))
  dir_data_out_proc <- file.path(dir_data_out_siteyear,"stacked_aop_data_processed")
  
  veg_plots_base_NEON <- readRDS(paste0(dir_data_out_proc,"/","veg_plots_base_neon_", site_code, "_", aop_data_year,"_", ndvi_tresh,"ndviTresh.rds"))
}


## ---- Open spectral diversity values for NEON baseplots -------------
Dir_data_out_proc <- file.path(Dir_data_out,paste(site_code,aop_data_year,sep="_"),"stacked_aop_data_processed")

if (file.exists(paste0(Dir_data_out_proc,"/",paste0("veg_plots_base_sdiv_",site_code, "_", aop_data_year,"_",ndvi_tresh,"ndviTresh.rds")))){
  aop_veg_plots_base_sdiv <- readRDS(file.path(Dir_data_out_proc,paste0("veg_plots_base_sdiv_",site_code, "_", aop_data_year,"_",ndvi_tresh,"ndviTresh.rds")))
  # aop_veg_plots_base_sdiv <- tibble::rownames_to_column(aop_veg_plots_base_sdiv, "id")
} else {
  dir_data_out_proc <- file.path(dir_data_out,paste(site_code,aop_data_year,sep="_"),"stacked_aop_data_processed")
  aop_veg_plots_base_sdiv <- readRDS(file.path(dir_data_out_proc,paste0("veg_plots_base_sdiv_",site_code, "_", aop_data_year,"_",ndvi_tresh,"ndviTresh.rds")))
  # aop_veg_plots_base_sdiv <- tibble::rownames_to_column(aop_veg_plots_base_sdiv, "id")
}

## ---- Concatenate all information for NEON baseplots -------------
veg_plots_base <- veg_plots_base_coords$veg_utm %>%
  dplyr::select(id, namedLocation, domainID, siteID, plotID, nlcdClass,
                decimalLatitude, decimalLongitude, elevation, easting, northing, 
                utmHemisphere, utmZoneNumber, utmZone, 
                coordinateSource, namedLocationCoordUncertainty, namedLocationElevUncertainty,
                geodeticDatum, slopeAspect, slopeGradient, soilTypeOrder) %>%
  mutate(id = plotID) %>%
  inner_join(aop_veg_plots_base_sdiv, by = c("id")) %>%
  inner_join(veg_plots_base_stat, by = c("id")) %>%
  inner_join(veg_plots_base_NEON, by = c("id"))

## ---- Save the concatenated dataframes -------------
dir_data_site_year <- file.path(dir_data_out,paste(site_code,aop_data_year,sep="_"))
check_create_dir(dir_data_site_year)
dir_data_out_proc <- file.path(dir_data_site_year,"stacked_aop_data_processed")
check_create_dir(dir_data_out_proc)

baseplot_filename = file.path(dir_data_out_proc,paste0("veg_plots_base_all_", 
                                                       site_code, "_", aop_data_year, "_",
                                                       ndvi_tresh,"ndviTresh.rds"))
saveRDS(veg_plots_base, file = baseplot_filename)  

### END of PART 8

##################################################
