###################################################################
###                                                             ###
###   Download, preprocess and visualize NEON airborne imagery  ###
###                                                             ###
###################################################################


# This script contains functions to download, stack, preporocess and visualize NEON AOP remote sensing imagery
# Before stacking, hyperspectral signatures are preprocessed as follows: 
# Removal of bad bands, application Savitzky-Golay filter

# Parts of this script are based on scripts provided on www.neonscience.org, written by Victoria Scholl and Etienne Laliberté

# by Elisa Van Cleemput, 2023
#######################################################
# Function to download 1 NEON hyperspectral tile of the site and year under investigation to extract the crs from
# Script based on: https://www.neonscience.org/resources/learning-hub/tutorials/neon-api-usage

# Used in PART 1

if (!require('neonUtilities')) install.packages('neonUtilities'); library('neonUtilities')
if (!require('httr')) install.packages('httr'); library('httr')
if (!require('jsonlite')) install.packages('jsonlite'); library('jsonlite')
if (!require('downloader')) install.packages('downloader'); library('downloader')
if (!require('rhdf5')) install.packages('rhdf5'); library('rhdf5')
if (!require('magrittr')) install.packages('magrittr'); library('magrittr')

# Ir rhdf5 does not install with the previous line
# if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
# BiocManager::install("rhdf5")

download_neon_refl_crsexample <- function(Dir_out,dir_out,site_code,aop_data_year){
  
  # check if an example tile was already downloaded for the current site and year
  if (file.exists(paste0(Dir_out, "/", site_code, aop_data_year, "_exampleforcrs.h5"))){
    # if it exists, read that instead of re-downloading the image.
    message("Example tile already downloaded for current site and year")
    f <- paste0(Dir_out, "/", site_code, aop_data_year, "_exampleforcrs.h5")
    return(f)
    
  } else {
    # Download 1 tile
    # Request reflectance mosaic data availability info
    req.aop <- GET("http://data.neonscience.org/api/v0/products/DP3.30006.001")
    # make this JSON readable
    avail.aop <- fromJSON(content(req.aop, as="text"),
                          simplifyDataFrame=T, flatten=T)
    # get data availability list for the product
    avi.urls <- unlist(avail.aop$data$siteCodes$availableDataUrls)
    
    # get data availability from location/date of interest
    avi_site_year <- GET(avi.urls[intersect(grep(site_code, avi.urls),
                                            grep(aop_data_year, avi.urls))])
    avi.files <- fromJSON(content(avi_site_year, as="text"))
    
    # # Optional: Check the first/last ten files
    # head(avi.files$data$files$name, 10)
    # tail(avi.files$data$files$name, 10)
    # head(avi.files$data$files$size, 10)
    # head(avi.files$data$files$url, 10)
    
    # define filepath to the hyperspectral dataset
    f <- paste0(dir_out, "/", site_code, aop_data_year, "_exampleforcrs.h5")
    download(avi.files$data$files$url[1],f,mode="wb")
    
    return(f)
  }
}

#######################################################
# Functions to download all NEON AOP remote sensing tiles in a dataframe with tile coordinates 

# Used in PART 2

# Each data product is downloaded to a folder named by its data product ID. 
# This script also moves the downloaded files into folders that are easier accessible 

# This script is an adapted version of the code written by Victoria Scholl:
# https://github.com/earthlab/neon-veg.git: 00-supporting-functions.R and 04-download_aop_imagery.R

# Function to summarize list of coordinates to a tibble list of the 1km x 1km tile coordinates covering them
list_tiles_covering_poly_tibble <- function(veg_df){
  
  print("Generating list of tiles covering polygon(s)...")
  
  # get easting, northing coordinates
  e <- NA
  n <- NA
  for (i in 1:dim(veg_df)[1]){
    easting <- floor(as.numeric(veg_df$eastings[i])/1000)*1000
    northing <- floor(as.numeric(veg_df$northings[i])/1000)*1000
    e[i] <- easting
    n[i] <- northing
  }
  
  # find unique rows repesenting tiles
  easting_northings <- data.frame(as.character(e),as.character(n))
  colnames(easting_northings) <- c('eastings','northings')
  tiles <- unique(easting_northings[,c('eastings','northings')])
  
  tiles_tibble <- as_tibble(tiles)
  
  options(scipen = 100, digits = 6)
  tiles_tibble$eastings <- as.numeric(tiles_tibble$eastings)
  tiles_tibble$northings <- as.numeric(tiles_tibble$northings)
  
  return(tiles_tibble)
}

# Functions to download and move NEON imagery to a temporary map before stacking
move_downloaded_files <- function(dir_out, dp_id, dp_name, file_pattern, delete_orig = FALSE, unzip = FALSE){
  # moves AOP files downloaded using the neonUtilities::byTileAOP
  # function into a folder with an intuitive name for each 
  # remote sensing data type. 
  # 
  # Args
  #   dir_out - character string path to the main directory where
  #             files were downloaded.
  #             example: "data/data_raw
  #
  #   dp_id - character string with the data product ID 
  #          example: "DP3.30015.001" for the Canopy Height Model
  #
  #   dp_name - character string with a short name describing the data 
  #            this is the folder name where the files will be moved. 
  #            this folder will be within the dir_out folder. 
  #            example: "chm" for Canopy Height Model 
  # 
  #   file_pattern - character string used to identify downloaded files 
  #                  that will be moved by this function. 
  #                  example: "*CHM.tif$" for Canopy Height Model tiles
  #                           that include CHM and end with a .tif extension
  # 
  #   delete_orig - optional logical parameter to delete the original 
  #                 downloaded files. Default value of FALSE 
  # 
  #   unzip - optional logical parameter to unxip the downloaded files
  #           after they have been moved. Deafult value of FALSE. 
  # 
  # Example function call: 
  # move_downloaded_files(dir_out = "data/data_raw/NIWO_2017"
  #                       ,dp_id = "DP3.30026.001"
  #                       ,dp_name = "veg_indices"
  #                       ,file_pattern = "*VegIndices.zip$"
  #                       ,delete_orig = TRUE
  #                       ,unzip = TRUE)
  
  
  # list all filenames 
  download_path <- file.path(dir_out, dp_id)
  list_files <- list.files(path = download_path
                           ,pattern = file_pattern
                           ,recursive = TRUE
                           ,full.names = TRUE)
  
  if (length(list_files) == 0) {
    print( "No tiles were downloaded")
    
  } else {
    # move files into a new folder with an intuitive name
    move_dir <- file.path(dir_out, dp_name)
    check_create_dir(move_dir)
    files_from <- list_files
    files_to <- paste(move_dir
                      ,sapply(stringr::str_split(list_files
                                                 ,.Platform$file.sep)
                              ,tail, 1)
                      ,sep = .Platform$file.sep)
    # copy the downloaded files into the simplified directory
    file.copy(from = files_from
              ,to = files_to
              ,overwrite = TRUE)
    
    # If the copy was successful, delete original files in the nested directories
    if(delete_orig){
      if(all(file.exists(files_to))){
        unlink(download_path
               ,recursive = TRUE)
      }
    }
    
    # Unzip compressed files 
    if(unzip){
      for(zipFile in files_to){
        utils::unzip(zipFile, exdir = move_dir)
      }
    }
    
    print(paste("Downloaded files moved from", download_path, "to",move_dir))
  }
}

download_neon_aop <- function(dp_chm,dp_aspect_slope,dp_elevation, 
                              dp_veg_indices,dp_lai,dp_fpar,
                              dp_hs_refl,
                              site_code,aop_data_year,
                              veg_coordinates,
                              buffer_val,
                              dir_data_raw,
                              check_size=T,delete_originals=T){
  # CHM ---------------------------------------------------------------------
  if (!is.na(dp_chm)){
    # dp_chm <- "DP3.30015.001"
    
    print("Downloading chm")
    
    neonUtilities::byTileAOP(
      dpID = dp_chm
      ,site = site_code
      ,year = aop_data_year
      ,savepath = dir_data_raw
      ,easting = veg_coordinates$eastings
      ,northing = veg_coordinates$northings
      ,check.size = check_size
      ,buffer = buffer_val)
    
    # move CHM files
    move_downloaded_files(dir_out = dir_data_raw, dp_id = dp_chm
                          ,dp_name = "chm", file_pattern = "*CHM.tif$"
                          ,delete_orig = delete_originals)
  }
  
  # aspect_slope ------------------------------------------------------------
  if(!is.na(dp_aspect_slope)){
    # dp_aspect_slope <- "DP3.30025.001"
    
    print("Downloading aspect and slope")
    
    neonUtilities::byTileAOP(
      dpID = dp_aspect_slope
      ,site = site_code
      ,year = aop_data_year
      ,savepath = dir_data_raw
      ,easting = veg_coordinates$eastings
      ,northing = veg_coordinates$northings
      ,check.size = check_size
      ,buffer = buffer_val)
    
    # move aspect files
    move_downloaded_files(dir_out = dir_data_raw, dp_id = dp_aspect_slope
                          ,dp_name = "aspect", file_pattern = "*aspect.tif$"
                          ,delete_orig = delete_originals)
    # move slope files; delete originally downloaded aspect and slope files 
    move_downloaded_files(dir_out = dir_data_raw, dp_id = dp_aspect_slope
                          ,dp_name = "slope", file_pattern = "*slope.tif$"
                          ,delete_orig = delete_originals)
  }
  
  # elevation ---------------------------------------------------------------------
  if(!is.na(dp_elevation)){
    # dp_elevation <- "DP3.30024.001" 
    
    print("Downloading elevation (DSM and DTM)")
    
    neonUtilities::byTileAOP(
      dpID = dp_elevation
      ,site = site_code
      ,year = aop_data_year
      ,savepath = dir_data_raw
      ,easting = veg_coordinates$eastings
      ,northing = veg_coordinates$northings
      ,check.size = check_size
      ,buffer = buffer_val)
    
    # move DSM files
    move_downloaded_files(dir_out = dir_data_raw, dp_id = dp_elevation
                          ,dp_name = "dsm", file_pattern = "*DSM.tif$"
                          ,delete_orig = delete_originals)
    # move DTM files; delete originally downloaded DSM and DTM files 
    move_downloaded_files(dir_out = dir_data_raw, dp_id = dp_elevation
                          ,dp_name = "dtm", file_pattern = "*DTM.tif$"
                          ,delete_orig = delete_originals)
    
  }  
  
  # Vegetation indices ------------------------------------------------------
  if(!is.na(dp_veg_indices)){
    # dp_veg_indices <- "DP3.30026.001"
    
    print("Downloading vegetation indices")
    
    neonUtilities::byTileAOP(
      dpID = dp_veg_indices
      ,site = site_code
      ,year = aop_data_year
      ,savepath = dir_data_raw
      ,easting = veg_coordinates$eastings
      ,northing = veg_coordinates$northings
      ,check.size = check_size
      ,buffer = buffer_val)
    
    move_downloaded_files(dir_out = dir_data_raw, dp_id = dp_veg_indices
                          ,dp_name = "veg_indices"
                          ,file_pattern = "*VegIndices.zip$"
                          ,delete_orig = delete_originals, unzip = TRUE)
  }
  
  
  # Leaf area index ------------------------------------------------------
  if(!is.na(dp_lai)){
    # dp_lai <- "DP3.30012.001"
    
    print("Downloading LAI")
    
    neonUtilities::byTileAOP(
      dpID = dp_lai
      ,site = site_code
      ,year = aop_data_year
      ,savepath = dir_data_raw
      ,easting = veg_coordinates$eastings
      ,northing = veg_coordinates$northings
      ,check.size = check_size
      ,buffer = buffer_val)
    
    move_downloaded_files(dir_out = dir_data_raw, dp_id = dp_lai
                          ,dp_name = "lai",file_pattern = "*LAI.tif$"
                          ,delete_orig = delete_originals)
  }
  
  # fPAR ------------------------------------------------------
  if(!is.na(dp_fpar)){
    # dp_fpar <- "DP3.30014.001"
    
    print("Downloading fPAR")
    
    neonUtilities::byTileAOP(
      dpID = dp_fpar
      ,site = site_code
      ,year = aop_data_year
      ,savepath = dir_data_raw
      ,easting = veg_coordinates$eastings
      ,northing = veg_coordinates$northings
      ,check.size = check_size
      ,buffer = buffer_val)
    
    move_downloaded_files(dir_out = dir_data_raw, dp_id = dp_fpar
                          ,dp_name = "fpar",file_pattern = "*fPAR.tif$"
                          ,delete_orig = delete_originals)
  }
  
  # Hyperspectral reflectance -----------------------------------------------
  if (!is.na(dp_hs_refl)) {
    # dp_hs_refl <- "DP3.30006.001"
    
    print("Downloading hyperspectral images")
    
    neonUtilities::byTileAOP(
      dpID = dp_hs_refl
      ,site = site_code
      ,year = aop_data_year
      ,savepath = dir_data_raw
      ,easting = veg_coordinates$eastings
      ,northing = veg_coordinates$northings
      ,check.size = check_size
      ,buffer = buffer_val)
    
    
    move_downloaded_files(dir_out = dir_data_raw, dp_id = dp_hs_refl
                          ,dp_name = "hyperspectral"
                          ,file_pattern = "*reflectance.h5$"
                          ,delete_orig = delete_originals)
  }
}  

#######################################################
# Functions for hyperspectral stacking and pre-processing

# Used in PART 2

# This script is an adapted version of the code written by Victoria Scholl:
# https://github.com/earthlab/neon-veg.git: 00-supporting-functions.R and 05-prep_aop_imagery.R

# and Laliberté et al. 2020:
# https://github.com/elaliberte/specdiv/blob/master/bartlett/scripts/bartlett.R
# --> masking and smoothing code

if (!require('hsdar')) install.packages('hsdar'); library('hsdar')
if (!require('stringr')) install.packages('stringr'); library('stringr')

stack_hyperspectral <- function(h5){
  # This function creates a rasterstack object for the specified HDF5 
  # filename. 
  #
  # Args: 
  # h5
  #   character string filename of HDF5 file 
  # out_dir
  #   directory for output files where the wavelengths will be written to
  #   a text file for further analysis
  #
  # Returns: 
  # s
  #   RasterStack (collection of Raster layers with the same spatial extent
  #   and resolution) containing nrows x ncols x nbands, based on the 
  #   resolution and number of bands in the HDF5 file. This RasterStack
  #   can then be clipped using Spatial vector layers. 
  
  # list the contents of HDF5 file
  h5_struct <- rhdf5::h5ls(h5, all=T)
  
  # construct the string using "/Reflectance/Metadata/Coordinate_System",
  # without explicitly using a site code 
  crs_tag <- h5_struct$group[grepl("/Reflectance/Metadata/Coordinate_System", 
                                   h5_struct$group)][1] 
  
  # read coordinate reference system data
  crs_info <- rhdf5::h5read(h5, crs_tag)
  
  # convert "UTM" to lowercase "utm" for proper usage later
  crs_info$Proj4 <- CRS(chartr("UTM", "utm", crs_info$Proj4))
  
  # get attributes for the Reflectance dataset.
  # construct the string using "/Reflectance/Reflectance_Data"" 
  refl_tag <- paste0(h5_struct$group[grepl("/Reflectance", 
                                           h5_struct$group)][1],
                     "/Reflectance_Data")
  
  # read the reflectance metadata
  refl_info <- rhdf5::h5readAttributes(h5,refl_tag)
  
  # get the dimensions of the reflectance data
  n_rows <- refl_info$Dimensions[1]
  n_cols <- refl_info$Dimensions[2]
  n_bands <- refl_info$Dimensions[3]
  
  # print dimensions 
  print(paste0("# Rows: ", as.character(n_rows)))
  print(paste0("# Columns: ", as.character(n_cols)))
  print(paste0("# Bands: ", as.character(n_bands)))
  
  # read the wavelengths of the hyperspectral image bands
  wavelength_tag <- paste0(h5_struct$group[grepl("/Reflectance/Metadata/Spectral_Data", 
                                                 h5_struct$group)][1],
                           "/Wavelength")
  wavelengths <- rhdf5::h5read(h5,
                               wavelength_tag)
  
  # define spatial extent: extract resolution and origin coordinates
  map_info <- unlist(strsplit(crs_info$Map_Info, 
                              split = ", "))
  res_x <- as.numeric(map_info[6])
  res_y <- as.numeric(map_info[7])
  x_min <- as.numeric(map_info[4])
  y_max <- as.numeric(map_info[5])
  
  # calculate the maximum X and minimum Y values 
  x_max <- (x_min + (n_cols * res_x))
  y_min <- (y_max - (n_rows * res_y))
  tile_extent <- raster::extent(x_min, x_max, y_min, y_max)
  print("tile extent")
  print(tile_extent)
  
  # read reflectance data for all bands
  refl <- rhdf5::h5read(h5, refl_tag,
                        index = list(1:n_bands, 1:n_cols, 1:n_rows))
  
  # view and apply scale factor to convert integer values to reflectance [0,1]
  # and data ignore value
  scale_factor <- refl_info$Scale_Factor
  data_ignore <- refl_info$Data_Ignore_Value
  refl[refl == data_ignore] <- NA 
  refl_scaled <- refl / scale_factor
  
  # create georeferenced raster using band 1 
  r1 <- (refl_scaled[1,,]) # convert first band to matrix
  # transpose the image pixels for proper orientation to match
  # the other layers. create a raster for this band and assign
  # the CRS.
  print("Transposing reflectance data for proper orientation")
  r1 <- raster::t(raster::raster(r1, crs = crs_info$Proj4))
  extent(r1) <- tile_extent
  
  # start the raster stack with first band 
  s <- raster::stack(r1)
  
  # loop through bands and create a giant rasterstack with 426 (n_bands) bands
  print("Stacking hyperspectral bands")
  for(b in 2:n_bands){
    # print(b)
    
    # create raster with current band
    r <- (refl_scaled[b,,]) # convert to matrix
    r <- raster::t(raster::raster(r, crs = crs_info$Proj4))
    extent(r) <- tile_extent
    
    # add additional band to the stack with the addLayer function
    s <- raster::addLayer(s, r)
    
  }
  
  # adjust the names for each layer in raster stack to correspond to wavelength
  names(s) <- round(wavelengths)
  
  # return the stacked hyperspectral data to clip with vector files 
  hs_stack <- list(s = s,
                   wavelengths = wavelengths)
  return(hs_stack)
}

mask_sg_hyperspectral <- function(hs_stack, out_dir, n_sg){
  cube_stack <- hs_stack$s
  wavelengths <- as.numeric(hs_stack$wavelengths)
  
  print("Creating hyperspectral RasterBrick and Speclib")
  cube_brick <- raster::brick(cube_stack)
  cube_speclib <- speclib(spectra = cube_brick, wavelength = wavelengths)
  # plot(cube_speclib)
  
  print("Masking bad bands")
  band_mask <- data.frame(lb = c(350, 1340, 1790, 2400),
                          ub = c(400, 1445, 1955, 2550))
  hsdar::mask(cube_speclib) <- band_mask
  # plot(cube_speclib)
  
  print("Applying Savitzky-Golay filter (Warning: takes some minutes)")
  system.time(cube_sg <- noiseFiltering(cube_speclib, method = 'sgolay', n = n_sg))
  # plot(cube_sg)
  
  # Transform to RasterBrick
  cube_sg_brick <- cube_sg@spectra@spectra_ra
  
  # Make sure to read and save the entire content, not just links to the Temp folder
  cube_sg_brick <- readAll(cube_sg_brick)
  
  # Identify bad bands that were removed
  wvls <- wavelengths
  window1 <- c(1340, 1445)
  window2 <- c(1790, 1955)
  window1_wvls <- which(wvls > window1[1] & wvls < window1[2])
  window2_wvls <- which(wvls > window2[1] & wvls < window2[2])
  below_400 <- which(wvls < 400)
  above_2400 <- which(wvls > 2400)
  
  # Get bad bands together
  bad_bands <- c(below_400, window1_wvls, window2_wvls, above_2400)
  
  wvl <- cbind(as.data.frame(wavelengths), band = seq(1,length(wavelengths)))
  wvl_bands <- wvl %>%
    dplyr::filter(!band %in% bad_bands)
  
  # Add band names to cube_sg_brick
  names(cube_sg_brick) <- round(wvl_bands$wavelengths)
  
  # write the exact wavelengths to a text file for future use
  print("Save wavelengths vector")
  write.table(data.frame(wavelengths = wvl_bands$wavelengths),
              file.path(out_dir,"wavelengths.txt"),
              row.names=FALSE)
  
  return(cube_sg_brick)
}

stack_neon_aop <- function(veg_coordinates, dir_data_raw, dir_data_out, n_sg){
  
  # Specify the paths for each data directory
  h5_dir <- file.path(dir_data_raw, "hyperspectral")
  chm_dir <- file.path(dir_data_raw, "chm")
  veg_indices_dir <- file.path(dir_data_raw, "veg_indices")
  
  # output directory for stacked AOP data
  stacked_aop_data_dir <- file.path(dir_data_out, "stacked_aop_data")
  check_create_dir(stacked_aop_data_dir)
  
  # list the files in each data directory; filter results based on file type
  # hyperspectral data - list the .h5 files 
  h5_list <- list.files(path = h5_dir, full.names = TRUE) 
  chm_list <- list.files(path = chm_dir, full.names = TRUE)
  veg_indices_list <- list.files(path = veg_indices_dir, full.names = TRUE)
  # remove any entries from the list that are .zip files
  veg_indices_list<- veg_indices_list[grepl("*.tif$", veg_indices_list)]
  veg_indices_names <- c("ARVI","EVI","NDLI","NDNI","NDVI","PRI","SAVI")
  
  # create data cubes with AOP-derived features 
  start_time <- Sys.time() # start the timer 
  
  # loop through the tiles; build up a data cube for each tile 
  # for (h5 in h5_list) {
  for (i in 1:nrow(veg_coordinates)) {
    
    # each current hyperspectral tile must be read and stacked into a 
    # georeferenced rasterstack object (so it can be clipped with point / polygon
    # shapefiles). The process of creating a rasterstack takes a while for 
    # each tile, so after creating each rasterstack once, each object gets 
    # written to a file. 
    
    # Build up the rasterstack filename by combining the easting/northing coordinates,
    # and use this to find corresponding tiles of various remote sensing data
    easting <- veg_coordinates[i,"eastings"]
    northing <- veg_coordinates[i,"northings"]
    east_north_string <- paste(easting, northing, sep="_")

    # generate a filename for the stacked AOP data
    stacked_aop_data_filename = file.path(stacked_aop_data_dir,
                                          paste0("stacked_aop_data_",
                                                 east_north_string, ".rds"))
    
    # Before proceeding, first check if the polygon coordinates lie within the flight area
    h5 <- grep(east_north_string, h5_list, value=TRUE)
    chm <- grep(east_north_string, chm_list, value=TRUE)
    veg_indices <- grep(east_north_string, veg_indices_list, value=TRUE)
    
    if(length(h5) == 0 | length(chm) == 0 | length(veg_indices) == 0){
      print(paste0("No (matching) tiles available for coordinates ", east_north_string))
    } else {
      
      # check if a .rds file already exists for the current feature data cube
      if (file.exists(stacked_aop_data_filename)){
        # if it exists, read that instead of re-generating the same rasterstack.
        message("stacked_aop_data already created for current tile.")
        
        # restore / read the rasterstack from file
        stacked_aop_data <- readRDS(file = stacked_aop_data_filename)
        
      } else { # this else statement runs until the end of this for-loop
        # if it doesn't exist, create the features from the aop data to file 
        message("creating stacked_aop_data for current tile...")
        
        # hyperspectral data ----------------------------------------
        # Build up the h5 rasterstack filename
        rasterstack_filename <- paste0(h5_dir, "rasterstack_",
                                       east_north_string, ".rds")
        
        print(paste("rasterstack filename: ", rasterstack_filename))
        
        # check to see if a .rds file already exists for the current tile
        # within the hyperspectral directory - just HS reflectance bands so far 
        
        if (file.exists(rasterstack_filename)){
          
          # if it exists, read that instead of re-generating the same rasterstack.
          message("reading hyperspectral rasterstack (already created for current tile)...")
          # restore / read the rasterstack from file
          s <- readRDS(file = rasterstack_filename)
          
        } else {
          # if it doesn't exist, generate the rasterstack. 
          message("creating hyperspectral rasterstack for current tile...")
          # create a georeferenced rasterstack using the current hyperspectral tile
          h5 <- grep(east_north_string, h5_list, value=TRUE)
          hs_stack <- stack_hyperspectral(h5)
          s <- mask_sg_hyperspectral(hs_stack, out_dir = dir_data_out, n_sg)
          
          # save the rasterstack to file 
          saveRDS(s, file = rasterstack_filename)
        }
        
        # It seems that in some cases tiles were not delivered in a 1000 m x 1000 m format
        # This causes errors when stacking. Therefore, if any of the tile coorinates 
        # does not match the filename coordinates, we will extend the tile with NA values
        tile_extent <- extent(as.numeric(veg_coordinates[i,"eastings"]), 
                              as.numeric(veg_coordinates[i,"eastings"])+1000, 
                              as.numeric(veg_coordinates[i,"northings"]),
                              as.numeric(veg_coordinates[i,"northings"])+1000)
        
        
        # read the corresponding remote sensing data layers for current tile
        # Lidar data ------------------------------------------------------------
        if (length(chm_list)==0) {
          print("NO lidar-derived Canopy Height Model layer to be added to the data cube...")
        } else {
          print("Reading lidar-derived Canopy Height Model layer for the data cube...")
          chm <- raster::raster(grep(east_north_string, chm_list, value=TRUE))
          # Make sure to read and save the entire content, not just links to the Temp folder
          chm <- readAll(chm)
          if (extent(chm) == tile_extent){chm <- chm } else {chm <- extend(chm,tile_extent)}
          # set the raster name for each layer to be simply the name of the data 
          names(chm) <- "chm"
          stacked_aop_data <- raster::addLayer(s, chm)
        }
        
        # Vegetation indices ------------------------------------------------------------
        if (length(veg_indices_list)==0) {
          print("NO vegetation indices to be added to the data cube...")
        } else {
          print("Reading vegetation indices for the training data cube...")
          # for the vegetation indices, go into the corresponding folder for current tile
          # and get a list of all the vegIndex geotiffs. then read all of those geotiffs 
          # into a single raster stack.
          veg_indices <- raster::stack(grep(east_north_string, 
                                            veg_indices_list, 
                                            value=TRUE))
          # Make sure to read and save the entire content, not just links to the Temp folder
          veg_indices <- readAll(veg_indices)
          if (extent(veg_indices) == tile_extent){veg_indices <- veg_indices } else {veg_indices <- extend(veg_indices,tile_extent)}
          # name each of the vegetation index layers based on the last piece of each 
          # respective filename, e.g. "NDVI" and 
          names(veg_indices) <- sapply(stringr::str_split(names(veg_indices),"_"),tail,1)
          stacked_aop_data <- raster::addLayer(stacked_aop_data, veg_indices)
          
        }
        
        # Other ------------------------------------------------------------

        # Create the pixel number grid as a layer to add to the data cube. 
        # this one keeps track of individual pixel ID's
        # to avoid duplicate spectra being extracted. Basically, assign an integer ID
        # to each pixel in the 1000x1000 raster. This raster needs to have the same 
        # dimensions, extent, crs as the other layers so they can be stacked together. 
        # create a vector of IDs from 1 to the number of pixels in one band (#rows x #cols)
        s <- stacked_aop_data
        pixelID <- 1:(nrow(s) * ncol(s))
        # add tile east, north coordinates - how to do this if raster values must be numeric? 
        #pixelID <- paste(pixelID, east_north_string, sep="_") 
        # reshape this 1D vector into a 2D matrix 
        dim(pixelID) <- c(nrow(s),ncol(s))
        # create a raster layer of pixel numbers 
        pixelNumbers <- raster::raster(pixelID, crs = crs(s))
        extent(pixelNumbers) <- extent(s)
        names(pixelNumbers) <- "pixelNumber"
        
        # Create similar layers to keep track of the tile where the pixel is located
        print("Creating layers for easting and northing per pixel...")
        eastingID <- rep(as.numeric(easting), times = (nrow(s) * ncol(s)))
        northingID <- rep(as.numeric(northing), times = (nrow(s) * ncol(s)))
        # reshape to be two-dimensional
        dim(eastingID) <- c(nrow(s),ncol(s))
        dim(northingID) <- c(nrow(s),ncol(s))
        # create rasters to contain the easting and northing values
        eastingIDs <- raster::raster(eastingID, crs = crs(s))
        northingIDs <- raster::raster(northingID, crs = crs(s))
        # assign extent and CRS to match the other layers in the stack
        extent(eastingIDs) <- extent(s)
        extent(northingIDs) <- extent(s)
        names(eastingIDs) <- "eastingIDs"
        names(northingIDs) <- "northingIDs"
        
        # now, all of the remote sensing data files have been read in for the current
        # tile. add each one to the hyperspectral data stack along with the 
        # layer to keep track pixel number within the tile. 
        stacked_aop_data <- raster::addLayer(stacked_aop_data, pixelNumbers, 
                                             eastingIDs, northingIDs)
        
        print("Stacked AOP data for current tile. ")
        
        # save the stacked AOP data to file for easy clipping later
        saveRDS(stacked_aop_data, file = stacked_aop_data_filename)
      }    }
    
  }
  
  end_time <- Sys.time()
  elapsed <- end_time - start_time
  print("Elapsed time: ")
  print(elapsed)
}

#######################################################
# Functions for stacking ancillary data that accompanies hyperspectral data collection

# Used in PART 2

stack_hyperspectral_ancillary <- function(h5){
  # list the contents of HDF5 file
  h5_struct <- rhdf5::h5ls(h5, all=T)
  
  ### Info on CRS
  # construct the string using "/Reflectance/Metadata/Coordinate_System",
  # without explicitly using a site code 
  crs_tag <- h5_struct$group[grepl("/Reflectance/Metadata/Coordinate_System", 
                                   h5_struct$group)][1] 
  
  # read coordinate reference system data
  crs_info <- rhdf5::h5read(h5, crs_tag)
  
  # convert "UTM" to lowercase "utm" for proper usage later
  crs_info$Proj4 <- CRS(chartr("UTM", "utm", crs_info$Proj4))
  
  ### Info on ancillary data
  # construct the string using "/Reflectance/Metadata/Ancillary_Imagery",
  # without explicitly using a site code 
  anc_tag <- h5_struct$group[grepl("/Reflectance/Metadata/Ancillary_Imagery", 
                                   h5_struct$group)][1] 
  
  # read the data
  anc_info <- rhdf5::h5read(h5, anc_tag)
  anc_info_weather <- anc_info[["Weather_Quality_Indicator"]]
  
  # get the dimensions of the ancillary weather data
  n_bands <- dim(anc_info_weather)[1]
  n_rows <- dim(anc_info_weather)[2]
  n_cols <- dim(anc_info_weather)[3]
  
  ### Info on sensor angles
  # construct the string using "/Reflectance/Metadata/",
  # without explicitly using a site code
  azimuth_tag <- paste0(h5_struct$group[grepl("/Reflectance/Metadata",
                                              h5_struct$group)][1],"/to-sensor_azimuth_angle")
  zenith_tag <- paste0(h5_struct$group[grepl("/Reflectance/Metadata",
                                             h5_struct$group)][1],"/to-sensor_zenith_angle")
  
  # read the data
  azimuth_info <- rhdf5::h5read(h5, azimuth_tag)
  zenith_info <- rhdf5::h5read(h5, zenith_tag)
  
  
  ### define spatial extent: extract resolution and origin coordinates
  map_info <- unlist(strsplit(crs_info$Map_Info, 
                              split = ", "))
  res_x <- as.numeric(map_info[6])
  res_y <- as.numeric(map_info[7])
  x_min <- as.numeric(map_info[4])
  y_max <- as.numeric(map_info[5])
  
  # calculate the maximum X and minimum Y values 
  x_max <- (x_min + (n_cols * res_x))
  y_min <- (y_max - (n_rows * res_y))
  tile_extent <- raster::extent(x_min, x_max, y_min, y_max)
  print("tile extent")
  print(tile_extent)
  
  ### Create rasterstack with all information
  # create georeferenced raster using band 1 
  r1 <- (anc_info_weather[1,,]) # convert first band to matrix
  # transpose the image pixels for proper orientation to match
  # the other layers. create a raster for this band and assign
  # the CRS.
  print("Transposing ancillary data for proper orientation")
  r1 <- raster::t(raster::raster(r1, crs = crs_info$Proj4))
  extent(r1) <- tile_extent
  
  # start the raster stack with first band 
  s <- raster::stack(r1)
  
  # loop through bands and create a rasterstack with 3 (n_bands) bands
  print("Stacking weather bands")
  for(b in 2:n_bands){
    # print(b)
    # create raster with current band
    r <- (anc_info_weather[b,,]) # convert to matrix
    r <- raster::t(raster::raster(r, crs = crs_info$Proj4))
    extent(r) <- tile_extent
    
    # add additional band to the stack with the addLayer function
    s <- raster::addLayer(s, r)
  }
  
  names(s) <- paste0("Weather_Quality_Indicator",seq(1:3))
  
  # add other ancillary information as additional layers
  anc_info_withoutWeather <- within(anc_info, rm("Weather_Quality_Indicator"))
  n_layers <- length(anc_info_withoutWeather) 
  print("Stacking other ancillary bands")
  for(l in 1:n_layers){
    # print(l)
    # create raster with current band
    r <- (anc_info_withoutWeather[[l]]) # convert to matrix
    r <- raster::t(raster::raster(r, crs = crs_info$Proj4))
    extent(r) <- tile_extent
    
    # add additional band to the stack with the addLayer function
    s <- raster::addLayer(s, r)
    names(s)[3+l] <- names(anc_info_withoutWeather[l])
  }
  
  # add other sensor angles information as additional layers
  azimuth_info_t <- raster::t(raster::raster(azimuth_info, crs = crs_info$Proj4))
  zenith_info_t <- raster::t(raster::raster(zenith_info, crs = crs_info$Proj4))
  extent(azimuth_info_t) <- tile_extent
  extent(zenith_info_t) <- tile_extent
  s <- raster::addLayer(s, azimuth_info_t)
  s <- raster::addLayer(s, zenith_info_t)
  names(s)[c(length(names(s))-1,length(names(s)))] <- c("to-sensor_azimuth_angle", "to-sensor_zenith_angle")
  
  return(s)
}

stack_neon_aop_ancillary <- function(veg_coordinates, dir_data_raw, dir_data_out){
  # Specify the paths for the data directory
  h5_dir <- file.path(dir_data_raw, "hyperspectral")
  
  # output directory for stacked AOP data
  stacked_aop_data_dir <- file.path(dir_data_out, "stacked_aop_data_ancillary")
  check_create_dir(stacked_aop_data_dir)
  
  # list the files in the data directory; filter results based on file type
  # hyperspectral data - list the .h5 files 
  h5_list <- list.files(path = h5_dir, full.names = TRUE) 
  
  # create data cubes with AOP-derived features 
  start_time <- Sys.time() # start the timer 
  
  # loop through the tiles; extract ancillary data for each tile 
  # for (h5 in h5_list) {
  for (i in 1:nrow(veg_coordinates)) {
    # each current hyperspectral tile must be read and stacked into a 
    # georeferenced rasterstack object (so it can be clipped with point / polygon
    # shapefiles). The process of creating a rasterstack takes a while for 
    # each tile, so after creating each rasterstack once, each object gets 
    # written to a file. 
    
    # Build up the rasterstack filename by combining the easting/northing coordinates,
    # and use this to find corresponding tiles of various remote sensing data
    easting <- veg_coordinates[i,"eastings"]
    northing <- veg_coordinates[i,"northings"]
    east_north_string <- paste(easting, northing, sep="_")
    # east_north_string <- paste(veg_coordinates[i,"eastings"], veg_coordinates[i,"northings"], sep="_")
    
    # generate a filename for the stacked AOP data
    stacked_aop_data_filename = file.path(stacked_aop_data_dir,
                                          paste0("stacked_aop_data_ancillary_",
                                                 east_north_string, ".rds"))
    
    # Before proceeding, first check if the polygon coordinates lie within the flight area
    h5 <- grep(east_north_string, h5_list, value=TRUE)
    if(length(h5) == 0){
      print(paste0("No (matching) tiles available for coordinates ", east_north_string))
    } else {
      # check if a .rds file already exists for the current feature data cube
      if (file.exists(stacked_aop_data_filename)){
        # if it exists, read that instead of re-generating the same rasterstack.
        message("stacked_aop_data already created for current tile.")
        
        # restore / read the rasterstack from file
        stacked_aop_data <- readRDS(file = stacked_aop_data_filename)
        
      } else { # this else statement runs until the end of this for-loop
        # if it doesn't exist, create the features from the aop data to file 
        message("creating stacked_aop_data for current tile...")
        
        # hyperspectral data ----------------------------------------
        # Build up the h5 rasterstack filename
        rasterstack_filename <- paste0(h5_dir, "rasterstack_ancillary_",
                                       east_north_string, ".rds")
        
        print(paste("rasterstack filename: ", rasterstack_filename))
        
        # check to see if a .rds file already exists for the current tile
        # within the hyperspectral directory - just HS reflectance bands so far 
        
        if (file.exists(rasterstack_filename)){
          
          # if it exists, read that instead of re-generating the same rasterstack.
          message("reading hyperspectral rasterstack (already created for current tile)...")
          # restore / read the rasterstack from file
          s <- readRDS(file = rasterstack_filename)
          
        } else {
          # if it doesn't exist, generate the rasterstack. 
          message("creating hyperspectral rasterstack for current tile...")
          # create a georeferenced rasterstack using the current hyperspectral tile
          h5 <- grep(east_north_string, h5_list, value=TRUE)
          s <- stack_hyperspectral_ancillary(h5)
          
          # save the rasterstack to file 
          saveRDS(s, file = rasterstack_filename)
        }
        
        # Create the pixel number grid as a layer to add to the data cube. 
        # this one keeps track of individual pixel ID's
        # to avoid duplicate spectra being extracted. Basically, assign an integer ID
        # to each pixel in the 1000x1000 raster. This raster needs to have the same 
        # dimensions, extent, crs as the other layers so they can be stacked together. 
        # create a vector of IDs from 1 to the number of pixels in one band (#rows x #cols)
        pixelID <- 1:(nrow(s) * ncol(s))
        # add tile east, north coordinates - how to do this if raster values must be numeric? 
        #pixelID <- paste(pixelID, east_north_string, sep="_") 
        # reshape this 1D vector into a 2D matrix 
        dim(pixelID) <- c(nrow(s),ncol(s))
        
        # create a raster layer of pixel numbers 
        pixelNumbers <- raster::raster(pixelID, crs = crs(s))
        extent(pixelNumbers) <- extent(s)
        names(pixelNumbers) <- "pixelNumber"
        
        # Create similar layers to keep track of the tile where the pixel is located
        print("Creating layers for easting and northing per pixel...")
        eastingID <- rep(as.numeric(easting), times = (nrow(s) * ncol(s)))
        northingID <- rep(as.numeric(northing), times = (nrow(s) * ncol(s)))
        # reshape to be two-dimensional
        dim(eastingID) <- c(nrow(s),ncol(s))
        dim(northingID) <- c(nrow(s),ncol(s))
        # create rasters to contain the easting and northing values
        eastingIDs <- raster::raster(eastingID, crs = crs(s))
        northingIDs <- raster::raster(northingID, crs = crs(s))
        # assign extent and CRS to match the other layers in the stack
        extent(eastingIDs) <- extent(s)
        extent(northingIDs) <- extent(s)
        names(eastingIDs) <- "eastingIDs"
        names(northingIDs) <- "northingIDs"
        
        stacked_aop_data <- raster::addLayer(s, pixelNumbers, 
                                             eastingIDs, northingIDs)
        
        print("Stacked AOP ancillary data for current tile. ")
        
        # save the stacked AOP data to file for easy clipping later
        saveRDS(stacked_aop_data, file = stacked_aop_data_filename)
        
      }
    }
  }   
  
  end_time <- Sys.time()
  elapsed <- end_time - start_time
  print("Elapsed time: ")
  print(elapsed)
}

#######################################################
# Function for stacking neon derivative products

# Used in PART 2

stack_neon_aop_derivatives <- function(veg_coordinates, dir_data_raw, dir_data_out){
  
  # Specify the paths for each data directory
  chm_dir <- file.path(dir_data_raw, "chm")
  slope_dir <- file.path(dir_data_raw, "slope")
  aspect_dir <-file.path(dir_data_raw, "aspect")
  dsm_dir <-file.path(dir_data_raw, "dsm")
  dtm_dir <-file.path(dir_data_raw, "dtm")
  # veg_indices_dir <- file.path(dir_data_raw, "veg_indices")
  veg_indices_dir <- file.path(move_dir, "veg_indices")
  lai_dir <-file.path(dir_data_raw, "lai")
  fpar_dir <-file.path(dir_data_raw, "fpar")
  
  # output directory for stacked AOP data
  stacked_aop_data_dir <- file.path(dir_data_out, "stacked_aop_data_derivatives")
  check_create_dir(stacked_aop_data_dir)
  
  # list the files in each data directory; filter results based on file type
  # hyperspectral data - list the .h5 files 
  chm_list <- list.files(path = chm_dir, full.names = TRUE)
  slope_list <- list.files(path = slope_dir, full.names = TRUE)
  aspect_list <- list.files(path = aspect_dir, full.names = TRUE)
  dsm_list <- list.files(path = dsm_dir, full.names = TRUE)
  dtm_list <- list.files(path = dtm_dir, full.names = TRUE)
  veg_indices_list <- list.files(path = veg_indices_dir, full.names = TRUE)
  lai_list <- list.files(path = lai_dir, full.names = TRUE)
  fpar_list <- list.files(path = fpar_dir, full.names = TRUE)
  # remove any entries from the list that are .zip files
  veg_indices_list<- veg_indices_list[grepl("*.tif$", veg_indices_list)]
  veg_indices_names <- c("ARVI","EVI","NDLI","NDNI","NDVI","PRI","SAVI")
  
  # create data cubes with AOP-derived features 
  start_time <- Sys.time() # start the timer 
  
  # loop through the tiles; build up a data cube for each tile 
  for (i in 1:nrow(veg_coordinates)) {
    
    # each current tile must be read and stacked into a 
    # georeferenced rasterstack object (so it can be clipped with point / polygon
    # shapefiles). The process of creating a rasterstack takes a while for 
    # each tile, so after creating each rasterstack once, each object gets 
    # written to a file. 
    
    # Build up the rasterstack filename by combining the easting/northing coordinates,
    # and use this to find corresponding tiles of various remote sensing data
    easting <- veg_coordinates[i,"eastings"]
    northing <- veg_coordinates[i,"northings"]
    east_north_string <- paste(easting, northing, sep="_")
    # east_north_string <- paste(veg_coordinates[i,"eastings"], veg_coordinates[i,"northings"], sep="_")
    
    # generate a filename for the stacked AOP data
    stacked_aop_data_filename = file.path(stacked_aop_data_dir,
                                          paste0("stacked_aop_data_derivatives_",
                                                 east_north_string, ".rds"))
    
    # Before proceeding, first check if the polygon coordinates lie within the flight area
    chm <- grep(east_north_string, chm_list, value=TRUE)
    veg_indices <- grep(east_north_string, veg_indices_list, value=TRUE)
    if(length(chm) == 0 | length(veg_indices) == 0){
      print(paste0("No tiles available for coordinates ", east_north_string))
    } else {
      
      # check if a .rds file already exists for the current feature data cube
      if (file.exists(stacked_aop_data_filename)){
        # if it exists, read that instead of re-generating the same rasterstack.
        message(paste0("stacked_aop_data already created for tile ",east_north_string))
        
        # restore / read the rasterstack from file
        stacked_aop_data <- readRDS(file = stacked_aop_data_filename)
        
      } else { # this else statement runs until the end of this for-loop
        # if it doesn't exist, create the features from the aop data to file 
        message(paste0("creating stacked_aop_data for tile ",east_north_string))
        
        # It seems that in some cases tiles were not delivered in a 1000 m x 1000 m format
        # This causes errors when stacking. Therefore, if any of the tile coorinates 
        # does not match the filename coordinates, we will extend the tile with NA values
        tile_extent <- extent(as.numeric(veg_coordinates[i,"eastings"]), 
                              as.numeric(veg_coordinates[i,"eastings"])+1000, 
                              as.numeric(veg_coordinates[i,"northings"]),
                              as.numeric(veg_coordinates[i,"northings"])+1000)
        
        # We build up the rasterstack starting from the chm raster, so this product needs to be present!
        
        # read the corresponding remote sensing data layers for current tile
        # Lidar data ------------------------------------------------------------
        print("Reading lidar-derived Canopy Height Model layer for the data cube...")
        chm <- raster::raster(grep(east_north_string, chm_list, value=TRUE))
        # Make sure to read and save the entire content, not just links to the Temp folder
        chm <- readAll(chm)
        if (extent(chm) == tile_extent){chm <- chm } else {chm <- extend(chm,tile_extent)}
        # set the raster name for each layer to be simply the name of the data 
        names(chm) <- "chm"
        stacked_aop_data <- chm
        
        if (length(slope_list)==0) {
          print("NO lidar-derived Aspect and Slope layers to be added to the data cube...")
        } else {
          print("Reading lidar-derived Aspect and Slope layers for the data cube...")
          east_north_string <- paste(veg_coordinates[i,"eastings"], veg_coordinates[i,"northings"], sep="_")
          slope <- raster::raster(grep(east_north_string, slope_list, value=TRUE))
          aspect <- raster::raster(grep(east_north_string, aspect_list, value=TRUE))
          # Make sure to read and save the entire content, not just links to the Temp folder
          slope <- readAll(slope)
          aspect <- readAll(aspect)
          if (extent(slope) == tile_extent){slope <- slope } else {
            slope <- extend(slope,tile_extent)
            slope <- crop(slope,tile_extent)}
          if (extent(aspect) == tile_extent){aspect <- aspect } else {
            aspect <- extend(aspect,tile_extent)
            aspect <- crop(aspect,tile_extent)
          }
          # set the raster name for each layer to be simply the name of the data 
          names(slope) <- "slope"
          names(aspect) <- "aspect"
          stacked_aop_data <- raster::addLayer(stacked_aop_data, slope, aspect)
        }
        
        if (length(dsm_list)==0) {
          print("NO lidar-derived DSM and DTM layers to be added to the data cube...")
        } else {
          print("Reading lidar-derived Aspect and Slope layers for the data cube...")
          dsm <- raster::raster(grep(east_north_string, dsm_list, value=TRUE))
          dtm <- raster::raster(grep(east_north_string, dtm_list, value=TRUE))
          # Make sure to read and save the entire content, not just links to the Temp folder
          dsm <- readAll(dsm)
          dtm <- readAll(dtm)
          if (extent(dsm) == tile_extent){dsm <- dsm } else {
            dsm <- extend(dsm,tile_extent)
            dsm <- crop(dsm,tile_extent)}
          if (extent(dtm) == tile_extent){dtm <- dtm } else {
            dtm <- extend(dtm,tile_extent)
            dtm <- crop(dtm,tile_extent)}
          # set the raster name for each layer to be simply the name of the data 
          names(dsm) <- "dsm"
          names(dtm) <- "dtm"
          stacked_aop_data <- raster::addLayer(stacked_aop_data, dsm, dtm)
        }
        
        
        # Vegetation indices ------------------------------------------------------------
        if (length(veg_indices_list)==0) {
          print("NO vegetation indices to be added to the data cube...")
        } else {
          print("Reading vegetation indices for the training data cube...")
          # for the vegetation indices, go into the corresponding folder for current tile
          # and get a list of all the vegIndex geotiffs. then read all of those geotiffs 
          # into a single raster stack.
          veg_indices <- raster::stack(grep(east_north_string, 
                                            veg_indices_list, 
                                            value=TRUE))
          # Make sure to read and save the entire content, not just links to the Temp folder
          veg_indices <- readAll(veg_indices)
          if (extent(veg_indices) == tile_extent){veg_indices <- veg_indices } else {veg_indices <- extend(veg_indices,tile_extent)}
          # name each of the vegetation index layers based on the last piece of each 
          # respective filename, e.g. "NDVI" and 
          names(veg_indices) <- sapply(stringr::str_split(names(veg_indices),"_"),tail,1)
          stacked_aop_data <- raster::addLayer(stacked_aop_data, veg_indices)
          
        }
        
        # LAI ------------------------------------------------------------
        if (length(lai_list)==0) {
          print("NO lidar-derived LAI layer to be added to the data cube...")
        } else {
          print("Reading LAI layer for the data cube...")
          lai <- raster::raster(grep(east_north_string, lai_list, value=TRUE))
          # Make sure to read and save the entire content, not just links to the Temp folder
          lai <- readAll(lai)
          if (extent(lai) == tile_extent){lai <- lai } else {lai <- extend(lai,tile_extent)}
          # set the raster name for each layer to be simply the name of the data 
          names(lai) <- "lai"
          stacked_aop_data <- raster::addLayer(stacked_aop_data, lai)
        }
        
        # fPAR ------------------------------------------------------------
        if (length(fpar_list)==0) {
          print("NO lidar-derived fPAR layer to be added to the data cube...")
        } else {
          print("Reading fPAR layer for the data cube...")
          fpar <- raster::raster(grep(east_north_string, fpar_list, value=TRUE))
          # Make sure to read and save the entire content, not just links to the Temp folder
          fpar <- readAll(fpar)
          if (extent(fpar) == tile_extent){fpar <- fpar } else {fpar <- extend(fpar,tile_extent)}
          # set the raster name for each layer to be simply the name of the data 
          names(fpar) <- "fpar"
          stacked_aop_data <- raster::addLayer(stacked_aop_data, fpar)
        }
        # Create the pixel number grid as a layer to add to the data cube. 
        # this one keeps track of individual pixel ID's
        # to avoid duplicate spectra being extracted. Basically, assign an integer ID
        # to each pixel in the 1000x1000 raster. This raster needs to have the same 
        # dimensions, extent, crs as the other layers so they can be stacked together. 
        # create a vector of IDs from 1 to the number of pixels in one band (#rows x #cols)
        s <- stacked_aop_data
        pixelID <- 1:(nrow(s) * ncol(s))
        # add tile east, north coordinates - how to do this if raster values must be numeric? 
        #pixelID <- paste(pixelID, east_north_string, sep="_") 
        # reshape this 1D vector into a 2D matrix 
        dim(pixelID) <- c(nrow(s),ncol(s))
        # create a raster layer of pixel numbers 
        pixelNumbers <- raster::raster(pixelID, crs = crs(s))
        extent(pixelNumbers) <- extent(s)
        names(pixelNumbers) <- "pixelNumber"
        
        # Create similar layers to keep track of the tile where the pixel is located
        print("Creating layers for easting and northing per pixel...")
        eastingID <- rep(as.numeric(easting), times = (nrow(s) * ncol(s)))
        northingID <- rep(as.numeric(northing), times = (nrow(s) * ncol(s)))
        # reshape to be two-dimensional
        dim(eastingID) <- c(nrow(s),ncol(s))
        dim(northingID) <- c(nrow(s),ncol(s))
        # create rasters to contain the easting and northing values
        eastingIDs <- raster::raster(eastingID, crs = crs(s))
        northingIDs <- raster::raster(northingID, crs = crs(s))
        # assign extent and CRS to match the other layers in the stack
        extent(eastingIDs) <- extent(s)
        extent(northingIDs) <- extent(s)
        names(eastingIDs) <- "eastingIDs"
        names(northingIDs) <- "northingIDs"
        
        
        
        # now, all of the remote sensing data files have been read in for the current
        # tile. add each one to the hyperspectral data stack along with the 
        # layer to keep track pixel number within the tile. 
        stacked_aop_data <- raster::addLayer(stacked_aop_data, pixelNumbers, 
                                             eastingIDs, northingIDs)
        
        # stacked_aop_data <- raster::addLayer(s, chm, slope, aspect, veg_indices, 
        #                                      rgb_features, pixelNumbers, 
        #                                      eastingIDs, northingIDs)
        print("Stacked AOP data for current tile. ")
        
        # save the stacked AOP data to file for easy clipping later
        saveRDS(stacked_aop_data, file = stacked_aop_data_filename)
      }    }
    
  }
  
  end_time <- Sys.time()
  elapsed <- end_time - start_time
  print("Elapsed time: ")
  print(elapsed)
}

#######################################################
# Functions to visualize AOP imagery for a specified polygon
# A polygon might cover more than 1 tile, and therefore a mosaicking process might be conducted

# Used in PART 3

# This script is an adapted version of the code written by Victoria Scholl:
# https://github.com/earthlab/neon-veg.git: 00-supporting-functions.R and 06-plot_aop_imagery.R

# Check this site for tips on improving plots:
# https://ourcodingclub.github.io/tutorials/spatial/

# Function to summarize list of coordinates to a list of the 1km x 1km tile coordinates covering them
list_tiles_covering_poly <- function(veg_df, out_dir){
  
  print("Generating list of tiles covering polygon(s)...")
  
  # get easting, northing coordinates
  e <- NA
  n <- NA
  for (i in 1:dim(veg_df)[1]){
    easting <- floor(as.numeric(veg_df$eastings[i])/1000)*1000
    northing <- floor(as.numeric(veg_df$northings[i])/1000)*1000
    e[i] <- easting
    n[i] <- northing
  }
  
  # find unique rows repesenting tiles
  easting_northings <- data.frame(as.character(e),as.character(n))
  colnames(easting_northings) <- c('e','n')
  tiles <- unique(easting_northings[,c('e','n')])
  
  # order by ascending tile coordinates 
  tiles <- tiles %>%
    arrange(e)
  
  tile_names <- paste(tiles$e, tiles$n, sep="_")
  
  # Optional: write list to a text file in the output directory provided
  # tiles_file <- file(file.path(out_dir,"list_tiles.txt"))
  # writeLines(tile_names, tiles_file)
  # close(tiles_file)
  
  # return(tiles)
  return(tile_names)
}

# Function to mosaic tiles from a list
mosaic_rasters <- function(raster_list){
  # x <- lapply(a, raster)
  names(raster_list)[1:2] <- c('x', 'y')
  raster_list$fun <- mean
  raster_list$na.rm <- TRUE
  mosaicked_rasters <- do.call(mosaic, raster_list) 
  return(mosaicked_rasters)
}

plot_aop_imagery <- function(filename, polyname, wavelengths, 
                             save_plots = FALSE, dir_data_out_maps, 
                             save_sign = FALSE, dir_data_out_sign, 
                             plot_rgb_hs = T, plot_chm = F, 
                             plot_slope = F, plot_aspect = F,plot_indices = T,
                             poly_sp_proj, add_poly=T){
  
  # read in the stacked AOP data for the specified polygon and plot the wanted feature maps
  
  ### Only 1 tile to be plotted
  if (length(filename) == 1){ 
    print("Reading single tile")
    stacked_aop_data <- readRDS(file = filename)
    
    # plot RGB composite using selected bands from the hyperspectral imagery.
    if (plot_rgb_hs == TRUE){
      # determine the index of each band closest to specific R,G,B wavelengths. 
      r_idx <- as.integer(which.min(abs(wavelengths - 669))) # Red 669 nm 
      g_idx <- as.integer(which.min(abs(wavelengths - 549))) # Green 549 nm 
      b_idx <- as.integer(which.min(abs(wavelengths - 474))) # Blue 474 nm 
      
      hs_rgb_stack <- raster::stack(stacked_aop_data@layers[[r_idx]],
                                    stacked_aop_data@layers[[g_idx]],
                                    stacked_aop_data@layers[[b_idx]])
      hs_rgb_brick <- raster::brick(hs_rgb_stack)
      raster::plotRGB(hs_rgb_brick,
                      r = 1, g = 2, b = 3,
                      stretch = "lin",
                      axes = TRUE,
                      main="Hyperspectrally derived \nRGB Composite using bands 669 nm, 549 nm and 474 nm",
                      #xlab="Easting (m)",
                      #ylab="Northing (m)", 
                      cex.main=1)
      if (add_poly == TRUE){
        plot(poly_sp_proj, add=T, border="red")
      }
      
      if(save_plots == TRUE){
        # write the RGB composite to an image file for visualizing in QGIS
        writeRaster(hs_rgb_brick,
                    file.path(dir_data_out_maps,paste0("hs_rg_", polyname, ".tif")), 
                    format="GTiff",
                    overwrite=TRUE)
      }
      
      if (save_sign == T){
        # generate 300 point samples from the polygons
        ptsamp <- spsample(poly_sp_proj, 300, type='regular')
        # extract values with points
        df <- raster::extract(stacked_aop_data,ptsamp,layer=1, nl=length(wavelengths))
        # To see some of the reflectance values
        # Create matrix with spectral statistics
        refl_stat <- rbind("max" = apply(df,2,max,na.rm=T), 
                           "mean" = apply(df,2,mean,na.rm=T),
                           "min" = apply(df,2,min,na.rm=T))
        colnames(refl_stat) <- wavelengths
        # Create a vector of color for the refl statistics for use in plotting
        library(colorspace)
        mycolor <- terrain_hcl(nrow(refl_stat))
        png(file= file.path(dir_data_out_sign,paste0("Spectral_signatures_", polyname, ".png")), 
            width=500, height=350)
        # First create an empty plot
        plot(0, ylim=c(0,0.8), xlim = c(400,2400), type='n', xlab="Wavelenghts (nm)", ylab = "Reflectance")
        # add the different statistics
        for (i in 1:nrow(refl_stat)){
          lines(x = wavelengths, y = refl_stat[i,], type = "l", lwd = 3, lty = 1, col = mycolor[i])
        }
        # Cover the atmospheric water absorption windows
        window1 <- c(1340, 1445)
        window2 <- c(1790, 1955)
        # We will also cover the 3 bands immediately at the edges of the water absoprtion windows (-3 and +3)
        window1_x1 <- wavelengths[as.integer(which.min(abs(wavelengths - window1[1])))-3]
        window1_x2 <- wavelengths[as.integer(which.min(abs(wavelengths - window1[2])))+3]
        window2_x1 <- wavelengths[as.integer(which.min(abs(wavelengths - window2[1])))-3]
        window2_x2 <- wavelengths[as.integer(which.min(abs(wavelengths - window2[2])))+3]
        # abline(v=c(window1_x1, window1_x2, window2_x1, window2_x2))
        polygon(x = c(window1_x1, window1_x1, window1_x2, window1_x2),                          
                y = c(0.01, 0.59, 0.59, 0.01),                             
                col = "white", border="white")                                    
        polygon(x = c(window2_x1, window2_x1, window2_x2, window2_x2),                          
                y = c(0.01, 0.59, 0.59, 0.01),                             
                col = "white", border="white")                                    
        title(main=paste0("Spectral Profile of ", sum(complete.cases(df)), " regularly sampled pixels"),
              font.main = 2)
        legend("topright", rownames(refl_stat),
               cex=0.8, col=mycolor, lty = 1, lwd =3, bty = "n")
        
        abline(v=c(669,549, 474), lty="dashed", col="gray")
        dev.off()
      }
    }
    
    # CHM
    if (plot_chm == TRUE){
      raster::plot(stacked_aop_data$chm,
                   col=grey(1:100/100),
                   axes = TRUE,
                   main="Lidar-Derived Canopy Height Model",
                   xlab="Easting (m)",
                   ylab="Northing (m)",
                   cex.main=1, # title text size 
                   legend.args=list(text='Height above ground [m]',side=2, font=2, 
                                    line=0.5, cex=0.8)) # legend on color bar
      if (add_poly == TRUE){
        plot(poly_sp_proj, add=T, border="red")
      }
    }
    
    # slope
    if (plot_slope == TRUE){
      raster::plot(stacked_aop_data$slope,
                   col=grey(1:100/100),
                   axes = TRUE,
                   main="Lidar-Derived Slope",
                   xlab="Easting (m)",
                   ylab="Northing (m)",
                   cex.main=1,
                   legend.args=list(text='Ratio of rise over run [degrees]',
                                    side=2, font=2, 
                                    line=0.5, cex=0.8)) # legend on color bar)
      if (add_poly == TRUE){
        plot(poly_sp_proj, add=T, border="red")
      }
    }
    
    # aspect 
    if (plot_aspect == TRUE){
      raster::plot(stacked_aop_data$aspect,
                   col=grey(1:100/100),
                   axes = TRUE,
                   main="Lidar-Derived Aspect",
                   xlab="Easting (m)",
                   ylab="Northing (m)",
                   cex.main=1,
                   legend.args=list(text='Direction of steepest slope \n [degrees from North]',
                                    side=2, font=2, 
                                    line=0.5, cex=0.8)) # legend on color bar)
      if (add_poly == TRUE){
        plot(poly_sp_proj, add=T, border="red")
      }
    } 
    
    # Vegetation indices
    if (plot_indices == TRUE){
      # ARVI vegetation index
      raster::plot(stacked_aop_data$NDVI,
                   col=grey(1:100/100),
                   axes = TRUE,
                   main="Hyperspectrally derived ARVI",
                   xlab="Easting (m)",
                   ylab="Northing (m)",
                   cex.main=1,
                   legend.args=list(text='ARVI [-1,1]',side=2, font=2, 
                                    line=0.5, cex=0.8)) # legend on color bar))
      if (add_poly == TRUE){
        plot(poly_sp_proj, add=T, border="red")
      }
      
      # EVI vegetation index
      raster::plot(stacked_aop_data$EVI,
                   col=grey(1:100/100),
                   axes = TRUE,
                   main="Hyperspectrally derived EVI",
                   xlab="Easting (m)",
                   ylab="Northing (m)",
                   cex.main=1,
                   legend.args=list(text='EVI [-1,1]',side=2, font=2, 
                                    line=0.5, cex=0.8)) # legend on color bar))
      if (add_poly == TRUE){
        plot(poly_sp_proj, add=T, border="red")
      }
      
      # NDLI vegetation index
      raster::plot(stacked_aop_data$NDLI,
                   col=grey(1:100/100),
                   axes = TRUE,
                   main="Hyperspectrally derived NDLI",
                   xlab="Easting (m)",
                   ylab="Northing (m)",
                   cex.main=1,
                   legend.args=list(text='NDLI [0,1]',side=2, font=2, 
                                    line=0.5, cex=0.8)) # legend on color bar))
      if (add_poly == TRUE){
        plot(poly_sp_proj, add=T, border="red")
      }
      
      # NDVI vegetation index
      raster::plot(stacked_aop_data$NDVI,
                   # col=grey(1:100/100),
                   col = rev(terrain.colors(10)),
                   axes = TRUE,
                   main="Hyperspectrally derived NDVI",
                   xlab="Easting (m)",
                   ylab="Northing (m)",
                   cex.main=1,
                   legend.args=list(text='NDVI [0,1]',side=2, font=2, 
                                    line=0.5, cex=0.8)) # legend on color bar))
      if (add_poly == TRUE){
        plot(poly_sp_proj, add=T, border="red")
      }
      
      # PRI vegetation index
      raster::plot(stacked_aop_data$PRI,
                   col=grey(1:100/100),
                   axes = TRUE,
                   main="Hyperspectrally derived PRI",
                   xlab="Easting (m)",
                   ylab="Northing (m)",
                   cex.main=1,
                   legend.args=list(text='PRI [-1,1]',side=2, font=2, 
                                    line=0.5, cex=0.8)) # legend on color bar))
      if (add_poly == TRUE){
        plot(poly_sp_proj, add=T, border="red")
      }
      
      # SAVI vegetation index
      raster::plot(stacked_aop_data$SAVI,
                   col=grey(1:100/100),
                   axes = TRUE,
                   main="Hyperspectrally derived SAVI",
                   xlab="Easting (m)",
                   ylab="Northing (m)",
                   cex.main=1,
                   legend.args=list(text='SAVI [-1,1]',side=2, font=2, 
                                    line=0.5, cex=0.8)) # legend on color bar))
      if (add_poly == TRUE){
        plot(poly_sp_proj, add=T, border="red")
      }
    }
    
    ### Multiple tiles need to be mosaicked  
  } else {
    aop_data_list <- vector(mode = "list", length = length(filename))
    names(aop_data_list) <- filename
    system.time(
      for (f in 1:length(filename)){
        print(paste0("Reading tile No. ", f, " of ", length(filename)))
        aop_data_list[[f]] <- readRDS(file = filename[f])
      }) 
    
    if (plot_rgb_hs == TRUE){
      # determine the index of each band closest to specific R,G,B wavelengths. 
      r_idx <- as.integer(which.min(abs(wavelengths - 669))) # Red 669 nm 
      g_idx <- as.integer(which.min(abs(wavelengths - 549))) # Green 549 nm 
      b_idx <- as.integer(which.min(abs(wavelengths - 474))) # Blue 474 nm 
      aop_data_list_rgb_hs <- lapply(aop_data_list, function(x) raster::subset(x, c(r_idx, g_idx, b_idx), drop=T))
      aop_data_mosaicked_rgb_hs <- mosaic_rasters(aop_data_list_rgb_hs)
      
      raster::plotRGB(aop_data_mosaicked_rgb_hs,
                      r = 1, g = 2, b = 3,
                      stretch = "lin",
                      axes = TRUE,
                      main="Hyperspectrally derived \nRGB Composite using bands 669 nm, 549 nm and 474 nm",
                      #xlab="Easting (m)",
                      #ylab="Northing (m)", 
                      cex.main=1)
      if (add_poly == TRUE){
        plot(poly_sp_proj, add=T, border="red")
      }
      
      if(save_plots == TRUE){
        # write the RGB composite to an image file for visualizing in QGIS
        writeRaster(aop_data_mosaicked_rgb_hs,
                    file.path(dir_data_out_maps,paste0("hs_rgb_", polyname, ".tif")), 
                    format="GTiff",
                    overwrite=TRUE)
      }
      
      if (save_sign == T){
        # generate 300 point samples from the polygons
        ptsamp <- spsample(poly_sp_proj, 300, type='regular')
        raster::plotRGB(aop_data_mosaicked_rgb_hs,
                        r = 1, g = 2, b = 3,
                        stretch = "lin",
                        axes = TRUE,
                        main="Hyperspectral-derived \nRGB Composite using bands 669nm, 549nm, 474nm",
                        #xlab="Easting (m)",
                        #ylab="Northing (m)", 
                        cex.main=1)
        plot(poly_sp_proj, add=T, border="red")
        plot(ptsamp, add=T)
        # extract values with points
        aop_sampled_list <- lapply(aop_data_list, function(x) raster::extract(x,ptsamp,layer=1, nl=length(wavelengths)))
        df <- do.call(rbind, aop_sampled_list)
        
        # To see some of the reflectance values
        # Create matrix with spectral statistics
        refl_stat <- rbind("max" = apply(df,2,max,na.rm=T), 
                           "mean" = apply(df,2,mean,na.rm=T),
                           "min" = apply(df,2,min,na.rm=T))
        colnames(refl_stat) <- wavelengths
        # Create a vector of color for the refl statistics for use in plotting
        library(colorspace)
        mycolor <- terrain_hcl(nrow(refl_stat))
        png(file= file.path(dir_data_out_sign,paste0("Spectral_signatures_", polyname, ".png")), 
            width=500, height=350)
        # First create an empty plot
        plot(0, ylim=c(0,0.8), xlim = c(400,2400), type='n', xlab="Wavelenghts (nm)", ylab = "Reflectance")
        # add the different statistics
        for (i in 1:nrow(refl_stat)){
          lines(x = wavelengths, y = refl_stat[i,], type = "l", lwd = 3, lty = 1, col = mycolor[i])
        }
        # Cover the atmospheric water absorption windows
        window1 <- c(1340, 1445)
        window2 <- c(1790, 1955)
        window1_x1 <- wavelengths[as.integer(which.min(abs(wavelengths - window1[1])))-3]
        window1_x2 <- wavelengths[as.integer(which.min(abs(wavelengths - window1[2])))+3]
        window2_x1 <- wavelengths[as.integer(which.min(abs(wavelengths - window2[1])))-3]
        window2_x2 <- wavelengths[as.integer(which.min(abs(wavelengths - window2[2])))+3]
        # abline(v=c(window1_x1, window1_x2, window2_x1, window2_x2))
        polygon(x = c(window1_x1, window1_x1, window1_x2, window1_x2),                          
                y = c(0.01, 0.59, 0.59, 0.01),                             
                col = "white", border="white")                                    
        polygon(x = c(window2_x1, window2_x1, window2_x2, window2_x2),                          
                y = c(0.01, 0.59, 0.59, 0.01),                             
                col = "white", border="white")                                    
        title(main=paste0("Spectral Profile of ", sum(complete.cases(df)), " regularly sampled pixels"),
              font.main = 2)
        legend("topright", rownames(refl_stat),
               cex=0.8, col=mycolor, lty = 1, lwd =3, bty = "n")
        
        abline(v=c(669,549, 474), lty="dashed", col="gray")
        dev.off()
      }
    }
    
    # CHM 
    if (plot_chm == TRUE){
      aop_data_list_chm <- lapply(aop_data_list, function(x) raster::subset(x, "chm", drop=T))
      aop_data_mosaicked_chm <- mosaic_rasters(aop_data_list_chm)
      raster::plot(aop_data_mosaicked_chm,
                   col=grey(1:100/100),
                   axes = TRUE,
                   main="Lidar-Derived Canopy Height Model",
                   xlab="Easting (m)",
                   ylab="Northing (m)",
                   cex.main=1, # title text size 
                   legend.args=list(text='Height above ground [m]',side=2, font=2, 
                                    line=0.5, cex=0.8)) # legend on color bar
      if (add_poly == TRUE){
        plot(poly_sp_proj, add=T, border="red")
      }
    }
    
    # slope
    if (plot_slope == TRUE){
      aop_data_list_slope <- lapply(aop_data_list, function(x) raster::subset(x, "slope", drop=T))
      aop_data_mosaicked_slope <- mosaic_rasters(aop_data_list_slope)
      raster::plot(aop_data_mosaicked_slope,
                   col=grey(1:100/100),
                   axes = TRUE,
                   main="Lidar-Derived Slope",
                   xlab="Easting (m)",
                   ylab="Northing (m)",
                   cex.main=1,
                   legend.args=list(text='Ratio of rise over run [degrees]',
                                    side=2, font=2, 
                                    line=0.5, cex=0.8)) # legend on color bar)
      if (add_poly == TRUE){
        plot(poly_sp_proj, add=T, border="red")
      }
    }
    
    # aspect 
    if (plot_aspect == TRUE){
      aop_data_list_aspect <- lapply(aop_data_list, function(x) raster::subset(x, "aspect", drop=T))
      aop_data_mosaicked_aspect <- mosaic_rasters(aop_data_list_aspect)
      raster::plot(aop_data_mosaicked_aspect,
                   col=grey(1:100/100),
                   axes = TRUE,
                   main="Lidar-Derived Aspect",
                   xlab="Easting (m)",
                   ylab="Northing (m)",
                   cex.main=1,
                   legend.args=list(text='Direction of steepest slope \n [degrees from North]',
                                    side=2, font=2, 
                                    line=0.5, cex=0.8)) # legend on color bar)
      if (add_poly == TRUE){
        plot(poly_sp_proj, add=T, border="red")
      }
    } 
    
    # Vegetation indices
    if (plot_indices == TRUE){
      
      # ARVI vegetation index
      aop_data_list_ARVI <- lapply(aop_data_list, function(x) raster::subset(x, "ARVI", drop=T))
      aop_data_mosaicked_ARVI <- mosaic_rasters(aop_data_list_ARVI)
      raster::plot(aop_data_mosaicked_ARVI,
                   # col=grey(1:100/100),
                   col = rev(terrain.colors(10)),
                   axes = TRUE,
                   main="Hyperspectrally derived ARVI",
                   xlab="Easting (m)",
                   ylab="Northing (m)",
                   cex.main=1,
                   legend.args=list(text='ARVI [-1,1]',side=2, font=2, 
                                    line=0.5, cex=0.8)) # legend on color bar))
      if (add_poly == TRUE){
        plot(poly_sp_proj, add=T, border="red")
      }
      
      # EVI vegetation index
      aop_data_list_EVI <- lapply(aop_data_list, function(x) raster::subset(x, "EVI", drop=T))
      aop_data_mosaicked_EVI <- mosaic_rasters(aop_data_list_EVI)
      raster::plot(aop_data_mosaicked_EVI,
                   # col=grey(1:100/100),
                   col = rev(terrain.colors(10)),
                   axes = TRUE,
                   main="Hyperspectrally derived EVI",
                   xlab="Easting (m)",
                   ylab="Northing (m)",
                   cex.main=1,
                   legend.args=list(text='EVI [-1,1]',side=2, font=2, 
                                    line=0.5, cex=0.8)) # legend on color bar))
      if (add_poly == TRUE){
        plot(poly_sp_proj, add=T, border="red")
      }
      
      # NDLI vegetation index
      aop_data_list_NDLI <- lapply(aop_data_list, function(x) raster::subset(x, "NDLI", drop=T))
      aop_data_mosaicked_NDLI <- mosaic_rasters(aop_data_list_NDLI)
      raster::plot(aop_data_mosaicked_NDLI,
                   # col=grey(1:100/100),
                   col = rev(terrain.colors(10)),
                   axes = TRUE,
                   main="Hyperspectrally derived NDLI",
                   xlab="Easting (m)",
                   ylab="Northing (m)",
                   cex.main=1,
                   legend.args=list(text='NDLI [0,1]',side=2, font=2, 
                                    line=0.5, cex=0.8)) # legend on color bar))
      if (add_poly == TRUE){
        plot(poly_sp_proj, add=T, border="red")
      }
      
      # NDVI vegetation index
      aop_data_list_NDVI <- lapply(aop_data_list, function(x) raster::subset(x, "NDVI", drop=T))
      aop_data_mosaicked_NDVI <- mosaic_rasters(aop_data_list_NDVI)
      raster::plot(aop_data_mosaicked_NDVI,
                   # col=grey(1:100/100),
                   col = rev(terrain.colors(10)),
                   axes = TRUE,
                   main="Hyperspectrally derived NDVI",
                   xlab="Easting (m)",
                   ylab="Northing (m)",
                   cex.main=1,
                   legend.args=list(text='NDVI [-1,1]',side=2, font=2, 
                                    line=0.5, cex=0.8)) # legend on color bar))
      if (add_poly == TRUE){
        plot(poly_sp_proj, add=T, border="red")
      }
      
      # PRI vegetation index
      aop_data_list_PRI <- lapply(aop_data_list, function(x) raster::subset(x, "PRI", drop=T))
      aop_data_mosaicked_PRI <- mosaic_rasters(aop_data_list_PRI)
      raster::plot(aop_data_mosaicked_PRI,
                   # col=grey(1:100/100),
                   col = rev(terrain.colors(10)),
                   axes = TRUE,
                   main="Hyperspectrally derived PRI",
                   xlab="Easting (m)",
                   ylab="Northing (m)",
                   cex.main=1,
                   legend.args=list(text='PRI [-1,1]',side=2, font=2, 
                                    line=0.5, cex=0.8)) # legend on color bar))
      if (add_poly == TRUE){
        plot(poly_sp_proj, add=T, border="red")
      }
      
      # SAVI vegetation index
      aop_data_list_SAVI <- lapply(aop_data_list, function(x) raster::subset(x, "SAVI", drop=T))
      aop_data_mosaicked_SAVI <- mosaic_rasters(aop_data_list_SAVI)
      raster::plot(aop_data_mosaicked_PRI,
                   # col=grey(1:100/100),
                   col = rev(terrain.colors(10)),
                   axes = TRUE,
                   main="Hyperspectrally derived SAVI",
                   xlab="Easting (m)",
                   ylab="Northing (m)",
                   cex.main=1,
                   legend.args=list(text='SAVI [-1,1]',side=2, font=2, 
                                    line=0.5, cex=0.8)) # legend on color bar))
      if (add_poly == TRUE){
        plot(poly_sp_proj, add=T, border="red")
      }
    }
  }
  
}

#######################################################
# Functions to preprocess and mosaic aop derivative tiles 

# Used in PART 5

## ---- Mosaic raster tiles to the extent of a polygon/plot/site --------------
mosaic_tiles <- function(aop_polymasked_list){
  if (length(aop_polymasked_list) == 1) {
    # access the first and only layer of the list
    aop_data_mosaicked <- aop_polymasked_list[[1]]
  } else {
    # Mosaic the layers in the list
    aop_data_mosaicked <- mosaic_rasters(aop_polymasked_list)
  }
  names(aop_data_mosaicked) <-names(aop_polymasked_list[[1]])
  
  # Make sure to read and save the entire content, not just links to the Temp folder
  # aop_data_mosaicked <- readAll(aop_data_mosaicked)
  
  return(aop_data_mosaicked)
}

## ---- Mask unwanted pixels in a tile (this can also be a mosaic)  --------------
mask_tile <- function(aop_data, ndvi_tresh){
  ## ---- Apply masks --------------
  ndvi <- subset(aop_data, "NDVI")
  # Show NDVI >= ndvi_tresh 
  # ndvi_mask <- ndvi >= ndvi_tresh
  ndvi_mask <- ndvi 
  ndvi_mask[ndvi_mask < ndvi_tresh] <- NA
  # plot(ndvi_mask)
  
  chm <- subset(aop_data, "chm")
  # Show chm < 2 (This is the treshold defined by NEON to identify trees)
  # chm_2 <- chm < 2
  chm_mask_2 <- chm
  chm_mask_2[chm >= 2] <- NA
  # plot(chm_2)
  
  combined_mask <- raster::mask(ndvi_mask, chm_mask_2, maskvalue = NA)
  
  aop_data_masked <- raster::mask(aop_data, combined_mask, maskvalue = NA)
  # updatevalue = NA, updateNA = T)
  
  freq_NA <- freq(aop_data_masked$chm)
  
  
  # Make sure to read and save the entire content, not just links to the Temp folder
  # aop_data_masked <- readAll(aop_data_masked)
  
  # The individual masks on the original data tiles, may be interesting for visualizing the workflow
  
  output <- list("aop_data_masked" = aop_data_masked,
                 "combined_ndvi_chm_mask" = combined_mask,
                 "freq_NA" = freq_NA)
  
  return(output)
}

## ---- Preprocess and mosaic all tiles of a site  --------------
derivatives_preprocessing <- function (aop_data_list, spdf, ndvi_tresh){
  
  ## Crop and mask raster to the extent of the feature
  print("Cropping tiles to the extent of the polygon")
  aop_cropped_list <- lapply(aop_data_list, function(x) raster::crop(x, extent(spdf)))
  print("Masking polygons from tiles")
  # system.time(aop_polymasked_list <- lapply(aop_cropped_list, function(x) raster::mask(x, rasterize(poly_sp_proj_ind, x))))
  system.time(aop_plotmasked_list <- lapply(aop_cropped_list, function(x) raster::mask(x, spdf)))
  
  ## mosaic tiles if a plot covers > 1 tile
  print("Mosaicking tiles")
  aop_data_mosaicked <- mosaic_tiles(aop_plotmasked_list)
  
  ## mask based on ndvi and chm treshold
  print("Masking mosaic based on NDVI and chm tresholds")
  aop_data_masked <- mask_tile(aop_data_mosaicked, ndvi_tresh = ndvi_tresh)
  
  output <- list("aop_data_mosaicked" = aop_data_mosaicked,
                 "aop_data_masked" = aop_data_masked$aop_data_masked,
                 "combined_ndvi_chm_mask"  = aop_data_masked$combined_ndvi_chm_mask,
                 "freq_NA" = aop_data_masked$freq_NA)
  return(output)
  
}

#######################################################
# Function to mosaic aop ancillary tiles

# Used in PART 5

ancillary_mosaicking <- function(aop_data_list, spdf){
  
  ## Crop and mask raster to the extent of the feature
  print("Cropping tiles to the extent of the polygon")
  aop_cropped_list <- lapply(aop_data_list, function(x) raster::crop(x, extent(spdf)))
  print("Masking polygons from tiles")
  system.time(aop_plotmasked_list <- lapply(aop_cropped_list, function(x) raster::mask(x, spdf)))
  
  ## mosaic tiles if a plot covers > 1 polygon
  print("Mosaicking tiles")
  aop_data_mosaicked <- mosaic_tiles(aop_plotmasked_list)
  
  return(aop_data_mosaicked)
  
}

#######################################################
# Function to spatially visualizate aop derivative image with plot locations

if (!require('RColorBrewer')) install.packages('RColorBrewer'); library('RColorBrewer')

plot_deriv_products_plot_loc <- function(deriv_prod, title, poly_sp_proj, veg_plots_coords, palette){
  raster::plot(deriv_prod,
               # col=grey(1:100/100),
               col=palette,
               axes = TRUE,
               main=paste0(title),
               xlab="Easting (m)",
               ylab="Northing (m)",
               cex.main=1, # title text size 
               legend.args=list(text=title,side=2, font=2, 
                                line=0.5, cex=0.8)) # legend on color bar
  plot(poly_sp_proj, add=T, border="black")
  plot(veg_plots_coords$veg_spdf_poly, add=T, border="blue")
  plot(veg_plots_coords$veg_spdf_point, add=T, col="blue", pch=18, cex=0.8)
} 

#######################################################
