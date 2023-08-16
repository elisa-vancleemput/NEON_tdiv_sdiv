###########################################################
###                                                     ###
###              Calculate spectral diversity           ###
###                                                     ###
###########################################################

# This script contains functions to first pre-processes spectra, and then calculate spectral diversity

# by Elisa Van Cleemput, 2023
#######################################################
# Base functions for organizing data per plot and pre-processing its data

## ---- General functions to mosaic and mask tiles --------------
source("01_download_preprocess_aop_imagery.R") 
# functions: 
# - list_tiles_covering_poly to select tiles that cover a specific individual plot 
# - mosaic_tiles to mosaic raster tiles to the extent of a site/polygon/plot
# - mask_tile to mask unwanted pixels in a tile (this can also be a mosaic)

## ---- Clean spectral bands --------------
# As indicated in the main script, with this function we will remove the 3 bands immediately adjacent to the atmospheric water absorption windows
clean_spectra <- function(aop_data_tile, wavelengths) {
  # Select hyperspectral bands from the stacked image
  aop_data_hs <- raster::subset(aop_data_tile, 1:length(wavelengths))
  
  window1 <- c(1340, 1445)
  window2 <- c(1790, 1955)
  window1_left <- c(which.min(abs(wavelengths - window1[1])),which.min(abs(wavelengths - window1[1]))-1, which.min(abs(wavelengths - window1[1]))-2)
  window1_right <- c(which.min(abs(wavelengths - window1[2])),which.min(abs(wavelengths - window1[2]))+1, which.min(abs(wavelengths - window1[2]))+2)
  window2_left <- c(which.min(abs(wavelengths - window2[1])),which.min(abs(wavelengths - window2[1]))-1, which.min(abs(wavelengths - window2[1]))-2)
  window2_right <- c(which.min(abs(wavelengths - window2[2])),which.min(abs(wavelengths - window2[2]))+1, which.min(abs(wavelengths - window2[2]))+2)
  bad_bands <- c(window1_left, window1_right, window2_left, window2_right)
  aop_data_hs_clean <- dropLayer(aop_data_hs, bad_bands)
  
  # Make sure to read and save the entire content, not just links to the Temp folder
  # aop_data_hs_clean <- readAll(aop_data_hs_clean)
  
  # list of clean bands (n = 333)
  wvl <- cbind(as.data.frame(wavelengths), band = seq(1,length(wavelengths)))
  wvl_bands <- wvl %>%
    dplyr::filter(!band %in% bad_bands)
  wavelengths_clean <- wvl_bands$wavelengths
  
  output <- list("aop_data_hs_clean" = aop_data_hs_clean,
                 "wavelengths_clean"  = wavelengths_clean)
  return(output)
}

## ---- PCA --------------
# Dimension reduction and feature extraction: Principal Component Analysis
# Use the function written by Laliberté et al. 2020: pca
# https://github.com/elaliberte/specdiv/blob/master/bartlett/scripts/bartlett.R
library(devtools)
source_url("https://github.com/elaliberte/specdiv/blob/master/functions/specdiv_functions.R?raw=TRUE")

########################################################### 
# Wrapping function that applies the above base functions

spectral_preprocessing_per_plot <- function (aop_data_list, plot_ind, ndvi_tresh, scaling_nb, wavelengths){
  ### Crop and mask raster to the extent of the polygons
  print("Cropping tiles to the extent of the polygon")
  aop_cropped_list <- lapply(aop_data_list, function(x) raster::crop(x, extent(plot_ind)))
  print("Masking polygons from tiles")
  system.time(aop_plotmasked_list <- lapply(aop_cropped_list, function(x) raster::mask(x, plot_ind)))
  
  ### mosaic tiles if a plot covers > 1 tile
  # aop_data_mosaicked <- readAll(mosaic_tiles(aop_plotmasked_list))
  aop_data_mosaicked <- mosaic_tiles(aop_plotmasked_list)
  
  ## mask based on ndvi and chm treshold
  aop_data_masked <- mask_tile(aop_data_mosaicked, ndvi_tresh = ndvi_tresh)
  
  ## spectral cleaning
  cleaned_spectra <- clean_spectra(aop_data_mosaicked, wavelengths)
  cleaned_masked_spectra <- clean_spectra(aop_data_masked$aop_data_masked, wavelengths)
  
  # Run PCA (Warning: takes several minutes)
  # print("Performing PCA")
  system.time(aop_data_hs_pca <- pca(cleaned_spectra$aop_data_hs_clean, scaling = scaling_nb, p=1)) # include all PC's
 
  
  nb_not_na <- ncell(cleaned_masked_spectra$aop_data_hs_clean[[1]]) - freq(cleaned_masked_spectra$aop_data_hs_clean[[1]], value=NA)
  if (nb_not_na >= 3){
    system.time(aop_data_hs_masked_pca <- pca(cleaned_masked_spectra$aop_data_hs_clean, scaling = scaling_nb, p=1))
  } else {
    aop_data_hs_masked_pca <- NA
  }

  output <- list("aop_data_mosaicked" = aop_data_mosaicked,
                "cleaned_spectra" = cleaned_spectra,
                "cleaned_masked_spectra" = cleaned_masked_spectra,
                "aop_data_hs_pca" = aop_data_hs_pca, 
                "aop_data_hs_masked_pca" = aop_data_hs_masked_pca)
                 
  return(output)
  
}

########################################################### 
# Main function in which loop is initiated running the above functions on each plot

get_preprocessed_plots <- function(plot_set, ndvi_tresh, scaling_nb, site_code, aop_data_year, Dir_data_out, dir_data_out) {
  
  # Loop through the polygons dataset
  aop_plot_list <- list()
  
  for (i in 1:(dim(plot_set)[1])){
   
    print(paste0("Processing plot No. ", i, " of ", dim(plot_set)[1]))
    
    ## ---- Make a list of unique easting, northing coordinates (500 m steps) representing a specific plot --------------
    eastings <- c(seq(plot_set[i,]@bbox[1,"min"],plot_set[i,]@bbox[1,"max"],by=500),plot_set[i,]@bbox[1,"max"]) # by adding the max we ensure that the full polygon is covered for download
    northings <- c(seq(plot_set[i,]@bbox[2,"min"],plot_set[i,]@bbox[2,"max"],by=500),plot_set[i,]@bbox[2,"max"])
    poly_coordinates <- tidyr::expand_grid(eastings, northings)
    
    ## ---- Generate a list of the 1km x 1km tile coordinates covering this specific plot --------------
    tile_coordinates <- list_tiles_covering_poly(poly_coordinates)
    # tile_coordinates
    
    ## ---- Specify which tiles (.rds stacked AOP data file) to plot --------------
    # If the tiles were already downloaded and stored previously, access them via the Data folder,
    # if they were downloaded during the current R sesseion, access them via the data folder
    # Attention: this if-loop only makes sense when entire sites were downloaded and stacked at once (see above)
    if (length(list.files(file.path(Dir_data_out,paste(site_code,aop_data_year,sep="_"),"stacked_aop_data"),all.files = TRUE)) != 0){
      Dir_data_out_site_year <- file.path(Dir_data_out,paste(site_code,aop_data_year,sep="_"))
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
    
    ## ---- Load the raster tiles covering the polygon(s) --------------
    filename = stacked_aop_data_filename
    aop_data_list <- vector(mode = "list", length = length(filename))
    names(aop_data_list) <- filename
    system.time(
      for (f in 1:length(filename)){
        # print(paste0("Reading tile No. ", f, " of ", length(filename)))
        aop_data_list[[f]] <- readRDS(file = filename[f])
      })
    
    ## ---- Apply spectral preprocessing --------------
    aop_plot_list[[i]] <- spectral_preprocessing_per_plot(aop_data_list, plot_set[i,], ndvi_tresh, scaling_nb, wavelengths)
  }
  
  names(aop_plot_list) <- plot_set$id

  
  return(aop_plot_list)
}

# Parallel processing version
get_preprocessed_plots_parallel <- function(plot_set, ndvi_tresh, scaling_nb, site_code, aop_data_year, Dir_data_out, dir_data_out, cores_nb) {
  
  # Loop through the polygons dataset, using parallel computing
  cores <- cores_nb
  cl <- makeCluster(cores)
  registerDoSNOW(cl)
  iterations <- dim(plot_set)[1]
  pb <- txtProgressBar(max = iterations, style = 3)
  progress <- function(n) setTxtProgressBar(pb, n)
  opts <- list(progress = progress)
  
  system.time(aop_plot_list <- foreach(i = 1:(dim(plot_set)[1]),
                                       .final = function(x) setNames(x, plot_set$id),
                                       .inorder = TRUE,
                                       .packages=c("raster", "tidyr", "magrittr", "dplyr", "tidyverse", "sp"),
                                       .export=c("list_tiles_covering_poly",
                                                 "mosaic_tiles", "mosaic_rasters",
                                                 "mask_tile", "clean_spectra", "pca",
                                                 "spectral_preprocessing_per_plot_unmasked"),
                                       .options.snow = opts,
                                       .verbose = TRUE) %dopar% {
                                         
                                         ## ---- Make a list of unique easting, northing coordinates (500 m steps) representing a specific plot --------------
                                         eastings <- c(seq(plot_set[i,]@bbox[1,"min"],plot_set[i,]@bbox[1,"max"],by=500),plot_set[i,]@bbox[1,"max"]) # by adding the max we ensure that the full polygon is covered for download
                                         northings <- c(seq(plot_set[i,]@bbox[2,"min"],plot_set[i,]@bbox[2,"max"],by=500),plot_set[i,]@bbox[2,"max"])
                                         # poly_coordinates <- rbind(poly_coordinates,tidyr::expand_grid(eastings, northings))
                                         poly_coordinates <- tidyr::expand_grid(eastings, northings)
                                         
                                         ## ----  generate a list of the 1km x 1km tile coordinates covering this specific plot --------------
                                         tile_coordinates <- list_tiles_covering_poly(poly_coordinates)
                                         
                                         ## ---- Specify which tiles (.rds stacked AOP data file) to plot --------------
                                         # If the tiles were already downloaded and stored previously, access them via the Data folder,
                                         # if they were downloaded during the current R sesseion, access them via the data folder
                                         # Attention: this if-loop only makes sense when entire sites were downloaded and stacked at once (see above)
                                         if (length(list.files(file.path(Dir_data_out,paste(site_code,aop_data_year,sep="_"),"stacked_aop_data"),all.files = TRUE)) != 0){
                                           Dir_data_out_site_year <- file.path(Dir_data_out,paste(site_code,aop_data_year,sep="_"))
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
                                         
                                         ## ---- Load the raster tiles covering the polygon(s) --------------
                                         filename = stacked_aop_data_filename
                                         aop_data_list <- vector(mode = "list", length = length(filename))
                                         names(aop_data_list) <- filename
                                         system.time(
                                           for (f in 1:length(filename)){
                                             print(paste0("Reading tile No. ", f, " of ", length(filename)))
                                             aop_data_list[[f]] <- readRDS(file = filename[f])
                                           })
                                         
                                         ## ---- Apply spectral preprocessing --------------
                                         spectral_preprocessing_per_plot(aop_data_list, plot_set[i,], ndvi_tresh, scaling_nb, wavelengths)
                                         
                                       })
  close(pb)
  stopCluster(cl)
  
  return(aop_plot_list)
}

########################################################### 
# Calculate spectral diversity

# Coefficient of variation
get_cv_plot <- function(aop_poly_list){
  cv_output <- data.frame(matrix(NA, nrow=length(aop_poly_list), ncol=4))
  rownames(cv_output) <- names(aop_poly_list)
  colnames(cv_output) <- c("cv", "cv_npixels", "cv_masked", "cv_masked_npixels")
  
  for(i in 1:length(aop_poly_list)){
    cube_hs_brick <- aop_poly_list[[i]]$cleaned_spectra$aop_data_hs_clean
    cube_hs_ma <- as.matrix(as.data.frame(cube_hs_brick))
    mean.spec <- t(as.matrix(apply(cube_hs_ma,2,mean)))
    sd.spec <- apply(cube_hs_ma,2,sd)
    cv_output[i,"cv"] <- mean(sd.spec/mean.spec)
    cv_output[i,"cv_npixels"] <- nrow(cube_hs_ma)
    
    cube_masked_hs_brick <- aop_poly_list[[i]]$cleaned_masked_spectra$aop_data_hs_clean
    cube_masked_hs_ma <- as.matrix(as.data.frame(cube_masked_hs_brick))
    cube_masked_hs_ma <- cube_masked_hs_ma[complete.cases(cube_masked_hs_ma),]
    if(is.null(dim(cube_masked_hs_ma))) { # If no more than one pixel left after masking
      cv_output[i,"cv_masked"] <- NA
      cv_output[i,"cv_masked_npixels"] <- NA
    } else {
      mean.spec_masked <- t(as.matrix(apply(cube_masked_hs_ma,2,mean)))
      sd.spec_masked <- apply(cube_masked_hs_ma,2,sd)
      cv_output[i,"cv_masked"] <- mean(sd.spec_masked/mean.spec_masked)
      cv_output[i,"cv_masked_npixels"] <- nrow(cube_masked_hs_ma)
    }

  }
  
  cv_output <- tibble::rownames_to_column(cv_output, "id")
  
  return(cv_output)
}

# Convex hull volume
if (!require('geometry')) install.packages('geometry'); library('geometry')
get_chv_plot <- function(aop_poly_list) {
  chv_output <- data.frame(matrix(NA, nrow=length(aop_poly_list), ncol=6))
  rownames(chv_output) <- names(aop_poly_list)
  colnames(chv_output) <- c("chv", "chv_npixels", "nPC", "chv_masked", "chv_masked_npixels", "nPC_masked")
  
  for (i in 1:length(aop_poly_list)){
    # 1) CHV on original image
    cube_pc_brick <- aop_poly_list[[i]]$aop_data_hs_pca$cube_pc 
    if (nlayers(cube_pc_brick) > 3){
      cube_pc_brick <- subset(cube_pc_brick, 1:3)
    } else {
      cube_pc_brick <- cube_pc_brick
    }
    cube_pc_ma <- as.matrix(as.data.frame(cube_pc_brick))
    chv_output[i,"chv"] <- convhulln(cube_pc_ma,"FA")$vol
    chv_output[i,"chv_npixels"] <- nrow(cube_pc_ma)
    chv_output[i,"nPC"] <- 3
    # chv_output[i,"nPC"] <- nlayers(aop_poly_list[[i]]$aop_data_hs_pca$cube_pc)
    
    # 2) CHV on masked image
    if (is.na(cube_masked_pc_brick <- aop_poly_list[[i]]$aop_data_hs_masked_pca[1])){
      # not enough pixels left after masking to perform PCA
      chv_output[i,"chv_masked"] <- NA
      chv_output[i,"chv_masked_npixels"] <- NA 
      chv_output[i,"nPC_masked"] <- NA
    } else {
      cube_masked_pc_brick <- aop_poly_list[[i]]$aop_data_hs_masked_pca$cube_pc 
      # if (nlayers(cube_masked_pc_brick) > 3){
        cube_masked_pc_brick <- subset(cube_masked_pc_brick, 1:3) # we will use 3 PCs
      # } else {
      #   cube_masked_pc_brick <- cube_masked_pc_brick
      # }
      cube_masked_pc_ma <- as.matrix(as.data.frame(cube_masked_pc_brick))
      cube_masked_pc_ma <- cube_masked_pc_ma[complete.cases(cube_masked_pc_ma),]
      chv_output[i,"chv_masked"] <- convhulln(cube_masked_pc_ma,"FA")$vol
      chv_output[i,"chv_masked_npixels"] <- nrow(cube_masked_pc_ma)
      chv_output[i,"nPC_masked"] <- 3
      # chv_output[i,"nPC_masked"] <- nlayers(aop_poly_list[[i]]$aop_data_hs_masked_pca$cube_pc)
    }
    
  }  
  
  chv_output <- tibble::rownames_to_column(chv_output, "id")
  
  return(chv_output)
}

# Mean spectral angle measure
# Function "spectral_angle" copied from Henry Frye:
# https://github.com/henryf6/GCFRSpectralSurrogacy/blob/main/spectral_div_funcs.R
#based on Kruse et al., 1993
rad2deg <- function(rad) {(rad * 180) / (pi)}
spectral_angle <- function(spec, ref){
  rad2deg(acos( (sum(spec*ref)/
                   (sqrt(sum(spec*spec))*sqrt(sum(ref*ref))))
                
  ))
}
get_sam_plot <- function(aop_poly_list){
  sam_output <- data.frame(matrix(NA, nrow=length(aop_poly_list), ncol=4))
  rownames(sam_output) <- names(aop_poly_list)
  colnames(sam_output) <- c("sam", "sam_npixels", "sam_masked", "sam_masked_npixels")
  
  for(i in 1:length(aop_poly_list)){
    cube_hs_brick <- aop_poly_list[[i]]$cleaned_spectra$aop_data_hs_clean
    cube_hs_ma <- as.matrix(as.data.frame(cube_hs_brick))
    mean.spec <- t(as.matrix(apply(cube_hs_ma,2,mean)))
    sam_output[i,"sam"] <- mean(apply(cube_hs_ma, 1, spectral_angle, ref = mean.spec)) 
    sam_output[i,"sam_npixels"] <- nrow(cube_hs_ma)
    
    cube_masked_hs_brick <- aop_poly_list[[i]]$cleaned_masked_spectra$aop_data_hs_clean
    cube_masked_hs_ma <- as.matrix(as.data.frame(cube_masked_hs_brick))
    cube_masked_hs_ma <- cube_masked_hs_ma[complete.cases(cube_masked_hs_ma),]
    if(is.null(dim(cube_masked_hs_ma))) { # If no more than one pixel left after masking 
      sam_output[i,"sam_masked"] <- NA
      sam_output[i,"sam_masked_npixels"] <- NA
      
    } else {
      mean.spec_masked <- t(as.matrix(apply(cube_masked_hs_ma,2,mean)))
      sam_output[i,"sam_masked"] <- mean(apply(cube_masked_hs_ma, 1, spectral_angle, ref = mean.spec_masked)) 
      sam_output[i,"sam_masked_npixels"] <- nrow(cube_masked_hs_ma)
    }
  }
  
  sam_output <- tibble::rownames_to_column(sam_output, "id")
  
  return(sam_output)
}

########################################################### 
