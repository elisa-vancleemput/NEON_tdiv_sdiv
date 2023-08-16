########################################################
###                                                  ###
###       Extract plot-level info from rasters       ###
###                                                  ###
########################################################

# This script overlaps plot boundaries with raster data obtained from 
#  - the Rangeland Analysis Platform
#  - NEON
# with the goal of obtaining plot-level values of NPP and abiotic conditions

# Using the GEE, we extracted raster products from the Rangeland Analysis Platform
# - cPFT cover %: annual forbs and grasses, perennial forbs and grasses, shrubs, trees, litter and bare ground 
# - cPFT NPP:  annual forbs and grasses, perennial forbs and grasses, shrubs, trees,  
# https://rangelands.app/products/  

# by Elisa Van Cleemput, 2023
#######################################################
# Functions to extract and visualize RAP data

# Function that extract values from RAP rasters, at locations in a spatial polygon dataframe or a spatial point dataframe
get_RAP_product_values_plot <- function(product, product_list, deriv_mosaic, veg_spdf, ID){
  
  cover_product_site <- product_list[[grep("cover", names(product_list))]]
  npp_product_site <- product_list[[grep("npp", names(product_list))]]

  # veg_spdf can be a spatial polygon dataframe or a spatial point dataframe
  
  ## ---- Alligning crs of spatial data frame and raster images --------------
  # For the sake of complementarity with previous analyses, 
  # we will apply the crs of veg_spdf to the new raster images, 
  # because this is the crs in which the NEON aop data were provided.
  cover_product_site_proj <- projectRaster(cover_product_site, crs = crs(veg_spdf))
  npp_product_site_proj <- projectRaster(npp_product_site, crs = crs(veg_spdf))

  ## ---- mask based on ndvi and chm treshold --------------
  # The mask might have a different extent than the RAP images, because
  # - the mask involves pre-processing (e.g., polygons removed)
  # - they have a different resolution (1 m vs. 30 m)
  # We will therefore resample the mask to the same resolution and extent as the RAP images.
  # We use nearest neighbor resampling
  mask_raster <- deriv_mosaic$combined_ndvi_chm_mask
  mask_raster_res_RAP <- raster::resample(mask_raster, cover_product_site_proj, method="ngb")
  cover_product_site_masked <- raster::mask(cover_product_site_proj, mask_raster_res_RAP, maskvalue = NA)
  npp_product_site_masked <- raster::mask(npp_product_site_proj, mask_raster_res_RAP, maskvalue = NA)

  ## ---- Extract information for each plot after masking --------------
  veg_plots_RAP1 <- bind_cols(as_tibble(raster::extract(cover_product_site_masked,veg_spdf,fun=median,na.rm = TRUE)),
                              as_tibble(raster::extract(npp_product_site_masked,veg_spdf,fun=median,na.rm = TRUE))) %>%
    dplyr::select(-QC) %>% # a column in npp_product_site_proj
    rename_all(~ paste0(product, "_", .x, "_masked")) %>%
    bind_cols(id = veg_spdf[[ID]])
  
  ## ---- Extract information for each plot without masking --------------
  veg_plots_RAP2 <- bind_cols(as_tibble(raster::extract(cover_product_site_proj,veg_spdf,fun=median,na.rm = TRUE)),
                              as_tibble(raster::extract(npp_product_site,veg_spdf,fun=median,na.rm = TRUE))) %>%
    dplyr::select(-QC) %>% # a column in npp_product_site_proj
    rename_all(~ paste(product, .x, sep="_")) %>%
    bind_cols(id = veg_spdf[[ID]])
  
  ## ---- Combine extracted RAP values --------------
  veg_plots_RAP <- veg_plots_RAP2 %>%
    full_join(veg_plots_RAP1, by=c("id"))
  
  return(veg_plots_RAP)
}

## ---- Visualize RAP products -------------
plot_RAP_products_plot_loc <- function(product_list, veg_plots_coords,
                              title, poly_sp_proj, palette){

  cover_product_site <- product_list[[grep("cover", names(product_list))]]
  npp_product_site <- product_list[[grep("npp", names(product_list))]]

  spdf <- veg_plots_coords

  ## ---- Alligning crs of polygons and raster images --------------
  # For the sake of complementarity with previous analyses, 
  # we will apply the crs of veg_spdf_poly to the new raster images, 
  # because this is the crs in which the NEON aop data were provided.
  cover_prod_site_proj <- projectRaster(cover_product_site, crs = crs(spdf))
  npp_prod_site_proj <- projectRaster(npp_product_site, crs = crs(spdf))

  ## ---- Visualize --------------
  raster::plot(cover_prod_site_proj$herbaceous,
               # col=grey(1:100/100),
               col=palette,
               axes = TRUE,
               main=paste0("Herbaceous cPFT cover ", title),
               xlab="Easting (m)",
               ylab="Northing (m)",
               cex.main=1, # title text size 
               legend.args=list(text=title,side=2, font=2, 
                                line=0.5, cex=0.8)) # legend on color bar
  plot(poly_sp_proj, add=T, border="black")
  if(class(spdf) == "SpatialPolygonsDataFrame"){
    plot(spdf, add=T, border ="blue")
  } else if (class(spdf) == "SpatialPointsDataFrame"){
    plot(spdf, add=T, col="blue", pch=18, cex=0.8)
    pointLabel(coordinates(spdf),labels=spdf$id, cex=0.5)
  }
  library(maptools)
  
  raster::plot(npp_prod_site_proj$herbaceousNPP,
               # col=grey(1:100/100),
               col=palette,
               axes = TRUE,
               main=paste0("Herbaceous NPP cover ", title),
               xlab="Easting (m)",
               ylab="Northing (m)",
               cex.main=1, # title text size 
               legend.args=list(text=title,side=2, font=2, 
                                line=0.5, cex=0.8)) # legend on color bar
  plot(poly_sp_proj, add=T, border="black")
  if(class(spdf) == "SpatialPolygonsDataFrame"){
    plot(spdf, add=T, border="blue")
  } else if (class(spdf) == "SpatialPointsDataFrame"){
    plot(spdf, add=T, col="blue", pch=18, cex=0.8)
  }
}

#######################################################
# Function to extract NEON data
get_NEON_product_values_plot <- function(deriv_mosaic, veg_spdf, ID){
  
  ## ---- Extract information for each plot after masking --------------
  # The mosaic_raster is already masked for ndvi and chm tresholds
  # and has the same coordinate system as veg_spdf
  veg_plots_NEON1 <- raster::extract(deriv_mosaic$aop_data_masked,veg_spdf,fun=median, na.rm = TRUE) %>%
    as.data.frame() %>%
    dplyr::select(dtm, slope, aspect, ARVI, EVI, NDVI, PRI, SAVI, lai) %>% 
    rename_all(~ paste0(.x, "_masked")) %>%
    bind_cols(id = veg_spdf[[ID]])
  
  ## ---- Extract information for each plot without masking --------------
  veg_plots_NEON2 <- raster::extract(deriv_mosaic$aop_data_mosaicked,veg_spdf,fun=median, na.rm = TRUE) %>%
    as.data.frame() %>%
    dplyr::select(dtm, slope, aspect, ARVI, EVI, NDVI, PRI, SAVI, lai) %>% 
    bind_cols(id = veg_spdf[[ID]])
  
  ## ---- Combine extracted values --------------
  veg_plots_NEON <- veg_plots_NEON2 %>%
    full_join(veg_plots_NEON1, by=c("id"))
  
  return(veg_plots_NEON)
}

#######################################################
