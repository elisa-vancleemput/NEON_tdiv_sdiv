########################################################
###                                                  ###
###         Download NEON plot level data            ###
###                                                  ###
########################################################

# This script contains functions to download NEON plot level data on species composition and cover
# The script also provides some functions for performing basic visualizations for exploratory purposes

# by Elisa Van Cleemput, 2023
#######################################################
## ---- Download vegetation plots information --------------

library(dplyr, quietly=T)
if (!require('neonUtilities')) install.packages('neonUtilities'); library('neonUtilities')
if (!require('tidyverse')) install.packages('tidyverse'); library('tidyverse')
if (!require('devtools')) install.packages('devtools'); library('devtools')
# install_github('NEONScience/NEON-geolocation/geoNEON', dependencies=TRUE)
# library(geoNEON)
if (!require('geoNEON')) install_github('NEONScience/NEON-geolocation/geoNEON', dependencies=TRUE); library('geoNEON')
if (!require('rgeos')) install.packages('rgeos'); library('rgeos')

#######################################################
# Function to NEON species composition survey data

# Used in PART 4

Extract_plot_spcomp <- function(site_code, aop_data_year, dp_spcomp){
  ## ----munging-and-organizing-------------------------------
  # Code adapted from https://www.neonscience.org/resources/learning-hub/tutorials/aquatic-diversity-macroinvertebrates
  Sp_comp_site_year <- neonUtilities::loadByProduct(
    dpID = dp_spcomp,
    site = site_code,
    startdate = paste0(aop_data_year,"-01"),
    enddate = paste0(aop_data_year,"-12"),
    # token = NEON_TOKEN, #Uncomment to use your token
    check.size = F)
  Sp_comp_site_year %>% list2env(.GlobalEnv)
  
  # fieldData = div_1m2Data
  # fieldData = div_10m2Data100m2Data
  
  fieldData <- div_1m2Data %>%
    plyr::rbind.fill(div_10m2Data100m2Data)
  
  # Create unique IDs for every subplot x date observation
  fieldData$month <- substr(fieldData$endDate, 6, 7)
  fieldData$sampleID <- paste0(fieldData$plotID, "_",fieldData$subplotID, "_", fieldData$month)
  
  # extract date and location information into a separate table
  table_date_location <- fieldData %>%
    # keep only the columns listed below
    dplyr::select(endDate,month,
                  siteID, 
                  domainID,
                  plotID, subplotID,
                  nlcdClass,
                  namedLocation,
                  decimalLatitude, 
                  decimalLongitude, 
                  elevation) %>%
    # keep rows with unique combinations of values, 
    # i.e., no duplicate records
    distinct()
  
  
  # create a taxon table, which describes each taxonID that appears in the data set
  table_taxon <- fieldData %>%
    # keep only the columns listed below
    dplyr::select(divDataType, 
                  taxonID, scientificName,taxonRank,
                  family,nativeStatusCode, 
                  otherVariables,
                  identificationQualifier,
                  identificationReferences) %>%
    # remove rows with duplicate information
    distinct()
  
  # Create observation table.
  table_observation <- fieldData %>% 
    # dplyr::select a subset of columns from
    dplyr::select(uid,
                  domainID,
                  siteID,
                  plotID,subplotID,sampleID,
                  endDate,month,
                  percentCover,
                  taxonID,
                  family, 
                  scientificName,
                  taxonRank,
                  otherVariables)
  
  # extract sample metadata
  table_sample_info <- fieldData %>%
    dplyr::select(sampleID, 
                  siteID, 
                  domainID,
                  plotID, subplotID,
                  namedLocation, 
                  decimalLatitude, 
                  decimalLongitude, 
                  elevation,
                  endDate,month,
                  nlcdClass,
                  measuredBy, recordedBy) %>%
    distinct()
  # View(table_sample_info)
  
  # dplyr::select records with species cover info only (= subplot surveys) and remove duplicate records
  table_observation_species_subplots <- table_observation %>%
    dplyr::select(plotID, subplotID, sampleID, endDate, month, taxonID, percentCover) %>%
    distinct() %>%     # keep rows with unique combinations of values, i.e., no duplicate records
    dplyr::filter(!is.na(percentCover)) %>%  # dplyr::select records with species cover info only
    dplyr::filter(!is.na(taxonID)) # remove non-species records
  
  # dplyr::select records with species presence info (= subplot + baseplot surveys) and remove duplicate records
  table_observation_species_baseplots <- table_observation %>%
    dplyr::select(plotID, subplotID, sampleID, endDate, month, taxonID, percentCover) %>%
    distinct() %>%     # keep rows with unique combinations of values, i.e., no duplicate records
    dplyr::filter(!is.na(taxonID)) %>% # remove non-species records
    dplyr::select(-percentCover)
  
  # Create table to retrieve plot coordinates from
  table_observation_coord <- table_date_location %>%
    dplyr::select(!c(endDate, month)) %>% # some plots may be sampled > 1 time
    mutate(namedLocationSubplot = paste0(namedLocation, ".", subplotID)) %>%
    distinct()
  
  # location of base plot centroid points only: reference point position = 41 = middle of the square plot
  veg_loc_base <- geoNEON::getLocByName(table_observation_coord, "namedLocation", locOnly = T) 
  
  # locations of all subplot centroids
  # Calculate mapped UTM locations based on distance and azimuth
  table_observation_coord_sub <- table_observation_coord %>%
    dplyr::filter(subplotID %in% table_observation_species_subplots[,"subplotID"])
  veg_loc_sub <- geoNEON::def.calc.geo.os(table_observation_coord_sub, "div_1m2Data")
  # veg_loc_sub2 <- geoNEON::getLocTOS(table_observation_coord, "div_1m2Data")
  
  # in the case of veg_loc_sub, rename the adjEasting and adjNorthing columns 
  veg_loc_sub_utm <- veg_loc_sub %>% dplyr::rename(easting = adjEasting, northing = adjNorthing)
  
  # Remove all rows with NA coordinates (they had missing inputs for easting, 
  # northing, or UTM zone and were not converted)
  veg_loc_base_utm <- veg_loc_base[complete.cases(veg_loc_base[, c("northing","easting")]),]
  veg_loc_sub_utm <- veg_loc_sub_utm[complete.cases(veg_loc_sub_utm[, c("northing","easting")]),]
  
  veg_loc_base_utm <- veg_loc_base_utm %>%
    mutate(northing = as.numeric(northing)) %>%
    mutate(easting = as.numeric(easting))
  veg_loc_sub_utm <- veg_loc_sub_utm %>%
    mutate(northing = as.numeric(northing)) %>%
    mutate(easting = as.numeric(easting))
  
  # Get the Coordinate Reference System (CRS) for these coordinates 
  datum_base <- as.character(veg_loc_base_utm$geodeticDatum[1])
  zone_base <- gsub("[^0-9\\.]", "", veg_loc_base_utm$utmZoneNumber[1])
  coord_ref_base <- paste("+proj=utm +zone=",zone_base," +datum=",datum_base," +units=m"," +no_defs",sep="")
  
  
  output = list("table_date_location" = table_date_location,
                "table_taxon" = table_taxon,
                "table_observation" = table_observation,
                "table_sample_info" = table_sample_info,
                "table_observation_species_baseplots" = table_observation_species_baseplots,
                "table_observation_species_subplots" = table_observation_species_subplots,
                "veg_loc_base_utm" = veg_loc_base_utm,
                "veg_loc_sub_utm" = veg_loc_sub_utm,
                "coord_ref" = coord_ref_base)
  return(output)
}

#######################################################
# Functions to convert plot centroids to squares

# Used in PART 4

# Creates a SpatialPolygonsDataFrame of square plots for a list of northing-easting coordinates
create_plot_squares <- function(veg_utm, radius_m, coord_ref, ID){
  ## ---- Create square buffer around plot centroids  --------------
  # Code adapted from:
  # https://www.neonscience.org/resources/learning-hub/tutorials/field-data-polygons-centroids
  
  # set the radius for the plots
  radius <- radius_m # radius in meters
  # A radius of 10 corresponds with a plot of 20 m x 20 m
  
  centroids = veg_utm
  
  # define the plot edges based upon the plot radius. 
  yPlus <- centroids$northing+radius
  xPlus <- centroids$easting+radius
  yMinus <- centroids$northing-radius
  xMinus <- centroids$easting-radius
  
  # calculate polygon coordinates for each plot centroid. 
  square <- cbind(xMinus,yPlus,  # NW corner
                  xPlus, yPlus,  # NE corner
                  xPlus,yMinus,  # SE corner
                  xMinus,yMinus, # SW corner
                  xMinus,yPlus)  # NW corner again - close polygon
  
  # Extract the plot ID information
  ID <- centroids[[ID]]
  
  # create spatial polygons from coordinates
  veg_sp_poly <- SpatialPolygons(mapply(function(poly, id) {
    xy <- matrix(poly, ncol=2, byrow=TRUE)
    Polygons(list(Polygon(xy)), ID=id)
  }, 
  split(square, row(square)), ID),
  proj4string = CRS(coord_ref))
  # proj4string=CRS(as.character("+proj=utm +zone=11 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0")))
  
  # Create SpatialPolygonDataFrame -- this step is required to output multiple polygons.
  veg_spdf_poly <- SpatialPolygonsDataFrame(veg_sp_poly, data.frame(id=ID, row.names=ID))
  
  return(veg_spdf_poly)
}

# Returns spatial point (plot centroids) and polygon (plot squares) data frame of NEON base- or subplots
get_plot_coords <- function(veg_plots_spcomp_site_year, plot_type, radius_m, poly_sp_proj){
  
  if (plot_type == "base"){
    veg_utm <- veg_plots_spcomp_site_year$veg_loc_base_utm
    ID <- "plotID"
  } else if (plot_type == "sub") {
    veg_utm <- veg_plots_spcomp_site_year$veg_loc_sub_utm
    ID <- "namedLocationSubplot"
  } else {
    print("Plot type does not match with NEON plot types")
  }
  
  coord_ref = veg_plots_spcomp_site_year$coord_ref
  
  # Created squared buffer around centroids
  veg_spdf_poly <- create_plot_squares(veg_utm, radius_m, coord_ref, ID)
  
  # Convert the point data frame to an sp object 
  veg_spdf_point <- SpatialPointsDataFrame(coords = veg_utm[,c("easting", "northing")], 
                                           data = veg_utm,
                                           proj4string = CRS(coord_ref))
  # head(veg_utm)
  # head(veg_spdf_point)
  crs(veg_spdf_point)
  crs(poly_sp_proj)
  
  
  ## ---- Alligning crs of points, polygons and raster images --------------
  # Normally the plot coordinates are provided in the same projection as the aop data, 
  # but to be sure, we will check this and adapt the projection if this is not the case.
  if (compareCRS(veg_spdf_point, poly_sp_proj) == FALSE){
    veg_spdf_point <- spTransform(veg_spdf_point, crs(poly_sp_proj))
  }
  if (compareCRS(veg_spdf_poly, poly_sp_proj) == FALSE){
    veg_spdf_poly <- spTransform(veg_spdf_poly, crs(poly_sp_proj))
  }
  
  
  output <- list("veg_utm" = veg_utm,
                 "veg_spdf_poly" = veg_spdf_poly,
                 "veg_spdf_point" = veg_spdf_point)
  return(output)
}

# Returns spatial point (plot centroids) and polygon (plot squares) data frame of randomly sampled plots
get_plot_coords_fromSPDF <- function(SPDF, radius_m) {
  veg_spdf_point <- SPDF
  
  # Created squared buffer around centroids
  veg_utm <- as.data.frame(SPDF@coords)
  colnames(veg_utm) <- c("easting","northing")
  veg_utm$plotID <- SPDF[["plotID"]]
  coord_ref <- as.character(SPDF@proj4string)
  veg_spdf_poly <- create_plot_squares(veg_utm, radius_m, coord_ref, "plotID")
  
  ## ---- Alligning crs of points, polygons and raster images --------------
  # This is not necessary here, because we defined the point locations based on a raster image
  
  output <- list("veg_utm" = veg_utm,
                 "veg_spdf_poly" = veg_spdf_poly,
                 "veg_spdf_point" = veg_spdf_point)
  return(output)
}


#######################################################
# Functions to remove plots that are not located within field boundaries or that are not herbaceous/grassland

# Used in PART 4

# Removes plots that are not located in grassland and/or that do not lie within the boundaries of poly_sp_proj
remove_unwanted_neon_plots <- function(veg_plots_coords, poly_sp_proj, veg_plots_spcomp_site_year){
  
  veg_spdf_point <- veg_plots_coords$veg_spdf_point
  
  # Remove subplots that are not located in grassland
  veg_spdf_point <- veg_spdf_point[veg_spdf_point$nlcdClass == "grasslandHerbaceous",]
  
  # Remove plots that do not lie within the boundaries of poly_sp_proj
  intersects <- which(colSums(gIntersects(veg_spdf_point, poly_sp_proj, byid = TRUE)) > 0)
  # Warning saying that proj4 strings are different, but the crs are actually the same:
  if (compareCRS(veg_spdf_point, poly_sp_proj) == TRUE){
    print("CRS of plots and polygon is the same")
  }
  veg_spdf_point <- veg_spdf_point[intersects,]
  
  # Also subset the polgyon dataset and veg_utm data frame
  veg_spdf_poly <- veg_plots_coords$veg_spdf_poly
  veg_spdf_poly <- veg_spdf_poly[veg_spdf_poly$id %in% veg_spdf_point$id,]
  veg_utm <- veg_plots_coords$veg_utm
  veg_utm <- veg_utm[veg_utm$id %in% veg_spdf_point$id,]
  
  output <- list("veg_utm" = veg_utm,
                 "veg_spdf_poly" = veg_spdf_poly,
                 "veg_spdf_point" = veg_spdf_point)
  return(output)
}

#######################################################



