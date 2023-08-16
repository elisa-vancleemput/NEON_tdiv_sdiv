##############################################################################
###                                                                        ###
###             Main script for base plot-level statistical analysis       ###
###                                                                        ###
##############################################################################
# Code accompanying
# Van Cleemput, E., Adler, P., Suding, K. in Global Ecology and Geography
# Making remote sense of biodiversity: What grassland characteristics make spectral diversity a good proxy for taxonomic diversity?

# Prior to this script, 00_Main_script_downloads should be run to download and prepare the data
# Script supporting this code
# - 01_Download_preprocess_aop_imagery.R: only for the functions "list_tiles_covering_poly_tibble" and "list_tiles_covering_poly"
# - 05_Analyses_functions.R

# This script contains 10 parts, parts 2-10 can be run independently after running general settings and part 1

### GENERAL LINES
# General settings and Data properties specifications
# PART 1: Read and prepare data
# PART 2: Explore data

### ANALYSES
# PART 3: Within and across sites - linear models: Taxonomic vs. spectral diversity
#  --> Output in Tables 1, S2 and S4, Fig. 4
# PART 4: Across sites - linear model: Taxonomic vs. spectral diversity with site interaction
#  --> Output in Tables 1 and S3
# PART 5: Calculate species turnover
#  --> Creation of Fig. S4
# PART 6: Evaluate relationship between site-level characteristics and strength of tdiv-sdiv relationship
#  --> Output in Tables 2, S5, S6 and S7, Fig. 5
# PART 7: Paired correlation plots and boxplots
#  --> Creation of Figs. 3 and S2

### VISUALIZATIONS FOR ILLUSTRATIVE PURPOSES ONLY
## PART 8: Visualize RGB of baseplots
#  --> Creation of Fig. S3
## PART 9: Visualize 1st spectral PC and hyperspectral signatures of baseplots
#  --> Creation of Fig. S3
## PART 10: Create illustrative hyperspectral signatures  
#  --> Creation of Fig. 1

# by Elisa Van Cleemput, 2023
#################################################################################
## General settings -------
#################################################################################
rm(list=ls())

# load general packages
library(magrittr)
library(dplyr)
library(raster)
library(plyr)
library(dplyr)
library(ggplot2)
if (!require('reshape2')) install.packages('reshape2'); library('reshape2')
if (!require('ggsci')) install.packages('ggsci'); library('ggsci')
# https://cran.r-project.org/web/packages/ggsci/vignettes/ggsci.html

# Set global option to NOT convert all character variables to factors
options(stringsAsFactors=F)

# ---- General function to create folders ---- 
check_create_dir <- function(new_dir){
  # check if directory exists. If it doesn't, create it. 
  if (!dir.exists(new_dir)){
    dir.create(new_dir)
  }
}

# ---- Specify input data folders and create output data folders ---- 
# Specify main directories where scripts and data are stored
main_dir_local <- # FILL IN
main_dir_external <- # FILL IN

# 1) Create folder to store raw data
dir_data_raw <- file.path(main_dir_external,"Data","data_raw")
dir_data_out_temp <- file.path(dir_data_raw,"Temp_AOPtiles")
check_create_dir(dir_data_out_temp)

# 2) Create folders to store processed data
main_dir_output <- main_dir_external
dir_data_out <- file.path(main_dir_output, "Data", "data_output")
check_create_dir(dir_data_out)

# Create more specific output folders
dir_data_out_results <- file.path(dir_data_out, "Across_sites_results")
check_create_dir(dir_data_out_results)

dir_results_explorative <- file.path(dir_data_out_results,"explorative")
check_create_dir(dir_results_explorative)

dir_results_tdiv_vs_sdiv_general <- file.path(dir_data_out_results,"tdiv_vs_sdiv_general")
check_create_dir(dir_results_tdiv_vs_sdiv_general)

dir_results_tdiv_vs_sdiv_general_Sitescombined <- file.path(dir_data_out_results,"tdiv_vs_sdiv_general_Sitescombined")
check_create_dir(dir_results_tdiv_vs_sdiv_general_Sitescombined)

dir_results_patchiness <- file.path(dir_data_out_results,"patchiness")
check_create_dir(dir_results_patchiness)

dir_results_expl_general <- file.path(dir_data_out_results,"tdiv_vs_sdiv_general_expl")
check_create_dir(dir_results_expl_general)

# 3) Folder where R scripts are stored
dir_scripts <- file.path(main_dir_local, "NEON_tdiv_sdiv")
setwd(dir_scripts)

# Source script with functions that will be used in this script
source("05_Analyses_functions.R")

#################################################################################
## Specify data properties -------
#################################################################################
# 1) NDVI treshold used for masking pixels during pre-processing of raw data
ndvi_tresh <- 0.2         

# 2) Year and days of flights
site_codes <- c("ABBY", "CLBJ", "CPER", "KONZ", "NIWO", "NOGP", "OAES", "SJER", "WOOD")
aop_data_year <- "2017"     # Four-digit year in character string "YYYY" for AOP imagery

# First and last flight date at every site
flight_time_ABBY1 <- as.POSIXct("2017-06-07", tz="GMT")
flight_time_ABBY2 <- as.POSIXct("2017-06-22", tz="GMT")
flight_time_CLBJ1 <- as.POSIXct("2017-05-06", tz="GMT")
flight_time_CLBJ2 <- as.POSIXct("2017-05-25", tz="GMT")
flight_time_CPER1 <- as.POSIXct("2017-05-24", tz="GMT")
flight_time_CPER2 <- as.POSIXct("2017-05-24", tz="GMT")
flight_time_KONZ1 <- as.POSIXct("2017-06-05", tz="GMT") 
flight_time_KONZ2 <- as.POSIXct("2017-06-09", tz="GMT") 
flight_time_NIWO1 <- as.POSIXct("2017-09-04", tz="GMT")
flight_time_NIWO2 <- as.POSIXct("2017-09-21", tz="GMT")
flight_time_NOGP1 <- as.POSIXct("2017-06-19", tz="GMT")
flight_time_NOGP2 <- as.POSIXct("2017-06-21", tz="GMT")
flight_time_OAES1 <- as.POSIXct("2017-05-01", tz="GMT")
flight_time_OAES2 <- as.POSIXct("2017-05-05", tz="GMT")
flight_time_SJER1 <- as.POSIXct("2017-03-28", tz="GMT")
flight_time_SJER2 <- as.POSIXct("2017-03-28", tz="GMT")
flight_time_WOOD1 <- as.POSIXct("2017-06-19", tz="GMT")
flight_time_WOOD2 <- as.POSIXct("2017-06-26", tz="GMT")
flight_times1 <- c(flight_time_ABBY1, flight_time_CLBJ1, flight_time_CPER1, 
                   flight_time_KONZ1, flight_time_NIWO1, flight_time_NOGP1, flight_time_OAES1, 
                   flight_time_SJER1, flight_time_WOOD1)
flight_times2 <- c(flight_time_ABBY2, flight_time_CLBJ2, flight_time_CPER2, 
                   flight_time_KONZ2, flight_time_NIWO2, flight_time_NOGP2, flight_time_OAES2, 
                   flight_time_SJER2, flight_time_WOOD2)

# 3) Specify time filter to discard plots that are not surveyed around the time of the flights
filter_name <- "61DaysFilter"
times_start <- flight_times1 - rep(61* 24*60*60, length(flight_times1))
times_end <- flight_times2 + rep(61* 24*60*60, length(flight_times2))

# 4) Vectors with variable names that will be used during analyses
vars <- c("cv_masked_log", "chv_masked_log", "sam_masked_log")
xlab_names <- c("Spectral Coefficient of Variation (log)", "Convex hull volume (log)", "Spectral angle measure (log)")

yvars <- c("SR","Shannon_exp", "Simpson_inv")
ylab_names <- c("Species Richness", "effective Shannon entropy", "inverse Simpson concentration")

# 5) Create base name for saving results
filename_base_all <- paste0("All_sites_",aop_data_year,"_",ndvi_tresh,"ndviTresh_",filter_name)


#################################################################################
# PART 1: Read and prepare data -------
#################################################################################
## ---- Read in plot level spcomp data and spectral diversity, and calculate taxonomic diversity --------------
# Concatenate info of different sites
All_veg_plots_base_sdiv_tdiv_list <- list()
list_species <- as.data.frame(matrix(nrow=0, ncol=5))
names(list_species) <- c("siteID", "taxonID", "family", "nativeStatusCode", "scientificName")

library(betapart)
library(textshape)

for (i in 1:length(site_codes)){
  print(paste0("Reading data from ", site_codes[i]))
  
  ## ---- Open baseplot level data  --------------
  veg_plots_base_data <- open_baseplot_data(dir_data_out, site_codes[i], aop_data_year, ndvi_tresh)
  
  ## ---- Prepare the baseplot level data  --------------
  # This includes concatenating all information, but also dealing with duplicate field surveys
  veg_plots_base_data_prep <- preparing_baseplot_data(veg_plots_base_data, times_start[i], times_end[i])
  veg_plots_base_spComp <- veg_plots_base_data_prep$veg_plots_base_species_abund_wide
  
  ## ---- Calculate abundance-based taxonomic diversity for each base plot --------------
  veg_plots_base_tdiv <- Calculate_plot_taxon_div(veg_plots_base_spComp)
  
  ## ---- Calculate species richness for each base plot --------------
  # Total species richness of base plot = species presence in subplots + rest of the plot
  plot_SR <- veg_plots_base_data_prep$veg_plots_base_species_pres %>%
    dplyr::group_by(plotID) %>%
    dplyr::summarise(SR = n_distinct(taxonID))
  
  # Species richness calculated as the total richness of all 1 m² subplots only
  plot_SR_incomplete <- veg_plots_base_data$veg_plots_spcomp_site_year$table_observation_species_subplots %>%
    dplyr::group_by(plotID) %>%
    dplyr::summarise(SR_incomplete = n_distinct(taxonID))

  
  veg_plots_base_tdiv <- veg_plots_base_tdiv %>%
    as.data.frame() %>%
    inner_join(plot_SR, by=c("plotID")) %>%
    inner_join(plot_SR_incomplete, by=c("plotID"))

  ## ---- Calculate species evenness for each base plot --------------
  veg_plots_base_tdiv <- veg_plots_base_tdiv %>%
    mutate(Shannon_exp = exp(Shannon),
           Evenness_Pielou = Shannon/log(SR_incomplete)) # I used SR_incomplete here, to match the dataset used to calculate Shannon
  
  # Pielou's evenness = Shannon/Shannon max (https://www.rpubs.com/roalle/mres_2019)
  
  ## ---- Calculate subplot dissimilarity for each base plot  --------------
  dissim <- as.data.frame(matrix(nrow=length(unique(veg_plots_base_spComp$plotID)), ncol=2))
  names(dissim) <- c("plotID", "dissim_Bray")
  
  for (j in 1:length(unique(veg_plots_base_spComp$plotID))){
    table_observation_data_baseplot <- veg_plots_base_spComp %>%
      dplyr::filter(plotID == unique(veg_plots_base_spComp$plotID)[j]) %>%
      dplyr::select(-siteID, -plotID, -subplotID, -namedLocation, -elevation, -decimalLatitude, -decimalLongitude, -month) %>%
      dplyr::group_by(plot_subplotID) %>%
      summarize_all(sum) %>%
      dplyr::select(where(~ any(. != 0))) %>%
      tibble::column_to_rownames(var = "plot_subplotID")
    Bray_1m <- vegan::vegdist(table_observation_data_baseplot, index="bray")
    
    dissim[j,"plotID"] <- unique(veg_plots_base_spComp$plotID)[j]
    dissim[j,"dissim_Bray"] <- mean(Bray_1m)
  }
  
  veg_plots_base_tdiv <- veg_plots_base_tdiv %>%
    inner_join(dissim, by=c("plotID"))
  
  ## ---- Calculate alien species richness (based on presence data) for each base plot --------------
  plot_SR_I <- veg_plots_base_data_prep$veg_plots_base_species_pres %>%
    dplyr::filter(nativeStatusCode == "I") %>%
    dplyr::group_by(plotID) %>%
    dplyr::summarise(SR_I = n_distinct(taxonID))
  
  plot_SR_NI <- veg_plots_base_data_prep$veg_plots_base_species_pres %>%
    dplyr::filter(nativeStatusCode == "NI") %>%
    dplyr::group_by(plotID) %>%
    dplyr::summarise(SR_NI = n_distinct(taxonID))
  
  plot_SR_N <- veg_plots_base_data_prep$veg_plots_base_species_pres %>%
    dplyr::filter(nativeStatusCode == "N") %>%
    dplyr::group_by(plotID) %>%
    dplyr::summarise(SR_N = n_distinct(taxonID))
  
  plot_SR_UNK <- veg_plots_base_data_prep$veg_plots_base_species_pres %>%
    dplyr::filter(nativeStatusCode == "UNK") %>%
    dplyr::group_by(plotID) %>%
    dplyr::summarise(SR_UNK = n_distinct(taxonID))
  
  veg_plots_base_tdiv <- veg_plots_base_tdiv %>%
    left_join(plot_SR_I, by=c("plotID")) %>%
    left_join(plot_SR_NI, by=c("plotID")) %>%
    left_join(plot_SR_N, by=c("plotID")) %>%
    left_join(plot_SR_UNK, by=c("plotID")) %>%
    mutate(SR_I_NI = rowSums(.[c("SR_I", "SR_NI")], na.rm = T)) %>%
    mutate(SR_sum = rowSums(.[c("SR_I", "SR_NI", "SR_N","SR_UNK")], na.rm = T)) %>% # checking if sum corresponds to SR
    mutate(SR_I_NI_perc = SR_I_NI / SR)
  
  veg_plots_base_tdiv[is.na(veg_plots_base_tdiv)] <- 0
  
  ## ---- Calculate IAS vegetation cover for each base plot --------------
  # Cover of all species
  PlotBySpecies <- veg_plots_base_spComp %>% 
    dplyr::select(-siteID, -subplotID, -plot_subplotID, -namedLocation, -elevation, -decimalLatitude, -decimalLongitude, -month) %>%
    dplyr::group_by(plotID) %>%
    summarize_all(sum) %>%
    tibble::column_to_rownames(var = "plotID")
  totalSpeciesCover <- apply(PlotBySpecies, 1, sum)
  
  # Cover of all invasive alien species
  Invasives_pres <- veg_plots_base_data_prep$veg_plots_base_species_pres %>%
    dplyr::filter(nativeStatusCode == "I" | nativeStatusCode == "NI") %>%
    dplyr::select(taxonID) %>%
    unique()
  print(paste0("number of invasive species (I + NI) in presence dataset is ", dim(Invasives_pres)[1]))
  
  Invasives_abund <- Invasives_pres %>%
    dplyr::filter(taxonID %in% names(PlotBySpecies))
  print(paste0("number of invasive species (I + NI) in abundance dataset is ", dim(Invasives_abund)[1]))
  
  Invasives_abund_vector <- paste(Invasives_abund[[1]],sep=",")
  PlotBySpecies_Invasives <- PlotBySpecies %>% dplyr::select(all_of(Invasives_abund_vector))
  
  InvasiveSpeciesCover <- apply(PlotBySpecies_Invasives, 1, sum)
  InvasiveSpeciesCover_perc <- InvasiveSpeciesCover/totalSpeciesCover*100  
  InvasiveSpeciesCover_perc_df <- InvasiveSpeciesCover_perc %>%
    as.data.frame() %>%
    dplyr::rename(cover_I_NI_perc = ".") %>%
    tibble::rownames_to_column(var = "plotID")
  
  veg_plots_base_tdiv <- veg_plots_base_tdiv %>%
    left_join(InvasiveSpeciesCover_perc_df, by=c("plotID"))
  
  ## ---- Calculate soil cover for each base plot --------------
  veg_plots_base_soil_abund <- veg_plots_base_data_prep$veg_plots_base_dead %>%
    dplyr::filter(otherVariables  == "soil") %>%
    dplyr::select(-siteID, -subplotID, -plot_subplotID, -namedLocation, -elevation, -decimalLatitude, -decimalLongitude, -month,
                  -sampleID, -endDate, -nlcdClass, -otherVariables) %>%
    group_by(plotID) %>%
    summarize_all(mean, na.rm = TRUE) %>%
    dplyr::rename("soil_perc" = "percentCover")
  
  veg_plots_base_tdiv <- veg_plots_base_tdiv %>%
    left_join(veg_plots_base_soil_abund, by=c("plotID")) %>%
    mutate(soil_perc = tidyr::replace_na(soil_perc,0)) # NA in this column actually means zero

  ## ---- Concatenate information on tdiv and sdiv --------------
  veg_plots_base_sdiv_tdiv <- veg_plots_base_data_prep$veg_plots_base_sdiv %>%
    inner_join(veg_plots_base_tdiv, c("plotID"))
  
  ## ---- Store the results in the list --------------
  All_veg_plots_base_sdiv_tdiv_list[[i]] <- veg_plots_base_sdiv_tdiv
  
  ## ---- Add species to species list --------------
  list_species_site <- veg_plots_base_data_prep$veg_plots_base_species_pres %>%
    dplyr::select(siteID, taxonID, family, nativeStatusCode) %>%
    unique() %>%
    left_join(veg_plots_base_data$veg_plots_spcomp_site_year$table_taxon %>%
                dplyr::select(taxonID, scientificName) %>%
                unique(), by=c("taxonID"))
  
  list_species <- rbind(list_species, list_species_site)
  
}

All_veg_plots_base_sdiv_tdiv_df <- ldply(All_veg_plots_base_sdiv_tdiv_list, rbind)

list_species
list_species_unique <- list_species %>% 
  dplyr::filter(taxonID != "2PLANT") %>%
  dplyr::select(-siteID, -nativeStatusCode)
dim(list_species_unique)

## ---- Transform data --------------
df <- All_veg_plots_base_sdiv_tdiv_df

df <- df %>%
  mutate(siteID = as.factor(siteID)) %>%
  mutate(
    cv_masked_log = log(cv_masked),
    chv_masked_log = log(chv_masked),
    sam_masked_log = log(sam_masked),
  
    stat_herbaceousNPP = 0.1 * stat_herbaceousNPP, # convert to g m-2 year-1
    # The rationale for this conversion is based on GEE code in
    # https://code.earthengine.google.com/312894322045efc8d0b9bb89af166816
    # They use a scalar of 0.0001 to convert raw values into kg m ² year-1
    # This implies that a scalar of 0.1 is needed to convert raw values into g m-2 year-1
  )

## ---- Remove unwanted data points --------------
# Most baseplots have 8 subsamples, except for some that we will discard
graphics.off()
hist(df$n_subplots)
discarded_plots <- df %>%
  dplyr::filter(n_subplots != 8)
df <- df %>%
  dplyr::filter(n_subplots == 8)

# Only retain base plots that are listed as "grasslandHerbaceous"
discarded_plots2 <- df %>%
  dplyr::filter(nlcdClass != "grasslandHerbaceous")
df <- df %>%
  dplyr::filter(nlcdClass == "grasslandHerbaceous")

df %>%
  group_by(siteID) %>%
  dplyr::summarise(count = n())

#################################################################################
## PART 2: Explore data -------
#################################################################################
## ---- Check data distribution --------------
par(mfrow=c(3,2))
hist(df$cv_masked, xlab = 'cv_masked', main = NA)
hist(df$chv_masked, xlab = 'chv_masked', main = NA)
hist(df$sam_masked, xlab = 'sam_masked', main = NA)
hist(df$cv_masked_log, xlab = 'cv_masked_log', main = NA)
hist(df$chv_masked_log, xlab = 'chv_masked_log', main = NA)
hist(df$sam_masked_log, xlab = 'sam_masked_log', main = NA)

par(mfrow=c(1,3))
hist(df$SR, xlab = 'SR', main = NA)
hist(df$Shannon_exp, xlab = 'Shannon_exp', main = NA)
hist(df$Simpson_inv, xlab = 'Simpson_inv', main = NA)

par(mfrow=c(3,2))
hist(df$lai, xlab = 'LAI', main = NA)
hist(df$NDVI, xlab = 'NDVI', main = NA)
hist(df$stab_herbaceousNPP, xlab = 'Herbaceous spp NPP stability', main = NA)
hist(df$soil_perc, xlab = 'Soil cover', main = NA)
hist(df$SR_I_NI, xlab = 'SR_I_NI', main = NA)
hist(df$SR_I_NI_perc, xlab = 'percentage SR_I_NI', main = NA)

graphics.off()


## ---- Create and save explorative plots for xvars and yvars --------------
vis_vars <- c(xvars, yvars)
vis_names <- c(xlab_names, ylab_names)

for (j in 1:length(vis_vars)){
  f1 <- ggplot(df, aes_string(vis_vars[j], fill = "siteID")) + 
    geom_histogram(binwidth = 1) + 
    facet_grid(siteID ~., margins = TRUE, scales = "free") +
    labs(x =vis_names[j]) +
    ggsci::scale_fill_igv() +
    theme_bw() +
    theme(legend.position = "none") 
  f1
  # If the variances within each level are higher than the means within each level,
  # this suggests that over-dispersion is present and that a Negative Binomial model would be appropriate.
  
  Fig_name_f1 <- paste0(dir_results_explorative, "/", filename_base_all, "_", vis_vars[j], "_histogram", ".png")
  ggsave(Fig_name_f1, f1, width = 4, height = 6, dpi = 600)
  
  f2 <- ggplot(df, aes_string(vis_vars[j], "siteID", fill = "siteID")) + 
    geom_boxplot() +
    coord_flip() +
    labs(x = vis_names[j], y="") +
    ggsci::scale_fill_igv() +
    theme_bw() +
    theme(legend.position = "none") 
  f2
  
  Fig_name_f2 <- paste0(dir_results_explorative, "/", filename_base_all, "_", vis_vars[j] , "_boxplots", ".png")
  ggsave(Fig_name_f2, f2, width = 4, height = 4, dpi = 600)
}

## ---- Generate data summaries for xvars and yvars --------------
print("Mean and SD of species richness per site")
with(df, tapply(SR, siteID, function(x) {
  sprintf("M (SD) = %1.2f (%1.2f)", mean(x), sd(x))
}))

print("Mean and SD of SR_I_NI_Perc per site")
with(df, tapply(SR_I_NI_perc, siteID, function(x) {
  sprintf("M (SD) = %1.2f (%1.2f)", mean(x), sd(x))
}))

summary(df$stat_herbaceousNPP)


## ---- Calculate total number of species per site, based on plots used in our study --------------
for (i in 1:length(site_codes)){
  print(paste0("Reading data from ", site_codes[i]))
  
  ## ---- Open baseplot level data  --------------
  veg_plots_base_data <- open_baseplot_data(dir_data_out, site_codes[i], aop_data_year, ndvi_tresh)
  
  ## ---- Prepare the baseplot level data  --------------
  # This includes concatenating all information, but also dealing with duplicate field surveys
  veg_plots_base_data_prep <- preparing_baseplot_data(veg_plots_base_data, times_start[i], times_end[i])
  veg_plots_base_spComp <- veg_plots_base_data_prep$veg_plots_base_species_abund_wide
  
  ## ---- Remove plots that were also removed earlier  --------------
  `%!in%` <- Negate(`%in%`)

  veg_plots_base_species_pres_filtered <- veg_plots_base_data_prep$veg_plots_base_species_pres %>%
    dplyr::filter(plotID %!in% discarded_plots[,"plotID"]) %>%
    dplyr::filter(plotID %!in% discarded_plots2[,"plotID"]) 
  
  ## ---- Calculate number of species across all plots  --------------
  present_species <- veg_plots_base_species_pres_filtered %>%
    dplyr::select(taxonID) %>%
    unique()
  print(paste0("Number of species in grassland/herbaceous plots at ", site_codes[i], " is ", nrow(present_species)))
  
  ## ---- Calculate number of invasive species across all plots  --------------
  species_I_NI <- veg_plots_base_data$veg_plots_spcomp_site_year$table_taxon %>%
    dplyr::filter(nativeStatusCode == "I" | nativeStatusCode == "NI") %>%
    dplyr::select(taxonID) %>%
    unique()
  
  present_species_I_NI <- present_species %>%
    dplyr::filter(taxonID %in% species_I_NI[,"taxonID"])
  
  # Same result as:
  # species_I_NI_pres <- veg_plots_base_data_prep$veg_plots_base_species_pres %>%
  #   dplyr::filter(nativeStatusCode == "I" | nativeStatusCode == "NI") %>%
  #   dplyr::select(taxonID) %>%
  #   unique()
  # present_species_I_NI <- present_species %>%
  #   dplyr::filter(taxonID %in% species_I_NI_pres[,"taxonID"])
  
  print(paste0("Number of I/NI species in grassland/herbaceous plots at ", site_codes[i], " is ", nrow(present_species_I_NI)))
  
}

#################################################################################
## PART 3: Within and across sites - linear models: Taxonomic vs. spectral diversity -------
#################################################################################
# Fit linear model for each site, linking tdiv to sdiv
# AND
# Fit linear model for all sites combined (denoted as "global" in the code), linking tdiv to sdiv

# Create empty dataframes to store the results in
tdiv_sdiv_intercept <- data.frame(matrix(NA, nrow=length(xvars)*length(yvars), ncol=length(site_codes)))
variable_combinations <- expand.grid(yvars,"_",xvars)
rownames(tdiv_sdiv_intercept) <- do.call(paste0,variable_combinations[order(variable_combinations$Var1),])
colnames(tdiv_sdiv_intercept) <- site_codes
tdiv_sdiv_slope <- tdiv_sdiv_intercept
tdiv_sdiv_pvalue <- tdiv_sdiv_intercept
tdiv_sdiv_n <- tdiv_sdiv_intercept
tdiv_sdiv_R2 <- tdiv_sdiv_intercept
tdiv_sdiv_RMSE <- tdiv_sdiv_intercept

tdiv_sdiv_global <- data.frame(matrix(NA, nrow=length(xvars)*length(yvars), ncol=4))
variable_combinations <- expand.grid(yvars,"_",xvars)
rownames(tdiv_sdiv_global) <- do.call(paste0,variable_combinations[order(variable_combinations$Var1),])
colnames(tdiv_sdiv_global) <- c("R2", "RMSE", "slope", "p")

# Function used during visualization of site-specific slopes
gather_data_for_slopes <- function(tdiv_sdiv_models, site_codes){
  y_predicted_list <- lapply(tdiv_sdiv_models, function(x) predict(x))
  x_values_list <- lapply(tdiv_sdiv_models, function(x) x$model[[xvars[i]]])
  
  y_predicted_df = melt(setNames(y_predicted_list, site_codes))
  names(y_predicted_df) <- c("y_pred", "siteID")
  x_values_df = melt(setNames(x_values_list, site_codes))
  names(x_values_df) <- c("x_value", "siteID2")
  
  tdiv_sdiv_predictions <- cbind(y_predicted_df, x_values_df)
}

# Fit models and visualize results
for (j in 1:length(yvars)){ 
  for (i in 1:length(xvars)){
    
    variable_combination <- paste0(yvars[j],"_",xvars[i])
    
    ### ------------- Fit site-specific tdiv-sdiv relationships and extract information on model fits  --------------
    tdiv_sdiv_formula <- as.formula(paste0(yvars[j]," ~ ",xvars[i]))
    tdiv_sdiv_models = lapply(split(df,df$siteID),function(s){lm(tdiv_sdiv_formula, data=s)})
    
    tdiv_sdiv_slope[variable_combination,] <- sapply(tdiv_sdiv_models, function(x) round(summary(x)$coefficients[2,"Estimate"],2))
    tdiv_sdiv_pvalue[variable_combination,] <- sapply(tdiv_sdiv_models, function(x) round(summary(x)$coefficients[2,"Pr(>|t|)"],3))
    tdiv_sdiv_n[variable_combination,] <- sapply(tdiv_sdiv_models, function(x) nrow(x$model))
    tdiv_sdiv_R2[variable_combination,] <- sapply(tdiv_sdiv_models, function(x) round(summary(x)$r.squared,3))
    tdiv_sdiv_RMSE[variable_combination,] <- sapply(tdiv_sdiv_models, function(x) sqrt(mean(x$residuals^2)))

    ### ------------- Fit global tdiv-sdiv relationship and extract information on --------------
    global_model = lm(as.formula(paste0(yvars[j]," ~ ",xvars[i])), data=df)
    tdiv_sdiv_global[variable_combination, "R2"] <- round(summary(global_model)$r.squared, 3)
    tdiv_sdiv_global[variable_combination, "RMSE"] <- round(sqrt(mean(global_model$residuals^2)),2)
    tdiv_sdiv_global[variable_combination, "slope"] <- round(summary(global_model)$coefficients[2,"Estimate"],3)
    tdiv_sdiv_global[variable_combination, "p"] <- round(summary(global_model)$coefficients[2,"Pr(>|t|)"],5)
    
    ### ------------- Visualize results --------------
    colors =  data.frame("siteID" = as.factor(levels(df[,"siteID"])),
                         "site_color" = as.factor(ggsci::pal_igv("default")(nlevels(df[,"siteID"]))))
    df2 <- df %>%
      left_join(colors, by=c("siteID"))
    
    # 1) Plot global model
    g <- ggplot() +
      geom_smooth(aes_string(x=xvars[i], y=yvars[j]),data=df2,
                  method="lm", size = 1, alpha = .25, color="black") +
      geom_point(aes_string(x=xvars[i], y=yvars[j],color="siteID"), data=df2, size=2) +
      ggsci::scale_colour_igv() +
      xlab("Spectral diversity")+
      ylab(ylab_names[j])+
      theme_classic() +
      theme(axis.text=element_text(size=12),
            axis.title=element_text(size=14),
            legend.title = element_blank(),
            legend.text = element_text(size=12))
    # print(g)
    Fig_name_g <- paste0(dir_results_tdiv_vs_sdiv_general, "/", filename_base_all, "_", yvars[j], "_", xvars[i], "_global.png")
    ggsave(Fig_name_g, g, width = 5, height = 4, dpi = 600)
    
    # 2) plot site-specific slopes
    p <- ggplot() +
      geom_point(aes_string(x=xvars[i], y=yvars[j],color="siteID"), data=df2, size=2) +
      ggsci::scale_colour_igv() +
      xlab("Spectral diversity")+
      ylab(ylab_names[j])+
      theme_classic() +
      theme(axis.text=element_text(size=12),
            axis.title=element_text(size=14),
            legend.title = element_blank(),
            legend.text = element_text(size=12))
    
    # Visualize significant slopes with a solid line and non significant slopes with a dashed line
    notSignSites <- which(tdiv_sdiv_pvalue[variable_combination, ] > 0.1)
    SignSites <- which(tdiv_sdiv_pvalue[variable_combination, ] <= 0.1)
    tdiv_sdiv_models_notSign <- tdiv_sdiv_models[notSignSites]
    tdiv_sdiv_models_Sign <- tdiv_sdiv_models[SignSites]
    site_codes_notSign <- site_codes[notSignSites]
    site_codes_Sign <- site_codes[SignSites]
    
    if (length(site_codes_notSign) > 1){
      tdiv_sdiv_predictions_notSign <- gather_data_for_slopes(tdiv_sdiv_models_notSign, site_codes_notSign)
      p1 <- p +
        geom_line(data = tdiv_sdiv_predictions_notSign, aes(x=x_value, y=y_pred, group=siteID, color=factor(siteID)),
                  size = 1, linetype = "dashed")
    } else {
      p1 <- p
    }
    
    if (length(site_codes_Sign) > 1){
      tdiv_sdiv_predictions_Sign <- gather_data_for_slopes(tdiv_sdiv_models_Sign, site_codes_Sign)
      p2 <- p1 +
        geom_line(data = tdiv_sdiv_predictions_Sign, aes(x=x_value, y=y_pred, group=siteID, color=factor(siteID)),
                  size = 1, linetype = "solid")
    } else {
      p2 <- p1
    }
    print(p2)
    
    Fig_name_p2 <- paste0(dir_results_tdiv_vs_sdiv_general, "/", filename_base_all, "_", yvars[j], "_", xvars[i], "_siteSpecific_solidAndDashed.png")
    ggsave(Fig_name_p2, p2, width = 5, height = 4, dpi = 600)
    
    
    
  }
}

# Save regression models
results_global_name <- paste0(dir_results_tdiv_vs_sdiv_general, "/", filename_base_all, "_tdiv_vs_sdiv_global_results", ".csv")
write.csv(tdiv_sdiv_global, results_global_name)

regression_results_siteSpecific <- list("slope" = tdiv_sdiv_slope,
                                        "pvalue" = tdiv_sdiv_pvalue,
                                        "n_plots" = tdiv_sdiv_n,
                                        "R2" = tdiv_sdiv_R2, 
                                        "RMSE" = tdiv_sdiv_RMSE)
                                        
results_name_siteSpecific <- paste0(dir_results_tdiv_vs_sdiv_general, "/", filename_base_all, "_tdiv_vs_sdiv_results_siteSpecific", ".rds")
saveRDS(regression_results_siteSpecific, results_name_siteSpecific)


# Concatenate results into 1 table for Supplementary Information
if (!require('purrr')) install.packages('purrr'); library('purrr')
if (!require('stringr')) install.packages('stringr'); library('stringr')
part1 <- map2_dfc("R² = ",
                  round(regression_results_siteSpecific$R2, 2),
                  ~ str_c(.x, .y, sep = ""))
part2 <- map2_dfc(" (slope = ",
                  regression_results_siteSpecific$slope,
                  ~ str_c(.x, .y, sep = ""))
part3 <- map2_dfc(", p = ",
                  round(regression_results_siteSpecific$pvalue,2),
                  ~ str_c(.x, .y, sep = ""))
part12 <- map2_dfc(part1, part2,  ~ str_c(.x, .y, sep = ""))
part123 <- map2_dfc(part12, part3,  ~ str_c(.x, .y, sep = ""))
final <- map2_dfc(part123, ")",  ~ str_c(.x, .y, sep = "")) %>%
  as.data.frame()
rownames(final) <- rownames(regression_results_siteSpecific$R2)
colnames(final) <- colnames(regression_results_siteSpecific$R2)


results_together_name <- paste0(dir_results_tdiv_vs_sdiv_general, "/", filename_base_all, "_tdiv_vs_sdiv_results_siteSpecific", ".csv")
write.csv(final, results_together_name)

#################################################################################
## PART 4: Across sites - linear model: Taxonomic vs. spectral diversity with site interaction-------
#################################################################################
# Fit linear model for all sites combined, linking tdiv to sdiv whereby slopes can vary with site (interaction)

if (!require('jtools')) install.packages('jtools'); library('jtools')
# https://cran.r-project.org/web/packages/interactions/vignettes/interactions.html

# Create empty dataframes to store the results in
regression_results <- data.frame(matrix(NA, nrow=length(xvars)*length(yvars), ncol=5))
variable_combinations <- expand.grid(yvars,"_",xvars)
rownames(regression_results) <- do.call(paste0,variable_combinations[order(variable_combinations$Var1),])
colnames(regression_results) <- c("R2", "RMSE", "p xvar","p siteID", "p interaction")

# Fit models and visualize results
for (j in 1:length(yvars)){ 
  for (i in 1:length(xvars)){
    
    lm_yvar_xvar <- lm(create_lm_formula(yvars[j], xvars[i], "siteID"), data=df)
    model_name <- paste0(dir_results_tdiv_vs_sdiv_general, "/", filename_base_all, "_", yvars[j], "_", xvars[i], ".rds")
    saveRDS(lm_yvar_xvar, model_name)
    
    variables <- paste0(yvars[j],"_",xvars[i])
    regression_results[variables,"R2"] <- round(summary(lm_yvar_xvar)$r.squared, 2)
    regression_results[variables,"RMSE"] <- round(sqrt(mean(lm_yvar_xvar$residuals^2)), 2) # sjstats::rmse(lm_yvar_xvar)
    regression_results[variables,"p xvar"] <- round(as.data.frame(anova(lm_yvar_xvar))[xvars[i],"Pr(>F)"], 3)
    regression_results[variables,"p siteID"] <- round(as.data.frame(anova(lm_yvar_xvar))["siteID","Pr(>F)"], 3)
    regression_results[variables,"p interaction"] <- round(as.data.frame(anova(lm_yvar_xvar))[paste0(xvars[i],":siteID"),"Pr(>F)"], 3)
    
    summary(lm_yvar_xvar)
    print("slope varies according to the site")

    colors =  data.frame("siteID" = as.factor(levels(df[,"siteID"])), 
                         "site_color" = as.factor(ggsci::pal_igv("default")(nlevels(df[,"siteID"]))))
    df2 <- df %>%
      left_join(colors, by=c("siteID")) %>%
      dplyr::select(siteID, xvars[i], yvars[j]) %>%
      na.omit()
    
    p <- ggplot(df2, aes_string(x=xvars[i], y=yvars[j], group="siteID", color="siteID", fill="siteID")) +
      geom_smooth(method="lm", size = 2, alpha = .25) +
      geom_line(aes(y=predict(lm_yvar_xvar), group=df2[,"siteID"]), size=2) +
      geom_point(size=1) +
      ggsci::scale_colour_igv() +
      ggsci::scale_fill_igv() +
      xlab(xlab_names[i])+
      ylab(ylab_names[j])+
      theme_classic()
    print(p)
    
    Fig_name_p <- paste0(dir_results_tdiv_vs_sdiv_general_Sitescombined, "/", filename_base_all, "_", yvars[j], "_", xvars[i], ".png")
    ggsave(Fig_name_p, p, width = 5, height = 4, dpi = 600)
    
    Fig_name_p_assump <- paste0(dir_results_tdiv_vs_sdiv_general_Sitescombined, "/", filename_base_all, "_", yvars[j], "_", xvars[i], "_assump.png")
    png(Fig_name_p_assump, width=700, height=700)
    par(mfrow=c(2,2))
    plot(lm_yvar_xvar)
    dev.off()
    
  }
}  

results_name <- paste0(dir_results_tdiv_vs_sdiv_general_Sitescombined, "/", filename_base_all, "_tdiv_vs_sdiv_results", ".csv")
write.csv(regression_results, results_name)

#################################################################################
## PART 5: Calculate species turnover -------
#################################################################################
# Calculate species turnover from species-area relationships

if (!require('utils')) install.packages('utils'); library('utils')
if (!require('sars')) install.packages('sars'); library('sars')

# Create empty lists to store the results in
All_veg_plots_base_sar_list_1m <- list()

## ---- Calculate area and number of species for all subplot combinations --------------
for (i in 1:length(site_codes)){
  ## ---- Open baseplot level data  --------------
  veg_plots_base_data <- open_baseplot_data(dir_data_out, site_codes[i], aop_data_year, ndvi_tresh)
  
  ## ---- Prepare the baseplot level data  --------------
  # This includes concatenating all information, but also dealing with duplicate field surveys
  veg_plots_base_data_prep <- preparing_baseplot_data(veg_plots_base_data, times_start[i], times_end[i])
  veg_plots_base_spPres <- veg_plots_base_data_prep$veg_plots_base_species_pres
  
  # Organize information per 1m² subplot
  veg_plots_base_spPres_1m <- veg_plots_base_spPres %>%
    dplyr::filter(subplotID != "31" & subplotID != "32" & subplotID != "40" & subplotID != "41" &
                    subplotID != "31.1.1" & subplotID != "31.4.1" & subplotID != "32.2.1" & subplotID != "32.4.1" & 
                    subplotID != "40.1.1"& subplotID != "40.3.1"& subplotID != "41.1.1" & subplotID != "41.4.1") %>% # remove areas sampled between subplots
    dplyr::select(-sampleID, -plot_subplotID) %>%
    unique()
  
  ## ---- Calculate species richness for different combinations of subplots  --------------
  veg_plots_base_sar_1m <- data.frame()
  
  # Loop through all combinations of subplots
  for (s in 1:length(unique(veg_plots_base_spPres_1m$subplotID))) { 
    subplot_combinations <- utils::combn(unique(veg_plots_base_spPres_1m$subplotID), m=s)
    dim(subplot_combinations)[2]
    
    for (c in 1:dim(subplot_combinations)[2]){
      veg_plots_base_sar_combi <- veg_plots_base_spPres_1m %>%
        dplyr::filter(subplotID %in% subplot_combinations[,c]) %>%
        dplyr::group_by(plotID) %>%
        dplyr::summarise(SR = n_distinct(taxonID)) %>%
        mutate(n_plots = s,
               area = s, # 1 subplot is 1m²
               plots = paste(subplot_combinations[,c], collapse="_"),
               siteID = site_codes[i])  
      
      veg_plots_base_sar_1m <- rbind(veg_plots_base_sar_1m, veg_plots_base_sar_combi)
    }
  }

  ## ---- Store the results in the list --------------
  All_veg_plots_base_sar_list_1m[[i]] <- veg_plots_base_sar_1m

}

All_veg_plots_base_sar_df_1m <- ldply(All_veg_plots_base_sar_list_1m, rbind)

## ---- Remove plots that were also removed earlier (PART 1) --------------
`%!in%` <- Negate(`%in%`)
All_veg_plots_base_sar_df_1m <- All_veg_plots_base_sar_df_1m %>%
  dplyr::filter(plotID %!in% discarded_plots[,"plotID"]) %>%
  dplyr::filter(plotID %!in% discarded_plots2[,"plotID"]) 

## ---- Fit and save Species richness - area relationships (SAR) for each base plot --------------
print("Fit Species richness - area relationships (SAR) for each base plot")
sar_fit_results_1m <- list()
sar_models_1m <- list()

for (i in 1:length(site_codes)){
  All_veg_plots_base_sar_site <- All_veg_plots_base_sar_df_1m %>%
    dplyr::filter(siteID == site_codes[i]) 
  
  sar_fit_results_site <- data.frame(matrix(ncol = 5, nrow = 0))
  colnames(sar_fit_results_site) <- c("siteID", "plotID", 
                                      "log_c", "log_z", 
                                      "log_R2")
  sar_models_site <- list()
  
  for (j in 1:length(unique(All_veg_plots_base_sar_site$plotID))){
    All_veg_plots_base_sar_site_plot <- All_veg_plots_base_sar_site %>%
      dplyr::filter(plotID == unique(All_veg_plots_base_sar_site$plotID)[j]) %>%
      dplyr::select(area, SR)
    
    sar_fit_results_site[j,"siteID"] <- site_codes[i]
    sar_fit_results_site[j,"plotID"] <- unique(All_veg_plots_base_sar_site$plotID)[j]
    
    sar_fit_log <-  sar_loga(All_veg_plots_base_sar_site_plot) # S == c + z * log(A)
    sar_fit_results_site[j, c("log_c", "log_z")] <-  sar_fit_log$par
    sar_fit_results_site[j, "log_R2"] <- sar_fit_log$R2
   
    site_plot_sar_models <- list("sar_fit_log" = sar_fit_log)
    sar_models_site[[j]] <- site_plot_sar_models
    
  }
  
  sar_fit_results_1m[[i]] <- sar_fit_results_site
  names(sar_models_site) <- unique(All_veg_plots_base_sar_site$plotID)
  sar_models_1m[[i]] <- sar_models_site
  
} 
names(sar_fit_results_1m) <- site_codes
names(sar_models_1m) <- site_codes

patchiness_baseplots_results_name_1m <- paste0(dir_results_patchiness, "/", filename_base_all, "_SAR_baseplots_1m.rds")
saveRDS(sar_fit_results_1m, patchiness_baseplots_results_name_1m)

## ---- Visualize and save base plot-specific Species richness - area relationships (SAR)--------------
print("Plot base plot-specific Species richness - area relationships")

data_log_all_sites_1m <- data.frame(matrix(nrow=0, ncol=5))
names(data_log_all_sites_1m) <- c("Area", "SR", "SR_pred", "siteID", "plotID")

for(i in 1:length(site_codes)){
  
  sar_models_site <- sar_models_1m[[i]]
  
  data_log_all_basesplots <- data_log_all_sites_1m

  for (j in 1:length(names(sar_models_site))){
    sar_models_site[[j]]
    
    data_log <- sar_models_site[[j]]$sar_fit_log$data %>%
      mutate(pred = sar_models_site[[j]]$sar_fit_log$calculated) %>%
      mutate(siteID = site_codes[i]) %>%
      mutate(plotID = names(sar_models_site)[j])
    names(data_log) <- c("Area", "SR", "SR_pred", "siteID", "plotID")
    data_log_all_basesplots <- rbind(data_log_all_basesplots, data_log)
  }
  
  data_log_all_sites_1m <- data_log_all_basesplots

}  

p_log_1m <- ggplot(data_log_all_sites_1m, aes(x=Area, y=SR, col=siteID)) +
  geom_point(size=0.5) +
  facet_wrap(~siteID) +
  geom_line(aes(x=Area, y=SR_pred, group=plotID), col="black") +
  labs(y = "Species richness", x = "Area (m²)") + 
  ggsci::scale_colour_igv() +
  theme_classic() +
  theme(legend.position = "none")
p_log_1m


fig_name_log_1m <- paste0(dir_results_patchiness, "/", filename_base_all, "_SAR_log_baseplots_1m.png")
ggsave(fig_name_log_1m, p_log_1m, width = 6, height = 5, dpi = 600)


# Calculate number of subplot combinations per site, this should be the same across sites
subplot_combis_nb <- All_veg_plots_base_sar_df_1m %>%
  dplyr::select(siteID, plots) %>%
  dplyr::group_by(siteID) %>%
  unique() %>%
  dplyr::group_by(siteID) %>%
  dplyr::summarise(count=n())
subplot_combis_nb

subplot_combis <- All_veg_plots_base_sar_df_1m %>%
  dplyr::select(siteID, plots) %>%
  dplyr::group_by(siteID) %>%
  unique()
# subplot_combis$plots

#################################################################################
## PART 6: Evaluate relationship between site-level characteristics and 
#          strength of tdiv-sdiv relationship 
#################################################################################
# This section can be run independently of PARTS 2-5 (because intermediate results are stored), 
# but make sure to read and prepare base data, by running PART 1

if (!require('ggrepel')) install.packages('ggrepel'); library('ggrepel')

# Specify explanatory variables: vegetation density, spatial species turnover, and invasion variables 
expl_variables <- c("lai", "NDVI", "stat_herbaceousNPP","soil_perc",
                    "log_z_mean_1m",
                    "SR_I_NI_perc", "cover_I_NI_perc",
                    "Evenness_Pielou", "dissim_Bray")
expl_variables_names <- c("Leaf area index", "NDVI", "Annual and perennial grasses and forbs NPP", "% plot cover soil",
                          "Species turnover",
                          "Non-native proportion", "Non-native proportion (cover)",
                          "Pielou's evenness", "Bray-Curtis dissimilarity")
strength_variables <- c("slope", "R2")

## ---- Read site-specific tdiv-sdiv relationships (generated in PART 3) --------------
regression_results_siteSpecific <- readRDS(paste0(dir_results_tdiv_vs_sdiv_general, "/", filename_base_all, "_tdiv_vs_sdiv_results_siteSpecific", ".rds"))

# Compile metrics of tdiv-sdviv strength in one data frame
df_slopes <- as.data.frame(t(regression_results_siteSpecific$slope))
names(df_slopes) <- paste0(names(df_slopes), ".slope")
df_slopes <- df_slopes %>% tibble::rownames_to_column("siteID")

df_R2 <- as.data.frame(t(regression_results_siteSpecific$R2))
names(df_R2) <- paste0(names(df_R2), ".R2")
df_R2 <- df_R2 %>% tibble::rownames_to_column("siteID")

df_strength <- df_slopes %>%
  inner_join(df_R2, by=c("siteID"))

## ---- Create empty dataframes to store results --------------
expl_regression_R2 <- data.frame(matrix(NA, nrow=ncol(df_strength)-1, ncol=0))
rownames(expl_regression_R2) <- colnames(df_strength)[-c(1)]
expl_regression_slope <- expl_regression_R2
expl_regression_pslope <- expl_regression_R2

# Add columns for explanatory variables to output tables
expl_regression_R2[, expl_variables] <- NA
expl_regression_slope[, expl_variables] <- NA
expl_regression_pslope[, expl_variables] <- NA

## ---- Read plot-level species turnover (SAR slopes; generated in PART 5) and calculate site-specific SAR --------------
sar_fit_results_baseplots_1m <- readRDS(paste0(dir_results_patchiness, "/", filename_base_all, "_SAR_baseplots_1m.rds"))

# Calculate mean SAR slope per site
calc_mean_sar <- function(x){
  apply(x[,-c(1,2)], 2, mean)
}
sar_baseplots_mean_1m <- sapply(sar_fit_results_baseplots_1m, calc_mean_sar) %>%
  t() %>%
  as.data.frame() %>%
  setNames(paste0(names(.),"_mean_1m")) %>%
  tibble::rownames_to_column(var = "siteID")

turnover <- sar_baseplots_mean_1m %>%
  dplyr::select(siteID, log_z_mean_1m)

## ---- Create dataset for analyses --------------
# add explanatory variables to the tdiv-sdiv strength dataframe
veg_density_inv <- df %>%
  dplyr::select(siteID, 
                lai, NDVI, stat_herbaceousNPP, soil_perc,
                SR_I_NI_perc, cover_I_NI_perc,
                Evenness_Pielou, dissim_Bray) %>%
  dplyr::group_by(siteID) %>%
  summarise_if(is.numeric, mean, na.rm = TRUE)

df_strength <- df_strength %>%
  inner_join(turnover, by = "siteID") %>%
  inner_join(veg_density_inv, by = "siteID")

## ---- Evaluate relationship between biotic variables and strength of tdiv-sdiv relationship --------------

for (k in 1:length(expl_variables)){
  for (l in 1:length(strength_variables)){
    for (j in 1:length(yvars)){
      for (i in 1:length(xvars)){
        
        strength_variable <- paste0(yvars[j],"_",xvars[i],".", strength_variables[l])
        model <- lm(as.formula(paste0(strength_variable, " ~ ", expl_variables[k])), data = df_strength)
        summary(model)
        anova(model)
        
        print(paste0(expl_variables[k], " effect on ", strength_variable))
        
        expl_regression_R2[strength_variable,expl_variables[k]] <- round(summary(model)$r.squared,6)
        expl_regression_slope[strength_variable,expl_variables[k]] <- summary(model)$coefficients[expl_variables[k],"Estimate"]
        expl_regression_pslope[strength_variable,expl_variables[k]] <- round(summary(model)$coefficients[expl_variables[k],"Pr(>|t|)"],4)
        
        # Create figure
        if (expl_regression_pslope[strength_variable,expl_variables[k]] <= 0.05) {
          h1 <- ggplot(df_strength, aes_string(x=expl_variables[k], y=strength_variable)) +
            geom_smooth(formula = y ~ x, method="lm", col="black",cex=1) +
            geom_point(aes(col=siteID), cex=4) +
            geom_text_repel(aes_string(label = "siteID"), cex=5.5) +
            # labs(y = paste0("Slope of relation between ", ylab_names[j] , "\n and ", xlab_names[i]),
            #      x = expl_variables_names[k]) +
            labs(y = "Strength of tdiv~sdiv relationship", x = expl_variables_names[k]) +
            ggsci::scale_colour_igv() +
            theme_classic() +
            theme(legend.position = "none",
                  axis.text = element_text(size=14),
                  axis.title = element_text(size=16))
          h1
          
          fig_name1 <- paste0(dir_results_expl_general, "/", filename_base_all, "_", expl_variables[k], "_", strength_variable, "_siteSpecific.png")
          ggsave(fig_name1, h1, width = 5, height = 4, dpi = 600)
        } else {
          h2 <- ggplot(df_strength, aes_string(x=expl_variables[k], y=strength_variable)) +
            # geom_smooth(formula = y ~ x, method="lm", col="black") +
            geom_point(aes(col=siteID),cex=4) +
            geom_text_repel(aes_string(label = "siteID"),cex=5.5) +
            # labs(y = paste0("Strength of relation between ", ylab_names[j] , "\n and ", xlab_names[i]),
            #      x = expl_variables_names[k]) +
            labs(y = "Strength of tdiv~sdiv relationship", x = expl_variables_names[k]) +
            ggsci::scale_colour_igv() +
            theme_classic() +
            theme(legend.position = "none",
                  axis.text = element_text(size=14),
                  axis.title = element_text(size=16))
          h2
          
          fig_name2 <- paste0(dir_results_expl_general, "/", filename_base_all, "_", expl_variables[k], "_",strength_variable, "_notSign_siteSpecific.png")
          ggsave(fig_name2, h2, width = 5, height = 4, dpi = 600)
        }
      } 
    }
  } 
}

# Save regression results
write.csv(expl_regression_R2, paste0(dir_results_expl_general, "/", filename_base_all, "_tdiv_vs_sdiv_expl_R2s", "_siteSpecific.csv"))
write.csv(expl_regression_slope, paste0(dir_results_expl_general, "/", filename_base_all, "_tdiv_vs_sdiv_expl_slope", "_siteSpecific.csv"))
write.csv(expl_regression_pslope, paste0(dir_results_expl_general, "/", filename_base_all, "_tdiv_vs_sdiv_expl_pslope", "_siteSpecific.csv"))

# Save df_strength dataframe
write.csv(df_strength, paste0(dir_results_expl_general, "/", filename_base_all, "_tdiv_vs_sdiv_expl_data_siteSpecific", ".csv"))

#################################################################################
## PART 7: Paired correlation plots and boxplots - 
#################################################################################
# This section can be run independently of PARTS 2-5 (because intermediate results are stored), 
# but make sure to read and prepare base data, by running PART 1

if (!require('corrplot')) install.packages('corrplot'); library('corrplot')

# Specify explanatory variables to be visualized
expl_variables <- c("lai", "NDVI", "stat_herbaceousNPP", "soil_perc",
                    "log_z", 
                    "SR_I_NI_perc", "cover_I_NI_perc",
                    "Evenness_Pielou", "dissim_Bray")
names_expl_variables <- c("LAI", "NDVI", "Herbaceous NPP", "% plot cover soil", 
                          "Species turnover",
                          "% SR non-native", "% cover non-native",
                          "Pielou's evenness",  "Bray-Curtis dissimilarity")

yvars_names_short <- c("Species richness", "Shannon", "Simpson") # short version of ylab_names

variables <- c(yvars,
               xvars,
               expl_variables)
variables_names <- c(yvars_names_short, xlab_names, names_expl_variables)


## ---- Read plot-level species turnover (SAR slopes; generated in PART 5) and calculate site-specific SAR --------------
sar_fit_results_baseplots_1m <- readRDS(paste0(dir_results_patchiness, "/", filename_base_all, "_SAR_baseplots_1m.rds"))

turnover <- do.call("rbind", sar_fit_results_baseplots_1m) 

## ---- Create and save paired correlation plot for ALL base plots combined --------------
vars_baseplot_for_cor <- df %>%
  inner_join(turnover, by=c("plotID")) %>%
  dplyr::select(one_of(c("plotID", variables))) %>%
  filter(complete.cases(.)) # from 157 to 156 plots -->  #SJER_007 removed because all pixels masked for the calculation of sdiv

vars_baseplot_for_cor <- vars_baseplot_for_cor[,variables]
labs <- names(vars_baseplot_for_cor)
names(vars_baseplot_for_cor) <- variables_names

Fig_name_corr_all <- paste0(dir_results_explorative, "/", filename_base_all, "_corrgram", ".png")
png(Fig_name_corr_all, width=1900, height=1900, res=200)
corrplot.mixed(cor(vars_baseplot_for_cor),
               lower = "number",
               upper = "ellipse",
               diag = "u" ,
               tl.pos = "lt",
               tl.col = "black", tl.srt = 50)
dev.off()


## ---- Visualize boxplots with Dunn test  --------------
vars_baseplot_for_Dunn <- df %>%
  inner_join(turnover, by=c("plotID", "siteID")) %>%
  filter(complete.cases(.))

# Perform Kruskal-Wallis test for each variable
# does not assume normality, nor homoscedasticity
# http://userweb.eng.gla.ac.uk/umer.ijaz/bioinformatics/ecological.html
ks.test <- kruskal.test(vars_baseplot_for_Dunn[,"SR"], g=vars_baseplot_for_Dunn[,"siteID"])
ks.test

# https://rcompanion.org/rcompanion/d_06.html
# If the Kruskal-Wallis test is significant, a post-hoc analysis can be performed 
# to determine which levels of the independent variable differ from each other level.
# options:
#   - Dunn test for multiple comparisons
#   - Nemenyi test for multiple comparisons
#   - Pairwise Mann-Whitney U-tests = Wilcoxon rank test
#  Zar (2010) states that the Dunn test is appropriate for groups with unequal numbers of observations.

if (!require('FSA')) install.packages('FSA'); library('FSA') # dunnTest
# install.packages("remotes")
# remotes::install_github("GegznaV/BioStat")
if (!require('biostat')) install.packages('biostat'); library('biostat') # make_cld


variables <- c("SR", "cv_masked_log", "lai", "log_z", "SR_I_NI_perc")
variables_names <- c("(a) Species richness", "(b) Spectral diversity", "(c) Leaf area index", "(d) Species turnover", "(e) Non-native proportion")
fig_ncol = 5

vars_baseplot_for_Dunn$siteID <- as.factor(vars_baseplot_for_Dunn$siteID)
boxplots_Dunn <- Dunn_test(vars_baseplot_for_Dunn, variables, variables_names, group_ID = "siteID", 
                           adj_method="bh", dunn.alpha="0.05", fig_ncol=fig_ncol)

Fig_name_Dunn <- paste0(dir_results_explorative, "/", filename_base_all, "_boxplots_Dunn_", fig_ncol, "cols", ".png")
ggsave(Fig_name_Dunn, boxplots_Dunn, width = 8, height = 4, dpi = 300)
Fig_name_Dunn_pdf <- paste0(dir_results_explorative, "/", filename_base_all, "_boxplots_Dunn_", fig_ncol, "cols", ".pdf")
ggsave(Fig_name_Dunn_pdf, boxplots_Dunn, width = 8.2, height = 4, dpi = 300)

## ----  Calculate min and max values at each site --------------
vars_baseplot_for_Dunn %>% 
  dplyr::select(siteID, lai, NDVI, stat_herbaceousNPP, soil_perc) %>%
  group_by(siteID) %>% 
  dplyr::summarise(across(everything(), list(mean)))


#################################################################################
## PART 8: Visualize RGB of baseplots
#################################################################################
# This section can be run independently of PARTS 1-5 (because intermediate results are stored), 

# The below code can be applied to any site

# USER-DEFINED-INPUTS
# Define parameters for the AOP remote sensing data
site_code <- "OAES"         # Four-digit NEON site code, character string type.
aop_data_year <- "2017"     # Four-digit year in character string "YYYY" for AOP imagery
dp_rgb <- "DP3.30010.001"   # High-resolution multispectral camera RGB
buffer_val <- 0 #[m]        # integer buffer size around coordinates to download image

## ---- Open plot coordinates --------------
filename <- paste0(dir_data_out,"/",site_code,"_",aop_data_year,"/",paste0(site_code, "_", aop_data_year,"_veg_plots_spcomp.rda"))
veg_plots_spcomp_site_year <- readRDS(filename)
folder_coords <- paste0(dir_data_out,"/",site_code,"_",aop_data_year,"/plot_locations/")
veg_plots_base_coords <- readRDS(paste0(folder_coords,paste0(site_code,"_veg_plots_base.rds")))

plots_herb_descr <- veg_plots_base_coords$veg_utm  %>%
  dplyr::filter(nlcdClass == "grasslandHerbaceous")
plots_herb_poly <- veg_plots_base_coords$veg_spdf_poly[veg_plots_base_coords$veg_spdf_poly$id %in% plots_herb_descr$plotID, ]

# Define plot coordinates
poly_coordinates <- data.frame()
for (i in 1:dim(plots_herb_poly)[1]){
  eastings <- c(seq(plots_herb_poly[i,]@bbox[1,"min"],plots_herb_poly[i,]@bbox[1,"max"],by=500),plots_herb_poly[i,]@bbox[1,"max"]) # by adding the max we ensure that the full polygon is covered for download
  northings <- c(seq(plots_herb_poly[i,]@bbox[2,"min"],plots_herb_poly[i,]@bbox[2,"max"],by=500),plots_herb_poly[i,]@bbox[2,"max"])
  poly_coordinates_plot <- tidyr::expand_grid(eastings, northings)
  poly_coordinates <- rbind(poly_coordinates, poly_coordinates_plot)
}

## ---- Convert list of plot coordinates to a list of 1km x 1km tile coordinates --------------
source("01_Download_preprocess_aop_imagery.R") # contains the function "list_tiles_covering_poly_tibble"
tile_coordinates <- list_tiles_covering_poly_tibble(poly_coordinates)

## ---- Download RGB tiles and save in raw data folder --------------
print("Downloading RGB")
neonUtilities::byTileAOP(
  dpID = dp_rgb
  ,site = site_code
  ,year = aop_data_year
  ,savepath = dir_data_out_temp
  ,easting = tile_coordinates$eastings
  ,northing = tile_coordinates$northings
  ,check.size = F
  ,buffer = buffer_val)

# move rgb files
move_downloaded_files(dir_out = dir_data_out_temp, dp_id = dp_rgb
                      ,dp_name = "rgb", file_pattern = "*image.tif$"
                      ,delete_orig = T)

## ---- Create site-specific output folder --------------
dir_data_out_siteyear <- file.path(dir_data_out,paste(site_code,aop_data_year,sep="_"))
check_create_dir(dir_data_out_siteyear)
dir_data_out_proc <- file.path(dir_data_out_siteyear,"rgb_plots")
check_create_dir(dir_data_out_proc)

## ---- Cut tiles to the extent of a plot and save the results as an .rds file --------------
veg_plots_base_rgb <- get_plot_rgb(plots_herb_poly, site_code, aop_data_year, dir_data_out_temp)
saveRDS(veg_plots_base_rgb, file=file.path(dir_data_out_proc,
                                           paste0("veg_plots_base_rgb_", site_code, "_", aop_data_year, ".rds")))

## ---- Visualize and save plot-level RGB images --------------
for (i in 1:length(veg_plots_base_rgb)){
  tiff(file.path(dir_data_out_proc, 
                 paste0("rgb_", aop_data_year, "_", names(veg_plots_base_rgb)[i], ".tif")))
  plotRGB(veg_plots_base_rgb[[i]])
  plot(plots_herb_poly, add=T, lwd=2)
  dev.off()
}


#################################################################################
## PART 9: Visualize 1st spectral PC and hyperspectral signatures of baseplots
#################################################################################
# This section can be run independently of PARTS 2-5 (because intermediate results are stored), 
# but make sure to read and prepare base data, by running PART 1

# The below code can be applied to any site

if (!require('rasterVis')) install.packages('rasterVis'); library('rasterVis')
if (!require('RColorBrewer')) install.packages('RColorBrewer'); library('RColorBrewer')

# USER-DEFINED-INPUTS
# Define parameters for the AOP remote sensing data
site_code <- "OAES"         # Four-digit NEON site code, character string type.
aop_data_year <- "2017"     # Four-digit year in character string "YYYY" for AOP imagery
ndvi_tresh <- 0.2           # NDVI treshold used for masking pixels (0.3 is the treshold applied in Lopatin et al. 2017)  
scaling_nb <-1              # PCA scaling type (i.e. type I or II, more info in Appendix S4 of Laliberté et al. 2020)


## ---- Open plot coordinates --------------
filename <- paste0(dir_data_out,"/",site_code,"_",aop_data_year,"/",paste0(site_code, "_", aop_data_year,"_veg_plots_spcomp.rda"))
veg_plots_spcomp_site_year <- readRDS(filename)
folder_coords <- paste0(dir_data_out,"/",site_code,"_",aop_data_year,"/plot_locations/")
veg_plots_base_coords <- readRDS(paste0(folder_coords,paste0(site_code,"_veg_plots_base.rds")))

plots_herb_descr <- veg_plots_base_coords$veg_utm  %>%
  dplyr::filter(nlcdClass == "grasslandHerbaceous")
plots_herb_poly <- veg_plots_base_coords$veg_spdf_poly[veg_plots_base_coords$veg_spdf_poly$id %in% plots_herb_descr$plotID, ]

## ---- Open baseplot spectral signatures --------------
dir_data_out_proc <- file.path(dir_data_out,paste(site_code,aop_data_year,sep="_"),"stacked_aop_data_processed")
aop_veg_plots_base_filename = file.path(dir_data_out_proc,paste0("aop_data_hs_veg_plots_base_",site_code, "_", aop_data_year, "_", 
                                                                 ndvi_tresh,"ndviTresh.rds"))
aop_veg_plots_base <- readRDS(aop_veg_plots_base_filename)

# Only retain base plots that were used in the analyses
df_site <- df %>%
  dplyr::filter(siteID == site_code)
aop_veg_plots_base_sel <- aop_veg_plots_base[df_site[,"plotID"]]

## ---- Create site-specific output folder --------------
dir_data_out_siteyear <- file.path(dir_data_out,paste(site_code,aop_data_year,sep="_"))
check_create_dir(dir_data_out_siteyear)
dir_data_out_aop_plots <- file.path(dir_data_out_siteyear,"aop_plots")
check_create_dir(dir_data_out_aop_plots)

## ---- Visualize and save 1st PC of baseplots --------------
palette <- colorRampPalette(rev(brewer.pal(n=9,name="Spectral")))

for (i in 1:length(aop_veg_plots_base_sel)){
  
  # Variation in PC1 before masking based on NDVI and canopy height
  tiff(file.path(dir_data_out_aop_plots,
                 paste0("pc1_", aop_data_year, "_", names(aop_veg_plots_base_sel)[i], ".tif")))
  r <- rasterVis::levelplot(aop_veg_plots_base_sel[[i]]$aop_data_hs_pca$cube_pc$PC1, 
                            colorkey=F,
                            col.regions = palette,
                            xlab=NULL, ylab=NULL,
                            scales=list(draw=FALSE),
                            margin = F) + 
    latticeExtra::layer(sp.polygons(plots_herb_poly, col="black"))
  print(r)
  dev.off()
  
  # Variation in PC1 after masking based on NDVI and canopy height
  if (is.na(aop_veg_plots_base_sel[[i]]$aop_data_hs_masked_pca)[[1]] == TRUE){
    print(paste0("masked image of ", names(aop_veg_plots_base_sel)[i], " contains NA only"))
  } else {
    tiff(file.path(dir_data_out_aop_plots,
                   paste0("pc1_", aop_data_year, "_", names(aop_veg_plots_base_sel)[i], 
                          "_masked_", ndvi_tresh, "ndviTresh.tif")))
    r2 <- rasterVis::levelplot(aop_veg_plots_base_sel[[i]]$aop_data_hs_masked_pca$cube_pc$PC1, 
                               # at=seq(ndvi_min, ndvi_max, length.out=500),
                               # colorkey = list(space = "bottom"),
                               colorkey=F,
                               col.regions = palette,
                               xlab=NULL, ylab=NULL,
                               scales=list(draw=FALSE),
                               margin = F) + 
      latticeExtra::layer(sp.polygons(plots_herb_poly, col="black"))
    print(r2)
    dev.off()
  }
  
  graphics.off()
}

## ---- Visualize and save spectral signatures of baseplots --------------
for (i in 1:length(aop_veg_plots_base_sel)){
  print(paste0("Plotting spectral signatures for ", names(aop_veg_plots_base_sel)[i]))
  
  # Signatures before masking based on NDVI and canopy height
  plots_hs_df <- as.data.frame(aop_veg_plots_base_sel[[i]]$cleaned_spectra$aop_data_hs_clean)
  wl <- aop_veg_plots_base[[i]]$cleaned_spectra$wavelengths_clean
  
  p1 <- plot_sign(plots_hs_df, wl)
  ggsave(file.path(dir_data_out_aop_plots, paste0("spectra_", aop_data_year, "_", names(aop_veg_plots_base_sel)[i],
                                                  ".png")),
         p1, width=4, height=3)
  
  # Signatures after masking based on NDVI and canopy height
  plots_hs_df <- as.data.frame(aop_veg_plots_base_sel[[i]]$cleaned_masked_spectra$aop_data_hs_clean) %>%
    tidyr::drop_na()
  wl <- aop_veg_plots_base[[i]]$cleaned_masked_spectra$wavelengths_clean
  
  p2 <- plot_sign(plots_hs_df, wl)
  ggsave(file.path(dir_data_out_aop_plots, paste0("spectra_", aop_data_year, "_", names(aop_veg_plots_base_sel)[i],
                                                  "_masked_", ndvi_tresh, "ndviTresh.png")),
         p2, width=4, height=3)
  
}

#################################################################################
## PART 10: Create illustrative hyperspectral signatures  
#################################################################################
# This section can be run independently of PARTS 2-5 (because intermediate results are stored), 
# but make sure to read and prepare base data, by running PART 1

# The spectral signatures chosen are purely illustrative, and not meant to represent the truth

# USER-DEFINED-INPUTS
# Define parameters for the AOP remote sensing data
site_code <- "SJER"         # Four-digit NEON site code, character string type.
aop_data_year <- "2017"     # Four-digit year in character string "YYYY" for AOP imagery
ndvi_tresh <- 0.2           # NDVI treshold used for masking pixels (0.3 is the treshold applied in Lopatin et al. 2017)  
scaling_nb <-1              # PCA scaling type (i.e. type I or II, more info in Appendix S4 of Laliberté et al. 2020)

## ---- Open plot coordinates --------------
filename <- paste0(dir_data_out,"/",site_code,"_",aop_data_year,"/",paste0(site_code, "_", aop_data_year,"_veg_plots_spcomp.rda"))
veg_plots_spcomp_site_year <- readRDS(filename)
folder_coords <- paste0(dir_data_out,"/",site_code,"_",aop_data_year,"/plot_locations/")
veg_plots_base_coords <- readRDS(paste0(folder_coords,paste0(site_code,"_veg_plots_base.rds")))

plots_herb_descr <- veg_plots_base_coords$veg_utm  %>%
  dplyr::filter(nlcdClass == "grasslandHerbaceous")
plots_herb_poly <- veg_plots_base_coords$veg_spdf_poly[veg_plots_base_coords$veg_spdf_poly$id %in% plots_herb_descr$plotID, ]

## ---- Open baseplot spectral signatures --------------
dir_data_out_proc <- file.path(dir_data_out,paste(site_code,aop_data_year,sep="_"),"stacked_aop_data_processed")
aop_veg_plots_base_filename = file.path(dir_data_out_proc,paste0("aop_data_hs_veg_plots_base_",site_code, "_", aop_data_year, "_", 
                                                                 ndvi_tresh,"ndviTresh.rds"))
aop_veg_plots_base <- readRDS(aop_veg_plots_base_filename)

# Only retain base plots that were used in the analyses
df_site <- df %>%
  dplyr::filter(siteID == site_code)
aop_veg_plots_base_sel <- aop_veg_plots_base[df_site[,"plotID"]]

## ---- Create output folder --------------
dir_data_out_conc <- file.path(dir_data_out_results,"conceptual_fig")
check_create_dir(dir_data_out_conc)

## ---- Visualize spectral signatures for plot and pixels of choice --------------
colors_gray <- c("gray5", "gray25", "gray45", "gray70")

# Low density figure
i = 5 # SJER_002
plot_hs_df_dens_low <- as.data.frame(aop_veg_plots_base_sel[[i]]$cleaned_masked_spectra$aop_data_hs_clean)
wl_dens_low <- aop_veg_plots_base[[i]]$cleaned_masked_spectra$wavelengths_clean
hs_df_dens_low <- plot_hs_df_dens_low[c(40, 397, 100, 200),]
colors_dens_low <- colors_gray
p_dens_low <- plot_sign_choice(hs_df_dens_low, wl_dens_low, colors_dens_low)
ggsave(file.path(dir_data_out_conc, paste0("spectra_density_low_", aop_data_year, "_", names(aop_veg_plots_base_sel)[i], "_gray.png")),
       p_dens_low, width=3, height=3)

# Low species turnover figure
i = 14  # SJER_019
plot_hs_df_turn_low <- as.data.frame(aop_veg_plots_base_sel[[i]]$cleaned_masked_spectra$aop_data_hs_clean)
wl_turn_low <- aop_veg_plots_base[[i]]$cleaned_masked_spectra$wavelengths_clean
hs_df_turn_low <- plot_hs_df_turn_low[c(40, 397, 100, 200),]
colors_turn_low <- colors_gray
p_turn_low <- plot_sign_choice(hs_df_turn_low, wl_turn_low, colors_turn_low)
ggsave(file.path(dir_data_out_conc, paste0("spectra_turnover_low_", aop_data_year, "_", names(aop_veg_plots_base_sel)[i], "_gray.png")),
       p_turn_low, width=3, height=3)

# High invasion figure
i = 16 # SJER_005
plot_hs_df_inva_high <- as.data.frame(aop_veg_plots_base_sel[[i]]$cleaned_masked_spectra$aop_data_hs_clean)
wl_inva_high <- aop_veg_plots_base[[i]]$cleaned_masked_spectra$wavelengths_clean
hs_df_inva_high <- plot_hs_df_inva_high[c(316, 360, 322,  325),]
colors_inva_high <- colors_gray
p_inva_low <- plot_sign_choice(hs_df_inva_high, wl_inva_high, colors_inva_high)
ggsave(file.path(dir_data_out_conc, paste0("spectra_invasion_high_", aop_data_year, "_", names(aop_veg_plots_base_sel)[i], "_gray.png")),
       p_inva_low, width=3, height=3)

# Good correlation figure: high density, high turnover and low invasion
i = 13 # SJER_001
plot_hs_df_dens_high <- as.data.frame(aop_veg_plots_base_sel[[i]]$cleaned_masked_spectra$aop_data_hs_clean)
wl_dens_high <- aop_veg_plots_base[[i]]$cleaned_masked_spectra$wavelengths_clean
# hs_df_dens_high <- plot_hs_df_dens_high[c(40, 397, 100, 200),]
hs_df_dens_high <- plot_hs_df_dens_high[c(43, 56, 182, 240),]
colors_dens_high <- colors_gray
p_dens_high <- plot_sign_choice(hs_df_dens_high, wl_dens_high, colors_dens_high)
ggsave(file.path(dir_data_out_conc, paste0("spectra_density_high_", aop_data_year, "_", names(aop_veg_plots_base_sel)[i], "_gray.png")),
       p_dens_high, width=3, height=3)

#################################################################################
