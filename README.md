# NEON_tdiv_sdiv
Code supporting 

Van Cleemput, E., Adler, P., Suding, K. in Global Ecology and Geography
Making remote sense of biodiversity: What grassland characteristics make spectral diversity a good proxy for taxonomic diversity?

**GEE_Download_RAP_NPP.txt**
Contains code to download herbaceous NPP from the Rangelands Analysis Platform in Google Earth Engine
The code is also directly accessible here:
https://code.earthengine.google.com/?accept_repo=users/elisavancleemput/NEON_tdiv_sdiv

**OO_Main_script_downloads.R**
Contains code to download, preprocess and combine imagery and vegetation plot information from NEON.
The code is supported by functions that are stored in 
 - 01_Download_preprocess_aop_imagery.R
 - 02_Download_veg_plots_data.R
 - 03_Extract_plot_level_raster_values.R
 - 04_Calculate_spectral_diversity.R

**00_Main_script_analyses.R**
Contains code to run analyses on the extracted data
The code is supported by functions that are stored in 
 - 01_Download_preprocess_aop_imagery.R: only for the functions "list_tiles_covering_poly_tibble" and "list_tiles_covering_poly"
 - 05_Analyses_functions.R
