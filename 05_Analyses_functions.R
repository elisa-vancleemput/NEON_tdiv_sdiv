###########################################################
###                                                     ###
###                  NEON baseplot analyses             ###
###                                                     ###
###########################################################

# This script contains functions to read and concatenate all information at the NEON baseplot level (20 m x 20 m)
# This information is used to 
# - calculate taxonomic diversity,
#   which can then be linked to spectral diversity by means of
# - linear regression models

# by Elisa Van Cleemput, 2023
###########################################################
# Open and concatenate baseplot datasets
# Used in PARTS 1, 2, 5

open_baseplot_data <- function(dir_data_out, site_code, aop_data_year, ndvi_tresh) {
  

  ## ---- Open species composition --------------
  dir_data_out_siteyear <- file.path(dir_data_out,paste(site_code,aop_data_year,sep="_"))
  filename <- paste0(dir_data_out_siteyear, "/",paste0(site_code, "_", aop_data_year,"_veg_plots_spcomp.rda"))
  veg_plots_spcomp_site_year <- readRDS(filename)

  
  # ---- Open baseplot information (sdiv)  -------------
  dir_data_out_proc <- file.path(dir_data_out,paste(site_code,aop_data_year,sep="_"),"stacked_aop_data_processed")
  baseplot_filename = file.path(dir_data_out_proc,paste0("veg_plots_base_all_", 
                                                         site_code, "_", aop_data_year, "_",
                                                         ndvi_tresh,"ndviTresh.rds"))
  veg_plots_base <- readRDS(baseplot_filename)    
  
  
  out <- list("veg_plots_spcomp_site_year" = veg_plots_spcomp_site_year,
              "veg_plots_base" = veg_plots_base)
  return(out)
}

preparing_baseplot_data <- function(veg_plots_base_data, time_start, time_end){
  
  ## ---- Linking ecological information to sdiv data --------------
  # Keep only those baseplots for which we calculated sdiv
  # Data on species presence (recorded for entire 400 m² plots)
  veg_plots_base_species_pres <- veg_plots_base_data$veg_plots_spcomp_site_year$table_observation_species_baseplots %>%
    inner_join(veg_plots_base_data$veg_plots_spcomp_site_year$table_sample_info %>% 
                 dplyr::select(siteID, sampleID, namedLocation, nlcdClass, elevation, decimalLatitude, decimalLongitude), c("sampleID")) %>%
    mutate(plot_subplotID = paste0(plotID, "_",subplotID)) %>%
    unique() %>%
    dplyr::filter(namedLocation %in% veg_plots_base_data$veg_plots_base$namedLocation)
  
  # Data on species abundance (recorded in 1 m² subplots only)
  veg_plots_base_species_abund <- veg_plots_base_data$veg_plots_spcomp_site_year$table_observation_species_subplots %>%
    inner_join(veg_plots_base_data$veg_plots_spcomp_site_year$table_sample_info %>% 
                 dplyr::select(siteID, sampleID, namedLocation, nlcdClass, elevation, decimalLatitude, decimalLongitude), c("sampleID")) %>%
    mutate(plot_subplotID = paste0(plotID, "_",subplotID)) %>%
    unique() %>%
    dplyr::filter(namedLocation %in% veg_plots_base_data$veg_plots_base$namedLocation)
  
  print(paste0(unique(veg_plots_base_species_pres$siteID), ": Total no. of species presence records is ",dim(veg_plots_base_species_pres)[1]))
  print(paste0(unique(veg_plots_base_species_abund$siteID), ": Total no. of species abundance records is ",dim(veg_plots_base_species_abund)[1]))
  
  # Data on abiotic elements and non photosynthetic active vegetation
  veg_plots_base_dead <- veg_plots_base_data$veg_plots_spcomp_site_year$table_observation %>%
    dplyr::select(plotID, subplotID, sampleID, endDate, month, otherVariables, percentCover) %>%
    # dplyr::filter(!is.na(otherVariables)) %>%  # dplyr::select records that are not live species
    dplyr::filter(otherVariables %in% c("soil", "standingdead", "rock", "litter", "wood", "water")) %>%  # dplyr::select records that are not live species
    distinct() %>%
    inner_join(veg_plots_base_data$veg_plots_spcomp_site_year$table_sample_info %>% 
                 dplyr::select(siteID, sampleID, namedLocation, nlcdClass, elevation, decimalLatitude, decimalLongitude), c("sampleID")) %>%
    mutate(plot_subplotID = paste0(plotID, "_",subplotID)) %>%
    unique() %>%
    dplyr::filter(namedLocation %in% veg_plots_base_data$veg_plots_base$namedLocation)
  
  
  ## ---- Filter on date --------------
  if(class(time_start)[1] == "POSIXct"){
    veg_plots_base_species_pres <- veg_plots_base_species_pres %>%
      dplyr::filter(endDate <= time_end) %>%
      dplyr::filter(endDate >= time_start)
    print(paste0("No. of species presence records after filtering based on time is ",dim(veg_plots_base_species_pres)[1]))
    veg_plots_base_species_abund <- veg_plots_base_species_abund %>%
      dplyr::filter(endDate <= time_end) %>%
      dplyr::filter(endDate >= time_start)
    print(paste0("No. of species abundance records after filtering based on time is ",dim(veg_plots_base_species_abund)[1]))
    veg_plots_base_dead <- veg_plots_base_dead %>%
      dplyr::filter(endDate <= time_end) %>%
      dplyr::filter(endDate >= time_start)
  } else {
    veg_plots_base_species_pres <- veg_plots_base_species_pres
    veg_plots_base_species_abund <- veg_plots_base_species_abund
    veg_plots_base_dead <- veg_plots_base_dead
  }
  
  ## ---- Add information on native/alien status --------------
  veg_plots_base_species_pres <- veg_plots_base_species_pres %>% 
    left_join(veg_plots_base_data$veg_plots_spcomp_site_year$table_taxon %>%
                dplyr::select(taxonID, nativeStatusCode, family) %>% unique(),
              by =c("taxonID"))
  
  ## ---- Linking spectral diversity calculations to ecological information --------------
  veg_plots_base_filtered <- veg_plots_base_data$veg_plots_base %>%
    dplyr::filter(namedLocation %in% unique(veg_plots_base_species_abund$namedLocation)) 
  print(paste0("No. of baseplots with sdiv data is ",dim(veg_plots_base_filtered)[1]))
  
  ## ----Dealing with duplicates --------------
  # pivot to wide format
  veg_plots_base_species_abund_wide <- veg_plots_base_species_abund %>%
    dplyr::group_by(siteID, plotID, subplotID, plot_subplotID, sampleID, namedLocation,
                    month, endDate,
                    elevation, decimalLatitude, decimalLongitude,
                    taxonID) %>%
    dplyr::summarise(percentCover = mean(percentCover), .groups = "keep") %>% # here, we average species with multiple records
    tidyr::pivot_wider(id_cols = c(siteID, plotID, subplotID, plot_subplotID, sampleID, namedLocation,
                                   month, endDate,
                                   elevation, decimalLatitude, decimalLongitude), 
                       names_from = taxonID,
                       values_from = percentCover,
                       values_fill = list(percentCover = 0))
  
  
  print(paste0("No. of subplots (incl. duplicates) surveyed for species abundance is ",dim(veg_plots_base_species_abund_wide)[1]))
  
  # Find out which plots occur more than once
  n_occur <- data.frame(table(veg_plots_base_species_abund_wide$plot_subplotID))
  n_occur[n_occur$Freq > 1,]
  duplicate_observation_species <- veg_plots_base_species_abund_wide[veg_plots_base_species_abund_wide$plot_subplotID %in% 
                                                                       n_occur$Var1[n_occur$Freq > 1],1:6]
  
  print(paste0("No. of duplicate (sampled multiple times) baseplots is ", length(unique(duplicate_observation_species$plot_subplotID)))) 
  
  # combine survey results whereby each species is attributed it's maximal abundance
  veg_plots_base_species_abund_wide_duplicMax <- veg_plots_base_species_abund_wide %>%
    group_by(siteID, plotID, subplotID, plot_subplotID, namedLocation,
             elevation, decimalLatitude, decimalLongitude) %>%
    summarize_all(max) %>%
    # summarize(across(everything(), max))
    as.data.frame() %>%
    dplyr::select(-endDate, -sampleID) 
  
  print(paste0("No. of baseplots without duplicates is ",dim(veg_plots_base_species_abund_wide_duplicMax)[1]))
  
  # Remark: In the previous lines month and enDate were also "maximized", but this is not meaningful.
  # Here, we will replace the month value of the duplicate baseplots
  veg_plots_base_species_abund_wide_duplicMax[veg_plots_base_species_abund_wide_duplicMax$plot_subplotID %in% 
                                                duplicate_observation_species$plot_subplotID,"month"] <- 
    rep("Multiple",length(unique(duplicate_observation_species$plot_subplotID)))
  
  ## ---- Combine everything in 1 output list --------------
  output <- list("veg_plots_base_species_pres" = veg_plots_base_species_pres,
                 "veg_plots_base_species_abund_wide" = veg_plots_base_species_abund_wide_duplicMax,
                 "veg_plots_base_sdiv" = veg_plots_base_filtered,
                 "veg_plots_base_dead" = veg_plots_base_dead)
  
  return(output)
}

###########################################################
# Calculate abundance-based taxonomic diversity metrics
# Used in PART 1

if (!require('vegan')) install.packages('vegan'); library('vegan')

Calculate_plot_taxon_div <- function(table_observation_species_wide){
  table_observation_data <- table_observation_species_wide %>% 
    dplyr::select(-siteID, -subplotID, -plot_subplotID, -namedLocation, -elevation, -decimalLatitude, -decimalLongitude, -month) %>%
    dplyr::group_by(plotID) %>%
    summarize_all(sum) %>%
    tibble::column_to_rownames(var = "plotID")
  
  Shannon <- vegan::diversity(table_observation_data, index="shannon")
  Simpson <- vegan::diversity(table_observation_data, index="simpson")
  Simpson_inv <- vegan::diversity(table_observation_data, index="invsimpson")
  n_subplots <- ddply(veg_plots_base_spComp, .(plotID), nrow)
  names(n_subplots) <- c("plotID", "n_subplots")
  taxon_diversity <- cbind(Shannon, Simpson, Simpson_inv) %>%
    as.data.frame() %>%
    tibble::rownames_to_column(var = "plotID") %>%
    inner_join(n_subplots, by="plotID")
  
  return(taxon_diversity)
}

###########################################################
# Create linear regression formula with categorical variable interacting with continuous explanatory variable
# Used in PART 4

create_lm_formula <- function(yvar, xvar_cont, xvars_categ){
  f <- as.formula(paste0(yvar, "~", xvar_cont, "+", 
                         paste(xvars_categ, collapse = "+"),"+",
                         paste(paste(xvar_cont, xvars_categ, sep="*"), collapse = "+")))
  return(f)
}

###########################################################
# Function to create boxplots with Dunn test
# code adapted from Van Cleemput et al. 2021. in Ecological Indicators
# "Spectrally defined plant functional types adequately capture multidimensional trait variation in herbaceous communities"
# https://doi.org/10.1016/j.ecolind.2020.106970

# Used in PART 7
boxplots_Dunn <- Dunn_test(vars_baseplot_for_Dunn, variables, variables_names, group_ID = "siteID", 
                           adj_method="bh", dunn.alpha="0.05", fig_ncol=fig_ncol)

Dunn_test <- function(dataset_full, variables, variables_names, group_ID, group_order_nb, adj_method, dunn.alpha, fig_ncol){
  
  dataset <- dataset_full[,c(variables)]
  names(dataset) <- variables_names
  Profile <- dataset_full[,c(group_ID, variables)]
  names(Profile) <- c(group_ID, variables_names)
  Profile.long <- melt(Profile,id.vars=c(group_ID))
  
  dunn.table <- data.frame()
  dunn.table_letters <- data.frame()
  dunn.table.adj <- data.frame()
  dunn.table.adj_letters <- data.frame()
  
  ## ---- 1) Perform Dunn test for each variables  --------------
  for (i in 1:length(variables_names)){
    print(paste0("Dunn test on ", variables_names[i]))
    d.test <- dunnTest(dataset[,variables_names[i]], g=dataset_full[,group_ID],method=adj_method)
    
    sign.pairs.adj <- d.test$res[which(d.test$res$P.adj <= dunn.alpha),]
    nb.sign.pairs.adj <- length(which(d.test$res$P.adj <= dunn.alpha))
    if (nb.sign.pairs.adj != 0){
      sign.pairs.adj <- cbind(sign.pairs.adj, variable = names(dataset)[i])
      dunn.table.adj <- rbind(dunn.table.adj, sign.pairs.adj)
      
      letters <- make_cld(P.adj ~ Comparison , data=d.test$res, alpha <= dunn.alpha)
      sign.pairs.adj_letters <- cbind(letters, variable = names(dataset)[i])
      dunn.table.adj_letters <- rbind(dunn.table.adj_letters, sign.pairs.adj_letters)
    } else {
      sign.pairs.adj <- 0
      dunn.table.adj <- dunn.table.adj
      dunn.table.adj_letters <- dunn.table.adj_letters
      
    }
    
    sign.pairs <- d.test$res[which(d.test$res$P.unadj < dunn.alpha),]
    nb.sign.pairs <- length(which(d.test$res$P.unadj < dunn.alpha))
    if (nb.sign.pairs != 0){
      sign.pairs <- cbind(sign.pairs, variable = names(dataset)[i])
      dunn.table <- rbind(dunn.table, sign.pairs)
      
      letters <- make_cld(P.unadj ~ Comparison , data=d.test$res, alpha <= dunn.alpha)
      sign.pairs.letters <- cbind(letters, variable = names(dataset)[i])
      dunn.table_letters <- rbind(dunn.table_letters, sign.pairs.letters)
    } else {
      sign.pairs <- 0
      dunn.table <- dunn.table
      dunn.table_letters <- dunn.table_letters
    }
  }
  
  ## ---- 2) Create figure --------------
  # letters from make_cld are the labels we need (cld = compact letter display)
  
  dunn.table.adj_letters$variable <- as.factor(dunn.table.adj_letters$variable)

  # Determine positions for labels
  abs_max <- apply(dataset, 2, max, na.rm=T)
  abs_max_rep <- rep(abs_max, nlevels(Profile$siteID))
  maxs <- Profile.long %>%
    dplyr::group_by(siteID, variable) %>%
    dplyr::summarise(max=max(value)) %>%
    mutate(value= abs_max_rep + 0.1 * abs_max_rep) %>%
    dplyr::rename ("group" = siteID) %>%
    arrange(factor(variable, levels=variables_names)) %>%
    left_join(dunn.table.adj_letters, by=c("group","variable"))

  # Make figure
  p<-ggplot(Profile.long)+ #,fill=groups
    geom_boxplot(aes(forcats::fct_rev(siteID),value, fill=siteID))+
    ggsci::scale_fill_igv() +
    geom_text(data=maxs, aes(x = group, y = value, label=cld), hjust="inward", color = "gray32")+
    coord_flip()+
    facet_wrap(~variable, scales="free_x", ncol = fig_ncol) +
    theme_bw() +
    theme(legend.position="none",
          axis.title.x = element_blank(), axis.title.y = element_blank(),)
  print(p)
  
  return(p)
}

###########################################################
# Functions to retrieve plot-level RGB data

# Used in PART 8

source("01_Download_preprocess_aop_imagery.R") # function: list_tiles_covering_poly

# mosaic function copied from "06_Locate_plot_samples.R"
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

get_plot_rgb <- function(plot_set, site_code, aop_data_year, dir_data_out_temp) {
  
  # Loop through the polygons dataset
  rgb_plot_list <- list()
  
  for (i in 1:(dim(plot_set)[1])){
    
    print(paste0("Processing plot No. ", i, " of ", dim(plot_set)[1]))
    
    ## ---- Make a list of unique easting, northing coordinates (500 m steps) representing a specific plot --------------
    eastings <- c(seq(plot_set[i,]@bbox[1,"min"],plot_set[i,]@bbox[1,"max"],by=500),plot_set[i,]@bbox[1,"max"]) # by adding the max we ensure that the full polygon is covered for download
    northings <- c(seq(plot_set[i,]@bbox[2,"min"],plot_set[i,]@bbox[2,"max"],by=500),plot_set[i,]@bbox[2,"max"])
    # poly_coordinates <- rbind(poly_coordinates,tidyr::expand_grid(eastings, northings))
    poly_coordinates <- tidyr::expand_grid(eastings, northings)
    
    ## ----  generate a list of the 1km x 1km tile coordinates covering this specific plot --------------
    tile_coordinates <- list_tiles_covering_poly(poly_coordinates)
    
    ## ---- Specify which tiles (.rds stacked AOP data file) to open --------------
    rgb_filename_nb <- grep(tile_coordinates, list.files(file.path(dir_data_out_temp, "rgb")))
    rgb_filename <- file.path(dir_data_out_temp, "rgb", 
                              list.files(file.path(dir_data_out_temp, "rgb"))[rgb_filename_nb])
    
    ## ---- Load the raster tiles covering the polygon(s) --------------
    filename = rgb_filename
    rgb_data_list <- vector(mode = "list", length = length(filename))
    names(rgb_data_list) <- filename
    system.time(
      for (f in 1:length(filename)){
        # print(paste0("Reading tile No. ", f, " of ", length(filename)))
        rgb_data_list[[f]] <- brick(filename[f])
      })
    
    ## ---- Extract plot level data --------------
    ### Crop raster(s) to the extent of the polygons
    print("Cropping tiles to the extent of the polygon")
    rgb_cropped_list <- lapply(rgb_data_list, function(x) raster::crop(x, extent(plot_set[i,])))
    
    ### mosaic tiles if a plot covers > 1 tile
    rgb_data_mosaicked <- mosaic_tiles(rgb_cropped_list)
    
    ## ---- Add results to list of plots --------------
    rgb_plot_list[[i]] <- rgb_data_mosaicked
  }
  
  names(rgb_plot_list) <- plot_set$id
  
  return(rgb_plot_list)
}


###########################################################
# Functions to plot hyperspectral signatures of one plot

# Used in PARTS 9 and 10

# Plot all signatures of one plot
plot_sign <- function(plots_hs_df, wl){
  names(plots_hs_df) <- wl
  plots_hs_df$pixelNb <- seq(1:nrow(plots_hs_df))
  
  plots_hs_df_long <- reshape2::melt(plots_hs_df, id.vars="pixelNb", variable.name="wl", value.name="refl")
  plots_hs_df_long$wl <- as.numeric(as.character(plots_hs_df_long$wl))
  
  window1 <- c(1340, 1445)
  window2 <- c(1790, 1955)
  # bad bands include 2 band left and right of these windows
  window1_left <- wl[which.min(abs(wl - window1[1]))-2]
  window1_right <- wl[which.min(abs(wl - window1[2]))+2]
  window2_left <- wl[which.min(abs(wl - window2[1]))-2]
  window2_right <- wl[which.min(abs(wl - window2[2]))+2]
  
  g <- ggplot()+
    geom_line(data=plots_hs_df_long,aes(x=wl,y=refl, col=pixelNb,group=pixelNb), lwd=0.5) +
    geom_rect(aes(xmin=window1_left, xmax=window1_right, ymin=0.0001, ymax=round(max(plots_hs_df_long[,"refl"]),1)), fill="white") +
    geom_rect(aes(xmin=window2_left, xmax=window2_right, ymin=0.0001, ymax=round(max(plots_hs_df_long[,"refl"]),1)), fill="white") +
    scale_x_continuous(breaks=c(400, 2400)) +
    labs(x = "Wavelength (nm)", y="Reflectance [0-1]") +
    theme_classic() +
    theme(legend.position="none") +
    theme(axis.title = element_text(size=14),
          axis.text = element_text(size=14))
  
  # Add number of pixels to plot
  if (!require('grid')) install.packages('grid'); library('grid')
  grob <- grobTree(textGrob(paste0("n pixels = ", nrow(plots_hs_df)), 
                            # x=max(plots_hs_df_long$refl),  y=2000, hjust=0,
                            x=0.55,  y=0.96, hjust=0,
                            gp=gpar(col="black", fontsize=16, fontface="italic")))
  g2 <- g + annotation_custom(grob)
  
  return(g2)
}

# Plot signatures with custom color pallette 
plot_sign_choice <- function(plots_hs_df, wl, colors){
  names(plots_hs_df) <- wl
  plots_hs_df$pixelNb <- seq(1:nrow(plots_hs_df))
  
  plots_hs_df_long <- reshape2::melt(plots_hs_df, id.vars="pixelNb", variable.name="wl", value.name="refl")
  plots_hs_df_long$wl <- as.numeric(as.character(plots_hs_df_long$wl))
  
  window1 <- c(1340, 1445)
  window2 <- c(1790, 1955)
  # bad bands include 2 band left and right of these windows
  window1_left <- wl[which.min(abs(wl - window1[1]))-2]
  window1_right <- wl[which.min(abs(wl - window1[2]))+2]
  window2_left <- wl[which.min(abs(wl - window2[1]))-2]
  window2_right <- wl[which.min(abs(wl - window2[2]))+2]
  
  
  g <- ggplot(plots_hs_df_long,aes(x=wl,y=refl))+
    geom_line(aes(col=factor(pixelNb), group=pixelNb), lwd=2) +
    scale_color_manual(values = colors) +
    geom_rect(aes(xmin=window1_left, xmax=window1_right, ymin=0.0001, ymax=max(plots_hs_df_long[,"refl"])), fill="white") +
    geom_rect(aes(xmin=window2_left, xmax=window2_right, ymin=0.0001, ymax=max(plots_hs_df_long[,"refl"])), fill="white") +
    scale_x_continuous(breaks=c(400, 2400)) +
    labs(x = "Wavelength (nm)") +
    theme_classic() +
    theme(legend.position = "none") +
    theme(axis.title = element_blank(),
          axis.text = element_blank(),
          axis.ticks =element_blank())
 
  return(g)
}


###########################################################

