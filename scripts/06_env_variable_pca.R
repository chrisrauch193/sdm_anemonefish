# scripts/06_env_variable_pca.R
# This script performs PCA on the selected environmental variables for anemones
# and anemonefish separately, saving the results to species-specific files.

library(ggplot2) # For plotting (if needed)
library(dplyr) # For data manipulation
library(vegan) #for the rda analysis
library(tools) #To add the function
source("helpers/extract_env.R") # Load function

# Define global variables
env_folder      <- "data/old_env/current"
save_location   <- "data/log"
coral_shapefile <- "data/shapefiles/WCMC008_CoralReef2018_Py_v4_1.shp"
occurrence_crs  <- "EPSG:4326"
env_crs         <- "EPSG:4326"

# Set the directory and load function
for (species in c("anemone", "anemonefish")) {
  # Set the occurrence folder and selected variables file path
  occurrence_folder <- file.path("data/occurrence", species)
  selected_vars_file <- file.path(save_location, paste0(species, "_selected_variables.txt"))
  # Check if the selected variables file exists
  if (!file.exists(selected_vars_file)) {
    cat("\nSelected variables file not found for", species, ". Skipping PCA.\n")
    next # Skip to the next species if the file doesn't exist
  }
  #Read the select file into names
  selected_variables =read.table(selected_vars_file)
  # Load and prepare environmental data:
  env_rast <- load_and_prepare_env_data(env_folder, coral_shapefile)
  if (is.null(env_rast)) {
    stop("Failed to load and prepare environmental data.")
  }
  
  #  Load occurrence data and extract environmental values:
  env_values_df <- load_occurrence_and_extract_env(occurrence_folder, env_rast, occurrence_crs)
  
  ## Subsetting Data and handling this error
  
  #Select  with select data - Check out and find variables from here -
  good_data = env_values_df%>%select(all_of(selected_variables[[1]]))
  
  # If no data is found , then return NA. This might help with future debugging
  if (nrow(good_data)<3||ncol(good_data) <3) {
    
    cat(number, "is empty")
    next
  }
  ## Load up code now
  #Run PCA analys and function
  # Run PCA - only if there are two factors to analysis for,
  
  data_stand <- decostand(good_data, method = "standardize")
  
  test =  metaMDS(data_stand , distance = "bray" ,k = 2, autotransform = FALSE)
  
  if (species == "anemone"){
    
    filename = "anemoneMDSplots.rds"
  }
  
  if (species == "anemonefish"){
    
    filename = "anenomefishMDSplots.rds"
  }
  #Create plot and call
  
  print(stressplot(test,  main =paste(species,"MDS")))
  plot= saveWidget(as_widget(stressplot(test,  main =paste(species,"MDS"))),  data_name2 )
  
  png_name <- paste0(save_location, "/",   "OrdinationGraph_",species,".png")
  png(filename= png_name )
  orditkplot(test,  main =paste(species,"MDS"))
  
  dev.off()
  
  show()
  show(test)
  cat("Code complete, now running next PCA: ",species,sep ="\n")
  next
}