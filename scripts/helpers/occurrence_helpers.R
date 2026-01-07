# scripts/helpers/occurrence_helpers.R
download_occurrence_data <- function(species_list_file, output_dir) {
  # Load the species list from the CSV file
  species_list <- read.csv(species_list_file)
  
  # Ensure the output directory exists
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # Define common OBIS fields
  obis_fields <- c("scientificName", "decimalLongitude", "decimalLatitude", "date_year", "month", "day", "eventDate", "dataset_id", "occurrenceStatus", "depth")
  
  # Loop through each species in the species_list using `seq_along`
  for (i in seq_along(species_list$AphiaID)) {
    aphia_id <- species_list$AphiaID[i]
    
    # Fetch occurrence data from OBIS
    print(paste("Processing OBIS", i, "/", length(species_list$AphiaID)))
    obis_results <- robis::occurrence(taxonid = aphia_id, fields = obis_fields)
    
    # Fetch occurrence data from GBIF
    print(paste("Processing GBIF", i, "/", length(species_list$AphiaID)))
    gbif_results <- obissdm::occurrence_gbif(
      taxonid = as.numeric(aphia_id),
      absence = "exclude",
      exclude = c("FOSSIL_SPECIMEN", "LIVING_SPECIMEN")
    )
    
    # Skip to the next iteration if either OBIS or GBIF data is NULL
    if (is.null(obis_results) || is.null(gbif_results)) {
      warning(paste("Skipping iteration", i, "due to missing OBIS or GBIF data."))
      next
    }
    
    # Select and rename columns from OBIS data
    obis_select <- obis_results %>%
      dplyr::select(all_of(obis_fields)) %>% # Use all_of to ensure columns exist
      dplyr::mutate(source = "obis.org")
    
    # Select and rename columns from GBIF data
    gbif_select <- gbif_results %>%
      dplyr::select(scientificName, decimalLongitude, decimalLatitude, year, month, day, eventDate, datasetKey, occurrenceStatus, depth)
    
    # Rename GBIF columns to match OBIS
    names(gbif_select) <- obis_fields
    
    # Add a source column to GBIF data
    gbif_select <- gbif_select %>%
      dplyr::mutate(source = "gbif.org")
    
    # Combine OBIS and GBIF data
    data_combined <- tryCatch({
      rbind(obis_select, gbif_select)
    }, error = function(e) {
      warning(paste("rbind failed for AphiaID", aphia_id, ":", e$message))
      return(NULL) # or an empty dataframe if that makes more sense for your workflow
    })
    
    # If combining failed, skip this iteration
    if (is.null(data_combined)) {
      next
    }
    
    # Convert eventDate to Date format
    data_combined$eventDate <- as.Date(data_combined$eventDate, format = "%Y-%m-%d")
    
    # Remove duplicate rows based on longitude, latitude, and eventDate
    data_deduplicated <- data_combined[!duplicated(cbind(data_combined$decimalLongitude, data_combined$decimalLatitude, data_combined$eventDate)), ]
    
    # Create a filename based on AphiaID with an "SP" prefix
    filename <- paste0(aphia_id, ".csv")
    filepath <- file.path(output_dir, filename)
    
    # Save the deduplicated data as a CSV file
    write_csv(data_deduplicated, file = filepath)
    
    # Print a message indicating the file has been saved
    print(paste("Saved data to", filepath))
  }
}