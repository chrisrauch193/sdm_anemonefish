#> Step 1: Obtain species list
final_anemone_list <- read.csv("final_anemone_species_list.csv")

# Settings
options(gbif_user = "chrisrauch193jp")
options(gbif_pwd = "pepsi100")
options(gbif_email = "chrisrauch193jp@gmail.com")


start_year <- 1950

# Loop through each anemone in the final_anemone_list using `seq_along`
for (i in seq_along(final_anemone_list$AphiaID)) {  # Loop based on the AphiaID column, assuming it's consistent
  
  # Print progress messages to the console
  print(paste("Processing OBIS", i, "/", length(final_anemone_list$AphiaID)))
  
  
  # scientificName, decimalLongitude, decimalLatitude, date_year, month, day, eventDate, dataset_id, bathymetry, occurrenceStatus, sst, depth
  fields = c("scientificName", "decimalLongitude", "decimalLatitude", "date_year", "month", "day", "eventDate", "dataset_id", "occurrenceStatus", "depth")
  
  # Fetch occurrence data from OBIS
  obis_results <- tryCatch({
    robis::occurrence(taxonid = final_anemone_list$AphiaID[i], fields = fields, startdate = as.Date(paste(start_year, 01, 01, sep = "-")))
  },
  error = function(e) {
    warning(paste("OBIS query failed for AphiaID", final_anemone_list$AphiaID[i], ":", e$message))
    return(NULL)  # Return NULL if there's an error
  }
  )
  
  
  print(paste("Processing GBIF", i, "/", length(final_anemone_list$AphiaID)))
  
  # # Fetch occurrence data from GBIF
  # gbif_results <- tryCatch({
  #   rgbif::occ_data(taxonKey = final_anemone_list$gbif_speciesKey[i], hasCoordinate = TRUE, limit=10000)$data
  # },
  # error = function(e) {
  #   warning(paste("GBIF query failed for gbif_speciesKey", final_anemone_list$gbif_speciesKey[i], ":", e$message))
  #   return(NULL)  # Return NULL if there's an error
  # }
  # )
  

  gbif_results <- obissdm::occurrence_gbif(
    taxonid = final_anemone_list$AphiaID[i],
    startdate = start_year,
    absence = "include",
    exclude = c("FOSSIL_SPECIMEN","LIVING_SPECIMEN")
  )


  # Skip to the next iteration if either OBIS or GBIF data is NULL
  if (is.null(obis_results) || is.null(gbif_results)) {
    warning(paste("Skipping iteration", i, "due to missing OBIS or GBIF data."))
    next  # Skips the rest of the current iteration
  }


  # Select and rename columns from OBIS data
  obis_select <- obis_results %>%
    dplyr::select(scientificName, decimalLongitude, decimalLatitude, date_year, month, day, eventDate, dataset_id, occurrenceStatus, depth) %>%
    dplyr::mutate(source = "obis.org")
  
  
  # Select and rename columns from GBIF data
  gbif_select <- gbif_results %>%
    dplyr::select(scientificName, decimalLongitude, decimalLatitude, year, month, day, eventDate, datasetKey, occurrenceStatus, depth)
  
  
  # Rename GBIF columns to match OBIS
  names(gbif_select) <- c("scientificName", "decimalLongitude", "decimalLatitude", "date_year", "month", "day", "eventDate", "dataset_id", "occurrenceStatus")
  
  
  # Add a source column to GBIF data
  gbif_select <- gbif_select %>%
    dplyr::mutate(source = "gbif.org")
  
  # Create a filename based on AphiaID with an "SP" prefix
  filename1 <- paste0("SP", final_anemone_list$AphiaID[i], "OBIS.csv")
  filepath1 <- file.path(output_dir, filename1)
  
  # Save the deduplicated data as a CSV file
  write_csv(obis_select, file = filepath1)
  
  # Print a message indicating the file has been saved
  print(paste("Saved data to", filepath1))
  
  # Create a filename based on AphiaID with an "SP" prefix
  filename2 <- paste0("SP", final_anemone_list$AphiaID[i], "GBIF.csv")
  filepath2 <- file.path(output_dir, filename2)
  
  # Save the deduplicated data as a CSV file
  write_csv(gbif_select, file = filepath2)
  
  # Print a message indicating the file has been saved
  print(paste("Saved data to", filepath2))
}
