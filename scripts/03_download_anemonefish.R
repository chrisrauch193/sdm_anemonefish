# Download all required occurrence data from OBIS, GBIF, and FishBase.
# Based on the following guide: https://www.gbif.us/post/2024/searching-with-aphiaids/

# Generalist
# Original:
# Amphiprion clarkii
# Amphiprion ocellaris
# Amphiprion perideraion

# Extra:
# Amphiprion bicinctus
# Amphiprion chrysopterus
# Amphiprion akindynos
# Amphiprion melanopus
# Amphiprion leucokranos
# Amphiprion percula
# Amphiprion omanensis
# Amphiprion tricinctus
# Amphiprion allardi
# Amphiprion chrysogaster
# Amphiprion latifasciatus (double check, it might be some sort of specialist)


# Specialist
# Original:
# Amphiprion frenatus
# Amphiprion sandaracinos
# Amphiprion polymnus

# Extra:
# Amphiprion akallopisos
# Amphiprion ephippium
# Amphiprion sebae
# Amphiprion chagosensis
# Amphiprion nigripes
# Amphiprion latezonatus
# Amphiprion biaculeatus
# Amphiprion mccullochi
# Amphiprion rubrocinctus
# Amphiprion barberi
# Amphiprion fuscocaudatus


# STEP 1: Get final accepted list of target species
# Read in in original species list from xlsx file
original_anemonefish_list <- readxl::read_excel("Anemonefish List.xlsx")
colnames(original_anemonefish_list) <- c(
  "vernacularName",
  "scientificNameOriginal",
  "rank"
)

# Match original Scientific Names with correct WoRMS ACCEPTED Scientific Names
new_names <- obistools::match_taxa(original_anemonefish_list$scientificNameOriginal)

# Bind new names
final_anemonefish_list <- bind_cols(new_names, original_anemonefish_list)

final_anemonefish_list <- final_anemonefish_list %>%
  mutate(AphiaID = acceptedNameUsageID) %>%
  relocate(AphiaID, scientificName)

print(final_anemonefish_list)

# Print the number of species
print(paste("Number of anemonefish species:", nrow(final_anemonefish_list)))

# Add GBIF information
gbif_names <- lapply(final_anemonefish_list$scientificName, rgbif::name_backbone)
gbif_names <- bind_rows(gbif_names)
gbif_names <- gbif_names %>%
  select(gbif_speciesKey = usageKey, gbif_scientificName = scientificName, gbif_matchType = matchType)

final_anemonefish_list <- final_anemonefish_list %>% bind_cols(gbif_names)

# Display the first few rows of the combined data frame
print(head(final_anemonefish_list))

# Add taxonID column (same as AphiaID)
final_anemonefish_list$taxonID <- final_anemonefish_list$AphiaID

# Save the anemonefish species list to a CSV file
write.csv(final_anemonefish_list, "final_anemonefish_species_list.csv", row.names = FALSE)





# STEP 2: Download all data from OBIS, GBIF and other sources
# Ensure the output directory exists
output_dir <- here("data", "occurrence", "anemonefish")
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Loop through each anemonefish in the final_anemonefish_list using `seq_along`
for (i in seq_along(final_anemonefish_list$AphiaID)) {  # Loop based on the AphiaID column, assuming it's consistent
  # Fetch occurrence data from OBIS
  print(paste("Processing OBIS", i, "/", length(final_anemonefish_list$AphiaID)))
  obis_fields = c("scientificName", "decimalLongitude", "decimalLatitude", "date_year", "month", "day", "eventDate", "dataset_id", "occurrenceStatus", "depth")
  obis_results <- tryCatch({
    robis::occurrence(taxonid = final_anemonefish_list$AphiaID[i], fields = obis_fields)
  },
  error = function(e) {
    warning(paste("OBIS query failed for AphiaID", final_anemonefish_list$AphiaID[i], ":", e$message))
    return(NULL) 
  }
  )
  
  # Fetch occurrence data from GBIF
  print(paste("Processing GBIF", i, "/", length(final_anemonefish_list$AphiaID)))
  gbif_results <- tryCatch({
    obissdm::occurrence_gbif(
      taxonid = final_anemonefish_list$AphiaID[i],
      absence = "include",
      exclude = c("FOSSIL_SPECIMEN","LIVING_SPECIMEN")
    )
  },
  error = function(e) {
    warning(paste("GBIF query failed for AphiaID", final_anemonefish_list$AphiaID[i], ":", e$message))
    return(NULL)
  }
  )
  
  
  # Skip to the next iteration if either OBIS or GBIF data is NULL
  if (is.null(obis_results) || is.null(gbif_results)) {
    warning(paste("Skipping iteration", i, "due to missing OBIS or GBIF data."))
    next
  }
  
  # Select and rename columns from OBIS data
  obis_select <- obis_results %>%
    dplyr::select(scientificName, decimalLongitude, decimalLatitude, date_year, month, day, eventDate, dataset_id, occurrenceStatus, depth) %>%
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
  data_combined <- rbind(obis_select, gbif_select)
  
  # Convert eventDate to Date format
  data_combined$eventDate <- as.Date(data_combined$eventDate, format = "%Y-%m-%d")
  
  
  # Remove duplicate rows based on longitude, latitude, and eventDate
  data_deduplicated <- data_combined[!duplicated(cbind(data_combined$decimalLongitude, data_combined$decimalLatitude, data_combined$eventDate)), ]
  
  # Create a filename based on AphiaID with an "SP" prefix
  filename <- paste0(final_anemonefish_list$AphiaID[i], ".csv")
  filepath <- file.path(output_dir, filename)
  
  # Save the deduplicated data as a CSV file
  write_csv(data_deduplicated, file = filepath)
  
  # Print a message indicating the file has been saved
  print(paste("Saved data to", filepath))
}