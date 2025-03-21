# Download all required occurrence data from OBIS, GBIF, and FishBase.
# Based on the following guide: https://www.gbif.us/post/2024/searching-with-aphiaids/

# Original:
# Entacmaea quadricolor
# Radianthus crispa (formerly Heteractis crispa)
# Radianthus magnifica (formerly Heteractis magnifica)
# Stichodactyla gigantea
# Cryptodendrum adhaesivum
# Macrodactyla doreensis
# Stichodactyla mertensii

# Extra:
# Stichodactyla haddoni
# Radianthus malu (formerly Heteractis malu)
# Radianthus aurora (formerly Heteractis aurora)


# STEP 1: Get final accepted list of target species
# Read in in original species list from xlsx file
original_anemone_list <- readxl::read_excel("Anemone List.xlsx")
colnames(original_anemone_list) <- c(
  "vernacularName",
  "scientificNameOriginal",
  "rank"
)

# Match original Scientific Names with correct WoRMS ACCEPTED Scientific Names
new_names <- obistools::match_taxa(original_anemone_list$scientificNameOriginal)

# Bind new names
final_anemone_list <- bind_cols(new_names, original_anemone_list)

final_anemone_list <- final_anemone_list %>%
  mutate(AphiaID = acceptedNameUsageID) %>%
  relocate(AphiaID, scientificName)

print(final_anemone_list)

# Print the number of species
print(paste("Number of anemone species:", nrow(final_anemone_list)))

# Add GBIF information
gbif_names <- lapply(final_anemone_list$scientificName, rgbif::name_backbone)
gbif_names <- bind_rows(gbif_names)
gbif_names <- gbif_names %>%
  select(gbif_speciesKey = usageKey, gbif_scientificName = scientificName, gbif_matchType = matchType)

final_anemone_list <- final_anemone_list %>% bind_cols(gbif_names)

# Display the first few rows of the combined data frame
print(head(final_anemone_list))

# Add taxonID column (same as AphiaID)
final_anemone_list$taxonID <- final_anemone_list$AphiaID


# Save the anemone species list to a CSV file
write.csv(final_anemone_list, "final_anemone_species_list.csv", row.names = FALSE)




# STEP 2: Download all data from OBIS, GBIF and other sources
# Ensure the output directory exists
output_dir <- here("data", "occurrence")
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Settings
options(gbif_user = "chrisrauch193jp")
options(gbif_pwd = "pepsi100")
options(gbif_email = "chrisrauch193jp@gmail.com")





# Loop through each anemone in the final_anemone_list using `seq_along`
for (i in seq_along(final_anemone_list$AphiaID)) {  # Loop based on the AphiaID column, assuming it's consistent
  
  # Print progress messages to the console
  print(paste("Processing OBIS", i, "/", length(final_anemone_list$AphiaID)))
  
  # Fetch occurrence data from OBIS
  obis_results <- tryCatch({
    robis::occurrence(taxonid = final_anemone_list$AphiaID[i], absence="include")
  },
  error = function(e) {
    warning(paste("OBIS query failed for AphiaID", final_anemone_list$AphiaID[i], ":", e$message))
    return(NULL)  # Return NULL if there's an error
  }
  )
  
  
  print(paste("Processing GBIF", i, "/", length(final_anemone_list$AphiaID)))
  
  # Fetch occurrence data from GBIF
  gbif_results <- tryCatch({
    rgbif::occ_data(taxonKey = final_anemone_list$gbif_speciesKey[i], hasCoordinate = TRUE, limit=10)$data
  },
  error = function(e) {
    warning(paste("GBIF query failed for gbif_speciesKey", final_anemone_list$gbif_speciesKey[i], ":", e$message))
    return(NULL)  # Return NULL if there's an error
  }
  )
  
  # Skip to the next iteration if either OBIS or GBIF data is NULL
  if (is.null(obis_results) || is.null(gbif_results)) {
    warning(paste("Skipping iteration", i, "due to missing OBIS or GBIF data."))
    next  # Skips the rest of the current iteration
  }
  
  
  # Select and rename columns from OBIS data
  obis_select <- obis_results %>%
    dplyr::select(scientificName, decimalLongitude, decimalLatitude, date_year, month, day, eventDate, dataset_id, bathymetry, occurrenceStatus, sst) %>%
    dplyr::mutate(source = "obis.org")
  
  
  # Select and rename columns from GBIF data
  gbif_select <- gbif_results %>%
    dplyr::select(scientificName, decimalLongitude, decimalLatitude, year, month, day, eventDate, datasetKey)
  
  
  # Add missing columns to GBIF data to match OBIS
  gbif_select$bathymetry <- NA
  gbif_select$occurrenceStatus <- NA
  gbif_select$sst <- NA
  
  
  # Rename GBIF columns to match OBIS
  names(gbif_select) <- c("scientificName", "decimalLongitude", "decimalLatitude", "date_year", "month", "day", "eventDate", "dataset_id", "bathymetry", "occurrenceStatus", "sst")
  
  
  # Add a source column to GBIF data
  gbif_select <- gbif_select %>%
    dplyr::mutate(source = "gbif.org")
  
  
  # Combine OBIS and GBIF data
  data_combined <- rbind(obis_select, gbif_select)
  
  
  # Rename dataset_id column to datasetKey for consistency
  data_combined <- data_combined %>%
    dplyr::rename(datasetKey = dataset_id)
  
  
  # Convert eventDate to Date format
  data_combined$eventDate <- as.Date(data_combined$eventDate, format = "%Y-%m-%d")
  
  
  # Remove duplicate rows based on longitude, latitude, and eventDate
  data_deduplicated <- data_combined[!duplicated(cbind(data_combined$decimalLongitude, data_combined$decimalLatitude, data_combined$eventDate)), ]
  
  
  # Create a filename based on AphiaID with an "SP" prefix
  filename <- paste0("SP", final_anemone_list$AphiaID[i], ".csv")
  filepath <- file.path(output_dir, filename)
  
  # Save the deduplicated data as a CSV file
  write_csv(data_deduplicated, file = filepath)
  
  # Print a message indicating the file has been saved
  print(paste("Saved data to", filepath))
}






# 
# 
# # Step 1. Query WoRMS to identify your AphiaID(s).
# AphiaID <- worrms::wm_name2id(name = "Coryphaena hippurus")
# 
# 
# # Step 2. Query WoRMS to get the equivalent Fishbase ID(s).
# fishbaseID <- request(base_url = 'https://www.marinespecies.org/rest/AphiaExternalIDByAphiaID/') %>% 
#   req_url_path_append(AphiaID) %>%
#   req_url_query(`type` = 'fishbase') %>% 
#   req_perform() %>%
#   resp_body_json() %>%
#   unlist()
# 
# fishbaseID <- wm_external(id = AphiaID, type = "fishbase")
# 
# 
# # Step 3. Query GBIF to get the equivalent identifier(s) from GBIF
# sourceId <- paste0('urn:lsid:marinespecies.org:taxname:', AphiaID)
# 
# response <- request(base_url = 'https://api.gbif.org/v1/species') %>% 
#   req_url_query(`datasetKey` = '2d59e5db-57ad-41ff-97d6-11f5fb264527', 
#                 `sourceId` = 'urn:lsid:marinespecies.org:taxname:159222') %>% 
#   req_perform() %>%
#   resp_body_json()
# 
# GBIF_backboneID <- response$results[[1]]$nubKey
# 
# 
# # Step 4. Query OBIS, GBIF, and FishBase for your taxa of interest.
# #Make convex hull of PR EEZ from MarineRegions.org.
# PR_EEZ <- mregions2::gaz_geometry(x = 33179) %>% 
#   sf::st_convex_hull() %>% 
#   sf::st_as_text() %>% 
#   wk::wkt() %>% 
#   wk::wk_orient()
# 
# obis_results <- robis::occurrence(taxonid = AphiaID, 
#                                   geometry = PR_EEZ)
# 
# # CITATIONS
# OBIS_metadata <- obis_results$dataset_id %>% unique() %>% 
#   robis::dataset(datasetid = .)
# 
# #generate citations
# OBIS_citations <- list()
# 
# for(i in 1:nrow(OBIS_metadata)){
#   
#   OBIS_citations[[i]] <- robis::generate_citation(title = OBIS_metadata[i,]$title,
#                                                   published = OBIS_metadata[i,]$published,
#                                                   url = OBIS_metadata[i,]$url,
#                                                   contacts = OBIS_metadata[i,]$contacts %>% as.data.frame())
#   
# } %>% unlist()
# 
# # Make citation for OBIS itself:
# date_accessed <- Sys.Date()
# query_title <- "Occurrence records of Coryphaena hippurus (Linnaeus, 1758) in the Puerto Rico EEZ "
# 
# OBIS_citation <- paste0("OBIS (2024) ", 
#                         query_title, 
#                         '(Available: Ocean Biodiversity Information System. Intergovernmental Oceanographic Commission of UNESCO. www.obis.org. Accessed:', 
#                         date_accessed, 
#                         ')'
# )
# 
# 
# # download without DOI, to explore the data.
# gbif_results <- rgbif::occ_data(taxonKey = GBIF_backboneID,
#                                 geometry = PR_EEZ) %>%
#   .[["data"]]
# 
# #download with DOI, so I can cite the data.
# rgbif::occ_download(
#   user = 'chrisrauch193jp',
#   email = 'chrisrauch193jp@gmail.com',
#   pwd = rstudioapi::askForPassword(prompt = "GBIF Password"),
#   pred_and(pred("taxonKey", GBIF_backboneID),
#            pred("geometry", PR_EEZ)
#   ))
# 
# GBIF_citation <- rgbif::gbif_citation(x = occ_download_meta(key = '0006266-250310093411724'))[['download']]
# 
# # FishBase citation
# repro_table <- rfishbase::species(SpecCode = fishbaseID) %>%
#   reproduction() %>%
#   select(Species,
#          ReproMode,
#          Fertilization,
#          Spawning, 
#          RepGuild1,
#          RepGuild2,
#          ParentalCare)
# 
# 
# # Step 5. Visualise Results
# #Get an outline of Puerto Rico for mapping
# PR <- rnaturalearth::ne_countries(country = 'Puerto Rico', 
#                                   returnclass = 'sf', 
#                                   scale = 'large')
# 
# # Select needed columns
# obis_select <- obis_results %>% 
#   select(occurrenceID, 
#          decimalLatitude, 
#          decimalLongitude) %>% 
#   mutate(Source = 'OBIS')
# 
# gbif_select <- gbif_results %>% 
#   select(occurrenceID, 
#          decimalLatitude, 
#          decimalLongitude) %>% 
#   mutate(Source = 'GBIF')
# 
# # Join Data from GBIF and OBIS
# mahi_joined <- rbind(obis_select,
#                      gbif_select)
# 
# map_plot <- ggplot(PR) +
#   geom_sf() +
#   geom_point(data = mahi_joined,
#              inherit.aes = FALSE,
#              aes(x = decimalLongitude, 
#                  y = decimalLatitude,
#                  color = Source)) +
#   
#   #everything below here only serves to stylize the plot
#   scale_color_manual(values = c('orange', 'skyblue')) +
#   theme_bw(base_size=14) +
#   theme(plot.title = ggtext::element_markdown(hjust = 0.5)) +
#   guides(color = guide_legend(override.aes = list(size = 5))) +
#   coord_sf(xlim = c(-64, -70)) +
#   
#   # This looks really complicated, but it's only style. Just making appropriate italics and line breaks.
#   labs(title = "Map of _Coryphaena hippurus_ occurrences <br>in the Puerto Rico EEZ. <br>Data sourced from OBIS and GBIF.")
# 
# map_plot
# 
# repro_table %>% 
#   tidyr::pivot_longer(cols = everything(),
#                       names_to = "Term", 
#                       values_to = "Value") %>% 
#   knitr::kable(caption = "Select reproductive traits for <i>Coryphaena hippurus</i> from FishBase.")
# 
# sessionInfo() %>% print(locale = FALSE)