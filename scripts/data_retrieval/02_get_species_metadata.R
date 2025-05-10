# scripts/data_retrieval/02_get_species_metadata.R
#-------------------------------------------------------------------------------
# Download species metadata from WoRMS and GBIF for anemones and anemonefish.
# Reads species names from Excel files, matches them to accepted names/IDs,
# and saves the final lists as CSV files.
# Based on the following guide: https://www.gbif.us/post/2024/searching-with-aphiaids/
#-------------------------------------------------------------------------------

# --- Species Lists (From User Comments) ---
# Anemone species list:
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

# Anemonefish species list:
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
#-------------------------------------------------------------------------------


cat("--- Running Script 02: Get Species Metadata (User's Original Function) ---\n")

# Ensure config is loaded if running standalone
if (!exists("config")) {
  # Source the configuration file to define the 'config' list
  source("scripts/config.R")
  if (!exists("config")) {
    stop("Failed to load config object from scripts/config.R")
  }
}

# Ensure required packages are loaded (based on original function code)
pacman::p_load(readxl, dplyr, obistools, rgbif, readr)


# --- User's Original Function Definition ---
process_species_list <- function(file_path, output_file_path) {
  # Read in species list from xlsx file
  species_list <- readxl::read_excel(file_path)
  colnames(species_list) <- c(
    "vernacularName",
    "scientificNameOriginal",
    "rank"
  )
  
  # Match original Scientific Names with correct WoRMS ACCEPTED Scientific Names
  new_names <- obistools::match_taxa(species_list$scientificNameOriginal)
  
  # Bind new names
  final_list <- bind_cols(new_names, species_list)
  
  final_list <- final_list %>%
    mutate(AphiaID = acceptedNameUsageID) %>%
    relocate(AphiaID, scientificName)
  
  print(final_list)
  
  # Print the number of species
  print(paste("Number of species:", nrow(final_list)))
  
  # Add GBIF information
  gbif_names <- lapply(final_list$scientificName, rgbif::name_backbone)
  gbif_names <- bind_rows(gbif_names)
  gbif_names <- gbif_names %>%
    select(gbif_speciesKey = usageKey, gbif_scientificName = scientificName, gbif_matchType = matchType)
  
  final_list <- final_list %>% bind_cols(gbif_names)
  
  # Display the first few rows of the combined data frame
  print(head(final_list))
  
  # Add taxonID column (same as AphiaID)
  final_list$taxonID <- final_list$AphiaID
  
  # Save the species list to a CSV file
  write.csv(final_list, output_file_path, row.names = FALSE)
  
  return(final_list) # Return the processed dataframe
}
# --- End of Original Function Definition ---


# --- Process Anemone Species List ---
# Call the function using paths from the config object
cat("\n--- Processing Anemone List ---\n")
tryCatch({
  # config$anemone_list_xlsx should point to "Anemone List.xlsx"
  # config$anemone_species_list_file should point to "final_anemone_species_list.csv"
  final_anemone_list <- process_species_list(
    file_path = config$anemone_list_xlsx,
    output_file_path = config$anemone_species_list_file
  )
}, error = function(e) {
  cat("ERROR processing anemone list:", e$message, "\n")
})

# --- Process Anemonefish Species List ---
# Call the function using paths from the config object
cat("\n--- Processing Anemonefish List ---\n")
tryCatch({
  # config$anemonefish_list_xlsx should point to "Anemonefish List.xlsx"
  # config$anemonefish_species_list_file should point to "final_anemonefish_species_list.csv"
  final_anemonefish_list <- process_species_list(
    file_path = config$anemonefish_list_xlsx,
    output_file_path = config$anemonefish_species_list_file
  )
}, error = function(e) {
  cat("ERROR processing anemonefish list:", e$message, "\n")
})


cat("\n--- Script 02 finished. ---\n")
#-------------------------------------------------------------------------------