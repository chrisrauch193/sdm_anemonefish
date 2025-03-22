# Download species metadata from WoRMS and GBIF for anemones and anemonefish.
# Based on the following guide: https://www.gbif.us/post/2024/searching-with-aphiaids/

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

# Function to process species list
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

# Process anemone species list
final_anemone_list <- process_species_list("Anemone List.xlsx", "final_anemone_species_list.csv")

# Process anemonefish species list
final_anemonefish_list <- process_species_list("Anemonefish List.xlsx", "final_anemonefish_species_list.csv")