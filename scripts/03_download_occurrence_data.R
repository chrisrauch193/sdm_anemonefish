# scripts/03_download_occurrence_data.R
# Download occurrence data from OBIS and GBIF for anemones and anemonefish.
# Based on the following guide: https://www.gbif.us/post/2024/searching-with-aphiaids/

source("helpers/occurrence_helpers.R")

# Run the function for anemones and anemonefish
download_occurrence_data("final_anemone_species_list.csv", "data/occurrence/anemone")
download_occurrence_data("final_anemonefish_species_list.csv", "data/occurrence/anemonefish")