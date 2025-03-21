#> Step 1: Obtain species list
final_anemone_list <- read.csv("final_anemone_species_list.csv")

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

# # Download from OBIS
# obissdm::mp_get_obis(sci_names = final_anemone_list$taxonID, mode = "full")
# 
# # Download from GBIF
# obissdm::mp_get_gbif(sci_names = final_anemone_list$gbif_speciesKey)

# Split datasets into individual species
# Select columns to keep on files
gbif_columns <- c("gbifid", "datasetkey", "occurrenceid", "species", "decimallatitude",
                  "decimallongitude", "depth", "depthaccuracy", "day", "month",
                  "year", "taxonkey", "specieskey", "basisofrecord", "catalognumber")

obis_columns <- c("occurrenceID", "catalogNumber", "fieldNumber", "aphiaID",
                  "materialSampleID", "institutionID", "collectionID", "datasetID",
                  "collectionCode", "institutionCode", "datasetName",
                  "eventID", "decimalLatitude", "decimalLongitude", 
                  "species", "eventDate", "date_year", "day", "month", "year",
                  "occurrenceStatus", "flags", "depth", "maximumDepthInMeters", 
                  "minimumDepthInMeters")

# Split GBIF dataset

# GBIF species are identified by a taxonKey, instead of AphiaID from  WoRMS used by OBIS.
gbif_new_keys <- final_anemone_list[,c("gbif_speciesKey", "AphiaID")]
colnames(gbif_new_keys) <- c("key", "new_key")

obissdm::split_dataset("data/raw/gbif_full_20250320.parquet",
              database_name = "gbif",
              grouping_key = "specieskey",
              sel_keys = gbif_new_keys$new_key,
              change_key = gbif_new_keys,
              sel_columns = gbif_columns,
              run_in_batches = T,
              batch_size = 100)

# # Split OBIS dataset
# split_dataset("data/raw/obis_full_20250320.parquet",
#               database_name = "obis",
#               grouping_key = "aphiaID",
#               sel_keys = final_anemone_list$taxonID,
#               sel_columns = obis_columns,
#               run_in_batches = T)

# Check
sum(grepl("ftype=gbif", list.files("data/species", recursive = T)))
sum(grepl("ftype=obis", list.files("data/species", recursive = T)))

