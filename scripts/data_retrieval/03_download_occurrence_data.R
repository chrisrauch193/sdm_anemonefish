# scripts/data_retrieval/03_download_occurrence_data.R
#-------------------------------------------------------------------------------
# Download occurrence data for all species listed in the final CSV files.
# Uses the helper function `download_occurrence_data`.
#-------------------------------------------------------------------------------
cat("--- Running Script 03: Download Occurrence Data ---\n")

# Ensure config is loaded if running standalone
if (!exists("config")) {
  source("scripts/config.R")
}

# Source the helper function
helper_path <- file.path(config$helpers_dir, "occurrence_helpers.R")
if (!file.exists(helper_path)) {
  stop("Helper file not found: ", helper_path)
}
source(helper_path)

# Check if species list files exist (should be created by script 02)
if (!file.exists(config$anemone_species_list_file)) {
  stop("Anemone species list CSV not found. Run script 02 first: ", config$anemone_species_list_file)
}
if (!file.exists(config$anemonefish_species_list_file)) {
  stop("Anemonefish species list CSV not found. Run script 02 first: ", config$anemonefish_species_list_file)
}


# --- Download for Anemones ---
cat("\n--- Starting Anemone Occurrence Download ---\n")
tryCatch({
  download_occurrence_data(
    species_list_file = config$anemone_species_list_file,
    output_dir = config$anemone_occurrence_dir
  )
}, error = function(e) {
  cat("ERROR during anemone occurrence download:", e$message, "\n")
})


# --- Download for Anemonefish ---
cat("\n--- Starting Anemonefish Occurrence Download ---\n")
tryCatch({
  download_occurrence_data(
    species_list_file = config$anemonefish_species_list_file,
    output_dir = config$anemonefish_occurrence_dir
  )
}, error = function(e) {
  cat("ERROR during anemonefish occurrence download:", e$message, "\n")
})


cat("\n--- Script 03 finished. ---\n")
#-------------------------------------------------------------------------------