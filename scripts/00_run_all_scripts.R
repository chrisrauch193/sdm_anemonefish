# scripts/00_run_all_scripts.R

# Force data retrieve if already existing?
force <- FALSE

source("scripts/01_install_requirements.R")

if (force | !any(grepl("final_anemone_species_list", list.files(".")))  | !any(grepl("final_anemonefish_species_list", list.files(".")))) {
  source("scripts/02_get_species_metadata.R")
}

if (force | !any(grepl(".csv", list.files("data/occurrence/anemone")))  | !any(grepl(".csv", list.files("data/occurrence/anemonefish")))) {
  source("scripts/03_download_occurrence_data.R")
}

# Don't need to check if files exist for this code
# source("scripts/04_download_env_data.R")

# Run the SDM scripts
print("Doing env variable selection")
source("scripts/05_env_variable_selection.R")

print("Doing PCA reduction")
source("scripts/06_env_variable_pca.R")

print("Doing the rest")
source("scripts/07_data_cleaning.R")
source("scripts/08_background_generation.R")
source("scripts/09_run_sdm.R")
source("scripts/run_anemonefish_sdm.R")
source("scripts/run_anemone_sdm.R")
