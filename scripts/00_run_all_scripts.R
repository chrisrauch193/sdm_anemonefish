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

source("scripts/04_download_env_data.R")