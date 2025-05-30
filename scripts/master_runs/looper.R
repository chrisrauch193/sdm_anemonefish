source("~/a0236995/sdm_anemonefish/scripts/config.R", echo=TRUE)

# Loop from 1 to 50
for (i in 1:1) {
  # Save the loop counter to config$global_seed
  config$global_seed <- i
  
  # Print a message to track progress
  cat(paste0("\n--- Starting iteration ", i, " with global_seed = ", config$global_seed, " ---\n"))
  
  source("~/a0236995/sdm_anemonefish/scripts/sdm_runs/06a_run_sdm_anemone.R", echo=TRUE)
  source("~/a0236995/sdm_anemonefish/scripts/sdm_runs/06b_run_sdm_anemonefish_env.R", echo=TRUE)
  source("~/a0236995/sdm_anemonefish/scripts/sdm_runs/06c_run_sdm_anemonefish_biotic_only.R", echo=TRUE)
  source("~/a0236995/sdm_anemonefish/scripts/sdm_runs/06d_run_sdm_anemonefish_combined.R", echo=TRUE)

  cat(paste0("--- Finished iteration ", i, " ---\n"))
}