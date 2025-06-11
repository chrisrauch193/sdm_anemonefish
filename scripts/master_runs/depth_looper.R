source("~/a0236995/sdm_anemonefish/scripts/config.R", echo=TRUE)


# depths <- c(50, 300, 500, 750, 1000)
depths <- c(500)

# Loop from 1 to 50
for (depth in depths) {
  # Print a message to track progress
  cat(paste0("\n[", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "] --- Starting iteration depth: ", depth, " with global_seed = ", config$global_seed, " ---\n"))  
  
  # Save the loop counter to config$global_seed
  config$depth_min <- -depth
  config$do_final_prediction <- TRUE
  
  
  for (i in 10:50) {
    # Print a message to track progress
    cat(paste0("\n[", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "] --- Run ", i, " for depth: ", depth, " with global_seed = ", config$global_seed, " ---\n"))
    
    config$global_seed <- depth + i
    
    
    source("~/a0236995/sdm_anemonefish/scripts/sdm_runs/06a_run_sdm_anemone.R", echo=TRUE)
    source("~/a0236995/sdm_anemonefish/scripts/sdm_runs/06b_run_sdm_anemonefish_env.R", echo=TRUE)
    source("~/a0236995/sdm_anemonefish/scripts/sdm_runs/06c_run_sdm_anemonefish_biotic_only.R", echo=TRUE)
    source("~/a0236995/sdm_anemonefish/scripts/sdm_runs/06d_run_sdm_anemonefish_combined.R", echo=TRUE)
  }
  
  cat(paste0("[", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "] --- Finished iteration depth: ", depth, " ---\n"))
}
