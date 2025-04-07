# scripts/helpers/logging_setup.R
#-------------------------------------------------------------------------------
# Setup for the logger package
#-------------------------------------------------------------------------------
if (!require("pacman")) install.packages("pacman")
pacman::p_load(logger, tools) # Ensure logger is loaded

setup_sdm_logging <- function(config, log_file_base = "sdm_run") {
  
  # Create log directory if it doesn't exist
  log_dir <- config$log_dir
  if (!dir.exists(log_dir)) {
    dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)
  }
  
  # Create a unique log file name with timestamp
  timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
  log_file_path <- file.path(log_dir, paste0(log_file_base, "_", timestamp, ".log"))
  
  # Configure the logger
  log_threshold(INFO) # Set minimum level to log (e.g., INFO, DEBUG, WARN, ERROR)
  
  # Log to both console and file
  log_appender(appender_tee(log_file_path))
  
  # Define a layout function (includes timestamp, level, message)
  log_layout(layout_glue_colors) # Or layout_glue, layout_simple, etc.
  # Example custom layout:
  # log_layout(layout_glue('{time} [{level}] {msg}'))
  
  log_info("Logging setup complete. Log file: {log_file_path}")
  
  return(log_file_path) # Return the path for reference if needed
}

#-------------------------------------------------------------------------------