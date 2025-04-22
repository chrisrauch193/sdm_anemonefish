# --- Preamble: Cleanup and GC ---

# Ensure the fs package is available
if (!requireNamespace("fs", quietly = TRUE)) {
  install.packages("fs")
}
library(fs)

message("--- Starting Pre-Script Cleanup ---")

# 1. Define parameters for cleanup
cleanup_dir <- "."              # Current working directory
log_pattern <- "\\.log$"       # Regex for files ending in .log (case-insensitive by default in fs)
# Use glob = "*.log" if you prefer globbing patterns instead of regex

# 2. Find and delete log files
message("Searching for and deleting log files (pattern: '", log_pattern, "') in tree: ", path_abs(cleanup_dir))
tryCatch({
  # Find files matching the pattern recursively
  log_files <- fs::dir_ls(
    path = cleanup_dir,
    recurse = TRUE,        # Search recursively
    type = "file",         # Only find files (-type f)
    regexp = log_pattern,  # Filter by regex (-name '*.log')
    # glob = "*.log",      # Alternative using globbing
    fail = FALSE           # Don't error if some directories are unreadable
  )
  
  if (length(log_files) > 0) {
    message("Found ", length(log_files), " log file(s). Attempting deletion...")
    
    # Delete the files found
    # file_delete returns paths of files it FAILED to delete
    failed_deletions <- fs::file_delete(log_files)
    
    if (length(failed_deletions) == 0) {
      message("Successfully deleted ", length(log_files), " log file(s).")
    } else {
      warning("Failed to delete the following ", length(failed_deletions), " files:\n",
              paste(failed_deletions, collapse = "\n"))
    }
  } else {
    message("No log files matching the pattern found to delete.")
  }
  
}, error = function(e) {
  # Catch potential errors during file operations
  warning("An error occurred during log file cleanup: ", e$message)
}) # End tryCatch

# 3. Run Garbage Collection
message("Running garbage collection (gc())...")
gc_start_time <- Sys.time()
gc(verbose = FALSE) # Set verbose = TRUE if you want detailed GC output
gc_end_time <- Sys.time()
message("Garbage collection complete. Time taken: ", format(gc_end_time - gc_start_time))

message("--- Pre-Script Cleanup Finished ---")
message(rep("-", 40), collapse="") # Separator