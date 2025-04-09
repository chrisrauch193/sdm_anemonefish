# scripts/helpers/logging_setup.R
#-------------------------------------------------------------------------------
# Helper function to set up the log4r logger
#-------------------------------------------------------------------------------
if (!require("pacman")) install.packages("pacman")
pacman::p_load(log4r)

#' Setup log4r Logger
#'
#' Creates and configures a log4r logger instance.
#'
#' @param log_file Path to the log file.
#' @param log_level Character string log level (e.g., "INFO", "DEBUG"). Case-insensitive.
#' @param append Logical, whether to append to the log file.
#' @param log_to_console Logical, whether to also log to console.
#' @param console_level Character string log level for console output. Case-insensitive.
#' @return A log4r logger object.
setup_logger <- function(log_file, log_level = "INFO", append = TRUE,
                         log_to_console = TRUE, console_level = "INFO") {
  
  dir.create(dirname(log_file), showWarnings = FALSE, recursive = TRUE)
  
  # Define appenders
  file_appender <- log4r::file_appender(log_file, append = append)
  appenders <- list(file_appender)
  
  if (log_to_console) {
    console_appender <- log4r::console_appender()
    appenders <- c(appenders, list(console_appender))
  }
  
  # Get the numeric level for the logger threshold safely
  logger_threshold_level <- tryCatch({
    # Convert input string to uppercase to match log4r constants
    level_const_name <- toupper(log_level)
    # Access the constant value directly using get()
    get(level_const_name, envir = asNamespace("log4r"))
  }, error = function(e) {
    # Fallback to INFO if the level string is invalid
    warning("Invalid log_level '", log_level, "' specified. Defaulting to INFO.", call. = FALSE)
    log4r::INFO # Use the actual constant here as fallback
  })
  
  # Create logger
  logger <- log4r::logger(threshold = logger_threshold_level, appenders = appenders)
  
  log4r::info(logger, paste("Logger initialized. Logging to file:", log_file))
  if (log_to_console) log4r::info(logger, paste("Logging to console with minimum level:", console_level))
  
  return(logger)
}

#-------------------------------------------------------------------------------