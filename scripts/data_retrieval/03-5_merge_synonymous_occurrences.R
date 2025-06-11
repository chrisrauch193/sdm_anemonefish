# scripts/data_processing/03-5_merge_synonymous_occurrences.R
#-------------------------------------------------------------------------------
# Merge occurrence data for synonymous species:
# - Heteractis aurora (from 290088.csv) occurrences will be merged into
#   Radianthus aurora (206704.csv).
# - Scientific names will be standardized to "Radianthus aurora".
# - Duplicates based on scientific name, location (rounded coordinates), and
#   eventDate will be removed.
# - The original 206704.csv file will be overwritten with the merged,
#   standardized, and deduplicated data.
#-------------------------------------------------------------------------------

cat("--- Running Script 04: Merge Synonymous Occurrences ---\n")

# --- Configuration and Setup ---
# Load necessary packages
if (!requireNamespace("readr", quietly = TRUE)) install.packages("readr")
if (!requireNamespace("dplyr", quietly = TRUE)) install.packages("dplyr")

library(readr)
library(dplyr)

# Define file paths
# Assuming the script is run from a project root where 'scripts/' and 'data/' directories exist.
base_dir <- "." # Assumes current working directory is the project root
occurrence_dir <- file.path(base_dir, "data", "occurrence", "anemone")

file_radianthus_aurora_source <- file.path(occurrence_dir, "206704.csv")
file_heteractis_aurora_source <- file.path(occurrence_dir, "290088.csv")
output_file_radianthus_aurora <- file.path(occurrence_dir, "206704.csv") # This file will be overwritten

# Define target scientific name
target_scientific_name <- "Radianthus aurora"

# Define expected column types to ensure consistent reading
# Most are character; numerics and dates will be converted after reading
col_spec <- cols(
  scientificName = col_character(),
  decimalLongitude = col_character(),
  decimalLatitude = col_character(),
  date_year = col_character(),
  month = col_character(),
  day = col_character(),
  eventDate = col_character(),
  dataset_id = col_character(),
  occurrenceStatus = col_character(),
  depth = col_character(),
  source = col_character()
)

# --- Read and Prepare Data ---
cat("Reading and preparing occurrence data...\n")

# Function to read and preprocess a single CSV
read_and_prep_occurrences <- function(file_path, desired_scientific_name) {
  if (!file.exists(file_path)) {
    stop("Source file not found: ", file_path, "\nMake sure this file exists and the path is correct.")
  }
  
  # read_csv will interpret "NA" as NA by default if not quoted.
  # Using na = c("", "NA", "NA ", "na", "N/A") for robustness.
  df <- read_csv(file_path, col_types = col_spec, na = c("", "NA", "NA ", "na", "N/A"))
  
  df <- df %>%
    mutate(
      # Standardize to the target scientific name
      scientificName = desired_scientific_name,
      # Convert types
      decimalLongitude = as.numeric(decimalLongitude),
      decimalLatitude = as.numeric(decimalLatitude),
      date_year = as.integer(date_year),
      month = as.integer(month),
      day = as.integer(day),
      # Parse eventDate carefully: readr::parse_date returns NA for unparseable strings with a warning.
      # This handles cases where eventDate might be empty or already NA from read_csv.
      eventDate = readr::parse_date(eventDate, format = "%Y-%m-%d", na = c("", "NA")),
      depth = as.numeric(depth)
    )
  return(df)
}

# Read Radianthus aurora data (from 206704.csv)
# This file might contain various names initially, so we standardize it.
df_radianthus <- read_and_prep_occurrences(file_radianthus_aurora_source, target_scientific_name)
cat("  Read", nrow(df_radianthus), "rows from", basename(file_radianthus_aurora_source), "\n")

# Read Heteractis aurora data (from 290088.csv)
df_heteractis <- read_and_prep_occurrences(file_heteractis_aurora_source, target_scientific_name)
cat("  Read", nrow(df_heteractis), "rows from", basename(file_heteractis_aurora_source), "\n")


# --- Combine Dataframes ---
cat("Combining dataframes...\n")
# The column names and types should now be consistent due to read_and_prep_occurrences
combined_df <- bind_rows(df_radianthus, df_heteractis)
cat("Total rows after combining: ", nrow(combined_df), "\n")

# --- Deduplication ---
cat("Deduplicating occurrences...\n")
# Round coordinates for robust matching (e.g., 5 decimal places ≈ 1.1m precision)
# This helps catch duplicates that only differ by very minor floating point variations.
deduplicated_df <- combined_df %>%
  mutate(
    # Create temporary rounded columns for deduplication key.
    # If coordinates are NA, they remain NA. distinct() treats NAs as actual values for grouping.
    rounded_longitude = if_else(is.na(decimalLongitude), NA_real_, round(decimalLongitude, 5)),
    rounded_latitude  = if_else(is.na(decimalLatitude), NA_real_, round(decimalLatitude, 5))
  ) %>%
  # Deduplicate based on scientificName, rounded location, and eventDate.
  # .keep_all = TRUE retains all columns of the first unique row encountered.
  # The order of rows in combined_df (df_radianthus first) means records from
  # the original 206704.csv are prioritized if an exact duplicate exists.
  dplyr::distinct(scientificName, rounded_longitude, rounded_latitude, eventDate, .keep_all = TRUE) %>%
  # Remove the temporary rounded coordinate columns
  dplyr::select(-rounded_longitude, -rounded_latitude)

cat("Total rows after deduplication: ", nrow(deduplicated_df), "\n")
cat(nrow(combined_df) - nrow(deduplicated_df), "duplicate rows removed.\n")


# --- Final Checks (Optional but good practice) ---
# Ensure all scientificNames in the final dataframe are the target scientific name.
if (nrow(deduplicated_df) > 0 && any(deduplicated_df$scientificName != target_scientific_name, na.rm = TRUE)) {
  warning("Warning: Some records still do not have the target scientific name '", target_scientific_name, "'. This is unexpected.")
  # Example to show differing names:
  # print(filter(deduplicated_df, scientificName != target_scientific_name) %>% select(scientificName) %>% distinct())
}

# --- Write Output ---
# Ensure output columns are in the original order and format eventDate back to character.
# Get the original column order from one of the input spec (they are the same)
original_col_order <- names(col_spec$cols)

deduplicated_df_final <- deduplicated_df %>%
  # Reorder columns to match original input
  dplyr::select(all_of(original_col_order)) %>%
  # Format eventDate back to character "YYYY-MM-DD" for CSV, "NA" for missing dates.
  mutate(eventDate = if_else(is.na(eventDate), NA_character_, format(eventDate, "%Y-%m-%d")))


cat("Writing merged and deduplicated data to: ", output_file_radianthus_aurora, "\n")
# write_csv uses "" for NA by default. To write "NA" string for NAs, matching input:
write_csv(deduplicated_df_final, output_file_radianthus_aurora, na = "NA")

cat("--- Script 04 finished. ---\n")