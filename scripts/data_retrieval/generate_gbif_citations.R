# scripts/citations/01_build_gbif_citations_from_occurrence_csvs.R
#-------------------------------------------------------------------------------
# Build GBIF citations from existing occurrence CSVs (no new downloads).
# - Scans occurrence CSVs (anemone + anemonefish)
# - Extracts unique GBIF dataset keys (UUIDs)
# - Fetches official dataset citation text from GBIF
# - Writes CSV + Markdown list + a single "derived data" citation block
#
# Outputs (created under data/citations/):
#   - gbif_dataset_citations.csv
#   - gbif_dataset_citations.md
#   - gbif_citation_full.txt
#-------------------------------------------------------------------------------

cat("--- Running: Build GBIF citations from existing occurrence CSVs ---\n")

# Ensure config is loaded if running standalone
if (!exists("config")) {
  source("scripts/config.R")
  if (!exists("config")) stop("Failed to load config from scripts/config.R")
}

# Packages
suppressPackageStartupMessages({
  if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman")
})
pacman::p_load(readr, dplyr, purrr, stringr, tibble, fs, rgbif)

# --- Helper: safe directory defaults ------------------------------------------------------------
# (Falls back to conventional locations if not present in config)
safe_get <- function(lst, name, default) if (!is.null(lst[[name]])) lst[[name]] else default

occ_dir_anemone    <- safe_get(config, "anemone_occurrence_dir",     "data/occurrence/anemone")
occ_dir_fish       <- safe_get(config, "anemonefish_occurrence_dir", "data/occurrence/anemonefish")
citations_dir      <- safe_get(config, "citations_dir",              "data/citations")

dir_create_safe <- function(path) {
  if (!fs::dir_exists(path)) fs::dir_create(path, recurse = TRUE)
  path
}
dir_create_safe(citations_dir)

# --- Helper: read minimal columns robustly ------------------------------------------------------
# We try to avoid reading full CSVs. If we can't select, we read and select after.
maybe_read_min <- function(file) {
  # Try fast: read, then select only a few likely columns
  suppressMessages({
    df <- tryCatch(
      readr::read_csv(file, show_col_types = FALSE, progress = FALSE, guess_max = 1000),
      error = function(e) NULL
    )
  })
  if (is.null(df) || nrow(df) == 0) return(tibble())
  
  # Normalize column names (lower snake for matching)
  nm <- names(df)
  names(df) <- nm
  
  # Keep only columns that might contain dataset keys + source/provider info
  keep <- intersect(
    c("datasetKey", "dataset_key", "dataset_id", "gbif_datasetKey", "source", "provider", "origin"),
    names(df)
  )
  if (length(keep) == 0) return(tibble())
  dplyr::select(df, dplyr::all_of(keep))
}

# --- Gather candidate files ---------------------------------------------------------------------
csvs <- c(
  fs::dir_ls(occ_dir_anemone,  regexp = "\\.csv$", type = "file", fail = FALSE),
  fs::dir_ls(occ_dir_fish,     regexp = "\\.csv$", type = "file", fail = FALSE)
)
csvs <- unique(csvs)
cat("Found", length(csvs), "occurrence CSV(s) to scan.\n")

# --- Extract dataset keys (UUIDs) from all files ------------------------------------------------
uuid_re <- "^[0-9a-fA-F]{8}-[0-9a-fA-F]{4}-[0-9a-fA-F]{4}-[0-9a-fA-F]{4}-[0-9a-fA-F]{12}$"

extract_keys_one <- function(df) {
  if (nrow(df) == 0) return(character())
  
  # Prefer GBIF rows if a 'source' column exists
  if ("source" %in% names(df)) {
    df <- df %>% filter(tolower(source) %in% c("gbif", "gbif.org"))
  }
  
  # Coalesce anything that looks like a dataset key column
  key <- dplyr::coalesce(
    df$datasetKey,
    df$dataset_key,
    df$gbif_datasetKey,
    df$dataset_id
  )
  
  key <- key[!is.na(key)]
  key <- key[str_detect(key, uuid_re)]
  unique(key)
}

all_keys <- csvs %>%
  map(maybe_read_min) %>%
  map(extract_keys_one) %>%
  unlist(use.names = FALSE) %>%
  unique()

cat("Unique GBIF dataset keys found:", length(all_keys), "\n")
if (length(all_keys) == 0) {
  cat("No GBIF dataset keys detected. Nothing to cite.\n")
  quit(save = "no", status = 0)
}

# --- Fetch dataset metadata & citation text from GBIF -------------------------------------------
fetch_one_ds <- function(key) {
  # Robust fetch with tryCatch; sometimes GBIF returns empty
  out <- tryCatch(rgbif::datasets(uuid = key), error = function(e) NULL)
  if (is.null(out) || is.null(out$data) || nrow(out$data) == 0) {
    return(tibble(
      datasetKey = key,
      title = NA_character_,
      publishingOrganizationTitle = NA_character_,
      doi = NA_character_,
      license = NA_character_,
      type = NA_character_,
      homepage = NA_character_,
      citation_text = NA_character_
    ))
  }
  
  # Only first row is relevant for a UUID query
  d <- out$data[1, , drop = FALSE]
  
  # Citation can appear as nested list; try multiple shapes
  # Prefer full citation text if present
  citation_text <- NA_character_
  if ("citation" %in% names(d)) {
    cit <- d$citation[[1]]
    if (is.list(cit)) {
      citation_text <- cit$text %||% cit$citation %||% NA_character_
    } else if (is.character(cit)) {
      citation_text <- cit[1]
    }
  }
  
  tibble(
    datasetKey = key,
    title = d$title %||% NA_character_,
    publishingOrganizationTitle = d$publishingOrganizationTitle %||% NA_character_,
    doi = d$doi %||% NA_character_,
    license = d$license %||% NA_character_,
    type = d$type %||% NA_character_,
    homepage = d$homepage %||% NA_character_,
    citation_text = citation_text
  )
}

`%||%` <- function(a, b) if (!is.null(a)) a else b

cat("Contacting GBIF for dataset metadata & citations...\n")
gbif_meta <- all_keys %>%
  map(fetch_one_ds) %>%
  bind_rows()

n_missing <- sum(is.na(gbif_meta$citation_text) | gbif_meta$citation_text == "")
if (n_missing > 0) {
  cat("Warning:", n_missing, "dataset(s) returned without a citation text.\n")
}

# --- Write CSV output ---------------------------------------------------------------------------
csv_out <- fs::path(citations_dir, "gbif_dataset_citations.csv")
readr::write_csv(gbif_meta, csv_out)
cat("Wrote dataset table:", csv_out, "\n")

# --- Write Markdown references list -------------------------------------------------------------
md_out <- fs::path(citations_dir, "gbif_dataset_citations.md")
access_date <- as.character(Sys.Date())

md_lines <- c(
  "# GBIF dataset citations (derived data)",
  paste0("_Accessed via GBIF.org on ", access_date, "._"),
  "",
  "If you also have a GBIF **download DOI**, prefer citing that single DOI. Otherwise, cite each dataset below:",
  ""
)

# Each bullet prefers the official citation text; falls back to "Title — Publisher. DOI"
md_bullets <- gbif_meta %>%
  mutate(
    fallback = paste0(
      ifelse(is.na(title), "(Untitled dataset)", title),
      ifelse(is.na(publishingOrganizationTitle), "", paste0(" — ", publishingOrganizationTitle)),
      ifelse(is.na(doi) | doi == "", "", paste0(". DOI: https://doi.org/", doi)),
      ". (GBIF dataset key: ", datasetKey, ")"
    ),
    cite_line = ifelse(is.na(citation_text) | citation_text == "", fallback, citation_text)
  ) %>%
  pull(cite_line) %>%
  paste0("- ", .)

readr::write_lines(c(c(md_lines, md_bullets, "")), md_out)
cat("Wrote Markdown references:", md_out, "\n")

# --- Write single derived-citation block --------------------------------------------------------
txt_out <- fs::path(citations_dir, "gbif_citation_full.txt")

header <- paste0(
  "Derived from datasets accessed via GBIF.org on ", access_date, ".\n",
  "Cite the GBIF download DOI if available; otherwise include the datasets below:\n\n"
)

body <- gbif_meta %>%
  mutate(
    one_line = ifelse(is.na(citation_text) | citation_text == "",
                      paste0(title %||% "(Untitled dataset)",
                             ifelse(is.na(publishingOrganizationTitle), "", paste0(" — ", publishingOrganizationTitle)),
                             ifelse(is.na(doi) | doi == "", "", paste0(". DOI: https://doi.org/", doi)),
                             ". GBIF dataset key: ", datasetKey, "."),
                      citation_text)
  ) %>%
  pull(one_line) %>%
  paste("*", ., collapse = "\n")

readr::write_file(paste0(header, body, "\n"), txt_out)
cat("Wrote single derived-citation block:", txt_out, "\n")

cat("--- Finished: GBIF citations built with no new downloads. ---\n")
