# scripts/citations/build_gbif_citations_from_existing_csvs.R
# ------------------------------------------------------------------------------
# Build GBIF citations from EXISTING occurrence CSVs (no new occurrence downloads).
# Scans anemone + anemonefish CSVs, extracts unique GBIF dataset keys (UUIDs),
# fetches official GBIF dataset citation text/DOI via rgbif::datasets(), and writes:
#   - data/citations/gbif_dataset_citations.csv
#   - data/citations/gbif_dataset_citations.md
#   - data/citations/gbif_citation_full.txt
# ------------------------------------------------------------------------------

suppressPackageStartupMessages({
  if (!requireNamespace("readr", quietly = TRUE)) install.packages("readr")
  if (!requireNamespace("dplyr", quietly = TRUE)) install.packages("dplyr")
  if (!requireNamespace("stringr", quietly = TRUE)) install.packages("stringr")
  if (!requireNamespace("purrr", quietly = TRUE)) install.packages("purrr")
  if (!requireNamespace("tibble", quietly = TRUE)) install.packages("tibble")
  if (!requireNamespace("fs", quietly = TRUE)) install.packages("fs")
  if (!requireNamespace("rgbif", quietly = TRUE)) install.packages("rgbif")
})

library(readr); library(dplyr); library(stringr); library(purrr)
library(tibble); library(fs); library(rgbif)

# ---- CONFIG (auto-detect config.R, else defaults) ----------------------------------------------
ANEMONE_DIR     <- "data/occurrence/anemone"
ANEMONEFISH_DIR <- "data/occurrence/anemonefish"
OUTDIR          <- "data/citations"

if (file.exists("scripts/config.R")) {
  source("scripts/config.R")
  if (exists("config")) {
    if (!is.null(config$anemone_occurrence_dir))     ANEMONE_DIR     <- config$anemone_occurrence_dir
    if (!is.null(config$anemonefish_occurrence_dir)) ANEMONEFISH_DIR <- config$anemonefish_occurrence_dir
    if (!is.null(config$citations_dir))              OUTDIR          <- config$citations_dir
  }
}

if (!fs::dir_exists(OUTDIR)) fs::dir_create(OUTDIR, recurse = TRUE)

# ---- HELPERS -----------------------------------------------------------------------------------
uuid_re <- "^[0-9a-fA-F]{8}-[0-9a-fA-F]{4}-[0-9a-fA-F]{4}-[0-9a-fA-F]{4}-[0-9a-fA-F]{12}$"
`%||%` <- function(a,b) if (!is.null(a)) a else b

list_csvs <- function(paths) {
  files <- character(0)
  for (p in paths) {
    if (fs::is_dir(p)) {
      files <- c(files, fs::dir_ls(p, regexp = "\\.csv$", type = "file", fail = FALSE))
    } else if (fs::file_exists(p) && grepl("\\.csv$", p, ignore.case = TRUE)) {
      files <- c(files, p)
    }
  }
  unique(files)
}

maybe_read_min <- function(file) {
  df <- tryCatch(readr::read_csv(file, show_col_types = FALSE, progress = FALSE, guess_max = 1000),
                 error = function(e) NULL)
  if (is.null(df) || !nrow(df)) return(tibble())
  keep <- intersect(
    c("datasetKey","dataset_key","gbif_datasetKey","dataset_id","source","origin"),
    names(df)
  )
  if (!length(keep)) return(tibble())
  dplyr::select(df, dplyr::all_of(keep))
}

extract_keys <- function(df) {
  if (!nrow(df)) return(character(0))
  # Prefer rows marked as GBIF
  if ("source" %in% names(df)) {
    df <- df %>% filter(tolower(source) %in% c("gbif","gbif.org"))
  }
  # Coalesce likely key columns
  key <- dplyr::coalesce(df$datasetKey, df$dataset_key, df$gbif_datasetKey, df$dataset_id)
  key <- key[!is.na(key) & stringr::str_detect(key, uuid_re)]
  base::unique(as.character(key))
}

fetch_dataset_meta <- function(key) {
  out <- tryCatch(rgbif::datasets(uuid = key), error = function(e) NULL)
  if (is.null(out) || is.null(out$data) || !nrow(out$data)) {
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
  d <- out$data[1, , drop = FALSE]
  cit_txt <- NA_character_
  if ("citation" %in% names(d)) {
    c <- d$citation[[1]]
    if (is.list(c)) cit_txt <- c$text %||% c$citation %||% NA_character_
    else if (is.character(c)) cit_txt <- c[1]
  }
  tibble(
    datasetKey = key,
    title = d$title %||% NA_character_,
    publishingOrganizationTitle = d$publishingOrganizationTitle %||% NA_character_,
    doi = d$doi %||% NA_character_,
    license = d$license %||% NA_character_,
    type = d$type %||% NA_character_,
    homepage = d$homepage %||% NA_character_,
    citation_text = cit_txt
  )
}

# ---- RUN ---------------------------------------------------------------------------------------
csvs <- list_csvs(c(ANEMONE_DIR, ANEMONEFISH_DIR))
message("Scanning ", length(csvs), " CSV file(s).")

keys <- csvs %>%
  map(maybe_read_min) %>%
  map(extract_keys) %>%
  list_flatten_chr() %>%
  base::unique()

message("Unique GBIF dataset keys found: ", length(keys))
if (!length(keys)) {
  message("No GBIF dataset keys detected. Exiting.")
  quit(save = "no", status = 0)
}

meta <- keys %>% map(fetch_dataset_meta) %>% list_rbind()

# CSV
csv_out <- fs::path(OUTDIR, "gbif_dataset_citations.csv")
readr::write_csv(meta, csv_out)
message("Wrote: ", csv_out)

# Markdown list
access_date <- as.character(Sys.Date())
fallback_line <- function(r) {
  paste0(
    ifelse(is.na(r$title), "(Untitled dataset)", r$title),
    ifelse(is.na(r$publishingOrganizationTitle),"", paste0(" — ", r$publishingOrganizationTitle)),
    ifelse(is.na(r$doi) || r$doi=="", "", paste0(". DOI: https://doi.org/", r$doi)),
    ". (GBIF dataset key: ", r$datasetKey, ")"
  )
}
md_out <- fs::path(OUTDIR, "gbif_dataset_citations.md")
md_header <- c(
  "# GBIF dataset citations (derived data)",
  paste0("_Accessed via GBIF.org on ", access_date, "._"),
  "",
  "If you created a **GBIF download DOI**, prefer citing that single DOI. Otherwise include the per-dataset citations below:",
  ""
)
md_bullets <- apply(meta, 1, function(row) {
  r <- as.list(row)
  line <- ifelse(is.na(r$citation_text) || r$citation_text=="", fallback_line(r), r$citation_text)
  paste0("- ", line)
})
readr::write_lines(c(md_header, md_bullets, ""), md_out)
message("Wrote: ", md_out)

# Single derived-citation block (TXT)
txt_out <- fs::path(OUTDIR, "gbif_citation_full.txt")
txt_header <- paste0(
  "Derived from datasets accessed via GBIF.org on ", access_date, ".\n",
  "Cite the GBIF download DOI if you generated one; otherwise include the dataset citations below:\n\n"
)
txt_body <- apply(meta, 1, function(row) {
  r <- as.list(row)
  if (!is.na(r$citation_text) && r$citation_text!="") return(paste("*", r$citation_text))
  paste("*", fallback_line(r))
})
readr::write_file(paste0(txt_header, paste(txt_body, collapse = "\n"), "\n"), txt_out)
message("Wrote: ", txt_out)
