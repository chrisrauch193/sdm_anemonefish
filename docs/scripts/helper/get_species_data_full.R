########################### Download species data ##############################


# Load packages ----
library(robis)
library(rgbif)
library(obissdm)
library(arrow)
library(dplyr)
fs::dir_create("data/raw/")

# TODO: Download OBIS subset data with api
# Get data from OBIS ----
# api_call <- httr::content(httr::GET("https://api.obis.org/export?complete=true"), as = "parsed")
# api_call <- api_call$results[[1]]
# latest_export <- paste0("https://obis-datasets.ams3.digitaloceanspaces.com/", api_call$s3path)
# 
# options(timeout = 999999999)
# download.file(
#     url = latest_export,
#     destfile = paste0("data/raw/", gsub("exports/", "", api_call$s3path)),
#     method = "wget"
# )


# Get data from GBIF ----
# sp_list_file <- recent_file("data", "all_splist")
# sp_list <- read.csv(sp_list_file)
sp_list <- read.csv("data/all_splist_20241213.csv")

# Download through the GBIF API ----
# Because we already have the GBIF taxonKey we can supply this to the function
# Download and save
gbif_res <- mp_get_gbif(sci_names = na.omit(sp_list$gbif_speciesKey))

if (inherits(gbif_res, "occ_download")) {
    done <- FALSE
    trials <- 0
    while (!done & trials < 50) {
        if (interactive()) {
            if (!exists("timetry")) timetry <- as.numeric(readline("Chose time to try again, in minutes\n")) * 60
        } else {
            timetry <- 60 * 60
        }
        cat("Waiting... \n")
        Sys.sleep(timetry)
        cat("Trying to download\n")
        gbif_status <- try(rgbif::occ_download_meta(gbif_res)$status)
        if (gbif_status == "SUCCEEDED") {
            dl <- rgbif::occ_download_get(gbif_res, path = "data/raw", curlopts = list(timeout_ms = 10000))
            unzip(paste0("data/raw/", gbif_res, ".zip"),
                exdir = "data/raw/"
            )

            file.rename(
                "data/raw/occurrence.parquet",
                glue::glue("data/raw/gbif_{save_acro}_{format(Sys.Date(), '%Y%m%d')}.parquet")
            )
            fs::file_delete(paste0(
                "data/raw/", gbif_res,
                ".zip"
            ))

            write.table(
                data.frame(
                    date = attr(
                        gbif_res,
                        "created"
                    ), download_code = as.character(gbif_res),
                    doi = attr(gbif_res, "doi"), citation = attr(
                        gbif_res,
                        "citation"
                    )
                ), glue::glue("data/gbif_full_download_log.txt"),
                row.names = F
            )

            done <- TRUE
        } else if (gbif_status == "FAILED") {
            stop("Problem in the download.\n")
        } else {
            trials <- trials + 1
        }
    }
}

### END