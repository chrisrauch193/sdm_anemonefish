# Download all required occurrence data from OBIS, GBIF, and FishBase.
# Based on the following guide: https://www.gbif.us/post/2024/searching-with-aphiaids/

library(dplyr)
library(mregions2)
library(sf)
library(wk)
library(worrms)
library(httr2)
library(tidyr)
library(ggplot2)
library(rnaturalearth)
library(rfishbase)
library(robis)
library(rgbif)
library(stringr)
library(ggtext)


# Step 1. Query WoRMS to identify your AphiaID(s).
AphiaID <- worrms::wm_name2id(name = "Coryphaena hippurus")


# Step 2. Query WoRMS to get the equivalent Fishbase ID(s).
fishbaseID <- request(base_url = 'https://www.marinespecies.org/rest/AphiaExternalIDByAphiaID/') %>% 
  req_url_path_append(AphiaID) %>%
  req_url_query(`type` = 'fishbase') %>% 
  req_perform() %>%
  resp_body_json() %>%
  unlist()

fishbaseID <- wm_external(id = AphiaID, type = "fishbase")


# Step 3. Query GBIF to get the equivalent identifier(s) from GBIF
sourceId <- paste0('urn:lsid:marinespecies.org:taxname:', AphiaID)

response <- request(base_url = 'https://api.gbif.org/v1/species') %>% 
  req_url_query(`datasetKey` = '2d59e5db-57ad-41ff-97d6-11f5fb264527', 
                `sourceId` = 'urn:lsid:marinespecies.org:taxname:159222') %>% 
  req_perform() %>%
  resp_body_json()

GBIF_backboneID <- response$results[[1]]$nubKey


# Step 4. Query OBIS, GBIF, and FishBase for your taxa of interest.
#Make convex hull of PR EEZ from MarineRegions.org.
PR_EEZ <- mregions2::gaz_geometry(x = 33179) %>% 
  sf::st_convex_hull() %>% 
  sf::st_as_text() %>% 
  wk::wkt() %>% 
  wk::wk_orient()

obis_results <- robis::occurrence(taxonid = AphiaID, 
                                  geometry = PR_EEZ)

# CITATIONS
OBIS_metadata <- obis_results$dataset_id %>% unique() %>% 
  robis::dataset(datasetid = .)

#generate citations
OBIS_citations <- list()

for(i in 1:nrow(OBIS_metadata)){
  
  OBIS_citations[[i]] <- robis::generate_citation(title = OBIS_metadata[i,]$title,
                                                  published = OBIS_metadata[i,]$published,
                                                  url = OBIS_metadata[i,]$url,
                                                  contacts = OBIS_metadata[i,]$contacts %>% as.data.frame())
  
} %>% unlist()

# Make citation for OBIS itself:
date_accessed <- Sys.Date()
query_title <- "Occurrence records of Coryphaena hippurus (Linnaeus, 1758) in the Puerto Rico EEZ "

OBIS_citation <- paste0("OBIS (2024) ", 
                        query_title, 
                        '(Available: Ocean Biodiversity Information System. Intergovernmental Oceanographic Commission of UNESCO. www.obis.org. Accessed:', 
                        date_accessed, 
                        ')'
)


# download without DOI, to explore the data.
gbif_results <- rgbif::occ_data(taxonKey = GBIF_backboneID,
                                geometry = PR_EEZ) %>%
  .[["data"]]

#download with DOI, so I can cite the data.
rgbif::occ_download(
  user = 'chrisrauch193jp',
  email = 'chrisrauch193jp@gmail.com',
  pwd = rstudioapi::askForPassword(prompt = "GBIF Password"),
  pred_and(pred("taxonKey", GBIF_backboneID),
           pred("geometry", PR_EEZ)
  ))

GBIF_citation <- rgbif::gbif_citation(x = occ_download_meta(key = '0006266-250310093411724'))[['download']]

# FishBase citation
repro_table <- rfishbase::species(SpecCode = fishbaseID) %>%
  reproduction() %>%
  select(Species,
         ReproMode,
         Fertilization,
         Spawning, 
         RepGuild1,
         RepGuild2,
         ParentalCare)


# Step 5. Visualise Results
#Get an outline of Puerto Rico for mapping
PR <- rnaturalearth::ne_countries(country = 'Puerto Rico', 
                                  returnclass = 'sf', 
                                  scale = 'large')

# Select needed columns
obis_select <- obis_results %>% 
  select(occurrenceID, 
         decimalLatitude, 
         decimalLongitude) %>% 
  mutate(Source = 'OBIS')

gbif_select <- gbif_results %>% 
  select(occurrenceID, 
         decimalLatitude, 
         decimalLongitude) %>% 
  mutate(Source = 'GBIF')

# Join Data from GBIF and OBIS
mahi_joined <- rbind(obis_select,
                     gbif_select)

map_plot <- ggplot(PR) +
  geom_sf() +
  geom_point(data = mahi_joined,
             inherit.aes = FALSE,
             aes(x = decimalLongitude, 
                 y = decimalLatitude,
                 color = Source)) +
  
  #everything below here only serves to stylize the plot
  scale_color_manual(values = c('orange', 'skyblue')) +
  theme_bw(base_size=14) +
  theme(plot.title = ggtext::element_markdown(hjust = 0.5)) +
  guides(color = guide_legend(override.aes = list(size = 5))) +
  coord_sf(xlim = c(-64, -70)) +
  
  # This looks really complicated, but it's only style. Just making appropriate italics and line breaks.
  labs(title = "Map of _Coryphaena hippurus_ occurrences <br>in the Puerto Rico EEZ. <br>Data sourced from OBIS and GBIF.")

map_plot

repro_table %>% 
  tidyr::pivot_longer(cols = everything(),
                      names_to = "Term", 
                      values_to = "Value") %>% 
  knitr::kable(caption = "Select reproductive traits for <i>Coryphaena hippurus</i> from FishBase.")

sessionInfo() %>% print(locale = FALSE)