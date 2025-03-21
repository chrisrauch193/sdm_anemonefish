############################# Get species lists ################################
# Based on the new gridded product: https://github.com/iobis/speciesgrids

# Load packages ----
library(readr)
library(h3jsr)
library(sf)
library(duckdb)
library(DBI)
library(dplyr)
library(worrms)
library(purrr)
library(rgbif)

sf_use_s2(FALSE)

print("ALright hwgg")

# Load study area
starea <- st_read("data/shapefiles/mpa_asia_starea_v1.shp")
# Resolution Guide: https://h3geo.org/docs/core-library/restable/
cells <- data.frame(cell = polygon_to_cells(starea, 7)[[1]]) # takes some time

# Set up DuckDB connection
con <- dbConnect(duckdb())
dbSendQuery(con, "install httpfs; load httpfs;")
duckdb_register(con, "cells", cells)

species <- dbGetQuery(con, "
  select species, AphiaID, source_obis, source_gbif
  from cells
  inner join read_parquet('s3://obis-products/speciesgrids/h3_7/*') h3 on cells.cell = h3.cell
  group by species, AphiaID, source_obis, source_gbif
")

head(species)

# Define the target species (Anemones and Anemonefish)
target_species <- c(
  # Anemone Species
  "Cryptodendrum adhaesivum",
  "Entacmaea quadricolor",
  "Heteractis crispa",
  "Heteractis magnifica",
  "Macrodactyla doreensis",
  "Stichodactyla gigantea",
  "Stichodactyla mertensii",
  # Anemonefish Species
  "Amphiprion clarkii",
  "Amphiprion frenatus",
  "Amphiprion ocellaris",
  "Amphiprion perideraion",
  "Amphiprion polymnus",
  "Amphiprion sandaracinos"
)
# Filter species to only the target taxa
species <- species %>% filter(species %in% target_species)

# Add WoRMS taxonomy
id_batches <- split(species$AphiaID, ceiling(seq_along(species$AphiaID) / 50))
taxa_batches <- purrr::map(id_batches, worrms::wm_record, .progress = T)

species_list <- bind_rows(taxa_batches) %>% 
  select(AphiaID, scientificname, kingdom, phylum, class, order, family, genus, scientificName = scientificname, rank, authority,
         status, valid_AphiaID, valid_name, isMarine, isBrackish, isTerrestrial, isFreshwater)

species_list <- species_list %>%
  filter(rank == "Species") %>%
  distinct(AphiaID, .keep_all = T) %>%
  # filter(!kingdom %in% c("Archaea", "Bacteria", "Fungi", "Protozoa", "Viruses", "Biota incertae sedis")) %>%
  # filter(!is.na(isMarine) | !is.na(isBrackish) | !is.na(isTerrestrial) | !is.na(isFreshwater)) %>%
  # mutate(isFreshwater = ifelse(is.na(isFreshwater), 0, isFreshwater)) %>%
  # rowwise() %>%
  # mutate(isAnyMBT = ifelse(all(is.na(c(isBrackish, isTerrestrial, isMarine))), 0,
  #                          ifelse(sum(na.omit(c(isBrackish, isTerrestrial, isMarine))) > 0, 1, 0))) %>%
  # filter(isAnyMBT == 1 | isFreshwater == 0)

head(species_list)
nrow(species_list)
length(unique(species_list$AphiaID))

# Join OBIS/GBIF information
species_source <- species %>% 
  select(AphiaID, source_obis, source_gbif) %>%
  distinct(AphiaID, .keep_all = T)

species_list_or <- left_join(species_list, species_source, by = "AphiaID")

species_list_or$source_both <- ifelse(species_list_or$source_obis + species_list_or$source_gbif == 2,
                                      TRUE, FALSE)

sum(species_list_or$source_obis)
sum(species_list_or$source_gbif)
sum(species_list_or$source_both)

# Add GBIF keys
gbif_keys <- rgbif::name_backbone_checklist(species_list_or, strict = T)

gbif_keys <- gbif_keys %>%
  select(gbif_scientificName = canonicalName, gbif_speciesKey = usageKey,
         gbif_order = order)

species_list_final <- bind_cols(species_list_or, gbif_keys) %>%
  select(-rank)
head(species_list_final)

# Add taxonID for compatibility with other parts
species_list_final <- species_list_final %>%
  mutate(taxonID = AphiaID)

# Save final list
write.csv(species_list_final, 
          paste0("data/all_splist_", format(Sys.Date(), "%Y%m%d"), ".csv"),
          row.names = F)

DBI::dbDisconnect(con)
