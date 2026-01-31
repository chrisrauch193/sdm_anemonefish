# SDM Anemonefish

Species distribution modeling for anemonefish and their host anemones under current and future climate scenarios.

This project models the distributions of 30 anemonefish species and 10 host anemone species, accounting for host-dependence relationships and climate change projections (SSP 1.1.9 and SSP 5.8.5 for 2050 and 2100).

## Repository Structure

```
sdm_anemonefish/
├── scripts/                    # R scripts for the modeling pipeline
├── data/                       # Input data (occurrence, environmental, interaction matrix)
├── outputs/                    # Model outputs, predictions, figures
├── index.qmd                   # Quarto report document
└── index.html                  # Rendered analysis report
```

### scripts/

Run in this order for entire start to finish pipeline

| File | Description |
|------|-------------|
| `00_config.R` | Global configuration (run ID, species selection, hyperparameters) |
| `data_retrieval/` | Scripts for downloading raw occurrence and environmental data |
| `00_download_and_clean_occurrences.R` | Download occurrences from OBIS/GBIF |
| `00_prepare_final_env_dataset.R` | Environmental data preparation and PCA |
| `0a_marine_regions.R` | Marine region spatial data prep |
| `0b_env_datapreparation.R` | Environmental variable preparation |
| `0c_env_var_selection.R` | Variable selection and justification |
| `0d_occ_datapreparation.R` | Occurrence data preparation |
| `01_functions_core.R` | Utility functions (thinning, spatial folding, Boyce index) |
| `02_functions_model.R` | Model fitting, tuning, bootstrap, predictions |
| `03_pipeline_runner.R` | Main pipeline executor |
| `04_post_analysis_helpers.R` | Post-hoc analysis and visualization |


### data/

- `anem_occ_env_final_dataset.csv` - Anemone occurrences with environmental variables
- `amph_occ_env_final_dataset.csv` - Fish occurrences with environmental variables
- `interaction_matrix.csv` - Host preference matrix (30 fish x 10 anemones)
- `selected_environmental_variables.csv` - Environmental variables at occurrence locations
- `final_env_stack.tif` - Environmental raster stack (PC1-5, rugosity)
- `shapefiles/` - Spatial reference data
  - `MarineRealms_BO.*` - Marine realms boundaries
  - `meow_ecos.*` - Marine Ecoregions of the World (MEOW)
  - Coral reef shapefile - **not tracked** (see External Data Requirements below)
- `env/pca_model.rds` - Fitted PCA model used to transform environmental variables in SDM runs
- `env/` - Environmental rasters (current, future, terrain) - not tracked in git
- `occurrence/` - Raw occurrence data by species - not tracked in git

### outputs/

```
outputs/final_run/
├── models/                     # Fitted maxnet model objects
├── models_stats/               # Model evaluation metrics
├── models_tuning/              # ENMevaluate tuning results
├── predictions/
│   ├── current/
│   │   ├── hosts/              # Anemone predictions
│   │   ├── fish_env_only/      # Fish models (environment only)
│   │   ├── fish_host_only/     # Fish models (host dependency only)
│   │   └── fish_combined/      # Fish models (combined)
│   └── future/                 # SSP scenario predictions
├── occurrences/                # Thinned occurrence datasets
├── stage_cache/                # Bootstrap intermediate results
├── analysis_tables/            # Analysis outputs
└── figures_final_index/        # Generated figures
```

## External Data Requirements

Due to size constraints and licensing, the raw coral reef shapefile is not hosted in this repository. You must download it manually to run the pipeline.

**Dataset:** Global Distribution of Warm-water Coral Reefs (v4.1)
**Source:** UNEP-WCMC Ocean Data Viewer
**Download Link:** [http://data.unep-wcmc.org/datasets/1](http://data.unep-wcmc.org/datasets/1)

**Setup Instructions:**
1. Download the dataset (Format: Shapefile).
2. Unzip the archive.
3. Place the files (`WCMC008_CoralReef2018_Py_v4_1.*`) in `data/shapefiles/`.

## Running the Pipeline

### Requirements

**System:**
- R 4.5.2 (tested on Ubuntu 24.04.1 LTS)
- Quarto 1.7+ (for rendering report)
- Pandoc 3.1+

**Key R packages:**
| Package | Version |
|---------|---------|
| terra | 1.8-80 |
| sf | 1.0-22 |
| dplyr | 1.1.4 |
| ggplot2 | 4.0.0 |
| tidyr | 1.3.1 |
| purrr | 1.2.0 |
| readr | 2.1.5 |
| raster | 3.6-32 |
| lme4 | 1.1-37 |
| rnaturalearth | 1.1.0 |

Additional packages: `ENMeval`, `maxnet`, `dismo`, `doParallel`, `foreach`, `factoextra`, `patchwork`, `tidyterra`

### Configuration

Edit `scripts/00_config.R` to set:
- `RUN_ID` - output directory name (default: "final_run")
- `N_CORES` - number of parallel workers
- `N_HOST_BOOT` - bootstrap iterations for hosts (default: 10)
- `N_FISH_BOOT` - bootstrap iterations for fish (default: 40)
- `TARGET_HOSTS`, `TARGET_FISH` - species to model (NULL for all)

### Execution

Full pipeline:
```r
source("scripts/03_pipeline_runner.R")
```

Data preparation only:
```r
source("scripts/00_download_and_clean_occurrences.R")
source("scripts/00_prepare_final_env_dataset.R")
```

Render analysis report:
```r
quarto::quarto_render("index.qmd")
```

## Model Types

The pipeline fits four types of models:

1. **Host models** - Anemone distributions using environmental predictors
2. **Fish env-only** - Fish distributions using environment only
3. **Fish host-only** - Fish distributions constrained by host suitability
4. **Fish combined** - Fish distributions using both environment and host suitability

Environmental predictors are PC1-4 from a PCA on climate variables (SST, salinity, oxygen, nitrate, chlorophyll, pH) plus terrain rugosity.

## Outputs

- Raster predictions (mean and SD) for current and future climate scenarios
- Model evaluation metrics (AUC, TSS, Boyce index)
- Bootstrap ensemble statistics
- HTML report with maps and analysis figures

## License

MIT License - see LICENSE file.
