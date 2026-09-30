# Install packages if not already installed
options(repos = c(CRAN = "https://cran.r-project.org"))

## -------------------------------------------------
## Package management
## -------------------------------------------------

required_packages <- c(
  # Shiny / Dashboards
  "shiny", "shinydashboard", "shinyFiles", "shinythemes",
  "shinyWidgets", "flexdashboard",
  
  # Data manipulation
  "tidyverse", "data.table", "lubridate", "reshape2", "terra",
  
  # Visualization
  "ggplot2", "cowplot", "plotly", "DT", "leaflet",
  "ggalluvial", "ggforce", "treemapify", "scales",
  "RColorBrewer", "colorspace", "gridExtra", "ellipse",
  "grDevices",
  
  # Mapping / spatial
  "sf", "rnaturalearth", "rnaturalearthdata",
  
  # Trees / widgets
  "collapsibleTree", "htmlwidgets",
  
  # I/O / external data
  "here", "readr", "jsonlite", "rentrez", "config",
  
  # Explicit tidyverse components you loaded
  "stringr"
)

# Load packages
invisible(lapply(required_packages, library, character.only = TRUE))
# Source imports
source("functions.R")


# Load (main) data
here::i_am("global.R")
config <- config::get()
metadata_path <- config$METADATA_FILE
dfs <- open_metadata(str_glue("{metadata_path}", sep=""), "Reference")
metadata <- dfs[[1]]
parent_df <- dfs[[2]]
metadata_full <- dfs[[3]]

# Sort mge_cluster order by frequency
value_counts <- table(metadata$mge_cluster)
sorted_categories <- names(sort(value_counts, decreasing = TRUE))
metadata$mge_cluster <- factor(metadata$mge_cluster, levels = sorted_categories, ordered = TRUE)

## Set up geo-data
nl_map_path <- config$NL_MAP_FILE
nl_map <- load_map(nl_map_path, config)

region_col <- config$MAP_COLUMNS$region_type %||% "regio_soort"

nl_municiple_map <- filter_map_type(
  nl_map,
  c(
    config$MAP_TYPES$municipalities,
    config$MAP_TYPES$extras
  )
)

nl_provinces <- filter_map_type(
  nl_map,
  c(
    config$MAP_TYPES$provinces,
    config$MAP_TYPES$extras
  )
)

# count samples per province
sample_counts <- metadata %>%
  group_by(submitter_province) %>%
  summarise(n_plasmids = n())
provinces <- left_join(nl_provinces, sample_counts, by = c("regio_naam" = "submitter_province"))


selectable_cols <- c(
  "Species","Genus", "ST", "replicon","replicon_family","mobility",
  "submitter_province","submitter_municipality","foreign_hospitalisation_history",
  "amr","amr_classes","AMR_plasmid","carba_allele", "CP_plasmid", "metal", "virulence",
  "biocide", "heat", "acid", "sampling_source", "DataSource"
)
