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
  "tidyverse", "data.table", "lubridate", "reshape2",
  
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
  "here", "readr", "jsonlite", "rentrez", "cbsodataR",
  
  # Explicit tidyverse components you loaded
  "stringr"
)

# Install missing packages
missing_packages <- required_packages[
  !required_packages %in% installed.packages()[, "Package"]
]

if (length(missing_packages) > 0) {
  install.packages(missing_packages)
}

if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes")
}
if (!requireNamespace("DARAvis", quietly = TRUE)) {
  remotes::install_gitlab("dara/DARAvis@main", host = "https://gitlab.rivm.nl", build = FALSE)
}

# Load packages
invisible(lapply(required_packages, library, character.only = TRUE))
library(DARAvis)
# Source imports
source("functions.R")

if (.Platform$OS.type == "windows") {
  windowsFonts(Verdana=windowsFont("Verdana"))
}



# Load (main) data
here::i_am("global.R")
config <- config::get()
metadata_path <- config$METADATA_FILE
print(metadata_path)
dfs <- open_metadata(str_glue("{metadata_path}", sep=""), "Reference")
metadata <- dfs[[1]]
parent_df <- dfs[[2]]
metadata_full <- dfs[[3]]

# Sort cluster order by frequency
value_counts <- table(metadata$cluster)
sorted_categories <- names(sort(value_counts, decreasing = TRUE))
metadata$cluster <- factor(metadata$cluster, levels = sorted_categories, ordered = TRUE)

## Set up geo-data
nl_map_path <- config$NL_MAP_FILE
nl_map <- load_map(nl_map_path, config)

region_col <- config$MAP_COLUMNS$region_type %||% "regio_soort"
nl_municiple_map <- nl_map %>%
  filter(.data[[region_col]] %in% c(
    config$MAP_TYPES$municipalities,
    config$MAP_TYPES$extras
  ))

nl_provinces <- nl_map %>%
  filter(.data[[region_col]] %in% c(
    config$MAP_TYPES$provinces,
    config$MAP_TYPES$extras
  ))

# count samples per province
sample_counts <- metadata %>%
  group_by(submitter_province) %>%
  summarise(n_plasmids = n())
provinces <- left_join(nl_provinces, sample_counts, by = c("regio_naam" = "submitter_province"))


selectable_cols <- c(
  "Species","Genus", "ST", "replicon","replicon_family","mobility",
  "submitter_province","submitter_municipality","foreign_hospitalisation_history",
  "amr_genes","amr_classes","AMR_plasmid","carba_allele", 
  "CP_plasmid", "metal_genes", "virulence_genes", "sampling_source", "sampling_OH_domain", "DataSource", "Institute"
)

