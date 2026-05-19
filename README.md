# PlasmidNL Shiny App

## Overview

This project provides an interactive, open-source tool built in R (Shiny) for rapid analysis of plasmids and antimicrobial resistance (AMR) patterns.
The application enables users to:
- Explore genomic patterns and relationships within AMR datasets
- Perform comparative analyses across multiple variables
- Analyse spatio-temporal trends to better understand the spread of AMR plasmids through integrated visualiswations and mapping tools

A hosted graphical user interface (GUI) version is available via the Dutch [National Institute for Public Health and the Environment (RIVM)](https://apps.rivm.nl/bsr-ids-ienv/plasmidnl/),
featuring our plasmid dataset obtained from [carbapenem-producing organisms (CPO)](http://dx.doi.org/10.5281/ZENODO.18920264).
This repository provides the standalone, open-source version of the tool, allowing users to run the application locally and perform analyses on their own datasets.

A template csv to upload your own data to compare, as well as a pipeline to type plasmid genomes in the same way as the metadata set is available at the [PlasmidNL_typing repository](https://github.com/gidjes/PlasmidNL_typing)

## Installation

The following steps should be taken to install our project.

Clone the repository and change directory

```bash
    git clone https://gitlab.com/gidjes/plasmidnl.git
    cd PlasmidNL
```

## Running the Application

### 1. Open R / RStudio and install shiny

Open Rstudio (recommended) and ensure your working directory is set to project root:

```R
setwd("path/to/PlasmidNL")
install.packages(c("shiny"))
```

### 2. Open global.R

The global.R script contains all the necessary code to run the application. In addition
it will tell your R / RStudio to install all required packages / dependencies. No manual
installation is required. Simply select 'Run App' in the top corner of your RStudio or
run either:

#### From Rstudio
```R
shiny::runApp()
```

#### From R console

In the project root:
```R
shiny::runApp(.)
```

## Application Structure

The application contains several modules / Rscripts to run the application:

- **global.R**      -->     Initialises the app as well installs and loads the required packages
- **ui.R**          -->     Defines the application user interface
- **server.R**      -->     Defines the server plotting and data manipulation
- **functions.R**   -->     Additional helper functions used in the application

## Notes

The initial or first run of the application may take some extra time in order to properly install all packages.
If package installation fails, try to install all dependencies manually:

```R
install.packages(c(
    # Shiny / Dashboards
    "shiny", "shinydashboard", "shinyFiles",
    "shinythemes", "shinyWidgets", "flexdashboard",
    
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
    ))
```

# Workflow

![PlasmidNL workflow](flowchart/PlasmidNL_flowchart.png)