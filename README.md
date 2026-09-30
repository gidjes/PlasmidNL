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

Clone the repository and change to the project directory:

```bash
    git clone https://gitlab.com/gidjes/plasmidnl.git
    cd PlasmidNL
```

The project uses ```renv``` to manage its R package dependencies and versions.

Open the project in R/RStudio. If ```renv``` is not already installed, install it once and then restore the project environment from the lockfile::

```R
install.packages("renv")
library(renv)
renv::restore()
```

This will install the package versions specified in renv.lock and recreate the project's R package environment.


## Running the Application

### 1. Open R / RStudio and install shiny

Once the renv environment has been restored, the application can be run as usual.
Open Rstudio (recommended) and ensure your working directory is set to project root:

```R
setwd("path/to/PlasmidNL")
install.packages(c("shiny"))
```

### 2. In Rstudio
Open the projects ```global.R``` script in RStudio and select Run App.

#### From R console
In the project root:
```R
shiny::runApp()
```

## Application Structure

The application contains several modules / Rscripts to run the application:

- **global.R**      -->     Initialises the app as well installs and loads the required packages
- **ui.R**          -->     Defines the application user interface
- **server.R**      -->     Defines the server plotting and data manipulation
- **functions.R**   -->     Additional helper functions used in the application

In addition, two files are used as data input in the application:
- **metadata.csv**              --> contains the data that is visualised
- **NLenBESenCAS_2024.json**    --> contains the map coordinate data

Both these files are placed in the shiny-data directory can be replaced with your own data to adapt the application to your circumstances.
See the instruction below about further details.

```
shiny-data
├── metadata.csv
└── NLenBASenCAS_2024.json
```

## Working with your own data
### Uploading your own file
You can use the [CSV template](upload_template) provided in this repository to upload and compare your own data in the dashboard.

If you are working with your own plasmid genome sequences, the [PlasmidNL_typing repository](https://github.com/gidjes/PlasmidNL_typing) provides a pipeline for typing and annotating plasmid genomes in the same consistent format as the genomic data used in this dashboard. The resulting annotations can then be combined with the metadata in the [CSV template](upload_template.csv) provided in this repository.

For help filling out the template, a [pre-filled example](upload_example_records.csv) containing a few reference records is provided in the repository. A [detailed description](upload_inputs_detailed.csv) of each metadata field, including the expected format and how each field is interpreted by the dashboard, is also available and is shown below.

| Column name                       | Data description                                                                                                                                   | Data type (R) | Example record                 |
| --------------------------------- | -------------------------------------------------------------------------------------------------------------------------------------------------- | ------------- | ------------------------------ |
| `Plasmid`                         | Unique identifier for the plasmid record.                                                                                                          | `character`   | `PLASMID_00123`                |
| `replicon`                        | Replicon type identified in the plasmid.                                                                                                           | `character`   | `IncFII`                       |
| `replicon_family`                 | Replicon family or group to which the replicon belongs.                                                                                            | `character`   | `IncF`                         |
| `mobility`                        | Predicted plasmid mobility category.                                                                                                               | `character`   | `conjugative`                  |
| `mge_cluster`                     | Cluster identifier from mge-cluster scheme.                                                                                           | `character`   | `26`                     |
| `tsne1D`                          | Coordinate of the plasmid on the first t-SNE dimension.                                                                                            | `numeric`     | `12.47`                        |
| `tsne2D`                          | Coordinate of the plasmid on the second t-SNE dimension.                                                                                           | `numeric`     | `-3.82`                        |
| `amr`                             | Antibiotic resistance genes detected in the plasmid. Enter multiple genes as a comma-separated list. Leave empty if no genes were detected.        | `character`   | `blaKPC-2,aac(6')-Ib-cr`    |
| `amr_classes`                     | Antibiotic classes for which resistance genes were detected. Enter multiple classes as a comma-separated list. Leave empty if none were detected.  | `character`   | `beta-lactam,aminoglycoside` |
| `AMR_plasmid`                     | Indicates whether the plasmid is classified as an antimicrobial resistance (AMR) plasmid.                                                          | `numeric`     | `1`                         |
| `carba_allele`                    | Carbapenem resistance genes or alleles detected in the plasmid. Enter multiple genes as a comma-separated list. Leave empty if none were detected. | `character`   | `blaKPC-2`                     |
| `CP_plasmid`                      | Indicates whether the plasmid is classified as a carbapenem-resistance plasmid.                                                                    | `numeric`     | `1`                         |
| `virulence`                       | Virulence-associated genes detected in the plasmid. Enter multiple genes as a comma-separated list. Leave empty if none were detected.             | `character`   | `iucC,iutA`                    |
| `metal`                           | Metal resistance genes detected in the plasmid. Enter multiple genes as a comma-separated list. Leave empty if none were detected.                 | `character`   | `merC`                    |
| `metal_classes`                   | Metals for which resistance was detected. Enter multiple metals as a comma-separated list. Leave empty if none were detected.                      | `character`   | `mercury`                |
| `biocide`                         | Biocide resistance genes detected in the plasmid. Enter multiple genes as a comma-separated list. Leave empty if none were detected.               | `character`   | `qacE`                         |
| `heat`                            | Heat resistance genes detected in the plasmid. Enter multiple genes as a comma-separated list. Leave empty if none were detected.                  | `character`   | `hsp20`                        |
| `acid`                            | Acid resistance genes detected in the plasmid. Enter multiple genes as a comma-separated list. Leave empty if none were detected.                  | `character`   | `gadA,gadB`                    |
| `GC_perc`                         | G+C content of the plasmid sequence, expressed as a percentage.                                                                                    | `numeric`     | `51.23`                        |
| `bp_length`                       | Length of the plasmid nucleotide sequence in base pairs.                                                                                           | `integer`     | `85432`                        |
| `Parent`                          | Identifier of the isolate record associated with the plasmid.                                                                                      | `character`   | `ISO_00456`                    |
| `ST`                              | Sequence type (ST) of the associated isolate.                                                                                                      | `character`   | `131`                        |
| `Species`                         | Species of the associated isolate.                                                                                                                 | `character`   | `Escherichia coli`             |
| `sampling_date`                   | Date on which the associated isolate was sampled. Use the format `YYYY-MM-DD`.                                                                     | `Date`        | `2024-03-15`                   |
| `submitter_municipality`          | Municipality associated with the submitting record, corresponding to the lower NUTS/local administrative level.                                    | `character`   | `Bilthoven`                      |
| `submitter_province`              | Province associated with the submitting record, corresponding to the higher NUTS administrative level.                                             | `character`   | `Utrecht`                      |
| `foreign_hospitalisation_history` | Country in which the patient was hospitalised within the 6 months preceding sampling. Leave empty if there was no foreign hospitalisation.         | `character`   | `Denmark`                      |
| `healthcare_employee`             | Indicates whether the patient works in healthcare.                                                                                                 | `logical`     | `1`                         |


### Using your own map
The app uses an `sf` spatial file (GeoJSON, GeoPackage, shapefile, etc.) to draw the map layers.
You can replace the default map of the Netherlands with your own regional map data by editing the configuration file.

#### 1. Add your map file
Place your map file in the shiny-data directory so it accessible to the app. It is not required to delete the Dutch map:

```
shiny-data
├── metadata.csv
├── NLenBASenCAS_2024.json
└── my_map.geojson
```

Supported formats include:
- .geojson
- .json
- .gpkg
- shapefiles (.shp)

**The file must be readable by the sf package.**

#### 2. Update the config file
If using a multi-layer map file, make sure `region_type` matches the column determining the layer. If using a single-layer map, simple omit this line entirely.

edit `config.yml`:

```yaml
NL_MAP_FILE: shiny-data/my_map.geojson

MAP_COLUMNS:
    region_type: column_region_type
    region_name: column_region_name

MAP_TYPES:
    municipalities: municipality
    provinces: province
    extras: boundary
```

- Replace `column_region_type` with the name of the column defining the layer types.
- Replace `column_region_name` with the name of the column defining the region names matching those in `submitter_municipaliy` in metadata.csv
- Replace `municipality` with variable name of the layer type you want to plot.
- Replace `province` with the variable name defining the provinces.
- Replace `boundary` with the variable name definig boundries. Used to plot inset borders.

##### Example file
| regio_naam | regio_soort | geometry |
| -------- | ------- | ------- |
| Utrecht | province | POLYGON(...) |
| Noord-Brabant | province | POLYGON(...) |
| Bilthoven | municipality | POLYGON(...) |
| Meierijstad | municipality | POLYGON(...) |

#### 3. Notes
- Coordinate reference systems (CRS) are handled automatically by sf.
- Large map files may increase startup time.
- MultiPolygon geometries are supported.
- If your dataset does not contain province/municipality distinctions, you may use the same value for multiple categories.

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

## Workflow

![PlasmidNL workflow](flowchart/PlasmidNL_flowchart.png)