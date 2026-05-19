### UI
ui <- tagList(
  tags$html(lang="en", #language attribute for accessibility
  dashboardPage(
  
  ## Define Header
  dashboardHeader(title = "PlasmidNL"),
  
  ## Define sidebar
  dashboardSidebar(
    
    # Upload field
    dashboardSidebar(
      fileInput(
        "userfile",
        label = tags$span("Upload CSV file", `aria-label` = "Upload CSV file"),
        accept = c(
          "text/csv",
          "text/comma-separated-values,text/plain",
          ".csv"
        )
      ),
      checkboxInput("header", "Header", TRUE),
      
      # Menu items
      sidebarMenu(
        menuItem("Home",
                 tabName = "home", icon = icon("home")
        ),
        menuItem("All Clusters Overview",
                 tabName = "genomes", icon = icon("dna"),
                 menuSubItem("Overview", tabName = "genomes_overview"),
                 menuSubItem("Plasmid types", tabName = "plasmid_type"),
                 menuSubItem("AMR data", tabName = "AMR_data"),
                 menuSubItem("Gene data", tabName = "genomic_gene_data"),
                 menuSubItem("Time series", tabName = "time_series"),
                 menuSubItem("mge tsne scatter", tabName = "tnse_scatter"),
                 menuSubItem("interactive mge tsne", tabName = "tsne2")
        ),
        menuItem("Cluster in-depth",
                 tabName = "clusters", icon = icon("microscope"),
                 menuSubItem("Isolate-level Overview", tabName = "isolate_overview"),
                 menuSubItem("Genomic-level Overiew", tabName = "gen_overview"),
                 menuSubItem("Metadata Sankey", tabName = "cluster_tracing"),
                 menuSubItem("Co-occurance", tabName = "cluster_co_occurance"),
                 menuSubItem("Cluster Table", tabName = "cluster_data")
        )
      ))),
  
  ## Define the actual content
  dashboardBody(
    # Add a tag list for accessibility attributes
    tags$head(
      # add RIVM blue color as nav bar color, to prevent too-low contrast for accessibility
      tags$style(HTML('
      .main-header .navbar {
        background-color: #007bc7 !important;
        }
      .main-header .logo {
          background-color: #007bc7 !important;
          color: #fff !important;
        }
      ')
      ),
      # fonts of the menu for bigger contrast
      tags$style(HTML('
        /* Sidebar menu item font color */
        .skin-blue .main-sidebar .sidebar .sidebar-menu a {
          color: #ffffff !important; 
        }
        /* Optional: active/selected menu in bold */
        .main-sidebar .sidebar .sidebar-menu .active a {
          font-weight: bold !important;
        }
    '))
    ),
    # Homepage
    tabItems(
      tabItem(tabName = "home", fluidRow(column(
        width = 12, tags$h1("PlasmidNL", style = "color: #0000FF;"),
        h2("Antimicrobial Resistance Plasmids in the Netherlands"),
        p("Carbapenem resistance spreads largely due to horizontal gene transfer, which is mostly driven by plasmids."),
        p("This Shiny app provides insight into antimicrobial resistance plasmids in the Netherlands through interactive maps and visualizations, allowing users to analyze resistance patterns and trends, and better understand plasmid spread."),
        p("Genomic typing of plasmids is performed through the PlasmidNL_typing pipeline, available through https://gitlab.com/mmb-umcu/plasmidnl_typing"),
        tags$hr(),
        plotlyOutput("voorblad_map"),
        tags$hr(),
        tags$h3("An interactive open source software tool for rapid analysis of plasmids and antimicrobial resistance patterns in R (AMR/R)."),
        tags$hr(),
        h3("Key Features"),
        HTML("<ul>
                        <li>Gain insight into genomic patterns.</li>
                        <li>Conduct comparative analyses across different variables.</li>
                        <li>Interact with geographical data to understand the spread of AMR plasmids through visualizations and maps.</li>
                        <li>Upload and analyze plasmid data: Upload metadata from plasmid isolates to compare and understand wider patterns of antimicrobial resistance in plasmids. Templete metadata file for uploads is also available at the PlasmidNL_typing Github https://gitlab.com/mmb-umcu/plasmidnl_typing</li>
                      </ul>"),
        tags$hr(),
        h3("Data Sources"),
        p("The data used in this application is sourced from the Dutch national CPE/CPPA/CRAb surveillance performed by the RIVM (Rijksinstituut voor Volksgezondheid
          en Milieu) (National Institute for Public Health and the Environment). Isolates are provided by participating healthcare institutes.
          It includes anonimysed epidemiological data, demographic information, and microbial data of resistant isolates."),
        tags$hr(),
        h3("Disclaimer"),
        p("This application is for research purposes only and should not be used for making clinical decisions. Always consult with healthcare professionals for medical advice and treatment."),
        tags$hr(),
      ))),
      
      ## Dataset tabs
      # Genomics
      tabItem(
        tabName = "genomes_overview",
        fluidRow(
          column(12, plotlyOutput("Overview", height = "90vh", width = "100%"))
        ),
      ),
      tabItem(
        tabName = "plasmid_type",
        fluidRow(
          column(12, plotlyOutput("ClusterSpecies", height = "30vh", width = "100%"))
        ),
        fluidRow(
          column(12, plotlyOutput("ClusterReplicon", height = "30vh", width = "100%"))
        ),
        fluidRow(
          column(12, plotlyOutput("ClusterMob", height = "30vh", width = "100%"))
        )
      ),
      tabItem(
        tabName = "AMR_data",
        fluidRow(
          column(12, plotlyOutput("ClusterGeneCount", height = "30vh", width = "100%"))
        ),
        fluidRow(
          column(12, plotlyOutput("ClusterClass", height = "30vh", width = "100%"))
        ),
        fluidRow(
          column(12, plotlyOutput("CarbaAllele", height = "30vh", width = "100%"))
        )
      ),
      tabItem(
        tabName = "genomic_gene_data",
        fluidRow(
          box(
            width = 12, # or any other appropriate value
            selectInput("gene_type",
                        "Select column",
                        choices = c(
                          "amr_genes",
                          "amr_classes",
                          "carba_allele",
                          "virulence_genes",
                          "metal_genes",
                          "metal_classes",
                          "biocide_genes",
                          "heat_genes",
                          "acid_genes"
                        ),
                        selected = "amr_classes", multiple = FALSE),
            column(12, plotlyOutput("ClusterGeneMap", height = "95vh", width = "100%"))
          )
        )
      ),
      tabItem(
        tabName = "tnse_scatter",
        fluidRow(
          column(12, plotlyOutput("tnse_scatter", height = "95vh", width = "100%"))
        )
      ),
      tabItem(
        tabName = "tsne2",
        fluidRow(
          box(
            width = 12, # or any other appropriate value
            selectInput("highlight", "Select column", choices = c("cluster", selectable_cols), selected = "cluster", multiple = FALSE),
            column(12, plotlyOutput("tsne2", height = "80vh", width = "100%"))
          )
        )
      ),
      
      # Epidiomology
      tabItem(
        tabName = "time_series",
        
        fluidPage(
          
          # Full-width instructions
          fluidRow(
            column(
              width = 12,
              tags$h3("Instructions"),
              tags$p("Select a grouping and date range for plasmid samples")
            )
          ),
          
          # Inputs side by side
          fluidRow(
            column(
              width = 6,
              selectInput(
                "subgrouping",
                label = "Subgrouping",
                choices = list(
                  "cluster",
                  "Species",
                  "Genus",
                  "mobility",
                  "replicon",
                  "replicon_family",
                  "submitter_province",
                  "submitter_municipality",
                  "AMR_plasmid",
                  "CP_plasmid"
                ),
                multiple = FALSE,
                selected = "cluster"
              )
            ),
            column(
              width = 6,
              dateRangeInput(
                "date_range",
                label = "Date Range",
                start = "2012-11-01",
                end = "2023-12-31",
                format = "yyyy"
              )
            )
          ),
          
          # Plot below inputs
          fluidRow(
            column(
              width = 12,
              plotlyOutput("TimeSeries", height = "70vh")
            )
          )
          
        )
      ),
      
      ## Cluster insight tabs
      # Group 1
      tabItem(
        tabName = "isolate_overview",
        fluidRow(
          box(
            width = 12,
            selectInput("cluster_over",
                        label="Select cluster",
                        choices=unique(metadata$cluster),
                        selected = "13",
                        multiple = FALSE
            )
          )
        ),
        fluidRow(
          splitLayout(cellWidths = c("50%", "50%"),
                      plotlyOutput("PlasmidMap"),
                      plotlyOutput("PopulationMap")
          )
        ),
        fluidRow(column(12, box("Distribution of cluster plasmid species origins",
                                width = 12,
                                plotlyOutput("SpeciesDistribution")
        )
        )
        ),
        fluidRow(
          column(
            12,
            box(
              "Proportion resistance to antimicrobial classes associated with detected antimicrobial resistance genes",
              width = 12,
              plotlyOutput("ResistenceProfile", height="10%", width="100%")
            )
          )
        ),
        fluidRow(
          box(
            width = 12,
            selectInput("group_col",
                        label="Select histogram subgroups to display",
                        choices=selectable_cols,
                        selected = "Genus",
                        multiple = FALSE
            )
          )
        ),
        fluidRow(column(12, box("Counts of cluster plasmids detected over time",
                                width = 12,
                                plotlyOutput("ClusterTimeSeries")
        )
        )
        )
      ),
      
      # Group 2
      tabItem(
        tabName = "gen_overview",
        fluidRow(
          box(
            width = 12,
            selectInput("cluster_over",
                        label="Select cluster",
                        choices=unique(metadata$cluster),
                        selected = "13",
                        multiple = FALSE
            )
          )
        ),
        fluidRow(column(12, box("Distribution known replicons for cluster plasmids",
                                width = 12,
                                plotlyOutput("RepliconDistribution")
        )
        )
        ),
        fluidRow(
          column(
            12,
            box(
              "",
              width = 12,
              plotlyOutput("ClusterGeneProfile", height="10%", width="100%")
            )
          )
        ),
        fluidRow(
          column(
            12,
            box(
              "",
              width = 12,
              plotlyOutput("ClusterMetalProfile", height="10%", width="100%")
            )
          )
        )
      ),
      
      # Sankey
      tabItem(
        tabName = "cluster_tracing",
        fluidRow(
          box(
            width = 12,
            selectInput("cluster_trace",
                        label="Select cluster",
                        choices=unique(metadata$cluster),
                        selected = "37",
                        multiple = FALSE
            ),
            selectInput("col1",
                        label="Select 1st column",
                        choices=selectable_cols,
                        selected = "Species",
                        multiple = FALSE
            ),
            selectInput("col2",
                        label="Select 2nd column",
                        choices=selectable_cols,
                        selected = "carba_allele",
                        multiple = FALSE
            ),
            selectInput("col3",
                        label="Select 3rd column",
                        choices=selectable_cols,
                        selected = "foreign_hospitalisation_history",
                        multiple = FALSE
            )
          )
        ),
        fluidRow(column(12, plotlyOutput("ClusterTracing", height = "95vh", width = "90%")))
      ),
      
      ## Table tab
      tabItem(
        tabName = "cluster_data",
        fluidPage(
          box(
            width = 12, # or any other appropriate value
            selectInput("cluster_select", "Select Cluster", choices = unique(metadata$cluster), multiple = TRUE),
            div(DTOutput("cluster_table"), style = "width: 100%; overflow-x: auto;")
          )
        )
      ),
      
      ## Co-occurance
      tabItem(
        tabName = "cluster_co_occurance",
        fluidRow(
          box(
            width = 12,
            selectInput("cluster_co_oc",
                        label="Select cluster",
                        choices=unique(metadata$cluster),
                        selected = "13",
                        multiple = FALSE
            ),
            selectInput("subdivision", "Select cluster subdivision", choices = selectable_cols, selected = "replicon", multiple = FALSE),
            selectInput("co_occur", "Select co-occurance group", choices = c("cluster", selectable_cols), selected = "cluster", multiple = FALSE),
          )
        ),
        fluidRow(column(12,
                        plotlyOutput("ClusterCoCluster", height = "95vh", width = "95%")
        )
        )
      )
    )
  )
)
)
)