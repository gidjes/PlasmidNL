server <- function(input, output, session) {
  # Reactive dataset used globally
  expand_data <- reactive({
    if (is.null(input$userfile)) {
      return(metadata)
    } else {
      df <- open_metadata(input$userfile$datapath, "UserUpload")[[1]]
      
      # Optional: check headers
      expected_headers <- colnames(metadata)
      missing_cols <- setdiff(expected_headers, colnames(df))
      if (!all(expected_headers %in% colnames(df))) {
        showNotification("Uploaded file headers do not match expected format", type = "error")
        print(missing_cols)
        return(metadata)
      }
      
      # Combine datasets for plots
      rbind(metadata, df)
    }
  })
  
  observeEvent({
    list(input$col1, input$col2, input$col3)
  }, {
    
    picked <- c(input$col1, input$col2, input$col3)
    
    for (id in c("col1", "col2", "col3")) {
      others <- picked[picked != input[[id]]]
      available <- setdiff(selectable_cols, others)  ## <-- Change cols to selectable_cols
      
      # keep current selection if still valid
      sel <- input[[id]]
      if (!sel %in% available) {
        sel <- available[1]
      }
      
      updateSelectInput(
        session, id,
        choices = available,
        selected = sel
      )
    }
    
  }, ignoreInit = FALSE)
  
  # cluster data
  output$cluster_table <- renderDT({
    selected_cluster <- input$cluster_select
    if (!is.null(selected_cluster)) {
      
      filtered_data <- metadata[metadata$cluster == selected_cluster, ]
      
      filtered_data[] <- lapply(filtered_data, function(x) {
        if (is.character(x)) {
          gsub("<[^>]+>", "", x)
        } else {
          x
        }
      })
      
      datatable(filtered_data, options = list(pageLength = 10))
    }
  })
  
  
  # Render the specific plots
  # Homepage
  output$voorblad_map <- renderPlotly({geo_plot_ly(parent_df, nl_provinces, "submitter_province", "Number of isolates included from each province", "Isolates")})
  
  ## Genomic Overview tabs
  # Cluster plasmid counts
  ClusterFreq <- metadata %>%
    ggplot(aes(x = cluster)) +
    geom_bar(fill = "#154273", colour="#b4b4b4", linewidth=0.2, alpha=0.9) +
    labs(y = "Count (n)", x = "Cluster") +
    theme_ggrivm() +
    theme(
      legend.position = "bottom",
      axis.text.x = element_text(angle = 45, hjust = 1),
      plot.title = element_text(size= 12, color="black", face = "plain"),
    )
  
  ClusterFreqPlotly <- ggplotly(ClusterFreq)
  
  
  # GC percentage per cluster
  ClusterGC <- metadata %>%
    ggplot(aes(x = cluster, y = gc_percentage)) +
    geom_boxplot(width = 0.5,
                 fill = "#154273",
                 colour="#b4b4b4",
                 linewidth=0.2,
                 alpha=0.9,
    ) +
    labs(y = "GC Content (%)", x = "Cluster") +
    theme_ggrivm() +
    theme(
      legend.position = "bottom",
      axis.text.x = element_text(angle = 45, hjust = 0.5),
      plot.title = element_text(size= 12, color="black", face = "plain"),
    )
  
  ClusterGCPlotly <- ggplotly(ClusterGC)
  
  # BP length per cluster
  ClusterBP <- metadata %>%
    ggplot(aes(x = cluster, y = bp_length)) +
    geom_boxplot(width = 0.5,
                 fill = "#154273",
                 colour="#b4b4b4",
                 linewidth=0.2,
                 alpha=0.9,
    ) +
    labs( y = "Sequence length (bp)", x = "Cluster") +
    scale_y_log10(labels = scales::label_number()) +
    theme_ggrivm() +
    theme(
      legend.position = "bottom",
      axis.text.x = element_text(angle = 45, hjust = 1),
      plot.title = element_text(size= 12, color="black", face = "plain"),
    )
  
  ClusterBPPlotly <- ggplotly(ClusterBP)
  
  annotations <- list(
    list(
      x=-0.05,
      y=1,
      text = "Plasmid count per cluster",  
      xref = "paper",  
      yref = "paper",  
      xanchor = "left",  
      yanchor = "bottom",  
      showarrow = FALSE
    ),
    list(
      x=-0.05,
      y=0.65,
      text = "Distribution of GC-content (%) of plasmids in each cluster",  
      xref = "paper",  
      yref = "paper",  
      xanchor = "left",  
      yanchor = "bottom",  
      showarrow = FALSE
    ),
    list(
      x=-0.05,
      y=0.31,
      text = "Distribution of plasmid lengths (bp) in each cluster",  
      xref = "paper",  
      yref = "paper",  
      xanchor = "left",  
      yanchor = "bottom",  
      showarrow = FALSE
    )
  )
  
  OverviewPlots <- subplot(ClusterFreqPlotly,
                           ClusterGCPlotly,
                           ClusterBPPlotly,
                           nrows=3,
                           margin=0.04,
                           shareX = TRUE,
                           titleX=TRUE,
                           titleY=TRUE) %>%
    layout(annotations=annotations)
  
  output$Overview <- renderPlotly({OverviewPlots})
  
  ## Plasmid Types
  # Mobility per cluster
  output$ClusterMob <- renderPlotly({
    df <- expand_data()
    categorical_bar(df, "mobility", categorical_colors[0:3], "Distribution of predicted plasmid mobility", "Cluster") %>%
      layout(
        autosize = TRUE,
        margin = list(l = 60, r = 20, b = 80, t = 28),
        xaxis = list(tickangle = -45),
        legend = list(
          orientation = "h",
          y = -0.5,
          font=list(size=10)
        )
      ) %>%
      config(responsive = TRUE)
  })
  output$ClusterSpecies <- renderPlotly({
    df <- expand_data()
    categorical_bar(df, "Species", species_palette, "Distribution of plasmid species origin", "") %>%
      layout(
        autosize = TRUE,
        margin = list(l = 60, r = 20, b = 80, t = 28),
        xaxis = list(tickangle = -45),
        legend = list(
          orientation = "h",
          y = -0.5,
          font=list(size=10)
        )
      ) %>%
      config(responsive = TRUE)
  })
  output$ClusterReplicon <- renderPlotly({
    df <- expand_data()
    categorical_bar(df, "replicon", rep_palette, "Distribution of known plasmid replicons", "") %>%
      layout(
        autosize = TRUE,
        margin = list(l = 60, r = 20, b = 80, t = 28),
        xaxis = list(tickangle = -45),
        legend = list(
          orientation = "h",
          y = -0.5,
          font=list(size=10)
        )
      ) %>%
      config(responsive = TRUE)
  })
  
  
  ## AMR Data 
  output$ClusterGeneCount <- renderPlotly({
    df <- expand_data()
    count_bar(df, "amr_genes", "Frequency of antimicrobial resistance genes detected", "") %>%
      layout(
        autosize = TRUE,
        margin = list(l = 60, r = 20, b = 80, t = 28),
        xaxis = list(tickangle = -45),
        legend = list(
          orientation = "h",
          y = -0.3
        )
      ) %>%
      config(responsive = TRUE)
  })
  output$ClusterClass <- renderPlotly({
    df <- expand_data()
    count_bar(df, "amr_classes", "Distribution of resistances to associated antibiotic classes", "") %>%
      layout(
        autosize = TRUE,
        margin = list(l = 60, r = 20, b = 80, t = 28),
        xaxis = list(tickangle = -45),
        legend = list(
          orientation = "h",
          y = -0.5
        )
      ) %>%
      config(responsive = TRUE)
  })
  output$CarbaAllele <- renderPlotly({
    df <- expand_data()
    categorical_bar(df, "carba_allele", carba_palette, "Number of carbapenemase-encoding alleles",  "") %>%
      layout(
        autosize = TRUE,
        margin = list(l = 60, r = 20, b = 80, t = 28),
        xaxis = list(tickangle = -45),
        legend = list(
          orientation = "h",
          y = -0.3
        )
      ) %>%
      config(responsive = TRUE)
  })
  
  
  
  output$ClusterGeneMap <- renderPlotly({frac_heatmap(metadata, input$gene_type, "cluster", str_glue("Ratio {input$gene_type}"))})
  
  
  output$tnse_scatter <- renderPlotly({
    df <- expand_data()
    ellips_scatter(df, "cluster")
  })
  output$tsne2 <- renderPlotly({
    df <- expand_data()
    ellips_scatter(df, input$highlight)
  })
  
  
  # Epi plots
  output$TimeSeries <- renderPlotly({
    metadata <- metadata %>% 
      mutate(year = format(as.Date(sampling_date, format = "%d-%m-%Y"), "%Y")) %>%
      filter(year >= format(input$date_range[1], "%Y")) %>%
      filter(year <= format(input$date_range[2], "%Y"))
    # creates frequency table
    # Create count table using the selected column from input$subgrouping
    count.table <- metadata %>%
      group_by(year, !!sym(input$subgrouping)) %>%
      summarise(count = n(), .groups = "drop")
    
    # Plot using ggplot2
    time_series <- ggplot(
      count.table, aes(
        x = year,
        y = count,
        color = !!sym(input$subgrouping),
        group = !!sym(input$subgrouping),
        text=sprintf(
          str_glue("%s<br>Plasmids (<i>n</i>): %s<br>Year: %s"),
          !!sym(input$subgrouping),
          count,
          year
        )
      )
    ) +
      geom_line() +
      geom_point() +
      labs(
        title = paste("Count of", input$subgrouping, "by Year"),
        x = "Year",
        y = "Count",
        color = input$subgrouping
      ) +
      theme_ggrivm()
    
    time_series <- add_palette_colour(time_series, input$subgrouping, metadata)
    
    interactive_time_series <- ggplotly(time_series, tooltip="text") %>%
      layout(
        autosize = TRUE,
        margin = list(l = 60, r = 20, b = 80, t = 40),
        xaxis = list(tickangle = -45),
        legend = list(
          orientation = "h",
          y = -0.3,
          title = list(text = "")
        )
      ) %>%
      config(responsive = TRUE)
  })
  
  output$ClusterTracing <- renderPlotly({
    
    req(input$col1, input$col2, input$col3, input$cluster_trace)
    
    col1 <- as.character(input$col1)
    col2 <- as.character(input$col2)
    col3 <- as.character(input$col3)
    
    col1_palette <- palette_by_stage[[col1]]
    
    # Filter and select relevant columns
    sankey_df <- metadata %>%
      filter(cluster %in% input$cluster_trace) %>%
      select(all_of(c(col1, col2, col3))) %>%
      drop_na()
    
    counts <- sankey_df %>%
      count(.data[[col1]], .data[[col2]], .data[[col3]], name = "value")
    
    # Build nodes (ordered by axis)
    nodes <- tibble(
      name = c(
        as.character(unique(counts[[col1]])),
        as.character(unique(counts[[col2]])),
        as.character(unique(counts[[col3]]))
      ),
      stage = c(
        rep(col1, length(unique(counts[[col1]]))),
        rep(col2, length(unique(counts[[col2]]))),
        rep(col3, length(unique(counts[[col3]])))
      )
    )
    
    # Helper to map names to indices
    node_index <- function(x) match(x, nodes$name) - 1
    
    # Build links: col1 → col2
    links_1 <- counts %>%
      mutate(
        source = node_index(.data[[col1]]),
        target = node_index(.data[[col2]]),
        link_group = .data[[col1]],
        from = .data[[col1]],
        via = .data[[col2]],
        to = .data[[col3]]
      )
    
    # Build links: col2 → col3
    links_2 <- counts %>%
      mutate(
        source = node_index(.data[[col2]]),
        target = node_index(.data[[col3]]),
        link_group = .data[[col1]],
        from = .data[[col1]],
        via = .data[[col2]],
        to = .data[[col3]]
      )
    
    links <- bind_rows(
      links_1 %>% select(source, target, value, link_group, from, via, to),
      links_2 %>% select(source, target, value, link_group, from, via, to)
    ) %>%
      mutate(
        color = col1_palette[link_group],
        hover = paste0(
          "<b>", col1, ":</b> ", from, "<br>",
          "<b>", col2, ":</b> ", via, "<br>",
          "<b>", col3, ":</b> ", to, "<br>",
          "<b>Count:</b> ", value
        )
      )
    links$color <- add_alpha_vec(links$color, 0.4)

    node_colors <- mapply(
      get_node_color,
      label = nodes$name,
      stage = nodes$stage,
      SIMPLIFY = TRUE
    )
    
    # Render plotly sankey
    plot_ly(
      type = "sankey",
      orientation = "h",
      arrangement = "fixed",
      node = list(
        label = nodes$name,
        color = node_colors,
        pad = 15,
        thickness = 18,
        x = nodes$x,
        y = nodes$y
      ),
      link = list(
        source = links$source,
        target = links$target,
        value = links$value,
        color = links$color,
        customdata = links$hover,
        hovertemplate = "%{customdata}<extra></extra>"
      )
    )
  })
  
  output$ResistenceProfile <- renderPlotly({
    frac_heatmap(metadata,
                 "amr_classes",
                 "cluster",
                 str_glue("Ratio of Occurance of Resistence to Antibiotic Classes "),
                 "-",
                 cluster_filter=input$cluster_over
    )
  })
  
  output$ClusterGeneProfile <- renderPlotly({
    frac_heatmap(metadata,
                 "amr_genes",
                 "cluster",
                 str_glue("Ratio of ARG Occurance in Cluster"),
                 "-",
                 cluster_filter=input$cluster_over
    )
  })
  
  output$ClusterMetalProfile <- renderPlotly({
    frac_heatmap(metadata,
                 "metal_genes",
                 "cluster",
                 str_glue("Ratio of Metal Resistance Genes in Cluster"),
                 "-",
                 cluster_filter=input$cluster_over
    )
  })
  
  output$PlasmidMap <- renderPlotly({geo_plot_ly(metadata, nl_municiple_map, "submitter_municipality", "Number of plasmids received\n of this cluster by municipality", "Plasmids", input$cluster_over)})
  output$PopulationMap <- renderPlotly({geo_plot_ly(metadata, nl_municiple_map, "submitter_municipality", "Fraction of isolates with plasmids\n of this cluster by municipality", "Plasmids", input$cluster_over, TRUE)})
  output$ClusterTimeSeries <- renderPlotly({categorical_time_series(metadata, input$cluster_over, input$group_col)})
  output$SpeciesDistribution <- renderPlotly({treemap(metadata, "Species", "Genus", input$cluster_over, c(genus_palette, species_palette))})
  output$RepliconDistribution <- renderPlotly({treemap(metadata, "replicon", "replicon_family", input$cluster_over, c(rep_family_palette, rep_palette))})
  
  
  output$ClusterCoCluster <- renderPlotly({
    cluster_df <- rbind(
      create_normalised_co_occurance(metadata, input$cluster_co_oc, "None", input$co_occur),
      create_normalised_co_occurance(metadata, input$cluster_co_oc, input$subdivision, input$co_occur)
    )
    
    y_categories <- cluster_df %>%
      separate_rows(cluster2, sep = ",") %>%
      filter(cluster2 != "") %>%
      pull(cluster2) %>%
      unique()
    
    row_height <- 20
    plot_height <- max(300, length(y_categories) * row_height)
    
    cluster_cluster_heat <- cluster_df %>%
      ggplot(
        aes(
          x = cluster1,
          y = cluster2
        )
      ) +
      geom_tile(
        aes(
          fill = ratio,
          text = sprintf(
            str_glue(
              "%s, <br>co-occurs with: %s<br>in %.1f%% of cases"
            ),
            cluster1,
            cluster2,
            ratio * 100
          )
        ),
        color="#535353",
        linewidth=0.1,
        alpha=0.2,
        width=0.1
      ) +
      scale_fill_gradientn(colors=numerical_palette,
                           limits=c(0, 1)
      ) +
      labs(
        title = "Co-occurrence Heatmap of clusters",
        x = input$subdivision,
        y = input$co_occur
      ) +
      theme_ggrivm()  +
      theme(
        legend.position = "bottom",
        legend.text = element_text(size = 12),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text.y = element_text(angle = 45, hjust = 0, vjust = 0, size=8, margin = margin(r = 5)),
        legend.key.size = unit(0.7, "cm"),
        plot.title = element_text(size= 12, color="black", face = "plain")
      )
    plot_ly <- ggplotly(cluster_cluster_heat, height = plot_height, tooltip="text") %>%
      layout(
        showlegend = FALSE
      )
  })
  
  
}

