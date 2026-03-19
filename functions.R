create_palettes <- function(metadata, categorical_colors) {
  numerical_palette <- c(
    "#FFFFFF",
    palette_rivm("full")[1],
    palette_rivm("categorical")[2],
    palette_rivm("categorical")[3]
  )
  
  cluster_palette <- custom_hierarchical_palette(metadata, "cluster", "cluster", categorical_colors)
  species_palette <- custom_hierarchical_palette(metadata, "Genus", "Species", categorical_colors)
  ST_palette <- custom_hierarchical_palette(metadata, "Genus", "ST", categorical_colors)
  geo_palette <- custom_hierarchical_palette(metadata, "submitter_province", "submitter_municipality", categorical_colors)
  rep_palette <- custom_hierarchical_palette(metadata, "replicon_family", "replicon", categorical_colors)
  
  carba_palette <- metadata %>%
    mutate(carba_family = factor(
      ifelse(
        carba_allele == "-",
        "",
        ifelse(
          grepl(",", carba_allele),
          "Mixed",
          str_split_fixed(carba_allele, "-", 2)[,1])),
      levels = c("", "blaIMP", "blaKPC", "blaNDM", "blaOXA", "blaVIM", "Mixed"),
      ordered = TRUE, exclude = NULL)) %>%
    custom_hierarchical_palette("carba_family", "carba_allele", c("#FFFFFF", categorical_colors))
  
  list(
    numerical_palette = numerical_palette,
    cluster_palette = cluster_palette,
    species_palette = species_palette,
    ST_palette = ST_palette,
    geo_palette = geo_palette,
    rep_palette = rep_palette,
    carba_palette = carba_palette
  )
}


# Determine the order of Standard_Cluster based on frequency
cluster_order <- function() {
  names(sort(table(metadata$cluster), decreasing = TRUE))
}

# Functions: 
open_metadata <- function(path, source) {
  metadata_full <- read.csv(path, header = TRUE, sep = ";", fileEncoding = "ISO-8859-1")
  metadata_full$DataSource <- source
  metadata <- metadata_full <- metadata_full[!(metadata_full$cluster %in% c("-", "-1")), ]
  
  # Sort cluster order by frequency
  value_counts <- table(metadata$cluster)
  sorted_categories <- names(sort(value_counts, decreasing = TRUE))
  
  metadata <- metadata %>%
    mutate(cluster = ifelse(cluster %in% c("", "-"), "-", as.character(as.integer(cluster))),
           cluster = factor(cluster, levels = sorted_categories, ordered = TRUE),
           sampling_date = as.Date(sampling_date, "%d-%m-%Y"),
           tsne1D = ifelse(tsne1D %in% c("", "-"), NA, as.numeric(tsne1D)),
           tsne2D = ifelse(tsne2D %in% c("", "-"), NA, as.numeric(tsne2D)),
           Genus = sub(" .*$", "", Species),
           submitter_province = str_replace(submitter_province, "eilanden", "islands")
    ) %>%
    ## Update the Italic columns
    mutate(
      Species := factor(
        Species,
        labels = paste0("<i>", levels(factor(Species)), "</i>")
      ),
      Genus := factor(
        Genus,
        labels = paste0("<i>", levels(factor(Genus)), "</i>")
      ),
      carba_allele := ifelse(
        is.na(carba_allele) | carba_allele == "" | carba_allele == "-",
        "-",
        map_chr(
          str_split(as.character(carba_allele), "\\s*,\\s*"),
          ~ paste0("<i>", .x, "</i>", collapse = ",")
        )
      ),
      amr_genes := ifelse(
        is.na(amr_genes) | amr_genes == "" | amr_genes == "-",
        "-",
        map_chr(
          str_split(as.character(amr_genes), "\\s*,\\s*"),
          ~ paste0("<i>", .x, "</i>", collapse = ",")
        )
      ),
      metal_genes := ifelse(
        is.na(metal_genes) | metal_genes == "" | metal_genes == "-",
        "-",
        map_chr(
          str_split(as.character(metal_genes), "\\s*,\\s*"),
          ~ paste0("<i>", .x, "</i>", collapse = ",")
        )
      ),
      virulence_genes := ifelse(
        is.na(virulence_genes) | virulence_genes == "" | virulence_genes == "-",
        "-",
        map_chr(
          str_split(as.character(virulence_genes), "\\s*,\\s*"),
          ~ paste0("<i>", .x, "</i>", collapse = ",")
        )
      ),
      heat_genes := ifelse(
        is.na(heat_genes) | heat_genes == "" | heat_genes == "-",
        "-",
        map_chr(
          str_split(as.character(heat_genes), "\\s*,\\s*"),
          ~ paste0("<i>", .x, "</i>", collapse = ",")
        )
      ),
      biocide_genes := ifelse(
        is.na(biocide_genes) | biocide_genes == "" | biocide_genes == "-",
        "-",
        map_chr(
          str_split(as.character(biocide_genes), "\\s*,\\s*"),
          ~ paste0("<i>", .x, "</i>", collapse = ",")
        )
      ),
      acid_genes := ifelse(
        is.na(acid_genes) | acid_genes == "" | acid_genes == "-",
        "-",
        map_chr(
          str_split(as.character(acid_genes), "\\s*,\\s*"),
          ~ paste0("<i>", .x, "</i>", collapse = ",")
        )
      )
    ) %>%
    mutate(
      Species = case_when(
        Species == "<i>Enterobacter cloacae complex</i>" ~ "<i>Enterobacter cloacae</i> complex",
        TRUE ~ Species
      )
    )
  
  metadata$person_traveled_to <- factor(metadata$person_traveled_to, levels = c(setdiff(unique(metadata$person_traveled_to), "No known travel history"), "No known travel history"))

  # Generate a table from metadata$submitter_province
  parent_df <- metadata[c("Parent", "submitter_province", "submitter_municipality", "DataSource")]
  parent_df <- parent_df[!duplicated(parent_df$Parent),]
  
  return(list(metadata, parent_df, metadata_full))
}


create_normalised_co_occurance <- function(df, subset_val, breakdown_col, alt_meta) {
  if (!(alt_meta == "None")) {
    df <- df %>%
      mutate(cluster = case_when(
        cluster %in% subset_val ~ cluster,
        TRUE ~ as.character(!!sym(alt_meta))
      ))
  }
  if (!(breakdown_col == "None")) {
    df <- df %>%
      mutate(cluster = case_when(
        cluster %in% subset_val ~ paste(as.character(cluster), !! rlang::ensym(breakdown_col), sep="_"),
        TRUE ~ cluster
      ))
  }
  subset_co_occurance <- df %>%
    count(Parent, cluster) %>%
    pivot_wider(names_from = cluster, values_from=n, values_fill=0) %>%
    select(!"Parent")
  binary_df <- subset_co_occurance > 0
  # Compute the co-occurrence matrix
  co_occurrence_matrix <- t(binary_df) %*% binary_df
  
  # Calculate the number of times each object was brought
  subset_counts <- colSums(binary_df)
  
  # Create an empty matrix to store the normalized values
  normalized_matrix <- matrix(0, nrow = ncol(subset_co_occurance), ncol = ncol(subset_co_occurance))
  
  # Normalize the co-occurrence matrix
  for (i in 1:ncol(subset_co_occurance)) {
    for (j in 1:ncol(subset_co_occurance)) {
      normalized_matrix[i, j] <- co_occurrence_matrix[i, j] / subset_counts[i]
    }
  }
  rownames(normalized_matrix) <- colnames(subset_co_occurance)
  colnames(normalized_matrix) <- colnames(subset_co_occurance)
  # Convert the normalized matrix to a dataframe for plotting
  normalized_df <- as.data.frame(as.table(normalized_matrix))
  names(normalized_df) <- c("cluster1", "cluster2", "ratio")
  
  if (!(breakdown_col == "None")) {
    cluster_df <- normalized_df %>%
      filter(str_detect(cluster1, paste("^", as.character(subset_val), "_", sep=""))) %>%
      filter(!str_detect(cluster2,paste("^", as.character(subset_val), "_", sep="")))
  } else {
    cluster_df <- normalized_df %>%
      filter(cluster1 == subset_val) %>%
      filter(!(cluster2 == subset_val))
  }
  return(cluster_df)
}


## Palette functions
# Function to generate 'n' discrete colors along the gradient
discrete_palette <- function(n, color_list) {
  colors <- colorRampPalette(color_list)(n)  # Sample 'n' colors from the gradient
  return(colors)
}

# Custom discrete scale functions for ggplot2
scale_color_custom_discrete <- function(n, color_list) {
  discrete_scale("color", "custom", palette = function(n) discrete_palette(n, color_list))
}

scale_fill_custom_discrete <- function(n, color_list) {
  discrete_scale("fill", "custom", palette = function(n) discrete_palette(n, color_list))
}

custom_hierarchical_palette <- function(df, main_categories_name, subcategories_name, colors) {
  # df: input dataframe to get data from
  # main_categories_name: column name of top in hierarchy
  # subcategories_name: column name for subcategories in hierarchy
  # colors: base colors for main categories

  # Split subcategories by main category
  families <- split(df[, subcategories_name], df[, main_categories_name])
  subcategories <- lapply(families, unique)

  # Get main categories
  main_categories <- sort(unique(df[, main_categories_name]))

  # Generate main category palette
  main_palette <- colorRampPalette(colors)(length(main_categories))
  names(main_palette) <- main_categories

  # Generate subcategory colors
  category_colors <- c()
  for (i in seq_along(main_categories)) {
    base_color <- main_palette[i]
    subs <- subcategories[[main_categories[i]]]
    n_subs <- length(subs)
    
    if (n_subs > 1) {
      shade_vals <- seq(0, 0.5, length.out = n_subs)
      sub_palette <- sapply(shade_vals, function(x) lighten(base_color, x))
    } else {
      sub_palette <- base_color
    }
    
    category_colors <- c(category_colors, setNames(sub_palette, subs))
  }

  # Return a list with both main and subcategory palettes
  return(list(
    main_palette = main_palette,
    subcategory_palette = category_colors
  ))
}

update_palette <- function(df, main_categories_name, subcategories_name, colors, existing_palette = NULL) {
  # df: new dataframe
  # main_categories_name: column name of top in hierarchy
  # subcategories_name: column name for subcategories in hierarchy
  # colors: base colors for main categories
  # existing_palette: named vector of colors (output of custom_hierarchical_palette)
  
  # Get all subcategories in new dataframe
  new_subcategories <- unique(df[[subcategories_name]])
  
  if (is.null(existing_palette)) {
    return(custom_hierarchical_palette(df, main_categories_name, subcategories_name, colors))
  }
  
  existing_subcategories <- names(existing_palette)
  
  # Check if all new subcategories are already present
  missing_subs <- setdiff(new_subcategories, existing_subcategories)
  
  if (length(missing_subs) == 0) {
    return(existing_palette)
  } else {
    return(custom_hierarchical_palette(df, main_categories_name, subcategories_name, colors))
  }
}

## Determine and add palette
add_palette <- function(plot_in, column, df = metadata) {
  
  pal <- palette_by_stage[[column]]
  
  if (!is.null(pal)) {
    # Use the named palette
    plot_out <- plot_in + scale_fill_manual(values = pal)
  } else {
    # Use categorical fallback
    n_unique <- length(unique(df[[column]]))
    plot_out <- plot_in + scale_fill_custom_discrete(n = n_unique, color_list = categorical_colors)
  }
  
  plot_out
}

add_palette_colour <- function(plot_in, column, df = metadata) {
  
  pal <- palette_by_stage[[column]]
  
  if (!is.null(pal)) {
    # Use the named palette
    plot_out <- plot_in + scale_color_manual(values = pal)
  } else {
    # Use categorical fallback
    n_unique <- length(unique(df[[column]]))
    plot_out <- plot_in + scale_color_custom_discrete(n = n_unique, color_list = categorical_colors)
  }
  
  plot_out
}


## Plotting functions
geo_plot_ly <- function(df, geo_df, geo_level, title, name, cluster_in = "", fractionalise = FALSE) {
  if (fractionalise==TRUE) {
    sample_count_population <- df %>%
      select(Parent, !!sym(geo_level)) %>%
      distinct() %>%
      group_by(!!sym(geo_level)) %>%
      summarise(n_plasmids = n())
    number_display = "%.1f%%"
  } else {
    number_display = "%s"
  }
  # count samples per province
  if (cluster_in != "") {
    df <- df %>%
      filter(cluster %in% cluster_in)
    style_name = "points"
  } else {
    style_name = "fills"
  }
  sample_counts <- df %>%
    group_by(!!sym(geo_level)) %>%
    summarise(n_isolates = n())
  
  if (fractionalise==TRUE) {
    sample_counts <- sample_counts %>%
      merge(., sample_count_population, by.x = "submitter_municipality", by.y = "submitter_municipality") %>%
      mutate(n_isolates = ((n_isolates / n_plasmids) * 100))
    rescale_factor = 0.001
  } else {
    rescale_factor = 0.999
  }
  
  geo_counts <- left_join(geo_df, sample_counts, by = c("regio_naam" = geo_level)) %>%
    mutate(
      n_isolates = case_when(
        regio_soort == "rand" ~ NA_real_,   # always NA → transparent
        is.na(n_isolates)     ~ 0,          # missing counts → 0
        TRUE                  ~ n_isolates  # keep actual counts
      )
    )
  max_val <- max(geo_counts$n_isolates, na.rm = TRUE)

  gg_plot <- geo_counts %>%
    ggplot() +
    geom_sf(data = filter(geo_counts, regio_soort == "rand"),
            fill = NA, color = "black", linewidth=0.1) +
    geom_sf(data = filter(geo_counts, regio_soort != "rand"), 
            aes(fill=n_isolates,
                text=sprintf(str_glue("Location: %s<br>{name} (<i>n</i>): {number_display}"), regio_naam, n_isolates)
                ),
            linewidth=0.1) +
    theme_ggrivm() +
    scale_fill_gradientn(colours=numerical_palette,
                         values = rescale(c(0, rescale_factor, max_val * 0.5, max_val)),  # note: second value is just above 0
                         limits = c(0, max_val), # ensures 0 is on the scale
                         na.value="transparent",
    ) +
    #geom_text(
    #  data = subset(geo_df, regio_soort %in% c("rand")),
    #  aes(x = x, y = y-1010, label = regio_naam),
    #  size = 4,
    #  vjust = -0.5
    #) +
    ggtitle(title) +
    #labs(title="Number of samples per municipality") +
    theme(
      panel.background = element_blank(),
      panel.grid = element_blank(),
      legend.text=element_text(size=10, angle = 45),
      axis.title = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      axis.line = element_blank(),
      plot.title = element_text(margin = margin(b = -50)),
    ) +
    guides(fill = guide_colourbar(title.position="top", title.hjust = 0.5)) +
    labs(fill=expression("Samples ("~italic("N")~")"))+
    coord_sf(clip = "off")
  
  bb <- sf::st_bbox(geo_df)
  plot_ly <- ggplotly(gg_plot, tooltip = "text") %>%
    style(hoveron=style_name) %>%
    config(responsive = TRUE)
  rand_labels <- subset(geo_df, regio_soort == "rand")
  
  for(i in seq_len(nrow(rand_labels))){
    plot_ly <- plot_ly %>% add_annotations(
      x = rand_labels$x[i],
      y = rand_labels$y[i],
      text = rand_labels$regio_naam[i],
      showarrow = FALSE,
      xanchor = "center",
      yanchor = "bottom",
      textfont = list(size=6)
    )
  }
  plot_ly <- plot_ly %>%
    layout(
      xaxis = list(autorange = TRUE),
      yaxis = list(autorange = TRUE)
    )
  return(plot_ly)
}

categorical_bar <- function(df, column_name, colourlist, title, xlab) {
  df <- df %>%
    group_by(cluster) %>%
    mutate(HasUploaded = ifelse(any(DataSource == "UserUpload"), "Contains upload", "Reference only")) %>%
    ungroup()
  
  # Compute fraction per cluster and category
  df_fraction <- df %>%
    group_by(cluster, !!sym(column_name)) %>%
    summarise(Count = n(), .groups = "drop") %>%
    group_by(cluster) %>%
    mutate(Fraction = Count / sum(Count)) %>%
    ungroup()
  
  # Join fraction back to original df for hover
  df <- df %>%
    left_join(df_fraction, by = c("cluster", column_name))
  
  gg_plot <- df %>%
    #mutate(safe_category = ifelse(!!sym(column_name) %in% c("", "-"), NA, !!sym(column_name))) %>%
    ggplot(aes(x=cluster)) +
    geom_bar(aes(
      fill = !!sym(column_name),
      colour = HasUploaded,
      linewidth = HasUploaded,
      text=sprintf(
        str_glue("Cluster: %s<br>{column_name}: %s<br>Fraction: %.1f%%"),
        cluster,
        !!sym(column_name),
        Fraction * 100
      )
    ),
    position = "fill",
    #colour="#b4b4b4",
    #linewidth=ifelse(df$DataSource=="Reference", 0, 0.5),
    alpha=0.95,
    ) +
    theme_ggrivm() +
    labs(title = str_glue("{title} by Cluster"), y = "Ratio", x = xlab) +
    scale_fill_manual(values=colourlist, na.value="#FFFFFF") +
    scale_color_manual(
      values=c("Reference only" = "#b4b4b4", "Contains upload" = "grey17"),
      guide="none"
    ) +
    scale_linewidth_manual(values=c("Reference only" = 0.1, "Contains upload" = 1)) +
    theme(
      legend.text = element_text(size = 10),
      legend.key.size = unit(0.1, "cm"),
      axis.text.x = element_text(angle = 45, size= 10, hjust = 0.5, vjust = 0),
      axis.title.x = element_text(vjust=-1),
      plot.title = element_text(size= 12, color="black")
    ) +
    coord_cartesian(clip = 'off') 
  
  plot_ly <- ggplotly(gg_plot, tooltip="text") %>%
    layout(
      showlegend = TRUE,
      legend = list(
        orientation = "h",
        x = 0.5,
        xanchor = "center",
        y = -0.2,
        yanchor="top",
        width = 200,
        itemwidth = 1,
        itemheight = 1,
        scrollbar = TRUE ,
        fill=NA,
        title = list(text = "")
      )
    )
  
  df_legend <- data.frame(id = seq_along(plot_ly$x$data), legend_entries = unlist(lapply(plot_ly$x$data, `[[`, "name")))
  keep_groups <- unique(df[[column_name]])
  # Split each legend entry into components
  df_legend$components <- strsplit(df_legend$legend_entries, ",")
  
  # Clean parentheses and whitespace
  df_legend$components <- lapply(df_legend$components, function(x) gsub("^\\(|\\)$", "", trimws(x)))
  
  # Determine which legend entries to keep based on the column values
  df_legend$is_keep <- sapply(df_legend$components, function(x) any(x %in% keep_groups))
  
  # Determine the label to show in the legend (e.g., the intersection with keep_groups)
  df_legend$legend_group <- sapply(df_legend$components, function(x) paste(x[x %in% keep_groups], collapse = ","))
  df_legend$is_first <- !duplicated(df_legend$legend_group)
  
  # Update plotly traces
  for (i in df_legend$id) {
    group <- df_legend$legend_group[[i]]
    is_first <- df_legend$is_first[[i]]
    is_keep <- df_legend$is_keep[[i]]
    
    plot_ly$x$data[[i]]$name <- group
    plot_ly$x$data[[i]]$legendgroup <- group
    plot_ly$x$data[[i]]$showlegend <- is_first && is_keep
  }
  return(plot_ly)
}

count_bar <- function(df, column_name, title, xlab, none_string="-") {
  df <- df %>%
    group_by(cluster) %>%
    mutate(HasUploaded = ifelse(any(DataSource == "UserUpload"), "Contains upload", "Reference only")) %>%
    ungroup()
  
  
  # Compute fraction per cluster and category
  df_fraction <- df %>%
    group_by(cluster, !!sym(column_name)) %>%
    summarise(Count = n(), .groups = "drop") %>%
    group_by(cluster) %>%
    mutate(Fraction = Count / sum(Count)) %>%
    ungroup()
  
  # Join fraction back to original df for hover
  df <- df %>%
    left_join(df_fraction, by = c("cluster", column_name))
  
  gg_plot <- df %>%
    mutate(count = sapply(strsplit(!!sym(column_name) , ","), function(x) length(unique(x)))) %>%
    mutate(count = ifelse(!!sym(column_name) == none_string, 0, count)) %>%
    ggplot(aes(x = cluster)) +
    geom_bar(aes(
      fill = factor(count),
      colour = HasUploaded,
      linewidth = HasUploaded,
      text=sprintf(
        str_glue(
          "Cluster: %s, <br>ARGs (<i>n</i>): %s<br>Fraction: %.1f%%"
        ),
        cluster,
        count,
        Fraction * 100
      )
    ),
    position = "fill",
    #colour="#b4b4b4",
    alpha=0.9,
    ) +
    labs(title = str_glue("{title} by Cluster"), y = "Ratio", x = xlab) +
    scale_fill_custom_discrete(n = length(levels(count)), color_list=numerical_palette) +
    scale_color_manual(values=c("Reference only" = "#b4b4b4", "Contains upload" = "grey17")) +
    scale_linewidth_manual(values=c("Reference only" = 0.1, "Contains upload" = 1)) +
    theme_ggrivm() +
    theme(
      legend.position = "bottom",
      legend.text = element_text(size = 12),
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.key.size = unit(0.7, "cm"),
      plot.title = element_text(size= 12, color="black", face = "plain")
    )
  
  plot_ly <- ggplotly(gg_plot, tooltip="text") %>%
    layout(
      showlegend = TRUE,
      legend = list(
        orientation = "h",
        x = 0.5,
        xanchor = "center",
        y = -0.15,
        yanchor="top",
        width = 200,
        itemwidth = 30,
        itemheight = 20,
        scrollbar = TRUE ,
        title = list(text = "") 
      )
    ) 
  df_legend <- data.frame(id = seq_along(plot_ly$x$data), legend_entries = unlist(lapply(plot_ly$x$data, `[[`, "name")))
  keep_groups <- df %>%
    mutate(count = sapply(strsplit(!!sym(column_name), ","), function(x) length(unique(x)))) %>%
    mutate(count = ifelse(!!sym(column_name) == none_string, 0, count)) %>%
    pull(count) %>%
    unique()
  # Split each legend entry into components
  df_legend$components <- strsplit(df_legend$legend_entries, ",")
  
  # Clean parentheses and whitespace
  df_legend$components <- lapply(df_legend$components, function(x) gsub("^\\(|\\)$", "", trimws(x)))
  
  # Determine which legend entries to keep based on the column values
  df_legend$is_keep <- sapply(df_legend$components, function(x) any(x %in% keep_groups))
  
  # Determine the label to show in the legend (e.g., the intersection with keep_groups)
  df_legend$legend_group <- sapply(df_legend$components, function(x) paste(x[x %in% keep_groups], collapse = ","))
  df_legend$is_first <- !duplicated(df_legend$legend_group)
  
  # Update plotly traces
  for (i in df_legend$id) {
    group <- df_legend$legend_group[[i]]
    is_first <- df_legend$is_first[[i]]
    is_keep <- df_legend$is_keep[[i]]
    
    plot_ly$x$data[[i]]$name <- group
    plot_ly$x$data[[i]]$legendgroup <- group
    plot_ly$x$data[[i]]$showlegend <- is_first && is_keep
  }
  return(plot_ly)
}

frac_heatmap <- function(df, column_name_y, column_name_x, title, none_string = "-", cluster_filter = "") {
  # Generate a table from metadata$cluster
  category_table <- table(df$cluster)
  category_table_df <- as.data.frame(category_table)
  names(category_table_df) <- c("cluster", "total")
  
  if (cluster_filter != "") {
    df <- df %>%
      filter(cluster %in% cluster_filter)
    category_table_df <- category_table_df %>%
      filter(cluster %in% cluster_filter)
    x_val = column_name_y
    y_val = column_name_x
    plot_height = 300
  } else {
    x_val = column_name_x
    y_val = column_name_y
    y_categories <- df %>%
      separate_rows(!!sym(column_name_y), sep = ",") %>%
      filter(!!sym(column_name_y) != "") %>%
      pull(!!sym(column_name_y)) %>%
      unique()
    
    row_height <- 20
    plot_height <- max(300, length(y_categories) * row_height)
  }
  

  
  # Resistences genes per cluster
  # Create the heatmap using ggplot2
  gg_plot <- df %>%
    mutate(!!sym(column_name_y) := ifelse(!!sym(column_name_y) %in% c(""), "-", !!sym(column_name_y))) %>%
    separate_rows(!!sym(column_name_y), sep = ",") %>%
    distinct() %>%
    group_by(!!sym(column_name_y), !!sym(column_name_x)) %>%
    summarise(count = n()) %>%
    merge(., category_table_df, by = column_name_x) %>%
    mutate(ratio = (count / total),
           !!sym(column_name_y) := factor(!!sym(column_name_y), levels = sort(unique(!!sym(column_name_y)), decreasing=TRUE))
    ) %>%
    #complete(., !!sym(column_name_y), nesting(!!sym(column_name_x)), fill = list(ratio = 0, count = 0)) %>%
    ggplot(aes(
      x = !!sym(x_val),
      y = !!sym(y_val))) +
    geom_tile(
      aes(
        fill = ratio,
        text = sprintf(
          str_glue(
            "Cluster: %s, <br>Gene name: %s<br>Fraction: %.1f%%"
          ),
          !!sym(x_val),
          !!sym(y_val),
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
    labs(title = str_glue("{title} by Cluster"), y = y_val, x = x_val) +
    theme_ggrivm()  +
    theme(
      legend.position = "bottom",
      legend.text = element_text(size = 12),
      axis.text.x = element_text(angle = 45, hjust = 1),
      axis.text.y = element_text(angle = 45, hjust = 0, vjust = 0, size=8, margin = margin(r = 5)),
      legend.key.size = unit(0.7, "cm"),
      plot.title = element_text(size= 12, color="black", face = "plain")
    )
  plot_ly <- ggplotly(gg_plot, height = plot_height, tooltip="text") %>%
    layout(
      showlegend = FALSE
    )
  return(plot_ly)
}

ellips_scatter <- function(df, column_name="cluster") {
  df <- df %>%
    filter(cluster != "-") %>%
    mutate(tsne1D = ifelse(tsne1D %in% c("", "-"), NA, as.numeric(tsne1D))) %>%
    mutate(tsne2D = ifelse(tsne2D %in% c("", "-"), NA, as.numeric(tsne2D)))
  df_clusters <- df %>%
    filter(cluster != "-1")
  # mge tsne scatterplot
  gg_plot <- df %>%
    ggplot() +
    stat_ellipse(
      data = df_clusters,
      geom = "polygon",
      type = "norm",
      fill = NA,
      level = 0.999,
      lwd = 0.8,
      alpha = 0.8,
      inherit.aes = FALSE,
      show.legend = FALSE,
      aes(
        x = tsne1D,
        y = tsne2D,
        color = cluster
      ),
    ) +
    geom_point(
      aes(
        x = tsne1D,
        y = tsne2D,
        fill = !!sym(column_name),
        shape = DataSource,
      ),
      color = ifelse(df$DataSource=="Reference", "grey100", "black"),
      stroke = ifelse(df$DataSource=="Reference", 0, 0.5),
      alpha=0.65
    ) +
    scale_color_manual(
      values = cluster_palette,  # this will control the ellipse outlines
      #aesthetics = "color"      # specifically for the stat_ellipse layer
      guide="none"
    ) +
    scale_shape_manual(
      values = c("Reference" = 21, "UserUpload" = 24),
      guide = "none"
    ) +
    labs(title = str_glue("tSNE-coordinate scatterplot of clustered plasmids\nColoured by {column_name}"), y = "tsne2D", x = "tsne1D") +
    theme_ggrivm() +
    theme(
      legend.position = "bottom",
      legend.text = element_text(size = 12),
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.key.size = unit(0.7, "cm"),
      plot.title = element_text(size= 12, color="black", face = "plain")
    ) +
    guides(fill = guide_legend(title = column_name))
  
  gg_plot <- add_palette(gg_plot, column_name, df)

  plot_ly <- ggplotly(gg_plot) %>%
    layout(
      showlegend = TRUE
    )
  
  df_legend <- data.frame(
    id = seq_along(plot_ly$x$data),
    legend_entries = unlist(lapply(plot_ly$x$data, `[[`, "name")),
    stringsAsFactors = FALSE
  )
  keep_groups <- unique(df[[column_name]])
  # Split each legend entry into components
  df_legend$components <- strsplit(df_legend$legend_entries, ",")
  
  # Clean parentheses and whitespace
  df_legend$components <- lapply(df_legend$components, function(x) gsub("^\\(|\\)$", "", trimws(x)))
  
  if (column_name == "cluster") {
    df_legend$components <- lapply(df_legend$components, function(x) gsub("^1$|^, 1$", "", x[x != ""]))
  }
  
  # Determine which legend entries to keep based on the column values
  df_legend$is_keep <- sapply(df_legend$components, function(x) any(x %in% keep_groups))
  
  # Determine the label to show in the legend (e.g., the intersection with keep_groups)
  df_legend$legend_group <- sapply(df_legend$components, function(x) paste(x[x %in% keep_groups], collapse = ","))
  df_legend$is_first <- !duplicated(df_legend$legend_group)
  
  # Update plotly traces
  for (i in df_legend$id) {
    group <- df_legend$legend_group[[i]]
    is_first <- df_legend$is_first[[i]]
    is_keep <- df_legend$is_keep[[i]]
    
    plot_ly$x$data[[i]]$name <- group
    plot_ly$x$data[[i]]$legendgroup <- group
    plot_ly$x$data[[i]]$showlegend <- is_first && is_keep
  }
  return(plot_ly)
}

categorical_time_series <- function(df, cluster_filter, column) {
  
  time_plot <- df %>%
    filter(cluster %in% cluster_filter) %>%
    mutate(month = floor_date(as.Date(sampling_date), "month")) %>%
    count(month, !!sym(column), name = "count") %>%
    ggplot(
      aes(
        x = month,
        y = count,
        fill = !!sym(column),
        text = paste0(
          !!sym(column), ": ", count
        )
      )
    ) +
    geom_col(alpha = 0.8) +
    theme_ggrivm() +
    theme(
      legend.position = "bottom",
      legend.text = element_text(size = 12),
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.key.size = unit(0.7, "cm"),
      plot.title = element_text(size = 12, color = "black", face = "plain")
    )
  
  time_plot <- add_palette(time_plot, column, df)
  
  time_plotly <- ggplotly(time_plot, tooltip = "text") %>%
    layout(
      hovermode = "x unified",
      xaxis = list(
        type = "date",
        hoverformat = "%Y-%m",
        tickformat = "%Y-%m"
      )
    )
  return(time_plotly)
}

treemap <- function(df, column_name, grouping_column, display_cluster, palette) {
  count_distribution <- df %>%
    filter(cluster %in% display_cluster) %>%
    group_by(.data[[column_name]], .data[[grouping_column]]) %>%
    summarise(n = n(), .groups = "drop")
  
  child_nodes <- count_distribution %>%
    mutate(
      label  = .data[[column_name]],
      parent = .data[[grouping_column]],
      value  = n,
      color  = palette[as.character(.data[[column_name]])]   # map by name
    ) %>%
    select(label, parent, value, color)
  
  parent_nodes <- child_nodes %>%
    group_by(parent) %>%
    summarise(value = sum(value), .groups = "drop") %>%
    mutate(
      label  = parent,
      parent = "",
      color  = palette[as.character(label)]   # map by name
    ) %>%
    select(label, parent, value, color)
  
  treemap_data <- bind_rows(parent_nodes, child_nodes) %>%
    mutate(id = ifelse(parent == "", label, paste(parent, label, sep = "-")))
  
  tree_plot <- plot_ly(
    data = treemap_data,
    type = "treemap",
    ids = ~id,
    labels  = ~label,
    parents = ~parent,
    values  = ~value,
    branchvalues = "total",
    textinfo = "label+value+percent entry",
    hovertemplate = paste(
      "%{label}<br>",
      "<i>N</i>: %{value}<br>",
      "Percentage: %{percentEntry:.2%}<extra></extra>"
    ),
    marker = list(colors = ~color)
  )
  
  return(tree_plot)
}


get_node_color <- function(label, stage) {
  pal <- palette_by_stage[[stage]]
  
  # If a predefined palette exists for this stage
  if (!is.null(pal)) {
    col <- pal[as.character(label)]
    if (!is.na(col)) return(col)
  }
  
  # Fallback: unnamed categorical palette (per column)
  stage_levels <- unique(nodes$name[nodes$stage == stage])
  
  categorical_colors[
    (match(label, stage_levels) - 1) %% length(categorical_colors) + 1
  ]
}


##### Create necessary variables
# Load (main) data
here::i_am("global.R")
config <- config::get()
metadata_path <- config$METADATA_FILE
print(metadata_path)
dfs <- open_metadata(str_glue("{metadata_path}", sep=""), "Reference")
metadata <- dfs[[1]]
parent_df <- dfs[[2]]
metadata_full <- dfs[[3]]


## RIVM-associated palettes
blues = c("#154273", "#01689b", "#007bc7", "#8fcae7")
pinks = c("#42145f", "#a90061", "#f092cd")
reds = c("#ca005d", "#d52b1e")
greens = c("#275937", "#39870c", "#777b00", "#76d2b6")
yellows = c("#673327", "#94710a", "#e17000", "#ffb612", "#f9e11e")
greys = c("#f3f3f3", "#e6e6e6", "#cccccc", "#b4b4b4", "#999999", "#696969", "#535353")
categorical_colors = c("#007bc7", "#ffb612", "#ca005d","#39870c", "#76d2b6", "#552c6f", "#e17000",  "#673327") 
basics = c("#007bc7", "#ca005d", "#552c6f")

numerical_palette <- c(
  "#FFFFFF",
  palette_rivm("full")[1],
  palette_rivm("categorical")[2],
  palette_rivm("categorical")[3]
)

# Cluster
cluster_palette <- custom_hierarchical_palette(metadata, "cluster", "cluster", categorical_colors)$main_palette

# Species
origin_palettes <- custom_hierarchical_palette(metadata, "Genus", "Species", categorical_colors)
genus_palette <- origin_palettes$main_palette
species_palette <- origin_palettes$subcategory_palette
ST_palette <- custom_hierarchical_palette(metadata, "Genus", "ST", categorical_colors)$subcategory_palette

# Replicons
rep_palettes <- custom_hierarchical_palette(metadata, "replicon_family", "replicon", categorical_colors)
rep_family_palette <- rep_palettes$main_palette
rep_palette <- rep_palettes$subcategory_palette

# Geographic
geo_palettes <- custom_hierarchical_palette(metadata, "submitter_province", "submitter_municipality", categorical_colors)
province_palette <- geo_palettes$main_palette
geo_palette <- geo_palettes$subcategory_palette
travel_palette <- custom_hierarchical_palette(metadata, "person_traveled_to", "person_traveled_to", categorical_colors)$main_palette

# AMR
amr_palettes <- custom_hierarchical_palette(metadata, "amr_classes", "amr_genes", categorical_colors)
amr_gene_palette <- amr_palettes$subcategory_palette
amr_class_palette <- amr_palettes$main_palette
amr_palette <- custom_hierarchical_palette(metadata, "AMR_plasmid", "AMR_plasmid", c("azure4", "darkred"))$main_palette

# Single
mobility_palette <- custom_hierarchical_palette(metadata, "mobility", "mobility", categorical_colors[0:3])$main_palette
CP_palette <- custom_hierarchical_palette(metadata, "CP_plasmid", "CP_plasmid", c("azure4", "darkred"))$main_palette
metal_palette <- custom_hierarchical_palette(metadata, "metal_genes", "metal_genes", categorical_colors)$main_palette
virulence_palette <- custom_hierarchical_palette(metadata, "virulence_genes", "virulence_genes", categorical_colors)$main_palette

# Carba alleles
carba_palettes <- metadata %>%
  mutate(carba_family = factor(
    ifelse(
      carba_allele == "-",
      "",
      ifelse(
        grepl(",", carba_allele),
        "Mixed",
        str_split_fixed(carba_allele, "-", 2)[,1])),
    levels=c("", "<i>blaIMP", "<i>blaKPC", "<i>blaNDM", "<i>blaOXA", "<i>blaVIM", "Mixed"),
    ordered = TRUE, exclude = NULL)) %>%
  custom_hierarchical_palette(., "carba_family", "carba_allele", c("#FFFFFF", categorical_colors))
carba_palette <- carba_palettes$subcategory_palette


# Obtain the palettes
palette_by_stage <- list(
  Species = species_palette,
  Genus = genus_palette,
  ST = ST_palette,
  replicon_family = rep_family_palette,
  replicon = rep_palette,
  mobility = mobility_palette,
  carba_allele = carba_palette,
  CP_plasmid = CP_palette,
  metal_genes = metal_palette,
  virulence_genes = virulence_palette,
  amr_genes = amr_gene_palette,
  amr_classes = amr_class_palette,
  AMR_plasmid = amr_palette,
  submitter_province = province_palette,
  submitter_municipality = geo_palette,
  person_traveled_to = travel_palette,
  cluster = cluster_palette
)

# Alpha parameter for Sankey
add_alpha_vec <- function(hex, alpha = 0.4) {
  rgb_val <- t(col2rgb(hex))
  sprintf(
    "rgba(%d,%d,%d,%.2f)",
    rgb_val[,1],
    rgb_val[,2],
    rgb_val[,3],
    alpha
  )
}