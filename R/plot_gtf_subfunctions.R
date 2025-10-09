# Plots annotated features over region of interest
# @param gtf A nine column annotation file
# @param reg_start an integer specifying the start location of region of interest
# @param reg_stop an integer specifying the end location of region of interest
# @return gtf_plot

.plot_exons_only <- function(gtf, reg_start, reg_stop) {
  
  roi_length <- reg_stop - reg_start + 1
  
  # Extract transcript ID from attribute field
  # Assign all transcripts to gene_id for easy matching
  for (i in 1:nrow(gtf)) {
    # Gene and transcript ids should be present in all remaining observations
    tmp <- unlist(strsplit(gtf$attribute[i], ";"))
    gene_idx <- grep("gene_id", tmp)
    gene <- gsub("gene_id ", "", tmp[gene_idx])
    gtf$gene_id[i] <- gene
    
    transcript_idx <- grep(" transcript_id", tmp)
    transcript <- tmp[transcript_idx]
    gtf$transcript_id[i] <- gsub(" transcript_id ", "", transcript)
  }
  
  
  plot_df <- gtf %>%
    dplyr::mutate(
      transcript_idx = match(transcript_id, unique(transcript_id)),
      polygon_idx = dplyr::row_number()
    ) %>%
    dplyr::select(-c(score, frame, attribute))
  
  # Just in case, put rows in a good order for determining plot layout
  plot_df <- plot_df %>%
    #dplyr::arrange(strand, gene_id, transcript_id, start, end)
    dplyr::arrange(strand, transcript_idx, start, end)
  
  # Determine row order for plotting
  plot_df <- plot_df %>%
    dplyr::mutate(
      row_index = match(transcript_id, unique(transcript_id))
    )

  num_plot_rows <- max(plot_df$row_index)
  
  if (num_plot_rows > 20) {
    warning_msg <- glue::glue("Warning: This region contains more features than are plottable [{num_plot_rows}]. Selecting first 20.\n")
    cat(file = logfile, warning_msg, append = TRUE)
    plot_df <- plot_df %>%
      dplyr::filter(row_index <= 20)
  }
  
  # Add a column to track if a feature's positions were truncated on a pointed end
  # We're not adding the arrow if the feature extends past the region of interest
  plot_df <- plot_df %>%
    dplyr::mutate(truncated = dplyr::case_when(
      strand == "+" & end > reg_stop ~ TRUE,
      strand == "-" & start < reg_start ~ TRUE,
      .default = FALSE
    ))
  
  # calculate the x and y coordinates for the exons
  exons <- plot_df %>%
    dplyr::group_by(transcript_id) %>%
    dplyr::mutate(
      exon_idx = dplyr::row_number(),
      n_exons  = dplyr::n(),
      pointed  = dplyr::case_when(
        strand == "+" & exon_idx == n_exons & !truncated ~ TRUE,   # last exon for +
        strand == "-" & exon_idx == 1 & !truncated ~ TRUE,         # first exon for -
        strand == "." ~ FALSE,
        TRUE ~ FALSE
      ),
      next_start = dplyr::lead(start)
    ) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(
      start    = ifelse(start < reg_start, reg_start, start),
      end      = ifelse(end > reg_stop, reg_stop, end),
      midpoint = (end - start)/2 + start
    )
  
  # keep separate datasets
  exons_normal <- exons %>% dplyr::filter(!pointed)   # only normal exons
  exons_pointed <- exons %>% dplyr::filter(pointed)   # only the directional ones
  
  # Generate labels safely
  # One per line and offset above row
  exon_labels <- exons %>%
    dplyr::select(row_index, gene_id, transcript_id) %>%
    dplyr::distinct() %>%
    dplyr::mutate(
      x = reg_stop - ((reg_stop - reg_start) / 2),
      y = row_index - 0.4, # Will plot text above features when y is transformed in reverse
      label = glue::glue("{gene_id} - {transcript_id}")
    )
  
  # Generate coordinates for exon polygons
  exon_poly_data <- .get_exon_poly_data(exons_normal, exons_pointed, roi_length)
  
  TXT_SIZE <- 6
  
  # Generate Plot Object
  gtf_plot <- ggplot2::ggplot() +
    
    # Draw the line that connects exons
    ggplot2::geom_segment(
      data = exons %>% dplyr::filter(!is.na(next_start)),
      ggplot2::aes(
        x = end,
        xend = next_start,
        y = row_index,
        yend = row_index
      ),
      color = "black",
      linewidth = 0.3
    ) +
    
    # Exons
    ggplot2::geom_polygon(
      data = exon_poly_data,
      ggplot2::aes(
        x = x,
        y = y,
        group = polygon_idx,
        fill = strand
      ),
      color = "black"
    ) +
    
    # Exon Labels
    ggplot2::geom_text(
      data = exon_labels,
      ggplot2::aes(
        x = x,
        y = y,
        label = label
      ),
      hjust = 0.5,
      vjust = 0.5,
      size = TXT_SIZE
    ) +
    
    # Plot layout
    ggplot2::scale_y_continuous(
      breaks = NULL,
      labels = NULL,
      # The row indexes are currently in ascending order
      # but the intention is to plot row 1 on top
      # this transform of the y axis accomplishes that
      transform = scales::transform_reverse()
    ) +
    ggplot2::scale_x_continuous(
      limits = c(reg_start, reg_stop)
    ) +
    ggplot2::scale_fill_manual(
      values = c(
        "+" = "lightcoral",
        "-" = "skyblue",
        "." = "lightgray"
      )
    ) +
    ggplot2::labs(x = "Genomic coordinate") +
    ggplot2::theme_classic() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
      axis.title.y = ggplot2::element_blank()
    )
  
  return(gtf_plot)
}
  
# Plots the coverage over an interval
# @param df a dataframe
# @return a histogram plot
.plot_transcripts_only <- function(gtf, reg_start, reg_stop) {
  
  roi_length <- reg_stop - reg_start + 1
  
  # Extract Gene ID and transcript ID from attribute field
  for (i in 1:nrow(gtf)) {
    # Gene and transcript ids should be present in all remaining observations
    tmp <- unlist(strsplit(gtf$attribute[i], ";"))
    gene_idx <- grep("gene_id", tmp)
    gene <- gsub("gene_id ", "", tmp[gene_idx])
    gtf$gene_id[i] <- gene
    
    transcript_idx <- grep(" transcript_id", tmp)
    transcript <- tmp[transcript_idx]
    gtf$transcript_id[i] <- gsub(" transcript_id ", "", transcript)
  }
  
  plot_df <- gtf %>%
    dplyr::mutate(
      transcript_idx = match(transcript_id, unique(transcript_id)),
      polygon_idx = dplyr::row_number()
    ) %>%
    dplyr::select(-c(score, frame, attribute))
  
  # Arrange rows for plotting
  # Feature precedence: gene > exon
  plot_df <- plot_df %>%
    #dplyr::arrange(strand, gene_id, transcript_id, start, end)
    dplyr::arrange(strand, transcript_idx, start, end)
  
  # Determine row order for plotting
  plot_df <- plot_df %>%
    dplyr::mutate(row_index = dplyr::row_number())
  
  num_plot_rows <- max(plot_df$row_index)
  
  if (num_plot_rows > 20) {
    warning_msg <- glue::glue("Warning: This region contains more features than are plottable [{num_plot_rows}]. Selecting first 20.\n")
    cat(file = logfile, warning_msg, append = TRUE)
    plot_df <- plot_df %>%
      dplyr::filter(row_index <= 20)
  }
  
  # Add a column to track if a feature's positions were truncated on a pointed end
  # We're not adding the arrow if the feature extends past the region of interest
  # Also not adding arrows if no strand info is present, but we'll handle that elsewhere
  plot_df <- plot_df %>%
    dplyr::mutate(truncated = dplyr::case_when(
      strand == "+" & end > reg_stop ~ TRUE,
      strand == "-" & start < reg_start ~ TRUE,
      .default = FALSE
    ))
  
  # calculate the x and y coordinates for the genes
  transcripts <- plot_df %>%
    dplyr::filter(feature == "transcript") %>%
    dplyr::distinct(transcript_id, strand, .keep_all = TRUE) %>%
    dplyr::mutate(
      start = ifelse(start < reg_start, reg_start, start),
      end = ifelse(end > reg_stop, reg_stop, end),
      midpoint = (end - start)/2 + start,
    )
  
  transcript_labels <- transcripts %>%
    dplyr::select(row_index, gene_id, transcript_id) %>%
    dplyr::mutate(
      x = reg_stop - ((reg_stop - reg_start) / 2),
      y = row_index - 0.4,
      label = glue::glue("{gene_id} - {transcript_id} - (transcript)")
    )
  
  # transcripts -> make polygons with tip
  transcript_poly_data <- .get_pointed_poly_data(transcripts, roi_length)
  
  TXT_SIZE <- 6
  
  # Generate Plot Object
  gtf_plot <- ggplot2::ggplot() +
    
    # Transcripts
    ggplot2::geom_polygon(
      data = transcript_poly_data,
      ggplot2::aes(
        x = x,
        y = y,
        group = polygon_idx,
        fill = strand
      ),
      color = "black"
    ) +
    
    # Transcript labels
    ggplot2::geom_text(
      data = transcript_labels,
      ggplot2::aes(
        x = x,
        y = y,
        label = label
      ),
      hjust = 0.5,
      vjust = 0.5,
      size = TXT_SIZE
    ) +
    
    # Plot layout
    ggplot2::scale_y_continuous(
      breaks = NULL,
      labels = NULL,
      transform = scales::transform_reverse()
    ) +
    ggplot2::scale_x_continuous(
      limits = c(reg_start, reg_stop)
    ) +
    ggplot2::scale_fill_manual(
      values = c(
        "+" = "lightcoral",
        "-" = "skyblue",
        "." = "lightgray"
      )
    ) +
    ggplot2::labs(x = "Genomic coordinate") +
    ggplot2::theme_classic() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
      axis.title.y = ggplot2::element_blank()
    )
  
  return(gtf_plot)
}


# Plots the coverage over an interval
# @param df a dataframe
# @return a histogram plot
.plot_genes_only <- function(gtf, reg_start, reg_stop) {
  
  roi_length <- reg_stop - reg_start + 1
  
  # Extract Gene ID from attribute field
  for (i in 1:nrow(gtf)) {
    # Gene ids should be present in all remaining observations
    tmp <- unlist(strsplit(gtf$attribute[i], ";"))
    gene_idx <- grep("gene_id", tmp)
    gene <- gsub("gene_id ", "", tmp[gene_idx])
    gtf$gene_id[i] <- gene
  }
  
  plot_df <- gtf %>%
    dplyr::mutate(
      gene_idx = match(gene_id, unique(gene_id)),
      polygon_idx = dplyr::row_number()
    ) %>%
    dplyr::select(-c(score, frame, attribute))
  
  # Arrange rows for plotting
  # Feature precedence: gene > exon
  plot_df <- plot_df %>%
    dplyr::arrange(strand, gene_id, start, end)
  
  # Determine row order for plotting
  plot_df <- plot_df %>%
    dplyr::mutate(row_index = dplyr::row_number())
  
  num_plot_rows <- max(plot_df$row_index)
  
  if (num_plot_rows > 20) {
    warning_msg <- glue::glue("Warning: This region contains more features than are plottable [{num_plot_rows}]. Selecting first 20.\n")
    cat(file = logfile, warning_msg, append = TRUE)
    plot_df <- plot_df %>%
      dplyr::filter(row_index <= 20)
  }
  
  # Add a column to track if a feature's positions were truncated on a pointed end
  # We're not adding the arrow if the feature extends past the region of interest
  # Also not adding arrows if no strand info is present, but we'll handle that elsewhere
  plot_df <- plot_df %>%
    dplyr::mutate(truncated = dplyr::case_when(
      strand == "+" & end > reg_stop ~ TRUE,
      strand == "-" & start < reg_start ~ TRUE,
      .default = FALSE
    ))
  
  # calculate the x and y coordinates for the genes
  genes <- plot_df %>%
    dplyr::filter(feature == "gene") %>%
    dplyr::distinct(gene_id, strand, .keep_all = TRUE) %>%
    dplyr::mutate(
      start = ifelse(start < reg_start, reg_start, start),
      end = ifelse(end > reg_stop, reg_stop, end),
      midpoint = (end - start)/2 + start,
    )
  
  gene_labels <- genes %>%
    dplyr::select(row_index, gene_id) %>%
    dplyr::mutate(
      x = reg_stop - ((reg_stop - reg_start) / 2),
      y = row_index - 0.4,
      label = glue::glue("{gene_id} - (gene)")
    )
  
  # genes -> make polygons with tip
  gene_poly_data <- .get_pointed_poly_data(genes, roi_length)
  
  TXT_SIZE <- 6
  
  # Generate Plot Object
  gtf_plot <- ggplot2::ggplot() +

    # Genes
    ggplot2::geom_polygon(
      data = gene_poly_data,
      ggplot2::aes(
        x = x,
        y = y,
        group = polygon_idx,
        fill = strand
      ),
      color = "black"
    ) +
    
    # Gene labels
    ggplot2::geom_text(
      data = gene_labels,
      ggplot2::aes(
        x = x,
        y = y,
        label = label
      ),
      hjust = 0.5,
      vjust = 0.5,
      size = TXT_SIZE
    ) +
    
    # Plot layout
    ggplot2::scale_y_continuous(
      breaks = NULL,
      labels = NULL,
      transform = scales::transform_reverse()
    ) +
    ggplot2::scale_x_continuous(
      limits = c(reg_start, reg_stop)
    ) +
    ggplot2::scale_fill_manual(
      values = c(
        "+" = "lightcoral",
        "-" = "skyblue",
        "." = "lightgray"
      )
    ) +
    ggplot2::labs(x = "Genomic coordinate") +
    ggplot2::theme_classic() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
      axis.title.y = ggplot2::element_blank()
    )
  
  return(gtf_plot)
}

# Plots the coverage over an interval
# @param df a dataframe
# @return a histogram plot
.plot_transcripts_exons <- function(gtf, reg_start, reg_stop) {
  
  roi_length <- reg_stop - reg_start + 1
  
  # if transcript and exon information is both present
  # Extract transcript ID from attribute field
  # Assign all transcripts to gene_id for easy matching
  for (i in 1:nrow(gtf)) {
    # Gene and transcript ids should be present in all remaining observations
    tmp <- unlist(strsplit(gtf$attribute[i], ";"))
    gene_idx <- grep("gene_id", tmp)
    gene <- gsub("gene_id ", "", tmp[gene_idx])
    gtf$gene_id[i] <- gene
    
    transcript_idx <- grep(" transcript_id", tmp)
    transcript <- tmp[transcript_idx]
    gtf$transcript_id[i] <- gsub(" transcript_id ", "", transcript)
  }
  

  plot_df <- gtf %>%
    dplyr::mutate(
      transcript_idx = match(transcript_id, unique(transcript_id)),
      polygon_idx = dplyr::row_number()
    ) %>%
    dplyr::select(-c(score, frame, attribute))
  
  # Just in case, put rows in a good order for determining plot layout
  plot_df <- plot_df %>%
    #dplyr::arrange(strand, transcript_id, desc(feature), start, end)
    dplyr::arrange(strand, transcript_idx, desc(feature), start, end)

  # Determine row order for plotting
  plot_df <- plot_df %>%
    dplyr::mutate(
      first_exon = feature == "exon" & !duplicated(ifelse(feature == "exon", transcript_id, NA)),
      inc = (feature == "transcript") | first_exon,
      row_index = cumsum(inc)
    ) %>%
    dplyr::select(-first_exon, -inc)
  
  num_plot_rows <- max(plot_df$row_index)
  
  if (num_plot_rows > 20) {
    warning_msg <- glue::glue("Warning: This region contains more features than are plottable [{num_plot_rows}]. Selecting first 20.\n")
    cat(file = logfile, warning_msg, append = TRUE)
    plot_df <- plot_df %>%
      dplyr::filter(row_index <= 20)
  }
  
  # Add a column to track if a feature's positions were truncated on a pointed end
  # We're not adding the arrow if the feature extends past the region of interest
  plot_df <- plot_df %>%
    dplyr::mutate(truncated = dplyr::case_when(
      strand == "+" & end > reg_stop ~ TRUE,
      strand == "-" & start < reg_start ~ TRUE,
      .default = FALSE
    ))
    
  
  # calculate the x and y coordinates for the exons/transcripts/genes
  transcripts <- plot_df %>%
    dplyr::filter(!is.na(transcript_id) & feature == "transcript") %>%
    dplyr::distinct(transcript_id, gene_id, strand, .keep_all = TRUE) %>%
    dplyr::arrange(transcript_id) %>%
    dplyr::mutate(
      start = ifelse(start < reg_start, reg_start, start),
      end = ifelse(end > reg_stop, reg_stop, end),
      midpoint = (end - start)/2 + start
    )

  exons <- plot_df %>%
    dplyr::filter(feature == "exon") %>%
    dplyr::group_by(transcript_id) %>%
    dplyr::mutate(
      exon_idx = dplyr::row_number(),
      n_exons  = dplyr::n(),
      pointed  = dplyr::case_when(
        strand == "+" & exon_idx == n_exons & !truncated ~ TRUE,   # last exon for +
        strand == "-" & exon_idx == 1 & !truncated ~ TRUE,         # first exon for -
        strand == "." ~ FALSE,
        TRUE ~ FALSE
      ),
      next_start = dplyr::lead(start)
    ) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(
      start    = ifelse(start < reg_start, reg_start, start),
      end      = ifelse(end > reg_stop, reg_stop, end),
      midpoint = (end - start)/2 + start
    )

  # keep separate datasets
  exons_normal <- exons %>% dplyr::filter(!pointed)   # only normal exons
  exons_pointed <- exons %>% dplyr::filter(pointed)   # only the directional ones

  # Generate labels safely
  # One per line and offset above row
  exon_labels <- exons %>%
    dplyr::select(row_index, gene_id, transcript_id) %>%
    dplyr::distinct() %>%
    dplyr::mutate(
      x = reg_stop - ((reg_stop - reg_start) / 2),
      y = row_index - 0.4, # Will plot text above features when y is transformed in reverse
      label = glue::glue("{gene_id} - {transcript_id}")
    )
  
  transcript_labels <- transcripts %>%
    dplyr::select(row_index, gene_id, transcript_id) %>%
    dplyr::mutate(
      x = reg_stop - ((reg_stop - reg_start) / 2),
      y = row_index - 0.4,
      label = glue::glue("{gene_id} - {transcript_id} - (transcript)")
    )
  
  # Generate polygon coordinates for exons
  exon_poly_data <- .get_exon_poly_data(exons_normal, exons_pointed, roi_length)
  
  # transcripts -> make polygons with tip like transcripts
  transcript_poly_data <- .get_pointed_poly_data(transcripts, roi_length)
  
  TXT_SIZE <- 6
  
  # Generate Plot Object
  gtf_plot <- ggplot2::ggplot() +
    
    # Draw the line that connects exons
    ggplot2::geom_segment(
      data = exons %>% dplyr::filter(!is.na(next_start)),
      ggplot2::aes(
        x = end,
        xend = next_start,
        y = row_index,
        yend = row_index
      ),
      color = "black",
      linewidth = 0.3
    ) +
    
    # Exons
    ggplot2::geom_polygon(
      data = exon_poly_data,
      ggplot2::aes(
        x = x,
        y = y,
        group = polygon_idx,
        fill = strand
      ),
      color = "black"
    ) +
    
    # Exon Labels
    ggplot2::geom_text(
      data = exon_labels,
      ggplot2::aes(
        x = x,
        y = y,
        label = label
      ),
      hjust = 0.5,
      vjust = 0.5,
      size = TXT_SIZE
    ) +
    
    # Transcripts
    ggplot2::geom_polygon(
      data = transcript_poly_data,
      ggplot2::aes(
        x = x,
        y = y,
        group = polygon_idx,
        fill = strand
      ),
      color = "black"
    ) +
    
    # Transcript labels
    ggplot2::geom_text(
      data = transcript_labels,
      ggplot2::aes(
        x = x,
        y = y,
        label = label
      ),
      hjust = 0.5,
      vjust = 0.5,
      size = TXT_SIZE
    ) +
    
    # Plot layout
    ggplot2::scale_y_continuous(
      breaks = NULL,
      labels = NULL,
      # The row indexes are currently in ascending order
      # but the intention is to plot row 1 on top
      # this transform of the y axis accomplishes that
      transform = scales::transform_reverse()
    ) +
    ggplot2::scale_x_continuous(
      limits = c(reg_start, reg_stop)
    ) +
    ggplot2::scale_fill_manual(
      values = c(
        "+" = "lightcoral",
        "-" = "skyblue",
        "." = "lightgray"
      )
    ) +
    ggplot2::labs(x = "Genomic coordinate") +
    ggplot2::theme_classic() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
      axis.title.y = ggplot2::element_blank()
    )
  
  return(gtf_plot)
}

# Plots the coverage over an interval
# @param df a dataframe
# @return a histogram plot
.plot_genes_exons <- function(gtf, reg_start, reg_stop) {
  
  roi_length <- reg_stop - reg_start + 1

  # Extract Gene ID from attribute field
  # Assign gene_id to gene and exon observations for grouping
  for (i in 1:nrow(gtf)) {
    # Gene ids should be present in all remaining observations
    tmp <- unlist(strsplit(gtf$attribute[i], ";"))
    gene_idx <- grep("gene_id", tmp)
    gene <- gsub("gene_id ", "", tmp[gene_idx])
    gtf$gene_id[i] <- gene
    
    if (gtf$feature[i] == "exon") {
      transcript_idx <- grep(" transcript_id", tmp)
      transcript <- tmp[transcript_idx]
      gtf$transcript_id[i] <- gsub(" transcript_id ", "", transcript) 
    } else {
      gtf$transcript_id[i] <- NA
    }
    
  }
  
  plot_df <- gtf %>%
    dplyr::mutate(
      gene_idx = match(gene_id, unique(gene_id)),
      transcript_idx = match(transcript_id, unique(transcript_id)),
      polygon_idx = dplyr::row_number()
    ) %>%
    dplyr::select(-c(score, frame, attribute))
  
  
  # Arrange rows for plotting
  # Feature precedence: gene > exon
  plot_df <- plot_df %>%
    #dplyr::arrange(strand, gene_id, desc(feature), transcript_id, start, end)
    dplyr::arrange(strand, gene_id, desc(feature), transcript_idx, start, end)
  
  # Determine row order for plotting
  plot_df <- plot_df %>%
    dplyr::mutate(
      first_exon = feature == "exon" & !duplicated(ifelse(feature == "exon", transcript_id, NA)),
      inc = (feature == "gene") | first_exon,
      row_index = cumsum(inc)
    ) %>%
    dplyr::select(-first_exon, -inc)
  
  num_plot_rows <- max(plot_df$row_index)
  
  if (num_plot_rows > 20) {
    warning_msg <- glue::glue("Warning: This region contains more features than are plottable [{num_plot_rows}]. Selecting first 20.\n")
    cat(file = logfile, warning_msg, append = TRUE)
    plot_df <- plot_df %>%
      dplyr::filter(row_index <= 20)
  }
  
  # Add a column to track if a feature's positions were truncated on a pointed end
  # We're not adding the arrow if the feature extends past the region of interest
  # Also not adding arrows if no strand info is present, but we'll handle that elsewhere
  plot_df <- plot_df %>%
    dplyr::mutate(truncated = dplyr::case_when(
      strand == "+" & end > reg_stop ~ TRUE,
      strand == "-" & start < reg_start ~ TRUE,
      .default = FALSE
    ))
  
  # calculate the x and y coordinates for the exons/transcripts/genes
  genes <- plot_df %>%
    dplyr::filter(feature == "gene") %>%
    dplyr::distinct(gene_id, strand, .keep_all = TRUE) %>%
    #dplyr::arrange(gene_id) %>%
    dplyr::mutate(
      start = ifelse(start < reg_start, reg_start, start),
      end = ifelse(end > reg_stop, reg_stop, end),
      midpoint = (end - start)/2 + start,
    )

  exons <- plot_df %>%
    dplyr::filter(feature == "exon") %>%
    dplyr::group_by(transcript_id) %>%
    dplyr::mutate(
      exon_idx= dplyr::row_number(),
      n_exons = dplyr::n(),
      pointed = dplyr::case_when(
        strand == "+" & exon_idx == n_exons & !truncated ~ TRUE,   # last exon for +
        strand == "-" & exon_idx == 1 & !truncated ~ TRUE,   # first exon for -
        strand == "." ~ FALSE, # No strand given
        TRUE ~ FALSE
      ),
      next_start = dplyr::lead(start)
    ) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(
      start = ifelse(start < reg_start, reg_start, start),
      end = ifelse(end > reg_stop, reg_stop, end),
      midpoint = (end - start)/2 + start,
    )
  
  # keep separate datasets
  exons_normal <- exons %>% dplyr::filter(!pointed)   # only normal exons
  exons_pointed <- exons %>% dplyr::filter(pointed)   # only the directional ones
  
  # Generate labels safely
  # One per line and offset above row
  exon_labels <- exons %>%
    dplyr::select(row_index, gene_id, transcript_id) %>%
    dplyr::distinct() %>%
    dplyr::mutate(
      x = reg_stop - ((reg_stop - reg_start) / 2),
      y = row_index - 0.4, # Will plot text above features when y is transformed in reverse
      label = glue::glue("{gene_id} - {transcript_id}")
    )
  
  gene_labels <- genes %>%
    dplyr::select(row_index, gene_id) %>%
    dplyr::mutate(
      x = reg_stop - ((reg_stop - reg_start) / 2),
      y = row_index - 0.4,
      label = glue::glue("{gene_id} - (gene)")
    )

  # Generate polygon coordinates for exons
  exon_poly_data <- .get_exon_poly_data(exons_normal, exons_pointed, roi_length)
  
  # genes -> make polygons with tip like transcripts
  gene_poly_data <- .get_pointed_poly_data(genes, roi_length)
  
  TXT_SIZE <- 6
  
  # Generate Plot Object
  gtf_plot <- ggplot2::ggplot() +
    
    # Draw the line that connects exons
    ggplot2::geom_segment(
      data = exons %>% dplyr::filter(!is.na(next_start)),
      ggplot2::aes(
        x = end,
        xend = next_start,
        y = row_index,
        yend = row_index
      ),
      color = "black",
      linewidth = 0.3
    ) +
    
    # Exons
    ggplot2::geom_polygon(
      data = exon_poly_data,
      ggplot2::aes(
        x = x,
        y = y,
        group = polygon_idx,
        fill = strand
      ),
      color = "black"
    ) +
    
    # Exon Labels
    ggplot2::geom_text(
      data = exon_labels,
      ggplot2::aes(
        x = x,
        y = y,
        label = label
      ),
      hjust = 0.5,
      vjust = 0.5,
      size = TXT_SIZE
    ) +
  
    # Genes
    ggplot2::geom_polygon(
      data = gene_poly_data,
      ggplot2::aes(
        x = x,
        y = y,
        group = polygon_idx,
        fill = strand
      ),
      color = "black"
    ) +
    
    # Gene labels
    ggplot2::geom_text(
      data = gene_labels,
      ggplot2::aes(
        x = x,
        y = y,
        label = label
      ),
      hjust = 0.5,
      vjust = 0.5,
      size = TXT_SIZE
    ) +
    
    # Plot layout
    ggplot2::scale_y_continuous(
      breaks = NULL,
      labels = NULL,
      transform = scales::transform_reverse()
    ) +
    ggplot2::scale_x_continuous(
      limits = c(reg_start, reg_stop)
    ) +
    ggplot2::scale_fill_manual(
      values = c(
        "+" = "lightcoral",
        "-" = "skyblue",
        "." = "lightgray"
      )
    ) +
    ggplot2::labs(x = "Genomic coordinate") +
    ggplot2::theme_classic() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
      axis.title.y = ggplot2::element_blank()
    )
  
  return(gtf_plot)
}


.plot_genes_transcripts <- function(gtf, reg_start, reg_stop) {
  
  roi_length <- reg_stop - reg_start + 1
  
  # Extract Gene ID and Transcript Id from attribute field
  # Assign gene_id to all observations
  # Assign transcript_id to transcripts and exons for grouping
  for (i in 1:nrow(gtf)) {
    # Gene ids should be present in all remaining observations
    tmp <- unlist(strsplit(gtf$attribute[i], ";"))
    gene_idx <- grep("gene_id", tmp)
    gene <- gsub("gene_id ", "", tmp[gene_idx])
    gtf$gene_id[i] <- gene
    
    if (gtf$feature[i] == "exon" | gtf$feature[i] == "transcript") {
      transcript_idx <- grep(" transcript_id", tmp)
      transcript <- tmp[transcript_idx]
      gtf$transcript_id[i] <- gsub(" transcript_id ", "", transcript) 
    } else {
      gtf$transcript_id[i] <- NA
    }
  }
  
  plot_df <- gtf %>%
    dplyr::mutate(
      gene_idx = match(gene_id, unique(gene_id)),
      transcript_idx = match(transcript_id, unique(transcript_id)),
      polygon_idx = dplyr::row_number()
    ) %>%
    dplyr::select(-c(score, frame, attribute))
  
  # Arrange rows for plotting
  # Feature precedence: gene > transcript
  plot_df <- plot_df %>%
    dplyr::arrange(strand, transcript_idx)
  
  # Determine row order for plotting
  plot_df <- plot_df %>%
    dplyr::mutate(row_index = dplyr::row_number())
  
  num_plot_rows <- max(plot_df$row_index)
  
  if (num_plot_rows > 20) {
    warning_msg <- glue::glue("Warning: This region contains more features than are plottable [{num_plot_rows}]. Selecting first 20.\n")
    cat(file = logfile, warning_msg, append = TRUE)
    plot_df <- plot_df %>%
      dplyr::filter(row_index <= 20)
  }
  
  # Add a column to track if a feature's positions were truncated on a pointed end
  # We're not adding the arrow if the feature extends past the region of interest
  # Also not adding arrows if no strand info is present, but we'll handle that elsewhere
  plot_df <- plot_df %>%
    dplyr::mutate(truncated = dplyr::case_when(
      strand == "+" & end > reg_stop ~ TRUE,
      strand == "-" & start < reg_start ~ TRUE,
      .default = FALSE
    ))
  
  # calculate the x and y coordinates for the exons/transcripts/genes
  genes <- plot_df %>%
    dplyr::filter(feature == "gene") %>%
    dplyr::distinct(gene_id, strand, .keep_all = TRUE) %>%
    dplyr::mutate(
      start = ifelse(start < reg_start, reg_start, start),
      end = ifelse(end > reg_stop, reg_stop, end),
      midpoint = (end - start)/2 + start
    )
  
  transcripts <- plot_df %>%
    dplyr::filter(feature == "transcript") %>%
    dplyr::distinct(transcript_id, strand, .keep_all = TRUE) %>%
    dplyr::mutate(
      start = ifelse(start < reg_start, reg_start, start),
      end = ifelse(end > reg_stop, reg_stop, end),
      midpoint = (end - start)/2 + start
    )
  
  gene_labels <- genes %>%
    dplyr::select(row_index, gene_id) %>%
    dplyr::mutate(
      x = reg_stop - ((reg_stop - reg_start) / 2),
      y = row_index - 0.4,
      label = glue::glue("{gene_id} - (gene)")
    )
  
  transcript_labels <- transcripts %>%
    dplyr::select(row_index, gene_id, transcript_id) %>%
    dplyr::mutate(
      x = reg_stop - ((reg_stop - reg_start) / 2),
      y = row_index - 0.4,
      label = glue::glue("{gene_id} - {transcript_id} - (transcript)")
    )
  
  # genes -> make polygons with tip like transcripts
  gene_poly_data <- .get_pointed_poly_data(genes, roi_length)
  
  # transcripts -> make polygons with tip like transcripts
  transcript_poly_data <- .get_pointed_poly_data(transcripts, roi_length)
  
  TXT_SIZE <- 6
  
  # Generate Plot Object
  gtf_plot <- ggplot2::ggplot() +
    
    # Genes
    ggplot2::geom_polygon(
      data = gene_poly_data,
      ggplot2::aes(
        x = x,
        y = y,
        group = polygon_idx,
        fill = strand
      ),
      color = "black"
    ) +
    
    # Gene labels
    ggplot2::geom_text(
      data = gene_labels,
      ggplot2::aes(
        x = x,
        y = y,
        label = label
      ),
      hjust = 0.5,
      vjust = 0.5,
      size = TXT_SIZE
    ) +
    
    # Transcripts
    ggplot2::geom_polygon(
      data = transcript_poly_data,
      ggplot2::aes(
        x = x,
        y = y,
        group = polygon_idx,
        fill = strand
      ),
      color = "black"
    ) +
    
    # Transcript labels
    ggplot2::geom_text(
      data = transcript_labels,
      ggplot2::aes(
        x = x,
        y = y,
        label = label
      ),
      hjust = 0.5,
      vjust = 0.5,
      size = TXT_SIZE
    ) +
    
    # Plot layout
    ggplot2::scale_y_continuous(
      breaks = NULL,
      labels = NULL,
      transform = scales::transform_reverse()
    ) +
    ggplot2::scale_x_continuous(
      limits = c(reg_start, reg_stop)
    ) +
    ggplot2::scale_fill_manual(
      values = c(
        "+" = "lightcoral",
        "-" = "skyblue",
        "." = "lightgray"
      )
    ) +
    ggplot2::labs(x = "Genomic coordinate") +
    ggplot2::theme_classic() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
      axis.title.y = ggplot2::element_blank()
    )
  
  return(gtf_plot)
}

# Plots the coverage over an interval
# @param df a dataframe
# @return a histogram plot
.plot_genes_exon_transcripts <- function(gtf, reg_start, reg_stop) {
  
  roi_length <- reg_stop - reg_start + 1
  
  # Extract Gene ID and Transcript Id from attribute field
  # Assign gene_id to all observations
  # Assign transcript_id to transcripts and exons for grouping
  for (i in 1:nrow(gtf)) {
    # Gene ids should be present in all remaining observations
    tmp <- unlist(strsplit(gtf$attribute[i], ";"))
    gene_idx <- grep("gene_id", tmp)
    gene <- gsub("gene_id ", "", tmp[gene_idx])
    gtf$gene_id[i] <- gene
    
    if (gtf$feature[i] == "exon" | gtf$feature[i] == "transcript") {
      transcript_idx <- grep(" transcript_id", tmp)
      transcript <- tmp[transcript_idx]
      gtf$transcript_id[i] <- gsub(" transcript_id ", "", transcript) 
    } else {
      gtf$transcript_id[i] <- NA
    }
  }
  
  plot_df <- gtf %>%
    dplyr::mutate(
      gene_idx = match(gene_id, unique(gene_id)),
      transcript_idx = match(transcript_id, unique(transcript_id)),
      polygon_idx = dplyr::row_number()
    ) %>%
    dplyr::select(-c(score, frame, attribute))
    
  # Modify "gene" features to "xgene" for temporary sorting
  plot_df$feature[plot_df$feature == "gene"] <- "xgene"
    
  # Arrange rows for plotting
  # Feature precedence: gene > transcript > exon
  plot_df <- plot_df %>%
    dplyr::arrange(strand, transcript_idx, desc(feature), start, end)
  
  # Change xgene back to gene
  plot_df$feature[plot_df$feature == "xgene"] <- "gene"
  
  # Determine row order for plotting
  plot_df <- plot_df %>%
    dplyr::mutate(
      first_exon = feature == "exon" & !duplicated(ifelse(feature == "exon", transcript_id, NA)),
      inc = (feature == "gene") | (feature == "transcript") | first_exon,
      row_index = cumsum(inc)
    ) %>%
    dplyr::select(-first_exon, -inc)
  
  num_plot_rows <- max(plot_df$row_index)
  
  if (num_plot_rows > 20) {
    warning_msg <- glue::glue("Warning: This region contains more features than are plottable [{num_plot_rows}]. Selecting first 20.\n")
    cat(file = logfile, warning_msg, append = TRUE)
    plot_df <- plot_df %>%
      dplyr::filter(row_index <= 20)
  }

  # Add a column to track if a feature's positions were truncated on a pointed end
  # We're not adding the arrow if the feature extends past the region of interest
  # Also not adding arrows if no strand info is present, but we'll handle that elsewhere
  plot_df <- plot_df %>%
    dplyr::mutate(truncated = dplyr::case_when(
      strand == "+" & end > reg_stop ~ TRUE,
      strand == "-" & start < reg_start ~ TRUE,
      .default = FALSE
    ))
  
  # calculate the x and y coordinates for the exons/transcripts/genes
  genes <- plot_df %>%
    dplyr::filter(feature == "gene") %>%
    dplyr::distinct(gene_id, strand, .keep_all = TRUE) %>%
    dplyr::mutate(
      start = ifelse(start < reg_start, reg_start, start),
      end = ifelse(end > reg_stop, reg_stop, end),
      midpoint = (end - start)/2 + start
    )
  
  transcripts <- plot_df %>%
    dplyr::filter(feature == "transcript") %>%
    dplyr::distinct(transcript_id, strand, .keep_all = TRUE) %>%
    dplyr::mutate(
      start = ifelse(start < reg_start, reg_start, start),
      end = ifelse(end > reg_stop, reg_stop, end),
      midpoint = (end - start)/2 + start
    )
  
  exons <- plot_df %>%
    dplyr::filter(feature == "exon") %>%
    dplyr::group_by(transcript_id) %>%
    dplyr::mutate(
      exon_idx= dplyr::row_number(),
      n_exons = dplyr::n(),
      pointed = dplyr::case_when(
        strand == "+" & exon_idx == n_exons & !truncated ~ TRUE,   # last exon for +
        strand == "-" & exon_idx == 1 & !truncated ~ TRUE,   # first exon for -
        strand == "." ~ FALSE, # No strand given
        TRUE ~ FALSE
      ),
      next_start = dplyr::lead(start)
    ) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(
      start = ifelse(start < reg_start, reg_start, start),
      end = ifelse(end > reg_stop, reg_stop, end),
      midpoint = (end - start)/2 + start,
    )
  
  # keep separate datasets
  exons_normal <- exons %>% dplyr::filter(!pointed)   # only normal exons
  exons_pointed <- exons %>% dplyr::filter(pointed)   # only the directional ones
  
  # Generate labels safely
  # One per line and offset above row
  exon_labels <- exons %>%
    dplyr::select(row_index, gene_id, transcript_id) %>%
    dplyr::distinct() %>%
    dplyr::mutate(
      x = reg_stop - ((reg_stop - reg_start) / 2),
      y = row_index - 0.4, # Will plot text above features when y is transformed in reverse
      label = glue::glue("{gene_id} - {transcript_id}")
    )
  
  gene_labels <- genes %>%
    dplyr::select(row_index, gene_id) %>%
    dplyr::mutate(
      x = reg_stop - ((reg_stop - reg_start) / 2),
      y = row_index - 0.4,
      label = glue::glue("{gene_id} - (gene)")
    )
  
  transcript_labels <- transcripts %>%
    dplyr::select(row_index, gene_id, transcript_id) %>%
    dplyr::mutate(
      x = reg_stop - ((reg_stop - reg_start) / 2),
      y = row_index - 0.4,
      label = glue::glue("{gene_id} - {transcript_id} - (transcript)")
    )
  
  # Generate polygon coordinates for exons
  exon_poly_data <- .get_exon_poly_data(exons_normal, exons_pointed, roi_length)
  
  # genes -> make polygons with tip like transcripts
  gene_poly_data <- .get_pointed_poly_data(genes, roi_length)
  
  # transcripts -> make polygons with tip like transcripts
  transcript_poly_data <- .get_pointed_poly_data(transcripts, roi_length)
  
  TXT_SIZE <- 6
  
  # Generate Plot Object
  gtf_plot <- ggplot2::ggplot() +
    
    # Draw the line that connects exons
    ggplot2::geom_segment(
      data = exons %>% dplyr::filter(!is.na(next_start)),
      ggplot2::aes(
        x = end,
        xend = next_start,
        y = row_index,
        yend = row_index
      ),
      color = "black",
      linewidth = 0.3
    ) +
    
    # Exons
    ggplot2::geom_polygon(
      data = exon_poly_data,
      ggplot2::aes(
        x = x,
        y = y,
        group = polygon_idx,
        fill = strand
      ),
      color = "black"
    ) +
    
    # Exon Labels
    ggplot2::geom_text(
      data = exon_labels,
      ggplot2::aes(
        x = x,
        y = y,
        label = label
      ),
      hjust = 0.5,
      vjust = 0.5,
      size = TXT_SIZE
    ) +
    
    # Genes
    ggplot2::geom_polygon(
      data = gene_poly_data,
      ggplot2::aes(
        x = x,
        y = y,
        group = polygon_idx,
        fill = strand
      ),
      color = "black"
    ) +
    
    # Gene labels
    ggplot2::geom_text(
      data = gene_labels,
      ggplot2::aes(
        x = x,
        y = y,
        label = label
      ),
      hjust = 0.5,
      vjust = 0.5,
      size = TXT_SIZE
    ) +
    
    # Transcripts
    ggplot2::geom_polygon(
      data = transcript_poly_data,
      ggplot2::aes(
        x = x,
        y = y,
        group = polygon_idx,
        fill = strand
      ),
      color = "black"
    ) +
    
    # Transcript labels
    ggplot2::geom_text(
      data = transcript_labels,
      ggplot2::aes(
        x = x,
        y = y,
        label = label
      ),
      hjust = 0.5,
      vjust = 0.5,
      size = TXT_SIZE
    ) +
    
    # Plot layout
    ggplot2::scale_y_continuous(
      breaks = NULL,
      labels = NULL,
      transform = scales::transform_reverse()
    ) +
    ggplot2::scale_x_continuous(
      limits = c(reg_start, reg_stop)
    ) +
    ggplot2::scale_fill_manual(
      values = c(
        "+" = "lightcoral",
        "-" = "skyblue",
        "." = "lightgray"
      )
    ) +
    ggplot2::labs(x = "Genomic coordinate") +
    ggplot2::theme_classic() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
      axis.title.y = ggplot2::element_blank()
    )
  
  return(gtf_plot)
}

# Generate exon polygon coordinates for pointed and non-pointed exons
.get_exon_poly_data <- function(normal_exons, pointed_exons, roi_length) {
  # Empty data frame for safe row binding
  empty_poly_data <- data.frame(
    seqname = character(0),
    source = character(0),
    feature = character(0),
    start = integer(0),
    end = integer(0),
    strand = character(0),
    gene_id = character(0),
    transcript_id = character(0),
    transcript_idx = integer(0),
    row_index = integer(0),
    truncated = logical(0),
    exon_idx = integer(0),
    n_exons = integer(0),
    pointed = logical(0),
    next_start = integer(0),
    midpoint = numeric(0),
    x = integer(0),
    y = numeric(0)      
  )
  
  # In certain cases, there may be no pointed exons
  # - If the only exons that would be pointed are truncated, for example
  if (nrow(pointed_exons) == 0) {
    pointed_poly_data <- empty_poly_data
  } else {
    pointed_poly_data <- .get_pointed_poly_data(pointed_exons, roi_length)
  }
  
  # There might be a situation where there are no non-pointed exons
  # - If there is only 1 exon per track and that exon is pointed
  if (nrow(normal_exons) == 0) {
    unpointed_poly_data <- empty_poly_data
  } else {
    unpointed_poly_data <- .get_unpointed_poly_data(normal_exons)
  }
  
  combined_poly_data <- dplyr::bind_rows(unpointed_poly_data, pointed_poly_data)
  
  return(combined_poly_data)
}

# Generate coordinates for pointed exons, genes, and transcripts
.get_pointed_poly_data <- function(df, roi_length) {
  
  tip_size <- 0.1 * roi_length # Fraction of roi length being plotted to determine arrow size
  track_height <- 0.2
  
  
  pointed_poly_data <- df %>%
    dplyr::rowwise() %>%
    dplyr::mutate(coords = list({
      width <- (end - start)
      #tip <- max(width * tip_size, 50)  # at least 50 bp tip
      # Arrows should be 10% of the region in length unless the feature length * .5 is smaller.
      # That would make the features/arrows look odd, so fallback to half the feature width
      tip <- ifelse(tip_size > width * .5, width * .5, tip_size)
      end_mod <- end - tip
      start_mod <- start + tip
      
      if (strand == "+") {
        # Right-pointing trapezoid
        data.frame(
          x = c(start, end_mod, end, end_mod, start),
          y = c(
            row_index + track_height,
            row_index + track_height,
            row_index, row_index - track_height,
            row_index - track_height)
        )
      } else {
        # Left-pointing trapezoid
        data.frame(
          x = c(start, start_mod, end, end, start_mod),
          y = c(
            row_index,
            row_index + track_height,
            row_index + track_height,
            row_index - track_height,
            row_index - track_height
          )
        )
      }
    })) %>%
    tidyr::unnest(coords) %>%
    dplyr::distinct()
  
  return(pointed_poly_data)
}

# Generate coordinates for unpointed exons
.get_unpointed_poly_data <- function(df) {
  
  track_height <- 0.2
  
  unpointed_poly_data <- df %>%
    dplyr::rowwise() %>%
    dplyr::mutate(coords = list({
      data.frame(
        x = c(start, end, end, start),
        y = c(
          row_index + track_height,
          row_index + track_height,
          row_index - track_height,
          row_index - track_height
        )
      )
    })) %>%
    tidyr::unnest(coords) %>%
    dplyr::distinct()
  
  return(unpointed_poly_data)
}
