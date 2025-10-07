# Plots annotated features over region of interest
# @param gtf A nine column annotation file
# @param reg_start an integer specifying the start location of region of interest
# @param reg_stop an integer specifying the end location of region of interest
# @return gtf_plot

.plot_exons_only <- function(gtf, reg_start, reg_stop){
  
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
    dplyr::mutate(transcript_idx = match(transcript_id, unique(transcript_id))) %>%
    dplyr::select(-c(score, frame, attribute))
  
  # Just in case, put rows in a good order for determining plot layout
  plot_df <- plot_df %>%
    dplyr::arrange(strand, gene_id, transcript_id, start, end)
  
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
  
  # make polygons for pointed exons
  tip_size <- 0.2  # fraction of exon width to convert to tip
  exon_height <- 0.2
  poly_data <- exons_pointed %>%
    dplyr::rowwise() %>%
    dplyr::mutate(coords = list({
      width <- (end - start)
      tip <- max(width * tip_size, 50)  # at least 50 bp tip
      end_mod <- end - tip
      start_mod <- start + tip
      end_tip <- min(end + tip, reg_stop)
      start_tip <- max(start - tip, reg_start)
      
      if (strand == "+") {
        # Right-pointing trapezoid
        data.frame(
          x = c(start, end_mod, end, end_mod, start),
          y = c(
            row_index + exon_height,
            row_index + exon_height,
            row_index, row_index - exon_height,
            row_index - exon_height)
        )
      } else {
        # Left-pointing trapezoid
        data.frame(
          x = c(start, start_mod, end, end, start_mod),
          y = c(
            row_index,
            row_index + exon_height,
            row_index + exon_height,
            row_index - exon_height,
            row_index - exon_height
          )
        )
      }
    })) %>%
    tidyr::unnest(coords) %>%
    dplyr::distinct()
  
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
    
    # Non pointed exons
    ggplot2::geom_rect(
      data = exons_normal,
      ggplot2::aes(
        xmin = start,
        xmax = end,
        ymin = row_index - exon_height,
        ymax = row_index + exon_height,
        fill = strand
      ),
      color = "black"
    ) +
    
    # Pointed exons
    ggplot2::geom_polygon(
      data = poly_data,
      ggplot2::aes(
        x = x,
        y = y,
        group = interaction(transcript_id, exon_idx),
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
.plot_transcripts_only <- function(gtf, reg_start, reg_stop){
  
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
  
  # Arrange rows for plotting
  # Feature precedence: gene > exon
  plot_df <- plot_df %>%
    dplyr::arrange(strand, gene_id, transcript_id, start, end)
  
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
  
  # make polygons for pointed exons
  tip_size <- 0.2  # fraction of exon width to convert to tip
  exon_height <- 0.2
  
  # transcripts -> make polygons with tip
  transcript_poly <- transcripts %>%
    dplyr::rowwise() %>%
    dplyr::mutate(coords = list({
      width <- end - start
      tip <- max(width * tip_size, 50)  # at least 50 bp tip
      end_mod <- end - tip
      start_mod <- start + tip
      
      if (strand == "+") {
        # Right-pointing trapezoid
        data.frame(
          x = c(start, end_mod, end, end_mod, start),
          y = c(
            row_index + exon_height,
            row_index + exon_height,
            row_index,
            row_index - exon_height,
            row_index - exon_height
          )
        )
      } else {
        # Left-pointing trapezoid
        data.frame(
          x = c(start, start_mod, end, end, start_mod),
          y = c(
            row_index,
            row_index + exon_height,
            row_index + exon_height,
            row_index - exon_height,
            row_index - exon_height
          )
        )
      }
    })) %>%
    tidyr::unnest(coords)
  
  TXT_SIZE <- 6
  
  # Generate Plot Object
  gtf_plot <- ggplot2::ggplot() +
    
    # Transcripts
    ggplot2::geom_polygon(
      data = transcript_poly,
      ggplot2::aes(
        x = x,
        y = y,
        group = interaction(transcript_id, gene_id),
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
.plot_genes_only <- function(gtf, reg_start, reg_stop){
  
  # Extract Gene ID from attribute field
  for (i in 1:nrow(gtf)) {
    # Gene ids should be present in all remaining observations
    tmp <- unlist(strsplit(gtf$attribute[i], ";"))
    gene_idx <- grep("gene_id", tmp)
    gene <- gsub("gene_id ", "", tmp[gene_idx])
    gtf$gene_id[i] <- gene
  }
  
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
  
  # make polygons for pointed exons
  tip_size <- 0.2  # fraction of exon width to convert to tip
  exon_height <- 0.2
  
  # genes -> make polygons with tip
  gene_poly <- genes %>%
    dplyr::rowwise() %>%
    dplyr::mutate(coords = list({
      width <- end - start
      tip <- max(width * tip_size, 50)  # at least 50 bp tip
      end_mod <- end - tip
      start_mod <- start + tip
      
      if (strand == "+") {
        # Right-pointing trapezoid
        data.frame(
          x = c(start, end_mod, end, end_mod, start),
          y = c(
            row_index + exon_height,
            row_index + exon_height,
            row_index,
            row_index - exon_height,
            row_index - exon_height
          )
        )
      } else {
        # Left-pointing trapezoid
        data.frame(
          x = c(start, start_mod, end, end, start_mod),
          y = c(
            row_index,
            row_index + exon_height,
            row_index + exon_height,
            row_index - exon_height,
            row_index - exon_height
          )
        )
      }
    })) %>%
    tidyr::unnest(coords)
  
  TXT_SIZE <- 6
  
  # Generate Plot Object
  gtf_plot <- ggplot2::ggplot() +

    # Genes
    ggplot2::geom_polygon(
      data = gene_poly,
      ggplot2::aes(
        x = x,
        y = y,
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
.plot_transcripts_exons <- function(gtf, reg_start, reg_stop){
  
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
    dplyr::mutate(transcript_idx = match(transcript_id, unique(transcript_id))) %>%
    dplyr::select(-c(score, frame, attribute))
  
  # Just in case, put rows in a good order for determining plot layout
  plot_df <- plot_df %>%
    dplyr::arrange(strand, transcript_id, desc(feature), start, end)

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
  
  # make polygons for pointed exons
  tip_size <- 0.2  # fraction of exon width to convert to tip
  exon_height <- 0.2
  poly_data <- exons_pointed %>%
    dplyr::rowwise() %>%
    dplyr::mutate(coords = list({
      width <- (end - start)
      tip <- max(width * tip_size, 50)  # at least 50 bp tip
      end_mod <- end - tip
      start_mod <- start + tip
      end_tip <- min(end + tip, reg_stop)
      start_tip <- max(start - tip, reg_start)
      
      if (strand == "+") {
        # Right-pointing trapezoid
        data.frame(
          x = c(start, end_mod, end, end_mod, start),
          y = c(
            row_index + exon_height,
            row_index + exon_height,
            row_index, row_index - exon_height,
            row_index - exon_height)
          )
      } else {
        # Left-pointing trapezoid
        data.frame(
          x = c(start, start_mod, end, end, start_mod),
          y = c(
            row_index,
            row_index + exon_height,
            row_index + exon_height,
            row_index - exon_height,
            row_index - exon_height
          )
        )
      }
    })) %>%
    tidyr::unnest(coords) %>%
    dplyr::distinct()
  
  # transcripts -> make polygons with tip like transcripts
  transcript_poly <- transcripts %>%
    dplyr::rowwise() %>%
    dplyr::mutate(coords = list({
      width <- end - start
      tip <- max(width * tip_size, 50)  # at least 50 bp tip
      end_mod <- end - tip
      start_mod <- start + tip
      #end_tip <- min(end + tip, reg_stop)
      #start_tip <- max(start - tip, reg_start)
      
      if (strand == "+") {
        # Right-pointing trapezoid
        data.frame(
          # To ensure that tip of arrow is at 
          x = c(start, end_mod, end, end_mod, start),
          y = c(
            row_index + exon_height,
            row_index + exon_height,
            row_index,
            row_index - exon_height,
            row_index - exon_height
          )
        )
      } else {
        # Left-pointing trapezoid
        data.frame(
          x = c(start, start_mod, end, end, start_mod),
          y = c(
            row_index,
            row_index + exon_height,
            row_index + exon_height,
            row_index - exon_height,
            row_index - exon_height
          )
        )
      }
    })) %>%
    tidyr::unnest(coords)
  
  
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
    
    # Non pointed exons
    ggplot2::geom_rect(
      data = exons_normal,
      ggplot2::aes(
        xmin = start,
        xmax = end,
        ymin = row_index - exon_height,
        ymax = row_index + exon_height,
        fill = strand
      ),
      color = "black"
    ) +
    
    # Pointed exons
    ggplot2::geom_polygon(
      data = poly_data,
      ggplot2::aes(
        x = x,
        y = y,
        group = interaction(transcript_id, exon_idx),
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
      data = transcript_poly,
      ggplot2::aes(
        x = x,
        y = y,
        group = interaction(transcript_id, gene_id),
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
.plot_genes_exons <- function(gtf, reg_start, reg_stop){

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
      transcript_idx = match(transcript_id, unique(transcript_id))
    ) %>%
    dplyr::select(-c(score, frame, attribute))
  
  
  # Arrange rows for plotting
  # Feature precedence: gene > exon
  plot_df <- plot_df %>%
    dplyr::arrange(strand, gene_id, desc(feature), transcript_id, start, end)
  
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

  
  # make polygons for pointed exons
  tip_size <- 0.2  # fraction of exon width to convert to tip
  exon_height <- 0.2
  
  poly_data <- exons_pointed %>%
    dplyr::rowwise() %>%
    dplyr::mutate(coords = list({
      width <- (end - start)
      tip <- max(width * tip_size, 50)  # at least 50 bp tip
      end_mod <- end - tip
      start_mod <- start + tip
  
      if (strand == "+") {
        # Right-pointing trapezoid
        data.frame(
          x = c(start, end_mod, end, end_mod, start),
          y = c(
            row_index + exon_height,
            row_index + exon_height,
            row_index,
            row_index - exon_height,
            row_index - exon_height
          )
        )
      } else {
        # Left-pointing trapezoid
        data.frame(
          x = c(start, start_mod, end, end, start_mod),
          y = c(
            row_index,
            row_index + exon_height,
            row_index + exon_height,
            row_index - exon_height,
            row_index - exon_height
          )
        )
      }
    })) %>%
    tidyr::unnest(coords) %>%
    dplyr::distinct()
  
  # genes -> make polygons with tip like transcripts
  gene_poly <- genes %>%
    dplyr::rowwise() %>%
    dplyr::mutate(coords = list({
      width <- end - start
      tip <- max(width * tip_size, 50)  # at least 50 bp tip
      end_mod <- end - tip
      start_mod <- start + tip
      
      if (strand == "+") {
        # Right-pointing trapezoid
        data.frame(
          x = c(start, end_mod, end, end_mod, start),
          y = c(
            row_index + exon_height,
            row_index + exon_height,
            row_index,
            row_index - exon_height,
            row_index - exon_height
          )
        )
      } else {
        # Left-pointing trapezoid
        data.frame(
          x = c(start, start_mod, end, end, start_mod),
          y = c(
            row_index,
            row_index + exon_height,
            row_index + exon_height,
            row_index - exon_height,
            row_index - exon_height
          )
        )
      }
    })) %>%
    tidyr::unnest(coords)
  
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
    
    # Non Pointed Exons
    ggplot2::geom_rect(
      data = exons_normal,
      ggplot2::aes(
        xmin = start,
        xmax = end,
        ymin = row_index - exon_height,
        ymax = row_index + exon_height,
        fill = strand
      ),
      color = "black"
    ) +
    
    # Pointed Exons
    ggplot2::geom_polygon(
      data = poly_data,
      ggplot2::aes(
        x = x,
        y = y,
        group = interaction(transcript_id, exon_idx),
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
      data = gene_poly,
      ggplot2::aes(
        x = x,
        y = y,
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
.plot_genes_exon_transcripts <- function(gtf, reg_start, reg_stop){
  
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
      transcript_idx = match(transcript_id, unique(transcript_id))
    ) %>%
    dplyr::select(-c(score, frame, attribute))
    
  # Modify "gene" features to "xgene" for temporary sorting
  plot_df$feature[plot_df$feature == "gene"] <- "xgene"
    
  # Arrange rows for plotting
  # Feature precedence: gene > transcript > exon
  plot_df <- plot_df %>%
    dplyr::arrange(strand, gene_id, desc(feature), transcript_id, start, end)
  
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
    
  # make polygons for pointed exons
  tip_size <- 0.2  # fraction of exon width to convert to tip
  exon_height <- 0.2
  
  poly_data <- exons_pointed %>%
    dplyr::rowwise() %>%
    dplyr::mutate(coords = list({
      width <- (end - start)
      tip <- max(width * tip_size, 50)  # at least 50 bp tip
      end_mod <- end - tip
      start_mod <- start + tip
      
      if (strand == "+") {
        # Right-pointing trapezoid
        data.frame(
          x = c(start, end_mod, end, end_mod, start),
          y = c(
            row_index + exon_height,
            row_index + exon_height,
            row_index,
            row_index - exon_height,
            row_index - exon_height
          )
        )
      } else {
        # Left-pointing trapezoid
        data.frame(
          x = c(start, start_mod, end, end, start_mod),
          y = c(
            row_index,
            row_index + exon_height,
            row_index + exon_height,
            row_index - exon_height,
            row_index - exon_height
          )
        )
      }
    })) %>%
    tidyr::unnest(coords) %>%
    dplyr::distinct()
  
  # genes -> make polygons with tip like transcripts
  gene_poly <- genes %>%
    dplyr::rowwise() %>%
    dplyr::mutate(coords = list({
      width <- end - start
      tip <- max(width * tip_size, 50)  # at least 50 bp tip
      end_mod <- end - tip
      start_mod <- start + tip
      
      if (strand == "+") {
        # Right-pointing trapezoid
        data.frame(
          x = c(start, end_mod, end, end_mod, start),
          y = c(
            row_index + exon_height,
            row_index + exon_height,
            row_index,
            row_index - exon_height,
            row_index - exon_height
          )
        )
      } else {
        # Left-pointing trapezoid
        data.frame(
          x = c(start, start_mod, end, end, start_mod),
          y = c(
            row_index,
            row_index + exon_height,
            row_index + exon_height,
            row_index - exon_height,
            row_index - exon_height
          )
        )
      }
    })) %>%
    tidyr::unnest(coords)
  
  # transcripts -> make polygons with tip like transcripts
  transcript_poly <- transcripts %>%
    dplyr::rowwise() %>%
    dplyr::mutate(coords = list({
      width <- end - start
      tip <- max(width * tip_size, 50)  # at least 50 bp tip
      end_mod <- end - tip
      start_mod <- start + tip
      #end_tip <- min(end + tip, reg_stop)
      #start_tip <- max(start - tip, reg_start)
      
      if (strand == "+") {
        # Right-pointing trapezoid
        data.frame(
          # To ensure that tip of arrow is at 
          x = c(start, end_mod, end, end_mod, start),
          y = c(
            row_index + exon_height,
            row_index + exon_height,
            row_index,
            row_index - exon_height,
            row_index - exon_height
          )
        )
      } else {
        # Left-pointing trapezoid
        data.frame(
          x = c(start, start_mod, end, end, start_mod),
          y = c(
            row_index,
            row_index + exon_height,
            row_index + exon_height,
            row_index - exon_height,
            row_index - exon_height
          )
        )
      }
    })) %>%
    tidyr::unnest(coords)
  
  TXT_SIZE <- 6

  # Generate plot object
  
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
    
    # Non Pointed Exons
    ggplot2::geom_rect(
      data = exons_normal,
      ggplot2::aes(
        xmin = start,
        xmax = end,
        ymin = row_index - exon_height,
        ymax = row_index + exon_height,
        fill = strand
      ),
      color = "black"
    ) +
    
    # Pointed Exons
    ggplot2::geom_polygon(
      data = poly_data,
      ggplot2::aes(
        x = x,
        y = y,
        group = interaction(transcript_id, exon_idx),
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
      data = gene_poly,
      ggplot2::aes(
        x = x,
        y = y,
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
      data = transcript_poly,
      ggplot2::aes(
        x = x,
        y = y,
        group = interaction(transcript_id, gene_id),
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
