# Plots annotated features over region of interest
# @param gtf A nine column annotation file
# @param reg_start an integer specifying the start location of region of interest
# @param reg_stop an integer specifying the end location of region of interest
# @return gtf_plot

.plot_exons_only <- function(gtf, reg_start, reg_stop){
  # Extract transcript ID from attribute field
  # Assign all transcripts to gene_id for easy matching
  na_counter <- 1
  na_idx <- vector()
  for (i in 1:nrow(gtf)) {
    # Gene ids should be present in all remaining observations
    tmp <- unlist(strsplit(gtf$attribute[i], ";"))
    gene_idx <- grep("gene_id", tmp)
    gene <- gsub("gene_id ", "", tmp[gene_idx])
    gtf$gene_id[i] <- gene
    #gtf$transcript_id[i] <- NA # By default, set to NA, will change if row is transcript or exon
    
    if (gtf$feature[i] == "transcript" | gtf$feature[i] == "exon") {
      # Transcript ids should only be present in transcript and exon rows
      trans_idx <- grep(" transcript_id", tmp)
      trans <- tmp[trans_idx]
      gtf$transcript_id[i] <- gsub(" transcript_id ", "", trans)
    } else {
      transcript_id <- paste0("NA_", na_counter)
      gtf$transcript_id[i] <- transcript_id
      na_idx <- append(na_idx, i)
      na_counter <- na_counter + 1
    }
  }
  
  # Separate exons
  exons <- gtf %>%
    dplyr::filter(feature == "exon")
  
  #Group exons belonging to same transcript and assign transcript index
  plot_df <- gtf %>%
    dplyr::mutate(trans_idx = match(transcript_id, unique(transcript_id))) %>%
    dplyr::select(-c(score, frame, attribute))
  
  # Revert numbered NAs to actual NAs
  plot_df$transcript_id[na_idx] <- NA
  
  mx <- max(exons$end)
  mn <- min(exons$start)
  
  if(length(mn) > 1){
    mn <- mn[1]
  }
  
  if(length(mx) > 1){
    mx <- mx[1]
  }
  # calculate the x and y coordinates for the exons/transcripts/genes
  transcripts <- plot_df %>%
    dplyr::filter(!is.na(transcript_id)) %>%
    dplyr::distinct(transcript_id, gene_id, strand, .keep_all = TRUE) %>%
    dplyr::arrange(transcript_id) %>%
    dplyr::mutate(track = dplyr::row_number() + length(unique(plot_df$transcript_id)), midpoint = ((mx - mn)/2 + reg_start)) #make midpoint for exons for text label
  
  # exons with transcript info
  exons <- plot_df %>%
    dplyr::filter(feature == "exon") %>%
    dplyr::left_join(transcripts %>% dplyr::select(transcript_id, track, strand),
                     by = "transcript_id") %>%
    dplyr::arrange(track, start) %>%
    dplyr::group_by(transcript_id) %>%
    dplyr::rename("strand" = "strand.x") %>%
    dplyr::mutate(
      exon_idx  = dplyr::row_number(),
      n_exons   = dplyr::n(),
      pointed   = dplyr::case_when(
        strand == "+" & exon_idx == n_exons ~ TRUE,
        strand == "-" & exon_idx == 1 ~ TRUE,
        TRUE ~ FALSE),
      next_start = dplyr::lead(start)) %>%
    dplyr::ungroup()
  
  exons <- exons %>%
    dplyr::mutate(midpoint = (end - start)/2 + start,
    start = ifelse(start < reg_start, reg_start, start),
    end = ifelse(end > reg_stop, reg_stop, end),
    midpoint = pmin(pmax(midpoint, reg_start), reg_stop))

  # keep separate datasets
  exons_normal <- exons %>% dplyr::filter(!pointed)   # only normal exons
  exons_pointed <- exons %>% dplyr::filter(pointed)   # only the directional ones
  
  # make polygons for pointed exons
  tip_size <- 0.2  # fraction of exon width to convert to tip
  exon_height <- 0.2
  poly_data <- exons_pointed %>%
    dplyr::rowwise() %>%
    dplyr::mutate(coords = list({
      width <- (end - start)
      tip <- max(width * tip_size, 50)  # at least 50 bp tip
      
      if (strand == "+") {
        # Right-pointing trapezoid
        data.frame(
          x = c(start, end - tip, end, end - tip, start),
          y = c(track - exon_height, track - exon_height,
                track, track + exon_height, track + exon_height))
      } else {
        # Left-pointing trapezoid
        data.frame(
          x = c(end, start + tip, start, start + tip, end),
          y = c(track - exon_height, track - exon_height,
                track, track + exon_height, track + exon_height))
      }})) %>%
    tidyr::unnest(coords) %>%
    dplyr::distinct()
  
  # plot
  gtf_plot <- ggplot2::ggplot()+
    # draw the line that connects exons
    ggplot2::geom_segment(
      data = exons %>% dplyr::filter(!is.na(next_start)),
      ggplot2::aes(x = end, xend = next_start,
                   y = track, yend = track),
      color = "black", linewidth = 0.3) +
    #normal exons
    ggplot2::geom_rect(
      data = exons_normal,
      ggplot2::aes(xmin = start, xmax = end,
                   ymin = track - exon_height, ymax = track + exon_height, fill = strand),
      color = "black") +
    # pointed exons
    ggplot2::geom_polygon(
      data = poly_data,
      ggplot2::aes(x = x, y = y, group = interaction(transcript_id, exon_idx), fill = strand),
      color = "black") +
    
    # transcript labels inside plot
    ggplot2::geom_text(
      data = transcripts,
      ggplot2::aes(x = midpoint, y = track, label = transcript_id),
      hjust = 1, vjust = 0.5, size = 3) +
    
    
    ggplot2::scale_y_continuous(limits = c(0.5, max(transcripts$track) + 0.5),
                                breaks = NULL, labels = NULL) +
    ggplot2::scale_x_continuous(limits = c(reg_start, reg_stop))+
    
    ggplot2::scale_fill_manual(values = c("+" = "lightcoral", "-" = "skyblue")) +
    ggplot2::labs(x = "Genomic coordinate") +
    ggplot2::theme_classic() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
      axis.title.y = ggplot2::element_blank())

    return(gtf_plot)
}
  
# Plots the coverage over an interval
# @param df a dataframe
# @return a histogram plot
.plot_transcripts_only <- function(gtf, reg_start, reg_stop){
  
  # Separate transcripts and genes
  transcripts <- gtf %>%
    dplyr::filter(feature == "transcript") %>%
    dplyr::distinct()
  
  
  if (nrow(transcripts) > 30) {
    message("Warning: region contains more features than are plottable. Selecting only first 30. See log file for details.")
    cat(file = logfile, "Region contains more features than are plottable Selecting only first 30. 
     Please consider filtering input GTF file (e.g. select only genes or transcripts or select transcript isoform of interest). \n", append = TRUE)
    transcripts <- transcripts[1:30,]
  }
  
  # Extract transcript_id and gene_id from attribute field
  transcripts$transcript_id <- sub('.*transcript_id "([^"]+)".*', '\\1', transcripts$attribute)
  transcripts$gene_id <- sub('.*gene_id "([^"]+)".*', '\\1', transcripts$attribute)
  transcripts$gene_id <- stringr::str_extract(transcripts$gene_id, "(?<=gene_id )[^;]+")
  
  # Keep only transcripts (no polygon building here)
  transcripts <- transcripts %>%
    dplyr::filter(feature == "transcript") %>%
    dplyr::distinct(transcript_id, .keep_all = TRUE) %>%
    dplyr::mutate(track = dplyr::row_number()) %>%
    dplyr::select(transcript_id, gene_id, start, end, strand, track) %>%
    dplyr::mutate(
      start = ifelse(start < reg_start, reg_start, start),
      end   = ifelse(end > reg_stop, reg_stop, end)
    )
  
   # make polygons for pointed exons
  tip_size <- 0.2  # fraction of exon width to convert to tip
  exon_height <- 0.2
  
  transcripts <- transcripts %>%
    dplyr::mutate(midpoint = (end - start)/2 + start,
    midpoint = pmin(pmax(midpoint, reg_start), reg_stop))  

  # Build polygons
  poly_data <- transcripts %>%
    dplyr::rowwise() %>%
    dplyr::mutate(coords = list({
      width <- end - start
      tip <- max(width * tip_size, 50)  # at least 50 bp tip
      
      if (strand == "+") {
        # Right-pointing trapezoid
        data.frame(
          x = c(start, end - tip, end, end - tip, start, start),
          y = c(track - exon_height, track - exon_height,
                track, track + exon_height, track + exon_height, track - exon_height)
        )
      } else {
        # Left-pointing trapezoid
        data.frame(
          x = c(end, start + tip, start, start + tip, end, end),
          y = c(track - exon_height, track - exon_height,
                track, track + exon_height, track + exon_height, track - exon_height)
        )
      }
    })) %>%
    tidyr::unnest(coords)
  
  # Plot transcripts as pointed boxes
  gtf_plot <- ggplot2::ggplot() +
    ggplot2::geom_polygon(
      data = poly_data,
      ggplot2::aes(x = x, y = y, group = transcript_id, fill = strand),
      color = "black") +
    
    ggplot2::geom_text(
      data = transcripts %>% dplyr::mutate(midpoint = (end + start) / 2),
      ggplot2::aes(x = midpoint, y = track, label = gene_id),
      vjust = 0.5, size = 3
    ) +
    ggplot2::scale_y_continuous(
      limits = c(0.5, max(transcripts$track) + 0.5),# add padding
      breaks = NULL, labels = NULL) +
    ggplot2::scale_x_continuous(limits = c(reg_start, reg_stop))+
    
    ggplot2::scale_fill_manual(values = c("+" = "lightcoral", "-" = "skyblue")) +
    ggplot2::labs(x = "Genomic coordinate", y = "Transcript") +
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
  #Group exons belonging to same transcript and assign transcript index
  
  # Separate transcripts and genes
  plot_df <- gtf %>%
    dplyr::mutate(gene_idx = match(gene_id, unique(gene_id))) %>%
    dplyr::select(-c(score, frame, attribute))
  
  # calculate the x and y coordinates for the exons/transcripts/genes
  genes <- plot_df %>%
    dplyr::filter(!is.na(gene_id)) %>%
    dplyr::distinct(gene_id, strand, .keep_all = TRUE) %>%
    dplyr::arrange(gene_id) %>%
    dplyr::mutate(track = dplyr::row_number(), midpoint = ((end - start)/2)+reg_start) #make genes plot first
  

  genes <- genes %>%
    dplyr::mutate(start = dplyr::case_when(
      start < reg_start ~ reg_start, # if start is less than cutoff, set it to cutoff
      TRUE ~ start)) %>%
    dplyr::mutate(end = dplyr::case_when(
      end > reg_stop ~ reg_stop, # if stop is greater than cutoff, set it to reg_stop
      TRUE ~ end))
  
    # make polygons for pointed exons
  tip_size <- 0.2  # fraction of exon width to convert to tip
  exon_height <- 0.2

  # genes -> make polygons with tip like transcripts
  gene_poly <- genes %>%
    dplyr::rowwise() %>%
    dplyr::mutate(coords = list({
      width <- end - start
      tip <- max(width * tip_size, 50)  # at least 50 bp tip
      
      if (strand == "+") {
        # right-pointing
        data.frame(
          x = c(start, end - tip, end, end - tip, start),
          y = c(track - exon_height, track - exon_height,
                track, track + exon_height, track + exon_height)
        )
      } else {
        # left-pointing
        data.frame(
          x = c(end, start + tip, start, start + tip, end),
          y = c(track - exon_height, track - exon_height,
                track, track + exon_height, track + exon_height)
        )
      }
    })) %>%
    tidyr::unnest(coords)
  
  
  # plot
  gtf_plot <- ggplot2::ggplot() +
    ggplot2::geom_polygon(
      data = gene_poly,
      ggplot2::aes(x = x, y = y, group = gene_id, fill = strand),
      color = "black") +
    
    ggplot2::geom_text(
      data = genes, 
      ggplot2::aes(x = midpoint, y = track, label = gene_id),
      vjust = 0.5, size = 3) +
    
    ggplot2::scale_y_continuous(limits = c(0.5, max(genes$track) + 0.5),
                                breaks = NULL, labels = NULL) +
    ggplot2::scale_x_continuous(limits = c(reg_start, reg_stop))+
    
    ggplot2::scale_fill_manual(values = c("+" = "lightcoral", "-" = "skyblue")) +
    ggplot2::labs(x = "Genomic coordinate") +
    ggplot2::theme_classic() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
      axis.title.y = ggplot2::element_blank())

return(gtf_plot) 

}

# Plots the coverage over an interval
# @param df a dataframe
# @return a histogram plot
.plot_transcripts_exons <- function(gtf, reg_start, reg_stop){
  
  # if transcript and exon information is both present
  # Extract transcript ID from attribute field
  # Assign all transcripts to gene_id for easy matching
  na_counter <- 1
  na_idx <- vector()
  for (i in 1:nrow(gtf)) {
    # Gene ids should be present in all remaining observations
    tmp <- unlist(strsplit(gtf$attribute[i], ";"))
    gene_idx <- grep("gene_id", tmp)
    gene <- gsub("gene_id ", "", tmp[gene_idx])
    gtf$gene_id[i] <- gene
    #gtf$transcript_id[i] <- NA # By default, set to NA, will change if row is transcript or exon
    
    if (gtf$feature[i] == "transcript" | gtf$feature[i] == "exon") {
      # Transcript ids should only be present in transcript and exon rows
      trans_idx <- grep(" transcript_id", tmp)
      trans <- tmp[trans_idx]
      gtf$transcript_id[i] <- gsub(" transcript_id ", "", trans)
    } else {
      transcript_id <- paste0("NA_", na_counter)
      gtf$transcript_id[i] <- transcript_id
      na_idx <- append(na_idx, i)
      na_counter <- na_counter + 1
    }
  }
  
  
  # Separate exons
  exons <- gtf %>%
    dplyr::filter(feature == "exon")
  
  # Separate transcripts and genes
  transcripts <- gtf %>%
    dplyr::filter(feature == "transcript") %>%
    dplyr::distinct()
  
  if (nrow(transcripts) > 30) {
    message("Warning: region contains more features than are plottable. Selecting only first 30. See log file for details.")
    cat(file = logfile, "Region contains more features than are plottable Selecting only first 30. 
   Please consider filtering input GTF file (e.g. select only genes or transcripts or select transcript isoform of interest). \n", append = TRUE)
    transcripts <- transcripts[1:30,]
  }
  
  transcripts$transcript_id <- sub('.*transcript_id "([^"]+)".*', '\\1', transcripts$attribute)
  
  plot_df <- gtf %>%
    dplyr::mutate(trans_idx = match(transcript_id, unique(transcript_id))) %>%
    dplyr::select(-c(score, frame, attribute))
  
  # Revert numbered NAs to actual NAs
  plot_df$transcript_id[na_idx] <- NA

  
  # calculate the x and y coordinates for the exons/transcripts/genes
  transcripts <- plot_df %>%
    dplyr::filter(!is.na(transcript_id)) %>%
    dplyr::distinct(transcript_id, gene_id, strand, .keep_all = TRUE) %>%
    dplyr::arrange(transcript_id) %>%
    dplyr::mutate(track = dplyr::row_number() + nrow(transcripts)) %>%
    dplyr::mutate(
      midpoint = (end - start)/2 + start,
      start = ifelse(start < reg_start, reg_start, start),
      end = ifelse(end > reg_stop, reg_stop, end),
      midpoint = pmin(pmax(midpoint, reg_start), reg_stop))
  
 
  exons <- plot_df %>%
    dplyr::filter(feature == "exon") %>%
    dplyr::mutate(track = dplyr::row_number()) %>%
    dplyr::arrange(track, start) %>%
    dplyr::group_by(transcript_id) %>%
    dplyr::mutate(
      exon_idx= dplyr::row_number(),
      n_exons = dplyr::n(),
      pointed = dplyr::case_when(
        strand == "+" & exon_idx == n_exons ~ TRUE,   # last exon for +
        strand == "-" & exon_idx == 1 ~ TRUE,   # first exon for -
        TRUE ~ FALSE),
      next_start = dplyr::lead(start)) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(
      midpoint = (end - start)/2 + start,
      start = ifelse(start < reg_start, reg_start, start),
      end = ifelse(end > reg_stop, reg_stop, end),
      midpoint = pmin(pmax(midpoint, reg_start), reg_stop))

  # keep separate datasets
  exons_normal <- exons %>% dplyr::filter(!pointed)   # only normal exons
  exons_pointed <- exons %>% dplyr::filter(pointed)   # only the directional ones

  
  # make polygons for pointed exons
  tip_size <- 0.2  # fraction of exon width to convert to tip
  exon_height <- 0.2
  poly_data <- exons_pointed %>%
    dplyr::rowwise() %>%
    dplyr::mutate(coords = list({
      width <- (end - start)
      tip <- max(width * tip_size, 50)  # at least 50 bp tip
      end_tip <- min(end + tip, reg_stop)
      start_tip <- max(start - tip, reg_start)
      
      if (strand == "+") {
        # Right-pointing trapezoid
        data.frame(
          x = c(start, end, end_tip, end, start),
          y = c(track + exon_height, track + exon_height,
                track, track - exon_height, track - exon_height))
      } else {
        # Left-pointing trapezoid
        data.frame(
          x = c(start_tip, start, end, end, start),
          y = c(track, track + exon_height, track + exon_height,
                track - exon_height, track - exon_height))
      }})) %>%
    tidyr::unnest(coords) %>%
    dplyr::distinct()
  
  # transcripts -> make polygons with tip like transcripts
  transcript_poly <- transcripts %>%
    dplyr::rowwise() %>%
    dplyr::mutate(coords = list({
      width <- end - start
      tip <- max(width * tip_size, 50)  # at least 50 bp tip
      end_tip <- min(end + tip, reg_stop)
      start_tip <- max(start - tip, reg_start)
      
      if (strand == "+") {
        # Right-pointing trapezoid
        data.frame(
          x = c(start, end, end_tip, end, start),
          y = c(track + exon_height, track + exon_height,
                track, track - exon_height, track - exon_height))
      } else {
        # Left-pointing trapezoid
        data.frame(
          x = c(start_tip, start, end, end, start),
          y = c(track, track + exon_height, track + exon_height,
                track - exon_height, track - exon_height))
      }
    })) %>%
    tidyr::unnest(coords)
  
  # plot
  gtf_plot <- ggplot2::ggplot() +
    # draw the line that connects exons
    ggplot2::geom_segment(
      data = exons %>% dplyr::filter(!is.na(next_start)),
      ggplot2::aes(x = end, xend = next_start,
                   y = track, yend = track),
      color = "black", linewidth = 0.3) +
    #normal exons
    ggplot2::geom_rect(
      data = exons_normal,
      ggplot2::aes(xmin = start, xmax = end,
                   ymin = track - exon_height, ymax = track + exon_height, fill = strand),
      color = "black") +
    
    ggplot2::geom_text(
      data = exons_normal,
      ggplot2::aes(x = midpoint, y = track, label = gene_id),
      hjust = 1, vjust = 0.5, size = 3) +
    
   
    # pointed exons
    ggplot2::geom_polygon(
      data = poly_data,
      ggplot2::aes(x = x, y = y, group = interaction(transcript_id, exon_idx), fill = strand),
      color = "black") +
    
    ggplot2::geom_text(
      data = exons_pointed,
      ggplot2::aes(x = midpoint, y = track, label = gene_id),
      hjust = 1, vjust = 0.5, size = 3) +
    
    
    ggplot2::geom_polygon(
      data = transcript_poly,
      ggplot2::aes(x = x, y = y, group = interaction(transcript_id, gene_id), fill = strand),
      color = "black")+
    
    # transcript labels inside plot
    ggplot2::geom_text(
      data = transcripts,
      ggplot2::aes(x = midpoint, y = track, label = transcript_id),
      hjust = 1, vjust = 0.5, size = 3) +
    
    
    ggplot2::scale_y_continuous(limits = c(0.5, max(transcripts$track) + 0.5),
                                breaks = NULL, labels = NULL) +
    ggplot2::scale_x_continuous(limits = c(reg_start, reg_stop))+
    
    ggplot2::scale_fill_manual(values = c("+" = "lightcoral", "-" = "skyblue")) +
    ggplot2::labs(x = "Genomic coordinate") +
    ggplot2::theme_classic() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
      axis.title.y = ggplot2::element_blank())  
}

# Plots the coverage over an interval
# @param df a dataframe
# @return a histogram plot
.plot_genes_exons <- function(gtf, reg_start, reg_stop){
  
  # if transcript and exon information is both present
  # Extract transcript ID from attribute field
  # Assign all transcripts to gene_id for easy matching
  
  na_counter <- 1
  na_idx <- vector()
  for (i in 1:nrow(gtf)) {
    # Gene ids should be present in all remaining observations
    tmp <- unlist(strsplit(gtf$attribute[i], ";"))
    gene_idx <- grep("gene_id", tmp)
    gene <- gsub("gene_id ", "", tmp[gene_idx])
    
    trans_idx <- grep(" transcript_id ", tmp)
    if(!identical(trans_idx, integer(0))){
      trans <- gsub(" transcript_id ", "", tmp[trans_idx])
      gtf$trans_id[i] <- trans
    } else {
      gtf$trans_id[i] <- "NA_"
    }
    gtf$gene_id[i] <- gene
    #gtf$transcript_id[i] <- NA # By default, set to NA, will change if row is transcript or exon
    
    if (gtf$feature[i] == "gene" | gtf$feature[i] == "exon") {
      # Transcript ids should only be present in transcript and exon rows
      gene_idx <- grep("gene_id", tmp)
      gene <- tmp[gene_idx]
      gtf$gene_id[i] <- gsub("gene_id ", "", gene)
    } else {
      gene_id <- paste0("NA_", na_counter)
      gtf$gene_id[i] <- gene_id
      na_idx <- append(na_idx, i)
      na_counter <- na_counter + 1
    }
  }
  
  
  # Separate exons
  exons <- gtf %>%
    dplyr::filter(feature == "exon")
  
  # Separate transcripts and genes
  genes <- gtf %>%
    dplyr::filter(feature == "gene") %>%
    dplyr::distinct()
  
  if (nrow(genes) > 30) {
    message("Warning: region contains more features than are plottable. Selecting only first 30. See log file for details.")
    cat(file = logfile, "Region contains more features than are plottable Selecting only first 30. 
   Please consider filtering input GTF file (e.g. select only genes or transcripts or select transcript isoform of interest). \n", append = TRUE)
    genes <- genes[1:30,]
  }
  
  genes$gene_id <- sub('.*gene_id "([^"]+)".*', '\\1', genes$attribute)
  
  plot_df <- gtf %>%
    dplyr::mutate(gene_idx = match(gene_id, unique(gene_id))) %>%
    dplyr::select(-c(score, frame, attribute))
  
  # Revert numbered NAs to actual NAs
  plot_df$gene_id[na_idx] <- NA
  plot_df[is.na(plot_df$gene_id)] <- NA
  
  # calculate the x and y coordinates for the exons/transcripts/genes
  genes <- plot_df %>%
    dplyr::filter(!is.na(gene_id)) %>%
    dplyr::distinct(gene_id, strand, .keep_all = TRUE) %>%
    dplyr::arrange(gene_id) %>%
    dplyr::mutate(track = dplyr::row_number() + nrow(exons)) %>%
    dplyr::mutate(
      midpoint = (end - start)/2 + start,
      start = ifelse(start < reg_start, reg_start, start),
      end = ifelse(end > reg_stop, reg_stop, end),
      midpoint = pmin(pmax(midpoint, reg_start), reg_stop))

  
  exons <- plot_df %>%
    dplyr::filter(feature == "exon") %>%
    dplyr::mutate(track = dplyr::row_number()) %>%
    dplyr::arrange(track, start) %>%
    dplyr::group_by(gene_id, trans_id) %>%
    dplyr::mutate(
      exon_idx= dplyr::row_number(),
      n_exons = dplyr::n(),
      pointed = dplyr::case_when(
        strand == "+" & exon_idx == n_exons ~ TRUE,   # last exon for +
        strand == "-" & exon_idx == 1 ~ TRUE,   # first exon for -
        TRUE ~ FALSE),
      next_start = dplyr::lead(start)) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(
      midpoint = (end - start)/2 + start,
      start = ifelse(start < reg_start, reg_start, start),
      end = ifelse(end > reg_stop, reg_stop, end),
      midpoint = pmin(pmax(midpoint, reg_start), reg_stop))
  
  
  # keep separate datasets
  exons_normal <- exons %>% dplyr::filter(!pointed)   # only normal exons
  exons_pointed <- exons %>% dplyr::filter(pointed)   # only the directional ones

  
  # make polygons for pointed exons
  tip_size <- 0.2  # fraction of exon width to convert to tip
  exon_height <- 0.2
  poly_data <- exons_pointed %>%
    dplyr::rowwise() %>%
    dplyr::mutate(coords = list({
      width <- (end - start)
      
      tip <- max(width * tip_size, 50)  # at least 50 bp tip
      end_tip <- min(end + tip, reg_stop)
      start_tip <- max(start - tip, reg_start)
  
      if (strand == "+") {
        # Right-pointing trapezoid
        data.frame(
          x = c(start, end, end_tip, end, start),
          y = c(track + exon_height, track + exon_height,
                track, track - exon_height, track - exon_height))
      } else {
        # Left-pointing trapezoid
        data.frame(
          x = c(start_tip, start, end, end, start),
          y = c(track, track + exon_height, track + exon_height,
                track - exon_height, track - exon_height, track))
      }})) %>%
    tidyr::unnest(coords) %>%
    dplyr::distinct()
  
  # genes -> make polygons with tip like transcripts
  gene_poly <- genes %>%
    dplyr::rowwise() %>%
    dplyr::mutate(coords = list({
      width <- end - start
      tip <- max(width * tip_size, 50)  # at least 50 bp tip
      end_tip <- min(end + tip, reg_stop)
      start_tip <- max(start - tip, reg_start)
      
      if (strand == "+") {
        # Right-pointing trapezoid
        data.frame(
          x = c(start, end, end_tip, end, start),
          y = c(track + exon_height, track + exon_height,
                track, track - exon_height, track - exon_height))
      } else {
        # Left-pointing trapezoid
        data.frame(
          x = c(start_tip, start, end, end, start),
          y = c(track, track + exon_height, track + exon_height,
                track - exon_height, track - exon_height, track))
      }
    })) %>%
    tidyr::unnest(coords)
  
  
  
  # plot
  gtf_plot <- ggplot2::ggplot() +
    # draw the line that connects exons
    ggplot2::geom_segment(
      data = exons %>% dplyr::filter(!is.na(next_start)),
      ggplot2::aes(x = end, xend = next_start,
                   y = track, yend = track),
      color = "black", linewidth = 0.3) +
    #normal exons
    ggplot2::geom_rect(
      data = exons_normal,
      ggplot2::aes(xmin = start, xmax = end,
                   ymin = track - exon_height, ymax = track + exon_height, fill = strand),
      color = "black") +
    
    ggplot2::geom_text(
      data = exons_normal,
      ggplot2::aes(x = midpoint, y = track, label = gene_id),
      hjust = 1, vjust = 0.5, size = 3) +
    
    
    # pointed exons
    ggplot2::geom_polygon(
      data = poly_data,
      ggplot2::aes(x = x, y = y, group = interaction(gene_id, trans_id), fill = strand),
      color = "black") +
    
    ggplot2::geom_text(
      data = exons_pointed,
      ggplot2::aes(x = midpoint, y = track, label = gene_id),
      hjust = 1, vjust = 0.5, size = 3) +
    
    
    ggplot2::geom_polygon(
      data = gene_poly,
      ggplot2::aes(x = x, y = y, fill = strand),
      color = "black")+
    
    # transcript labels inside plot
    ggplot2::geom_text(
      data = genes,
      ggplot2::aes(x = midpoint, y = track, label = gene_id),
      hjust = 1, vjust = 0.5, size = 3) +
    
    
    ggplot2::scale_y_continuous(limits = c(0.5, max(genes$track, exons$track) + 0.5),
                                breaks = NULL, labels = NULL) +
    ggplot2::scale_x_continuous(limits = c(reg_start, reg_stop))+
    
    ggplot2::scale_fill_manual(values = c("+" = "lightcoral", "-" = "skyblue")) +
    ggplot2::labs(x = "Genomic coordinate") +
    ggplot2::theme_classic() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
      axis.title.y = ggplot2::element_blank())  
    
  
  
}

# Plots the coverage over an interval
# @param df a dataframe
# @return a histogram plot
.plot_genes_exon_transcripts <- function(gtf, reg_start, reg_stop){
  
  # Assign all transcripts to gene_id for easy matching
  na_counter <- 1
  na_idx <- vector()
  for (i in 1:nrow(gtf)) {
    # Gene ids should be present in all remaining observations
    tmp <- unlist(strsplit(gtf$attribute[i], ";"))
    gene_idx <- grep("gene_id", tmp)
    gene <- gsub("gene_id ", "", tmp[gene_idx])
    gtf$gene_id[i] <- gene
    #gtf$transcript_id[i] <- NA # By default, set to NA, will change if row is transcript or exon
    
    if (gtf$feature[i] == "transcript" | gtf$feature[i] == "exon") {
      # Transcript ids should only be present in transcript and exon rows
      trans_idx <- grep(" transcript_id", tmp)
      trans <- tmp[trans_idx]
      gtf$transcript_id[i] <- gsub(" transcript_id ", "", trans)
    } else {
      transcript_id <- paste0("NA_", na_counter)
      gtf$transcript_id[i] <- transcript_id
      na_idx <- append(na_idx, i)
      na_counter <- na_counter + 1
    }
  }
  
  # Separate exons
  exons <- gtf %>%
    dplyr::filter(feature == "exon")
  
  # Separate transcripts and genes
  genes_transcripts <- gtf %>%
    dplyr::filter(feature == "transcript" | feature == "gene") %>%
    dplyr::distinct()
  
  transcripts <- gtf %>% dplyr::filter(feature == "transcript") %>%
    dplyr::distinct()
  # TO DO in all subfunctions check that rows of gtf do not exceed 30
  if (nrow(genes_transcripts) > 30) {
    message("Warning: region contains more features than are plottable. Selecting only first 30. See log file for details.")
    cat(file = logfile, "Region contains more features than are plottable Selecting only first 30. 
     Please consider filtering input GTF file (e.g. select only genes or transcripts or select transcript isoform of interest). \n", append = TRUE)
    genes_transcripts <- genes_transcripts[1:30,]
  }

    #Group exons belonging to same transcript and assign transcript index
    plot_df <- gtf %>%
      dplyr::mutate(trans_idx = match(transcript_id, unique(transcript_id))) %>%
      dplyr::select(-c(score, frame, attribute))
    
    # Revert numbered NAs to actual NAs
    plot_df$transcript_id[na_idx] <- NA
    
    
    # calculate the x and y coordinates for the exons/transcripts/genes
    transcripts <- plot_df %>%
      dplyr::filter(!is.na(transcript_id)) %>%
      dplyr::distinct(transcript_id, gene_id, strand, .keep_all = TRUE) %>%
      dplyr::arrange(transcript_id) %>%
      #dplyr::mutate(track = dplyr::row_number() + nrow(transcripts)) %>%
      dplyr::mutate(track = nrow(genes_only) + length(unique(exons$transcript_id)) + dplyr::row_number()) %>%
      dplyr::mutate(
        midpoint = (end - start)/2 + start,
        start = ifelse(start < reg_start, reg_start, start),
        end = ifelse(end > reg_stop, reg_stop, end),
        midpoint = pmin(pmax(midpoint, reg_start), reg_stop))
   
    # genes with no transcript_id (gene-only features)
    genes_only <- plot_df %>%
      dplyr::filter(feature == "gene" & is.na(transcript_id)) %>%
      dplyr::mutate(track = min(transcripts$track, 0) + dplyr::row_number()) 
    
    genes_only <- genes_only %>%
      dplyr::mutate(midpoint = (end - start)/2 + start,
      start = ifelse(start < reg_start, reg_start, start),
      end = ifelse(end > reg_stop, reg_stop, end),
      midpoint = pmin(pmax(midpoint, reg_start), reg_stop))
    
    # exons with transcript info
    exons <- plot_df %>%
      dplyr::filter(feature == "exon") %>%
      dplyr::left_join(transcripts %>% dplyr::select(transcript_id, track, strand),
                       by = "transcript_id") %>%
      dplyr::arrange(track, start) %>%
      dplyr::group_by(transcript_id) %>%
      dplyr::rename("strand" = "strand.x") %>%
      dplyr::mutate(
        exon_idx  = dplyr::row_number(),
        track = nrow(genes_only) + length(unique(transcript_id)) + dplyr::cur_group_id(),
        n_exons   = dplyr::n(),
        pointed   = dplyr::case_when(
          strand == "+" & exon_idx == n_exons ~ TRUE,
          strand == "-" & exon_idx == 1 ~ TRUE,
          TRUE ~ FALSE),
        next_start = dplyr::lead(start)) %>%
      dplyr::ungroup() %>%
      dplyr::mutate(
        midpoint = (end - start)/2 + start,
        start = ifelse(start < reg_start, reg_start, start),
        end = ifelse(end > reg_stop, reg_stop, end),
        midpoint = pmin(pmax(midpoint, reg_start), reg_stop))
      
    
    exons_normal <- exons %>% dplyr::filter(!pointed)   # only normal exons
    exons_pointed <- exons %>% dplyr::filter(pointed)   # only the directional ones
    
    # make polygons for pointed exons
    tip_size <- 0.2  # fraction of exon width to convert to tip
    exon_height <- 0.2
    poly_data <- exons_pointed %>%
      dplyr::rowwise() %>%
      dplyr::mutate(coords = list({
        width <- (end - start)
        tip <- max(width * tip_size, 50)  # at least 50 bp tip
        end_tip <- min(end + tip, reg_stop)
        start_tip <- max(start - tip, reg_start)
        
        if (strand == "+") {
          # Right-pointing trapezoid
          data.frame(
            x = c(start, end, end_tip, end, start),
            y = c(track + exon_height, track + exon_height,
                  track, track - exon_height, track - exon_height))
        } else {
          # Left-pointing trapezoid
          data.frame(
            x = c(start_tip, start, end, end, start),
            y = c(track, track + exon_height, track + exon_height,
                  track - exon_height, track - exon_height))
        }})) %>%
      tidyr::unnest(coords) %>%
      dplyr::distinct()
    
    # genes -> make polygons with tip like transcripts
    gene_poly <- genes_only %>%
      dplyr::rowwise() %>%
      dplyr::mutate(coords = list({
        width <- end - start
        tip <- max(width * tip_size, 50)  # at least 50 bp tip
        end_tip <- min(end + tip, reg_stop)
        start_tip <- max(start - tip, reg_start)
        
        if (strand == "+") {
          # right-pointing
          data.frame(
            x = c(start, end, end_tip, end, start),
            y = c(track + exon_height, track + exon_height,
                  track, track - exon_height, track - exon_height)
          )
        } else {
          # left-pointing
          data.frame(
            x = c(start_tip, start, end, end, start),
            y = c(track, track + exon_height, track + exon_height,
                  track - exon_height, track - exon_height))
            #midpoint = (end - start)/2 + reg_start)
        }
      })) %>%
      tidyr::unnest(coords)
    
    
    # genes -> make polygons with tip like transcripts
    transcript_poly <- transcripts %>%
      dplyr::rowwise() %>%
      dplyr::mutate(coords = list({
        width <- end - start
        tip <- max(width * tip_size, 50)  # at least 50 bp tip
        end_tip <- min(end + tip, reg_stop)
        start_tip <- max(start - tip, reg_start)
        midpoint <- (end - start)/2 + start
    
        if (strand == "+") {
          # right-pointing
          data.frame(
            x = c(start, end, end_tip, end, start),
            y = c(track + exon_height, track + exon_height,
                  track, track - exon_height, track - exon_height))
        } else {
          # left-pointing
          data.frame(
            x = c(start_tip, start, end,  end, start),
            y = c(track, track + exon_height, track + exon_height,
                  track - exon_height, track - exon_height))
        }
      })) %>%
      tidyr::unnest(coords)
  
    # plot
    gtf_plot <- ggplot2::ggplot() +
      # draw the line that connects exons
      ggplot2::geom_polygon(
        data = gene_poly,
        ggplot2::aes(x = x, y = y, group = interaction(transcript_id, gene_id), fill = strand),
        color = "black")+
      
      ggplot2::geom_text(
        data = genes_only,
        ggplot2::aes(x= midpoint, y = track, label = gene_id),
        hjust = 1, vjust = 0.5, size = 3)+
      
      ggplot2::geom_segment(
        data = exons %>% dplyr::filter(!is.na(next_start)),
        ggplot2::aes(x = end, xend = next_start,
                     y = track, yend = track),
        color = "black", linewidth = 0.3) +
      #normal exons
      ggplot2::geom_rect(
        data = exons_normal,
        ggplot2::aes(xmin = start, xmax = end,
                     ymin = track - exon_height, ymax = track + exon_height, fill = strand),
        color = "black") +
      
      ggplot2::geom_text(
        data = exons_normal,
        ggplot2::aes(x = midpoint, y = track, label = transcript_id),
        hjust = 1, vjust = 0.5, size = 3) +
      
      
      # pointed exons
      ggplot2::geom_polygon(
        data = poly_data,
        ggplot2::aes(x = x, y = y, group = interaction(transcript_id, exon_idx), fill = strand),
        color = "black") +
      
      ggplot2::geom_text(
        data = exons_pointed,
        ggplot2::aes(x = midpoint, y = track, label = gene_id),
        hjust = 1, vjust = 0.5, size = 3) +
      
      
      ggplot2::geom_polygon(
        data = transcript_poly,
        ggplot2::aes(x = x, y = y, group = interaction(transcript_id, gene_id), fill = strand),
        color = "black")+
      
      # transcript labels inside plot
      ggplot2::geom_text(
        data = transcripts,
        ggplot2::aes(x = midpoint, y = track, label = transcript_id),
        hjust = 1, vjust = 0.5, size = 3) +
      
      
      ggplot2::scale_y_continuous(limits = c(0.5, max(transcripts$track) + 0.5),
                                  breaks = NULL, labels = NULL) +
      ggplot2::scale_x_continuous(limits = c(reg_start, reg_stop))+
      
      ggplot2::scale_fill_manual(values = c("+" = "lightcoral", "-" = "skyblue")) +
      ggplot2::labs(x = "Genomic coordinate") +
      ggplot2::theme_classic() +
      ggplot2::theme(
        axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
        axis.title.y = ggplot2::element_blank())  
    
    return(gtf_plot)
  
}
