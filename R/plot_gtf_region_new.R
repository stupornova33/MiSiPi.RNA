# takes a gtf file
# outputs a plot
# @param gtf_file a string path leading to the gtf file
# @param chrom_name a string
# @param reg_start an integer
# @param reg_stop an integer
# @param logfile
#
# @return plot

.plot_gtf_region_new <- function(gtf_file, chrom_name, reg_start, reg_stop, logfile) {
  # Load GTF
  gtf <- read.csv(gtf_file, header = FALSE, sep = "\t", fill = TRUE, quote="")
  colnames(gtf)[1:9] <- c("seqname","source","feature","start","end","score","strand","frame","attribute")
  
  # Filter for region
  gtf <- gtf %>% 
    dplyr::filter(seqname == chrom_name, 
           (start >= reg_start & start <= reg_stop) |
             (end   >= reg_start & end   <= reg_stop))  %>%
    dplyr::filter(feature != "start_codon" & feature != "CDS")
  
  
  if (nrow(gtf) == 0) {
    message("No features found in region.")
    return(NULL)
  }


  
  ### test with exons/introns
  # Separate exons and introns
  exons <- gtf %>% dplyr::filter(feature == "exon")
  
  # if only transcripts in gtf 
  if(nrow(exons) != 0){

    
  
  # Extract transcript ID from attribute field
  # Assign all transcripts to gene_id for easy matching
  for(i in 1:nrow(gtf)){
    if(gtf$feature[i] == "transcript" | gtf$feature[i] == "exon"){
      tmp <- unlist(strsplit(gtf$attribute[i], ";"))
      trans_idx <- grep(" transcript_id", tmp)
      trans <- tmp[trans_idx]
      gtf$transcript_id[i] <- gsub(" transcript_id ", "", trans)
      gene_idx <- grep("gene_id ", tmp)
      gene <- tmp[gene_idx]
      gtf$gene_id[i] <- gsub("gene_id ", "", gene)
      
    } else if(gtf$feature[i] == "gene") {
      tmp <- unlist(strsplit(gtf$attribute[i], ";"))
      gene_idx <- grep("gene_id", tmp)
      
      gene <- unlist(strsplit(gtf$attribute[i], ";"))[1]
      gtf$transcript_id[i] <- gsub("gene_id", "", gene) # if feature is gene, use gene name for plot
      gtf$gene_id[i] <- gsub("gene_id ", "", gene) 
      
    } 
  } 
  
  # Separate exons

  exons <- gtf %>% dplyr::filter(feature == "exon")
  genes_transcripts <- gtf %>% dplyr::filter(feature == "transcript" | feature == "gene") %>% dplyr::distinct()
   if(nrow(genes_transcripts) > 30){
     message("Warning: region contains more features than are plottable. Selecting only first 10. See log file for details.")
     cat(file = logfile, "Region contains more features than are plottable Selecting only first 10. 
     Please consider filtering input GTF file (e.g. select only genes or transcripts or select transcript isoform of interest). \n", append = TRUE)
     genes_transcripts <- genes_transcripts[1:30,]
   }
  
  #Group exons belonging to same transcript
  gene_list <- list()
  for(i in 1:nrow(genes_transcripts)){
    if(genes_transcripts$feature[i] == "transcript" | genes_transcripts$feature[i] == "exon"){
      idx <- which(exons$transcript_id == genes_transcripts$transcript_id[i])
      gene_list[[i]] <- c(exons[idx,])
    } else if(genes_transcripts$feature[i] == "gene"){
      gene_list[[i]] <- as.list(genes_transcripts[i,])
    }
  }
  
  #assign transcript number to exons for plotting
  plot_df <- data.frame()
  for(i in 1:length(gene_list)){
    df <- data.frame(gene_list[[i]])
    df$trans_idx <- i
    plot_df <- rbind(plot_df, df)
  }
  
  plot_df <- plot_df %>% dplyr::select(-c(score, frame, attribute))
  # calculate the x and y coordinates for the exons/transcripts/genes
  # Compute transcript spans

  # transcripts -> track positions
  transcripts <- plot_df %>%
    dplyr::filter(!is.na(transcript_id)) %>%
    dplyr::distinct(transcript_id, gene_id, strand, .keep_all = TRUE) %>%
    dplyr::arrange(transcript_id) %>%
    dplyr::mutate(track = dplyr::row_number())
  
  # exons -> add info on first/last per transcript
  exons <- plot_df %>%
    dplyr::filter(feature == "exon") %>%
    dplyr::left_join(transcripts %>% dplyr::select(transcript_id, track, strand), by = "transcript_id") %>%
    dplyr::rename("strand" = "strand.x") %>%
    dplyr::arrange(track, start) %>%
    dplyr::group_by(transcript_id) %>%
    dplyr::mutate(
      exon_idx  = dplyr::row_number(),
      n_exons   = dplyr::n(),
      pointed   = dplyr::case_when(
        strand == "+" & exon_idx == n_exons ~ TRUE,   # last exon for +
        strand == "-" & exon_idx == 1 ~ TRUE,   # first exon for -
        TRUE ~ FALSE),
      next_start = dplyr::lead(start)) %>%
    dplyr::ungroup()
  
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
      }
    })) %>%
    tidyr::unnest(coords)
  
  # plot
  gtf_plot <- ggplot2::ggplot() +
    # introns
    ggplot2::geom_segment(
      data = exons %>% dplyr::filter(!is.na(next_start)),
      ggplot2::aes(x = end, xend = next_start,
          y = track, yend = track),
      color = "black", linewidth = 0.3) +
    # normal exons
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
      ggplot2::aes(x = start, y = track, label = transcript_id),
      hjust = 1.1, vjust = 0.5, size = 3) +
    
    # remove y-axis labels since we now draw them manually
    ggplot2::scale_y_continuous(limits = c(0.5, max(transcripts$track) + 0.5),
                                breaks = NULL, labels = NULL) +
    
    ggplot2::scale_fill_manual(values = c("+" = "lightcoral", "-" = "skyblue")) +
    ggplot2::labs(x = "Genomic coordinate") +
    ggplot2::theme_classic() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
      axis.title.y = ggplot2::element_blank())
  
  } else {  #if exon information does not exist, only transcript
  
    
    genes_transcripts <- gtf %>% dplyr::filter(feature == "transcript" | feature == "gene") %>% dplyr::distinct()
     if(nrow(genes_transcripts) > 30){
       message("Warning: region contains more features than are plottable. Selecting only first 10. See log file for details.")
       cat(file = logfile, "Region contains more features than are plottable Selecting only first 10. 
     Please consider filtering input GTF file (e.g. select only genes or transcripts or select transcript isoform of interest). \n", append = TRUE)
       genes_transcripts <- genes_transcripts[1:30,]
     
     }
    
    # Extract transcript_id and gene_id from attribute field
    gtf$transcript_id <- sub('.*transcript_id "([^"]+)".*', '\\1', gtf$attribute)
    gtf$gene_id <- sub('.*gene_id "([^"]+)".*', '\\1', gtf$attribute)
    
    # Keep only transcript rows
    transcripts <- gtf %>%
      dplyr::filter(feature == "transcript") %>%
      dplyr::distinct(transcript_id, .keep_all = TRUE) %>%
      dplyr::mutate(track = dplyr::row_number()) %>%
      dplyr::select(transcript_id, gene_id, start, end, strand, track)
    
    # Parameters
    tip_size <- 0.2
    exon_height <- 0.2
    
    # Build polygons for transcript rectangles with tip
    poly_data <- transcripts %>%
      dplyr::rowwise() %>%
      dplyr::mutate(coords = list({
        width <- end - start
        tip <- max(width * tip_size, 50)  # at least 50 bp tip
        
        if (strand == "+") {
          # Rectangle with right-pointing tip
          data.frame(
            x = c(start, end - tip, end, end - tip, start),
            y = c(track - exon_height, track - exon_height,
                  track, track + exon_height, track + exon_height))
        } else {
          # Rectangle with left-pointing tip
          data.frame(
            x = c(end, start + tip, start, start + tip, end),
            y = c(track - exon_height, track - exon_height,
                  track, track + exon_height, track + exon_height))
        }
      })) %>%
      tidyr::unnest(coords)
    
    # Plot transcripts as pointed boxes
    gtf_plot <- ggplot2::ggplot() +
      ggplot2::geom_polygon(
        data = poly_data,
        ggplot2::aes(x = x, y = y, group = transcript_id, fill = strand),
        color = "black") +
      ggplot2::scale_y_continuous(
        limits = c(0.5, max(transcripts$track) + 0.5),# add padding
        breaks = NULL, labels = NULL) +   
      
      ggplot2::scale_fill_manual(values = c("+" = "lightcoral", "-" = "skyblue")) +
      ggplot2::labs(x = "Genomic coordinate", y = "Transcript") +
      ggplot2::theme_classic() +
      ggplot2::theme(
        axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
        axis.title.y = ggplot2::element_blank()
      )
  }
 
    
}
