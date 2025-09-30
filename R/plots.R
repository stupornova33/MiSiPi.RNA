# Plots the coverage over an interval
# @param df a dataframe
# @return a histogram plot

.plot_coverage <- function(df) {
  pos <- count <- NULL
  coverage_hist <- ggplot2::ggplot(df, ggplot2::aes(pos, count)) +
    ggplot2::geom_col(width = 1) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45)) +
    ggplot2::scale_x_continuous(breaks = round(seq(min(df$pos), max(df$pos), by = (max(df$pos) - min(df$pos)) / 5), 5)) +
    ggplot2::labs(x = "Bp position", y = "Depth") +
    ggplot2::theme_minimal()
  return(coverage_hist)
}

# Plot densities of different read sizes
# @param data a data frame
# @param reg_start a whole number
# @param reg_stop a whole number
# @return density

.plot_density <- function(data, reg_start, reg_stop) {
  density <- ggplot2::ggplot(data, ggplot2::aes(x = pos, y = count, group = size)) +
    ggplot2::geom_area(data = subset(data, strand == "pos"), ggplot2::aes(color = size, fill = size)) +
    ggplot2::geom_area(data = subset(data, strand == "neg"), ggplot2::aes(color = size, fill = size)) +
    # ggplot2::scale_x_continuous(limits = c(reg_start - 10, reg_stop + 10), expand = c(0, 0)) +
    # ggplot2::scale_y_continuous(expand = c(0, 0)) +
    # data is organized: pos_26_32, neg_26-32, pos_23_25, neg_23-25, pos_20-22, neg_20-22, pos_18-19, neg_18
    ggplot2::scale_fill_manual(values = c("yellow", "red", "blue", "black")) +
    ggplot2::scale_color_manual(values = c("yellow", "red", "blue", "black")) +
    ggplot2::ggtitle("Read Density") +
    ggplot2::theme_classic() +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5)) +
    # ggplot2::theme(legend.key.size = ggplot2::unit(0.3, 'cm'))+
    ggplot2::theme(axis.text.x = ggplot2::element_text(size = 10, angle = 60, hjust = 1))
  #              axis.ticks.x = ggplot2::element_line(),
  #              axis.line.x = ggplot2::element_line(),
  #              panel.grid = ggplot2::element_blank(),
  #              panel.border = ggplot2::element_blank())
  
  return(density)
}

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
    # filter features that 
    dplyr::filter(seqname == chrom_name, start >= reg_start & start <= reg_stop | end >= reg_start & end <= reg_stop)  %>%
    dplyr::filter(feature != "start_codon" & feature != "CDS" & feature != "stop_codon")
  
  if (nrow(gtf) == 0) {
    message("No features found in region.")
    return(NULL)
  }
  
  # Remove quotation marks
  gtf$attribute <- gsub("\"", "", gtf$attribute)
  
  exons_present <- any(gtf$feature == "exon")
  
  # if transcript and exon information is both present
  if (exons_present) {
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
    genes_transcripts <- gtf %>%
      dplyr::filter(feature == "transcript" | feature == "gene") %>%
      dplyr::distinct()
    
    if (nrow(genes_transcripts) > 30) {
      message("Warning: region contains more features than are plottable. Selecting only first 30. See log file for details.")
      cat(file = logfile, "Region contains more features than are plottable Selecting only first 30. 
     Please consider filtering input GTF file (e.g. select only genes or transcripts or select transcript isoform of interest). \n", append = TRUE)
      genes_transcripts <- genes_transcripts[1:30,]
    }
    
    if (nrow(genes_transcripts) > 0) {
      #Group exons belonging to same transcript and assign transcript index
      plot_df <- gtf %>%
        dplyr::mutate(trans_idx = match(transcript_id, unique(transcript_id))) %>%
        dplyr::select(-c(score, frame, attribute))
      
      # Revert numbered NAs to actual NAs
      plot_df$transcript_id[na_idx] <- NA
      
      
      # genes with no transcript_id (gene-only features)
      genes_only <- plot_df %>%
        dplyr::filter(feature == "gene" & is.na(transcript_id)) %>%
        dplyr::mutate(track = min(transcripts$track, 0) + dplyr::row_number())
      
      genes_only <- genes_only %>%
        dplyr::mutate(start = dplyr::case_when(
          start < reg_start ~ reg_start, # if start is less than cutoff, set it to cutoff
          TRUE ~ start)) %>%
        dplyr::mutate(end = dplyr::case_when(
          end > reg_stop ~ reg_stop, # if stop is greater than cutoff, set it to reg_stop
          TRUE ~ end))
      
      # calculate the x and y coordinates for the exons/transcripts/genes
      transcripts <- plot_df %>%
        dplyr::filter(!is.na(transcript_id)) %>%
        dplyr::distinct(transcript_id, gene_id, strand, .keep_all = TRUE) %>%
        dplyr::arrange(transcript_id) %>%
        dplyr::mutate(track = dplyr::row_number() + nrow(genes_only), midpoint = (end - start)/2) #make genes plot first
      
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
      
      
    } else if (nrow(genes_transcripts) == 0) {
      # exons -> add info on first/last per transcript
      plot_df <- exons
      exons <- plot_df %>%
        dplyr::filter(feature == "exon") %>%
        #dplyr::left_join(transcripts %>% dplyr::select(transcript_id, track, strand), by = "transcript_id") %>%
        dplyr::rename("strand" = "strand.x") %>%
        #dplyr::arrange(track, start) %>%
        dplyr::mutate(track = dplyr::row_number()) %>%
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
      transcripts <- exons
    }
    
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
    
    transcripts <- transcripts %>%
      dplyr::mutate(midpoint = (end + start) / 2)
    
    
    # genes -> make polygons with tip like transcripts
    gene_poly <- genes_only %>%
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
        color = "black"
      ) +
      
      ggplot2::geom_text(
        data = genes_only %>% dplyr::mutate(midpoint = (end + start) / 2),
        ggplot2::aes(x = midpoint, y = track, label = gene_id),
        vjust = 0.5, size = 3
      ) +
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
    
  } else {  #if exon information does not exist, only transcript
    
    
    genes_transcripts <- gtf %>% dplyr::filter(feature == "transcript" | feature == "gene") %>% dplyr::distinct()
    if(nrow(genes_transcripts) > 30){
      message("Warning: region contains more features than are plottable. Selecting only first 30. See log file for details.")
      cat(file = logfile, "Region contains more features than are plottable Selecting only first 30. 
     Please consider filtering input GTF file (e.g. select only genes or transcripts or select transcript isoform of interest). \n", append = TRUE)
      genes_transcripts <- genes_transcripts[1:30,]
      
    }
    
    # Extract transcript_id and gene_id from attribute field
    genes_transcripts$transcript_id <- sub('.*transcript_id "([^"]+)".*', '\\1', genes_transcripts$attribute)
    genes_transcripts$gene_id <- sub('.*gene_id "([^"]+)".*', '\\1', genes_transcripts$attribute)
    genes_transcripts$gene_id <- stringr::str_extract(genes_transcripts$gene_id, "(?<=gene_id )[^;]+")
    
    # Keep only transcripts (no polygon building here)
    transcripts <- genes_transcripts %>%
      dplyr::filter(feature == "transcript") %>%
      dplyr::distinct(transcript_id, .keep_all = TRUE) %>%
      dplyr::mutate(track = dplyr::row_number()) %>%
      dplyr::select(transcript_id, gene_id, start, end, strand, track) %>%
      dplyr::mutate(
        start = ifelse(start < reg_start, reg_start, start),
        end   = ifelse(end > reg_stop, reg_stop, end)
      )
    
    transcripts <- transcripts %>%
      dplyr::mutate(start = dplyr::case_when(
        start < reg_start ~ reg_start, # if start is less than cutoff, set it to cutoff
        TRUE ~ start)) %>%
      dplyr::mutate(end = dplyr::case_when(
        end > reg_stop ~ reg_stop, # if stop is greater than cutoff, set it to reg_stop
        TRUE ~ end))
    
    
    # Parameters
    tip_size <- 0.2
    exon_height <- 0.2
    
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
    
  }
  
  
  
}

# plot heatmaps for si or pi rna
# takes a matrix as input and a specified color palette
# outputs plots
# @param results a data table
# @param method "piRNA" or "siRNA"
# @param pal a string
# @return plots

.plot_heat <- function(results, method = c("piRNA", "siRNA"), pal = c("RdYlBl", "yelOrRed", "MagYel", "Greens", "BlYel")) {
  method <- match.arg(method)
  m_pal <- match.arg(pal)
  m_pal <- switch(m_pal,
                  "RdYlBl" = c("#4575b4", "#74add1", "#abd9e9", "#e0f3f8", "#fee090", "#fdae61", "#f46d43", "#d73027"),
                  "BlYel" = c("#0c2c84", "#225ea8", "#1d91c0", "#41b6c4", "#7fcdbb", "#c7e9b4", "#edf8b1"),
                  "yelOrRed" = c("#FFFFCC", "#FFFF66", "#FFCC00", "#FF9900", "#FF6600"),
                  "MagYel" = c("#330033", "#660066", "#990066", "#FF6600", "#FFCC00", "#FFFF66"),
                  "Greens" = c("#003333", "#006666", "#009966", "#00CC33", "#CCFF99", "#FFFF99", "#FFFF00")
  )
  
  if (method == "piRNA") {
    plot_title <- "Reads With Proper Overlaps By Size"
  } else {
    plot_title <- "Reads With Proper Overhangs By Size"
  }
  
  p <- pheatmap::pheatmap(results, main = plot_title, cluster_cols = FALSE, cluster_rows = FALSE, fontsize = 12, color = m_pal)
  
  # Wrap plot in ggplotify::as.grob since pheatmaps can't be coerced to grob by default
  p <- ggplotify::as.grob(p)
  
  return(p)
}

# plots the arc diagram
# @param filePath a string
# @return an arc plot

.plot_helix <- function(filePath) {
  R4RNA::plotHelix(helix = R4RNA::readHelix(filePath), line = TRUE, arrow = FALSE, lwd = 2.25, scale = FALSE)
  grDevices::dev.control("enable")
  temp <- grDevices::recordPlot()
  return(temp)
}

.plot_siRNA_hp_phasing_probability_combined <- function(plus_phased_table, minus_phased_table) {
  phased_dist <- phased_z <- NULL
  
  plus_color <- "red"
  minus_color <- "blue"
  plot_title <- "siRNA Phasing Probability"
  
  p <- ggplot2::ggplot(plus_phased_table, ggplot2::aes(x = phased_dist, y = phased_z)) +
    ggplot2::geom_line(linewidth = 1.25, color = plus_color) +
    ggplot2::geom_line(data = minus_phased_table, ggplot2::aes(x = phased_dist, y = phased_z), linewidth = 1.25, color = minus_color) +
    ggplot2::scale_x_continuous("3' to 5' Distance", labels = seq(1, 50, by = 5), breaks = seq(1, 50, by = 5)) +
    ggplot2::scale_y_continuous("Z-score") +
    ggplot2::ggtitle(plot_title) +
    ggplot2::theme_classic() +
    ggplot2::theme(axis.text = ggplot2::element_text(size = 12), plot.title = ggplot2::element_text(size = 14, hjust = 0.5)) +
    ggplot2::theme(text = ggplot2::element_text(size = 12))
  
  return(p)
}

# Plots read sizes over large loci
# Faster than geom_area but doesn't look as good
# @param data a dataframe consisting of three columns, "Position", "Count" and "Length".
# @param reg_start The start position of the region of interest.
# @param reg_stop The end position of a region of interest.
# @return Density Plot

.plot_large_density <- function(data, reg_start, reg_stop) {
  density <- ggplot2::ggplot(data, ggplot2::aes(x = pos, y = count, fill = size)) +
    ggplot2::geom_col() +
    ggplot2::scale_fill_manual(values = c("yellow", "red", "blue", "black", "yellow", "red", "blue", "black")) +
    ggplot2::ggtitle("Read Density") +
    ggplot2::theme_classic() +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5)) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(size = 10, angle = 60, hjust = 1))
  
  return(density)
}

# function to plot the zscore for Dicer overhangs
# @param overhang1 a dataframe
# @param strand a string, "+", "-", or "none"
# @return plot

.plot_overhangz <- function(overhang1, strand) {
  # calculates Dicer signature in pairs
  
  # shift <- value <- proper_count <- NULL
  shift <- c(-4, -3, -2, -1, 0, 1, 2, 3, 4)
  
  if (strand == "+") {
    col <- "red"
    title <- "Dicer Overhang Prob (plus strand)"
  } else if (strand == "-") {
    col <- "blue"
    title <- "Dicer Overhang Prob (minus strand)"
  } else {
    col <- "black"
    title <- "siRNA Dicer Overhang Probability (dual strand)"
  }
  
  p <- ggplot2::ggplot(overhang1, ggplot2::aes(x = shift, y = zscore)) +
    ggplot2::geom_line(color = col, linewidth = 2) +
    ggplot2::ggtitle(title) +
    ggplot2::scale_x_continuous("Shift") +
    ggplot2::scale_y_continuous("Z-score") +
    ggplot2::theme_classic() +
    ggplot2::theme(text = ggplot2::element_text(size = 15), plot.title = ggplot2::element_text(hjust = 0.5, size = 14)) +
    ggplot2::theme(plot.margin = ggplot2::unit(c(2, 0, 0, 0), "cm"))
  return(p)
}

.plot_siRNA_overhangs_combined <- function(plus_overhangs, minus_overhangs, dual_strand_overhangs) {
  shift <- c(-4, -3, -2, -1, 0, 1, 2, 3, 4)
  
  plus_color <- "red"
  minus_color <- "blue" 
  dual_color <- "black"
  plot_title <- "siRNA Dicer Overhang Probability"
  
  
  
  p <- ggplot2::ggplot(plus_overhangs, ggplot2::aes(x = shift, y = zscore)) +
    ggplot2::geom_line(color = plus_color, linewidth = 2) +
    ggplot2::geom_line(data = minus_overhangs, ggplot2::aes(x = shift, y = zscore), color = minus_color, linewidth = 2) +
    ggplot2::geom_line(data = dual_strand_overhangs, ggplot2::aes(x = shift, y = zscore), color = dual_color, linewidth = 2) +
    ggplot2::ggtitle(plot_title) +
    ggplot2::scale_x_continuous("Shift") +
    ggplot2::scale_y_continuous("Z-score") +
    ggplot2::theme_classic() +
    ggplot2::theme(text = ggplot2::element_text(size = 15), plot.title = ggplot2::element_text(hjust = 0.5, size = 14)) +
    ggplot2::theme(plot.margin = ggplot2::unit(c(2, 0, 0, 0), "cm"))
  return(p) 
  
}

.plot_miRNA_dicer_overhang_probability <- function(plus_overhangs, minus_overhangs) {
  # shift <- value <- proper_count <- NULL
  shift <- c(-4, -3, -2, -1, 0, 1, 2, 3, 4)
  
  p_overhangs <- data.frame(
    shift = plus_overhangs$shift,
    zscore = plus_overhangs$zscore,
    Strand = "Plus"
  )
  
  m_overhangs <- data.frame(
    shift = minus_overhangs$shift,
    zscore = minus_overhangs$zscore,
    Strand = "Minus"
  )
  
  all_overhangs <- dplyr::bind_rows(p_overhangs, m_overhangs)
  
  p <- ggplot2::ggplot(all_overhangs, ggplot2::aes(x = shift, y = zscore, color = Strand)) +
    ggplot2::geom_line(linewidth = 2) +
    ggplot2::scale_x_continuous("Shift") +
    ggplot2::scale_y_continuous("Z-score") +
    ggplot2::scale_color_manual(values = c("Plus" = "red", "Minus" = "blue")) +
    ggplot2::ggtitle("miRNA Dicer Overhang Probability") +
    ggplot2::theme_classic() +
    ggplot2::theme(text = ggplot2::element_text(size = 15), plot.title = ggplot2::element_text(hjust = 0.5, size = 14)) +
    ggplot2::theme(plot.margin = ggplot2::unit(c(2, 0, 0, 0), "cm"))
  return(p)
}

.plot_siRNA_overhangs_combined <- function(plus_overhangs, minus_overhangs, dual_strand_overhangs) {
  shift <- c(-4, -3, -2, -1, 0, 1, 2, 3, 4)
  
  plus_color <- "red"
  minus_color <- "blue"
  dual_color <- "black"
  plot_title <- "siRNA Dicer Overhang Probability"
  
  
  p <- ggplot2::ggplot(plus_overhangs, ggplot2::aes(x = shift, y = zscore)) +
    ggplot2::geom_line(color = plus_color, linewidth = 2) +
    ggplot2::geom_line(data = minus_overhangs, ggplot2::aes(x = shift, y = zscore), color = minus_color, linewidth = 2) +
    ggplot2::geom_line(data = dual_strand_overhangs, ggplot2::aes(x = shift, y = zscore), color = dual_color, linewidth = 2) +
    ggplot2::ggtitle(plot_title) +
    ggplot2::scale_x_continuous("Shift") +
    ggplot2::scale_y_continuous("Z-score") +
    ggplot2::theme_classic() +
    ggplot2::theme(text = ggplot2::element_text(size = 15), plot.title = ggplot2::element_text(hjust = 0.5, size = 14)) +
    ggplot2::theme(plot.margin = ggplot2::unit(c(2, 0, 0, 0), "cm"))
  
  return(p)
}

# plot the overlap probability
# @param z_df a data frame
# @return a ggplot object

# Previously named .plot_overlapz
.plot_piRNA_overlap_probability <- function(z_df) {
  Z_score <- Overlap <- NULL
  p <- ggplot2::ggplot(z_df, ggplot2::aes(x = Overlap, y = zscore)) +
    ggplot2::geom_line(color = "black", linewidth = 1.5) +
    ggplot2::theme_classic() +
    ggplot2::scale_x_continuous("Overlap", breaks = seq(0, 32, 5), labels = seq(0, 32, 5)) +
    ggplot2::theme(axis.title.x = ggplot2::element_text(size = 14), axis.title.y = ggplot2::element_text(size = 14)) +
    ggplot2::theme(plot.title = ggplot2::element_text(size = 14, hjust = 0.5)) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(size = 12), axis.text.y = ggplot2::element_text(size = 12)) +
    ggplot2::ggtitle("piRNA Overlap Probability")
  
  return(p)
}

.plot_miRNA_overlap_probability <- function(plus_z_df, minus_z_df) {
  zscore <- Overlap <- NULL
  
  
  p_overlaps <- data.frame(
    Overlap = plus_z_df$Overlap,
    zscore = plus_z_df$zscore,
    Strand = "Plus"
  )
  
  m_overlaps <- data.frame(
    Overlap = minus_z_df$Overlap,
    zscore = minus_z_df$zscore,
    Strand = "Minus"
  )
  
  all_overlaps <- dplyr::bind_rows(p_overlaps, m_overlaps)
  
  p <- ggplot2::ggplot(all_overlaps, ggplot2::aes(x = Overlap, y = zscore, color = Strand)) +
    ggplot2::geom_line(linewidth = 1.5) +
    ggplot2::scale_x_continuous("Overlap") +
    ggplot2::scale_y_continuous("Z-score") +
    ggplot2::scale_color_manual(values = c("Plus" = "red", "Minus" = "blue")) +
    ggplot2::ggtitle("miRNA Overlap Probability") +
    ggplot2::theme_classic() +
    ggplot2::theme(text = ggplot2::element_text(size = 15), plot.title = ggplot2::element_text(hjust = 0.5, size = 14)) +
    ggplot2::theme(plot.margin = ggplot2::unit(c(2, 0, 0, 0), "cm"))
  
  return(p)
}

# phasing function
# processes reads according to phased piRNA algorithm
# plots output
# takes data frame
# returns plots
#
# @param df1 a data frame
# @param strand a string, "+" or "-"
# @return plots

.plot_phasedz <- function(df1, strand = NULL) {
  phased_num <- phased_dist <- value <- phased_z <- phased26_dist <- phased26_z <- phased_dist1 <- phased_z1 <- phased_dist2 <- phased_z2 <- NULL
  
  if (!is.null(strand)) {
    if (strand == "+") {
      title <- "Phasing Prob (Plus Strand)"
      col <- "red"
    } else {
      title <- "Phasing Prob (Minus Strand)"
      col <- "blue"
    }
    
    if (ncol(df1) > 3 & colnames(df1)[1] != "phased_dist1") {
      p <- ggplot2::ggplot(df1, ggplot2::aes(x = phased_dist)) +
        ggplot2::scale_y_continuous("Z-score") +
        ggplot2::ggtitle(title) +
        ggplot2::scale_x_continuous("3' to 5' Distance", labels = seq(1, 65, by = 5), breaks = seq(1, 65, by = 5)) +
        ggplot2::geom_line(ggplot2::aes(y = phased_z), linewidth = 1.25, color = "black") +
        ggplot2::geom_line(ggplot2::aes(x = phased26_dist, y = phased26_z), linewidth = 1.5, color = col) +
        ggplot2::theme_classic() +
        ggplot2::theme(axis.text = ggplot2::element_text(size = 12), plot.title = ggplot2::element_text(size = 14, hjust = 0.5)) +
        ggplot2::theme(axis.text.x = ggplot2::element_text(size = 12, angle = 45)) +
        ggplot2::theme(text = ggplot2::element_text(size = 12))
    } else if (ncol(df1) > 3 & colnames(df1)[1] == "phased_dist1") {
      p <- ggplot2::ggplot(df1, ggplot2::aes(x = phased_dist1)) +
        ggplot2::scale_y_continuous("Z-score") +
        ggplot2::ggtitle(title) +
        ggplot2::scale_x_continuous("3' to 5' Distance", labels = seq(1, 65, by = 5), breaks = seq(1, 65, by = 5)) +
        ggplot2::geom_line(ggplot2::aes(y = phased_z1), linewidth = 1.25, color = "blue") +
        ggplot2::geom_line(ggplot2::aes(x = phased_dist2, y = phased_z2), linewidth = 1.5, color = "red") +
        ggplot2::theme_classic() +
        ggplot2::theme(axis.text = ggplot2::element_text(size = 12), plot.title = ggplot2::element_text(size = 14, hjust = 0.5)) +
        ggplot2::theme(axis.text.x = ggplot2::element_text(size = 12, angle = 45)) +
        ggplot2::theme(text = ggplot2::element_text(size = 12))
    } else {
      p <- ggplot2::ggplot(df1, ggplot2::aes(x = phased_dist)) +
        ggplot2::scale_y_continuous("Z-score") +
        ggplot2::ggtitle(title) +
        ggplot2::scale_x_continuous("3' to 5' Distance", labels = seq(1, 65, by = 5), breaks = seq(1, 65, by = 5)) +
        ggplot2::geom_line(ggplot2::aes(y = phased_z), linewidth = 1.25, color = col) +
        ggplot2::theme_classic() +
        ggplot2::theme(axis.text = ggplot2::element_text(size = 12), plot.title = ggplot2::element_text(size = 14, hjust = 0.5)) +
        ggplot2::theme(axis.text.x = ggplot2::element_text(size = 12, angle = 45)) +
        ggplot2::theme(text = ggplot2::element_text(size = 12))
    }
  }
  
  return(p)
}

.plot_piRNA_phasing_probability_combined <- function(plus_phased_table, minus_phased_table) {
  phased_num <- phased_dist <- value <- phased_z <- phased26_dist <- phased26_z <- phased_dist1 <- phased_z1 <- phased_dist2 <- phased_z2 <- NULL
  
  # Convert data frames to long and combine into a single data frame for plotting
  plus_df <- data.frame(
    phased_dist = plus_phased_table$phased_dist,
    value = plus_phased_table$phased_z,
    Type = "(+) All Sizes"
  )
  plus_df2 <- data.frame(
    phased_dist = plus_phased_table$phased_dist,
    value = plus_phased_table$phased26_z,
    Type = "(+) >= 26nt"
  )
  minus_df <- data.frame(
    phased_dist = minus_phased_table$phased_dist,
    value = minus_phased_table$phased_z,
    Type = "(-) All Sizes"
  )
  minus_df2 <- data.frame(
    phased_dist = minus_phased_table$phased_dist,
    value = minus_phased_table$phased26_z,
    Type = "(-) >= 26nt"
  )
  
  # Combine all into one data frame
  plot_df <- rbind(plus_df, plus_df2, minus_df, minus_df2)
  
  p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = phased_dist, y = value, color = Type)) +
    ggplot2::geom_line(linewidth = 1, alpha = 0.7) +
    ggplot2::scale_color_manual(values = c("(+) All Sizes" = "red", "(+) >= 26nt" = "red4", "(-) All Sizes" = "blue", "(-) >= 26nt" = "blue4")) +
    ggplot2::scale_y_continuous("Z-score") +
    ggplot2::scale_x_continuous("3' to 5' Distance", labels = seq(1, 65, by = 5), breaks = seq(1, 65, by = 5)) +
    ggplot2::ggtitle("piRNA Phasing Probability") +
    ggplot2::theme_classic() +
    ggplot2::theme(
      axis.text = ggplot2::element_text(size = 12),
      axis.text.x = ggplot2::element_text(size = 12, angle = 45),
      plot.title = ggplot2::element_text(size = 14, hjust = 0.5),
      text = ggplot2::element_text(size = 12)
    )
  
  return(p)
}


# Write the plots to a file
.print_miRNA_plots <- function(read_distribution_plot, read_density_plot, dicer_overhang_plot, overlap_probability_plot, out_type, prefix, wkdir, plot_details) {
  
  plot_title <- cowplot::ggdraw() +
    cowplot::draw_label(
      plot_details$title,
      x = 0.5,
      hjust = 0.5,
      size = 12
    ) +
    ggplot2::theme(plot.margin = ggplot2::margin(7, 0, 0, 0))
  
  plot_caption <- cowplot::ggdraw() +
    cowplot::draw_label(
      plot_details$caption,
      x = 0.5,
      hjust = 0.5,
      size = 12
    ) +
    ggplot2::theme(plot.margin = ggplot2::margin(0, 0, 5, 0))
  
  plot_body <- cowplot::plot_grid(
    read_distribution_plot, read_density_plot,
    dicer_overhang_plot, overlap_probability_plot,
    ncol = 2,
    align = "hv",
    axis = "lrtb"
  )
  
  all_plot <- cowplot::plot_grid(
    plot_title,
    plot_body,
    plot_caption,
    ncol = 1,
    rel_heights = c(0.07, 1, 0.07)
  )
  
  if (out_type == "png") {
    grDevices::png(file = file.path(wkdir, paste(prefix, "combined.png", sep = "_")), height = 8, width = 11, units = "in", res = 300)
  } else {
    grDevices::pdf(file = file.path(wkdir, paste(prefix, "combined.pdf", sep = "_")), height = 8, width = 11)
  }
  print(all_plot)
  grDevices::dev.off()
  return(NULL)
}

.plot_siRNA <- function(read_distribution_plot, density_plot, phasedz_plot, overhang_probability_plot, arc_plot, gtf_plot, heat_plot, out_type, prefix, wkdir, plot_details) {
  plot_title <- cowplot::ggdraw() +
    cowplot::draw_label(
      plot_details$title,
      x = 0.5,
      hjust = 0.5,
      size = 12
    ) +
    ggplot2::theme(plot.margin = ggplot2::margin(7, 0, 0, 0))
  
  plot_caption <- cowplot::ggdraw() +
    cowplot::draw_label(
      plot_details$caption,
      x = 0.5,
      hjust = 0.5,
      size = 12
    ) +
    ggplot2::theme(plot.margin = ggplot2::margin(0, 0, 5, 0))
  
  plot_body <- cowplot::plot_grid(
    read_distribution_plot, arc_plot,
    overhang_probability_plot, density_plot,
    phasedz_plot, gtf_plot,
    heat_plot,
    ncol = 2,
    #align = "hv",
    axis = "lrtb",
    rel_widths = c(1, 1,0.8,1,0.8,1,0.8)
    
  )
  
  all_plot <- cowplot::plot_grid(
    plot_title,
    plot_body,
    plot_caption,
    ncol = 1,
    rel_heights = c(0.07, 1, 0.07)
  )
  
  # TODO: Change dimensions to better suit new layout
  if (out_type == "png") {
    grDevices::png(file = file.path(wkdir, paste0(prefix, "_si_plot.png")), height = 20, width = 14, units = "in", res = 300)
  } else {
    grDevices::pdf(file = file.path(wkdir, paste0(prefix, "_si_plot.pdf")), height = 20, width = 14)
  }
  
  print(all_plot)
  grDevices::dev.off()
  return(NULL)
}

plot_piRNA <- function(read_distribution_plot, density_plot, overlap_probability_plot, phased_probability_plot, heat_plot, out_type, prefix, wkdir, plot_details) {
  # Modify density_plot to display better in this layout
  density_plot <- density_plot +
    ggplot2::theme(axis.text.x = ggplot2::element_text(size = 12, angle = 60, hjust = 1),
                   axis.text.y = ggplot2::element_text(size = 12),
                   axis.title.x = ggplot2::element_text(size = 12),
                   axis.title.y = ggplot2::element_text(size = 12))
  
  
  
  plot_title <- cowplot::ggdraw() +
    cowplot::draw_label(
      plot_details$title,
      x = 0.5,
      hjust = 0.5,
      size = 12
    ) +
    ggplot2::theme(plot.margin = ggplot2::margin(7, 0, 0, 0))
  
  plot_caption <- cowplot::ggdraw() +
    cowplot::draw_label(
      plot_details$caption,
      x = 0.5,
      hjust = 0.5,
      size = 12
    ) +
    ggplot2::theme(plot.margin = ggplot2::margin(0, 0, 5, 0))
  
  plot_body_top <- cowplot::plot_grid(
    read_distribution_plot, heat_plot,
    overlap_probability_plot, phased_probability_plot,
    ncol = 2,
    align = "hv",
    axis = "lrtb"
  )
  
  plot_body_bottom <- cowplot::plot_grid(
    NULL, density_plot, NULL,
    ncol = 3,
    rel_widths = c(0.1, 1, 0.1),
    align = "hv",
    axis = "lrtb"
  )
  
  all_plot <- cowplot::plot_grid(
    plot_title,
    plot_body_top,
    plot_body_bottom,
    plot_caption,
    ncol = 1,
    rel_heights = c(0.25, 2, 0.8, 0.15)
  )
  
  if (out_type == "png" || out_type == "PNG") {
    grDevices::png(file = file.path(wkdir, paste0(prefix, "_pi-zscore.png")), height = 17, width = 14, units = "in", res = 300)
  } else {
    grDevices::cairo_pdf(file = file.path(wkdir, paste0(prefix, "_pi-zscore.pdf")), height = 18, width = 14)
  }
  
  print(all_plot)
  grDevices::dev.off()
}

# Arrange the plots from .run_all() and write the combined plot to a file
plot_combined_plots <- function(p, out_type, prefix, output_dir, plot_details) {
  
  plot_title <- cowplot::ggdraw() +
    cowplot::draw_label(
      plot_details$title,
      x = 0.5,
      hjust = 0.5,
      size = 12
    ) +
    ggplot2::theme(plot.margin = ggplot2::margin(7, 0, 0, 0))
  
  plot_caption <- cowplot::ggdraw() +
    cowplot::draw_label(
      plot_details$caption,
      x = 0.5,
      hjust = 0.5,
      size = 12
    ) +
    ggplot2::theme(plot.margin = ggplot2::margin(0, 0, 5, 0))
  
  plot_body <- cowplot::plot_grid(
    p$read_distribution_plot, p$siRNA_arc_plot, p$siRNA_dicer_overhang_probability_plot, p$piRNA_overlap_probability_plot,
    p$miRNA_dicer_overhang_plot, p$read_density_plot, p$siRNA_phasing_probability_plot, p$piRNA_phasing_probability_plot,
    p$miRNA_overlap_probability_plot, p$siRNA_gtf_plot, p$siRNA_proper_overhangs_by_size_plot, p$piRNA_proper_overlaps_by_size_plot,

    ncol = 4,
    align = "hv",
    axis = "lrtb"
  )
  
  all_plot <- cowplot::plot_grid(
    plot_title,
    plot_body,
    plot_caption,
    ncol = 1,
    rel_heights = c(0.1, 1, 0.1)
  )
  
  if (out_type == "png") {
    grDevices::png(file = file.path(output_dir, "combined_plots", paste(prefix, "combined.png", sep = "_")), height = 15, width = 26, units = "in", res = 300)
  } else {
    grDevices::pdf(file = file.path(output_dir, "combined_plots", paste(prefix, "combined.pdf", sep = "_")), height = 15, width = 26)
  }
  print(all_plot)
  grDevices::dev.off()
}



plot_title <- function(bam_file, bed_file, genome_file, prefix, i) {
  now <- format(lubridate::now(), "%Y-%m-%d %H:%M:%S")
  
  misipi_version <- packageVersion("MiSiPi.RNA")
  
  iteration_str <- paste0(i, ")")
  
  # Drop full path of file and keep basename
  bam_file <- basename(bam_file)
  bed_file <- basename(bed_file)
  genome_file <- basename(genome_file)
  
  p_title <- paste("MiSiPi Results for locus:", prefix, "(Bed file line:", iteration_str)
  p_subtitle <- paste("Bam:", bam_file, "| Bed:", bed_file, "| Genome:", genome_file)
  p_caption <- paste("Run at:", now, "with MiSiPi.RNA Version:", misipi_version)
  
  plot_details <- list()
  plot_details$title <- paste0(p_title, "\n", p_subtitle)
  #plot_details$subtitle <- p_subtitle
  plot_details$caption <- p_caption
  return(plot_details)
}

null_plot <- function(type, reason) {
  p <- ggplot2::ggplot() +
    ggplot2::annotate("text", x = 0.5, y = 0.5, label = reason, size = 5) +
    ggplot2::ggtitle(type) +
    ggplot2::theme_void() +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))
  
  return(p)
}

# .plot_sizes_by_strand
# plots the size distribution of reads separately by strand
# @param stranded_size_dist
# @return A plot

.plot_sizes_by_strand <- function(stranded_size_dist) {
  options(scipen = 999)
  
  p_size_dist <- stranded_size_dist %>% dplyr::filter(strand == "plus")
  m_size_dist <- stranded_size_dist %>% dplyr::filter(strand == "minus")
  
  m_size_dist$count <- m_size_dist$count * -1
  
  empty <- data.frame(width = c(seq(18, 32, by = 1)), strand = "plus", count = numeric(15L))
  p_dist <- merge(empty, p_size_dist, by = "width", all.x = TRUE) %>%
    dplyr::select(-c(count.x, strand.y)) %>%
    dplyr::rename("total_count" = "count.y")
  p_dist[is.na(p_dist)] <- 0
  
  empty <- data.frame(width = c(seq(18, 32, by = 1)), strand = "minus", count = numeric(15L))
  m_dist <- merge(empty, m_size_dist, by = "width", all.x = TRUE) %>%
    dplyr::select(-c(count.x, strand.y)) %>%
    dplyr::rename("total_count" = "count.y")
  m_dist[is.na(m_dist)] <- 0
  m_dist$total_count <- m_dist$total_count * -1
  
  df <- rbind(p_dist, m_dist)
  
  df$first <- sub("T", "U", df$first)
  
  df$first <- factor(df$first, levels = c("A", "C", "G", "U", "N"))
  
  letter_levels <- df %>% dplyr::select(first) %>% dplyr::distinct()
  if(nrow(letter_levels) == 4){
    pal <- c( "darkgreen", "red","blue", "yellow")
  } else {
    pal <- c( "darkgreen", "red","blue", "darkgrey","yellow")
  }
  
  p <- ggplot2::ggplot(df, ggplot2::aes(x = factor(width), y = total_count, fill = first)) +
    ggplot2::geom_bar(stat = "identity", position = "stack") +
    ggplot2::scale_y_continuous(labels = abs) +
    ggplot2::scale_fill_manual(values = c( "darkgreen", "red","blue", "yellow", "darkgrey")) +
    ggplot2::labs(x = "Read size", y = "Read count", title = "Read Distribution") +
    ggplot2::theme_minimal() +
    ggplot2::theme(axis.text=ggplot2::element_text(size=14),
                   plot.title = ggplot2::element_text(size = 14, hjust = 0.5))
  
  
  return(p)
}

# plots the size distribution of reads
#
# @param dist a data table
# @return plots

.plot_sizes <- function(dist) {
  options(scipen = 999)
  
  num_nt <- dist %>%
    dplyr::select(c(first)) %>%
    dplyr::distinct()
  
  if (nrow(num_nt) == 4) {
    man_pal <- c("green4", "red3", "blue4", "yellow3")
  } else {
    man_pal <- c("green4", "red3", "blue4", "grey", "yellow3")
  }
  
  dist$first <- sub("T", "U", dist$first)
  
  width <- count <- first <- NULL
  p <- ggplot2::ggplot(dist, ggplot2::aes(x = width, y = count, fill = first)) +
    ggplot2::geom_bar(stat = "identity") +
    ggplot2::theme_classic() +
    ggplot2::labs(title = "Read Size Distribution") +
    ggplot2::scale_x_continuous(breaks = seq(18, 32, 2)) +
    ggplot2::scale_fill_manual(values = man_pal) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(size = 16, hjust = 0.5),
      # plot.subtitle = ggplot2::element_text(size = 12, hjust = 0.5),
      axis.text.x = ggplot2::element_text(size = 12, angle = 45, vjust = 0.5),
      axis.text.y = ggplot2::element_text(size = 12),
      axis.title.y = ggplot2::element_text(size = 12, margin = ggplot2::margin(r = 5)),
      axis.title.x = ggplot2::element_text(size = 12)
    )
  
  return(p)
}

# plots the dot bracket structure above the nucleotide sequence
# @param vienna a string
# @param sequence a string
# @return a plot

.plot_text <- function(vienna, sequence) {
  x <- y <- text <- NULL
  v <- as.vector(unlist(strsplit(vienna, "")))
  c <- as.vector(unlist(strsplit(sequence, "")))
  
  v_df <- data.frame("x" = c(1:length(v)), "y" = c(5), "text" = v)
  c_df <- data.frame("x" = c(1:length(c)), "y" = c(4.6), "text" = c)
  
  f_df <- rbind(v_df, c_df)
  
  g <- ggplot2::ggplot(f_df, ggplot2::aes(x, y, label = text)) +
    ggplot2::geom_text() +
    ggplot2::xlim(0, max(length(v), length(c) + 2)) +
    ggplot2::scale_y_continuous(breaks = seq(2, 5, by = 1), limits = c(2, 5)) +
    # ggbreak::scale_y_cut(breaks=c(4), which=c(1), scales=c(3))+
    ggplot2::theme(plot.margin = ggplot2::unit(c(0, 0, 0, 0), "cm")) +
    ggplot2::theme_void()
  
  return(g)
}
