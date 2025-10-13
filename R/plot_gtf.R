# takes a gtf file
# outputs a plot
# @param gtf_file a string path leading to the gtf file
# @param chrom_name a string
# @param reg_start an integer
# @param reg_stop an integer
# @param logfile
#
# @return plot

# Wrapper function for gtf sub functions called based on combination of features present
.plot_gtf <- function(gtf_file, chrom_name, reg_start, reg_stop, logfile) {
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
  transcripts_present <- any(gtf$feature == "transcript")
  genes_present <- any(gtf$feature == "gene")
  
  if (exons_present & !transcripts_present & !genes_present) { # exons only
    gtf_plot <- .plot_exons_only(gtf, reg_start, reg_stop, logfile)
  } else if (exons_present & transcripts_present & !genes_present) { # exons and transcripts
    gtf_plot <- .plot_transcripts_exons(gtf, reg_start, reg_stop, logfile)
  } else if (!exons_present & transcripts_present & genes_present) { # genes and transcripts
    gtf_plot <- .plot_genes_transcripts(gtf, reg_start, reg_stop, logfile)
  } else if (exons_present & transcripts_present & genes_present) { # exons and transcripts and genes
    gtf_plot <- .plot_genes_exon_transcripts(gtf, reg_start, reg_stop, logfile)
  } else if (!exons_present & transcripts_present & !genes_present) { # transcripts only
    gtf_plot <- .plot_transcripts_only(gtf, reg_start, reg_stop, logfile)
  } else if (!exons_present & !transcripts_present & genes_present) { # genes only
    gtf_plot <- .plot_genes_only(gtf, reg_start, reg_stop, logfile)
  } else { # exons and genes only
    gtf_plot <- .plot_genes_exons(gtf, reg_start, reg_stop, logfile)
  }
  
  return(gtf_plot)
}

# Plots annotated features over region of interest
# @param gtf A nine column annotation file
# @param reg_start an integer specifying the start location of region of interest
# @param reg_stop an integer specifying the end location of region of interest
# @return gtf_plot

.plot_exons_only <- function(gtf, reg_start, reg_stop, logfile) {
  
  roi_length <- reg_stop - reg_start + 1
  
  gtf <- gtf %>%
    dplyr::mutate(
      gene_id = .get_gene_id(attribute),
      transcript_id = .get_transcript_id(attribute)
    )
  
  #### SEPARATE STRANDS ####
  sense_df <- gtf %>%
    dplyr::filter(strand == "+")
  antisense_df <- gtf %>%
    dplyr::filter(strand == "-")
  nonsense_df <- gtf %>%
    dplyr::filter(strand == ".")
  
  #### SORT ####
  # Arrange the data frames for proper plotting
  arrange_rows <- function(plot_df) {
    arranged_df <- plot_df %>%
      dplyr::arrange(transcript_id, start, end)
    
    return(arranged_df)
  }
  
  sense_df <- arrange_rows(sense_df)
  antisense_df <- arrange_rows(antisense_df)
  nonsense_df <- arrange_rows(nonsense_df)
  
  #### COUNT ####
  # Count the number of plot rows for each strand
  
  num_s_rows <- nrow(
    sense_df %>%
      dplyr::select(transcript_id) %>%
      dplyr::distinct()
  )
  
  num_a_rows <- nrow(
    antisense_df %>%
      dplyr::select(transcript_id) %>%
      dplyr::distinct()
  )
  
  num_n_rows <- nrow(
    nonsense_df %>%
      dplyr::select(transcript_id) %>%
      dplyr::distinct()
  )
  
  #### CALCULATE PLOT ROW ORDER ####
  
  # Determine row order for plotting
  sense_df <- sense_df %>%
    dplyr::mutate(
      first_exon = feature == "exon" & !duplicated(ifelse(feature == "exon", transcript_id, NA)),
      inc = first_exon,
      row_index = cumsum(inc)
    ) %>%
    dplyr::select(-first_exon, -inc)
  
  antisense_df <- antisense_df %>%
    dplyr::mutate(
      first_exon = feature == "exon" & !duplicated(ifelse(feature == "exon", transcript_id, NA)),
      inc = first_exon,
      row_index = cumsum(inc)
    ) %>%
    dplyr::select(-first_exon, -inc)
  
  nonsense_df <- nonsense_df %>%
    dplyr::mutate(
      first_exon = feature == "exon" & !duplicated(ifelse(feature == "exon", transcript_id, NA)),
      inc = first_exon,
      row_index = cumsum(inc)
    ) %>%
    dplyr::select(-first_exon, -inc)
  
  #### SUBSET ####
  SUBSET <- FALSE
  
  # We can't really plot more than 18-20 rows without things starting to look bad, so we're limiting how many features here
  # 18 was chosen due to being divisible by 2 and 3 for cases where there are numerous sense, antisense, and nonstranded features present
  max_plot_rows <- 18
  num_plot_rows <- sum(num_s_rows, num_a_rows, num_n_rows)
  
  if (num_plot_rows > 18) {
    SUBSET <- TRUE
    warning_msg <- glue::glue("Warning: This region contains more features than are plottable [{num_plot_rows}]. Selecting first 18.\n")
    cat(file = logfile, warning_msg, append = TRUE)
    
    s_ratio <- dplyr::if_else(
      num_s_rows >= max_plot_rows,
      1,                           # TRUE
      num_s_rows / max_plot_rows   # FALSE
    )
    
    a_ratio <- dplyr::if_else(
      num_a_rows >= max_plot_rows,
      1,                           # TRUE
      num_a_rows / max_plot_rows   # FALSE
    )
    
    n_ratio <- dplyr::if_else(
      num_n_rows >= max_plot_rows,
      1,                           # TRUE
      num_n_rows / max_plot_rows   # FALSE
    )
    
    ratio_sum <- sum(s_ratio, a_ratio, n_ratio)
    
    s_factor <- s_ratio / ratio_sum
    a_factor <- a_ratio / ratio_sum
    n_factor <- n_ratio / ratio_sum
    
    # The number of rows contributed by each strands observations should be 18 +- 1 depending on rounding
    s_rows <- round(max_plot_rows * s_factor)
    a_rows <- round(max_plot_rows * a_factor)
    n_rows <- round(max_plot_rows * n_factor)
    
    sense_df <- sense_df %>%
      dplyr::filter(row_index <= s_rows)
    
    antisense_df <- antisense_df %>%
      dplyr::filter(row_index <= a_rows)
    
    nonsense_df <- nonsense_df %>%
      dplyr::filter(row_index <= n_rows)
  }
  
  #### COMBINE ####
  # Modify the row_index of antisense to start after sense and nonsense to start after antisense
  if (nrow(sense_df) == 0) {
    max_sense_row <- 0
  } else {
    max_sense_row <- max(sense_df$row_index)
  }
  
  if (nrow(antisense_df) == 0) {
    max_antisense_row <- 0
  } else {
    max_antisense_row <- max(antisense_df$row_index)
  }
  
  antisense_df$row_index <- antisense_df$row_index + max_sense_row
  nonsense_df$row_index <- nonsense_df$row_index + max_sense_row + max_antisense_row
  
  plot_df <- dplyr::bind_rows(sense_df, antisense_df, nonsense_df)
  
  # Add in a unique index for each row to allow for polygon grouping
  plot_df <- plot_df %>%
    dplyr::mutate(polygon_idx = dplyr::row_number())
  
  # Add a column to track if a feature's positions were truncated on a pointed end
  # We're not adding the arrow if the feature extends past the region of interest
  # Also not adding arrows if no strand info is present, but we'll handle that elsewhere
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
  
  # Calculating this again in case the data frame was subset
  num_plot_rows <- max(plot_df$row_index)
  
  if (num_plot_rows < 8) {
    TXT_SIZE <- 6
  } else if (num_plot_rows < 12) {
    TXT_SIZE <- 5
  } else {
    TXT_SIZE <- 4
  }
  
  if (SUBSET) {
    plot_title <- "Features filtered to fit on plot"
  } else {
    plot_title <- ""
  }
  
  # Generate Plot Object
  gtf_plot <- .get_gtf_plot(
    plot_title = plot_title,
    TXT_SIZE = TXT_SIZE,
    reg_start = reg_start,
    reg_stop = reg_stop,
    exons = exons,
    exon_poly_data = exon_poly_data,
    exon_labels = exon_labels
  )
  
  return(gtf_plot)
}
  
# Plots the coverage over an interval
# @param df a dataframe
# @return a histogram plot
.plot_transcripts_only <- function(gtf, reg_start, reg_stop, logfile) {
  
  roi_length <- reg_stop - reg_start + 1
  
  gtf <- gtf %>%
    dplyr::mutate(
      gene_id = .get_gene_id(attribute),
      transcript_id = .get_transcript_id(attribute)
    )
  
  #### SEPARATE STRANDS ####
  sense_df <- gtf %>%
    dplyr::filter(strand == "+")
  antisense_df <- gtf %>%
    dplyr::filter(strand == "-")
  nonsense_df <- gtf %>%
    dplyr::filter(strand == ".")
  
  #### SORT ####
  # Arrange the data frames for proper plotting
  arrange_rows <- function(plot_df) {
    arranged_df <- plot_df %>%
      dplyr::arrange(start, end)
    
    return(arranged_df)
  }
  
  sense_df <- arrange_rows(sense_df)
  antisense_df <- arrange_rows(antisense_df)
  nonsense_df <- arrange_rows(nonsense_df)
  
  #### COUNT ####
  # Count the number of plot rows for each strand
  
  num_s_rows <- nrow(sense_df)
  num_a_rows <- nrow(antisense_df)
  num_n_rows <- nrow(nonsense_df)
  
  #### CALCULATE PLOT ROW ORDER ####
  # Since we just have transcript features in this data frame
  #   the plot order is simply the row order
  sense_df <- sense_df %>%
    dplyr::mutate(row_index = dplyr::row_number())
  antisense_df <- antisense_df %>%
    dplyr::mutate(row_index = dplyr::row_number())
  nonsense_df <- nonsense_df %>%
    dplyr::mutate(row_index = dplyr::row_number())
  
  #### SUBSET ####
  SUBSET <- FALSE
  
  # We can't really plot more than 18-20 rows without things starting to look bad, so we're limiting how many features here
  # 18 was chosen due to being divisible by 2 and 3 for cases where there are numerous sense, antisense, and nonstranded features present
  max_plot_rows <- 18
  num_plot_rows <- sum(num_s_rows, num_a_rows, num_n_rows)
  
  if (num_plot_rows > 18) {
    SUBSET <- TRUE
    warning_msg <- glue::glue("Warning: This region contains more features than are plottable [{num_plot_rows}]. Selecting first 18.\n")
    cat(file = logfile, warning_msg, append = TRUE)
    
    s_ratio <- dplyr::if_else(
      num_s_rows >= max_plot_rows,
      1,                           # TRUE
      num_s_rows / max_plot_rows   # FALSE
    )
    
    a_ratio <- dplyr::if_else(
      num_a_rows >= max_plot_rows,
      1,                           # TRUE
      num_a_rows / max_plot_rows   # FALSE
    )
    
    n_ratio <- dplyr::if_else(
      num_n_rows >= max_plot_rows,
      1,                           # TRUE
      num_n_rows / max_plot_rows   # FALSE
    )
    
    ratio_sum <- sum(s_ratio, a_ratio, n_ratio)
    
    s_factor <- s_ratio / ratio_sum
    a_factor <- a_ratio / ratio_sum
    n_factor <- n_ratio / ratio_sum
    
    # The number of rows contributed by each strands observations should be 18 +- 1 depending on rounding
    s_rows <- round(max_plot_rows * s_factor)
    a_rows <- round(max_plot_rows * a_factor)
    n_rows <- round(max_plot_rows * n_factor)
    
    sense_df <- sense_df %>%
      dplyr::filter(row_index <= s_rows)
    
    antisense_df <- antisense_df %>%
      dplyr::filter(row_index <= a_rows)
    
    nonsense_df <- nonsense_df %>%
      dplyr::filter(row_index <= n_rows)
  }
  
  #### COMBINE ####
  # Modify the row_index of antisense to start after sense and nonsense to start after antisense
  if (nrow(sense_df) == 0) {
    max_sense_row <- 0
  } else {
    max_sense_row <- max(sense_df$row_index)
  }
  
  if (nrow(antisense_df) == 0) {
    max_antisense_row <- 0
  } else {
    max_antisense_row <- max(antisense_df$row_index)
  }
  
  antisense_df$row_index <- antisense_df$row_index + max_sense_row
  nonsense_df$row_index <- nonsense_df$row_index + max_sense_row + max_antisense_row
  
  plot_df <- dplyr::bind_rows(sense_df, antisense_df, nonsense_df)
  
  # Add in a unique index for each row to allow for polygon grouping
  plot_df <- plot_df %>%
    dplyr::mutate(polygon_idx = dplyr::row_number())
  
  # calculate the x and y coordinates for the genes
  transcripts <- plot_df %>%
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
  
  # Calculating this again in case the data frame was subset
  num_plot_rows <- max(plot_df$row_index)
  
  if (num_plot_rows < 8) {
    TXT_SIZE <- 6
  } else if (num_plot_rows < 12) {
    TXT_SIZE <- 5
  } else {
    TXT_SIZE <- 4
  }
  
  if (SUBSET) {
    plot_title <- "Features filtered to fit on plot"
  } else {
    plot_title <- ""
  }

  # Generate Plot Object
  gtf_plot <- .get_gtf_plot(
    plot_title = plot_title,
    TXT_SIZE = TXT_SIZE,
    reg_start = reg_start,
    reg_stop = reg_stop,
    transcript_poly_data = transcript_poly_data,
    transcript_labels = transcript_labels
  )
  
  return(gtf_plot)
}


# Plots the coverage over an interval
# @param df a dataframe
# @return a histogram plot
.plot_genes_only <- function(gtf, reg_start, reg_stop, logfile) {
  
  roi_length <- reg_stop - reg_start + 1
  
  gtf <- gtf %>%
    dplyr::mutate(
      gene_id = .get_gene_id(attribute)
    )
  
  #### SEPARATE STRANDS ####
  sense_df <- gtf %>%
    dplyr::filter(strand == "+")
  antisense_df <- gtf %>%
    dplyr::filter(strand == "-")
  nonsense_df <- gtf %>%
    dplyr::filter(strand == ".")
  
  #### SORT ####
  # Arrange the data frames for proper plotting
  arrange_rows <- function(plot_df) {
    # Join the transcript ranks back, then build ordering keys and arrange
    arranged_df <- plot_df %>%
      dplyr::arrange(start, end)
    
    return(arranged_df)
  }
  
  sense_df <- arrange_rows(sense_df)
  antisense_df <- arrange_rows(antisense_df)
  nonsense_df <- arrange_rows(nonsense_df)
  
  #### COUNT ####
  # Count the number of plot rows for each strand
  
  num_s_rows <- nrow(sense_df)
  num_a_rows <- nrow(antisense_df)
  num_n_rows <- nrow(nonsense_df)
  
  #### CALCULATE PLOT ROW ORDER ####
  # Since we only have genes in this sub-function,
  #   The plot order is simply the row order
  sense_df <- sense_df %>%
    dplyr::mutate(row_index = dplyr::row_number())
  antisense_df <- antisense_df %>%
    dplyr::mutate(row_index = dplyr::row_number())
  nonsense_df <- nonsense_df %>%
    dplyr::mutate(row_index = dplyr::row_number())
  
  #### SUBSET ####
  SUBSET <- FALSE
  
  # We can't really plot more than 18-20 rows without things starting to look bad, so we're limiting how many features here
  # 18 was chosen due to being divisible by 2 and 3 for cases where there are numerous sense, antisense, and nonstranded features present
  max_plot_rows <- 18
  num_plot_rows <- sum(num_s_rows, num_a_rows, num_n_rows)
  
  if (num_plot_rows > 18) {
    SUBSET <- TRUE
    warning_msg <- glue::glue("Warning: This region contains more features than are plottable [{num_plot_rows}]. Selecting first 18.\n")
    cat(file = logfile, warning_msg, append = TRUE)
    
    s_ratio <- dplyr::if_else(
      num_s_rows >= max_plot_rows,
      1,                           # TRUE
      num_s_rows / max_plot_rows   # FALSE
    )
    
    a_ratio <- dplyr::if_else(
      num_a_rows >= max_plot_rows,
      1,                           # TRUE
      num_a_rows / max_plot_rows   # FALSE
    )
    
    n_ratio <- dplyr::if_else(
      num_n_rows >= max_plot_rows,
      1,                           # TRUE
      num_n_rows / max_plot_rows   # FALSE
    )
    
    ratio_sum <- sum(s_ratio, a_ratio, n_ratio)
    
    s_factor <- s_ratio / ratio_sum
    a_factor <- a_ratio / ratio_sum
    n_factor <- n_ratio / ratio_sum
    
    # The number of rows contributed by each strands observations should be 18 +- 1 depending on rounding
    s_rows <- round(max_plot_rows * s_factor)
    a_rows <- round(max_plot_rows * a_factor)
    n_rows <- round(max_plot_rows * n_factor)
    
    sense_df <- sense_df %>%
      dplyr::filter(row_index <= s_rows)
    
    antisense_df <- antisense_df %>%
      dplyr::filter(row_index <= a_rows)
    
    nonsense_df <- nonsense_df %>%
      dplyr::filter(row_index <= n_rows)
  }
  
  #### COMBINE ####
  # Modify the row_index of antisense to start after sense and nonsense to start after antisense
  if (nrow(sense_df) == 0) {
    max_sense_row <- 0
  } else {
    max_sense_row <- max(sense_df$row_index)
  }
  
  if (nrow(antisense_df) == 0) {
    max_antisense_row <- 0
  } else {
    max_antisense_row <- max(antisense_df$row_index)
  }
  
  antisense_df$row_index <- antisense_df$row_index + max_sense_row
  nonsense_df$row_index <- nonsense_df$row_index + max_sense_row + max_antisense_row
  
  plot_df <- dplyr::bind_rows(sense_df, antisense_df, nonsense_df)
  
  # Add in a unique index for each row to allow for polygon grouping
  plot_df <- plot_df %>%
    dplyr::mutate(polygon_idx = dplyr::row_number())
  
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
  
  # Calculating this again in case the data frame was subset
  num_plot_rows <- max(plot_df$row_index)
  
  if (num_plot_rows < 8) {
    TXT_SIZE <- 6
  } else if (num_plot_rows < 12) {
    TXT_SIZE <- 5
  } else {
    TXT_SIZE <- 4
  }
  
  if (SUBSET) {
    plot_title <- "Features filtered to fit on plot"
  } else {
    plot_title <- ""
  }
  
  # Generate Plot Object
  gtf_plot <- .get_gtf_plot(
    plot_title = plot_title,
    TXT_SIZE = TXT_SIZE,
    reg_start = reg_start,
    reg_stop = reg_stop,
    gene_poly_data = gene_poly_data,
    gene_labels = gene_labels
  )
  
  return(gtf_plot)
}

# Plots the coverage over an interval
# @param df a dataframe
# @return a histogram plot
.plot_transcripts_exons <- function(gtf, reg_start, reg_stop, logfile) {
  
  roi_length <- reg_stop - reg_start + 1
  
  gtf <- gtf %>%
    dplyr::mutate(
      gene_id = .get_gene_id(attribute),
      transcript_id = .get_transcript_id(attribute)
    )
  
  #### SEPARATE STRANDS ####
  sense_df <- gtf %>%
    dplyr::filter(strand == "+")
  antisense_df <- gtf %>%
    dplyr::filter(strand == "-")
  nonsense_df <- gtf %>%
    dplyr::filter(strand == ".")
  
  #### SORT ####
  # Arrange the data frames for proper plotting
  arrange_rows <- function(plot_df) {
    
    # Rank transcripts within each gene by start, end
    tx_rank <- plot_df %>%
      dplyr::filter(feature == "transcript") %>%
      dplyr::arrange(gene_id, start, end) %>%
      dplyr::mutate(transcript_rank = dplyr::row_number()) %>%
      dplyr::select(gene_id, transcript_id, transcript_rank)
    
    # Rank exons within each transcript by start, end
    arranged_df <- plot_df %>%
      dplyr::group_by(transcript_id) %>%
      dplyr::arrange(start, end, .by_group = TRUE) %>%
      dplyr::mutate(exon_rank = dplyr::if_else(feature == "exon", cumsum(feature == "exon"), NA_integer_)) %>%
      dplyr::ungroup()
    
    # Join the transcript ranks back, then build ordering keys and arrange
    arranged_df <- arranged_df %>%
      dplyr::left_join(tx_rank, by = c("gene_id", "transcript_id")) %>%
      dplyr::mutate(
        # within a transcript group: the transcript row, then its exons
        within_tx = dplyr::case_when(
          feature == "transcript" ~ 0L,
          feature == "exon" ~ 1L,
          TRUE ~ 0L
        )
      ) %>%
      dplyr::arrange(
        transcript_rank,   # transcripts in (start, end) order; exons inherit their transcript's rank
        within_tx,         # transcript row before its exons
        exon_rank          # exons in (start, end) order
      ) %>%
      dplyr::select(-within_tx)
    
    return(arranged_df)
  }
  
  sense_df <- arrange_rows(sense_df)
  antisense_df <- arrange_rows(antisense_df)
  nonsense_df <- arrange_rows(nonsense_df)
  
  #### COUNT ####
  # Count the number of plot rows for each strand
  
  sense_tx_plot_rows <- nrow(
    sense_df %>%
      dplyr::filter(feature == "transcript")
  )
  
  # Note this is intended to count groups of exons with the same transcript_id instead of total observations
  sense_exon_rows <- nrow(
    sense_df %>%
      dplyr::filter(feature == "exon") %>%
      dplyr::select(transcript_id) %>%
      dplyr::distinct()
  )
  
  num_s_rows <- sum(sense_tx_plot_rows, sense_exon_rows)
  
  # Antisense rows
  antisense_tx_plot_rows <- nrow(
    antisense_df %>%
      dplyr::filter(feature == "transcript")
  )
  
  # Note this is intended to count groups of exons with the same transcript_id instead of total observations
  antisense_exon_rows <- nrow(
    antisense_df %>%
      dplyr::filter(feature == "exon") %>%
      dplyr::select(transcript_id) %>%
      dplyr::distinct()
  )
  
  num_a_rows <- sum(antisense_tx_plot_rows, antisense_exon_rows)
  
  # Nonsense rows
  nonsense_tx_plot_rows <- nrow(
    nonsense_df %>%
      dplyr::filter(feature == "transcript")
  )
  
  # Note this is intended to count groups of exons with the same transcript_id instead of total observations
  nonsense_exon_rows <- nrow(
    nonsense_df %>%
      dplyr::filter(feature == "exon") %>%
      dplyr::select(transcript_id) %>%
      dplyr::distinct()
  )
  
  num_n_rows <- sum(nonsense_tx_plot_rows, nonsense_exon_rows)
  
  #### CALCULATE PLOT ROW ORDER ####
  
  # Determine row order for plotting
  sense_df <- sense_df %>%
    dplyr::mutate(
      first_exon = feature == "exon" & !duplicated(ifelse(feature == "exon", transcript_id, NA)),
      inc = (feature == "transcript") | first_exon,
      row_index = cumsum(inc)
    ) %>%
    dplyr::select(-first_exon, -inc)
  
  antisense_df <- antisense_df %>%
    dplyr::mutate(
      first_exon = feature == "exon" & !duplicated(ifelse(feature == "exon", transcript_id, NA)),
      inc = (feature == "transcript") | first_exon,
      row_index = cumsum(inc)
    ) %>%
    dplyr::select(-first_exon, -inc)
  
  nonsense_df <- nonsense_df %>%
    dplyr::mutate(
      first_exon = feature == "exon" & !duplicated(ifelse(feature == "exon", transcript_id, NA)),
      inc = (feature == "transcript") | first_exon,
      row_index = cumsum(inc)
    ) %>%
    dplyr::select(-first_exon, -inc)
  
  #### SUBSET ####
  SUBSET <- FALSE
  
  # We can't really plot more than 18-20 rows without things starting to look bad, so we're limiting how many features here
  # 18 was chosen due to being divisible by 2 and 3 for cases where there are numerous sense, antisense, and nonstranded features present
  max_plot_rows <- 18
  num_plot_rows <- sum(num_s_rows, num_a_rows, num_n_rows)
  
  if (num_plot_rows > 18) {
    SUBSET <- TRUE
    warning_msg <- glue::glue("Warning: This region contains more features than are plottable [{num_plot_rows}]. Selecting first 18.\n")
    cat(file = logfile, warning_msg, append = TRUE)
    
    s_ratio <- dplyr::if_else(
      num_s_rows >= max_plot_rows,
      1,                           # TRUE
      num_s_rows / max_plot_rows   # FALSE
    )
    
    a_ratio <- dplyr::if_else(
      num_a_rows >= max_plot_rows,
      1,                           # TRUE
      num_a_rows / max_plot_rows   # FALSE
    )
    
    n_ratio <- dplyr::if_else(
      num_n_rows >= max_plot_rows,
      1,                           # TRUE
      num_n_rows / max_plot_rows   # FALSE
    )
    
    ratio_sum <- sum(s_ratio, a_ratio, n_ratio)
    
    s_factor <- s_ratio / ratio_sum
    a_factor <- a_ratio / ratio_sum
    n_factor <- n_ratio / ratio_sum
    
    # The number of rows contributed by each strands observations should be 18 +- 1 depending on rounding
    s_rows <- round(max_plot_rows * s_factor)
    a_rows <- round(max_plot_rows * a_factor)
    n_rows <- round(max_plot_rows * n_factor)
    
    sense_df <- sense_df %>%
      dplyr::filter(row_index <= s_rows)
    
    antisense_df <- antisense_df %>%
      dplyr::filter(row_index <= a_rows)
    
    nonsense_df <- nonsense_df %>%
      dplyr::filter(row_index <= n_rows)
  }
  
  #### COMBINE ####
  # Modify the row_index of antisense to start after sense and nonsense to start after antisense
  if (nrow(sense_df) == 0) {
    max_sense_row <- 0
  } else {
    max_sense_row <- max(sense_df$row_index)
  }
  
  if (nrow(antisense_df) == 0) {
    max_antisense_row <- 0
  } else {
    max_antisense_row <- max(antisense_df$row_index)
  }
  
  antisense_df$row_index <- antisense_df$row_index + max_sense_row
  nonsense_df$row_index <- nonsense_df$row_index + max_sense_row + max_antisense_row
  
  plot_df <- dplyr::bind_rows(sense_df, antisense_df, nonsense_df)
  
  # Add in a unique index for each row to allow for polygon grouping
  plot_df <- plot_df %>%
    dplyr::mutate(polygon_idx = dplyr::row_number())
  
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
  
  # Calculating this again in case the data frame was subset
  num_plot_rows <- max(plot_df$row_index)
  
  if (num_plot_rows < 8) {
    TXT_SIZE <- 6
  } else if (num_plot_rows < 12) {
    TXT_SIZE <- 5
  } else {
    TXT_SIZE <- 4
  }
  
  if (SUBSET) {
    plot_title <- "Features filtered to fit on plot"
  } else {
    plot_title <- ""
  }
  
  # Generate Plot Object
  gtf_plot <- .get_gtf_plot(
    plot_title = plot_title,
    TXT_SIZE = TXT_SIZE,
    reg_start = reg_start,
    reg_stop = reg_stop,
    exons = exons,
    exon_poly_data = exon_poly_data,
    exon_labels = exon_labels,
    transcript_poly_data = transcript_poly_data,
    transcript_labels = transcript_labels
  )
  
  return(gtf_plot)
}

# Plots the coverage over an interval
# @param df a dataframe
# @return a histogram plot
.plot_genes_exons <- function(gtf, reg_start, reg_stop, logfile) {
  
  roi_length <- reg_stop - reg_start + 1

  gtf <- gtf %>%
    dplyr::mutate(
      gene_id = .get_gene_id(attribute),
      transcript_id = .get_transcript_id(attribute)
    )

  #### SEPARATE STRANDS ####
  sense_df <- gtf %>%
    dplyr::filter(strand == "+")
  antisense_df <- gtf %>%
    dplyr::filter(strand == "-")
  nonsense_df <- gtf %>%
    dplyr::filter(strand == ".")
  
  #### SORT ####
  # Arrange the data frames for proper plotting
  arrange_rows <- function(plot_df) {
    if (nrow(plot_df) == 0) return(plot_df)
    
    # Rank genes by (start, end) using the gene rows
    gene_rank <- plot_df %>%
      dplyr::filter(feature == "gene") %>%
      dplyr::arrange(start, end) %>%
      dplyr::mutate(gene_rank = dplyr::row_number()) %>%
      dplyr::select(gene_id, gene_rank)
    
    # Infer transcript spans from exons, then rank transcripts within each gene
    # Use (tx_start = min exon start, tx_end = max exon end) for ordering
    inferred_tx_rank <- plot_df %>%
      dplyr::filter(feature == "exon") %>%
      dplyr::group_by(gene_id, transcript_id) %>%
      dplyr::summarise(
        tx_start = min(start, na.rm = TRUE),
        tx_end = max(end, na.rm = TRUE),
        .groups = "drop"
      ) %>%
      dplyr::arrange(gene_id, tx_start, tx_end) %>%
      dplyr::group_by(gene_id) %>%
      dplyr::mutate(transcript_rank = dplyr::row_number()) %>%
      dplyr::ungroup() %>%
      dplyr::select(gene_id, transcript_id, transcript_rank)
    
    # Rank exons within each transcript by (start, end)
    arranged_df <- plot_df %>%
      dplyr::group_by(transcript_id) %>%
      dplyr::arrange(start, end, .by_group = TRUE) %>%
      dplyr::mutate(exon_rank = dplyr::if_else(feature == "exon", dplyr::row_number(), NA_integer_)) %>%
      dplyr::ungroup()
    
    # 3) Join ranks and arrange: gene first, then transcripts (via inferred ranks), then exons
    arranged_df <- arranged_df %>%
      dplyr::left_join(gene_rank, by = "gene_id") %>%
      dplyr::left_join(inferred_tx_rank, by = c("gene_id", "transcript_id")) %>%
      dplyr::mutate(
        block = dplyr::if_else(feature == "gene", 0L, 1L) # genes before exons
      ) %>%
      dplyr::arrange(
        gene_rank,          # genes ordered by (start, end)
        block,              # gene row first
        transcript_rank,    # transcript order inferred from exon spans
        exon_rank           # exons ordered by (start, end)
      ) %>%
      dplyr::select(-block)
    
    return(arranged_df)
  }
  
  sense_df <- arrange_rows(sense_df)
  antisense_df <- arrange_rows(antisense_df)
  nonsense_df <- arrange_rows(nonsense_df)
  
  #### COUNT ####
  # Count the number of plot rows for each strand
  
  sense_gene_plot_rows <- nrow(
    sense_df %>%
      dplyr::filter(feature == "gene")
  )
  
  # Note this is intended to count groups of exons with the same transcript_id instead of total observations
  sense_exon_rows <- nrow(
    sense_df %>%
      dplyr::filter(feature == "exon") %>%
      dplyr::select(transcript_id) %>%
      dplyr::distinct()
  )
  
  num_s_rows <- sum(sense_gene_plot_rows, sense_exon_rows)
  
  antisense_gene_plot_rows <- nrow(
    antisense_df %>%
      dplyr::filter(feature == "gene")
  )
  
  antisense_exon_rows <- nrow(
    antisense_df %>%
      dplyr::filter(feature == "exon") %>%
      dplyr::select(transcript_id) %>%
      dplyr::distinct()
  )
  
  num_a_rows <- sum(antisense_gene_plot_rows, antisense_exon_rows)
  
  nonsense_gene_plot_rows <- nrow(
    nonsense_df %>%
      dplyr::filter(feature == "gene")
  )
  
  nonsense_exon_rows <- nrow(
    nonsense_df %>%
      dplyr::filter(feature == "exon") %>%
      dplyr::select(transcript_id) %>%
      dplyr::distinct()
  )
  
  num_n_rows <- sum(nonsense_gene_plot_rows, nonsense_exon_rows)
  
  #### CALCULATE PLOT ROW ORDER ####
  
  # Determine row order for plotting
  sense_df <- sense_df %>%
    dplyr::mutate(
      first_exon = feature == "exon" & !duplicated(ifelse(feature == "exon", transcript_id, NA)),
      inc = (feature == "gene") | first_exon,
      row_index = cumsum(inc)
    ) %>%
    dplyr::select(-first_exon, -inc)
  
  antisense_df <- antisense_df %>%
    dplyr::mutate(
      first_exon = feature == "exon" & !duplicated(ifelse(feature == "exon", transcript_id, NA)),
      inc = (feature == "gene") | first_exon,
      row_index = cumsum(inc)
    ) %>%
    dplyr::select(-first_exon, -inc)
  
  nonsense_df <- nonsense_df %>%
    dplyr::mutate(
      first_exon = feature == "exon" & !duplicated(ifelse(feature == "exon", transcript_id, NA)),
      inc = (feature == "gene") | first_exon,
      row_index = cumsum(inc)
    ) %>%
    dplyr::select(-first_exon, -inc)
  
  #### SUBSET ####
  SUBSET <- FALSE
  
  # We can't really plot more than 18-20 rows without things starting to look bad, so we're limiting how many features here
  # 18 was chosen due to being divisible by 2 and 3 for cases where there are numerous sense, antisense, and nonstranded features present
  max_plot_rows <- 18
  num_plot_rows <- sum(num_s_rows, num_a_rows, num_n_rows)
  
  if (num_plot_rows > 18) {
    SUBSET <- TRUE
    warning_msg <- glue::glue("Warning: This region contains more features than are plottable [{num_plot_rows}]. Selecting first 18.\n")
    cat(file = logfile, warning_msg, append = TRUE)
    
    s_ratio <- dplyr::if_else(
      num_s_rows >= max_plot_rows,
      1,                           # TRUE
      num_s_rows / max_plot_rows   # FALSE
    )
    
    a_ratio <- dplyr::if_else(
      num_a_rows >= max_plot_rows,
      1,                           # TRUE
      num_a_rows / max_plot_rows   # FALSE
    )
    
    n_ratio <- dplyr::if_else(
      num_n_rows >= max_plot_rows,
      1,                           # TRUE
      num_n_rows / max_plot_rows   # FALSE
    )
    
    ratio_sum <- sum(s_ratio, a_ratio, n_ratio)
    
    s_factor <- s_ratio / ratio_sum
    a_factor <- a_ratio / ratio_sum
    n_factor <- n_ratio / ratio_sum
    
    # The number of rows contributed by each strands observations should be 18 +- 1 depending on rounding
    s_rows <- round(max_plot_rows * s_factor)
    a_rows <- round(max_plot_rows * a_factor)
    n_rows <- round(max_plot_rows * n_factor)
    
    sense_df <- sense_df %>%
      dplyr::filter(row_index <= s_rows)
    
    antisense_df <- antisense_df %>%
      dplyr::filter(row_index <= a_rows)
    
    nonsense_df <- nonsense_df %>%
      dplyr::filter(row_index <= n_rows)
  }
  
  #### COMBINE ####
  # Modify the row_index of antisense to start after sense and nonsense to start after antisense
  if (nrow(sense_df) == 0) {
    max_sense_row <- 0
  } else {
    max_sense_row <- max(sense_df$row_index)
  }
  
  if (nrow(antisense_df) == 0) {
    max_antisense_row <- 0
  } else {
    max_antisense_row <- max(antisense_df$row_index)
  }
  
  antisense_df$row_index <- antisense_df$row_index + max_sense_row
  nonsense_df$row_index <- nonsense_df$row_index + max_sense_row + max_antisense_row
  
  plot_df <- dplyr::bind_rows(sense_df, antisense_df, nonsense_df)
  
  # Add in a unique index for each row to allow for polygon grouping
  plot_df <- plot_df %>%
    dplyr::mutate(polygon_idx = dplyr::row_number())
  
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
  
  # Calculating this again in case the data frame was subset
  num_plot_rows <- max(plot_df$row_index)
  
  if (num_plot_rows < 8) {
    TXT_SIZE <- 6
  } else if (num_plot_rows < 12) {
    TXT_SIZE <- 5
  } else {
    TXT_SIZE <- 4
  }
  
  if (SUBSET) {
    plot_title <- "Features filtered to fit on plot"
  } else {
    plot_title <- ""
  }

  # Generate Plot Object
  gtf_plot <- .get_gtf_plot(
    plot_title = plot_title,
    TXT_SIZE = TXT_SIZE,
    reg_start = reg_start,
    reg_stop = reg_stop,
    exons = exons,
    exon_poly_data = exon_poly_data,
    exon_labels = exon_labels,
    gene_poly_data = gene_poly_data,
    gene_labels = gene_labels
  )
  
  return(gtf_plot)
}


.plot_genes_transcripts <- function(gtf, reg_start, reg_stop, logfile) {
  
  roi_length <- reg_stop - reg_start + 1
  
  gtf <- gtf %>%
    dplyr::mutate(
      gene_id = .get_gene_id(attribute),
      transcript_id = .get_transcript_id(attribute)
    )
  
  # SEPARATE STRANDS ####
  sense_df <- gtf %>%
    dplyr::filter(strand == "+")
  antisense_df <- gtf %>%
    dplyr::filter(strand == "-")
  nonsense_df <- gtf %>%
    dplyr::filter(strand == ".")
  
  #### SORT ####
  # Arrange the data frames for proper plotting
  arrange_rows <- function(plot_df) {
    
    # Rank genes by their start and end positions
    gene_rank <- plot_df %>%
      dplyr::filter(feature == "gene") %>%
      dplyr::arrange(start, end) %>%
      dplyr::mutate(gene_rank = dplyr::row_number()) %>%
      dplyr::select(gene_id, gene_rank)
    
    # Rank transcripts within each gene by start, end
    tx_rank <- plot_df %>%
      dplyr::filter(feature == "transcript") %>%
      dplyr::arrange(gene_id, start, end) %>%
      dplyr::group_by(gene_id) %>%
      dplyr::mutate(transcript_rank = dplyr::row_number()) %>%
      dplyr::ungroup() %>%
      dplyr::select(gene_id, transcript_id, transcript_rank)
    
    # Join the transcript ranks back, then build ordering keys and arrange
    arranged_df <- plot_df %>%
      dplyr::left_join(gene_rank, by = c("gene_id")) %>%
      dplyr::left_join(tx_rank, by = c("gene_id", "transcript_id")) %>%
      dplyr::mutate(
        # gene first, then everything else
        block = dplyr::if_else(feature == "gene", 0L, 1L)
      ) %>%
      dplyr::arrange(
        gene_rank,         
        block,            # genes ordered by (start, end)
        transcript_rank   # transcripts in (start, end) order
      ) %>%
      dplyr::select(-block)
    
    return(arranged_df)
  }
  
  sense_df <- arrange_rows(sense_df)
  antisense_df <- arrange_rows(antisense_df)
  nonsense_df <- arrange_rows(nonsense_df)
  
  #### COUNT ####
  # Count the number of plot rows for each strand
  
  num_s_rows <- nrow(sense_df)
  num_a_rows <- nrow(antisense_df)
  num_n_rows <- nrow(nonsense_df)
  
  #### CALCULATE PLOT ROW ORDER ####
  
  # Determine row order for plotting
  # In this particular case with no exons present, the row order is already established
  # So just mutate in a row_index based on row_number
  sense_df <- sense_df %>%
    dplyr::mutate(row_index = dplyr::row_number())
  antisense_df <- antisense_df %>%
    dplyr::mutate(row_index = dplyr::row_number())
  nonsense_df <- nonsense_df %>%
    dplyr::mutate(row_index = dplyr::row_number())
  
  #### SUBSET ####
  SUBSET <- FALSE
  
  # We can't really plot more than 18-20 rows without things starting to look bad, so we're limiting how many features here
  # 18 was chosen due to being divisible by 2 and 3 for cases where there are numerous sense, antisense, and nonstranded features present
  max_plot_rows <- 18
  num_plot_rows <- sum(num_s_rows, num_a_rows, num_n_rows)
  
  if (num_plot_rows > 18) {
    SUBSET <- TRUE
    warning_msg <- glue::glue("Warning: This region contains more features than are plottable [{num_plot_rows}]. Selecting first 18.\n")
    cat(file = logfile, warning_msg, append = TRUE)

    # The purpose of this command is to set a max value of 1 on an data frame that has the max_plot_rows number of rows or more
    # Any data frame that has less will be a ratio of rows:max_plot_rows
    # This will let us determine how many rows each stranded data frame should contribute
    # This is agnostic to how many rows a data frame has if it is equal to or over the max_plot_rows value
    #   So in a case where sense_df has 4000 plottable rows, antisense_df has 18 plottable rows, and nonsense_df has 0 plottable rows,
    #   Both sense and antisense would contribute 9 rows.
    s_ratio <- dplyr::if_else(
      num_s_rows >= max_plot_rows,
      1,                           # TRUE
      num_s_rows / max_plot_rows   # FALSE
    )
    
    a_ratio <- dplyr::if_else(
      num_a_rows >= max_plot_rows,
      1,                           # TRUE
      num_a_rows / max_plot_rows   # FALSE
    )
    
    n_ratio <- dplyr::if_else(
      num_n_rows >= max_plot_rows,
      1,                           # TRUE
      num_n_rows / max_plot_rows   # FALSE
    )
    
    ratio_sum <- sum(s_ratio, a_ratio, n_ratio)
    
    s_factor <- s_ratio / ratio_sum
    a_factor <- a_ratio / ratio_sum
    n_factor <- n_ratio / ratio_sum
    
    # The number of rows contributed by each strands observations should be 18 +- 1 depending on rounding
    s_rows <- round(max_plot_rows * s_factor)
    a_rows <- round(max_plot_rows * a_factor)
    n_rows <- round(max_plot_rows * n_factor)
    
    sense_df <- sense_df %>%
      dplyr::filter(row_index <= s_rows)
    
    antisense_df <- antisense_df %>%
      dplyr::filter(row_index <= a_rows)
    
    nonsense_df <- nonsense_df %>%
      dplyr::filter(row_index <= n_rows)
  }
  
  #### COMBINE ####
  # Modify the row_index of antisense to start after sense and nonsense to start after antisense
  if (nrow(sense_df) == 0) {
    max_sense_row <- 0
  } else {
    max_sense_row <- max(sense_df$row_index)
  }
  
  if (nrow(antisense_df) == 0) {
    max_antisense_row <- 0
  } else {
    max_antisense_row <- max(antisense_df$row_index)
  }
  
  antisense_df$row_index <- antisense_df$row_index + max_sense_row
  nonsense_df$row_index <- nonsense_df$row_index + max_sense_row + max_antisense_row
  
  plot_df <- dplyr::bind_rows(sense_df, antisense_df, nonsense_df)
  
  # Add in a unique index for each row to allow for polygon grouping
  plot_df <- plot_df %>%
    dplyr::mutate(polygon_idx = dplyr::row_number())
    
  # calculate the x and y coordinates for the transcripts/genes
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
  
  # Calculating this again in case the data frame was subset
  num_plot_rows <- max(plot_df$row_index)
  
  if (num_plot_rows < 8) {
    TXT_SIZE <- 6
  } else if (num_plot_rows < 12) {
    TXT_SIZE <- 5
  } else {
    TXT_SIZE <- 4
  }
  
  if (SUBSET) {
    plot_title <- "Features filtered to fit on plot"
  } else {
    plot_title <- ""
  }
  
  # Generate Plot Object
  gtf_plot <- .get_gtf_plot(
    plot_title = plot_title,
    TXT_SIZE = TXT_SIZE,
    reg_start = reg_start,
    reg_stop = reg_stop,
    gene_poly_data = gene_poly_data,
    gene_labels = gene_labels,
    transcript_poly_data = transcript_poly_data,
    transcript_labels = transcript_labels
  )
  
  return(gtf_plot)
}

# Plots the coverage over an interval
# @param df a dataframe
# @return a histogram plot
.plot_genes_exon_transcripts <- function(gtf, reg_start, reg_stop, logfile) {
  
  roi_length <- reg_stop - reg_start + 1
  
  gtf <- gtf %>%
    dplyr::mutate(
      gene_id = .get_gene_id(attribute),
      transcript_id = .get_transcript_id(attribute)
    )

  #### SEPARATE STRANDS ####
  sense_df <- gtf %>%
    dplyr::filter(strand == "+")
  antisense_df <- gtf %>%
    dplyr::filter(strand == "-")
  nonsense_df <- gtf %>%
    dplyr::filter(strand == ".")
  
  #### SORT ####
  # Arrange the data frames for proper plotting
  arrange_rows <- function(plot_df) {
    
    # Rank genes by their start and end positions
    gene_rank <- plot_df %>%
      dplyr::filter(feature == "gene") %>%
      dplyr::arrange(start, end) %>%
      dplyr::mutate(gene_rank = dplyr::row_number()) %>%
      dplyr::select(gene_id, gene_rank)
    
    # Rank transcripts within each gene by start, end
    tx_rank <- plot_df %>%
      dplyr::filter(feature == "transcript") %>%
      dplyr::arrange(gene_id, start, end) %>%
      dplyr::group_by(gene_id) %>%
      dplyr::mutate(transcript_rank = dplyr::row_number()) %>%
      dplyr::ungroup() %>%
      dplyr::select(gene_id, transcript_id, transcript_rank)
    
    # Rank exons within each transcript by start, end
    arranged_df <- plot_df %>%
      dplyr::group_by(transcript_id) %>%
      dplyr::arrange(start, end, .by_group = TRUE) %>%
      dplyr::mutate(exon_rank = dplyr::if_else(feature == "exon", cumsum(feature == "exon"), NA_integer_)) %>%
      dplyr::ungroup()
    
    # Join the transcript ranks back, then build ordering keys and arrange
    arranged_df <- arranged_df %>%
      dplyr::left_join(gene_rank, by = c("gene_id")) %>%
      dplyr::left_join(tx_rank, by = c("gene_id", "transcript_id")) %>%
      dplyr::mutate(
        # gene first, then everything else
        block = dplyr::if_else(feature == "gene", 0L, 1L),
        # within a transcript group: the transcript row, then its exons
        within_tx = dplyr::case_when(
          feature == "transcript" ~ 0L,
          feature == "exon" ~ 1L,
          TRUE ~ 0L # gene rows (ignored by transcript_rank)
        )
      ) %>%
      dplyr::arrange(
        gene_rank,         
        block,             # genes ordered by (start, end)
        transcript_rank,   # transcripts in (start, end) order; exons inherit their transcript's rank
        within_tx,         # transcript row before its exons
        exon_rank,         # exons in (start, end) order
        start,             # final tiebreakers - probably not needed
        end
      ) %>%
      dplyr::select(-block, -within_tx)
    
    return(arranged_df)
  }
  
  sense_df <- arrange_rows(sense_df)
  antisense_df <- arrange_rows(antisense_df)
  nonsense_df <- arrange_rows(nonsense_df)
  
  #### COUNT ####
  # Count the number of plot rows for each strand
  
  sense_gene_tx_plot_rows <- nrow(sense_df %>%
    dplyr::filter(feature == "gene" | feature == "transcript"))
    
  # sense_tx_plot_rows <- nrow(sense_df %>%
  #   dplyr::filter(feature == "transcript"))
  
  # Note this is intended to count groups of exons with the same transcript_id instead of total observations
  sense_exon_rows <- nrow(sense_df %>%
    dplyr::filter(feature == "exon") %>%
    dplyr::select(transcript_id) %>%
    dplyr::distinct()
  )
  
  num_s_rows <- sum(sense_gene_tx_plot_rows, sense_exon_rows)
  
  # Antisense rows
  antisense_gene_tx_plot_rows <- nrow(
    antisense_df %>%
      dplyr::filter(feature == "gene" | feature == "transcript"))
  
  # Note this is intended to count groups of exons with the same transcript_id instead of total observations
  antisense_exon_rows <- nrow(
    antisense_df %>%
    dplyr::filter(feature == "exon") %>%
    dplyr::select(transcript_id) %>%
    dplyr::distinct()
  )
  
  num_a_rows <- sum(antisense_gene_tx_plot_rows, antisense_exon_rows)
  
  # Nonsense rows
  nonsense_gene_tx_plot_rows <- nrow(
    nonsense_df %>%
    dplyr::filter(feature == "gene" | feature == "transcript"))
  
  # Note this is intended to count groups of exons with the same transcript_id instead of total observations
  nonsense_exon_rows <- nrow(
    nonsense_df %>%
    dplyr::filter(feature == "exon") %>%
    dplyr::select(transcript_id) %>%
    dplyr::distinct()
  )
  
  num_n_rows <- sum(nonsense_gene_tx_plot_rows, nonsense_exon_rows)
  
  #### CALCULATE PLOT ROW ORDER ####
  
  # Determine row order for plotting
  sense_df <- sense_df %>%
    dplyr::mutate(
      first_exon = feature == "exon" & !duplicated(ifelse(feature == "exon", transcript_id, NA)),
      inc = (feature == "gene") | (feature == "transcript") | first_exon,
      row_index = cumsum(inc)
    ) %>%
    dplyr::select(-first_exon, -inc)
  
  antisense_df <- antisense_df %>%
    dplyr::mutate(
      first_exon = feature == "exon" & !duplicated(ifelse(feature == "exon", transcript_id, NA)),
      inc = (feature == "gene") | (feature == "transcript") | first_exon,
      row_index = cumsum(inc)
    ) %>%
    dplyr::select(-first_exon, -inc)
  
  nonsense_df <- nonsense_df %>%
    dplyr::mutate(
      first_exon = feature == "exon" & !duplicated(ifelse(feature == "exon", transcript_id, NA)),
      inc = (feature == "gene") | (feature == "transcript") | first_exon,
      row_index = cumsum(inc)
    ) %>%
    dplyr::select(-first_exon, -inc)
  
  #### SUBSET ####
  SUBSET <- FALSE
  
  # We can't really plot more than 18-20 rows without things starting to look bad, so we're limiting how many features here
  # 18 was chosen due to being divisible by 2 and 3 for cases where there are numerous sense, antisense, and nonstranded features present
  max_plot_rows <- 18
  num_plot_rows <- sum(num_s_rows, num_a_rows, num_n_rows)
  
  if (num_plot_rows > 18) {
    SUBSET <- TRUE
    warning_msg <- glue::glue("Warning: This region contains more features than are plottable [{num_plot_rows}]. Selecting first 18.\n")
    cat(file = logfile, warning_msg, append = TRUE)
    # Subset the data frames to only include up to 18 plotable rows combined
    # Need to make sure each strand get some features if present
    # Need to make sure that all 3 features are still present after subsetting
    # Need to indicate on the plot that this is a subsetted plot
    # Figuring how many rows each data frame should get is annoying
    # There has to be a better way to do this
    
    # Case 1 - Only Sense rows
    
    # Case 2 - Only Sense and Antisense rows
    
    # Case 3 - Only Sense and Nonsense rows
    
    # Case 4 - Sense, Antisense, and Nonsense rows
    
    # Case 5 - Only Antisense rows
    
    # Case 6 - Only Antisense and Nonsense rows
    
    # Case 7 - Only Nonsense rows
    
    s_ratio <- dplyr::if_else(
      num_s_rows >= max_plot_rows,
      1,                           # TRUE
      num_s_rows / max_plot_rows   # FALSE
    )
    
    a_ratio <- dplyr::if_else(
      num_a_rows >= max_plot_rows,
      1,                           # TRUE
      num_a_rows / max_plot_rows   # FALSE
    )
    
    n_ratio <- dplyr::if_else(
      num_n_rows >= max_plot_rows,
      1,                           # TRUE
      num_n_rows / max_plot_rows   # FALSE
    )
    
    ratio_sum <- sum(s_ratio, a_ratio, n_ratio)
    
    s_factor <- s_ratio / ratio_sum
    a_factor <- a_ratio / ratio_sum
    n_factor <- n_ratio / ratio_sum
    
    # The number of rows contributed by each strands observations should be 18 +- 1 depending on rounding
    s_rows <- round(max_plot_rows * s_factor)
    a_rows <- round(max_plot_rows * a_factor)
    n_rows <- round(max_plot_rows * n_factor)
    
    # There is a slight possibility that not all the features that are expected in this subfunction will remain after subsetting
    # TODO we need to add a check for that here, and either warn and return empty or...or I don't know. Randomize the subsetting until we get what we want?? That seems terrible.
    # This seems even more likely for strands with a very small ratio/factor and will be getting very few rows represented.
    # We might want to change this to include all rows from very small data frames. Bleh
    
    sense_df <- sense_df %>%
      dplyr::filter(row_index <= s_rows)
    
    antisense_df <- antisense_df %>%
      dplyr::filter(row_index <= a_rows)
    
    nonsense_df <- nonsense_df %>%
      dplyr::filter(row_index <= n_rows)
  }
  
  #### COMBINE ####
  # Modify the row_index of antisense to start after sense and nonsense to start after antisense
  if (nrow(sense_df) == 0) {
    max_sense_row <- 0
  } else {
    max_sense_row <- max(sense_df$row_index)
  }
  
  if (nrow(antisense_df) == 0) {
    max_antisense_row <- 0
  } else {
    max_antisense_row <- max(antisense_df$row_index)
  }
  
  antisense_df$row_index <- antisense_df$row_index + max_sense_row
  nonsense_df$row_index <- nonsense_df$row_index + max_sense_row + max_antisense_row
  
  plot_df <- dplyr::bind_rows(sense_df, antisense_df, nonsense_df)
  
  # Add in a unique index for each row to allow for polygon grouping
  plot_df <- plot_df %>%
    dplyr::mutate(polygon_idx = dplyr::row_number())
  
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
      exon_idx = dplyr::row_number(),
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
  
  # Calculating this again in case the data frame was subset
  num_plot_rows <- max(plot_df$row_index)
  
  if (num_plot_rows < 8) {
    TXT_SIZE <- 6
  } else if (num_plot_rows < 12) {
    TXT_SIZE <- 5
  } else {
    TXT_SIZE <- 4
  }
  
  if (SUBSET) {
    plot_title <- "Features filtered to fit on plot"
  } else {
    plot_title <- ""
  }
  
  # Generate Plot Object
  gtf_plot <- .get_gtf_plot(
    plot_title = plot_title,
    TXT_SIZE = TXT_SIZE,
    reg_start = reg_start,
    reg_stop = reg_stop,
    exons = exons,
    exon_poly_data = exon_poly_data,
    exon_labels = exon_labels,
    gene_poly_data = gene_poly_data,
    gene_labels = gene_labels,
    transcript_poly_data = transcript_poly_data,
    transcript_labels = transcript_labels
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

.get_gene_id <- Vectorize(function(x) {
  tmp <- unlist(strsplit(x, ";"))
  
  # Locate gene_id in the attributes column
  gene_idx <- grep("gene_id", tmp)
  gene_field <- tmp[gene_idx]
  
  # Strip out the text gene_id
  gene_id <- gsub("gene_id ", "", gene_field)
  
  # In case gene_id wasn't the first attribute, strip another space
  gene_id <- gsub(" ", "", gene_id)
  
  return(gene_id)
})

.get_transcript_id <- Vectorize(function(x) {
  tmp <- unlist(strsplit(x, ";"))
  
  # Locate transcript_id in the attributes column
  transcript_idx <- grep(" transcript_id", tmp)
  transcript_field <- tmp[transcript_idx]
  
  # Strip out the text transcript_id
  transcript_id <- gsub(" transcript_id ", "", transcript_field)
  # In case transcript_id wasn't the first attribute, strip another space
  # transcript_id <- gsub(" ", "", transcript_id)
  
  # If no actual transcript_id exists, set to NA
  if (transcript_id == "") {
    transcript_id <- NA
  }
  
  return(transcript_id)
})

.get_gtf_plot <- function(
    plot_title, TXT_SIZE, reg_start, reg_stop,
    exons = NULL, exon_poly_data = NULL, exon_labels = NULL,
    gene_poly_data = NULL, gene_labels = NULL,
    transcript_poly_data = NULL, transcript_labels = NULL
) {
  # Build plot in layers depending on what is present
  gtf_plot <- ggplot2::ggplot()
  
  # Exon Layers
  if (!is.null(exons)) {
    gtf_plot <- gtf_plot +
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
        )
  }
  
  # Gene Layers
  if (!is.null(gene_poly_data)) {
    gtf_plot <- gtf_plot +
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
      )
  }
  
  # Transcript Layers
  if (!is.null(transcript_poly_data)) {
    gtf_plot <- gtf_plot +
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
      )
  }
  
  # Plot layout
  gtf_plot <- gtf_plot +
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
      ggplot2::labs(x = "Genomic coordinate", title = plot_title) +
      ggplot2::theme_classic() +
      ggplot2::theme(
        axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
        axis.title.y = ggplot2::element_blank(),
        plot.title = ggplot2::element_text(hjust = 0.5)
      )
  
  return(gtf_plot)
}
