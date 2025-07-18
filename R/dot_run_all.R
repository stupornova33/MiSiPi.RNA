# run_all function
# run all modules of the package and store metrics for ML
# @param chrom_name a string
# @param reg_stop an integer
# @param reg_start an integer
# @param length an integer
# @param bam_file a string
# @param roi a string
# @param genome_file a string
# @param si_pal a string
# @param pi_pal a string
# @param plot_output a bool, TRUE or FALSE
# @param path_to_RNAfold a string
# @param path_to_RNAplot a string
# @param annotate_region a bool, TRUE or FALSE
# @param weight_reads a bool, TRUE or FALSE
# @param gtf_file a string
# @param write_fastas a bool, TRUE or FALSE. Default is FALSE
# @param out_type Specifies whether file types for plots are png or pdf. Default is pdf.
# @param output_dir The current output directory
# @param i
# @param i_total
# @return results

.run_all <- function(chrom_name, reg_start, reg_stop,
                     length, bam_file,
                     roi, genome_file,
                     si_pal, pi_pal, plot_output,
                     path_to_RNAfold, path_to_RNAplot,
                     annotate_region, weight_reads,
                     gtf_file, write_fastas, out_type,
                     output_dir, i, i_total) {
  width <- pos <- start <- end <- NULL

  .inform_iteration(i, i_total, chrom_name)
  
  calling_method <- "all"

  # create empty data table for results
  local_ml <- data.table::data.table(
    locus = numeric(1),
    locus_length = numeric(1),
    log_shap_p = numeric(1),
    auc = numeric(1),
    strand_bias = numeric(1),
    perc_GC = numeric(1),
    ave_size = numeric(1),
    perc_first_nucT = numeric(1),
    perc_A10 = numeric(1),
    highest_si_col = numeric(1),
    si_dicerz = numeric(1),
    num_si_dicer_reads = numeric(1),
    hp_perc_paired = numeric(1),
    hp_phasedz = numeric(1),
    hp_mfe = numeric(1),
    hp_dicerz = numeric(1),
    mi_perc_paired = numeric(1),
    mirna_dicerz = numeric(1),
    mirna_mfe = numeric(1),
    mirna_overlapz = numeric(1),
    pingpong_col = numeric(1),
    max_pi_count = numeric(1),
    max_piz_overlap = numeric(1),
    pi_phasedz = numeric(1),
    pi_phased26z = numeric(1)
  )

  local_ml$locus <- .get_region_string(chrom_name, reg_start, reg_stop)

  local_ml$locus_length <- reg_stop - reg_start + 1

  prefix <- .get_region_string(chrom_name, reg_start, reg_stop)

  ####################################################################### process bam input files #############################################################################

  logfile <- file.path(output_dir, "run_all_log.txt")

  bam_obj <- .open_bam(bam_file)
  bam_header <- Rsamtools::scanBamHeader(bam_obj)
  chr_name <- names(bam_header[["targets"]])
  chr_length <- unname(bam_header[["targets"]])
  bam_header <- NULL

  cat(file = logfile, paste0("chrom_name: ", chrom_name, " reg_start: ", reg_start - 1, " reg_stop: ", reg_stop - 1, "\n"), append = TRUE)
  cat(file = logfile, "Filtering forward and reverse reads by length\n", append = TRUE)

  chromP <- .get_chr(bam_obj, chrom_name, reg_start, reg_stop, strand = "plus")
  chromM <- .get_chr(bam_obj, chrom_name, reg_start, reg_stop, strand = "minus")

  forward_dt <- data.table::setDT(.make_si_BamDF(chromP)) %>%
    subset(width <= 32 & width >= 16) %>%
    dplyr::rename(start = pos) %>%
    dplyr::mutate(end = start + width - 1) %>%
    dplyr::group_by_all() %>%
    # get the number of times a read occurs
    dplyr::summarize(count = dplyr::n()) %>%
    na.omit()

  reverse_dt <- data.table::setDT(.make_si_BamDF(chromM)) %>%
    subset(width <= 32 & width >= 16) %>%
    dplyr::rename(start = pos) %>%
    dplyr::mutate(end = start + width - 1) %>%
    dplyr::group_by_all() %>%
    dplyr::summarize(count = dplyr::n()) %>%
    na.omit()

  stranded_read_dist <- .get_stranded_read_dist(bam_obj, chrom_name, reg_start, reg_stop)
  
  .plot_sizes_by_strand(wkdir, stranded_read_dist, chrom_name, reg_start, reg_stop)
  
  chromP <- NULL
  chromM <- NULL
  size_dist <- NULL

  ############################################################################ get extra metrics for ML ###################################################################

  # calculate
  # shap_p
  # auc
  # strand_bias
  # perc_GC
  # ave_size
  # perc_first_nucT
  # perc_A10

  sizes <- data.frame(width = c(forward_dt$width, reverse_dt$width))

  set.seed(1234)
  sample <- sizes[sample(1:nrow(sizes)), ]
  sample <- head(sample, 5000)
  sizes <- NULL

  if ((length(sample) > 3) && !(length(unique(sample)) == 1)) {
    local_ml$log_shap_p <- log10(as.numeric(unlist(unname(shapiro.test(sample)))[2]))
  } else {
    local_ml$log_shap_p <- 2
  }

  total_read_count <- sum(forward_dt$count) + sum(reverse_dt$count)

  unique_read_count <- nrow(forward_dt) + nrow(reverse_dt)

  local_ml$unique_read_bias <- unique_read_count / total_read_count

  if (nrow(forward_dt) > 0) {
    forward_dt <- .weight_reads(forward_dt, weight_reads = "none", 0L, 0L)
  } else {
    forward_dt <- forward_dt %>%
      dplyr::select(-c(count))
  }

  if (nrow(reverse_dt) > 0) {
    reverse_dt <- .weight_reads(reverse_dt, weight_reads = "none", 0L, 0L)
  } else {
    reverse_dt <- reverse_dt %>% dplyr::select(-c(count))
  }

  if (nrow(forward_dt) + nrow(reverse_dt) > 1) {
    all_widths <- c(forward_dt$width, reverse_dt$width)

    m <- mean(all_widths)
    std <- sd(all_widths)

    if (!(std == 0)) {
      bin_width <- KernSmooth::dpih(all_widths, scalest = "stdev")

      nbins <- seq(min(all_widths) - bin_width,
        max(all_widths) + bin_width,
        by = bin_width
      )

      hist(all_widths,
        density = 20, breaks = 5, prob = TRUE,
        xlab = "Read size", ylim = c(0, 2),
        main = "Normal curve over histogram"
      )

      curve <- curve(dnorm(x, mean = m, sd = std),
        col = "darkblue", lwd = 2, add = TRUE, yaxt = "n"
      )

      local_ml$auc <- sum(diff(curve$x) * (head(curve$y, -1) + tail(curve$y, -1))) / 2
    } else {
      local_ml$auc <- 0
    }
  } else {
    local_ml$auc <- -1
  }

  all_widths <- NULL

  # If no data present, set machine learning defaults
  if (nrow(forward_dt) == 0 && nrow(reverse_dt) == 0) {
    
    local_ml$strand_bias <- -33
    local_ml$perc_GC <- -33
    local_ml$ave_size <- -33
    local_ml$perc_first_nucT <- -33
    local_ml$perc_A10 <- -33
    
  } else {
    perc_plus <- nrow(forward_dt) / (nrow(forward_dt) + nrow(reverse_dt))
    perc_minus <- nrow(reverse_dt) / (nrow(reverse_dt) + nrow(forward_dt))
    
    # combine perc minus and perc plus into "strand bias"
    if (perc_plus > perc_minus) {
      local_ml$strand_bias <- perc_plus
    } else {
      local_ml$strand_bias <- perc_minus
    }
    
    all_seqs <- c(forward_dt$seq, reverse_dt$seq)
    local_ml$perc_GC <- .get_GC_content(all_seqs)
    
    read_dist <- .get_read_size_dist(forward_dt, reverse_dt)
    
    ave_size <- .highest_sizes(read_dist)
    
    local_ml$ave_size <- ave_size
    
    cat(file = logfile, "Creating size plots\n", append = TRUE)
    
    if (plot_output == TRUE) {
      size_dir <- file.path(getwd(), "run_all", "size_plots")
      
      if (!dir.exists(size_dir)) {
        dir.create(size_dir)
      }
      
      plot_context <- paste0(chrom_name, ": ", reg_start, "-", reg_stop)
      
      size_plot <- .plot_sizes(read_dist)
      plot_filename <- paste0(prefix, "_read_size_distribution.", out_type)
      plot_file <- file.path(size_dir, plot_filename)
      ggplot2::ggsave(
        plot = size_plot,
        filename = plot_file,
        device = out_type,
        height = 8,
        width = 11,
        units = "in"
      )
    }
    
    local_ml$perc_first_nucT <- .first_nuc_T(forward_dt, reverse_dt)
    
    all_nuc_10 <- all_seqs %>%
      stringr::str_sub(10, 10)
    all_nuc_10_A <- sum(all_nuc_10 == "A")
    local_ml$perc_A10 <- all_nuc_10_A / length(all_nuc_10)
    
    all_seqs <- NULL
    all_nuc_10 <- NULL
    all_nuc_10_A <- NULL
    max_sizes <- NULL
    read_dist <- NULL
  }

  ############################################################################ run siRNA function #######################################################################
  # calculate
  # highest_si_col
  # si_dicerz
  # num_si_dicer_reads
  # hp_perc_paired

  cat(file = logfile, "Beginning siRNA function\n", append = TRUE)

  si_dir <- file.path(output_dir, "siRNA")
  si_log <- file.path(si_dir, "siRNA_log.txt")

  si_res <- .siRNA(
    chrom_name, reg_start, reg_stop,
    length, genome_file,
    bam_file, si_log, si_dir,
    si_pal, plot_output, path_to_RNAfold,
    annotate_region, weight_reads, gtf_file,
    write_fastas, out_type, calling_method
  )
  
  # If plots object is not in si_res, it will return NULL
  siRNA_plots <- si_res$plots

  max_si_heat <- .get_max_si_heat(si_res)
  local_ml$highest_si_col <- max_si_heat$highest_si_col
  
  si_dicerz <- si_res$si_dicer$ml_zscore[5]
  local_ml$si_dicerz <- si_dicerz

  plus_perc_paired <- si_res$dsh$plus_res$perc_paired
  minus_perc_paired <- si_res$dsh$minus_res$perc_paired

  # changed 3/25 to be RPM
  local_ml$num_si_dicer_reads <- (si_res$si_dicer$proper_count[5] * 1000000) / total_read_count
  local_ml$hp_perc_paired <- max(plus_perc_paired, minus_perc_paired)

  ######################################################################### get hairpin-specific results ###############################################################

  plus_phasedz <- si_res$dsh$plus_res$phased_tbl.phased_mlz
  plus_mean <- mean(plus_phasedz[1:4])
  
  minus_phasedz <- si_res$dsh$minus_res$phased_tbl.phased_mlz
  minus_mean <- mean(minus_phasedz[1:4])

  local_ml$hp_phasedz <- max(plus_mean, minus_mean)
  local_ml$hp_mfe <- min(unlist(unname(si_res$dsh$minus_res$MFE)), unlist(unname(si_res$dsh$plus_res$MFE)))

  plus_dicerz <- si_res$dsh$plus_res$hp_overhang_mlz
  minus_dicerz <- si_res$dsh$minus_res$hp_overhang_mlz

  local_ml$hp_dicerz <- max(plus_dicerz, minus_dicerz)

  si_res <- NULL
  max_si_heat <- NULL

  ############################################################################# run miRNA function ####################################################################
  # mi_perc_paired
  # mirna_dicerz
  # mirna_mfe
  # mirna_overlapz

  cat(file = logfile, "Beginning miRNA function\n", append = TRUE)

  mi_dir <- file.path(output_dir, "miRNA")
  mi_log <- file.path(mi_dir, "miRNA_log.txt")

  mi_res <- .miRNA(
    chrom_name, reg_start, reg_stop,
    length, "+",
    genome_file, bam_file,
    mi_log, mi_dir,
    plot_output,
    path_to_RNAfold,
    path_to_RNAplot,
    weight_reads,
    write_fastas,
    out_type,
  )

  miRNA_plus_plots <- mi_res$plots
  miRNA_plus_overhangs <- mi_res$overhangs
  miRNA_plus_overlaps <- mi_res$overlaps

  # Look at first result
  # mi_res <- mi_res[[1]]
  mirnaMFE_plus <- mi_res$mfe

  pp_plus <- mi_res$perc_paired
  
  mirna_dicerz_plus <- mi_res$overhangs$ml_zscore[5]

  plus_overlapz <- mean(mi_res$overlaps$ml_zscore[17:19])

  mi_res <- .miRNA(
    chrom_name, reg_start, reg_stop,
    length, "-",
    genome_file, bam_file,
    mi_log, mi_dir,
    plot_output,
    path_to_RNAfold,
    path_to_RNAplot,
    weight_reads,
    write_fastas,
    out_type,
  )
  
  miRNA_minus_plots <- mi_res$plots
  miRNA_minus_overhangs <- mi_res$overhangs
  miRNA_minus_overlaps <- mi_res$overlaps
  
  miRNA_dicer_overhang_plot <- .plot_miRNA_dicer_overhang_probability(miRNA_plus_overhangs, miRNA_minus_overhangs)
  miRNA_overlap_probability_plot <- .plot_miRNA_overlap_probability(miRNA_plus_overlaps, miRNA_minus_overlaps)

  # mi_res <- mi_res[[1]]
  mirnaMFE_minus <- mi_res$mfe
  
  pp_minus <- mi_res$perc_paired
  
  mirna_dicerz_minus <- mi_res$overhangs$ml_zscore[5]
  
  minus_overlapz <- mean(mi_res$overlaps$ml_zscore[17:19])


  if (is.na(mirnaMFE_minus) && !is.na(mirnaMFE_plus)) {
    local_ml$mirna_mfe <- mirnaMFE_plus
  } else if (is.na(mirnaMFE_plus) && !is.na(mirnaMFE_minus)) {
    local_ml$mirna_mfe <- mirnaMFE_minus
  } else if (is.na(mirnaMFE_minus) && is.na(mirnaMFE_plus)) {
    local_ml$mirna_mfe <- 0
  } else {
    local_ml$mirna_mfe <- min(mirnaMFE_plus, mirnaMFE_minus)
  }
  
  local_ml$mi_perc_paired <- max(pp_plus, pp_minus)
  local_ml$mirna_dicerz <- max(mirna_dicerz_plus, mirna_dicerz_minus)

  if (minus_overlapz == -33 && plus_overlapz != -33) {
    local_ml$mirna_overlapz <- plus_overlapz
  } else if (minus_overlapz != -33 && plus_overlapz == -33) {
    local_ml$mirna_overlapz <- minus_overlapz
  } else if (minus_overlapz == -33 && plus_overlapz == -33) {
    local_ml$mirna_overlapz <- -33
  } else {
    local_ml$mirna_overlapz <- max(minus_overlapz, plus_overlapz)
  }
  

  ############################################################################# run piRNA function ####################################################################
  # calculates pingpong_col
  # max_pi_count
  # max_piz_overlap

  cat(file = logfile, "Begin piRNA function\n", append = TRUE)

  pi_dir <- file.path(output_dir, "piRNA")
  pi_log <- file.path(pi_dir, "piRNA_log.txt")

  pi_res <- .piRNA(chrom_name, reg_start, reg_stop,
    length, bam_file, genome_file,
    pi_log, pi_dir, pi_pal,
    plot_output = plot_output,
    weight_reads,
    write_fastas,
    out_type,
    calling_method
  )
  
  piRNA_plots <- pi_res$plots

  if (sum(pi_res$heat_results) != 0) {
    max_pi_heat <- .get_max_pi_heat(pi_res)
    local_ml$pingpong_col <- max_pi_heat$highest_pi_col
    # changed pi_count to CPM
    local_ml$max_pi_count <- ((max_pi_heat$highest_pi_count) * 1000000) / total_read_count
    local_ml$max_piz_overlap <- .get_max_zscore(unlist(pi_res$z_df$ml_zscore), unlist(pi_res$z_df$Overlap))[[1]]
  } else {
    local_ml$pingpong_col <- -33
    local_ml$max_pi_count <- -33
    local_ml$max_piz_overlap <- -33
  }


  max_pi_heat <- NULL

  ## extract phasing results
  phasedz_plus <- pi_res$phased_plus_z
  phasedz26_plus <- pi_res$phased_26plus_z
  phasedz_minus <- pi_res$phased_minus_z
  phasedz26_minus <- pi_res$phased_26minus_z


  if (phasedz_minus == -33 & phasedz_plus != -33) {
    local_ml$pi_phasedz <- phasedz_plus
  } else if (phasedz_minus != -33 & phasedz_plus == -33) {
    local_ml$pi_phasedz <- phasedz_minus
  } else if (phasedz_minus == -33 & phasedz_plus == -33) {
    local_ml$pi_phasedz <- -33
  } else {
    local_ml$pi_phasedz <- max(phasedz_plus, phasedz_minus)
  }
  
  if (phasedz26_minus == -33 & phasedz26_plus != -33) {
    local_ml$pi_phased26z <- phasedz26_plus
  } else if (phasedz26_minus != -33 & phasedz26_plus == -33) {
    local_ml$pi_phased26z <- phasedz26_minus
  } else if (phasedz26_minus == -33 & phasedz26_plus == -33) {
    local_ml$pi_phased26z <- -33
  } else {
    local_ml$pi_phased26z <- max(phasedz26_plus, phasedz26_minus)
  }

  pi_res <- NULL

  ####################################################################### add results to table ########################################################################

  tbl_pref <- strsplit(roi, "[.]")[[1]][1]
  tbl_pref <- unlist(strsplit(tbl_pref, "[/]"))
  tbl_pref <- tbl_pref[length(tbl_pref)]

  tmp <- unlist(strsplit(bam_file, "[/]"))
  input_pref <- tmp[length(tmp)]
  input_pref <- strsplit(input_pref, "[.]")[[1]][1]

  ml_file <- file.path(output_dir, paste0(tbl_pref, "_", input_pref, "_ml.txt"))

  local_ml <- as.matrix(local_ml)

  cat(file = logfile, "Writing machine learning results to table\n", append = TRUE)
  .write.quiet(local_ml, ml_file)
 
  
  #### Make combined plot for current locus #### 
  # Set all possible plots to NULL initially
  density_plot <- NULL
  distribution_plot <- NULL
  # miRNA specific plots
  # Generated above
  
  # piRNA specific plots
  piRNA_overlap_probability_plot <- NULL # piRNA_plots$z
  piRNA_proper_overlaps_by_size_plot <- NULL # piRNA_plots$heat_plot
  piRNA_phasing_probability_plot <- NULL # piRNA_plots$phased_plot
  # siRNA specific plots
  siRNA_arc_plot <- NULL # siRNA_plots$arc_plot
  siRNA_dicer_overhang_probability_plot <- NULL # siRNA_plots$overhang_probability_plot
  siRNA_phasing_probability_plot <- NULL # siRNA_plots$phasedz
  siRNA_proper_overhangs_by_size_plot <- NULL # siRNA_plots$heat_plot
  siRNA_gtf_plot <- NULL # siRNA_plots$gtf_plot
  
  # miRNA Plus Strand
  if (!is.null(miRNA_plus_plots)) {
    # Redundant plots
    density_plot <- miRNA_plus_plots$density
    distribution_plot <- miRNA_plus_plots$distribution
  }
  
  # miRNA Minus Strand
  if (!is.null(miRNA_minus_plots)) {
    # Redundant Plots
    if (is.null(density_plot)) {
      density_plot <- miRNA_minus_plots$density
    }
    if (is.null(distribution_plot)) {
      distribution_plot <- miRNA_minus_plots$distribution
    }
  }
  
  # piRNA
  if (!is.null(piRNA_plots)) {
    # Redundant Plots
    if (is.null(density_plot)) {
      density_plot <- piRNA_plots$density_plot
    }
    if (is.null(distribution_plot)) {
      distribution_plot <- piRNA_plots$dist_plot
    }
    # piRNA specific plots
    piRNA_overlap_probability_plot <- piRNA_plots$z
    piRNA_proper_overlaps_by_size_plot <-piRNA_plots$heat_plot
    piRNA_phasing_probability_plot <- piRNA_plots$phased_plot
  }
  
  # siRNA
  if (!is.null(siRNA_plots)) {
    # Redundant Plots
    if (is.null(density_plot)) {
      density_plot <- siRNA_plots$density_plot
    }
    if (is.null(distribution_plot)) {
      distribution_plot <- siRNA_plots$size_plot
    }
    # siRNA specific plots
    siRNA_arc_plot <- siRNA_plots$arc_plot
    siRNA_dicer_overhang_probability_plot <- siRNA_plots$overhang_probability_plot
    siRNA_phasing_probability_plot <- siRNA_plots$phasedz
    siRNA_proper_overhangs_by_size_plot <- siRNA_plots$heat_plot
    siRNA_gtf_plot <- siRNA_plots$gtf_plot
  }
  
  # Generate plot object
  # There are currently a total of 18 possible plots
  # Work will be done to considate these plots soon, but for now,
  # we'll test out making a plot that is 2 columns wide and up to
  # 9 columns long
  
  # miRNA_dicer_overhang_plot
  # miRNA_overlap_probability_plot
  
  # Check if each plot is null, if not place in all_plots list
  all_plots <- list()
  i <- 1
  
  if (!is.null(distribution_plot)) {
    all_plots[[i]] <- distribution_plot
    i <- i + 1
  }
  if (!is.null(density_plot)) {
    all_plots[[i]] <- density_plot
    i <- i + 1
  }
  
  if (!is.null(miRNA_dicer_overhang_plot)) {
    all_plots[[i]] <- miRNA_dicer_overhang_plot
    i <- i + 1
  }
  if (!is.null(miRNA_overlap_probability_plot)) {
    all_plots[[i]] <- miRNA_overlap_probability_plot
    i <- i + 1
  }
  if (!is.null(piRNA_overlap_probability_plot)) {
    all_plots[[i]] <- piRNA_overlap_probability_plot
    i <- i + 1
  }
  if (!is.null(piRNA_proper_overlaps_by_size_plot)) {
    all_plots[[i]] <- piRNA_proper_overlaps_by_size_plot
    i <- i + 1
  }
  if (!is.null(piRNA_phasing_probability_plot)) {
    all_plots[[i]] <- piRNA_phasing_probability_plot
    i <- i + 1
  }
  if (!is.null(siRNA_arc_plot)) {
    all_plots[[i]] <- siRNA_arc_plot
    i <- i + 1
  }
  if (!is.null(siRNA_dicer_overhang_probability_plot)) {
    all_plots[[i]] <- siRNA_dicer_overhang_probability_plot
    i <- i + 1
  }
  if (!is.null(siRNA_phasing_probability_plot)) {
    all_plots[[i]] <- siRNA_phasing_probability_plot
    i <- i + 1
  }
  if (!is.null(siRNA_proper_overhangs_by_size_plot)) {
    all_plots[[i]] <- siRNA_proper_overhangs_by_size_plot
    i <- i + 1
  }
  if (!is.null(siRNA_gtf_plot)) {
    all_plots[[i]] <- siRNA_gtf_plot
    i <- i + 1
  }
  
  # Undo the last increment
  i <- i - 1
  
  # Create plot row
  # Since there will be 2 plots per row, we'll just round up
  NUM_PLOT_ROWS <- ceiling(i / 2)
  
  plot_rows <- list() 
  
  for (i in seq(1, length(all_plots), 2)) {
    left_plot <- all_plots[[i]]
    
    # Last Row
    if ((length(all_plots) - i) < 2) {
      if (length(all_plots) %% 2 == 0) {
        right_plot <- all_plots[[i + 1]]
      } else {
        right_plot <- NULL
      }
    } else {
      right_plot <- all_plots[[i + 1]]
    }
    
    current_row <- cowplot::plot_grid(
      left_plot,
      NULL,
      right_plot,
      ncol = 3,
      rel_widths = c(1, 0.3, 1),
      rel_heights = c(1, 1, 1),
      align = "hv", # Align both vertically and horizontally
      axis = "lrtb"
    )
    
    current_row_number <- length(plot_rows) + 1
    plot_rows[[current_row_number]] <- current_row
  }
  
  print(paste0("NUM_PLOT_ROWS: ", NUM_PLOT_ROWS))
  # I hate that we currently have to do this series of if checks
  # TODO: revisit this later
  if (NUM_PLOT_ROWS == 1) {
    all_plot <- plot_rows[[1]]
  } else if (NUM_PLOT_ROWS == 2) {
    all_plot <- cowplot::plot_grid(
      plot_rows[[1]],
      plot_rows[[2]],
      ncol = 1,
      rel_heights = rep(1, NUM_PLOT_ROWS),
      rel_widths = rep(1, NUM_PLOT_ROWS),
      align = "hv",
      axis = "lrtb"
    )
  } else if (NUM_PLOT_ROWS == 3) {
    all_plot <- cowplot::plot_grid(
      plot_rows[[1]],
      plot_rows[[2]],
      plot_rows[[3]],
      ncol = 1,
      rel_heights = rep(1, NUM_PLOT_ROWS),
      rel_widths = rep(1, NUM_PLOT_ROWS),
      align = "hv",
      axis = "lrtb"
    )
  } else if (NUM_PLOT_ROWS == 4) {
    all_plot <- cowplot::plot_grid(
      plot_rows[[1]],
      plot_rows[[2]],
      plot_rows[[3]],
      plot_rows[[4]],
      ncol = 1,
      rel_heights = rep(1, NUM_PLOT_ROWS),
      rel_widths = rep(1, NUM_PLOT_ROWS),
      align = "hv",
      axis = "lrtb"
    )
  } else if (NUM_PLOT_ROWS == 5) {
    all_plot <- cowplot::plot_grid(
      plot_rows[[1]],
      plot_rows[[2]],
      plot_rows[[3]],
      plot_rows[[4]],
      plot_rows[[5]],
      ncol = 1,
      rel_heights = rep(1, NUM_PLOT_ROWS),
      rel_widths = rep(1, NUM_PLOT_ROWS),
      align = "hv",
      axis = "lrtb"
    )
  } else if (NUM_PLOT_ROWS == 6) {
    all_plot <- cowplot::plot_grid(
      plot_rows[[1]],
      plot_rows[[2]],
      plot_rows[[3]],
      plot_rows[[4]],
      plot_rows[[5]],
      plot_rows[[6]],
      ncol = 1,
      rel_heights = rep(1, NUM_PLOT_ROWS),
      rel_widths = rep(1, NUM_PLOT_ROWS),
      align = "hv",
      axis = "lrtb"
    )
  }
  
  height <- 4 * NUM_PLOT_ROWS
  
  if (out_type == "png") {
    grDevices::png(file = file.path(output_dir, "combined_plots", paste(prefix, "combined.png", sep = "_")), height = height, width = 11, units = "in", res = 300)
  } else {
    grDevices::pdf(file = file.path(output_dir, "combined_plots", paste(prefix, "combined.pdf", sep = "_")), height = height, width = 11)
  }
  print(all_plot)
  grDevices::dev.off()
  
}
