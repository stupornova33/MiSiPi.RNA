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
    subset(width <= 32 & width >= 18) %>%
    dplyr::rename(start = pos) %>%
    dplyr::mutate(end = start + width - 1) %>%
    dplyr::group_by_all() %>%
    # get the number of times a read occurs
    dplyr::summarize(count = dplyr::n()) %>%
    na.omit()

  reverse_dt <- data.table::setDT(.make_si_BamDF(chromM)) %>%
    subset(width <= 32 & width >= 18) %>%
    dplyr::rename(start = pos) %>%
    dplyr::mutate(end = start + width - 1) %>%
    dplyr::group_by_all() %>%
    dplyr::summarize(count = dplyr::n()) %>%
    na.omit()

  size_dist <- dplyr::bind_rows(forward_dt, reverse_dt) %>%
    dplyr::group_by(width) %>%
    dplyr::summarize(count = sum(count))
  # Append read size distribution table to output file
  .output_readsize_dist(size_dist, prefix, output_dir, strand = NULL, type = "all")
  
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
    
    # TODO
    # Move this down into the plot section
    # We are currently using the read size distribution plots returned from sub modules
    # However, there is a slight chances that all 3 sub modules could return null results
    # So we should perhaps generate the size dist plot here for some output to be displayed
    # In each region's combined plot
    #size_plot <- .plot_sizes(read_dist)
    
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
    bam_file, roi, si_log, si_dir,
    si_pal, plot_output, path_to_RNAfold,
    annotate_region, weight_reads, gtf_file,
    write_fastas, out_type, calling_method
  )
  
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

  max_si_heat <- NULL

  ############################################################################# run miRNA function ####################################################################
  # mi_perc_paired
  # mirna_dicerz
  # mirna_mfe
  # mirna_overlapz

  cat(file = logfile, "Beginning miRNA function\n", append = TRUE)

  mi_dir <- file.path(output_dir, "miRNA")
  mi_log <- file.path(mi_dir, "miRNA_log.txt")

  mi_res_plus <- .miRNA(
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

  

  mirnaMFE_plus <- mi_res_plus$mfe

  pp_plus <- mi_res_plus$perc_paired
  
  mirna_dicerz_plus <- mi_res_plus$overhangs$ml_zscore[5]

  plus_overlapz <- mean(mi_res_plus$overlaps$ml_zscore[17:19])

  mi_res_minus <- .miRNA(
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

  mirnaMFE_minus <- mi_res_minus$mfe
  
  pp_minus <- mi_res_minus$perc_paired
  
  mirna_dicerz_minus <- mi_res_minus$overhangs$ml_zscore[5]
  
  minus_overlapz <- mean(mi_res_minus$overlaps$ml_zscore[17:19])

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
    length, bam_file, genome_file, roi,
    pi_log, pi_dir, pi_pal,
    plot_output = plot_output,
    weight_reads,
    write_fastas,
    out_type,
    calling_method
  )

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
 
  
  if (plot_output == FALSE) {
    si_res <- NULL
    mi_res_plus <- NULL
    mi_res_minus <- NULL
    pi_res <- NULL
    return()
  }
  
  
  null_plot <- function(type, reason) {
    p <- ggplot2::ggplot() +
      ggplot2::annotate("text", x = 0.5, y = 0.5, label = reason, size = 5) +
      ggplot2::ggtitle(type) +
      ggplot2::theme_void() +
      ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))
    
    return(p)
  }
  
  plot_title <- function(bam_file, bed_file, genome_file, chrom_name, reg_start, reg_stop, i) {
    locus <- paste0(chrom_name, ":", reg_start, "-", reg_stop)
    now <- format(lubridate::now(), "%Y-%m-%d %H:%M:%S")
    
    misipi_version <- packageVersion("MiSiPi.RNA")
    
    iteration_str <- paste0(i, ")")
    
    p_title <- paste("MiSiPi Results for locus:", locus, "(Bed file line:", iteration_str)
    p_subtitle <- paste("Bam:", bam_file, "| Bed:", bed_file, "| Genome:", genome_file)
    p_caption <- paste("Run at:", now, "with MiSiPi.RNA Version:", misipi_version)

    plot_details <- list()
    plot_details$title <- paste0(p_title, "\n", p_subtitle)
    #plot_details$subtitle <- p_subtitle
    plot_details$caption <- p_caption
    return(plot_details)
  }

  
  
  
  #### Make combined plot for current locus ####
  # Generate General Read Plots
  
  density_data <- .read_densityBySize(chrom_name, reg_start, reg_stop, bam_file, wkdir)
  read_density_plot <- .plot_density(density_data, reg_start, reg_stop)
  
  stranded_read_dist <- .get_stranded_read_dist(bam_obj, chrom_name, reg_start, reg_stop)
  read_distribution_plot <- .plot_sizes_by_strand(stranded_read_dist, chrom_name, reg_start, reg_stop)
  
  # Generate miRNA Plots
  miRNA_plus_overhangs <- mi_res_plus$overhangs
  miRNA_plus_overlaps <- mi_res_plus$overlaps
  mi_res_plus <- NULL
  
  miRNA_minus_overhangs <- mi_res_minus$overhangs
  miRNA_minus_overlaps <- mi_res_minus$overlaps
  mi_res_minus <- NULL
  
  miRNA_dicer_overhang_plot <- .plot_miRNA_dicer_overhang_probability(miRNA_plus_overhangs, miRNA_minus_overhangs)
  miRNA_overlap_probability_plot <- .plot_miRNA_overlap_probability(miRNA_plus_overlaps, miRNA_minus_overlaps)
  
  # Gather siRNA Plots
  siRNA_plots <- si_res$plots
  
  if (is.null(siRNA_plots)) {
    # Generate plot placeholders for missing plots
    siRNA_arc_plot <- NULL
    siRNA_dicer_overhang_probability_plot <- NULL
    siRNA_phasing_probability_plot <- NULL
    siRNA_proper_overhangs_by_size_plot <- NULL
    siRNA_gtf_plot <- NULL
  } else {
    # Arc Plot
    siRNA_arc_plot <- siRNA_plots$arc_plot
    
    # Dicer Overhang Probability
    siRNA_dicer_overhang_probability_plot <- siRNA_plots$overhang_probability_plot
    
    # Phasing Probability
    if (length(siRNA_plots$phasedz) == 1) {
      siRNA_phasing_probability_plot <- NULL
    } else {
      siRNA_phasing_probability_plot <- siRNA_plots$phasedz
    }
    
    # Proper Overhangs by Size
    # NULL plot is now handled closer to plot generation
    siRNA_proper_overhangs_by_size_plot <- siRNA_plots$heat_plot
    
    # GTF Annotation
    siRNA_gtf_plot <- siRNA_plots$gtf_plot
  }
  
  si_res <- NULL
  
  
  
  # Gather piRNA Plots
  piRNA_plots <- pi_res$plots
  
  if (is.null(piRNA_plots)) {
    piRNA_overlap_probability_plot <- NULL
    piRNA_proper_overlaps_by_size_plot <- NULL
    piRNA_phasing_probability_plot <- NULL
  } else {
    # Overlap Probability
    piRNA_overlap_probability_plot <- piRNA_plots$overlap_probability_plot
    
    # Proper Overlaps by Size
    piRNA_proper_overlaps_by_size_plot <- piRNA_plots$heat_plot
    
    # Phasing Probability
    piRNA_phasing_probability_plot <- piRNA_plots$phased_probability_plot
  }

  pi_res <- NULL
  
  # Combined plot is 3 rows and 4 columns and will display in the order shown below
  #
  # Read Distribution Plot   -  Arc Plot      -  siRNA Overhangs by Size  -  piRNA Overlaps by Size
  # miRNA Dcr Overhang Prob  -  Read Density  -  siRNA Dcr Overhang Prob  -  piRNA Overlap Prob
  # miRNA Overlap Prob       -  GTF Plot      -  siRNA Phasing Prob       -  piRNA Phasing Prob
  
  plot_body <- cowplot::plot_grid(
    read_distribution_plot, siRNA_arc_plot, siRNA_dicer_overhang_probability_plot, piRNA_overlap_probability_plot,
    miRNA_dicer_overhang_plot, read_density_plot, siRNA_phasing_probability_plot, piRNA_phasing_probability_plot,
    miRNA_overlap_probability_plot, siRNA_gtf_plot, siRNA_proper_overhangs_by_size_plot, piRNA_proper_overlaps_by_size_plot,
    
    # read_distribution_plot, siRNA_arc_plot, siRNA_proper_overhangs_by_size_plot, piRNA_proper_overlaps_by_size_plot,
    # miRNA_dicer_overhang_plot, read_density_plot, siRNA_dicer_overhang_probability_plot, piRNA_overlap_probability_plot,
    # miRNA_overlap_probability_plot, siRNA_gtf_plot, siRNA_phasing_probability_plot, piRNA_phasing_probability_plot,
    ncol = 4,
    align = "hv",
    axis = "lrtb"
  )
  
  plot_details <- plot_title(bam_file, roi, genome_file, chrom_name, reg_start, reg_stop, i)
  
  plot_title <- cowplot::ggdraw() +
    cowplot::draw_label(
      plot_details$title,
      x = 0.5,
      hjust = 0.5,
    ) +
    ggplot2::theme(plot.margin = ggplot2::margin(0, 0, 0, 7))
  
  plot_caption <- cowplot::ggdraw() +
    cowplot::draw_label(
      plot_details$caption,
      x = 0.5,
      hjust = 0.5
    ) +
    ggplot2::theme(plot.margin = ggplot2::margin(0, 0, 0, 7))
  
  
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
