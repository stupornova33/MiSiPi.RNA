# run_all function
# run all modules of the package and store metrics for ML
# @param chrom_name a string
# @param reg_stop an integer
# @param reg_start an integer
# @param chromosome an integer representing the chromosome number
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
# @param i
# @param i_total
# @return results

.run_all <- function(chrom_name, reg_start, reg_stop,
                     chromosome, length, bam_file,
                     roi, genome_file,
                     si_pal, pi_pal, plot_output,
                     path_to_RNAfold, path_to_RNAplot,
                     annotate_region, weight_reads,
                     gtf_file, write_fastas, out_type,
                     i, i_total) {
  width <- pos <- start <- end <- NULL

  .inform_iteration(i, i_total, chrom_name)

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

  all_dir <- "run_all"
  logfile <- file.path(all_dir, "run_all_logfile.txt")

  bam_obj <- .open_bam(bam_file)
  bam_header <- Rsamtools::scanBamHeader(bam_obj)
  chr_name <- names(bam_header[["targets"]])
  chr_length <- unname(bam_header[["targets"]])
  bam_header <- NULL
  chromosome <- which(chr_name == chrom_name)

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

  size_dist <- dplyr::bind_rows(forward_dt, reverse_dt) %>%
    dplyr::group_by(width) %>%
    dplyr::summarize(count = sum(count))

  .output_readsize_dist(size_dist, prefix, all_dir, strand = NULL, type = "all")

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

  # TO DO: make a null_res to return
  if (nrow(forward_dt) == 0 && nrow(reverse_dt) == 0) {
    return()
  }

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

  ############################################################################ run siRNA function #######################################################################
  # calculate
  # highest_si_col
  # si_dicerz
  # num_si_dicer_reads
  # hp_perc_paired

  cat(file = logfile, "Beginning siRNA function\n", append = TRUE)

  si_dir <- file.path(all_dir, "siRNA_outputs")
  si_log <- file.path(si_dir, "siRNA_logfile.txt")

  si_res <- .siRNA(
    chrom_name, reg_start, reg_stop,
    length, genome_file,
    bam_file, si_log, si_dir,
    si_pal, plot_output, path_to_RNAfold,
    annotate_region, weight_reads, gtf_file,
    write_fastas, out_type
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

  si_res <- NULL
  max_si_heat <- NULL

  ############################################################################# run miRNA function ####################################################################
  # mi_perc_paired
  # mirna_dicerz
  # mirna_mfe
  # mirna_overlapz

  cat(file = logfile, "Beginning miRNA function\n", append = TRUE)

  mi_dir <- file.path(all_dir, "miRNA_outputs")
  mi_log <- file.path(mi_dir, "miRNA_logfile.txt")

  mi_res <- .miRNA(
    chrom_name, reg_start, reg_stop,
    chromosome, length, "+",
    genome_file, bam_file,
    mi_log, mi_dir,
    plot_output,
    path_to_RNAfold,
    path_to_RNAplot,
    weight_reads,
    write_fastas,
    out_type
  )

  # Look at first result
  # mi_res <- mi_res[[1]]
  mirnaMFE_plus <- mi_res$mfe

  pp_plus <- mi_res$perc_paired
  
  mirna_dicerz_plus <- mi_res$overhangs$ml_zscore[5]

  plus_overlapz <- mean(mi_res$overlaps$ml_zscore[17:19])

  mi_res <- .miRNA(
    chrom_name, reg_start, reg_stop,
    chromosome, length, "-",
    genome_file, bam_file,
    mi_log, mi_dir,
    plot_output,
    path_to_RNAfold,
    path_to_RNAplot,
    weight_reads,
    write_fastas,
    out_type
  )

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

  pi_dir <- file.path(all_dir, "piRNA_outputs")
  pi_log <- file.path(pi_dir, "piRNA_logfile.txt")

  pi_res <- .piRNA(chrom_name, reg_start, reg_stop,
    length, bam_file, genome_file,
    pi_log, pi_dir, pi_pal,
    plot_output = plot_output,
    weight_reads,
    write_fastas,
    out_type
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

  pi_res <- NULL

  ####################################################################### add results to table ########################################################################

  tbl_pref <- strsplit(roi, "[.]")[[1]][1]
  tbl_pref <- unlist(strsplit(tbl_pref, "[/]"))
  tbl_pref <- tbl_pref[length(tbl_pref)]

  tmp <- unlist(strsplit(bam_file, "[/]"))
  input_pref <- tmp[length(tmp)]
  input_pref <- strsplit(input_pref, "[.]")[[1]][1]

  ml_file <- file.path(all_dir, paste0(tbl_pref, "_", input_pref, "_ml.txt"))

  local_ml <- as.matrix(local_ml)

  cat(file = logfile, "Writing machine learning results to table\n", append = TRUE)
  .write.quiet(local_ml, ml_file)
}
