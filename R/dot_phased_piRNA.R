# phased piRNA function
# processes reads according to phased piRNA algorith
# plots output
# takes bam object and strand
# returns max_zscore and plots
#
# param strand a character passed in, either "+" or "-"
# param chrom_name a string
# param reg_start a whole number
# param reg_stop a whole number
# param input_file a string
# param logfile a string
# param wkdir a string
# param pal a string
# param plot_output a bool, TRUE or FALSE, default is TRUE
# return max_zscore, plots

.phased_piRNA <- function(strand, chrom_name, reg_start, reg_stop, input_file, logfile, wkdir, pal, plot_output) {
  # Calculates a probability based on observations, mean, and standard deviation

  # process bam input files
  bam_obj <- .open_bam(input_file, logfile)
  bam_header <- Rsamtools::scanBamHeader(bam_obj)
  chr_name <- names(bam_header[["targets"]])
  chr_length <- unname(bam_header[["targets"]])
  bam_header <- width <- pos <- first <- rname <- start <- phased_num.y <- phased_num.x <- phased_dist <- phased_z <- phased26_z <- phased_num <- NULL

  chromosome <- which(chr_name == chrom_name)

  cat(file = paste0(wkdir, logfile), paste0("chrom_name: ", chrom_name, " reg_start: ", reg_start - 1, " reg_stop: ", reg_stop - 1, "\n"), append = TRUE)
  cat(file = paste0(wkdir, logfile), "Filtering forward and reverse reads by length\n", append = TRUE)

  # for the read size dist plot

  read_dist <- .get_read_dist(bam_obj, chrom_name, reg_start, reg_stop)

  if (strand == "+") {
    chrom <- .get_chr(bam_obj, chrom_name, reg_start, reg_stop, strand = "plus")

    filter_dt <- data.table::setDT(.makeBamDF(chrom)) %>%
      base::subset(width <= 32 & width >= 18) %>%
      dplyr::mutate(start = pos, end = pos + (width - 1)) %>%
      dplyr::select(-c(rname, pos)) %>%
      dplyr::distinct()

    filter_dt <- filter_dt %>% dplyr::mutate(rname = chrom_name)

    filter_r1_dt <- filter_dt %>%
      dplyr::filter(first == "T") %>%
      dplyr::mutate(end = start + (width - 1) + 59) %>%
      dplyr::select(-c(first))

    over_26_dt <- filter_dt %>% subset(width >= 26)
    over_26_dt <- over_26_dt %>%
      dplyr::filter(first == "T") %>%
      dplyr::mutate(end = start + (width - 1) + 59) %>%
      dplyr::select(-c(first))

    filter_r2_dt <- filter_dt %>%
      dplyr::filter(first == "T") %>%
      dplyr::select(-c(first))
  } else {
    chrom <- .get_chr(bam_obj, chrom_name, reg_start, reg_stop, strand = "minus")

    filter_dt <- data.table::setDT(.makeBamDF(chrom)) %>%
      base::subset(width <= 32 & width >= 18) %>%
      dplyr::mutate(start = pos, end = pos + (width - 1)) %>%
      dplyr::select(-c(rname, pos)) %>%
      dplyr::distinct()

    filter_dt <- filter_dt %>% dplyr::mutate(rname = chrom_name)

    over_26_dt <- filter_dt %>% subset(width >= 26)

    over_26_dt <- over_26_dt %>%
      dplyr::filter(first == "T") %>%
      dplyr::mutate(end = start + (width - 1) + 59) %>%
      dplyr::select(-c(first))

    filter_r1_dt <- filter_dt %>%
      dplyr::filter(first == "T") %>%
      dplyr::mutate(end = start + (width - 1) + 59) %>%
      dplyr::select(-c(first))

    filter_r2_dt <- filter_dt %>%
      dplyr::filter(first == "T") %>%
      dplyr::select(-c(first))
  }

  all_table <- data.table::data.table(phased_dist = seq(0, 50), phased_num = rep(0, 51))


  if (!nrow(filter_r1_dt) == 0) {
    phased_counts <- .calc_phasing(filter_r1_dt, filter_r2_dt, 50)
  } else {
    phased_counts <- data.table::data.table(phased_dist = 1, phased_num = 0L, phased_z = 0)
    phased_counts <- data.table::setDT(dplyr::full_join(phased_counts, all_table, by = "phased_dist", "phased_num")) %>%
      dplyr::select(-c(phased_num.y)) %>%
      dplyr::rename("phased_num" = phased_num.x)
    phased_counts[is.na(phased_counts)] <- 0
  }

  if (!nrow(over_26_dt) == 0) {
    phased_26_counts <- .calc_phasing(over_26_dt, over_26_dt, 50)
  } else {
    phased_26_counts <- data.table::data.table(phased_dist = 1, phased_num = 0L, phased_z = 0)
    phased_26_counts <- data.table::setDT(dplyr::full_join(phased_26_counts, all_table, by = "phased_dist", "phased_num")) %>%
      dplyr::select(-c(phased_num.y)) %>%
      dplyr::rename("phased_num" = phased_num.x)
    phased_26_counts[is.na(phased_26_counts)] <- 0
  }

  # make the results data table

  phased_26_counts <- phased_26_counts %>%
    dplyr::rename(phased26_dist = phased_dist, phased26_num = phased_num, phased26_z = phased_z)


  df <- cbind(phased_counts, phased_26_counts)

  prefix <- .get_region_string(chrom_name, reg_start, reg_stop)


  phased_output <- phased_counts %>% dplyr::select(c(phased_z))
  phased_output <- t(c(prefix, t(phased_output)))

  phased26_output <- phased_26_counts %>% dplyr::select(c(phased26_z))
  phased26_output <- t(c(prefix, t(phased26_output)))

  phased_file <- file.path(wkdir, "phased_piRNA_zscores.txt")
  phased26_file <- file.path(wkdir, "phased26_piRNA_zscores.txt")

  .write.quiet(phased_output, phased_file)
  .write.quiet(phased26_output, phased26_file)

  if (plot_output == TRUE) {
    dist_plot <- .plot_sizes(read_dist)
    phased_plot <- .plot_phasedz(df)
    all_plot <- cowplot::plot_grid(phased_plot, dist_plot, rel_widths = c(1, 0.7), rel_heights = c(1, 1), ncol = 2)
    grDevices::pdf(file = paste0(wkdir, prefix, "_", strand, "_phased_zscore.pdf"), height = 6, width = 11)
    print(all_plot)
    grDevices::dev.off()
  }

  if (!is.na(sum(phased_counts$phased_z))) {
    ave_z <- mean(phased_counts$phased_z[1:4])
  } else {
    ave_z <- -30
  }


  if (!is.na(sum(phased_26_counts$phased26_z))) {

  }

  ave_26z <- mean(phased_26_counts$phased26_z[1:4])
  return(c(ave_z, ave_26z))
}
