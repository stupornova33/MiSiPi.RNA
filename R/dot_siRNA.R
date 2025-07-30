# siRNA function
# process forward and reverse dt and plot si heatmap
# @param chrom_name a string
# @param reg_start a whole number
# @param reg_stop a whole number
# @param length an integer
# @param genome_file a string
# @param bam_file a string
# @param logfile a string
# @param wkdir a string
# @param pal a string
# @param plot_output a bool, TRUE or FALSE. Default is TRUE
# @param path_to_RNAfold a string
# @param annotate_region a bool, TRUE or FALSE
# @param weight_reads Determines whether read counts will be weighted and with which method. Valid options are "weight_by_prop", "locus_norm", or a user-defined value. Default is none. See MiSiPi documentation for descriptions of the weighting methods.
# @param gtf_file a string
# @param write_fastas a bool, Determines whether siRNA pairs will be written to a fasta file. TRUE or FALSE expected. Default: FALSE
# @param out_type The type of file to write the plots to. Options are "png" or "pdf". Default is PDF.
# @return results

.siRNA <- function(chrom_name, reg_start, reg_stop, length,
                   genome_file, bam_file, logfile, wkdir, pal, plot_output,
                   path_to_RNAfold, annotate_region, weight_reads, gtf_file,
                   write_fastas, out_type, method = c("self", "all"),
                   i = NULL, i_total = NULL) {
  
  # i and i_total will be null if called from run_all
  if (!is.null(i)) {
    .inform_iteration(i, i_total, chrom_name)
  }
  
  prefix <- .get_region_string(chrom_name, reg_start, reg_stop)
  width <- pos <- phased_dist <- phased_num <- phased_z <- phased_dist2 <- plus_num2 <- phased_dist1 <- phased_num1 <- NULL

  # use Rsamtools to process the bam file
  bam_obj <- .open_bam(bam_file, logfile)
  bam_header <- Rsamtools::scanBamHeader(bam_obj)
  chr_name <- names(bam_header[["targets"]])
  chr_length <- unname(bam_header[["targets"]])
  bam_header <- NULL

  cat(file = logfile, paste0("chrom_name: ", chrom_name, " reg_start: ", reg_start - 1, " reg_stop: ", reg_stop - 1, "\n"), append = TRUE)
  cat(file = logfile, "Filtering forward and reverse reads by length\n", append = TRUE)

  # extract reads by strand
  # this creates a list object
  chromP <- .get_chr(bam_obj, chrom_name, reg_start, reg_stop, strand = "plus")
  chromM <- .get_chr(bam_obj, chrom_name, reg_start, reg_stop, strand = "minus")

  # turn the list object into a more useable data frame and filter reads by length,
  # bam only contains pos and width, need to add an end column
  cat(file = logfile, "Making Forward DT\n", append = TRUE)
  forward_dt <- data.table::setDT(.make_si_BamDF(chromP)) %>%
    subset(width <= 32 & width >= 18) %>%
    dplyr::rename(start = pos) %>%
    dplyr::mutate(end = start + width - 1) %>%
    dplyr::group_by_all() %>%
    dplyr::summarize(count = dplyr::n())

  cat(file = logfile, "Making Reverse DT\n", append = TRUE)
  reverse_dt <- data.table::setDT(.make_si_BamDF(chromM)) %>%
    subset(width <= 32 & width >= 18) %>%
    dplyr::rename(start = pos) %>%
    dplyr::mutate(end = start + width - 1) %>%
    dplyr::group_by_all() %>%
    dplyr::summarize(count = dplyr::n())

  if (method == "self") {
    size_dist <- dplyr::bind_rows(forward_dt, reverse_dt) %>%
      dplyr::group_by(width) %>%
      dplyr::summarise(count = sum(count))
    .output_readsize_dist(size_dist, prefix, wkdir, strand = NULL, "siRNA")
    size_dist <- NULL
  }
  
  #stranded_read_dist <- .get_stranded_read_dist(bam_obj, chrom_name, reg_start, reg_stop)
  #.plot_sizes_by_strand(wkdir, stranded_read_dist, chrom_name, reg_start, reg_stop)
  
  chromP <- NULL
  chromM <- NULL
  size_dist <- NULL

  # If the data frames are empty there are no reads, can't do siRNA calculations
  if (nrow(forward_dt) == 0 & nrow(reverse_dt) == 0) {
    cat(file = logfile, "No reads detected on one strand. \n", append = TRUE)
    
    # Setting the ml_zscore to -33 for use in later machine learning calculations
    # Setting zscore to 0 for plotting purposes
    dicer_overhangs <- data.frame(shift = seq(-4, 4),
                                  proper_count = c(rep(0, times = 9)),
                                  zscore = c(rep(0, times = 9)),
                                  ml_zscore = c(rep(-33, times = 9)))
    results <- rep(0, times = 324)
    
  } else {
    cat(file = logfile, "Calculating overhangs\n", append = TRUE)
    
    # Get expanded-weighted reads
    locus_length <- reg_stop - reg_start + 1
    
    forward_dt <- .weight_reads(forward_dt, weight_reads, locus_length, sum(forward_dt$count))
    reverse_dt <- .weight_reads(reverse_dt, weight_reads, locus_length, sum(reverse_dt$count))
    
    # Verify weighted data frames contain reads before proceeding
    if (nrow(forward_dt) == 0 & nrow(reverse_dt) == 0) {
      # Even though empty, results are being kept in case the "all" method is being used
      # In that case, they will be written to a table at the end of the that function
      
      cat(file = logfile, "No reads detected on one strand. \n", append = TRUE)
      
      # the data.frame should be modified if using calc_expand_overhangs
      dicer_overhangs <- data.frame(shift = seq(-4, 4),
                                    proper_count = c(rep(0, times = 9)),
                                    zscore = c(rep(0, times = 9)),
                                    ml_zscore = c(rep(-33, times = 9)))
      results <- rep(0, times = 324)
    } else {
      # Now that the DTs have been weighted and re-expanded,
      # Let's summarize them again and keep track of the duplicates with the column "n"
      # This will be crucial in keeping memory and cpu usage down during .find_overlaps()
      f_summarized <- forward_dt %>%
        dplyr::group_by_all() %>%
        dplyr::count()
      
      r_summarized <- reverse_dt %>%
        dplyr::group_by_all() %>%
        dplyr::count()
      
      # get overlapping reads
      overlaps <- .find_overlaps(f_summarized, r_summarized) %>%
        dplyr::mutate(
          p5_overhang = r1_start - r2_start,
          p3_overhang = r1_end - r2_end
        )
      
      # TODO This function runs very slowly on large loci
      # See if it can be run on the summarized dts
      if (write_fastas == TRUE) .write_proper_overhangs(forward_dt, reverse_dt, wkdir, prefix, overlaps, "")
      
      # calculate the number of dicer pairs for the zscore
      dicer_overhangs <- calc_overhangs(overlaps$r1_start, overlaps$r1_end,
                                        overlaps$r2_start, overlaps$r2_width,
                                        dupes_present = TRUE,
                                        overlaps$r1_dupes, overlaps$r2_dupes
      )
      
      dicer_overhangs$zscore <- .calc_zscore(dicer_overhangs$proper_count)
      dicer_overhangs$ml_zscore <- .calc_ml_zscore(dicer_overhangs$proper_count)
      
      cat(file = logfile, "get_si_overlaps\n", append = TRUE)
      
      # calculate the siRNA pairs for the heatmap
      # TODO See if this can be run on the summarized dts for cpu time improvement
      results <- new_get_si_overlaps(
        reverse_dt$start, reverse_dt$end, reverse_dt$width,
        forward_dt$start, forward_dt$end, forward_dt$width
      )
      
      row.names(results) <- c("18", "", "20", "", "22", "", "24", "", "26", "", "28", "", "30", "", "32")
      colnames(results) <- c("18", "", "20", "", "22", "", "24", "", "26", "", "28", "", "30", "", "32")
    }
  }
  
  # transform the data frame for writing to table by row
  # output is the locus followed by all ml_zscores
  overhang_output <- data.frame(t(dicer_overhangs$ml_zscore))
  colnames(overhang_output) <- dicer_overhangs$shift
  overhang_output <- overhang_output %>% dplyr::mutate(locus = prefix)
  overhang_output <- overhang_output[, c(10, 1:9)]

  # heat output needs to be a matrix, so transform
  heat_output <- t(c(prefix, as.vector(results)))
  
  dicerz_file <- file.path(wkdir, "siRNA_dicerz.txt")
  heatmap_file <- file.path(wkdir, "siRNA_heatmap.txt")
  .write.quiet(overhang_output, dicerz_file)
  .write.quiet(heat_output, heatmap_file)
  
  
  # 7/17/25 - passing dicer_overhangs to dual_strand_hairpin in order to make a combined
  #           plot with the individual strands
  #dicer_plot <- .plot_overhangz(dicer_overhangs, "none")
  
  # run the hairpin function on each strand separately
  dsh <- .dual_strand_hairpin(
    chrom_name, reg_start, reg_stop, length, genome_file, bam_file, logfile,
    wkdir, plot_output, path_to_RNAfold, annotate_region, weight_reads,
    gtf_file, write_fastas, out_type, dicer_overhangs
  )

  
  #### Plot Output ####
  if (plot_output) {
    results_present <- sum(results) != 0
    is_small_locus <- (reg_stop - reg_start + 1) <= 10000
    
    if (results_present) {
      heat_plot <- .plot_heat(results, chrom_name, reg_start, reg_stop, wkdir, "siRNA", pal = pal)
    } else {
      heat_plot <- NULL
    }
    
    #dist <- .get_weighted_read_dist(forward_dt, reverse_dt)
    #size_plot <- .plot_sizes(dist)
    stranded_read_dist <- .get_stranded_read_dist(bam_obj, chrom_name, reg_start, reg_stop)
    size_plot <- .plot_sizes_by_strand(wkdir, stranded_read_dist, chrom_name, reg_start, reg_stop)
    
    if (method == "all") {
      plots <- list()
      plots$prefix <- prefix
      plots$size_plot <- size_plot
      #plots$dicer_plot <- dicer_plot
      plots$overhang_probability_plot <- dsh$overhang_probability_plot
      
      # Wrap heat_plot in ggplotify::as.grob if not null since pheatmaps can't be coerced to grob by default
      if (!is.null(heat_plot)) {
        heat_plot <- ggplotify::as.grob(heat_plot)
      }
      
      plots$heat_plot <- heat_plot
      plots$density_plot <- dsh$density_plot
      #plots$plus_overhang_plot <- dsh$plus_overhang_plot
      #plots$minus_overhang_plot <- dsh$minus_overhang_plot
      plots$arc_plot <- dsh$arc_plot
      plots$phasedz <- dsh$phasedz
      #plots$plus_phasedz <- dsh$plus_phasedz
      #plots$minus_phasedz <- dsh$minus_phasedz
      plots$gtf_plot <- dsh$gtf_plot
    } else {
      plots <- NULL
      .plot_siRNA(dsh, is_small_locus, annotate_region, results_present, size_plot, heat_plot, out_type, prefix, wkdir)
    }
  } else {
    plots <- NULL
  }

  return(list(heat = results, si_dicer = dicer_overhangs, dsh = dsh, plots = plots))
}
