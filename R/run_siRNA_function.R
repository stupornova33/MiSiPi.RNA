#' run siRNA function
#' process forward and reverse dt and plot si heatmap
#' @param chrom_name a string
#' @param reg_start a whole number
#' @param reg_stop a whole number
#' @param length an integer
#' @param min_read_count an integer
#' @param genome_file a string
#' @param bam_file a string
#' @param logfile a string
#' @param wkdir a string
#' @param pal a string
#' @param plot_output a bool, TRUE or FALSE. Default is TRUE
#' @param path_to_RNAfold a string
#' @param annotate_region a bool, TRUE or FALSE
#' @param weight_reads Determines whether read counts will be weighted and with which method. Valid options are "weight_by_prop", "locus_norm", or a user-defined value. Default is none. See MiSiPi documentation for descriptions of the weighting methods.
#' @param gtf_file a string
#' @param write_fastas a bool, Determines whether siRNA pairs will be written to a fasta file. TRUE or FALSE expected. Default: FALSE
#' @param out_type The type of file to write the plots to. Options are "png" or "pdf". Default is PDF.
#' @return results
#' @export

run_siRNA_function <- function(chrom_name, reg_start, reg_stop, length, min_read_count, genome_file, bam_file, logfile, wkdir,
                           pal, plot_output, path_to_RNAfold, annotate_region, weight_reads, gtf_file, write_fastas, out_type) {
  prefix <- .get_region_string(chrom_name, reg_start, reg_stop)
  print(prefix)
  width <- pos <- phased_dist <- phased_num <- phased_z <- phased_dist2 <- plus_num2 <- phased_dist1 <- phased_num1 <- NULL

  # use Rsamtools to process the bam file
  bam_obj <- OpenBamFile(bam_file, logfile)
  bam_header <- Rsamtools::scanBamHeader(bam_obj)
  chr_name <- names(bam_header[['targets']])
  chr_length <- unname(bam_header[['targets']])
  bam_header <- NULL

  cat(file = paste0(wkdir, logfile), paste0("chrom_name: ", chrom_name, " reg_start: ", reg_start - 1, " reg_stop: ", reg_stop - 1, "\n"), append = TRUE)
  cat(file = paste0(wkdir, logfile), "Filtering forward and reverse reads by length\n", append = TRUE)

  # extract reads by strand
  # this creates a list object
  chromP <- .get_chr(bam_obj, chrom_name, reg_start, reg_stop, strand = "plus")
  chromM <- .get_chr(bam_obj, chrom_name, reg_start, reg_stop, strand = "minus")

  # turn the list object into a more useable data frame and filter reads by length,
  # bam only contains pos and width, need to add an end column
  cat(file = paste0(wkdir, logfile), "Making Forward DT\n", append = TRUE)
  forward_dt <- data.table::setDT(.make_si_BamDF(chromP)) %>%
    subset(width <= 32 & width >= 18) %>%
    dplyr::rename(start = pos) %>%
    dplyr::mutate(end = start + width - 1) %>%
    dplyr::group_by_all() %>%
    # get the number of times a read occurs
    dplyr::summarize(count = dplyr::n())

  cat(file = paste0(wkdir, logfile), "Making Reverse DT\n", append = TRUE)
  reverse_dt <- data.table::setDT(.make_si_BamDF(chromM)) %>%
    subset(width <= 32 & width >= 18) %>%
    dplyr::rename(start = pos) %>%
    dplyr::mutate(end = start + width - 1) %>%
    dplyr::group_by_all() %>%
    dplyr::summarize(count = dplyr::n())

  size_dist <- dplyr::bind_rows(forward_dt, reverse_dt) %>%
    dplyr::group_by(width) %>%
    dplyr::summarise(count = sum(count))

  output_readsize_dist(size_dist, prefix, wkdir, strand = NULL, "siRNA")

  chromP <- NULL
  chromM <- NULL
  size_dist <- NULL

  # If the data frames are empty there are no reads, can't do siRNA calculations
  if (nrow(forward_dt) > 0 & nrow(reverse_dt) > 0) {
    print("f_dt and r_dt are not empty")
    cat(file = paste0(wkdir, logfile), "Calc overhangs\n", append = TRUE)

    # Get expanded-weighted reads
    locus_length <- reg_stop - reg_start + 1
    
    forward_dt <- .weight_reads(forward_dt, weight_reads, locus_length, sum(forward_dt$count))
    reverse_dt <- .weight_reads(reverse_dt, weight_reads, locus_length, sum(reverse_dt$count))

    print("Completed getting weighted dataframes.")

    # check to see if subsetted dfs are empty
    # have to keep doing this at each step otherwise errors will happen
    if (nrow(forward_dt) > 0 & nrow(reverse_dt) > 0) {

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
      #overlaps <- .find_overlaps(f_summarized, r_summarized) %>%
      #  dplyr::mutate(p5_overhang = r2_end - r1_end,
      #                p3_overhang = r2_start - r1_start)# %>%
        #dplyr::filter(p5_overhang >= 0 & p3_overhang >= 0)
      overlaps <- .find_overlaps(f_summarized, r_summarized) %>%
       dplyr::mutate(p5_overhang = r1_start - r2_start,
                     p3_overhang = r1_end - r2_end) #%>%
       #dplyr::filter(p3_overhang >= 0 & p5_overhang >= 0)


      # TODO This function runs very slowly on large loci
      # See if it can be run on the summarized dts
      if (write_fastas == TRUE) write_proper_overhangs(forward_dt, reverse_dt, wkdir, prefix, overlaps, "")

      #calculate the number of dicer pairs for the zscore
      dicer_overhangs <- calc_overhangs(overlaps$r1_start, overlaps$r1_end,
                                        overlaps$r2_start, overlaps$r2_width,
                                        dupes_present = TRUE,
                                        overlaps$r1_dupes, overlaps$r2_dupes)

      dicer_overhangs$Z_score <- .calc_zscore(dicer_overhangs$proper_count)

      cat(file = paste0(wkdir, logfile), "get_si_overlaps\n", append = TRUE)
      # calculate the siRNA pairs for the heatmap

      #system.time(results <- get_si_overlaps(reverse_dt$start, reverse_dt$end, reverse_dt$width,
      #                          forward_dt$start, forward_dt$end, forward_dt$width))

      # TODO See if this can be run on the summarized dts for cpu time improvement
      results <- new_get_si_overlaps(reverse_dt$start, reverse_dt$end, reverse_dt$width,
                                     forward_dt$start, forward_dt$end, forward_dt$width)

      row.names(results) <- c('15','','17','','19','','21','','23','','25','','27','','29','','31','')
      colnames(results) <- c('15','','17','','19','','21','','23','','25','','27','','29','','31','')
    } else {
      # results are being stored also in case the run_all function is being used, at the end they will be written to a table
      cat(file = paste0(wkdir, logfile), "No reads detected on one strand. \n", append = TRUE)
      #the data.frame should be modified if using calc_expand_overhangs
      dicer_overhangs <- data.frame(shift = seq(-4,4), proper_count = c(rep(0, times = 9)), Z_score = c(rep(-33, times = 9)))
      #dicer_overhangs <- data.frame(shift = seq(-8,8), proper_count = c(rep(0, times = 17)), Z_score = c(rep(-33, times = 17)))
      results <- rep(0, times = 324)
    }

  } else {
    cat(file = paste0(wkdir, logfile), "No reads detected on one strand. \n", append = TRUE)
    #dicer_overhangs <- data.frame(shift = seq(-8,8), proper_count = c(rep(0, times = 17)), Z_score = c(rep(-33, times = 17)))
    dicer_overhangs <- data.frame(shift = seq(-4,4), proper_count = c(rep(0, times = 9)), Z_score = c(rep(-33, times = 9)))
    results <- rep(0, times = 324)
  }

  # transform the data frame for writing to table by row
  # output is the locus followed by all zscores
  overhang_output <- data.frame(t(dicer_overhangs$Z_score))
  colnames(overhang_output) <- dicer_overhangs$shift
  print(overhang_output)
  overhang_output <- overhang_output %>% dplyr::mutate(locus = prefix)
  overhang_output <- overhang_output[, c(10, 1:9)]

   suppressWarnings(
      if(!file.exists(paste0(wkdir, "siRNA_dicerz.txt"))){
         utils::write.table(overhang_output, file = paste0(wkdir, "siRNA_dicerz.txt"), sep = "\t", quote = FALSE, append = FALSE, col.names = TRUE, na = "NA", row.names = F)
      } else {
         utils::write.table(overhang_output, file = paste0(wkdir, "siRNA_dicerz.txt"), quote = FALSE, sep = "\t", col.names = FALSE, append = TRUE, na = "NA", row.names = F)
      }
   )


  #heat output nees to be a matrix, so transform
  heat_output <- t(c(prefix, as.vector(results)))

  suppressWarnings(
    if (!file.exists(paste0(wkdir, "siRNA_heatmap.txt"))) {
      utils::write.table(heat_output, file = paste0(wkdir, "siRNA_heatmap.txt"), sep = "\t", quote = FALSE, append = FALSE, col.names = TRUE, na = "NA", row.names = FALSE)
    } else {
      utils::write.table(heat_output, file = paste0(wkdir, "siRNA_heatmap.txt"), quote = FALSE, sep = "\t", col.names = FALSE, append = TRUE, na = "NA", row.names = FALSE)
    }
  )
  print("heatmap has been written.")

  print("Beginning hairpin function.")
  #run the hairpin function on each strand separately
  dsh <- .dual_strand_hairpin(chrom_name, reg_start, reg_stop, length, 1, genome_file, bam_file, logfile, wkdir, plot_output,
                              path_to_RNAfold, annotate_region, weight_reads, gtf_file, write_fastas, out_type)

  #if there are results then use the heat plot

  #user provides argument plot = T or plot = F
  if (plot_output == TRUE) {
    if (!sum(results) == 0) {
      cat(file = paste0(wkdir, logfile), "plot_si_heat\n", append = TRUE)
      print("Making heatmap.")
      heat_plot <- plot_si_heat(results, chrom_name, reg_start, reg_stop, wkdir, pal = pal)
    }
    #heat_plot <- plot_si_heat(results, chrom_name, reg_start, reg_stop, wkdir, pal = pal)
    cat(file = paste0(wkdir, logfile), "get_read_dist\n", append = TRUE)
    print("Making size_dist plot")
    dist <- .get_weighted_read_dist(forward_dt, reverse_dt)

    cat(file = paste0(wkdir, logfile), "plot_sizes\n", append = TRUE)
    size_plot <- plot_sizes(dist)
    cat(file = paste0(wkdir, logfile), "plot_overhangz\n", append = TRUE)

    dicer_overhangs$zscore <- .calc_zscore(dicer_overhangs$proper_count)
    print("Making dicer plot")
    dicer_plot <- plot_overhangz(dicer_overhangs, "none")

    #order of results in dsh obj =
    #plus_res, minus_res, plus_overhang_plot, minus_overhang_plot, density_plot,
    #arc_plot, gtf_plot, plus_phasedz, minus_phasedz

    ### combine siRNA and hpRNA plots

    if (reg_stop - reg_start < 10000) {
      print("Region is less than 10kb in length")

      plus_hp_overhangs <- dsh[[3]]
      minus_hp_overhangs <- dsh[[4]]

      density_plot <- dsh[[5]]
      arc_plot <- dsh[[6]]

      if (annotate_region == TRUE) {
        print("Annotate_region == TRUE")
        gtf_plot <- dsh[[7]]
        plus_phasedz <- dsh[[8]]
        minus_phasedz <- dsh[[9]]

        #gtf_plot <- plot_gtf(gtf_file, chrom_name, reg_start, reg_stop)

        #if there are results for the heatmap, plot, otherwise omit
        if (!sum(results) == 0) {
          print("Heat map contains results")
          left <- cowplot::plot_grid(arc_plot, gtf_plot, density_plot, size_plot, ggplotify::as.grob(heat_plot), rel_widths = c(0.6,1.1,0.9,0.9,0.4), rel_heights = c(0.7,0.7,0.7,0.7,1.4), ncol = 1, align = "vh", axis = "lrtb")
          right <- cowplot::plot_grid(plus_hp_overhangs, minus_hp_overhangs, plus_phasedz, minus_phasedz, dicer_plot, ncol = 1, align = "vh", axis = "l", rel_widths = c(1,1,1,1,1), rel_heights = c(1,1,1,1,1))
        } else {
          left <- cowplot::plot_grid(arc_plot, gtf_plot, density_plot, size_plot, rel_widths = c(0.6,1.1,0.9,0.9), rel_heights = c(0.7,0.7,0.7,0.7,1.4), ncol = 1, align = "vh", axis = "lrtb")
          right <- cowplot::plot_grid(plus_hp_overhangs, minus_hp_overhangs, plus_phasedz, minus_phasedz, dicer_plot, ncol = 1, align = "vh", axis = "l", rel_widths = c(1,1,1,1,1), rel_heights = c(1,1,1,1,1))

        }

        all_plot <- cowplot::plot_grid(left, NULL, right, ncol = 3, rel_widths = c(0.9, 0.01,0.7), align = "vh", axis = "lrtb")


      } else { #if annotate_region == FALSE
        print("Annotate_region == FALSE")
        plus_phasedz<- dsh[[7]]
        minus_phasedz <- dsh[[8]]

        if (!sum(results) == 0) {
          left <- cowplot::plot_grid(arc_plot, density_plot, size_plot, ggplotify::as.grob(heat_plot), ncol = 1, rel_widths = c(0.6,0.9,0.9,0.4), rel_heights = c(0.7,0.7,0.7,1.4), align = "vh", axis = "lrtb")
          right <- cowplot::plot_grid(plus_hp_overhangs, minus_hp_overhangs, plus_phasedz, minus_phasedz, dicer_plot, ncol = 1, align = "vh", axis = "l", rel_widths = c(1,1,1,1,1), rel_heights = c(1,1,1,1,1))
        } else {
          left <- cowplot::plot_grid(arc_plot, density_plot, size_plot, rel_widths = c(0.6,0.9,0.9), rel_heights = c(0.7,0.7,0.7,1.4), ncol = 1, align = "vh", axis = "lrtb")
          right <- cowplot::plot_grid(plus_hp_overhangs, minus_hp_overhangs, plus_phasedz, minus_phasedz, dicer_plot, ncol = 1, align = "vh", axis = "l", rel_widths = c(1,1,1,1,1), rel_heights = c(1,1,1,1,1))
        }

        all_plot <- cowplot::plot_grid(left, NULL, right, ncol = 3, rel_widths = c(0.9, 0.01,0.7), align = "vh", axis = "lrtb")
      }

      if (out_type == "png" || out_type == "PNG") {
        cat(file = paste0(wkdir, logfile), "Making png\n", append = TRUE)
        grDevices::png(file = file.path(wkdir, paste0(prefix, "_si_plot.png")), height = 16, width = 14, units = "in", res = 300)
        print(all_plot)
        grDevices::dev.off()
      } else {
        cat(file = paste0(wkdir, logfile), "Making pdf\n", append = TRUE)
        grDevices::pdf(file = file.path(wkdir, paste0(prefix, "_si_plot.pdf")), height = 16, width = 14)
        print(all_plot)
        grDevices::dev.off()
      }

    } else { #none of the hairpin plots were made because the region > 10kb
      print("Region is greater than 10kb in length")
      if (annotate_region == TRUE) {
        print("Annotate_region == TRUE")
        #gtf_plot <- dsh[[7]]
        #plus_phasedz <- dsh[[8]]
        #minus_phasedz <- dsh[[9]]

        # things are broken here
        gtf_plot <- plot_gtf(gtf_file, chrom_name, reg_start, reg_stop)


        if (!sum(results) == 0) {
          print("Heatplot results are not empty")
          data <- read_densityBySize(bam_obj, chrom_name, reg_start, reg_stop, bam_file, wkdir)

          #density_plot <- plot_density(data, reg_start, reg_stop)
          density_plot <- plot_large_density(data, reg_start, reg_stop)
          top <- cowplot::plot_grid(density_plot, size_plot, ncol = 2,  rel_widths = c(1,1), rel_heights = c(1,1), align = "vh", axis = "lrtb")
          bottom <- cowplot::plot_grid(ggplotify::as.grob(heat_plot), dicer_plot, ncol = 2, rel_widths = c(1,1), rel_heights = c(1,1), align = "vh", axis = "lrtb")
          all_plot <- cowplot::plot_grid(top, NULL, bottom, ncol = 1, rel_widths = c(1,1,1), rel_heights = c(1,0.1,1), align = "vh", axis = "lrtb")

        } else { #if no heat plot
           print("Heat plot results are empty")
           data <- read_densityBySize(bam_obj, chrom_name, reg_start, reg_stop, bam_file, wkdir)

           density_plot <- plot_large_density(data, reg_start, reg_stop)
           top <- cowplot::plot_grid(density_plot, size_plot, ncol = 2, rel_widths = c(1,1), rel_heights = c(1,1), align = "vh", axis = "lrtb")
           bottom <- cowplot::plot_grid(dicer_plot, rel_widths = c(1), rel_heights = c(1), align = "vh", axis = "lrtb")
           all_plot <- cowplot::plot_grid(top, NULL, bottom, ncol = 1, rel_widths = c(1,1,1), rel_heights = c(1,0.1,1), align = "vh", axis = "lrtb")
        }

        if (out_type == "png" || out_type == "PNG") {
          cat(file = paste0(wkdir, logfile), "Making png\n", append = TRUE)
          grDevices::png(file = file.path(wkdir, paste0(prefix, "_si_plot.png")), height = 9, width = 14, units = "in", res = 300)
          print(all_plot)
          grDevices::dev.off()
        } else {
          cat(file = paste0(wkdir, logfile), "Making pdf\n", append = TRUE)
          grDevices::pdf(file = file.path(wkdir, paste0(prefix, "_si_plot.pdf")), height = 9, width = 14)
          print(all_plot)
          grDevices::dev.off()
        }
      } else { # if annotate_region == FALSE
        #plus_phasedz<- dsh[[7]]
        #minus_phasedz <- dsh[[8]]
        data <- read_densityBySize(bam_obj, chrom_name, reg_start, reg_stop, bam_file, wkdir)
        density_plot <- plot_large_density(data, reg_start, reg_stop)

        if (!sum(results) == 0) {
          left <- cowplot::plot_grid(density_plot, size_plot, dicer_plot, ncol = 1, rel_widths = c(1,1), rel_heights = c(1,1), align = "vh", axis = "lrtb")
          right <- cowplot::plot_grid(NULL, ggplotify::as.grob(heat_plot), NULL, ncol = 1, align = "vh", axis = "l", rel_widths = c(0.4,1,0.4), rel_heights = c(1,1))
          all_plot <- cowplot::plot_grid(left, NULL, right, ncol = 3, rel_widths = c(1,0.1,0.8), rel_heights = c(1,1,1), align = "vh", axis = "lrtb")

        } else {
          bottom <- cowplot::plot_grid(dicer_plot, size_plot, rel_widths = c(1,1), rel_heights = c(1,1), ncol = 2, align = "vh", axis = "lrtb")
          top <- cowplot::plot_grid(density_plot, ncol = 1, align = "vh", axis = "l", rel_widths = c(1), rel_heights = c(1))
          all_plot <- cowplot::plot_grid(top, NULL, bottom, ncol = 1, rel_widths = c(1,1), rel_heights = c(1,0.1,1), align = "vh", axis = "lrtb")
        }
      }

      if (out_type == "png" || out_type == "PNG") {
        cat(file = paste0(wkdir, logfile), "Making png\n", append = TRUE)
        grDevices::png(file = file.path(wkdir, paste0(prefix, "_si_plot.png")), height = 13, width = 13, units = "in", res = 300)
        print(all_plot)
        grDevices::dev.off()
      } else {
        cat(file = paste0(wkdir, logfile), "Making pdf\n", append = TRUE)
        grDevices::pdf(file = file.path(wkdir, paste0(prefix, "_si_plot.pdf")), height = 13, width = 13)
        print(all_plot)
        grDevices::dev.off()
      }
    }
  }

                                 #1         #2       #3             #4         #5
      #left <- cowplot::plot_grid(arc_plot, gtf_plot, density_plot, size_plot, heat, rel_widths = c(0.6,1.1,0.9,0.9,0.4), rel_heights = c(0.7,0.7,0.7,0.7,1.4), ncol = 1, align = "vh", axis = "lrtb")
      #right <- cowplot::plot_grid(plus_hp_overhangs, minus_hp_overhangs, plus_phasedz, minus_phasedz, dicer_plot, ncol = 1, align = "vh", axis = "l", rel_widths = c(1,1,1,1,1), rel_heights = c(1,1,1,1,1))

   return(list(heat = results, si_dicer = dicer_overhangs, dsh))
}

