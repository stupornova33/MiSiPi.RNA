#' siRNA function
#' process forward and reverse dt and plot si heatmap
#' @param chrom_name a string
#' @param reg_start a whole number
#' @param reg_stop a whole number
#' @param length an integer
#' @param min_read_count an integer
#' @param genome_file a string
#' @param input_file a string
#' @param logfile a string
#' @param dir a string 
#' @param pal a string
#' @param plot_output a string, 'T' or 'F'. Default is 'T'
#' @param path_to_RNAfold a string
#' @param annotate_bed a string, "T" or "F"
#' @param gff_file a string
#' @return results

#' @export
siRNA_function <- function(chrom_name, reg_start, reg_stop, length, min_read_count, genome_file, input_file, logfile, dir, 
                           pal, plot_output, path_to_RNAfold, annotate_bed, gff_file){
   print(paste0(chrom_name, "_", reg_start, "_", reg_stop))
   print(paste0('Length: ', reg_stop - reg_start))
   bam_obj <- OpenBamFile(input_file, logfile)
   bam_header <- Rsamtools::scanBamHeader(bam_obj)
   chr_name <- names(bam_header[['targets']])
   chr_length <- unname(bam_header[['targets']])
   bam_header <- NULL
   
   cat(file = paste0(dir, logfile), paste0("chrom_name: ", chrom_name, " reg_start: ", reg_start, " reg_stop: ", reg_stop, "\n"), append = TRUE)
   cat(file = paste0(dir, logfile), "Filtering forward and reverse reads by length\n", append = TRUE)
   
   chromP <- getChrPlus(bam_obj, chrom_name, reg_start, reg_stop)
   chromM <- getChrMinus(bam_obj, chrom_name, reg_start, reg_stop)
   
   cat(file = paste0(dir, logfile), "Making Forward DT\n", append = TRUE)
   forward_dt <- data.table::setDT(makeBamDF(chromP)) %>%
      subset(width <= 32 & width >= 15) %>% 
      dplyr::mutate(start = pos, end = pos + width - 1) %>% dplyr::distinct() %>%
      dplyr::select(-c(pos))
   
   cat(file = paste0(dir, logfile), "Making Reverse DT\n", append = TRUE)
   reverse_dt <- data.table::setDT(makeBamDF(chromM)) %>%
      subset(width <= 32 & width >= 15) %>%
      dplyr::mutate(start = pos, end = pos + width - 1) %>% dplyr::distinct() %>%
      dplyr::select(-c(pos))
   
   #chromP <- NULL
   #chromM <- NULL
   

   if(!nrow(forward_dt) == 0 && !nrow(reverse_dt) == 0){
      cat(file = paste0(dir, logfile), "Calc overhangs\n", append = TRUE)
      
      ### added, testing
      overlaps <- find_overlaps(forward_dt, reverse_dt)
      dicer_overhangs <- calc_overhangs(overlaps$r1_start, overlaps$r1_end, overlaps$r2_start, overlaps$r2_width)
      ### 
      cat(file = paste0(dir, logfile), "get_si_overlaps\n", append = TRUE)
      results <- get_si_overlaps(forward_dt$start, forward_dt$end, forward_dt$width,
                              reverse_dt$start, reverse_dt$end, reverse_dt$width)
      
       row.names(results) <- c('15','','17','','19','','21','','23','','25','','27','','29','','31','')
       colnames(results) <- c('15','','17','','19','','21','','23','','25','','27','','29','','31','')

   } else {
      cat(file = paste0(dir, logfile), "No reads detected on one strand. \n", append = TRUE)
      dicer_overhangs <- data.frame(pos_shift = seq(-4,4), proper_count = c(rep(0, times = 9)), Z_score = c(rep(NA, times = 9)))
      results <- 0
   }
   

 if(!sum(results) == 0){
   
   
   ### phasing
   plus_phased_counts <- calc_phasing(forward_dt, forward_dt)
   minus_phased_counts <- calc_phasing(reverse_dt, reverse_dt)
   
   plus_phased_counts <- plus_phased_counts %>% dplyr::rename('phased_dist2' = phased_dist, 'plus_num2' = phased_num, 'phased_z2' = phased_z)
   minus_phased_counts <- minus_phased_counts %>% dplyr::rename('phased_dist1' = phased_dist, 'phased_num1' = phased_num, 'phased_z1' = phased_z)
   
   df <- cbind(minus_phased_counts, plus_phased_counts)
   
   ### hairpins
   if(plot_output == "T"){
      cat(file = paste0(dir, logfile), "plot_si_heat\n", append = TRUE)
      
      heat_plot <- plot_si_heat(results, chrom_name, reg_start, reg_stop, dir, pal = pal)
      cat(file = paste0(dir, logfile), "get_read_dist\n", append = TRUE)
      dist <- get_read_dist(chromM, chromP)
      #forward_dt <- NULL
      #reverse_dt <- NULL
      cat(file = paste0(dir, logfile), "plot_sizes\n", append = TRUE)
      size_plot <- plot_sizes(dist)
      cat(file = paste0(dir, logfile), "plot_overhangz\n", append = TRUE)
      
      dicer_plot <- plot_overhangz(dicer_overhangs)
      print('making all_plot')
         
      phased_plot <- plot_phasedz(df)   
      #will output the hairpin plots to file "_hairpin_fold.pdf"
      #mapply(dual_strand_hairpin, chrom_name, reg_start, reg_stop, length, 1, genome_file, input_file, logfile, dir, plot_output, 
      #       path_to_RNAfold, annotate_bed, gff_file)
      dsh <- dual_strand_hairpin(chrom_name, reg_start, reg_stop, length, 1, genome_file, input_file, logfile, dir, plot_output, 
                                 path_to_RNAfold, annotate_bed, gff_file)
     
      ### make siRNA plots
      #top <- cowplot::plot_grid(ggplotify::as.grob(heat_plot), NULL, phased_plot, ncol = 3, rel_widths = c(1,0.2, 1), align = "vh", axis = "lrtb")
      #bottom <- cowplot::plot_grid(size_plot, NULL, dicer_plot, rel_widths = c(1,0.1,1.3), nrow = 1, ncol = 3, align = "vh", axis = "lrtb")
      top <- cowplot::plot_grid(size_plot, NULL, ncol = 3,dicer_plot, rel_widths = c(1,0.2,1), align = "vh", axis = "lrtb")
      print("Make the bottom")
      bottom <- cowplot::plot_grid(ggplotify::as.grob(heat_plot), phased_plot, ncol = 2, align = "vh", axis = "lrtb")
      print("Make all plot")
      all_plot <- cowplot::plot_grid(top, bottom, ncol = 1, rel_widths = c(1, 1), rel_heights = c(1,1))
      print("Make pdf")
      cat(file = paste0(dir, logfile), "Making PDF\n", append = TRUE)
      pdf(file = paste0(dir, chrom_name, "_", reg_start, "_", reg_stop, "_si_plot.pdf"), height = 9, width = 9)
      print(all_plot)
      dev.off()
   }
      
   } else {
      forward_dt <- NULL
      reverse_dt <- NULL
      cat(file = paste0(dir, logfile), "Zero overlapping reads found.\n", append = TRUE)
      results <- NA
   }
   
     return(list(results, dicer_overhangs))
}
