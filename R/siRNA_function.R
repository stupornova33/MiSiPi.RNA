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
   prefix <- paste0(chrom_name, "_", reg_start, "_", reg_stop)
   width <- pos <- phased_dist <- phased_num <- phased_z <- phased_dist2 <- plus_num2 <- phased_dist1 <- phased_num1 <- NULL
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
      dplyr::mutate(start = pos, end = pos + width - 1) %>%
      dplyr::select(-c(pos)) %>%
      dplyr::group_by_all() %>%
      dplyr::reframe(count = dplyr::n())


   cat(file = paste0(dir, logfile), "Making Reverse DT\n", append = TRUE)
   reverse_dt <- data.table::setDT(makeBamDF(chromM)) %>%
       subset(width <= 32 & width >= 15) %>%
       dplyr::mutate(start = pos, end = pos + width - 1) %>%
       dplyr::select(-c(pos)) %>%
       dplyr::group_by_all() %>%
       dplyr::summarize(count = dplyr::n())

   #chromP <- NULL
   #chromM <- NULL


   if(nrow(forward_dt) > 0 & nrow(reverse_dt) > 0){
      print("f_dt and r_dt are not empty")
      cat(file = paste0(dir, logfile), "Calc overhangs\n", append = TRUE)
      forward_dt <- get_top_n_weighted(forward_dt, chrom_name, 10)
      reverse_dt <- get_top_n_weighted(reverse_dt, chrom_name, 10)

      #check to see if subsetted dfs are empty
      if(nrow(forward_dt) > 0 & nrow(reverse_dt) > 0){
      ### added, testing
      overlaps <- find_overlaps(forward_dt, reverse_dt)
      #dicer_overhangs <- calc_overhangs(overlaps$r1_start, overlaps$r1_end, overlaps$r2_start, overlaps$r2_width)
      dicer_overhangs <- calc_expand_overhangs(overlaps$r1_start, overlaps$r1_end, overlaps$r2_start, overlaps$r2_width)
      dicer_overhangs$Z_score <- calc_zscore(dicer_overhangs$proper_count)
      ###
      cat(file = paste0(dir, logfile), "get_si_overlaps\n", append = TRUE)
      results <- get_si_overlaps(forward_dt$start, forward_dt$end, forward_dt$width,
                              reverse_dt$start, reverse_dt$end, reverse_dt$width)

      row.names(results) <- c('15','','17','','19','','21','','23','','25','','27','','29','','31','')
      colnames(results) <- c('15','','17','','19','','21','','23','','25','','27','','29','','31','')
      } else {
        cat(file = paste0(dir, logfile), "No reads detected on one strand. \n", append = TRUE)
        dicer_overhangs <- data.frame(shift = seq(-8,8), proper_count = c(rep(0, times = 17)), Z_score = c(rep(-33, times = 17)))
        results <- 0
      }

   } else {
      cat(file = paste0(dir, logfile), "No reads detected on one strand. \n", append = TRUE)
      dicer_overhangs <- data.frame(shift = seq(-8,8), proper_count = c(rep(0, times = 17)), Z_score = c(rep(-33, times = 17)))
      results <- 0
   }


   overhang_output <- data.frame(t(dicer_overhangs$Z_score))
   colnames(overhang_output) <- dicer_overhangs$shift
   print(overhang_output)
   overhang_output <- overhang_output %>% dplyr::mutate(locus = prefix)
   overhang_output <- overhang_output[, c(18, 1:17)]

   suppressWarnings(
      if(!file.exists("siRNA_dicerz.txt")){
         utils::write.table(overhang_output, file = paste0(dir, "siRNA_dicerz.txt"), sep = "\t", quote = FALSE, append = T, col.names = T, na = "NA", row.names = F)
      } else {
         utils::write.table(overhang_output, file = paste0(dir, "siRNA_dicerz.txt"), quote = FALSE, sep = "\t", col.names = F, append = TRUE, na = "NA", row.names = F)
      }
   )


   heat_output <- t(c(prefix, as.vector(results)))

   suppressWarnings(
      if(!file.exists("siRNA_heatmap.txt")){
         utils::write.table(heat_output, file = paste0(dir, "siRNA_heatmap.txt"), sep = "\t", quote = FALSE, append = T, col.names = F, na = "NA", row.names = F)
      } else {
         utils::write.table(heat_output, file = paste0(dir, "siRNA_heatmap.txt"), quote = FALSE, sep = "\t", col.names = F, append = TRUE, na = "NA", row.names = F)
      }
   )


 if(!sum(results) == 0){
   ### phasing
   plus_phased_counts <- calc_phasing(forward_dt, forward_dt)
   minus_phased_counts <- calc_phasing(reverse_dt, reverse_dt)

   plus_phased_counts <- plus_phased_counts %>% dplyr::rename('phased_dist2' = phased_dist, 'plus_num2' = phased_num, 'phased_z2' = phased_z)
   minus_phased_counts <- minus_phased_counts %>% dplyr::rename('phased_dist1' = phased_dist, 'phased_num1' = phased_num, 'phased_z1' = phased_z)

   df <- cbind(minus_phased_counts, plus_phased_counts)

   ## write phased results to file
   plus_phased_output <- plus_phased_counts %>% dplyr::select(-c(phased_dist2, plus_num2))
   plus_phased_output <- t(c(prefix, t(plus_phased_output)))
   minus_phased_output <- minus_phased_counts %>% dplyr::select(-c(phased_dist1, phased_num1))
   minus_phased_output <- t(c(prefix, t(minus_phased_output)))

   suppressWarnings(
      if(!file.exists("siRNA_plus_phasedz.txt")){
         utils::write.table(plus_phased_output, file = paste0(dir, "siRNA_plus_phasedz.txt"), sep = "\t", quote = FALSE, append = T, col.names = F, na = "NA", row.names = F)
      } else {
         utils::write.table(plus_phased_output, file = paste0(dir, "siRNA_plus_phasedz.txt"), quote = FALSE, sep = "\t", col.names = F, append = TRUE, na = "NA", row.names = F)
      }
   )

   suppressWarnings(
      if(!file.exists("siRNA_minus_phasedz.txt")){
         utils::write.table(minus_phased_output, file = paste0(dir, "siRNA_minus_phasedz.txt"), sep = "\t", quote = FALSE, append = T, col.names = F, na = "NA", row.names = F)
      } else {
         utils::write.table(minus_phased_output, file = paste0(dir, "siRNA_minus_phasedz.txt"), quote = FALSE, sep = "\t", col.names = F, append = TRUE, na = "NA", row.names = F)
      }
   )

   dsh <- dual_strand_hairpin(chrom_name, reg_start, reg_stop, length, 1, genome_file, input_file, logfile, dir, plot_output,
                              path_to_RNAfold, annotate_bed, gff_file)
   ### hairpins
   if(plot_output == "T"){
      cat(file = paste0(dir, logfile), "plot_si_heat\n", append = TRUE)

      heat_plot <- plot_si_heat(results, chrom_name, reg_start, reg_stop, dir, pal = pal)
      cat(file = paste0(dir, logfile), "get_read_dist\n", append = TRUE)
      dist <- get_read_dist(bam_obj, chrom_name, reg_start, reg_stop)
      #forward_dt <- NULL
      #reverse_dt <- NULL
      cat(file = paste0(dir, logfile), "plot_sizes\n", append = TRUE)
      size_plot <- plot_sizes(dist)
      cat(file = paste0(dir, logfile), "plot_overhangz\n", append = TRUE)

      dicer_plot <- plot_overhangz(dicer_overhangs)

      phased_plot <- plot_phasedz(df)
      #will output the hairpin plots to file "_hairpin_fold.pdf"


      ### make siRNA plots

      top <- cowplot::plot_grid(size_plot, NULL, ncol = 3,dicer_plot, rel_widths = c(1,0.2,1), align = "vh", axis = "lrtb")

      bottom <- cowplot::plot_grid(ggplotify::as.grob(heat_plot), phased_plot, ncol = 2, align = "vh", axis = "lrtb")

      all_plot <- cowplot::plot_grid(top, bottom, ncol = 1, rel_widths = c(1, 1), rel_heights = c(1,1))

      cat(file = paste0(dir, logfile), "Making PDF\n", append = TRUE)
      grDevices::pdf(file = paste0(dir, chrom_name, "_", reg_start, "_", reg_stop, "_si_plot.pdf"), height = 9, width = 9)
      print(all_plot)
      grDevices::dev.off()
     }

   } else {
      dsh <- null_hp_res()
      forward_dt <- NULL
      reverse_dt <- NULL
      cat(file = paste0(dir, logfile), "Zero overlapping reads found.\n", append = TRUE)
      results <- NA
   }

     return(list(heat = results, si_dicer = dicer_overhangs, dsh))
}
