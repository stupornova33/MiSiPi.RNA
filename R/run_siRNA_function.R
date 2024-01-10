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
#' @param plot_output a string, 'T' or 'F'. Default is 'T'
#' @param path_to_RNAfold a string
#' @param annotate_bed a string, "T" or "F"
#' @param weight_reads a string, "T" or "F"
#' @param bed_file a string
#' @return results

#' @export

run_siRNA_function <- function(chrom_name, reg_start, reg_stop, length, min_read_count, genome_file, bam_file, logfile, wkdir,
                           pal, plot_output, path_to_RNAfold, annotate_bed, weight_reads, bed_file){
   print(paste0(chrom_name, "_", reg_start, "_", reg_stop))
   prefix <- paste0(chrom_name, "_", reg_start, "_", reg_stop)
   width <- pos <- phased_dist <- phased_num <- phased_z <- phased_dist2 <- plus_num2 <- phased_dist1 <- phased_num1 <- NULL

   # use Rsamtools to process the bam file
   bam_obj <- OpenBamFile(bam_file, logfile)
   bam_header <- Rsamtools::scanBamHeader(bam_obj)
   chr_name <- names(bam_header[['targets']])
   chr_length <- unname(bam_header[['targets']])
   bam_header <- NULL

   cat(file = paste0(wkdir, logfile), paste0("chrom_name: ", chrom_name, " reg_start: ", reg_start, " reg_stop: ", reg_stop, "\n"), append = TRUE)
   cat(file = paste0(wkdir, logfile), "Filtering forward and reverse reads by length\n", append = TRUE)

   # extract reads by strand
   # this creates a list object
   chromP <- getChrPlus(bam_obj, chrom_name, reg_start, reg_stop)
   chromM <- getChrMinus(bam_obj, chrom_name, reg_start, reg_stop)

   # turn the list object into a more useable data frame and filter reads by length,
   # bam only contains pos and width, need to add an end column
   cat(file = paste0(wkdir, logfile), "Making Forward DT\n", append = TRUE)
   forward_dt <- data.table::setDT(make_si_BamDF(chromP)) %>%
     subset(width <= 32 & width >= 18) %>%
      dplyr::mutate(start = pos, end = pos + width - 1) %>%
      dplyr::select(-c(pos)) %>%
      dplyr::group_by_all() %>%
     # get the number of times a read occurs
      dplyr::summarize(count = dplyr::n())


   cat(file = paste0(wkdir, logfile), "Making Reverse DT\n", append = TRUE)
   reverse_dt <- data.table::setDT(make_si_BamDF(chromM)) %>%
       subset(width <= 32 & width >= 18) %>%
       dplyr::mutate(start = pos, end = pos + width - 1) %>%
       dplyr::select(-c(pos)) %>%
       dplyr::group_by_all() %>%
       dplyr::summarize(count = dplyr::n())


   chromP <- NULL
   chromM <- NULL

   # If the data frames are empty there are no reads, can't do siRNA calculations
   if(nrow(forward_dt) > 0 & nrow(reverse_dt) > 0){
      print("f_dt and r_dt are not empty")
      cat(file = paste0(wkdir, logfile), "Calc overhangs\n", append = TRUE)

      #include "T" argument to return read sequences
      if(weight_reads == "Top"){
        forward_dt <- get_top_n_weighted(forward_dt, chrom_name, "T")
        reverse_dt <- get_top_n_weighted(reverse_dt, chrom_name, "T")

        print("Completed getting weighted dataframes.")
      } else if(weight_reads == "Locus_norm"){

        forward_dt <- locus_norm(forward_dt, sum(forward_dt$count, reverse_dt$count), "T")
        reverse_dt <- locus_norm(reverse_dt, sum(reverse_dt$count, reverse_dt$count), "T")

      } else {
        forward_dt <- get_top_n(forward_dt, chrom_name, "T")
        reverse_dt <- get_top_n(reverse_dt, chrom_name, "T")
      }
      # check to see if subsetted dfs are empty
      # have to keep doing this at each step otherwise errors will happen
      if(nrow(forward_dt) > 0 & nrow(reverse_dt) > 0){

      # get overlapping reads
      overlaps <- find_overlaps(forward_dt, reverse_dt) %>% dplyr::mutate(p5_overhang = r2_end - r1_end, p3_overhang = r2_start - r1_start) %>%
        dplyr::filter(p5_overhang >= 0 & p3_overhang >= 0)

      # create a GRanges object to extract the chromosome sequence from the genome file
      mygranges <- GenomicRanges::GRanges(
        seqnames = c(chrom_name),
        ranges = IRanges::IRanges(start=c(1), end=c(length)))

      geno_seq <- Rsamtools::scanFa(genome_file, mygranges)

      #get sequence for only the region of interest from the chromosome
      geno_seq <- as.character(unlist(Biostrings::subseq(geno_seq, start = 1, end = length)))

      #check to see how many overlaps have the 2nt overhang
      #currently this extracts pairs with overhang at 3p OR 5p
      #tried getting pairs with overhang at 3p AND 5p but few/no results returned in many cases...
      proper_overlaps <- overlaps %>% dplyr::filter(r2_start - r1_start == 2 | r2_end - r1_end == 2)

      #write the pairs to file
      if(nrow(proper_overlaps > 0)){
        paired_seqs <- proper_overlaps %>%
         dplyr::mutate(r1_seq = paste0(">",chrom_name, ":", proper_overlaps$r1_start, "-", proper_overlaps$r1_end, " ", substr(geno_seq, proper_overlaps$r1_start, proper_overlaps$r1_end)), r2_seq = paste0(">",chrom_name, ":", proper_overlaps$r2_start, "-", proper_overlaps$r2_end, " " , substr(geno_seq, proper_overlaps$r2_start, proper_overlaps$r2_end))) %>%
          dplyr::distinct()

        paired_seqs <- paired_seqs %>% dplyr::transmute(col1 = paste0(r1_seq, ",", r2_seq)) %>% tidyr::separate_rows(col1, sep = ",")
        fastas <- stringi::stri_split_regex(paired_seqs$col1, " ")

        write.table(unlist(fastas), paste0(wkdir, "siRNA_pairs.fa"), sep = " ", quote = FALSE, row.names = FALSE, col.names = FALSE)
      } else {
        cat(file = paste0(wkdir, logfile), "No proper siRNA pairs were found.\n", append = TRUE)
      }


      #calculate the number of dicer pairs for the zscore
      dicer_overhangs <- calc_overhangs(overlaps$r1_start, overlaps$r1_end, overlaps$r2_start, overlaps$r2_width)
      dicer_overhangs$Z_score <- calc_zscore(dicer_overhangs$proper_count)

      cat(file = paste0(wkdir, logfile), "get_si_overlaps\n", append = TRUE)
      # calculate the siRNA pairs for the heatmap
      results <- get_si_overlaps(forward_dt$start, forward_dt$end, forward_dt$width,
                              reverse_dt$start, reverse_dt$end, reverse_dt$width)

      row.names(results) <- c('15','','17','','19','','21','','23','','25','','27','','29','','31','')
      colnames(results) <- c('15','','17','','19','','21','','23','','25','','27','','29','','31','')
      } else {


        # results are being stored also in case the run_all function is being used, at the end they will be written to a table
        cat(file = paste0(wkdir, logfile), "No reads detected on one strand. \n", append = TRUE)
        #the data.frame should be modified if using calc_expand_overhangs
        dicer_overhangs <- data.frame(shift = seq(-4,4), proper_count = c(rep(0, times = 9)), Z_score = c(rep(-33, times = 9)))
        #dicer_overhangs <- data.frame(shift = seq(-8,8), proper_count = c(rep(0, times = 17)), Z_score = c(rep(-33, times = 17)))
        results <- 0
      }

   } else {
      cat(file = paste0(wkdir, logfile), "No reads detected on one strand. \n", append = TRUE)
      #dicer_overhangs <- data.frame(shift = seq(-8,8), proper_count = c(rep(0, times = 17)), Z_score = c(rep(-33, times = 17)))
      dicer_overhangs <- data.frame(shift = seq(-4,4), proper_count = c(rep(0, times = 9)), Z_score = c(rep(-33, times = 9)))
      results <- 0
   }


   # transform the data frame for writing to table by row
   # output is the locus followed by all zscores
   overhang_output <- data.frame(t(dicer_overhangs$Z_score))
   colnames(overhang_output) <- dicer_overhangs$shift
   print(overhang_output)
   overhang_output <- overhang_output %>% dplyr::mutate(locus = prefix)
   overhang_output <- overhang_output[, c(10, 1:9)]

   suppressWarnings(
      if(!file.exists("siRNA_dicerz.txt")){
         utils::write.table(overhang_output, file = paste0(wkdir, "siRNA_dicerz.txt"), sep = "\t", quote = FALSE, append = T, col.names = T, na = "NA", row.names = F)
      } else {
         utils::write.table(overhang_output, file = paste0(wkdir, "siRNA_dicerz.txt"), quote = FALSE, sep = "\t", col.names = F, append = TRUE, na = "NA", row.names = F)
      }
   )


   #heat output nees to be a matrix, so transform
   heat_output <- t(c(prefix, as.vector(results)))

   suppressWarnings(
      if(!file.exists("siRNA_heatmap.txt")){
         utils::write.table(heat_output, file = paste0(wkdir, "siRNA_heatmap.txt"), sep = "\t", quote = FALSE, append = T, col.names = F, na = "NA", row.names = F)
      } else {
         utils::write.table(heat_output, file = paste0(wkdir, "siRNA_heatmap.txt"), quote = FALSE, sep = "\t", col.names = F, append = TRUE, na = "NA", row.names = F)
      }
   )


 #if there are results then do the hairpin functions and plots
 if(!sum(results) == 0){
   ### hairpin function
   print("Beginning hairpin function.")


   #run the hairpin function on each strand separately
   dsh <- dual_strand_hairpin(chrom_name, reg_start, reg_stop, length, 1, genome_file, bam_file, logfile, wkdir, plot_output,
                              path_to_RNAfold, annotate_bed, weight_reads, bed_file)

   #user provides argument plot = T or plot = F
   if(plot_output == "T"){
      cat(file = paste0(wkdir, logfile), "plot_si_heat\n", append = TRUE)

      heat_plot <- plot_si_heat(results, chrom_name, reg_start, reg_stop, wkdir, pal = pal)
      cat(file = paste0(wkdir, logfile), "get_read_dist\n", append = TRUE)

      dist <- get_weighted_read_dist(forward_dt, reverse_dt)

      cat(file = paste0(wkdir, logfile), "plot_sizes\n", append = TRUE)
      size_plot <- plot_sizes(dist)
      cat(file = paste0(wkdir, logfile), "plot_overhangz\n", append = TRUE)

      dicer_overhangs$zscore <- calc_zscore(dicer_overhangs$proper_count)
      dicer_plot <- plot_overhangz(dicer_overhangs, "none")


      ### make siRNA plots

      top <- cowplot::plot_grid(size_plot, NULL, ncol = 3,dicer_plot, rel_widths = c(1,0.2,1), align = "vh", axis = "lrtb")

      bottom <- cowplot::plot_grid(ggplotify::as.grob(heat_plot), NULL, ncol = 2, align = "vh", axis = "lrtb", rel_widths = c(1,1))

      all_plot <- cowplot::plot_grid(top, bottom, ncol = 1, rel_widths = c(1, 1), rel_heights = c(1,1))

      cat(file = paste0(wkdir, logfile), "Making PDF\n", append = TRUE)
      grDevices::pdf(file = paste0(wkdir, chrom_name, "_", reg_start, "_", reg_stop, "_si_plot.pdf"), height = 9, width = 9)
      print(all_plot)
      grDevices::dev.off()
     }

   } else {

      dsh <- null_hp_res()
      forward_dt <- NULL
      reverse_dt <- NULL
      cat(file = paste0(wkdir, logfile), "Zero overlapping reads found.\n", append = TRUE)
      results <- NA
   }

     return(list(heat = results, si_dicer = dicer_overhangs, dsh))
}

