#' plot_sizes_only
#' A function which will process the reads at a locus but only plot the size distribution and will not run the miRNA, siRNA or piRNA modules.
#' @param bam_file A string - the name of the BAM file.
#' @param bed_file A string - the name of the 3-column BED file.
#' @return results

#' @export

plot_sizes_only <- function(bam_file, bed_file){
  
  roi <- read.csv(bed_file, header = FALSE, sep = '\t')
  
  for(i in 1:nrow(roi)){
    
    chrom_name <- roi$V1[i]
    reg_start <- roi$V2[i]
    reg_stop <- roi$V3[i]
    logfile <- "plot_sizes.log.txt"
    print(paste0(chrom_name, "_", reg_start, "_", reg_stop))
    prefix <- paste0(chrom_name, "_", reg_start, "_", reg_stop)

  # use Rsamtools to process the bam file
    bam_obj <- OpenBamFile(bam_file, logfile)
    bam_header <- Rsamtools::scanBamHeader(bam_obj)
    chr_name <- names(bam_header[['targets']])
    chr_length <- unname(bam_header[['targets']])
    bam_header <- NULL

  # extract reads by strand
  # this creates a list object
   chromP <- getChrPlus(bam_obj, chrom_name, reg_start, reg_stop)
   chromM <- getChrMinus(bam_obj, chrom_name, reg_start, reg_stop)
  
  
   forward_dt <- data.table::setDT(make_si_BamDF(chromP)) %>%
     subset(width <= 32 & width >= 18) %>%
     dplyr::mutate(start = pos, end = pos + width - 1) %>%
     dplyr::select(-c(pos)) %>%
     dplyr::group_by_all() %>%
     # get the number of times a read occurs
     dplyr::summarize(count = dplyr::n())
  
   reverse_dt <- data.table::setDT(make_si_BamDF(chromM)) %>%
     subset(width <= 32 & width >= 18) %>%
     dplyr::mutate(start = pos, end = pos + width - 1) %>%
     dplyr::select(-c(pos)) %>%
     dplyr::group_by_all() %>%
     dplyr::summarize(count = dplyr::n())
  
   dist <- get_weighted_read_dist(forward_dt, reverse_dt)
  
   size_plot <- plot_sizes(dist)
  
   grDevices::pdf(file = paste0(chrom_name, "_", reg_start, "_", reg_stop, "_sizes.pdf"), height = 4, width = 4)
   print(size_plot)
   grDevices::dev.off() 
  }
  
}
