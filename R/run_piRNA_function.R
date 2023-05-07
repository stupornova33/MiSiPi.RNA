#' run_piRNA_function
#' takes test_list and input_file
#' returns pi_res and plots
#' 
#' @param bed_file  a tab-delineated file with three columns. Chrom, start, stop
#' @param input_file a BAM file
#' @param pal a string
#' @param plot_output a string, 'T' or 'F', default = 'T
#' @return pi_res, plots

#' @export


run_piRNA_function <- function(bed_file, input_file, pal, plot_output){
   #`%>%` <- magrittr::`%>%`
   #test_list <- read.csv(bed_file, sep = "\t", header = FALSE) %>% 
   #   dplyr::mutate(V2 = V2 + 1, V3 = V3 + 1)
   # add siRNA heat for hairpins
   
   dir <- 'piRNA_outputs/'
   if(!dir.exists(dir) == TRUE) dir.create(dir)
   
   logfile = "piRNA_logfile.txt"
   if(!file.exists(logfile) == TRUE) file.create(paste0(dir, logfile))
   
   #bam_obj <- OpenBamFile(input_file)
   #bam_header <- Rsamtools::scanBamHeader(bam_obj)
   #chr_name <- names(bam_header[['targets']])
   #chr_length <- unname(bam_header[['targets']])
   #bam_header <- NULL
   
   mapply(piRNA_function, vars[[1]], vars[[2]], vars[[3]], vars[[10]], logfile, dir, vars[[12]], vars[[6]])
   #mapply(piRNA_function, test_list$V1, test_list$V2, test_list$V3, input_file, logfile, dir, pal, plot_output)


}