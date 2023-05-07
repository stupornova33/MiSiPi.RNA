#' run_phased_piRNA
#'
#' @param vars a list
#' @return plots

#' @export

run_phased_piRNA <- function(bed_file, input_file, min_read_count, plot_output = 'T'){
   #`%>%` <- magrittr::`%>%`
   #test_list <- read.csv(bed_file, sep = "\t", header = FALSE) %>% 
   #   dplyr::mutate(V2 = V2 + 1, V3 = V3 + 1)
   #input_file <- input_file
   #min_read_count <- min_read_count

   
   logfile = "phased_piRNA_logfile.txt"
   #bam_obj <- OpenBamFile(input_file)
   #bam_header <- Rsamtools::scanBamHeader(bam_obj)
   #chr_name <- names(bam_header[['targets']])
   #chr_length <- unname(bam_header[['targets']])
   #bam_header <- NULL
   
   
   dir <- 'phased_piRNA_outputs/'
   if(!dir.exists(dir) == TRUE) dir.create(dir)
   
   logfile = "phased_piRNA_logfile.txt"
   if(!file.exists(logfile) == TRUE) file.create(paste0(dir, logfile))
  
   mapply(phased_piRNA_function, "+", vars[[1]], vars[[2]], vars[[3]], vars[[10]], logfile, dir, vars[[6]])
   mapply(phased_piRNA_function, "-", vars[[1]], vars[[2]], vars[[3]], vars[[10]], logfile, dir, vars[[6]])
   #mapply(phased_piRNA_function, "+", test_list$V1, test_list$V2, test_list$V3, input_file, logfile, dir, plot_output)
   #mapply(phased_piRNA_function, "-", test_list$V1, test_list$V2, test_list$V3, input_file, logfile, dir, plot_output)

}

