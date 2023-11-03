#' run the hairpin function
#' runs hairpin function on both strands
#' takes a bed file, bam file, and genome file specified by user
#' takes min_read_count specified by user. Default = 1
#' 
#' @param vars a list

#' @return plots
#' @export

run_hairpin_function <- function(vars){
   dir <- 'hairpin_outputs/'
   if(!dir.exists(dir) == TRUE) dir.create(dir)

   
   logfile = "hairpin_logfile.txt"
   if(!file.exists(logfile) == TRUE) file.create(paste0(dir, logfile))
  
   mapply(dual_strand_hairpin, vars[[1]], vars[[2]], vars[[3]], vars[[4]], 1, vars[[9]], vars[[10]], logfile, dir, vars[[6]], vars[[7]], vars[[14]], vars[[15]])
   #mapply(hairpin_function, vars[[1]], vars[[2]], vars[[3]], "-", vars[[4]], 1, vars[[9]], vars[[10]], logfile, dir, vars[[6]], vars[[7]], vars[[14]], vars[[15]])
   #mapply(hairpin_function, test_list$V1, test_list$V2, test_list$V3, "+", test_list$length, min_read_count, genome_file, input_file, logfile, dir, plot_output, path_to_RNAfold)
   #mapply(hairpin_function, test_list$V1, test_list$V2, test_list$V3, "-", test_list$length, min_read_count, genome_file, input_file, logfile, dir, plot_output, path_to_RNAfold)
}
