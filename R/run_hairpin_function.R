#' run the hairpin function
#' runs hairpin function on both strands
#' takes a gtf file, bam file, and genome file specified by user
#' takes min_read_count specified by user. Default = 1
#'
#' @param vars a list

#' @return plots
#' @export

run_hairpin_function <- function(vars){
   wkdir <- 'hairpin_outputs/'
   if(!dir.exists(wkdir) == TRUE) dir.create(wkdir)


   logfile = "hairpin_logfile.txt"
   if(!file.exists(logfile) == TRUE) file.create(paste0(wkdir, logfile))

   mapply(dual_strand_hairpin, vars[[1]], vars[[2]], vars[[3]], vars[[4]], 1, vars[[9]], vars[[10]], logfile, wkdir, vars[[6]], vars[[7]], vars[[14]], vars[[15]])

}
