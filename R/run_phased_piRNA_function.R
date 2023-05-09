#' run_phased_piRNA
#' calls phased_piRNA function on both strand, creates output dir and logfile
#' @param vars a list
#' @return plots

#' @export

run_phased_piRNA <- function(vars){

   
   logfile = "phased_piRNA_logfile.txt"
   
   dir <- 'phased_piRNA_outputs/'
   if(!dir.exists(dir) == TRUE) dir.create(dir)
   
   logfile = "phased_piRNA_logfile.txt"
   if(!file.exists(logfile) == TRUE) file.create(paste0(dir, logfile))
  
   mapply(phased_piRNA_function, "+", vars[[1]], vars[[2]], vars[[3]], vars[[10]], logfile, dir, vars[[6]])
   mapply(phased_piRNA_function, "-", vars[[1]], vars[[2]], vars[[3]], vars[[10]], logfile, dir, vars[[6]])

}

