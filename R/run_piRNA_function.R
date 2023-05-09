#' runs piRNA_function and creates dir and logfile for outputs
#' 
#' @param vars
#' @return pi_res, plots
#' @export


run_piRNA_function <- function(){
   
   dir <- 'piRNA_outputs/'
   if(!dir.exists(dir) == TRUE) dir.create(dir)
   
   logfile = "piRNA_logfile.txt"
   if(!file.exists(logfile) == TRUE) file.create(paste0(dir, logfile))
   

   
   mapply(piRNA_function, vars[[1]], vars[[2]], vars[[3]], vars[[10]], logfile, dir, vars[[12]], vars[[6]])  


}