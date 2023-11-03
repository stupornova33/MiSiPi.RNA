#' runs piRNA_function and creates dir and logfile for outputs
#'
#' @param vars a list
#' @return pi_res, plots
#' @export


piRNA_function <- function(vars){

   wkdir <- 'piRNA_outputs/'
   if(!dir.exists(wkdir) == TRUE) dir.create(wkdir)

   logfile = "piRNA_logfile.txt"
   if(!file.exists(logfile) == TRUE) file.create(paste0(wkdir, logfile))



   mapply(run_piRNA_function, vars[[1]], vars[[2]], vars[[3]], vars[[10]], logfile, wkdir, vars[[12]], vars[[6]])


}
