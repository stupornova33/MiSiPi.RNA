#' run_phased_piRNA
#' calls phased_piRNA function on both strand, creates output wkdir and logfile
#' @param vars a list
#' @return plots

#' @export

run_phased_piRNA <- function(vars){


   logfile = "phased_piRNA_logfile.txt"

   wkdir <- 'phased_piRNA_outputs/'
   if(!dir.exists(wkdir) == TRUE) dir.create(wkdir)

   logfile = "phased_piRNA_logfile.txt"
   if(!file.exists(logfile) == TRUE) file.create(paste0(wkdir, logfile))

   mapply(phased_piRNA_function, "+", vars[[1]], vars[[2]], vars[[3]], vars[[10]], logfile, wkdir, vars[[6]])
   mapply(phased_piRNA_function, "-", vars[[1]], vars[[2]], vars[[3]], vars[[10]], logfile, wkdir, vars[[6]])

}

