#' run miRNA function
#' calls the miRNA function for both strands, creates miRNA logfile
#' @param vars a list of variables provided by user
#' @importFrom Rcpp sourceCpp
#' @return plots

#' @export

miRNA_function <- function(vars){

   logfile <- "miRNA_logfile.txt"

   wkdir <- 'miRNA_outputs/'
   if(!dir.exists(wkdir) == TRUE) dir.create(wkdir)

   if(!file.exists(logfile) == TRUE) file.create(paste0(wkdir, logfile))

   mapply(run_miRNA_function, vars[[1]], vars[[2]], vars[[3]], vars[[5]], vars[[4]], "+", 1, vars[[9]], vars[[10]], logfile,
          wkdir, vars[[6]], vars[[7]], vars[[15]], vars[[17]])
   mapply(run_miRNA_function, vars[[1]], vars[[2]], vars[[3]], vars[[5]], vars[[4]], "-", 1, vars[[9]], vars[[10]], logfile,
          wkdir, vars[[6]], vars[[7]], vars[[15]], vars[[17]])


}
