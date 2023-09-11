#' run miRNA function
#' calls the miRNA function for both strands, creates miRNA logfile
#' @param vars a list of variables provided by user
#' @importFrom Rcpp sourceCpp
#' @return plots

#' @export

run_miRNA_function <- function(vars){
   `%>%` <- magrittr::`%>%`

   logfile = "miRNA_logfile.txt"

   dir <- 'miRNA_outputs/'
   if(!dir.exists(dir) == TRUE) dir.create(dir)

   if(!file.exists(logfile) == TRUE) file.create(paste0(dir, logfile))

   mapply(miRNA_function, vars[[1]], vars[[2]], vars[[3]], vars[[5]], vars[[4]], "+", 1, vars[[9]], vars[[10]], logfile, dir, vars[[6]], vars[[7]], vars[[15]])
   mapply(miRNA_function, vars[[1]], vars[[2]], vars[[3]], vars[[5]], vars[[4]], "-", 1, vars[[9]], vars[[10]], logfile, dir, vars[[6]], vars[[7]], vars[[15]])


}
