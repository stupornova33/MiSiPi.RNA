#' run_siRNA_function
#'
#' @param vars a list
#' @return si_res

#' @export



run_siRNA_function <- function(vars){
   `%>%` <- magrittr::`%>%`
   dir <- "siRNA_outputs/"
   logfile <- "siRNA_logfile.txt"
   
   if(!dir.exists(dir) == TRUE) dir.create(dir)
   
   if(!file.exists(logfile) == TRUE) file.create(paste0(dir, logfile))

 
   si_res <- mapply(siRNA_function, vars[[1]], vars[[2]], vars[[3]], vars[[4]], 1, vars[[9]], vars[[10]], logfile, dir, vars[[13]], vars[[6]], vars[[7]], 
                    vars[[14]], vars[[15]])
   return(si_res) 
   
}