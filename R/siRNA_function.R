#' siRNA_function
#'
#' @param vars a list
#' @return si_res

#' @export



siRNA_function <- function(vars){
<<<<<<< Updated upstream
   dir <- "siRNA_outputs/"
=======
   `%>%` <- magrittr::`%>%`

   print("calling run_siRNA_function")
   wkdir <- "siRNA_outputs/"
>>>>>>> Stashed changes
   logfile <- "siRNA_logfile.txt"

   if(!dir.exists(wkdir) == TRUE) dir.create(wkdir)

   if(!file.exists(logfile) == TRUE) file.create(paste0(wkdir, logfile))


   si_res <- mapply(run_siRNA_function, vars[[1]], vars[[2]], vars[[3]], vars[[4]], 1, vars[[9]], vars[[10]], logfile, wkdir, vars[[13]], vars[[6]], vars[[7]],
                    vars[[14]], vars[[15]], vars[[16]])
   return(si_res)

}
