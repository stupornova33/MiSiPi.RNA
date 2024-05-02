#' run siRNA function
#' calls the siRNA function, creates siRNA logfile
#' @param vars a list of variables provided by user
#' @importFrom Rcpp sourceCpp
#' @return si_res

#' @export

siRNA_function <- function(vars){
   dir <- "siRNA_outputs/"

   `%>%` <- magrittr::`%>%`

   print("calling run_siRNA_function")
   wkdir <- "siRNA_outputs/"
   logfile <- "siRNA_logfile.txt"

   if(!dir.exists(wkdir) == TRUE) dir.create(wkdir)

   if(!file.exists(logfile) == TRUE) file.create(paste0(wkdir, logfile))
   print("Logfile created. Now running siRNA_function.")

   #si_res <- mapply(run_siRNA_function, vars[[1]], vars[[2]], vars[[3]], vars[[4]], 1, vars[[9]], vars[[10]], logfile, wkdir, vars[[13]], vars[[6]], vars[[7]],
  #                  vars[[14]], vars[[15]], vars[[16]], vars[[17]], vars[[18]])


   si_res <- mapply(run_siRNA_function, vars$chrom_name, vars$reg_start, vars$reg_stop, vars$length, 1, vars$genome, vars$bam_file, logfile,
                    wkdir, vars$si_pal, vars$plot_output, vars$path_to_RNAfold, vars$annotate_region, vars$weight_reads, vars$gtf_file,
                    vars$write_fastas, vars$out_type)
   return(si_res)

}
