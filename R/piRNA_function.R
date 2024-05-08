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

                            #chrom_name  reg_start  reg_stop   length     bam_file    genome     logfile  wikdir  pal
   mapply(run_piRNA_function, vars[[1]], vars[[2]], vars[[3]], vars[[4]], vars[[10]], vars[[9]], logfile, wkdir, vars[[12]],
            #plot_output weight_reads write_fastas out_type
           vars[[6]], vars[[15]], vars[[17]], vars[[18]])

}
