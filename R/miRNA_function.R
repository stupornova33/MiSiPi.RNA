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


   #mapply(run_miRNA_function, vars$chrom_name, vars$reg_start, vars$reg_stop, vars$chromosome, vars$length, "+", 1, vars$genome,
   #      vars$bam_file, logfile, wkdir, vars$plot_output, vars$path_to_RNAfold,vars$write_fastas, vars$weight_reads, vars$out_type)

   #mapply(run_miRNA_function, vars$chrom_name, vars$reg_start, vars$reg_stop, vars$chromosome, vars$length, "-", 1, vars$genome,
   #      vars$bam_file, logfile, wkdir, vars$plot_output, vars$path_to_RNAfold, vars$write_fastas, vars$weight_reads, vars$out_type)

   mapply(new_miRNA_function, vars$chrom_name, vars$reg_start, vars$reg_stop, vars$chromosome, vars$length, "+", 1, vars$genome,
          vars$bam_file, logfile, wkdir, vars$plot_output, vars$path_to_RNAfold,
          vars$path_to_RNAplot, vars$write_fastas, vars$weight_reads, vars$out_type)

   mapply(new_miRNA_function, vars$chrom_name, vars$reg_start, vars$reg_stop, vars$chromosome, vars$length, "-", 1, vars$genome,
          vars$bam_file, logfile, wkdir, vars$plot_output, vars$path_to_RNAfold,
          vars$path_to_RNAplot, vars$write_fastas, vars$weight_reads, vars$out_type)


}
