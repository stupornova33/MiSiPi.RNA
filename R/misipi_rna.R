#' run the run_all function
#' @param vars a list
#' @return plots

#' @export

misipi_rna <- function(vars){


mapply(new_run_all, vars$chrom_name, vars$reg_start, vars$reg_stop, vars$chromosome, vars$length, vars$bam_file,
       vars$roi,vars$genome, 1, vars$si_pal, vars$pi_pal, vars$plot_output, vars$path_to_RNAfold, vars$annotate_region,
       vars$weight_reads, vars$gtf_file, vars$write_fastas, vars$out_type)

}
