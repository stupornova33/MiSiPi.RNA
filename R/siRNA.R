#' run .siRNA function
#' calls the siRNA function, creates siRNA logfile
#' @param vars a list of variables provided by user
#' @importFrom Rcpp sourceCpp
#' @return si_res

#' @export

siRNA <- function(vars) {
  dir <- "siRNA_outputs/"

  wkdir <- "siRNA_outputs/"
  logfile <- "siRNA_logfile.txt"

  if (!dir.exists(wkdir) == TRUE) dir.create(wkdir)

  if (!file.exists(logfile) == TRUE) file.create(paste0(wkdir, logfile))

  si_res <- mapply(
    .siRNA, vars$chrom_name, vars$reg_start, vars$reg_stop, vars$length, 1, vars$genome, vars$bam_file, logfile,
    wkdir, vars$si_pal, vars$plot_output, vars$path_to_RNAfold, vars$annotate_region, vars$weight_reads, vars$gtf_file,
    vars$write_fastas, vars$out_type
  )
  return(si_res)
}
