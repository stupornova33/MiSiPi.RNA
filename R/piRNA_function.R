#' runs piRNA_function and creates dir and logfile for outputs
#'
#' @param vars a list
#' @return pi_res, plots
#' @export

piRNA_function <- function(vars) {
  wkdir <- "piRNA_outputs/"
  if (!dir.exists(wkdir) == TRUE) dir.create(wkdir)

  logfile <- "piRNA_logfile.txt"
  if (!file.exists(logfile) == TRUE) file.create(paste0(wkdir, logfile))

  mapply(
    run_piRNA_function, vars$chrom_name, vars$reg_start, vars$reg_stop, vars$length, vars$bam_file, vars$genome, logfile, wkdir, vars$pi_pal,
    vars$plot_output, vars$weight_reads, vars$write_fastas, vars$out_type
  )
}
