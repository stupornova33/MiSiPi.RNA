# run .siRNA function
# calls the siRNA function, creates siRNA logfile
# @param vars a list of variables provided by user

siRNA <- function(vars) {
  total_iterations <- length(vars$chrom_name)
  idx_vec <- 1:total_iterations
  
  .print_intro(
    roi = vars$roi,
    bam = vars$bam_file,
    genome = vars$genome,
    method = "siRNA"
  )
  
  dir <- "siRNA_outputs/"

  wkdir <- "siRNA_outputs/"
  logfile <- "siRNA_logfile.txt"

  if (!dir.exists(wkdir) == TRUE) dir.create(wkdir)

  if (!file.exists(logfile) == TRUE) file.create(paste0(wkdir, logfile))

  invisible(
    mapply(
      .siRNA,
      vars$chrom_name,
      vars$reg_start,
      vars$reg_stop,
      vars$length,
      1,
      vars$genome,
      vars$bam_file,
      logfile,
      wkdir,
      vars$si_pal,
      vars$plot_output,
      vars$path_to_RNAfold,
      vars$annotate_region,
      vars$weight_reads,
      vars$gtf_file,
      vars$write_fastas,
      vars$out_type,
      idx_vec,
      total_iterations
    )
  )
}
