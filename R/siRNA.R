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
  
  si_dir <- "siRNA_outputs"
  logfile <- file.path(si_dir, "siRNA_logfile.txt")
  
  if (dir.exists(si_dir)) {
    # If running interactively, give the option to delete old files if present,
    # or move them to a new timestamped directory if not
    if (!interactive()) {
      unlink(si_dir, recursive = TRUE)
    } else {
      .overwrite_warning(si_dir)
    }
  }
  
  dir.create(si_dir)
  file.create(logfile)

  invisible(
    mapply(
      .siRNA,
      vars$chrom_name,
      vars$reg_start,
      vars$reg_stop,
      vars$length,
      vars$genome,
      vars$bam_file,
      logfile,
      si_dir,
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
  
  .inform_complete(si_dir)
}
