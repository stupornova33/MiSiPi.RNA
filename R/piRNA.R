# runs .piRNA and creates dir and logfile for outputs
#
# @param vars a list

piRNA <- function(vars) {
  total_iterations <- length(vars$chrom_name)
  idx_vec <- 1:total_iterations
  
  .print_intro(
    roi = vars$roi,
    bam = vars$bam_file,
    genome = vars$genome,
    method = "piRNA"
  )

  pi_dir <- "piRNA_outputs/"
  logfile <- file.path(pi_dir, "piRNA_logfile.txt")
  
  if (dir.exists(pi_dir)) {
    # If running interactively, give the option to delete old files if present,
    # or move them to a new timestamped directory if not
    if (!interactive()) {
      unlink(pi_dir, recursive = TRUE)
    } else {
      .overwrite_warning(pi_dir)
    }
  }
  
  dir.create(pi_dir)
  file.create(logfile)
  
  invisible(
    mapply(
      .piRNA,
      vars$chrom_name,
      vars$reg_start,
      vars$reg_stop,
      vars$length,
      vars$bam_file,
      vars$genome,
      logfile,
      pi_dir,
      vars$pi_pal,
      vars$plot_output,
      vars$weight_reads,
      vars$write_fastas,
      vars$out_type,
      idx_vec,
      total_iterations
    )
  )
  
  .inform_complete(pi_dir)
}
