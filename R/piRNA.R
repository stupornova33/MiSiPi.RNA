# runs .piRNA and creates dir and logfile for outputs
#
# @param vars a list

piRNA <- function(vars, output_dir) {
  total_iterations <- length(vars$chrom_name)
  idx_vec <- 1:total_iterations
  
  .print_intro(
    roi = vars$roi,
    bam = vars$bam_file,
    genome = vars$genome,
    method = "piRNA"
  )

  pi_dir <- file.path(output_dir, "piRNA")
  logfile <- file.path(pi_dir, "piRNA_log.txt")
  
  dir.create(pi_dir)
  file.create(logfile)
  
  calling_method <- "self"
  
  invisible(
    mapply(
      .piRNA,
      vars$chrom_name,
      vars$reg_start,
      vars$reg_stop,
      vars$prefix,
      vars$bam_file,
      vars$genome,
      vars$roi,
      logfile,
      pi_dir,
      vars$pi_pal,
      vars$plot_output,
      vars$weight_reads,
      vars$write_fastas,
      vars$out_type,
      calling_method,
      idx_vec,
      total_iterations,
      vars$iteration_input
    )
  )
  
  .inform_complete(pi_dir)
}
