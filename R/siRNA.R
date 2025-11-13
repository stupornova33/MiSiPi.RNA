# run .siRNA function
# calls the siRNA function, creates siRNA logfile
# @param vars a list of variables provided by user

siRNA <- function(vars, output_dir) {
  total_iterations <- length(vars$chrom_name)
  idx_vec <- 1:total_iterations
  
  .print_intro(
    roi = vars$roi,
    bam = vars$bam_file,
    genome = vars$genome,
    method = "siRNA"
  )
  
  si_dir <- file.path(output_dir, "siRNA")
  logfile <- file.path(si_dir, "siRNA_log.txt")
  
  dir.create(si_dir)
  file.create(logfile)
  
  calling_method <- "self"

  invisible(
    mapply(
      .siRNA,
      vars$chrom_name,
      vars$reg_start,
      vars$reg_stop,
      vars$prefix,
      vars$genome,
      vars$bam_file,
      vars$roi,
      logfile,
      si_dir,
      vars$si_pal,
      vars$plot_output,
      vars$path_to_RNAfold,
      vars$annotate_region,
      vars$weight_reads,
      vars$write_fastas,
      vars$out_type,
      calling_method,
      idx_vec,
      total_iterations,
      vars$iteration_output,
      vars$density_timeout,
      MoreArgs = list(vars$gtf_df)
    )
  )
  
  .inform_complete(si_dir)
}
