# miRNA function
# calls the .miRNA function for both strands, creates miRNA logfile
# @param vars a list of variables provided by user
# @param output_dir the current iteration's output directory

miRNA <- function(vars, output_dir) {
  total_iterations <- length(vars$chrom_name)
  idx_vec <- 1:total_iterations

  .print_intro(
    roi = vars$roi,
    bam = vars$bam_file,
    genome = vars$genome,
    method = "miRNA"
  )

  mi_dir <- file.path(output_dir, "miRNA")
  logfile <- file.path(mi_dir, "miRNA_logfile.txt")
  
  dir.create(mi_dir)
  file.create(logfile)
  
  cli::cli_inform("Starting plus strand")
  # Process the positive strand
  invisible(
    mapply(
      .miRNA,
      vars$chrom_name,
      vars$reg_start,
      vars$reg_stop,
      vars$chromosome,
      vars$length,
      "+",
      vars$genome,
      vars$bam_file,
      logfile,
      mi_dir,
      vars$plot_output,
      vars$path_to_RNAfold,
      vars$path_to_RNAplot,
      vars$write_fastas,
      vars$weight_reads,
      vars$out_type,
      idx_vec,
      total_iterations
    )
  )

  cli::cli_inform(c("", "Starting minus strand"))
  # Process the minus strand
  invisible(
    mapply(
      .miRNA,
      vars$chrom_name,
      vars$reg_start,
      vars$reg_stop,
      vars$chromosome,
      vars$length,
      "-",
      vars$genome,
      vars$bam_file,
      logfile,
      mi_dir,
      vars$plot_output,
      vars$path_to_RNAfold,
      vars$path_to_RNAplot,
      vars$write_fastas,
      vars$weight_reads,
      vars$out_type,
      idx_vec,
      total_iterations
    )
  )
  
  .compare_miRNA_strands(vars$chrom_name, vars$reg_start, vars$reg_stop, output_dir)
  
  .inform_complete(mi_dir)
}
