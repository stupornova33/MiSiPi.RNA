# miRNA function
# calls the .miRNA function for both strands, creates miRNA logfile
# @param vars a list of variables provided by user

miRNA <- function(vars) {
  total_iterations <- length(vars$chrom_name)
  idx_vec <- 1:total_iterations

  .print_intro(
    roi = vars$roi,
    bam = vars$bam_file,
    genome = vars$genome,
    method = "miRNA"
  )

  logfile <- "miRNA_logfile.txt"

  wkdir <- "miRNA_outputs/"
  if (!dir.exists(wkdir) == TRUE) dir.create(wkdir)

  if (!file.exists(logfile) == TRUE) file.create(paste0(wkdir, logfile))

  print("Starting plus strand")
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
      1,
      vars$genome,
      vars$bam_file,
      logfile,
      wkdir,
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
  
  print("Starting minus strand")
  # Process the negative strand
  invisible(
    mapply(
      .miRNA,
      vars$chrom_name,
      vars$reg_start,
      vars$reg_stop,
      vars$chromosome,
      vars$length,
      "-",
      1,
      vars$genome,
      vars$bam_file,
      logfile,
      wkdir,
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
}
