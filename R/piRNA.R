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
  
  wkdir <- "piRNA_outputs/"
  if (!dir.exists(wkdir) == TRUE) dir.create(wkdir)

  logfile <- "piRNA_logfile.txt"
  if (!file.exists(logfile) == TRUE) file.create(paste0(wkdir, logfile))

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
      wkdir,
      vars$pi_pal,
      vars$plot_output,
      vars$weight_reads,
      vars$write_fastas,
      vars$out_type,
      idx_vec,
      total_iterations
    )
  )
}
