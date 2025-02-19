# Wrapper function that calls .run_all for each region of interest

run_all <- function(vars) {
  total_iterations <- length(vars$chrom_name)
  idx_vec <- 1:total_iterations
  
  .print_intro(
    roi = vars$roi,
    bam = vars$bam_file,
    genome = vars$genome,
    method = "all"
  )
  
  invisible(
    mapply(
      .run_all,
      vars$chrom_name,
      vars$reg_start,
      vars$reg_stop,
      vars$chromosome,
      vars$length,
      vars$bam_file,
      vars$roi,
      vars$genome,
      vars$si_pal,
      vars$pi_pal,
      vars$plot_output,
      vars$path_to_RNAfold,
      vars$path_to_RNAplot,
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