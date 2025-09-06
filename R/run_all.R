# Wrapper function that calls .run_all for each region of interest

run_all <- function(vars, output_dir) {
  total_iterations <- length(vars$chrom_name)
  idx_vec <- 1:total_iterations
  
  .print_intro(
    roi = vars$roi,
    bam = vars$bam_file,
    genome = vars$genome,
    method = "all"
  )
  
  all_log <- "run_all_log.txt"
  
  mi_dir <- file.path(output_dir, "miRNA")
  mi_log <- "miRNA_log.txt"
  
  pi_dir <- file.path(output_dir, "piRNA")
  pi_log <- "piRNA_log.txt"
  
  si_dir <- file.path(output_dir, "siRNA")
  si_log <- "siRNA_log.txt"
  
  plot_dir <- file.path(output_dir, "combined_plots")
  
  
  dir.create(mi_dir)
  dir.create(pi_dir)
  dir.create(si_dir)
  dir.create(plot_dir)
  
  file.create(file.path(output_dir, all_log))
  file.create(file.path(mi_dir, mi_log))
  file.create(file.path(pi_dir, pi_log))
  file.create(file.path(si_dir, si_log))
  
  invisible(
    mapply(
      .run_all,
      vars$chrom_name,
      vars$reg_start,
      vars$reg_stop,
      vars$prefix,
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
      output_dir,
      vars$use_bed_names,
      idx_vec,
      total_iterations,
      vars$iteration_output
    )
  )
  
  .inform_complete(output_dir)
}