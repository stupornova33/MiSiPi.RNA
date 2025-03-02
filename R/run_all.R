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
  
  all_dir <- "run_all"
  all_log <- "run_all_logfile.txt"
  
  mi_dir <- file.path(all_dir, "miRNA_outputs")
  mi_log <- "miRNA_logfile.txt"
  
  pi_dir <- file.path(all_dir, "piRNA_outputs")
  pi_log <- "piRNA_logfile.txt"
  
  si_dir <- file.path(all_dir, "siRNA_outputs")
  si_log <- "siRNA_logfile.txt"
  
  if (dir.exists(all_dir)) {
    # If running interactively, give the option to delete old files if present,
    # or move them to a new timestamped directory if not
    if (!interactive()) {
      unlink(all_dir, recursive = TRUE)
    } else {
      .overwrite_warning(all_dir)
    }
  }
  
  dir.create(all_dir)
  dir.create(mi_dir)
  dir.create(pi_dir)
  dir.create(si_dir)
  
  file.create(file.path(all_dir, all_log))
  file.create(file.path(mi_dir, mi_log))
  file.create(file.path(pi_dir, pi_log))
  file.create(file.path(si_dir, si_log))
  
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
  
  .inform_complete(all_dir)
}