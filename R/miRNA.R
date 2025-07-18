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
  
  genome <- vars$genome
  bam <- vars$bam_file
  plot_output <- vars$plot_output
  RNAfold_path <- vars$path_to_RNAfold
  RNAplot_path <- vars$path_to_RNAplot
  write_fastas <- vars$write_fastas
  weight_reads <- vars$weight_reads
  out_type <- vars$out_type
  
  for (i in idx_vec) {
    chrom <- vars$chrom_name[i]
    reg_start <- vars$reg_start[i]
    reg_stop <- vars$reg_stop[i]
    chr_length <- vars$length[i]
    
    cli::cli_inform("Starting plus strand")
    plus_results <- .miRNA(
      chrom_name = chrom,
      reg_start = reg_start,
      reg_stop = reg_stop,
      length = chr_length,
      strand = "+",
      genome_file = genome,
      bam_file = bam,
      logfile = logfile,
      wkdir = mi_dir,
      plot_output = plot_output,
      path_to_RNAfold = RNAfold_path,
      path_to_RNAplot = RNAplot_path,
      write_fastas = write_fastas,
      weight_reads = weight_reads,
      out_type = out_type,
      i = i,
      i_total = total_iterations
    )
    
    cli::cli_inform(c("", "Starting minus strand"))
    minus_results <- .miRNA(
      chrom_name = chrom,
      reg_start = reg_start,
      reg_stop = reg_stop,
      length = chr_length,
      strand = "-",
      genome_file = genome,
      bam_file = bam,
      logfile = logfile,
      wkdir = mi_dir,
      plot_output = plot_output,
      path_to_RNAfold = RNAfold_path,
      path_to_RNAplot = RNAplot_path,
      write_fastas = write_fastas,
      weight_reads = weight_reads,
      out_type = out_type,
      i = i,
      i_total = total_iterations
    )
    
    .compare_miRNA_strands(vars$chrom_name, vars$reg_start, vars$reg_stop, output_dir)
    
    # Generate plots
    plus_plots <- plus_results$plots
    plus_overhangs <- plus_results$overhangs
    plus_overlaps <- plus_results$overlaps
    
    minus_plots <- minus_results$plots
    minus_overhangs <- minus_results$overhangs
    minus_overlaps <- minus_results$overlaps
    
    dicer_overhang_plot <- .plot_miRNA_dicer_overhang_probability(plus_overhangs, minus_overhangs)
    overlap_probability_plot <- .plot_miRNA_overlap_probability(plus_overlaps, minus_overlaps)
    
    
    if (!is.null(plus_plots$distribution)) {
      read_distribution_plot <- plus_plots$distribution
    } else {
      read_distribution_plot <- minus_plots$distribution
    }
    
    if (!is.null(plus_plots$density)) {
      read_density_plot <- plus_plots$density
    } else {
      read_density_plot <- minus_plots$density
    }
    
    prefix <- .get_region_string(chrom, reg_start, reg_stop)
    
    # Print the plots to a file
    .print_miRNA_plots(read_distribution_plot, read_density_plot, dicer_overhang_plot, overlap_probability_plot, out_type, prefix, mi_dir)

  }
  
  .inform_complete(mi_dir)
}
