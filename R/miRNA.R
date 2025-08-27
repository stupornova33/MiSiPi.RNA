# miRNA function
# calls the .miRNA function for both strands, creates miRNA logfile
# @param vars a list of variables provided by user
# @param output_dir the current iteration's output directory

miRNA <- function(vars, output_dir) {
  total_iterations <- length(vars$chrom_name)
  idx_vec <- 1:total_iterations

  genome <- vars$genome
  bam <- vars$bam_file
  bed <- vars$roi
  plot_output <- vars$plot_output
  RNAfold_path <- vars$path_to_RNAfold
  RNAplot_path <- vars$path_to_RNAplot
  write_fastas <- vars$write_fastas
  weight_reads <- vars$weight_reads
  out_type <- vars$out_type
  
  .print_intro(
    roi = bed,
    bam = bam,
    genome = genome,
    method = "miRNA"
  )

  mi_dir <- file.path(output_dir, "miRNA")
  logfile <- file.path(mi_dir, "miRNA_logfile.txt")
  
  dir.create(mi_dir)
  file.create(logfile)
  
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
    
    .compare_miRNA_strands(output_dir)
    
    if (plot_output == TRUE) {
      # Generate plots
      plus_overhangs <- plus_results$overhangs
      plus_overlaps <- plus_results$overlaps
      
      minus_overhangs <- minus_results$overhangs
      minus_overlaps <- minus_results$overlaps
      
      dicer_overhang_plot <- .plot_miRNA_dicer_overhang_probability(plus_overhangs, minus_overhangs)
      overlap_probability_plot <- .plot_miRNA_overlap_probability(plus_overlaps, minus_overlaps)
      
      # Moved the density and distribution plots here so they wouldn't be called twice if both strands get processed
      # If miRNA is called from run_all, then density and distribution will be generated in that function
      bam_obj <- .open_bam(bam_file, logfile)
      plus_df <- .get_filtered_bam_df(bam_obj, chrom_name, reg_start, reg_stop, "plus", 18, 32, FALSE)
      minus_df <- .get_filtered_bam_df(bam_obj, chrom_name, reg_start, reg_stop, "plus", 18, 32, FALSE)

      stranded_size_dist <- .get_stranded_read_dist(plus_df, minus_df)
      read_distribution_plot <- .plot_sizes_by_strand(stranded_size_dist)
      .close_bam(bam_obj)
      
      density_data <- .read_densityBySize(chrom, reg_start, reg_stop, bam, mi_dir)
      read_density_plot <- .plot_density(density_data, reg_start, reg_stop)
      
      prefix <- .get_region_string(chrom, reg_start, reg_stop)
      
      plot_details <- plot_title(bam, bed, genome, chrom, reg_start, reg_stop, i)
      
      # Print the plots to a file
      .print_miRNA_plots(read_distribution_plot, read_density_plot, dicer_overhang_plot, overlap_probability_plot, out_type, prefix, mi_dir, plot_details)
    }
  }
  
  .inform_complete(mi_dir)
}
