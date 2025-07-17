.plot_miRNA <- function(chrom_name, reg_start, reg_stop, strand, bam_file, fold_list, overhangs, stranded_size_dist, z_df, out_type, prefix, wkdir, method = c("self", "all")) {
  
  method = match.arg(method)
  
  dicer_sig <- .plot_overhangz(overhangs, strand = strand)
  
  # make a table with the same positions but empty cov column for combining with pileups
  # necessary because the pileups table doesn't have all positions from the first nt to the last because
  # coverages of zero aren't reported
  # set these to zero below
  empty_table <- data.frame(pos = c(seq(fold_list$start, fold_list$stop)), count = c(0))
  
  density <- .read_densityBySize(chrom_name, reg_start, reg_stop, bam_file, wkdir)
  
  density_plot <- .plot_density(density, reg_start, reg_stop)
  
  dist_plot <- .plot_sizes_by_strand(wkdir, stranded_size_dist, chrom_name, reg_start, reg_stop)
  zplot <- .plot_overlapz(z_df)
  
  # If .miRNA was called by run_all, then don't write the plots yet
  # Return each individual plot as an object to combine in .run_all
  if (method == "all") {
    plots <- list()
    plots$prefix <- prefix
    plots$strand <- strand
    plots$density <- density_plot
    plots$distribution <- dist_plot
    plots$dicer_sig <- dicer_sig
    plots$zplot <- zplot
    return(plots)
  }
  
  
  left_top <- cowplot::plot_grid(
    dist_plot,
    dicer_sig,
    ncol = 1,
    rel_widths = c(1, 1),
    rel_heights = c(1, 1),
    align = "vh",
    axis = "lrtb"
  )
  
  right_top <- cowplot::plot_grid(
    NULL,
    density_plot,
    zplot,
    ncol = 1,
    rel_widths = c(1, 1, 1),
    rel_heights = c(0.4, 1, 1),
    align = "vh",
    axis = "lrtb"
  )
  
  all_plot <- cowplot::plot_grid(
    left_top,
    right_top,
    rel_heights = c(1, 1),
    rel_widths = c(1, 1),
    align = "vh",
    axis = "lrtb"
  )
  
  if (out_type == "png") {
    grDevices::png(file = file.path(wkdir, paste(prefix, strand, "combined.png", sep = "_")), height = 8, width = 11, units = "in", res = 300)
  } else {
    grDevices::pdf(file = file.path(wkdir, paste(prefix, strand, "combined.pdf", sep = "_")), height = 8, width = 11)
  }
  print(all_plot)
  grDevices::dev.off()
  return(NULL)
}



.plot_siRNA <- function(dsh, is_small_locus, annotate_region, results_present, dicer_plot, size_plot, heat_plot, out_type, prefix, wkdir) {
  ### combine siRNA and hpRNA plots
  
  method = match.arg(method)
  
  if (is_small_locus) {
    #### Regions Less Than or Equal to 10kb (RLT10K) ####

    plus_hp_overhangs <- dsh$plus_overhang_plot
    minus_hp_overhangs <- dsh$minus_overhang_plot
    density_plot <- dsh$density_plot
    arc_plot <- dsh$arc_plot
    # In the case that no arc plot get made in dual_strand_hairpins
    # In order to keep track of that fact, we set the arc_plot to NA
    # Setting it to NULL at the time wouldn't have allowed us to store that
    # In a list with the other results. But we do need it to be NULL to keep
    # cowplot from complaining, so check for NA here and set to NULL
    # if it is not NA, then we can't really check it for NA, since the length
    # of the object will be greater than 1, so check for length should be
    # sufficient
    if (length(arc_plot) == 1) {
      arc_plot <- NULL
    }
    
    plus_phasedz <- dsh$plus_phasedz
    minus_phasedz <- dsh$minus_phasedz
    
    #### RLT10K - Annotate True ####
    if (annotate_region) {
      gtf_plot <- dsh$gtf_plot
      
      # if there are results for the heatmap, plot, otherwise omit
      if (results_present) {
        left <- cowplot::plot_grid(
          arc_plot,
          gtf_plot,
          density_plot,
          size_plot,
          ggplotify::as.grob(heat_plot),
          rel_widths = c(0.6, 1.1, 0.9, 0.9, 0.4),
          rel_heights = c(0.7, 0.7, 0.7, 0.7, 1.4),
          ncol = 1,
          align = "vh",
          axis = "lrtb")
      } else {
        left <- cowplot::plot_grid(
          arc_plot,
          gtf_plot,
          density_plot,
          size_plot,
          rel_widths = c(0.6, 1.1, 0.9, 0.9),
          rel_heights = c(0.7, 0.7, 0.7, 0.7, 1.4),
          ncol = 1,
          align = "vh",
          axis = "lrtb")
      }
      
      #### RLT10K - Annotate False ####
    } else {
      if (results_present) {
        left <- cowplot::plot_grid(
          arc_plot,
          density_plot,
          size_plot,
          ggplotify::as.grob(heat_plot),
          ncol = 1,
          rel_widths = c(0.6, 0.9, 0.9, 0.4),
          rel_heights = c(0.7, 0.7, 0.7, 1.4),
          align = "vh",
          axis = "lrtb")
      } else {
        left <- cowplot::plot_grid(
          arc_plot,
          density_plot,
          size_plot,
          rel_widths = c(0.6, 0.9, 0.9),
          rel_heights = c(0.7, 0.7, 0.7, 1.4),
          ncol = 1,
          align = "vh",
          axis = "lrtb")
      }
    }
    
    right <- cowplot::plot_grid(
      plus_hp_overhangs,
      minus_hp_overhangs,
      plus_phasedz,
      minus_phasedz,
      dicer_plot,
      ncol = 1,
      align = "vh",
      axis = "l",
      rel_widths = c(1, 1, 1, 1, 1),
      rel_heights = c(1, 1, 1, 1, 1))
    
    all_plot <- cowplot::plot_grid(
      left,
      NULL,
      right,
      ncol = 3,
      rel_widths = c(0.9, 0.01, 0.7),
      align = "vh",
      axis = "lrtb")
    
    if (out_type == "png") {
      grDevices::png(file = file.path(wkdir, paste0(prefix, "_si_plot.png")), height = 16, width = 14, units = "in", res = 300)
    } else {
      grDevices::pdf(file = file.path(wkdir, paste0(prefix, "_si_plot.pdf")), height = 16, width = 14)
    }
    
    print(all_plot)
    grDevices::dev.off()
    
    
  } else {
    #### Regions Greater Than 10kb (RGT10k) ####
    plus_hp_overhangs <- dsh$plus_overhang_plot
    minus_hp_overhangs <- dsh$minus_overhang_plot

    # None of the hairpin plots were made because the region > 10kb
    density_plot <- dsh$density_plot
    
    #### RGT10k - Annotate True ####
    if (annotate_region) {
      gtf_plot <- dsh$gtf_plot
      
      if (results_present) {
        bottom <- cowplot::plot_grid(
          ggplotify::as.grob(heat_plot),
          dicer_plot,
          ncol = 2,
          rel_widths = c(1, 1),
          rel_heights = c(1, 1),
          align = "vh",
          axis = "lrtb")
      } else { # if no heat plot
        bottom <- cowplot::plot_grid(
          dicer_plot,
          rel_widths = c(1),
          rel_heights = c(1),
          align = "vh",
          axis = "lrtb")
      }
      
      top <- cowplot::plot_grid(
        density_plot,
        size_plot,
        ncol = 2,
        rel_widths = c(1, 1),
        rel_heights = c(1, 1),
        align = "vh",
        axis = "lrtb")
      
      all_plot <- cowplot::plot_grid(
        top,
        NULL,
        bottom,
        ncol = 1,
        rel_widths = c(1, 1, 1),
        rel_heights = c(1, 0.1, 1),
        align = "vh",
        axis = "lrtb")
      
      if (out_type == "png") {
        grDevices::png(file = file.path(wkdir, paste0(prefix, "_si_plot.png")), height = 9, width = 14, units = "in", res = 300)
      } else {
        grDevices::pdf(file = file.path(wkdir, paste0(prefix, "_si_plot.pdf")), height = 9, width = 14)
      }
      
      print(all_plot)
      grDevices::dev.off()
      
      #### RGT10k - Annotate False ####
    } else {
      if (results_present) {
        left <- cowplot::plot_grid(
          density_plot,
          size_plot,
          dicer_plot,
          ncol = 1,
          rel_widths = c(1, 1),
          rel_heights = c(1, 1),
          align = "vh",
          axis = "lrtb")
        
        right <- cowplot::plot_grid(
          NULL,
          ggplotify::as.grob(heat_plot),
          NULL,
          ncol = 1,
          align = "vh",
          axis = "l",
          rel_widths = c(0.4, 1, 0.4),
          rel_heights = c(1, 1))
        
        all_plot <- cowplot::plot_grid(
          left,
          NULL,
          right,
          ncol = 3,
          rel_widths = c(1, 0.1, 0.8),
          rel_heights = c(1, 1, 1),
          align = "vh",
          axis = "lrtb")
        
      } else {
        bottom <- cowplot::plot_grid(
          dicer_plot,
          size_plot,
          rel_widths = c(1, 1),
          rel_heights = c(1, 1),
          ncol = 2,
          align = "vh",
          axis = "lrtb")
        
        top <- cowplot::plot_grid(
          density_plot,
          ncol = 1,
          align = "vh",
          axis = "l",
          rel_widths = c(1),
          rel_heights = c(1))
        
        all_plot <- cowplot::plot_grid(
          top,
          NULL,
          bottom,
          ncol = 1,
          rel_widths = c(1, 1),
          rel_heights = c(1, 0.1, 1),
          align = "vh",
          axis = "lrtb")
      }
    }
    
    if (out_type == "png") {
      grDevices::png(file = file.path(wkdir, paste0(prefix, "_si_plot.png")), height = 13, width = 13, units = "in", res = 300)
    } else {
      grDevices::pdf(file = file.path(wkdir, paste0(prefix, "_si_plot.pdf")), height = 13, width = 13)
    }
    
    print(all_plot)
    grDevices::dev.off()
    
  }
}
