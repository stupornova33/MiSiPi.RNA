.plot_siRNA_results <- function(dsh, is_small_locus, annotate_region, results_present, dicer_plot, size_plot, heat_plot, out_type, prefix, wkdir) {
  ### combine siRNA and hpRNA plots
  
  
  if (is_small_locus) {
    #### Regions Less Than or Equal to 10kb (RLT10K) ####
    
    plus_hp_overhangs <- dsh$plus_overhang_plot
    minus_hp_overhangs <- dsh$minus_overhang_plot
    
    density_plot <- dsh$density_plot
    arc_plot <- dsh$arc_plot
    
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
