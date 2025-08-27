# Write the plots to a file
.print_miRNA_plots <- function(read_distribution_plot, read_density_plot, dicer_overhang_plot, overlap_probability_plot, out_type, prefix, wkdir, plot_details) {
  
  plot_title <- cowplot::ggdraw() +
    cowplot::draw_label(
      plot_details$title,
      x = 0.5,
      hjust = 0.5,
      size = 12
    ) +
    ggplot2::theme(plot.margin = ggplot2::margin(7, 0, 0, 0))
  
  plot_caption <- cowplot::ggdraw() +
    cowplot::draw_label(
      plot_details$caption,
      x = 0.5,
      hjust = 0.5,
      size = 12
    ) +
    ggplot2::theme(plot.margin = ggplot2::margin(0, 0, 5, 0))
  
  plot_body <- cowplot::plot_grid(
    read_distribution_plot, read_density_plot,
    dicer_overhang_plot, overlap_probability_plot,
    ncol = 2,
    align = "hv",
    axis = "lrtb"
  )
  
  all_plot <- cowplot::plot_grid(
    plot_title,
    plot_body,
    plot_caption,
    ncol = 1,
    rel_heights = c(0.07, 1, 0.07)
  )
  
  if (out_type == "png") {
    grDevices::png(file = file.path(wkdir, paste(prefix, "combined.png", sep = "_")), height = 8, width = 11, units = "in", res = 300)
  } else {
    grDevices::pdf(file = file.path(wkdir, paste(prefix, "combined.pdf", sep = "_")), height = 8, width = 11)
  }
  print(all_plot)
  grDevices::dev.off()
  return(NULL)
}

.plot_siRNA <- function(read_distribution_plot, density_plot, phasedz_plot, overhang_probability_plot, arc_plot, gtf_plot, heat_plot, out_type, prefix, wkdir, plot_details) {
  plot_title <- cowplot::ggdraw() +
    cowplot::draw_label(
      plot_details$title,
      x = 0.5,
      hjust = 0.5,
      size = 12
    ) +
    ggplot2::theme(plot.margin = ggplot2::margin(7, 0, 0, 0))
  
  plot_caption <- cowplot::ggdraw() +
    cowplot::draw_label(
      plot_details$caption,
      x = 0.5,
      hjust = 0.5,
      size = 12
    ) +
    ggplot2::theme(plot.margin = ggplot2::margin(0, 0, 5, 0))
  
  plot_body <- cowplot::plot_grid(
    read_distribution_plot, arc_plot,
    overhang_probability_plot, density_plot,
    phasedz_plot, gtf_plot,
    heat_plot,
    ncol = 2,
    #align = "hv",
    axis = "lrtb",
    rel_widths = c(1, 1,0.8,1,0.8,1,0.8)
    
  )
  
  all_plot <- cowplot::plot_grid(
    plot_title,
    plot_body,
    plot_caption,
    ncol = 1,
    rel_heights = c(0.07, 1, 0.07)
  )
  
  # TODO: Change dimensions to better suit new layout
  if (out_type == "png") {
    grDevices::png(file = file.path(wkdir, paste0(prefix, "_si_plot.png")), height = 20, width = 14, units = "in", res = 300)
  } else {
    grDevices::pdf(file = file.path(wkdir, paste0(prefix, "_si_plot.pdf")), height = 20, width = 14)
  }
  
  print(all_plot)
  grDevices::dev.off()
  return(NULL)
}

plot_piRNA <- function(read_distribution_plot, density_plot, overlap_probability_plot, phased_probability_plot, heat_plot, out_type, prefix, wkdir, plot_details) {
  # Modify density_plot to display better in this layout
  density_plot <- density_plot +
    ggplot2::theme(axis.text.x = ggplot2::element_text(size = 12, angle = 60, hjust = 1),
                   axis.text.y = ggplot2::element_text(size = 12),
                   axis.title.x = ggplot2::element_text(size = 12),
                   axis.title.y = ggplot2::element_text(size = 12))
  
  
  
  plot_title <- cowplot::ggdraw() +
    cowplot::draw_label(
      plot_details$title,
      x = 0.5,
      hjust = 0.5,
      size = 12
    ) +
    ggplot2::theme(plot.margin = ggplot2::margin(7, 0, 0, 0))
  
  plot_caption <- cowplot::ggdraw() +
    cowplot::draw_label(
      plot_details$caption,
      x = 0.5,
      hjust = 0.5,
      size = 12
    ) +
    ggplot2::theme(plot.margin = ggplot2::margin(0, 0, 5, 0))
  
  plot_body_top <- cowplot::plot_grid(
    read_distribution_plot, heat_plot,
    overlap_probability_plot, phased_probability_plot,
    ncol = 2,
    align = "hv",
    axis = "lrtb"
  )
  
  plot_body_bottom <- cowplot::plot_grid(
    NULL, density_plot, NULL,
    ncol = 3,
    rel_widths = c(0.1, 1, 0.1),
    align = "hv",
    axis = "lrtb"
  )
  
  all_plot <- cowplot::plot_grid(
    plot_title,
    plot_body_top,
    plot_body_bottom,
    plot_caption,
    ncol = 1,
    rel_heights = c(0.25, 2, 0.8, 0.15)
  )
  
  if (out_type == "png" || out_type == "PNG") {
    grDevices::png(file = file.path(wkdir, paste0(prefix, "_pi-zscore.png")), height = 17, width = 14, units = "in", res = 300)
  } else {
    grDevices::cairo_pdf(file = file.path(wkdir, paste0(prefix, "_pi-zscore.pdf")), height = 18, width = 14)
  }
  
  print(all_plot)
  grDevices::dev.off()
}





plot_title <- function(bam_file, bed_file, genome_file, prefix, i) {
  now <- format(lubridate::now(), "%Y-%m-%d %H:%M:%S")
  
  misipi_version <- packageVersion("MiSiPi.RNA")
  
  iteration_str <- paste0(i, ")")
  
  p_title <- paste("MiSiPi Results for locus:", prefix, "(Bed file line:", iteration_str)
  p_subtitle <- paste("Bam:", bam_file, "| Bed:", bed_file, "| Genome:", genome_file)
  p_caption <- paste("Run at:", now, "with MiSiPi.RNA Version:", misipi_version)
  
  plot_details <- list()
  plot_details$title <- paste0(p_title, "\n", p_subtitle)
  #plot_details$subtitle <- p_subtitle
  plot_details$caption <- p_caption
  return(plot_details)
}

null_plot <- function(type, reason) {
  p <- ggplot2::ggplot() +
    ggplot2::annotate("text", x = 0.5, y = 0.5, label = reason, size = 5) +
    ggplot2::ggtitle(type) +
    ggplot2::theme_void() +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))
  
  return(p)
}

