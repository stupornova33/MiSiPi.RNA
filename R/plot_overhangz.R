# function to plot the zscore for Dicer overhangs
# @param overhang1 a dataframe
# @param strand a string, "+", "-", or "none"
# @return plot

.plot_overhangz <- function(overhang1, strand) {
  # calculates Dicer signature in pairs

  # shift <- value <- proper_count <- NULL
  shift <- c(-4, -3, -2, -1, 0, 1, 2, 3, 4)

  if (strand == "+") {
    col <- "red"
    title <- "Dicer Overhang Prob (plus strand)"
  } else if (strand == "-") {
    col <- "blue"
    title <- "Dicer Overhang Prob (minus strand)"
  } else {
    col <- "black"
    title <- "siRNA Dicer Overhang Probability (dual strand)"
  }

  p <- ggplot2::ggplot(overhang1, ggplot2::aes(x = shift, y = zscore)) +
    ggplot2::geom_line(color = col, linewidth = 2) +
    ggplot2::ggtitle(title) +
    ggplot2::scale_x_continuous("Shift") +
    ggplot2::scale_y_continuous("Z-score") +
    ggplot2::theme_classic() +
    ggplot2::theme(text = ggplot2::element_text(size = 15), plot.title = ggplot2::element_text(hjust = 0.5, size = 14)) +
    ggplot2::theme(plot.margin = ggplot2::unit(c(2, 0, 0, 0), "cm"))
  return(p)
}

.plot_siRNA_overhangs_combined <- function(plus_overhangs, minus_overhangs, dual_strand_overhangs) {
  shift <- c(-4, -3, -2, -1, 0, 1, 2, 3, 4)
   
   plus_color <- "red"
   minus_color <- "blue" 
   dual_color <- "black"
   plot_title <- "siRNA Dicer Overhang Probability"
  

   
   p <- ggplot2::ggplot(plus_overhangs, ggplot2::aes(x = shift, y = zscore)) +
     ggplot2::geom_line(color = plus_color, linewidth = 2) +
     ggplot2::geom_line(data = minus_overhangs, ggplot2::aes(x = shift, y = zscore), color = minus_color, linewidth = 2) +
     ggplot2::geom_line(data = dual_strand_overhangs, ggplot2::aes(x = shift, y = zscore), color = dual_color, linewidth = 2) +
     ggplot2::ggtitle(plot_title) +
     ggplot2::scale_x_continuous("Shift") +
     ggplot2::scale_y_continuous("Z-score") +
     ggplot2::theme_classic() +
     ggplot2::theme(text = ggplot2::element_text(size = 15), plot.title = ggplot2::element_text(hjust = 0.5, size = 14)) +
     ggplot2::theme(plot.margin = ggplot2::unit(c(2, 0, 0, 0), "cm"))
  return(p) 

}

.plot_miRNA_dicer_overhang_probability <- function(plus_overhangs, minus_overhangs) {
  # shift <- value <- proper_count <- NULL
  shift <- c(-4, -3, -2, -1, 0, 1, 2, 3, 4)

  p_overhangs <- data.frame(
    shift = plus_overhangs$shift,
    zscore = plus_overhangs$zscore,
    Strand = "Plus"
  )
  
  m_overhangs <- data.frame(
    shift = minus_overhangs$shift,
    zscore = minus_overhangs$zscore,
    Strand = "Minus"
  )
  
  all_overhangs <- dplyr::bind_rows(p_overhangs, m_overhangs)
  
  p <- ggplot2::ggplot(all_overhangs, ggplot2::aes(x = shift, y = zscore, color = Strand)) +
    ggplot2::geom_line(linewidth = 2) +
    ggplot2::scale_x_continuous("Shift") +
    ggplot2::scale_y_continuous("Z-score") +
    ggplot2::scale_color_manual(values = c("Plus" = "red", "Minus" = "blue")) +
    ggplot2::ggtitle("miRNA Dicer Overhang Probability") +
    ggplot2::theme_classic() +
    ggplot2::theme(text = ggplot2::element_text(size = 15), plot.title = ggplot2::element_text(hjust = 0.5, size = 14)) +
    ggplot2::theme(plot.margin = ggplot2::unit(c(2, 0, 0, 0), "cm"))
  return(p)
}

.plot_siRNA_overhangs_combined <- function(plus_overhangs, minus_overhangs, dual_strand_overhangs) {
  shift <- c(-4, -3, -2, -1, 0, 1, 2, 3, 4)

  plus_color <- "red"
  minus_color <- "blue"
  dual_color <- "black"
  plot_title <- "siRNA Dicer Overhang Probability"
  
  
  p <- ggplot2::ggplot(plus_overhangs, ggplot2::aes(x = shift, y = zscore)) +
    ggplot2::geom_line(color = plus_color, linewidth = 2) +
    ggplot2::geom_line(data = minus_overhangs, ggplot2::aes(x = shift, y = zscore), color = minus_color, linewidth = 2) +
    ggplot2::geom_line(data = dual_strand_overhangs, ggplot2::aes(x = shift, y = zscore), color = dual_color, linewidth = 2) +
    ggplot2::ggtitle(plot_title) +
    ggplot2::scale_x_continuous("Shift") +
    ggplot2::scale_y_continuous("Z-score") +
    ggplot2::theme_classic() +
    ggplot2::theme(text = ggplot2::element_text(size = 15), plot.title = ggplot2::element_text(hjust = 0.5, size = 14)) +
    ggplot2::theme(plot.margin = ggplot2::unit(c(2, 0, 0, 0), "cm"))
  
  return(p)
}
