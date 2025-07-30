# Plots phasing probability
# Used in siRNA/dual_strand_hairpin
#
# Parameters:
# [df] - A data frame containing phasing distances and Z Scores
# [strand] - A character, eitehr "+" or "-"
#
# Returns a phasing probability plot

.plot_hp_phasedz <- function(df, strand) {
  phased_dist <- phased_z <- NULL

  if (strand == "+") {
    title <- "Phasing Prob. (plus_strand)"
    col <- "red"
  } else {
    title <- "Phasing Prob. (minus strand)"
    col <- "blue"
  }

  p <- ggplot2::ggplot(df, ggplot2::aes(x = phased_dist)) +
    ggplot2::geom_line(ggplot2::aes(y = phased_z), linewidth = 1.25, color = col) +
    ggplot2::scale_x_continuous("3' to 5' Distance", labels = seq(1, 50, by = 5), breaks = seq(1, 50, by = 5)) +
    ggplot2::scale_y_continuous("Z-score") +
    ggplot2::ggtitle(title) +
    ggplot2::theme_classic() +
    ggplot2::theme(axis.text = ggplot2::element_text(size = 12), plot.title = ggplot2::element_text(size = 14, hjust = 0.5)) +
    ggplot2::theme(text = ggplot2::element_text(size = 12))

  return(p)
}

.plot_siRNA_hp_phasing_probability_combined <- function(plus_phased_table, minus_phased_table) {
  phased_dist <- phased_z <- NULL

  plus_color <- "red"
  minus_color <- "blue"
  plot_title <- "siRNA Phasing Probability"

  p <- ggplot2::ggplot(plus_phased_table, ggplot2::aes(x = phased_dist, y = phased_z)) +
    ggplot2::geom_line(linewidth = 1.25, color = plus_color) +
    ggplot2::geom_line(data = minus_phased_table, ggplot2::aes(x = phased_dist, y = phased_z), linewidth = 1.25, color = minus_color) +
    ggplot2::scale_x_continuous("3' to 5' Distance", labels = seq(1, 50, by = 5), breaks = seq(1, 50, by = 5)) +
    ggplot2::scale_y_continuous("Z-score") +
    ggplot2::ggtitle(plot_title) +
    ggplot2::theme_classic() +
    ggplot2::theme(axis.text = ggplot2::element_text(size = 12), plot.title = ggplot2::element_text(size = 14, hjust = 0.5)) +
    ggplot2::theme(text = ggplot2::element_text(size = 12))

  return(p)
}
