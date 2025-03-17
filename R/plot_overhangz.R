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
    title <- "Dicer Overhang Prob (siRNA)."
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
