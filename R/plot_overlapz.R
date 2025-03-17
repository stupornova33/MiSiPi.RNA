# plot the overlap probability
# @param z_df a data frame
# @return a ggplot object

.plot_overlapz <- function(z_df) {
  Z_score <- Overlap <- NULL
  p <- ggplot2::ggplot(z_df, ggplot2::aes(x = Overlap, y = zscore)) +
    ggplot2::geom_line(color = "black", linewidth = 1.5) +
    ggplot2::theme_classic() +
    ggplot2::scale_x_continuous("Overlap", breaks = seq(0, 32, 5), labels = seq(0, 32, 5)) +
    ggplot2::theme(axis.title.x = ggplot2::element_text(size = 14), axis.title.y = ggplot2::element_text(size = 14)) +
    ggplot2::theme(plot.title = ggplot2::element_text(size = 14, hjust = 0.5)) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(size = 12), axis.text.y = ggplot2::element_text(size = 12)) +
    ggplot2::ggtitle("Overlap Probability")

  return(p)
}
