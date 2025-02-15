# plots the coverage over an interval
# @param dt a dataframe
# @return a histogram plot

.plot_coverage <- function(dt) {
  pos <- count <- NULL
  coverage_hist <- ggplot2::ggplot(dt, ggplot2::aes(pos, count)) +
    ggplot2::geom_col(width = 1) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45)) +
    ggplot2::scale_x_continuous(breaks = round(seq(min(dt$pos), max(dt$pos), by = (max(dt$pos) - min(dt$pos)) / 5), 5)) +
    ggplot2::labs(x = "Bp position", y = "Depth") +
    ggplot2::theme_minimal()
  return(coverage_hist)
}
