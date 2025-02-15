# function to plot densities of different read sizes
# @param data a data frame
# @param reg_start a whole number
# @param reg_stop a whole number
# @return density

.plot_density <- function(data, reg_start, reg_stop) {
  density <- ggplot2::ggplot(data, ggplot2::aes(x = pos, y = count, group = size)) +
    ggplot2::geom_area(data = subset(data, strand == "pos"), ggplot2::aes(color = size, fill = size)) +
    ggplot2::geom_area(data = subset(data, strand == "neg"), ggplot2::aes(color = size, fill = size)) +
    # ggplot2::scale_x_continuous(limits = c(reg_start - 10, reg_stop + 10), expand = c(0, 0)) +
    # ggplot2::scale_y_continuous(expand = c(0, 0)) +
    # data is organized: pos_26_32, neg_26-32, pos_23_25, neg_23-25, pos_20-22, neg_20-22, pos_18-19, neg_18
    ggplot2::scale_fill_manual(values = c("yellow", "red", "blue", "black")) +
    ggplot2::scale_color_manual(values = c("yellow", "red", "blue", "black")) +
    ggplot2::ggtitle("Read Density") +
    ggplot2::theme_classic() +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5)) +
    # ggplot2::theme(legend.key.size = ggplot2::unit(0.3, 'cm'))+
    ggplot2::theme(axis.text.x = ggplot2::element_text(size = 10, angle = 60, hjust = 1))
    #              axis.ticks.x = ggplot2::element_line(),
    #              axis.line.x = ggplot2::element_line(),
    #              panel.grid = ggplot2::element_blank(),
    #              panel.border = ggplot2::element_blank())

  return(density)
}
