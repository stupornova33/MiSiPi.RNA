#' plots the size distribution of reads
#'
#' @param dist a data table
#' @return plots
#' @export

plot_sizes <- function(dist) {
  options(scipen=999)

  num_nt <- dist %>% dplyr::select(c(first)) %>% dplyr::distinct()

  if (nrow(num_nt) == 4) {
    man_pal <- c('green4', 'red3', 'blue4', "yellow3")
  } else {
    man_pal <- c('green4', 'red3', 'blue4', 'grey', "yellow3")
  }

  width <- count <- first <- NULL
  p <- ggplot2::ggplot(dist, ggplot2::aes(x=width, y=count, fill=first)) +
    ggplot2::geom_bar(stat = "identity") +
    ggplot2::theme_classic() +
    ggplot2::labs(title = "Read Size Distribution") +
    ggplot2::scale_x_continuous(breaks = seq(18,32,2)) +
    ggplot2::scale_fill_manual(values = man_pal) +
    ggplot2::theme(plot.title = ggplot2::element_text(size = 16, hjust = 0.5),
                   #plot.subtitle = ggplot2::element_text(size = 12, hjust = 0.5),
                   axis.text.x = ggplot2::element_text(size=12, angle=45, vjust = 0.5),
                   axis.text.y = ggplot2::element_text(size=12),
                   axis.title.y = ggplot2::element_text(size = 12, margin = ggplot2::margin(r = 5)),
                   axis.title.x = ggplot2::element_text(size = 12))

  return(p)
}
