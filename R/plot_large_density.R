#' Plots read sizes over large loci
#' Faster than geom_area but doesn't look as good
#' @param data a dataframe consisting of three columns, "Position", "Count" and "Length".
#' @param reg_start The start position of the region of interest.
#' @param reg_stop The end position of a region of interest.
#' @return Density Plot
#' @export


plot_large_density <- function(data, reg_start, reg_stop){

  #density <- ggplot2::ggplot(data, ggplot2::aes(x = Position, y = Count,
  #                                  fill = forcats::fct_reorder(Lengths, Count),
  #                                  colour = forcats::fct_reorder(Lengths, Count))) +
  #  ggplot2::geom_col(position="identity")

   density <- ggplot2::ggplot(data, ggplot2::aes(x = Position, y = Count, fill = Lengths)) +
    ggplot2::geom_col() +
    ggplot2::scale_fill_manual(values = c("yellow","red","blue", "black","yellow", "red","blue","black" )) +
    ggplot2::ggtitle("Read Density") +
    ggplot2::theme_classic()+
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))+
    ggplot2::theme(axis.text.x = ggplot2::element_text(size = 10, angle = 60, hjust = 1))


  return(density)

}
