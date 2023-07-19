#' function to plot densities of different read sizes
#' @param all_df a data frame
#' @param reg_start a whole number
#' @param reg_stop a whole number
#' @return density

#' @export


plot_density <- function(all_df, reg_start, reg_stop){
   position <- pos_26_32 <- neg_26_32 <- pos_23_25 <- neg_23_25 <- pos_20_22 <- neg_20_22 <- pos_18_19 <- neg_18_19 <- NULL
   length <- reg_stop - reg_start

   density <- ggplot2::ggplot(all_df, ggplot2::aes(x = position)) +
      ggplot2::geom_area(ggplot2::aes(y=pos_26_32), fill = "black") +
      ggplot2::geom_area(ggplot2::aes(y=neg_26_32), fill = "black") +
      ggplot2::geom_area(ggplot2::aes(y=pos_23_25), fill = "blue") +
      ggplot2::geom_area(ggplot2::aes(y=neg_23_25), fill = "blue") +
      ggplot2::geom_area(ggplot2::aes(y=pos_20_22), fill = "red") +
      ggplot2::geom_area(ggplot2::aes(y=neg_20_22), fill = "red") +
      ggplot2::geom_area(ggplot2::aes(y=pos_18_19), fill = "yellow") +
      ggplot2::geom_area(ggplot2::aes(y=neg_18_19), fill = "yellow") +
      ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0,1))) +
      ggplot2::scale_x_continuous(breaks = seq(reg_start, reg_stop, by = round(length/10)), expand = ggplot2::expansion(mult = c(0,0.06))) +
      ggplot2::ggtitle("Read Density")+
      ggplot2::theme_void() +
      ggplot2::theme(axis.text.x= ggplot2::element_text(size = 10, angle = 60),
            axis.ticks.x= ggplot2::element_line(),
            axis.line.x = ggplot2::element_line(),
            plot.title = ggplot2::element_text(size=14, hjust = 0.5))+
      ggplot2::scale_fill_discrete(breaks = c("26-32", "23-25", "20-22", "18-19"))
   return(density)

}

