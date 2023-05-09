#' plots the zscore for the 2nt overhang
#' takes a table of overhang zscores
#' only output is overhang plots
#' 
#' @param overhangs a table
#' 
#' @return plot
#' @export

overhang_z_plot <- function(overhangs){
   proper_count <- NULL
   pos_shift <- c(-4, -3, -2, -1, 0, 1, 2, 3, 4)
   p <- ggplot2::ggplot(overhangs, ggplot2::aes(x = pos_shift, y = scale(proper_count))) +
      ggplot2::scale_x_continuous("Shift")+
      ggplot2::scale_y_continuous("Z-score")+
      ggplot2::theme_classic() +
      ggplot2::theme(text = ggplot2::element_text(size = 10))+
      ggplot2::geom_line(color= "black", linewidth = 1.25)
      
   p
   return(p)
   
   
}
