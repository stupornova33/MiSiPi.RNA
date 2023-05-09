#' plot_miRNA structure
#' @param f_df a data frame

#' @return plot
#' @export


plot_miRNA_struct <- function(f_df){

   text <- x <- y <- colorbins <- NULL
   pal_length <- length(unique(f_df$color_bins))
   palt <- viridis::plasma(pal_length)
   
   letters <- f_df %>% dplyr::filter(text == 'C' | text == 'U' | text == 'G' | text == 'A' )
   not_letters <- f_df %>% dplyr::filter(text == ' ' | text == '|' | text == '-')
      
   p <- ggplot2::ggplot(letters) + # ggplot2::aes(x,y, label = text)) +
      ggplot2::geom_label(ggplot2::aes(x, y, label = text, fill = color_bins), color = "white", show.legend = FALSE, label.r = ggplot2::unit(0.5, "lines"), 
                          label.padding = ggplot2::unit(0.3, "lines")) +
      ggplot2::geom_text(data = not_letters, mapping = ggplot2::aes(x, y, label = text), color = "black", fontface = "bold") +
      ggplot2::scale_colour_manual(values=c(palt))+
      ggplot2::scale_fill_manual(values = c(palt))+ 
      ggplot2::scale_y_continuous(breaks = seq(2.5,5.5, by = 1), limits = c(2.5,5.5))+ 
      ggplot2::theme(plot.margin = ggplot2::unit(c(0, 0, 0, 0), "cm")) +
      ggplot2::theme_void()

   
   return(p)
}
