#' plots the dot bracket structure above the nucleotide sequence
#' @param vienna a string
#' @param sequence a string
#' @return a plot
#' @export

plot_text <- function(vienna, sequence){
  
   x <- y <- text <- NULL
   v <- as.vector(unlist(strsplit(vienna, "")))
   c <- as.vector(unlist(strsplit(sequence, "")))
   
   v_df <- data.frame("x" = c(1:length(v)), "y" = c(5), "text" = v)
   c_df <- data.frame("x" = c(1:length(c)), "y" = c(4.6), "text" = c)
   
   f_df <- rbind(v_df, c_df)
   
   g <- ggplot2::ggplot(f_df, ggplot2::aes(x,y, label = text)) +
      ggplot2::geom_text() + 
      ggplot2::xlim(0, max(length(v), length(c) + 2)) + 
      
      ggplot2::scale_y_continuous(breaks = seq(2,5, by = 1), limits = c(2,5))+ 
      #ggbreak::scale_y_cut(breaks=c(4), which=c(1), scales=c(3))+ 
      ggplot2::theme(plot.margin = ggplot2::unit(c(0, 0, 0, 0), "cm")) +
      ggplot2::theme_void()
   
   return(g)
}
