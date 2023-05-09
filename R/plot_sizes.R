#' plots the size distribution of reads
#' 
#' @param dist a data table
#' @return plots
#' @export

plot_sizes <- function(dist){
   options(scipen=999)
   width <- count <- first <- NULL
   p <- ggplot2::ggplot(dist, ggplot2::aes(x=width, y=count, fill=first))+ 
      ggplot2::geom_bar(stat = "identity") +
      ggplot2::theme_classic()+
      ggplot2::scale_x_continuous(breaks = seq(18,32,2))+
      ggplot2::scale_fill_manual(values = c('green4', 'red3', 'blue4', 'yellow3'))+
      ggplot2::theme(axis.text.x = ggplot2::element_text(size=13, angle=45), axis.text.y = ggplot2::element_text(size=13))+
      ggplot2::theme(axis.title.x = ggplot2::element_text(size = 13))+
      ggplot2::theme(axis.title.y = ggplot2::element_text(size = 13))

   return(p)
}
