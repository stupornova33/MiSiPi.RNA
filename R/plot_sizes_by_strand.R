#' plots the size distribution of reads separately by strand
#'
#' @param dist A data frame made by get_stranded_read_dist function.
#' @return A plot
#' @export

plot_sizes_by_strand <- function(dist){
  options(scipen=999)
  

  p_size_dist <- dist %>% dplyr::filter(strand = "positive")
  m_size_dist <- dist %>% dplyr::filter(strand = "negative")
  
  g1 <- ggplot2::ggplot(p_size_dist, ggplot2::aes(x=width, y=count, fill=strand))+
    ggplot2::geom_bar(position = "stack", stat = "identity") +
    ggplot2::theme_classic()+
    ggplot2::labs(title = "Read Size Distribution By Strand", fill = "first nt")+
    ggplot2::ylab("sense count")+
    ggplot2::theme(plot.title = ggplot2::element_text(size = 12, hjust = 0.5))+
    ggplot2::scale_x_continuous(breaks = seq(18,32,2), expand = c(0, 0))+
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0, 0),
                                          add = c(0, 0))) +
    ggplot2::scale_fill_manual(values = "red")+
    ggplot2::theme(axis.text.x = ggplot2::element_blank(), axis.text.y = ggplot2::element_text(size=12))+
    ggplot2::theme(axis.title.x = ggplot2::element_blank(), axis.ticks.x = ggplot2::element_blank(),
                   axis.title.y = ggplot2::element_text(size = 12), legend.position = "none")
  
  
  g2 <- ggplot2::ggplot(m_size_dist, ggplot2::aes(x=width, y=count, fill=strand))+
    ggplot2::geom_bar(position = "stack", stat = "identity") +
    ggplot2::theme_classic()+
    ggplot2::ylab("antisense count")+
    #ggplot2::labs(title = "Read Size Distribution")+
    ggplot2::theme(plot.title = ggplot2::element_text(size = 12, hjust = 0.5))+
    ggplot2::scale_x_continuous(position = "top", breaks = seq(18,32,2)) +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0, 0),
                                                            add = c(0, 0))) +
    ggplot2::scale_fill_manual(values = "blue")+
    ggplot2::theme(axis.text.x = ggplot2::element_text(size=12, angle=45), axis.text.y = ggplot2::element_text(size=12))+
    ggplot2::theme(axis.title.x = ggplot2::element_blank())+
    ggplot2::theme(axis.title.y = ggplot2::element_text(size = 12), legend.position = "none")
  
  
  #gridExtra::grid.arrange(g1,g2,ncol=1,nrow = 2, heights = c(1,1), widths = c(1), lrt)
  
  p <- cowplot::plot_grid(g1, g2, ncol = 1, align = "vh", axis = "lrtb", rel_heights = c(1,1), rel_widths = c(1,1))
  
  
  
  
  
  
  
  return(p)
}
