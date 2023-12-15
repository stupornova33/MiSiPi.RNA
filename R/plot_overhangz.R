#' function to plot the zscore for Dicer overhangs
#' @param overhang1 a dataframe
#' @param strand a string, "+", "-", or "none"
#' @return plot

#' @export


plot_overhangz <- function(overhang1, strand){
   #calculates Dicer signature in pairs

   #shift <- value <- proper_count <- NULL
   shift <- c(-4, -3, -2, -1, 0, 1, 2, 3, 4)

   if(strand == "+"){
     col <- "red"
   } else if(strand == "-") {
     col <- "blue"
   } else {
     col <- "black"
   }
  # if(!is.null(overhang2)){
  #    overhang1 <- overhang1 %>% dplyr::mutate(strand = "plus")
  #    overhang2 <- overhang2 %>% dplyr::mutate(strand = "minus")

  #    all_dat <- rbind(overhang1, overhang2) #%>% dplyr::select(-c(improper_count, proper_count))
  #    data_long <- reshape2::melt(all_dat, id = c("strand", "shift"))

   #p <- ggplot2::ggplot(data_long, ggplot2::aes(x = shift, y = value, color = strand)) +
  #       ggplot2::geom_line(linewidth = 2) +
  #       ggplot2::scale_x_continuous("Shift")+
  #       ggplot2::scale_color_manual(values = c("red", "blue"))+
  #       ggplot2::scale_y_continuous("Z-score")+
  #       ggplot2::ggtitle("Dicer Overhang Probability By Strand")+
  #       ggplot2::theme_classic() +
  #       ggplot2::theme(text = ggplot2::element_text(size = 15), plot.title = ggplot2::element_text(hjust = 0.5, size = 14)) +
  #       ggplot2::theme(plot.margin = ggplot2::unit(c(2, 0, 0, 0), "cm"))
  # } else {

   p <- ggplot2::ggplot(overhang1, ggplot2::aes(x = shift, y = zscore)) +
      ggplot2::geom_line(color= col, linewidth = 2)+
      ggplot2::ggtitle("Dicer Overhang Probability (siRNA)")+
      ggplot2::scale_x_continuous("Shift")+
      ggplot2::scale_y_continuous("Z-score")+

      ggplot2::theme_classic() +
      ggplot2::theme(text = ggplot2::element_text(size = 15), plot.title = ggplot2::element_text(hjust = 0.5, size = 14)) +
      ggplot2::theme(plot.margin = ggplot2::unit(c(2, 0, 0, 0), "cm"))

   #}
   return(p)
}
