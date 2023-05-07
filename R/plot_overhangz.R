#' function to plot the zscore for Dicer overhangs
#' @param overhang1 a dataframe
#' @param overhang2 a dataframe
#' @return plot

#' @export


plot_overhangz <- function(overhang1, overhang2 = NULL){
   #calculates Dicer signature in pairs
   pos_shift <- c(-4, -3, -2, -1, 0, 1, 2, 3, 4)
   
   if(!is.null(overhang2)){
      overhang1 <- overhang1 %>% dplyr::mutate(type = "plus")
      overhang2 <- overhang2 %>% dplyr::mutate(type = "minus")

      all_dat <- rbind(overhang1, overhang2) %>% dplyr::select(-c(improper_count, proper_count))
      data_long <- reshape2::melt(all_dat, id = c("type", "shift"))
      
   p <- ggplot2::ggplot(data_long, ggplot2::aes(x = shift, y = value, color = type)) +  
         ggplot2::geom_line(size = 2) +
         ggplot2::scale_x_continuous("Shift")+
         ggplot2::scale_color_manual(values = c("red", "blue"))+
         ggplot2::scale_y_continuous("Z-score")+
         
         ggplot2::theme_classic() +
         ggplot2::theme(text = ggplot2::element_text(size = 15)) +
         ggplot2::theme(plot.margin = ggplot2::unit(c(2, 0, 0, 0), "cm"))
   } else {

   p <- ggplot2::ggplot(overhang1, ggplot2::aes(x = pos_shift, y = scale(proper_count))) +
      ggplot2::geom_line(color= "black", size = 2)+
      
      ggplot2::scale_x_continuous("Shift")+
      ggplot2::scale_y_continuous("Z-score")+
      
      ggplot2::theme_classic() +
      ggplot2::theme(text = ggplot2::element_text(size = 15)) +
      ggplot2::theme(plot.margin = ggplot2::unit(c(2, 0, 0, 0), "cm"))
      
   }
   return(p)
}
