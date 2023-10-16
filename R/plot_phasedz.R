#' phasing function
#' processes reads according to phased piRNA algorithm
#' plots output
#' takes data frame
#' returns plots
#'
#' @param df1 a data frame
#' @param strand a string, "+" or "-"
#' @param df2 a data frame
#' @return plots

#' @export


plot_phasedz <- function(df1, strand = NULL, df2 = NULL) {
  phased_num <- phased_dist <- value <- phased_z <- phased26_dist <- phased26_z <- phased_dist1 <- phased_z1 <- phased_dist2 <- phased_z2 <- NULL

  if(!is.null(strand)){
     if(strand == "+"){
        title <- "Phasing Prob (Plus Strand)"
     } else if(strand == "-"){
         title <- "Phasing Prob (Minus Strand)"
      }
  } else {
      title <- "Phasing Prob"
  }


  if(!is.null(df2)){
         df1 <- df1 %>% dplyr::mutate(strand = "plus")
         df2 <- df2 %>% dplyr::mutate(strand = "minus")

         all_dat <- rbind(df1, df2) %>% dplyr::select(-c(phased_num))
         data_long <- reshape2::melt(all_dat, id = c("strand", "phased_dist"))
         p <- ggplot2::ggplot(data_long, ggplot2::aes(x = phased_dist, y = value, color = strand)) +
            ggplot2::geom_line(linewidth = 1.25) +
            ggplot2::ggtitle(title)+
            ggplot2::scale_color_manual(values = c("blue", "red"))+
            ggplot2::scale_y_continuous("Z-score")+
            ggplot2::scale_x_continuous("3' to 5' Distance", labels = seq(1,50, by = 5), breaks = seq(1,50, by=5)) +
            ggplot2::theme_classic() +
            ggplot2::theme(text = ggplot2::element_text(size = 12), plot.title = ggplot2::element_text(size = 14, hjust = 0.5)) +
            ggplot2::theme(axis.text.x = ggplot2::element_text(size = 12, angle = 45))+
            ggplot2::theme(plot.margin = ggplot2::unit(c(2, 0, 0, 0), "cm"))

   } else {
      if(ncol(df1) > 3 & colnames(df1)[1] != "phased_dist1"){
          p <- ggplot2::ggplot(df1, ggplot2::aes(x = phased_dist))+
               ggplot2::scale_y_continuous("Z-score")+
               ggplot2::ggtitle(title)+
               ggplot2::scale_x_continuous("3' to 5' Distance", labels = seq(1,50, by = 5), breaks = seq(1,50, by=5)) +
               ggplot2::geom_line(ggplot2::aes(y = phased_z), linewidth = 1.25, color = "black") +
               ggplot2::geom_line(ggplot2::aes(x = phased26_dist, y = phased26_z), linewidth = 1.5, color = "red") +ggplot2::theme_classic() +
               ggplot2::theme(axis.text = ggplot2::element_text(size = 12), plot.title = ggplot2::element_text(size = 14, hjust = 0.5))+
               ggplot2::theme(axis.text.x = ggplot2::element_text(size = 12, angle = 45))+
               ggplot2::theme(text = ggplot2::element_text(size = 12))

      } else if(ncol(df1) > 3 & colnames(df1)[1] == "phased_dist1"){
         p <- ggplot2::ggplot(df1, ggplot2::aes(x = phased_dist1))+
            ggplot2::scale_y_continuous("Z-score")+
            ggplot2::ggtitle(title)+
            ggplot2::scale_x_continuous("3' to 5' Distance", labels = seq(1,50, by = 5), breaks = seq(1,50, by=5)) +
            ggplot2::geom_line(ggplot2::aes(y = phased_z1), linewidth = 1.25, color = "blue") +
            ggplot2::geom_line(ggplot2::aes(x = phased_dist2, y = phased_z2), linewidth = 1.5, color = "red") +ggplot2::theme_classic() +
            ggplot2::theme(axis.text = ggplot2::element_text(size = 12), plot.title = ggplot2::element_text(size = 14, hjust = 0.5))+
            ggplot2::theme(axis.text.x = ggplot2::element_text(size = 12, angle = 45))+
            ggplot2::theme(text = ggplot2::element_text(size = 12))

      } else {
      p <- ggplot2::ggplot(df1, ggplot2::aes(x = phased_dist))+

         ggplot2::scale_y_continuous("Z-score")+
         ggplot2::ggtitle(title)+
         ggplot2::scale_x_continuous("3' to 5' Distance", labels = seq(1,50, by = 5), breaks = seq(1,50, by=5)) +
         ggplot2::geom_line(ggplot2::aes(y = phased_z), linewidth = 1.25, color = "black") +
         ggplot2::theme_classic() +
         ggplot2::theme(axis.text = ggplot2::element_text(size = 12), plot.title = ggplot2::element_text(size = 14, hjust = 0.5))+
         ggplot2::theme(axis.text.x = ggplot2::element_text(size = 12, angle = 45))+
         ggplot2::theme(text = ggplot2::element_text(size = 12))
      }
   }

   return(p)
}

