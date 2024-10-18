#' takes a gtf file
#' outputs a plot
#' @param gtf_file a string path leading to the gtf file
#' @param chrom_name a string
#' @param reg_start an integer
#' @param reg_stop an integer
#'
#'
#' @return plot
#' @export



plot_gtf <- function(gtf_file, chrom_name, reg_start, reg_stop){
   gtf <- utils::read.table(gtf_file, header = FALSE, sep = '\t')
   #x coordinates will be gtf coordinates
   #y coordinates will be predefined
   options(scipen = 999)

   df <- data.frame(V1 = chrom_name, start = reg_start, end = reg_stop)
   idx <- vector()
   get_matches <- function(i){

      gtf <- gtf %>% dplyr::filter(gtf$V1 == df$V1)
      rng <- seq(df$start, df$end)
      idx <- which(gtf$V4 %in% rng | gtf$V5 %in% rng)
      return(idx)
   }

   idx <- unlist(lapply(seq(nrow(df)), get_matches))


   gtf <- gtf[idx,] %>% dplyr::distinct()


   ###########################################
   ## MAKE ENVIRONMENT TO STORE PLOT OBJECT ##
   ###########################################
   plot_env <- new.env(parent = emptyenv())

   plot_env$g <- ggplot2::ggplot()
   #ggplot2::coord_cartesian(ylim = c(0,6), xlim = c(reg_start, reg_stop))


   get_g <- function() {
      plot_env$g
   }

   set_g <- function(value) {
      old <- plot_env$g
      plot_env$g <- value
      invisible(old)
   }

   get_coordinates <- function() {
      plot_env$coordinates
   }
   length <- reg_stop - reg_start
   if (nrow(gtf) < 1) {
      print('nrow gtf < 1')
      current_plot <- get_g()
      print('making new_plot.')
      new_plot <- current_plot +
         ggplot2::coord_cartesian(ylim = c(0,nrow(gtf) + 3), xlim = c(reg_start, reg_stop))
         #ggplot2::scale_x_continuous(breaks = seq(reg_start, reg_stop, by = round(length/10)), expand = ggplot2::expansion(mult = c(0,0.06)))
      print('setting new_plot')
      set_g(new_plot)
   } else {
      current_plot <- get_g()
      new_plot <- current_plot +
         ggplot2::coord_cartesian(ylim = c(0,nrow(gtf) + 3), xlim = c(reg_start, reg_stop))
         #ggplot2::scale_x_continuous(breaks = seq(reg_start, reg_stop, by = round(length/10)), expand = ggplot2::expansion(mult = c(0,0.06)))
      set_g(new_plot)
      coord_df <- data.frame(ids = numeric(0), xvals = numeric(0), yvals = numeric(0), midx = numeric(0), midy = numeric(0), fill_col = numeric(0))
      calculate_coordinates <- function(gtf) {
         #take a data frame with multiple regions, determine coordinates of points, return as list
         start <- as.numeric(gtf[4])
         stop <- as.numeric(gtf[5])
         if (gtf[7] == "+") {
            ymin <- 3
            ymax <- 4
            perc_length <- (stop - start)*0.15
            xvals <- c(start, start, stop, stop + perc_length, stop)
            ids <- c(rep(gtf[3], times = length(xvals)))
            yvals <- c(ymin, ymax, ymax, ymin + 0.5, ymin)
            #midx <- start + 14
            midx <- ((stop - start)/2) + start
            midy <- ymin + 0.5
            fill_col <- rep("red", times = length(xvals))
         } else {
            ymin <- 1
            ymax <- 2
            perc_length <- (stop - start)*0.05
            xvals <- c(start, start - perc_length, start, stop, stop)
            yvals <- c(ymin, ymin + 0.5, ymax, ymax, ymin)
            ids <- c(rep(gtf[3], times = length(xvals)))
            midx <- ((stop - start)/2) + start
            #midx <- start + 14
            midy <- ymin + 0.5
            fill_col <- rep("blue", times = length(xvals))
         }
         rbind(coord_df, data.frame(ids, xvals, yvals, midx, midy, fill_col))
      }

      plot_env$coordinates <- apply(gtf, 1, calculate_coordinates)

      append_plots <- function(i) {
        xvals <- yvals <- midx <- midy <- ids <- NULL
         current_plot <- get_g()
         coordinates <- get_coordinates()
         new_plot <- current_plot +
            ggplot2::geom_polygon(data = coordinates[[i]], ggplot2::aes(x=xvals, y=yvals+(i-1)), fill = coordinates[[i]]['fill_col'], alpha = 0.6) +
            ggplot2::geom_text(data = coordinates[[i]], ggplot2::aes(midx, midy+(i-1), label = ids, size = 3), color = "black") +
            ggplot2::theme_classic() +
            ggplot2::theme(axis.text.y = ggplot2::element_blank(), axis.text.x = ggplot2::element_text(size = 12, angle = 45, hjust = 1)) +
            ggplot2::theme(axis.ticks.y = ggplot2::element_blank())+
            ggplot2::theme(legend.position = "none")+
            ggplot2::labs(y = "", x = "")
         set_g(new_plot)
      }
      lapply(seq(get_coordinates()), append_plots)


   }
   print('setting g')
   g <- get_g()
   return(g)
}

