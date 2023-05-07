#' takes a gff file 
#' only output is hairpin plots
#' @param gff a string path leading to the gff file
#' @param chrom_name a string
#' @param reg_start an integer
#' @param reg_stop an integer
#' 
#' 
#' @return plot
#' @export



plot_bed_annotation <- function(gff_path, chrom_name, reg_start, reg_stop){
#gff <- read.table("C:/Users/tmjar/Desktop/SharedUbuntu/hemiptera/GCF_000001215.4_Release_6_plus_ISO1_MT_genomic.gff", sep = "\t", header = FALSE) %>%
#  dplyr::filter(V3 != "exon" & V3 != 'region') %>%
#  dplyr::select(V1, V3, V4,V5, V7,V9) %>%
#  dplyr::distinct(V4,V5, .keep_all = TRUE)
  
#gff <- tidyr::separate(data = gff, col = V9, into = c('a', 'name'), sep = 'ID=')
#gff <- tidyr::separate(data = gff, col = name, into = c('feature', 'b'), sep = ';')
#gff <- gff %>% dplyr::select(-c(a, b))

#write.table(gff, file = "processed_dmel_annot.txt", sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)
#bed <- gff[2:3,]

gff <- read.table(gff_path, header = FALSE, sep = '\t')
#x coordinates will be bed coordinates
#y coordinates will be predefined
options(scipen = 999)
#reg_start <- 280100
#reg_stop <- 280300

bed <- gff %>% dplyr::filter(V1 == chrom_name) %>% 
  dplyr::filter(V3 >= reg_start & V4 <= reg_stop)

if(nrow(bed) < 1){
  g <- ggplot2::ggplot() +
    ggplot2::coord_cartesian(ylim = c(0,5), xlim = c(reg_start, reg_stop))+
    ggplot2::theme_classic() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(size = 12), axis.text.y = ggplot2::element_blank())
  
} else{ 


  ###################################################################
  #test outside of function
  coord_df <- data.frame(ids = numeric(0), xvals = numeric(0), yvals = numeric(0), midx = numeric(0), midy = numeric(0), fill_col = numeric(0))
  func <- function(bed){
    #take a data frame with multiple regions, determine coordinates of points, return as list
    start <- as.numeric(bed[1])
    stop <- as.numeric(bed[2])
    if(bed[3] == "+"){
      ymin <- 2
      ymax <- 3
      
      
      print(stop + 5)
      xvals <- c(start, start, stop, stop + 5, start)
      ids <- c(rep(bed[4], times = length(xvals)))
      yvals <- c(ymin, ymax, ymax, ymin + 0.5, ymin)
      midx <- stop + 14
      midy <- ymin + 0.5
      fill_col <- rep("red", times = length(xvals))
    } else {
      ymin <- 1
      ymax <- 2
      tipx <- start - 5
      xvals <- c(tipx, start, stop, stop, start)
      yvals <- c(ymin + 0.5, ymax, ymax, ymin, ymin)
      ids <- c(rep(bed[4], times = length(xvals)))
      midx <- stop + 14
      midy <- ymin + 0.5
      fill_col <- rep("blue", times = length(xvals))
    }
    rbind(coord_df, data.frame(ids, xvals, yvals,midx, midy, fill_col))
  }
  
  res <- apply(unname(bed), 1, func)
  plot_env <- new.env()
  
  g <- ggplot2::ggplot() +
    ggplot2::coord_cartesian(ylim = c(0,6), xlim = c(reg_start, reg_stop))
    
    my_env <- new.env(parent = emptyenv())
    
    # Package state variable for g with getter and setter functions
    my_env$g <- ggplot() +
      coord_cartesian(ylim = c(0,6), xlim = c(reg_start, reg_stop))
    
    get_g <- function() {
      my_env$g
    }
    
    set_g <- function(value) {
      old <- my_env$g
      my_env$g <- value
      invisible(old)
    }
    
    # Package state variable for res with getter function
    my_env$res <- res
    
    get_res <- function() {
      my_env$res
    }
    
    
    new_func <- function(i) {
      current_plot <- get_g()
      res_list <- get_res()
      
      new_plot <- current_plot +
        geom_polygon(data = res_list[[i]], aes(x=xvals, y=yvals), fill = res_list[[i]]['fill_col'], alpha = 0.6) +
        geom_text(data = res_list[[i]], aes(midx - 27, midy, label = ids, size = 4), color = "white") +
        theme_classic() +
        theme(axis.text.x = element_text(size = 12), axis.text.y = element_blank())
      
      set_g(new_plot)
    }
    
    
    
    lapply(seq(get_res()), new_func)
    get_g()
    
  }
 
  test <- lapply(seq(res), func)
 
  
  return(test[-1])
  
  #df("test_plot.pdf", width = 10, height = 8)
  #print(test5[-1])
  #dev.off()


}

