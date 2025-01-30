#' plots the size distribution of reads separately by strand
#' @param chrom_name a string
#' @param reg_start a numerical value
#' @param reg_stop a numerical value
#' @param bam_file a string
#' @param libsize a numeric value
#' @return A plot
#' @export

plot_sizes_by_strand <- function(chrom_name, reg_start, reg_stop, bam_file, libsize){
  options(scipen=999)
  prefix <- get_region_string(chrom_name, reg_start, reg_stop)
  print(prefix)

  # use Rsamtools to process the bam file
  bam_obj <- OpenBamFile(bam_file, "logfile.txt")
  bam_header <- Rsamtools::scanBamHeader(bam_obj)
  chr_name <- names(bam_header[['targets']])
  chr_length <- unname(bam_header[['targets']])
  bam_header <- NULL

  cat(file = paste0(wkdir, logfile), paste0("chrom_name: ", chrom_name, " reg_start: ", reg_start - 1, " reg_stop: ", reg_stop - 1, "\n"), append = TRUE)
  cat(file = paste0(wkdir, logfile), "Filtering forward and reverse reads by length\n", append = TRUE)

  # extract reads by strand
  # this creates a list object
  chromP <- getChrPlus(bam_obj, chrom_name, reg_start, reg_stop)
  chromM <- getChrMinus(bam_obj, chrom_name, reg_start, reg_stop)


  # turn the list object into a more useable data frame and filter reads by length,
  # bam only contains pos and width, need to add an end column
  cat(file = paste0(wkdir, logfile), "Making Forward DT\n", append = TRUE)
  forward_dt <- data.table::setDT(make_si_BamDF(chromP)) %>%
    subset(width <= 32 & width >= 18) %>%
    dplyr::rename(start = pos) %>%
    dplyr::mutate(end = start + width - 1) %>%
    dplyr::group_by_all() %>%
    # get the number of times a read occurs
    dplyr::summarize(count = dplyr::n())
  forward_dt$strand <- "positive"

  cat(file = paste0(wkdir, logfile), "Making Reverse DT\n", append = TRUE)
  reverse_dt <- data.table::setDT(make_si_BamDF(chromM)) %>%
    subset(width <= 32 & width >= 18) %>%
    dplyr::rename(start = pos) %>%
    dplyr::mutate(end = start + width - 1) %>%
    dplyr::group_by_all() %>%
    dplyr::summarize(count = dplyr::n())
  reverse_dt$strand <- "negative"

  size_dist <- dplyr::bind_rows(forward_dt, reverse_dt) %>%
    dplyr::group_by(width, strand) %>%
    dplyr::summarise(count = sum(count))

  #normalize counts to library size

  size_dist <- size_dist %>% dplyr::mutate(norm_count = count/(libsize/1000000))
  size_dist$norm_count <- ceiling(size_dist$norm_count)
  #output_readsize_dist(size_dist, prefix, wkdir, strand = NULL, "siRNA")

  #dist <- .get_read_size_dist(forward_dt, reverse_dt)

  p_size_dist <- size_dist %>% dplyr::filter(strand == "positive")
  m_size_dist <- size_dist %>% dplyr::filter(strand == "negative")

  m_size_dist$count <- m_size_dist$count*-1

  empty <- data.frame(width = c(seq(18,32, by = 1)), strand = "positive", count = numeric(15L))
  p_dist <- merge(empty, p_size_dist, by = "width", all.x = TRUE) %>% select(-c(count.x, strand.y)) %>%
    dplyr::rename( "total_count"= "count.y")
  p_dist[is.na(p_dist)] <- 0

  empty <- data.frame(width = c(seq(18,32, by = 1)), strand = "negative", count = numeric(15L))
  m_dist <- merge(empty, m_size_dist, by = "width", all.x = TRUE) %>% select(-c(count.x, strand.y)) %>%
    dplyr::rename("total_count"="count.y")
  m_dist[is.na(m_dist)] <- 0
  m_dist$norm_count <- m_dist$norm_count * -1


  m_dist <- m_dist %>% tidyr::pivot_longer(c(total_count, norm_count), names_to = "count_type")
  p_dist <- p_dist %>% tidyr::pivot_longer(c(total_count, norm_count), names_to = "count_type")



   g1 <- ggplot2::ggplot(p_dist, ggplot2::aes(x=width, y=value, fill=count_type))+
    ggplot2::geom_bar(position = "dodge", stat = "identity") +
    ggplot2::theme_classic()+
    ggplot2::labs(title = "Read Size Distribution By Strand")+
    ggplot2::ylab("sense count")+
    ggplot2::theme(plot.title = ggplot2::element_text(size = 12, hjust = 0.5))+
    ggplot2::scale_x_continuous(breaks = seq(18,32,2))+
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0, 0),
                                          add = c(0, 0))) +
    ggplot2::scale_fill_manual(values = c("red", "black"))+
    ggplot2::theme(axis.text.y = ggplot2::element_text(size=12)) +
    ggplot2::theme(axis.title.x = ggplot2::element_blank()) #axis.ticks.x = ggplot2::element_blank(),
    #               axis.title.y = ggplot2::element_text(size = 12), legend.position = "none")


  g2 <- ggplot2::ggplot(m_dist, ggplot2::aes(x=width, y=value, fill=count_type))+
    ggplot2::geom_bar(position = "dodge", stat = "identity") +
    ggplot2::theme_classic()+
    ggplot2::ylab("antisense count")+
    #ggplot2::labs(title = "Read Size Distribution")+
    ggplot2::theme(plot.title = ggplot2::element_text(size = 12, hjust = 0.5))+
    ggplot2::scale_x_continuous(position = "top", breaks = seq(18,32,2)) +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0, 0),
                                                            add = c(0, 0))) +
    ggplot2::scale_fill_manual(values = c("blue", "black"))+
    ggplot2::theme(axis.text.y = ggplot2::element_text(size=12))+
    ggplot2::theme(axis.text.x = ggplot2::element_blank(), axis.title.x = ggplot2::element_blank())+
    ggplot2::theme(axis.title.y = ggplot2::element_text(size = 12))


  #gridExtra::grid.arrange(g1,g2,ncol=1,nrow = 2, heights = c(1,1), widths = c(1), lrt)

  p <- cowplot::plot_grid(g1, g2, ncol = 1, align = "vh", axis = "lrtb", rel_heights = c(1,1), rel_widths = c(1,1))

  grDevices::png(file = paste0(prefix, "_sizes.png"), height = 13, width = 13, units = "in", res = 300)
  print(p)
  grDevices::dev.off()




  return(p)
}
