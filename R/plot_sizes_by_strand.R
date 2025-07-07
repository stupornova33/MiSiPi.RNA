# .plot_sizes_by_strand
# plots the size distribution of reads separately by strand
# @param chrom_name a string
# @param reg_start a numerical value
# @param reg_stop a numerical value
# @param bam_file a string
# @param libsize a numeric value
# @return A plot

.plot_sizes_by_strand <- function(wkdir, stranded_size_dist, chrom_name, reg_start, reg_stop) {
  options(scipen = 999)
  prefix <- .get_region_string(chrom_name, reg_start, reg_stop)

  
  p_size_dist <- stranded_size_dist %>% dplyr::filter(strand == "plus")
  m_size_dist <- stranded_size_dist %>% dplyr::filter(strand == "minus")

  m_size_dist$count <- m_size_dist$count * -1

  empty <- data.frame(width = c(seq(18, 32, by = 1)), strand = "plus", count = numeric(15L))
  p_dist <- merge(empty, p_size_dist, by = "width", all.x = TRUE) %>%
    dplyr::select(-c(count.x, strand.y)) %>%
    dplyr::rename("total_count" = "count.y")
  p_dist[is.na(p_dist)] <- 0

  empty <- data.frame(width = c(seq(18, 32, by = 1)), strand = "minus", count = numeric(15L))
  m_dist <- merge(empty, m_size_dist, by = "width", all.x = TRUE) %>%
    dplyr::select(-c(count.x, strand.y)) %>%
    dplyr::rename("total_count" = "count.y")
  m_dist[is.na(m_dist)] <- 0
  m_dist$total_count <- m_dist$total_count * -1

  df <- rbind(p_dist, m_dist)
  df$first <- factor(df$first, levels = c("A", "C", "G", "T", "N"))
  
  letter_levels <- df %>% dplyr::select(first) %>% dplyr::distinct()
  if(nrow(letter_levels) == 4){
    pal <- c( "darkgreen", "red","blue", "yellow")
  } else {
    pal <- c( "darkgreen", "red","blue", "darkgrey","yellow")
  }

  p <- ggplot2::ggplot(df, ggplot2::aes(x = factor(width), y = total_count, fill = first)) +
    ggplot2::geom_bar(stat = "identity", position = "stack") +
    ggplot2::scale_y_continuous(labels = abs) +
    ggplot2::scale_fill_manual(values = c( "darkgreen", "red","blue", "yellow", "darkgrey")) +
    ggplot2::labs(x = "Read size", y = "Read count") +
    ggplot2::theme(axis.text=ggplot2::element_text(size=14))+
    ggplot2::theme_minimal()

  #grDevices::png(file = file.path(wkdir, paste0(prefix, "_sizes_by_strand.png")), height = 10, width = 10, units = "in", res = 300)
  #print(p)
  #grDevices::dev.off()

  return(p)
}
