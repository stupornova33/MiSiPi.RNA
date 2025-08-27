# .plot_sizes_by_strand
# plots the size distribution of reads separately by strand
# @param stranded_size_dist
# @return A plot

.plot_sizes_by_strand <- function(stranded_size_dist) {
  options(scipen = 999)

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
  
  df$first <- sub("T", "U", df$first)
  
  df$first <- factor(df$first, levels = c("A", "C", "G", "U", "N"))
  
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
    ggplot2::labs(x = "Read size", y = "Read count", title = "Read Distribution") +
    ggplot2::theme_minimal() +
    ggplot2::theme(axis.text=ggplot2::element_text(size=14),
                   plot.title = ggplot2::element_text(size = 14, hjust = 0.5))
    

  return(p)
}
