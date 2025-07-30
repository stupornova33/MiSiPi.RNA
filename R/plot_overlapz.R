# plot the overlap probability
# @param z_df a data frame
# @return a ggplot object

# Previously named .plot_overlapz
.plot_piRNA_overlap_probability <- function(z_df) {
  Z_score <- Overlap <- NULL
  p <- ggplot2::ggplot(z_df, ggplot2::aes(x = Overlap, y = zscore)) +
    ggplot2::geom_line(color = "black", linewidth = 1.5) +
    ggplot2::theme_classic() +
    ggplot2::scale_x_continuous("Overlap", breaks = seq(0, 32, 5), labels = seq(0, 32, 5)) +
    ggplot2::theme(axis.title.x = ggplot2::element_text(size = 14), axis.title.y = ggplot2::element_text(size = 14)) +
    ggplot2::theme(plot.title = ggplot2::element_text(size = 14, hjust = 0.5)) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(size = 12), axis.text.y = ggplot2::element_text(size = 12)) +
    ggplot2::ggtitle("piRNA Overlap Probability")

  return(p)
}

.plot_miRNA_overlap_probability <- function(plus_z_df, minus_z_df) {
  zscore <- Overlap <- NULL
  
  
  p_overlaps <- data.frame(
    Overlap = plus_z_df$Overlap,
    zscore = plus_z_df$zscore,
    Strand = "Plus"
  )
  
  m_overlaps <- data.frame(
    Overlap = minus_z_df$Overlap,
    zscore = minus_z_df$zscore,
    Strand = "Minus"
  )
  
  all_overlaps <- dplyr::bind_rows(p_overlaps, m_overlaps)
  
  p <- ggplot2::ggplot(all_overlaps, ggplot2::aes(x = Overlap, y = zscore, color = Strand)) +
    ggplot2::geom_line(linewidth = 1.5) +
    ggplot2::scale_x_continuous("Overlap") +
    ggplot2::scale_y_continuous("Z-score") +
    ggplot2::scale_color_manual(values = c("Plus" = "red", "Minus" = "blue")) +
    ggplot2::ggtitle("miRNA Overlap Probability") +
    ggplot2::theme_classic() +
    ggplot2::theme(text = ggplot2::element_text(size = 15), plot.title = ggplot2::element_text(hjust = 0.5, size = 14)) +
    ggplot2::theme(plot.margin = ggplot2::unit(c(2, 0, 0, 0), "cm"))
  
  return(p)
}
