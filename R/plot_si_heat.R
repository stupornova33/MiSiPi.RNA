# plot the si heatmaps
# takes a matrix as input and a specified color palette
# outputs plots
# @param results a data table
# @param chrom_name a string
# @param reg_start a whole number
# @param reg_stop a whole number
# @param wkdir a string
# @param pal a string
# @return plots

# DEPRICATED - Use .plot_heat()
.plot_si_heat <- function(results, chrom_name, reg_start, reg_stop, wkdir, pal = c("RdYlBl", "yelOrRed", "MagYel", "Greens", "BlYel")) {
  if (!missing(pal)) {
    m_pal <- match.arg(pal)
    m_pal <- switch(m_pal,
      "RdYlBl" = c("#4575b4", "#74add1", "#abd9e9", "#e0f3f8", "#fee090", "#fdae61", "#f46d43", "#d73027"),
      "BlYel" = c("#0c2c84", "#225ea8", "#1d91c0", "#41b6c4", "#7fcdbb", "#c7e9b4", "#edf8b1"),
      "yelOrRed" = c("#FFFFCC", "#FFFF66", "#FFCC00", "#FF9900", "#FF6600"),
      "MagYel" = c("#330033", "#660066", "#990066", "#FF6600", "#FFCC00", "#FFFF66"),
      "Greens" = c("#003333", "#006666", "#009966", "#00CC33", "#CCFF99", "#FFFF99", "#FFFF00")
    )

    p <- pheatmap::pheatmap(results, main = "Reads With Proper Overhangs By Size", cluster_cols = FALSE, cluster_rows = FALSE, fontsize = 12, color = m_pal)
  } else {
    p <- pheatmap::pheatmap(results, main = "Reads With Proper Overhangs By Size", cluster_cols = FALSE, cluster_rows = FALSE, fontsize = 12, color = grDevices::colorRampPalette("RdYlBl")(length("RdYlBl")))
  }
  return(p)
}
