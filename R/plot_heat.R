# plot heatmaps for si or pi rna
# takes a matrix as input and a specified color palette
# outputs plots
# @param results a data table
# @param method "piRNA" or "siRNA"
# @param pal a string
# @return plots

.plot_heat <- function(results, method = c("piRNA", "siRNA"), pal = c("RdYlBl", "yelOrRed", "MagYel", "Greens", "BlYel")) {
  method <- match.arg(method)
  m_pal <- match.arg(pal)
  m_pal <- switch(m_pal,
                  "RdYlBl" = c("#4575b4", "#74add1", "#abd9e9", "#e0f3f8", "#fee090", "#fdae61", "#f46d43", "#d73027"),
                  "BlYel" = c("#0c2c84", "#225ea8", "#1d91c0", "#41b6c4", "#7fcdbb", "#c7e9b4", "#edf8b1"),
                  "yelOrRed" = c("#FFFFCC", "#FFFF66", "#FFCC00", "#FF9900", "#FF6600"),
                  "MagYel" = c("#330033", "#660066", "#990066", "#FF6600", "#FFCC00", "#FFFF66"),
                  "Greens" = c("#003333", "#006666", "#009966", "#00CC33", "#CCFF99", "#FFFF99", "#FFFF00")
  )
  
  if (method == "piRNA") {
    plot_title <- "Reads With Proper Overlaps By Size"
  } else {
    plot_title <- "Reads With Proper Overhangs By Size"
  }
  
  p <- pheatmap::pheatmap(results, main = plot_title, cluster_cols = FALSE, cluster_rows = FALSE, fontsize = 12, color = m_pal)
  
  return(p)
}
