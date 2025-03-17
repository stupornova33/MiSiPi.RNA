# plots the arc diagram
# @param filePath a string
# @return an arc plot

.plot_helix <- function(filePath) {
  R4RNA::plotHelix(helix = R4RNA::readHelix(filePath), line = TRUE, arrow = FALSE, lwd = 2.25, scale = FALSE)
  grDevices::dev.control("enable")
  temp <- grDevices::recordPlot()
  return(temp)
}
