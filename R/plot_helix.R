#' plots the arc diagram
#' @param filePath a string
#' @return an arc plot
#' @export

plot_helix <- function(filePath){
   grDevices::dev.control("enable")
   R4RNA::plotHelix(helix = R4RNA::readHelix(filePath), line = TRUE, arrow = FALSE, lwd = 2.25, scale = FALSE)
   temp <- grDevices::recordPlot()
   return(temp)
}
