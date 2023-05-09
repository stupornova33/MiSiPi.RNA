#' calc zscore
#' takes a data frame
#' returns zscore
#' 
#' @param allCounts
#' @return z_score

#' @export

calc_zscore <- function(allCounts){
  #function to calculate z-score
  sd <- NULL
  z_score <- (allCounts - mean(allCounts))/sd(allCounts)
  return(z_score)
}
