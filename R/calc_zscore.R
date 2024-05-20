#' calc z_score
#' takes a data frame
#' returns z_score
#' 
#' @param allCounts a data frame
#' @return z_score

#' @export

calc_zscore <- function(allCounts){
  #function to calculate z_score
  sd <- NULL
  z_score <- (allCounts - mean(allCounts))/sd(allCounts)
  return(z_score)
}
