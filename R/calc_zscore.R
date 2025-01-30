# Calculate zscore
#
# allCounts - Data Frame
# return z_score

.calc_zscore <- function(allCounts) {
  # function to calculate z-score
  sd <- NULL
  z_score <- (allCounts - mean(allCounts)) / sd(allCounts)
  return(z_score)
}
