# calc_zscore and calc_ml_zscore
#
# Description: Calculates Z Scores
# Valid Z Scores are calculated the same way in both functions
# They differ only in how they treat scenarios where they divide by 0
#
# calc_zscore() assigns 0 as the Z Score if the standard devition is 0 to prevent dividing by 0
# calc_ml_zscore() assigns -33 as the Z Score in that same situation
#
# Intended use:
# calc_zscore() returns zscores used in plotting
# calc_ml_zscore() return zscores that harshly weight divide by zero situations
#    for use in downstream machine learning
#
# Parameters:
# [counts] - A vector of integers usually passed in from as a column from a data frame
#
# Return: A vector of Z Scores

.calc_zscore <- function(counts) {
  sd <- NULL
  
  if (sd(counts) == 0) {
    return(rep(0, length(counts)))
  }
  return((counts - mean(counts)) / sd(counts))
}

.calc_ml_zscore <- function(counts) {
  sd <- NULL
  
  if (sd(counts) == 0) {
    return(rep(-33, length(counts)))
  }
  return((counts - mean(counts)) / sd(counts))
}
