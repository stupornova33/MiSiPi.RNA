# function to find max_pi heatmap cell values
# takes a matrix created from piRNA function
# finds max row, column, and count of highest cell
# returns highest_pi_row, highest_pi_col
#
# @param pi_res a data table
# @return max_pi_df

.get_max_pi_heat <- function(pi_res) {
  pi_table <- as.matrix(pi_res[[1]])

  max_results <- as.matrix(unlist(unname(which(pi_res[[1]] == max(pi_res[[1]]), arr.ind = TRUE))))

  # This offset will be used to make the matrix indicies human readable
    # Matrix of size 15 with indicies of 1-15 when offset would indicate
    # Read widths of 18-32
    # Minimum read width is 18, so 18 - 1 is the offset
    MINIMUM_READ_WIDTH <- 18
    MATRIX_INDEX_OFFSET <- MINIMUM_READ_WIDTH - 1

  if (nrow(max_results) == 0) {
    highest_pi_count <- 0
    highest_pi_row <- NA
    highest_pi_col <- NA
  } else if (nrow(max_results) == 1) {
    highest_pi_row <- max_results[1, 1]
    highest_pi_col <- max_results[1, 2]
    highest_pi_count <- pi_table[highest_pi_row, highest_pi_col]
  } else {
    counts <- pi_table[max_results]
    highest_pi_count <- mean(counts)
    highest_pi_row <- mean(max_results[, 1])
    highest_pi_col <- mean(max_results[, 2])
  }

  highest_pi_row <- highest_pi_row + MATRIX_INDEX_OFFSET
  highest_pi_col <- highest_pi_col + MATRIX_INDEX_OFFSET

  if (is.na(highest_pi_row) && is.na(highest_pi_col)) {
    highest_pi_count <- -33
    highest_pi_row <- -33
    highest_pi_col <- -33
  }
  max_pi_df <- data.frame(highest_pi_row = highest_pi_row, highest_pi_col = highest_pi_col, highest_pi_count = highest_pi_count)
  return(max_pi_df)
}
