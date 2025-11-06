# function to find max_si heatmap cell values
# takes a matrix created from siRNA function
# finds max row, column, and count of highest cell
# returns highest_si_row, highest_si_col
#
# @param si_res a data table
# @return max_si_df

.get_max_si_heat <- function(heat_matrix) {
  # If heat matrix is null, na, or all zeros, set default values
  if (!is.numeric(heat_matrix) | sum(heat_matrix) <= 0) {
    highest_si_row <- -33
    highest_si_col <- -33
    highest_si_count <- -33
  } else {
    # This offset will be used to make the matrix indexes human readable
    # Matrix of size 15 with indexes of 1-15 when offset would indicate
    # Read widths of 18-32
    # Minimum read width is 18, so 18 - 1 is the offset
    MINIMUM_READ_WIDTH <- 18
    MATRIX_INDEX_OFFSET <- MINIMUM_READ_WIDTH - 1

    max_results <- unname(which(heat_matrix == max(heat_matrix), arr.ind = TRUE))

    if (nrow(max_results) == 0) {
      highest_si_count <- -33
      highest_si_row <- -33
      highest_si_col <- -33
    } else if (nrow(max_results) == 1) {
      highest_si_row <- max_results[1, 1]
      highest_si_col <- max_results[1, 2]
      highest_si_count <- heat_matrix[highest_si_row, highest_si_col]
    } else {
      counts <- heat_matrix[max_results]
      highest_si_count <- mean(counts)
      highest_si_row <- mean(max_results[, 1])
      highest_si_col <- mean(max_results[, 2])
    }
    
    highest_si_row <- highest_si_row + MATRIX_INDEX_OFFSET
    highest_si_col <- highest_si_col + MATRIX_INDEX_OFFSET
  }

  max_si_df <- data.frame(
    highest_si_row = highest_si_row,
    highest_si_col = highest_si_col,
    highest_si_count = highest_si_count)
  
  return(max_si_df)
}
