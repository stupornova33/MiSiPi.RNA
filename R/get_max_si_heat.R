# function to find max_si heatmap cell values
# takes a matrix created from siRNA function
# finds max row, column, and count of highest cell
# returns highest_si_row, highest_si_col
#
# @param si_res a data table
# @return max_si_df

.get_max_si_heat <- function(si_res) {
  if (sum(si_res[[1]]) > 0) {

    # This offset will be used to make the matrix indicies human readable
    # Matrix of size 15 with indicies of 1-15 when offset would indicate
    # Read widths of 18-32
    # Minimum read width is 18, so 18 - 1 is the offset
    MINIMUM_READ_WIDTH <- 18
    MATRIX_INDEX_OFFSET <- MINIMUM_READ_WIDTH - 1

    max_results <- as.matrix(unlist(unname(which(si_res[[1]] == max(si_res[[1]]), arr.ind = TRUE))))
    si_table <- as.matrix(si_res[[1]])
    if (nrow(max_results) == 0) {
      highest_si_count <- 0
      highest_si_row <- NA
      highest_si_col <- NA
    } else if (nrow(max_results) == 1) {
      highest_si_row <- max_results[1, 1]
      highest_si_col <- max_results[1, 2]
      highest_si_count <- si_table[highest_si_row, highest_si_col]
    } else {
      counts <- si_table[max_results]
      highest_si_count <- mean(counts)
      highest_si_row <- mean(max_results[, 1])
      highest_si_col <- mean(max_results[, 2])
    }
    highest_si_row <- highest_si_row + MATRIX_INDEX_OFFSET
    highest_si_col <- highest_si_col + MATRIX_INDEX_OFFSET
  } else {
    highest_si_row <- NA
    highest_si_col <- NA
  }

  if (is.na(highest_si_row) && is.na(highest_si_col)) {
    highest_si_count <- highest_si_row <- highest_si_col <- -33
  }
  max_si_df <- data.frame(highest_si_row = highest_si_row, highest_si_col = highest_si_col, highest_si_count = highest_si_count)
  return(max_si_df)
}
