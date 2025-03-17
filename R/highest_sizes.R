# function to find the read size with the highest abundance
# takes a table of read sizes and counts
# returns highest and second highest read sizes
#
# @param read_dist a table
# @return a list

.highest_sizes <- function(read_dist) {
  max <- which(read_dist$count == max(read_dist$count))
  if (length(max) > 1) {
    # take the average of the first two results
    highest_size <- mean(read_dist$width[max[1]], read_dist$width[max[2]])
  } else {
    highest_size <- read_dist$width[max]
  }

  return(highest_size)
}
