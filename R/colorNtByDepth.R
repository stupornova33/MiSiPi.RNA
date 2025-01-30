# function to specify colors for each nt as a function of depth
#
# seq_cov_dat - Data Frame
# return - Data Frame

.colorNtByDepth <- function(seq_cov_dat) {
  count <- bin <- NULL
  # order table by count for binning
  ordered_tbl <- seq_cov_dat[order(seq_cov_dat$count), ]
  ordered_tbl <- ordered_tbl %>% dplyr::mutate(count = ceiling(count))

  max_count <- max(ordered_tbl$count)
  min_count <- min(ordered_tbl$count)
  print(paste0("min_count: ", min_count, " max_count: ", max_count))

  # divide the expression values into 20 arbitrary bins
  num_bins <- 20
  print(paste0("num_bins: ", num_bins))
  print(paste0(" num_bins: ", num_bins, "  (max_count/num_bins): ", (max_count / num_bins)))

  if ((max_count / num_bins) < 1) {
    print("less than 1")
    window <- 1
  } else {
    print("greater than 1")
    window <- ceiling(max_count / num_bins)
  }

  print(paste0("window: ", window))

  # make a vector from min_count - 1 to max_count + extra to encompass all the values
  vec <- seq((min_count - 1), max_count + window, by = window)

  ordered_tbl[ordered_tbl == 0] <- 0.00001 # ??

  # assign each expression value a bin
  binned_tbl <- ordered_tbl %>% dplyr::mutate(bin = cut(count, breaks = c(vec)))

  # binned_tbl$count <- round(binned_tbl$count)

  # remove those pesky parentheses
  binned_tbl$bin <- gsub("[()]", "", binned_tbl$bin)
  binned_tbl$bin <- gsub("[[]", "", binned_tbl$bin)
  binned_tbl$bin <- gsub("[]]", "", binned_tbl$bin)

  # right now the bin is a range e.g. (-1,200)
  # need to split the values, want to convert them to indexes

  get_bin_num <- function(x) {
    res <- strsplit(binned_tbl$bin, ",")[[x]][2]
    return(res)
  }

  bins <- lapply(1:nrow(binned_tbl), get_bin_num)
  bins <- unlist(bins)
  # remove the original bin from the table
  binned_tbl <- binned_tbl %>% dplyr::select(-c(bin))

  # assign the old bins to new bins (the length of the binned_tbl; so 200 becomes 1, next one becomes 2, etc.)
  new_table <- data.frame(old_bin = c(unique(bins)), new_bin = c(seq(1, length(unique(bins)))))
  new_bins <- vector()
  convert_bins <- function(x) {
    match <- which(new_table$old_bin == bins[[x]])
    match <- new_table$new_bin[[match]]
    new_bins <- append(new_bins, match)
    return(new_bins)
  }

  new_bins <- unlist(lapply(1:length(bins), convert_bins))

  # give the new bins to the old table
  binned_tbl$bin <- new_bins

  # need to put the new values in the same order as the original table
  ordered_color_dat <- binned_tbl[order(binned_tbl$pos), ]
  ordered_color_dat$bin <- format(as.numeric(ordered_color_dat$bin, scientific = F))

  return(ordered_color_dat)
}
