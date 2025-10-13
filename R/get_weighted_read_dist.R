# function to calculate the read distribution
# takes forward and reverse data frames
# summarizes read size by count
# returns table of size distribution
#
# @param forward_dt a minus strand object
# @param reverse_dt a plus strand object

.get_weighted_read_dist <- function(forward_dt, reverse_dt) {
  options(scipen = 999)
  width <- pos <- count.y <- count.x <- NULL
  all_dat <- rbind(forward_dt, reverse_dt)

  sizeDist <- all_dat %>%
    dplyr::group_by(width, first) %>%
    dplyr::summarize(count = dplyr::n())

  all_dat <- NULL

  widths <- c(
    rep(18, times = 4), rep(19, times = 4), rep(20, times = 4), rep(21, times = 4),
    rep(22, times = 4), rep(23, times = 4), rep(24, times = 4), rep(25, times = 4),
    rep(26, times = 4), rep(27, times = 4), rep(28, times = 4), rep(29, times = 4),
    rep(30, times = 4)
  )
  first <- c("A", "G", "C", "T")
  empty_dt <- data.frame(width = widths, first = rep(first, times = 13), count = 0)

  size_dist <- merge(sizeDist, empty_dt, by = c("width", "first"), all = TRUE) %>%
    dplyr::select(-c(count.y)) %>%
    dplyr::rename("count" = count.x)
  size_dist["count"][is.na(size_dist["count"])] <- 0

  sizeDist <- NULL

  return(size_dist)
}
