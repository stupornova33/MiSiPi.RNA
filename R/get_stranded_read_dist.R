# .get_stranded_read_dist
# takes forward and reverse data frames
# summarizes read size by count
# returns table of size distribution
# @param plus_df - A data frame of plus strand reads
# @param minus_df - A data frame of minus strand reads
# @return size_dist

.get_stranded_read_dist <- function(plus_df, minus_df) {
  options(scipen = 999)
  
  plus_dist <- plus_df %>%
    dplyr::group_by(width, first) %>%
    dplyr::summarize(count = dplyr::n()) %>%
    dplyr::mutate(strand = "plus")
  
  minus_dist <- minus_df %>%
    dplyr::group_by(width, first) %>%
    dplyr::summarize(count = dplyr::n()) %>%
    dplyr::mutate(strand = "minus") %>%
    dplyr::mutate(count = (count * -1))
  
  widths <- c(
    rep(18, times = 4), rep(19, times = 4), rep(20, times = 4), rep(21, times = 4),
    rep(22, times = 4), rep(23, times = 4), rep(24, times = 4), rep(25, times = 4),
    rep(26, times = 4), rep(27, times = 4), rep(28, times = 4), rep(29, times = 4),
    rep(30, times = 4), rep(31, times = 4), rep(32, times = 4)
  )
  first <- c("A", "G", "C", "T")
  empty_df <- data.frame(width = widths, first = rep(first, times = 15), count = 0)
  
  p_size_dist <- merge(plus_dist, empty_df, by = c("width", "first"), all = TRUE) %>%
    dplyr::select(-c(count.y)) %>%
    dplyr::rename("count" = count.x)
  p_size_dist["count"][is.na(p_size_dist["count"])] <- 0
  p_size_dist["strand"][is.na(p_size_dist["strand"])] <- "plus"
  plus_dist <- NULL
  
  m_size_dist <- merge(minus_dist, empty_df, by = c("width", "first"), all = TRUE) %>%
    dplyr::select(-c(count.y)) %>%
    dplyr::rename("count" = count.x)
  m_size_dist["count"][is.na(m_size_dist["count"])] <- 0
  m_size_dist["strand"][is.na(m_size_dist["strand"])] <- "minus"
  minus_dist <- NULL
  
  all_sizes <- rbind(p_size_dist, m_size_dist)
  
  p_size_dist <- NULL
  m_size_dist <- NULL
  
  return(all_sizes)
}
