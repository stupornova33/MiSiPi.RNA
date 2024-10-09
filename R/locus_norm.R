#' normalize read counts by total number of reads mapped to locus
#' @param dt a dataframe of reads summarized by count
#' @param locus_read_count an integer
#' @return expanded_dt a dataframe

#' @export

locus_norm <- function(dt, locus_read_count) {
  options(scipen = 999)

  #CPM-like
  expanded_dt <- dt %>%
    dplyr::mutate(weighted_count = round(count * (1/locus_read_count) * 10^6))

  expanded_dt <- rep_seq_reads(expanded_dt$weighted_count, expanded_dt$rname,
                               expanded_dt$start, expanded_dt$end,
                               expanded_dt$first, expanded_dt$seq) %>%
    dplyr::mutate(width = end - start + 1)

  return(expanded_dt)
}






