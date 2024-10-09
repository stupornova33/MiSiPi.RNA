#' Weight_by_uservalue
#' Takes a "collapsed" data.frame of reads summarized by count.
#' Replicates them according to a factor provided by the user.
#'
#' @param dt A dataframe of reads.
#' @param norm An integer chosen by the user which will be used to normalize the read counts.
#' @param locus_length An integer corresponding to the length of the region of interest
#' @return expanded_dt A dataframe of reads which have been "uncollapsed" and replicated by the normalization factor.
#' @export

weight_by_uservalue <- function(dt, norm, locus_length) {
  expanded_dt <- dt
  expanded_dt$rep_count <- round((dt$count * 10^9)/(norm * locus_length))

  expanded_dt <- expanded_dt %>%
    dplyr::arrange(rep_count)

  expanded_dt <- rep_seq_reads(expanded_dt$rep_count, expanded_dt$rname,
                               expanded_dt$start, expanded_dt$end,
                               expanded_dt$first, expanded_dt$seq) %>%
    dplyr::mutate(width = end - start + 1)

  return(expanded_dt)
}
