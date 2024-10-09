#' function takes a data table of reads
#' summarizes count of grouped reads
#' returns dataframe of top reads replicated to reflect their proportional weight to the top read
#'
#' @param dt a data table of reads
#' @param chrom_name a string
#' @return expanded_dt
#' @export

weight_by_prop <- function(dt, chrom_name) {
  width <- start <- end <- first <- count <- NULL

  if(nrow(dt) < 1) {
    blank_df <- data.frame(matrix(ncol = 3, nrow = 0))
    return(blank_df)
  } else {
    expanded_dt <- dt
    mn <- min(counts_dt$count)
    mx <- max(counts_dt$count)

    expanded_dt <- expanded_dt %>%
      dplyr::mutate(prop = count/mx)
    expanded_dt <- counts_dt %>%
      dplyr::mutate(weighted_count = prop*count)

    #remove values < 1?
    expanded_dt <- expanded_dt[expanded_dt$weighted_count > 1,]

    if(nrow(expanded_dt) > 0) {
      expanded_dt$weighted_count <- round(expanded_dt$weighted_count)

      expanded_dt <- rep_seq_reads(expanded_dt$weighted_count, expanded_dt$rname,
                                   expanded_dt$start, expanded_dt$end,
                                   expanded_dt$first, expanded_dt$seq) %>%
        dplyr::mutate(width = end - start + 1)
    }
  }
  return(expanded_dt)
}
