#' function takes a data table of reads
#' summarizes count of grouped reads
#' returns an unweighted top n% of reads
#'
#' @param filter_dt a data table or data frame of reads
#' @param chrom_name a string
#' @return expanded_dt a dataframe
#' @export

### this is the unweighted version.
no_weight <- function(filter_dt, chrom_name){
  width <- pos <- start <- end <- first <- NULL
  filter_dt$rname <- chrom_name

  expanded_dt <- filter_dt %>%
    dplyr::arrange(count)

  expanded_dt <- rep_seq_reads(expanded_dt$count, expanded_dt$rname,
                               expanded_dt$start, expanded_dt$end,
                               expanded_dt$first, expanded_dt$seq) %>%
    dplyr::mutate(width = end - start + 1)

  return(expanded_dt)
}
