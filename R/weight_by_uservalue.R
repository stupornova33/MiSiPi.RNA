#' Weight_by_uservalue
#' Takes a "collapsed" data.frame of reads summarized by count.
#' Replicates them according to a factor provided by the user.
#'
#' @param dt A dataframe of reads.
#' @param norm An integer chosen by the user which will be used to normalize the read counts.
#' @param locus_length An integer corresponding to the length of the region of interest
#' @return final_df A dataframe of reads which have been "uncollapsed" and replicated by the normalization factor.
#' @export

weight_by_uservalue <- function(dt, norm, locus_length){
  dt$rname <- chrom_name
  dt$rep_count <- round((dt$count * 10^9)/(norm * locus_length))

  dt <- dt %>% dplyr::arrange(rep_count)

  res_df <- rep_seq_reads(counts_dt$rep_count, counts_dt$rname, counts_dt$start, counts_dt$end, counts_dt$first, counts_dt$seq)

  #shuffle order randomly
  #res_df <- res_df[sample(1:nrow(res_df)), ] #%>%
  # dplyr::mutate(width = end - start + 1)


  #if(nrow(res_df) > 0){
  #select a random subset
  #final_df <- utils::head(res_df, 10000)
  #} else {
  final_df <- res_df
  #}

  return(final_df)





}
