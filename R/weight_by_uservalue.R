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

  rep_reads <- function(i) {
    rep_count <- dt$count[i]
    rname <- rep(dt$rname[i], rep_count)
    start <- rep(dt$start[i], rep_count)
    width <- rep(dt$width[i], rep_count)
    end <- rep(dt$end[i], rep_count)
    first <- rep(dt$first[i], rep_count)
    
    # if sequence is present, need to return it as well
    if(!is.null(seq)){
      seq <- rep(dt$seq[i], rep_count)
      current_df <- data.frame(rname = rname, start = start, width = width, end = end, first = first, seq = seq)
    } else {
      current_df <- data.frame(rname = rname, start = start, width = width, end = end, first = first)
    }
    return(current_df)
  }
  
  # (reads * 10^9)/(miRNA_reads * locus_length)
  dt$rep_count <- round((dt$count * 10^9)/(norm * locus_length))
  
  dt <- dt %>% dplyr::arrange(rep_count)
  
  res <- lapply(seq(nrow(dt)), rep_reads)
  res_df <- dplyr::bind_rows(res)
  
  #shuffle order randomly
  res_df <- res_df[sample(1:nrow(res_df)), ] #%>%
  # dplyr::mutate(width = end - start + 1)
  
  
  if(nrow(res_df) > 0){
  #select a random subset
  final_df <- utils::head(res_df, 10000)
  } else {
  final_df <- res_df
  }
  
  return(final_df)
  
  
    
    
  
}