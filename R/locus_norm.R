#' normalize read counts by total number of reads mapped to locus
#' @param dt a dataframe of reads summarized by count
#' @param locus_read_count an integer
#' @return final_df a dataframe

#' @export


locus_norm <- function(dt, locus_read_count){
  options(scipen = 999)
  
  rep_reads <- function(i) {
    rep_count <- counts_dt$weighted_count[i]
    rname <- rep(counts_dt$rname[i], rep_count)
    start <- rep(counts_dt$start[i], rep_count)
    end <- rep(counts_dt$end[i], rep_count)
    first <- rep(counts_dt$first[i], rep_count)
    
    current_df <- data.frame(rname = rname, start = start, end = end, first = first)
    
    return(current_df)
  }
  
  counts_dt <- dt %>% dplyr::mutate(weighted_count = round((count/locus_read_count)*1000000)) 

  res <- lapply(seq(nrow(counts_dt)), rep_reads)
  res_df <- dplyr::bind_rows(res)
  
  #shuffle order randomly
  res_df <- res_df[sample(1:nrow(res_df)), ] %>%
    dplyr::mutate(width = end - start + 1)
  
  
  
  if(nrow(res_df) > 0){
    #select a subset
    final_df <- utils::head(res_df, 10000)
  } else {
    final_df <- counts_dt %>% dplyr::mutate(width = end - start + 1)
  }
  
  
}





