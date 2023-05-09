#' replicates reads from counts dt
#' @param counts_dt a data table or data frame
#' @return res_df

#' @export


rep_reads <- function(counts_dt){
#need to create the reads that represent the counts
  end <- start <- NULL
  inner_func <- function(i) {
    rep_count <- counts_dt$count[i]
    rname <- rep(counts_dt$rname[i], rep_count)
    start <- rep(counts_dt$start[i], rep_count)
    end <- rep(counts_dt$end[i], rep_count)
    
    current_df <- data.frame(rname = rname, start = start, end = end)
    
    return(current_df)
  }
  
  if(nrow(counts_dt) > 0){
    res <- lapply(seq(nrow(counts_dt)), inner_func)
    res_df <- dplyr::bind_rows(res)
    
    #shuffle order randomly
    res_df <- res_df[sample(1:nrow(res_df)), ] %>% 
      dplyr::mutate(width = end - start + 1)
    
    #select a subset 
    final_df <- utils::head(res_df, 10000)
  } else {
    final_df <- counts_dt
  }

return(res_df)
}
