#' function takes a data table of reads
#' summarizes count of grouped reads
#' returns dataframe of top reads replicated to reflect their proportional weight to the top read
#'
#' @param dt a data table of reads
#' @param chrom_name a string
#' @return filter_dt
#' @export


weight_by_prop <- function(dt, chrom_name){
  width <- pos <- start <- end <- first <- count <- NULL
  dt$rname <- chrom_name

  if(nrow(dt) < 1 ){
    final_df <- data.frame(matrix(ncol = 3, nrow = 0))
    return(final_df)
  } else {
    counts_dt <- dt
     mn <- min(counts_dt$count)
     mx <- max(counts_dt$count)

     counts_dt <- counts_dt %>% dplyr::mutate(prop = count/mx)
     counts_dt <- counts_dt %>% dplyr::mutate(weighted_count = prop*count)

     #remove values < 1?
     counts_dt <- counts_dt[counts_dt$weighted_count > 1,]


     if(nrow(counts_dt) > 0){
     counts_dt$weighted_count <- round(counts_dt$weighted_count)


     final_df <- rep_seq_reads(counts_dt$weighted_count, counts_dt$rname, counts_dt$start, counts_dt$end, counts_dt$first, counts_dt$seq)

     # speed things up by selecting random subset.
     # find overlaps struggles with larger vectors

     #shuffle order randomly
     #res_df <- res_df[sample(1:nrow(res_df)), ] %>%
     #   dplyr::mutate(width = end - start + 1)

     #if(nrow(res_df) > 0){
       #select a subset
       #final_df <- utils::head(res_df, 10000)
     #} else {
     #    final_df <- counts_dt %>% dplyr::mutate(width = end - start + 1)
     #}

     #} else {
      # final_df <- data.frame(matrix(ncol = 3, nrow = 0))
  }

 }

  return(final_df)
}

