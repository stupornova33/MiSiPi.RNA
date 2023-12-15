#' function takes a data table of reads
#' summarizes count of grouped reads
#' returns dataframe of top reads replicated to reflect their weight
#'
#' @param dt a data table of reads
#' @param chrom_name a string
#' @param seq a string, "T" or "F"
#' @return filter_dt
#' @export


get_top_n_weighted <- function(dt, chrom_name, seq = NULL){
  width <- pos <- start <- end <- first <- count <- NULL
  dt$rname <- chrom_name

  rep_reads <- function(i) {
    rep_count <- counts_dt$weighted_count[i]
    rname <- rep(counts_dt$rname[i], rep_count)
    start <- rep(counts_dt$start[i], rep_count)
    end <- rep(counts_dt$end[i], rep_count)
    first <- rep(counts_dt$first[i], rep_count)

    #for getting siRNA pairs
    if(!is.null(seq)){
      seq <- rep(counts_dt$seq[i], rep_count)
      current_df <- data.frame(rname = rname, start = start, end = end, first = first, seq = seq)
    } else {
      current_df <- data.frame(rname = rname, start = start, end = end, first = first)
    }

    return(current_df)
  }


  if(nrow(dt) < 1 ){
    final_df <- data.frame(matrix(ncol = 3, nrow = 0))
    return(final_df)
  } else {
    counts_dt <- dt
     mn <- min(counts_dt$count)
     mx <- max(counts_dt$count)

     counts_dt <- counts_dt %>% dplyr::mutate(prop = count/mx)
     counts_dt <- counts_dt %>% dplyr::mutate(weighted_count = prop*count)

     #counts cant be fractions, set counts < 1 to 1
     #counts_dt$weighted_count[counts_dt$weighted_count < 1] <- 1

     #remove values < 1?
     counts_dt <- counts_dt[counts_dt$weighted_count > 1,]


     if(nrow(counts_dt) > 0){
     counts_dt$weighted_count <- round(counts_dt$weighted_count)

     res <- lapply(seq(nrow(counts_dt)), rep_reads)
     res_df <- dplyr::bind_rows(res)


     # speed things up by selecting random subset.
     # find overlaps struggles with larger vectors

     #shuffle order randomly
     res_df <- res_df[sample(1:nrow(res_df)), ] %>%
       dplyr::mutate(width = end - start + 1)

     if(nrow(res_df) > 0){
       #select a subset
       final_df <- utils::head(res_df, 10000)
     } else {
       final_df <- counts_dt %>% dplyr::mutate(width = end - start + 1)
     }

     } else {
       final_df <- data.frame(matrix(ncol = 3, nrow = 0))
  }

 }

  return(final_df)
}

