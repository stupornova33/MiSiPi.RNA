#' normalize read counts by total number of reads mapped to locus
#' @param dt a dataframe of reads summarized by count
#' @param locus_read_count an integer
#' @param seq a string, "T" or "F"
#' @return final_df a dataframe

#' @export


locus_norm <- function(dt, locus_read_count, seq = NULL){
  options(scipen = 999)

  #rep_reads <- function(i) {
  #  rep_count <- counts_dt$weighted_count[i]
  #  rname <- rep(counts_dt$rname[i], rep_count)
  #  start <- rep(counts_dt$start[i], rep_count)
  #  end <- rep(counts_dt$end[i], rep_count)
  #  first <- rep(counts_dt$first[i], rep_count)

    #for getting siRNA pairs
  #  if(!is.null(seq)){
  #    seq <- rep(counts_dt$seq[i], rep_count)
  #    current_df <- data.frame(rname = rname, start = start, end = end, first = first, seq = seq)
  #  } else {
  #    current_df <- data.frame(rname = rname, start = start, end = end, first = first)
  #  }

  #  return(current_df)
  #}

  #do an RPKM-like normalization
  #counts_dt <- dt %>% dplyr::mutate(weighted_count = count/(locus_length/1000*locus_read_count/100000))

  #CPM-like (thousand)
  counts_dt <- dt %>% dplyr::mutate(weighted_count = round(count * 1/locus_read_count * 10^3))

  #res <- lapply(seq(nrow(counts_dt)), rep_reads)
  #res_df <- dplyr::bind_rows(res)
  res <- rep_seq_reads(dt$weighted_count, dt$rname, dt$start, dt$end, dt$first, dt$seq)
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






