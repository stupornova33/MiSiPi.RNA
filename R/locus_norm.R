#' normalize read counts by total number of reads mapped to locus
#' @param dt a dataframe of reads summarized by count
#' @param locus_read_count an integer
#' @return final_df a dataframe

#' @export


locus_norm <- function(dt, locus_read_count){
  options(scipen = 999)
  dt$rname <- chrom_name
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

  #CPM-like
  counts_dt <- dt %>% dplyr::mutate(weighted_count = round(count * (1/locus_read_count) * 10^6))

  #res <- lapply(seq(nrow(counts_dt)), rep_reads)
  #res_df <- dplyr::bind_rows(res)
  res <- rep_seq_reads(counts_dt$weighted_count, counts_dt$rname, counts_dt$start, counts_dt$end, counts_dt$first, counts_dt$seq)
  #shuffle order randomly
  res_df <- res[sample(1:nrow(res)), ] %>%
    dplyr::mutate(width = end - start + 1)



  if(nrow(res_df) > 0){
    #select a subset
    final_df <- utils::head(res_df, 10000)
  } else {
    final_df <- counts_dt %>% dplyr::mutate(width = end - start + 1)
  }


}






