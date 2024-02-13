#' function takes a data table of reads
#' summarizes count of grouped reads
#' returns an unweighted top n% of reads
#'
#' @param chrom a chrom object
#' @param chrom_name a string
#' @return filter_dt
#' @export


### this is the unweighted version.
no_weight <- function(filter_dt, chrom_name){
   width <- pos <- start <- end <- first <- count <- NULL

   #rep_reads <- function(i) {
   #  rep_count <- filter_dt$count[i]
   #  rname <- rep(filter_dt$rname[i], rep_count)
   #  start <- rep(filter_dt$start[i], rep_count)
   #  width <- rep(filter_dt$width[i], rep_count)
   #  end <- rep(filter_dt$end[i], rep_count)
   #  first <- rep(filter_dt$first[i], rep_count)

      #if sequence is present, need to return it as well
   #    if(!is.null(seq)){
   #    seq <- rep(counts_dt$seq[i], rep_count)
   #    current_df <- data.frame(rname = rname, start = start, width = width, end = end, first = first, seq = seq)
   #  } else {
   #    current_df <- data.frame(rname = rname, start = start, width = width, end = end, first = first)
   #  }
   #  return(current_df)
   #}

   counts_dt <- filter_dt %>% dplyr::arrange(count)

   #if(!is.null(seq)){
  #   res <- rep_seq_reads(counts_dt$count, counts_dt$rname, counts_dt$start, counts_dt$end, counts_dt$first, counts_dt$seq)
  # } else {
  #   res <- rep_nonseq_reads(counts_dt$count, counts_dt$rname, counts_dt$start, counts_dt$end, counts_dt$first)
  # }

   #res <- lapply(seq(nrow(filter_dt)), rep_reads)
   #res_df <- dplyr::bind_rows(res)
   res_df <- rep_seq_reads(counts_dt$count, counts_dt$rname, counts_dt$start, counts_dt$end, counts_dt$first, counts_dt$seq)
   #shuffle order randomly
   res_df <- res_df[sample(1:nrow(res_df)), ] #%>%
    # dplyr::mutate(width = end - start + 1)


   if(nrow(counts_dt) > 0){
      #select a random subset
      final_df <- utils::head(res_df, 10000)
   } else {
      final_df <- counts_dt
   }

   return(final_df)
}
