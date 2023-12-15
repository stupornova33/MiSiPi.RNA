#' function takes a data table of reads
#' summarizes count of grouped reads
#' returns an unweighted top n% of reads
#'
#' @param chrom a chrom object
#' @param chrom_name a string
#' @param seq a string, "T" or "F"
#' @return filter_dt
#' @export


### this is the unweighted version.
# TODO change the name
get_top_n <- function(filter_dt, chrom_name, seq = NULL){
   width <- pos <- start <- end <- first <- count <- NULL

   rep_reads <- function(i) {
      rep_count <- filter_dt$count[i]
      rname <- rep(filter_dt$rname[i], rep_count)
      start <- rep(filter_dt$start[i], rep_count)
      end <- rep(filter_dt$end[i], rep_count)
      first <- rep(filter_dt$first[i], rep_count)

      #for getting siRNA pairs
      if(!is.null(seq)){
        seq <- rep(counts_dt$seq[i], rep_count)
        current_df <- data.frame(rname = rname, start = start, end = end, first = first, seq = seq)
      } else {
        current_df <- data.frame(rname = rname, start = start, end = end, first = first)
      }
      return(current_df)
   }

   counts_dt <- filter_dt %>% dplyr::arrange(count)

   res <- lapply(seq(nrow(filter_dt)), rep_reads)
   res_df <- dplyr::bind_rows(res)

   #shuffle order randomly
   res_df <- res_df[sample(1:nrow(res_df)), ] %>%
     dplyr::mutate(width = end - start + 1)


   if(nrow(counts_dt) > 0){
      #select a random subset
      final_df <- utils::head(res_df, 10000)
   } else {
      final_df <- counts_dt
   }

   return(final_df)
}
