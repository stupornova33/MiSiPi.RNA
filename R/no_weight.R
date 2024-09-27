#' function takes a data table of reads
#' summarizes count of grouped reads
#' returns an unweighted top n% of reads
#'
#' @param filter_dt a data table or data frame of reads
#' @param chrom_name a string
#' @return filter_dt
#' @export


### this is the unweighted version.
no_weight <- function(filter_dt, chrom_name){
   width <- pos <- start <- end <- first <- NULL
   filter_dt$rname <- chrom_name

   counts_dt <- filter_dt %>% dplyr::arrange(count)

   res <- rep_seq_reads(counts_dt$count, counts_dt$rname, counts_dt$start, counts_dt$end, counts_dt$first, counts_dt$seq)

   res_df <- res[sample(1:nrow(res)), ] %>%
     dplyr::mutate(width = end - start + 1)

   res_df <- na.omit(res_df)

   final_df <- res_df
   return(final_df)
}
