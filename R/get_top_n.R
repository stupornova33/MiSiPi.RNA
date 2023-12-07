#' function takes a data table of reads
#' summarizes count of grouped reads
#' returns an unweighted top n% of reads
#'
#' @param chrom a chrom object
#' @param chrom_name a string
#' @param n a decimal
#' @return filter_dt
#' @export

get_top_n <- function(filter_dt, chrom_name, n){
   width <- pos <- start <- end <- first <- count <- NULL
   #n <- n/100
   rep_reads <- function(i) {
      rep_count <- filter_dt$count[i]
      rname <- rep(filter_dt$rname[i], rep_count)
      start <- rep(filter_dt$start[i], rep_count)
      end <- rep(filter_dt$end[i], rep_count)
      first <- rep(filter_dt$first[i], rep_count)
      current_df <- data.frame(rname = rname, start = start, end = end, first= first)

      return(current_df)
   }

   counts_dt <- filter_dt %>% dplyr::arrange(count)


   #counts_dt <- counts_dt[counts_dt$count > stats::quantile(counts_dt$count,prob=n),]

   res <- lapply(seq(nrow(filter_dt)), rep_reads)
   res_df <- dplyr::bind_rows(res)

   #res_df$rname <- chrom_name
   #shuffle order randomly
   res_df <- res_df[sample(1:nrow(res_df)), ] %>%
     dplyr::mutate(width = end - start + 1)


   #counts_dt <- res_df %>%
  #    dplyr::group_by_all() %>% dplyr::summarize(count = dplyr::n())

###################################################################################################################################
   #need to create the reads that represent the counts

   if(nrow(counts_dt) > 0){
      #select a subset
      final_df <- utils::head(res_df, 10000)
   } else {
      final_df <- counts_dt
   }

   return(final_df)
}
