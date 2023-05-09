#' function to filter data table for short hairpin function
#' makes chrom obj into bam df and filters reads between 18-25 nt
#' summarizes count of grouped reads
#' returns top 10% of reads
#' 
#' @param chrom a chrom object
#' @param chrom_name a string
#' @param n a decimal 
#' @return filter_dt

#' @export


get_top_n <- function(chrom, chrom_name, n){
   width <- pos <- start <- end <- first <- count <- NULL
   n <- n/100
   filter_dt <- data.table::setDT(makeBamDF(chrom)) %>%
      base::subset(width <= 25 & width >= 18) %>% 
      dplyr::mutate(start = pos, end = pos + width - 1) %>%
      dplyr::group_by_at(dplyr::vars(start, end)) %>% 
      dplyr::select(-c(pos, seq, first))
   filter_dt$rname <- chrom_name
   counts_dt <- filter_dt %>% 
      dplyr::group_by_all() %>% dplyr::summarise(count = dplyr::n()) %>%
      dplyr::filter(count > 2)
   
   counts_dt <- counts_dt %>% dplyr::arrange(count) 
  
   counts_dt <- counts_dt[counts_dt$count > stats::quantile(counts_dt$count,prob=n),]
   
   #need to create the reads that represent the counts
   rep_reads <- function(i) {
      rep_count <- counts_dt$count[i]
      rname <- rep(counts_dt$rname[i], rep_count)
      start <- rep(counts_dt$start[i], rep_count)
      end <- rep(counts_dt$end[i], rep_count)
      
      current_df <- data.frame(rname = rname, start = start, end = end)
      
      return(current_df)
   }
   
   if(nrow(counts_dt) > 0){
   res <- lapply(seq(nrow(counts_dt)), rep_reads)
   res_df <- dplyr::bind_rows(res)
   
   #shuffle order randomly
   res_df <- res_df[sample(1:nrow(res_df)), ] %>% 
      dplyr::mutate(width = end - start + 1)
   
   #select a subset 
   final_df <- utils::head(res_df, 10000)
   } else {
      final_df <- counts_dt
   }
   
   return(final_df)
}
