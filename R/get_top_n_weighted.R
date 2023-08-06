#' function takes a data table of reads
#' summarizes count of grouped reads
#' returns a weighted top n% of reads
#'
#' @param dt a data table of reads
#' @param chrom_name a string
#' @param n a decimal representing percent of reads to return
#' @return filter_dt
#' @export


get_top_n_weighted <- function(dt, chrom_name, n){
  width <- pos <- start <- end <- first <- count <- NULL
  n <- n/100
  dt$rname <- chrom_name

  rep_reads <- function(i) {
    rep_count <- counts_dt$count[i]
    rname <- rep(counts_dt$rname[i], rep_count)
    start <- rep(counts_dt$start[i], rep_count)
    end <- rep(counts_dt$end[i], rep_count)

    current_df <- data.frame(rname = rname, start = start, end = end)

    return(current_df)
  }



  if(nrow(dt) < 1 ){
    final_df <- data.frame(matrix(ncol = 3, nrow = 0))
    return(final_df)
  } else {
   #testing #############
    #counts_dt <- dt %>%
   #  dplyr::group_by_all() %>% dplyr::summarize(count = dplyr::n())
    counts_dt <- dt
     mn <- min(counts_dt$count)
     mx <- max(counts_dt$count)
     window <- round(mx - mn)/10

     if(window == 0){
        binned_tbl <- counts_dt %>% dplyr::mutate(wt_count = count)
        counts_dt <- counts_dt %>% dplyr::arrange(count)
        final_df <- counts_dt %>% dplyr::mutate(width = end - start + 1)
     } else {
        vec <- seq((mn - 1), mx + window, by = window)
        counts_dt <- counts_dt[counts_dt$count > stats::quantile(counts_dt$count,prob=n),]
        binned_tbl <- counts_dt %>% dplyr::mutate(bin = cut(count, breaks = c(vec)))


        binned_tbl$bin <- gsub("[()]", "", binned_tbl$bin)
        binned_tbl$bin <- gsub("[[]", "", binned_tbl$bin)
        binned_tbl$bin <- gsub("[]]", "", binned_tbl$bin)


        levels <- unique(binned_tbl$bin)
        wt_df <- data.frame(level = levels)

        wt_df$weight <- seq(1, (2- (1/length(levels))), by = 1/length(levels))

        weight <- vector("integer", nrow(counts_dt))

        binned_tbl$weight <- weight

        for(i in 1:nrow(wt_df)){
         idx <- which(binned_tbl$bin == wt_df$level[i])
         binned_tbl$weight[idx] <- wt_df$weight[i]
        }

       binned_tbl <- binned_tbl %>% dplyr::mutate(wt_count = round(weight*count))

       res <- lapply(seq(nrow(binned_tbl)), rep_reads)
       res_df <- dplyr::bind_rows(res)

       #shuffle order randomly
       res_df <- res_df[sample(1:nrow(res_df)), ] %>%
         dplyr::mutate(width = end - start + 1)



      if(nrow(binned_tbl) > 0){
        #select a subset
        final_df <- utils::head(res_df, 10000)
      } else {
        final_df <- counts_dt %>% dplyr::mutate(width = end - start + 1)
      }



    }
  }

  #final_df <- final_df %>% dplyr::select(-c(count))
  return(final_df)
}

