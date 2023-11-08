#' function to find hairpin overlaps
#' takes two data tables of reads with start and stop positions
#' converts to gr objects and finds overlaps, transforms positions back to original
#' returns overlap table
#'
#' @param r1_dt a data frame containing chrom name, start, and stop
#' @param r2_dt a data frame containing chrom name, start, and stop
#' @return overlap table



#' @export


find_hp_overlaps <- function(r1_dt, r2_dt){
   gr1_obj <- makeGRobj(r1_dt, "rname", "start", "end")
   gr2_obj <- makeGRobj(r2_dt, "rname", "start", "end")
   r1_stop <- r1_start <- dist <- start_r1 <- widthx <- start_r2 <- widthy <- NULL

   res <- GenomicRanges::findOverlaps(gr1_obj, gr2_obj,
                                      maxgap=-1L, minoverlap=1L,
                                      type=c("any"),
                                      select=c("all"),
                                      ignore.strand=TRUE)

   #queryHits and subjectHits are functions to get the value at an index
   query_idx <- S4Vectors::queryHits(res)
   sub_idx <- S4Vectors::subjectHits(res)
   res <- NULL

   query_values <- data.frame(cbind(r1_start = IRanges::start(gr1_obj[query_idx]), r1_stop = IRanges::end(gr1_obj[query_idx]), r1_width = IRanges::width(gr1_obj[query_idx]))) %>%
      dplyr::mutate(r1_stop = r1_stop - 63, r1_width = r1_stop - r1_start + 1)
   query_idx <- NULL

   sub_values <- data.frame(cbind(r2_start = IRanges::start(gr2_obj[sub_idx]), r2_stop = IRanges::end(gr2_obj[sub_idx]), r2_width = IRanges::width(gr2_obj[sub_idx])))
   sub_idx <- NULL


   res_tbl <- data.frame(cbind(r1_start = query_values$r1_start, r1_width = query_values$r1_width, r1_end = query_values$r1_stop,
                               r2_start = sub_values$r2_start, r2_width = sub_values$r2_width, r2_end = sub_values$r2_stop)) %>%
     dplyr::mutate(dist = r2_start - r1_end) %>%
     dplyr::filter(dist >= 0 & dist < 64) %>%
     dplyr::mutate(r1_end = r1_start + r1_width - 1, r2_end = r2_start + r2_width - 1) %>%
     dplyr::filter(r1_end <= r2_start)
   #need to filter to remove self matches, get dist between end of R1 and start of R2
   #size <- nrow(res_tbl)
   #overlaps <- get_nearby(res_tbl$r1_start, res_tbl$r1_end, res_tbl$r2_start, res_tbl$r2_end, 60, size) %>%
   #dplyr::bind_rows() %>%
   #    dplyr::filter(dist >= 0) %>%
   #    dplyr::mutate(end_r1 = start_r1 + widthx - 1, end_r2 = start_r2 + widthy - 1)
   query_idx <- NULL
   sub_idx <- NULL
   query_values <- NULL
   sub_values <- NULL
   size <- NULL
   return(res_tbl)
}
