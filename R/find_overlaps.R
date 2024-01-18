#' function to find regular overlaps
#' takes two data tables of reads with start and stop positions
#' converts to gr objects and finds overlaps
#' returns overlap table
#'
#' @param r1_dt a data frame which contains a read_name, start, and stop
#' @param r2_dt a data frame which contains a read_name, start, and stop
#' @return overlap table



#' @export

find_overlaps <- function(r1_dt, r2_dt = NULL){
   if(!is.null(r2_dt)){
   gr1_obj <- makeGRobj(r1_dt, "rname", "start", "end")
   gr2_obj <- makeGRobj(r2_dt, "rname", "start", "end")

   res <- GenomicRanges::findOverlaps(gr1_obj, gr2_obj,
                                      maxgap=-1L, minoverlap=1L,
                                      type=c("any"),
                                      select=c("all"),
                                      ignore.strand=TRUE)


   } else {
      gr1_obj <- makeGRobj(r1_dt, "rname", "start", "end")
      res <- GenomicRanges::findOverlaps(gr1_obj,maxgap=-1L, minoverlap=1L,
                                         type=c("any"),select=c("all"),
                                         ignore.strand=TRUE, drop.self = TRUE, drop.redundant = TRUE)
      gr2_obj <- gr1_obj
   }
   #queryHits and subjectHits are functions to get the value at an index
   query_idx <- S4Vectors::queryHits(res)
   sub_idx <- S4Vectors::subjectHits(res)
   res <- NULL

   query_values <- data.frame(cbind(r1_start = IRanges::start(gr1_obj[query_idx]), r1_width = IRanges::width(gr1_obj[query_idx]), r1_end = IRanges::end(gr1_obj[query_idx]))) #%>%
      #dplyr::mutate(r1_end = r1_start + r1_width - 1)
   query_idx <- NULL

   sub_values <- data.frame(cbind(r2_start = IRanges::start(gr2_obj[sub_idx]), r2_width = IRanges::width(gr2_obj[sub_idx]), r2_end = IRanges::end(gr2_obj[sub_idx])))#%>%
      #dplyr::mutate(r2_end = r2_start + r2_width - 1)
   sub_idx <- NULL

   #size <- nrow(query_values)
   #reorganize so that 5' read comes first in the table
   #blank_tbl <- data.frame(r1_start = numeric(size), r1_width = numeric(size), r1_end = numeric(size),
    #                     r2_start = numeric(size), r2_width = numeric(size), r2_end = numeric(size))

   #tmp_tbl <- data.frame(cbind(r1_start = query_values$r1_start, r1_width = query_values$r1_width, r1_end = query_values$r1_end,
    #                           r2_start = sub_values$r2_start, r2_width = sub_values$r2_width, r2_end = sub_values$r2_end))
   #res_tbl <- dplyr::bind_rows(blank_tbl, tmp_tbl)

   res_tbl <- data.frame(cbind(r1_start = query_values$r1_start, r1_width = query_values$r1_width, r1_end = query_values$r1_end,
                               r2_start = sub_values$r2_start, r2_width = sub_values$r2_width, r2_end = sub_values$r2_end))

   query_idx <- NULL
   sub_idx <- NULL
   query_values <- NULL
   sub_values <- NULL
   size <- NULL
   return(res_tbl)
}
