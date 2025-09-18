# function to find hairpin overlaps
# takes two data tables of reads with start and stop positions
# converts to gr objects and finds overlaps, transforms positions back to original
# returns overlap table
#
# @param r1_dt a data frame containing chrom name, start, and stop
# @param r2_dt a data frame containing chrom name, start, and stop
# @param n the number of bases to transform the read by
# @return overlap table

find_overlaps_chunks  <- function(r1_dt, r2_dt, n){
  r1_stop <- r1_start <- dist <- start_r1 <- widthx <- start_r2 <- widthy <- NULL
  
  tot_rows <- nrow(r1_dt)
  
  max_val <- tot_rows
  #n_chunks <- ceiling(log2(max_val))
  #n_chunks <- ceiling(sqrt(max_val))
  #chunk_size <- ceiling(max_val/n_chunks)
  chunk_size <- 2000
  n_chunks <- max_val/chunk_size
  
  gr1_obj <- .makeGRobj(r1_dt, "rname", "start", "end")
  gr2_obj <- .makeGRobj(r2_dt, "rname", "start", "end")
  final_list <- list()
  
  for(i in 1:n_chunks){
    # use size of query obj
    if(i == 1){
      start <- 1
      end <- chunk_size
      old_end <- end
    } else {
      start <- old_end + 1
      end <- start + chunk_size
      old_end <- end
    }
    
    #print(paste0("start: ", start, " end: ", end))
    query_chunk <- gr1_obj[start:end]
    
    hits <- GenomicRanges::findOverlaps(query_chunk, gr2_obj, 
                                        maxgap = -1L, 
                                        minoverlap = 1L,
                                        type = c("any"),
                                        select = c("all"),
                                        ignore.strand = TRUE)
    if (length(hits) > 0) {
      
      # indices inside the chunk
      local_query_idx <- S4Vectors:: queryHits(hits)
      subject_idx <- S4Vectors::subjectHits(hits)
      
      # convert to indices in the full gr1_obj
      # should be n_chunks unique values in this
      global_query_idx <- local_query_idx + (start - 1)
      
      
      overlap_df <- data.frame(
        r1_start = IRanges::start(gr1_obj[global_query_idx]),
        r1_end   = IRanges::end(gr1_obj[global_query_idx]) - n,
        r1_dupes = query_chunk$n[global_query_idx],
        r2_start = IRanges::start(gr2_obj[subject_idx]),
        r2_width = IRanges::width(gr2_obj[subject_idx]),
        r2_end   = IRanges::end(gr2_obj[subject_idx]),
        r2_dupes = gr2_obj$n[subject_idx]) %>%
        
        dplyr::mutate(
          r1_width = r2_end - r2_start + 1,
          dist = r2_start - r1_end) %>%
        dplyr::filter(
          dist >= 0 & dist < 50,
          r1_end <= r2_start)
      
      #final_list[[i]] <- overlap_df
      final_list[[i]] <- overlap_df
    }
    # Combine results
    all_overlaps_df <- dplyr::bind_rows(final_list)
  }
  
}
