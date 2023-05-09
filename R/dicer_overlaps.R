#' calculates overlaps for Dicer signature 
#' @param filter_dt a data frame
#' @param helix a data frame
#' @param chrom_name a string
#' @param reg_start an integer
#' @return all_overlaps

#' @export

dicer_overlaps <- function(dicer_dt, helix_df, chrom_name, reg_start){
   # convert helix positions to chromosome positions
   bin <- j <- X.End <- X.Start <- Y.Start <- Y.End <- paired_pos <- start <- width <- NULL
  
   helix_df <- helix_df %>% dplyr::mutate(i = i + reg_start, j = j + reg_start)
   ### paired bases are several nucleotides diff. than read starts... pad the ends of the helix
   
   i_begin <- helix_df$i[1]
   i_end <- helix_df$i[nrow(helix_df)]
   j_begin <- helix_df$j[1]
   j_end <- helix_df$j[nrow(helix_df)]
   
   istartvals <- seq(i_begin - 8, i_begin - 1)
   iendvals <- seq(i_end + 1, i_end + 8)
   jstartvals <- seq(j_begin + 8,j_begin + 1 )
   jendvals <- seq(j_end - 1,j_end - 8)
   
   startdf <- data.frame(i = istartvals, j = jstartvals, length = c(rep(1, times = 4)), value = c(rep(NA, times = 4)))
   enddf <- data.frame(i = iendvals, j = jendvals, length = c(rep(1, times = 4)), value = c(rep(NA, times = 4)))
   
   tmp_df <- rbind(startdf, helix_df)
   final_helix_df <- rbind(tmp_df, enddf)
   final_helix_df <- final_helix_df %>% dplyr::filter(j - i > 15)
   
   if(nrow(final_helix_df) == 0){
      i_j_overlaps <- data.frame(r1_start = 0, r1_width = 0, r1_end = 0,
                                 r2_start = 0, r2_width = 0, r2_end = 0)
      return(i_j_overlaps)
   }
   
   grouped_helix <- group_helix_res(final_helix_df$i, final_helix_df$j)
   filter_helix <- grouped_helix %>% dplyr::filter(X.End - X.Start > 17 | Y.Start - Y.End > 17)
  
   i_dat <- data.frame(rname = numeric(0), start = numeric(0), end = numeric(0), width = numeric(0))
   j_dat <- data.frame(rname = numeric(0), start = numeric(0), end = numeric(0), width = numeric(0))
   
   if(nrow(filter_helix) > 0){
      for(i in 1:nrow(filter_helix)){
         x_rng <- seq(filter_helix$X.Start[i], filter_helix$X.End[i])
         y_rng <- seq(filter_helix$Y.End[i], filter_helix$Y.Start[i])
        
         i_idx <- which(dicer_dt$start %in% x_rng | dicer_dt$end %in% x_rng)
         j_idx <- which(dicer_dt$start %in% y_rng | dicer_dt$end %in% y_rng)
         
         i_dat <- rbind(i_dat, dicer_dt[i_idx,])
         j_dat <- rbind(j_dat, dicer_dt[j_idx,])
      }
       i_idx <- match(i_dat$start, i_dat$start)
       j_idx <- match(j_dat$start, j_dat$start)
       
      #get reads that start or end at a paired position
      i_dat <- i_dat %>% dplyr::mutate(paired_pos = final_helix_df$j[i_idx]) %>%
         dplyr::mutate(start = paired_pos, end = start + width, rname = chrom_name) %>%
         #dplyr::select(-c(r1_start, r1_width, r1_end, paired_pos)) %>%
         stats::na.omit(i_dat)
    
      j_idx <- match(j_dat$start, j_dat$start)
      
      j_dat <- j_dat %>% dplyr::mutate(paired_pos = final_helix_df$j[j_idx]) %>%
         dplyr::mutate(start = paired_pos, end = paired_pos + width, rname = chrom_name) %>%
         #dplyr::select(-c(r1_start, r1_width, r1_end, paired_pos)) %>%
         na.omit(j_dat)
      
      if(nrow(i_dat) > 2000 | nrow(j_dat) > 2000){
         #shuffle order of dts and randomly select sample
         #because find_overlaps can't handle large vector sizes....
         i_dat <- i_dat[sample(1:nrow(i_dat)),]
         j_dat <- j_dat[sample(1:nrow(j_dat)),]
         i_dat <- utils::head(i_dat, 2000)
         j_dat <- utils::head(j_dat, 2000)
      }
   
      i_j_overlaps <- find_overlaps(i_dat, j_dat)
   
   } else {
      i_j_overlaps <- data.frame(r1_start = 0, r1_width = 0, r1_end = 0,
                                 r2_start = 0, r2_width = 0, r2_end = 0)
   }

   #return table of overlapping read pairs
   return(i_j_overlaps)
   
}
