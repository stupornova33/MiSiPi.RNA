#' calculates overlaps for Dicer signature
#' @param dicer_dt a data frame
#' @param helix_df a data frame
#' @param chrom_name a string
#' @param reg_start an integer
#' @return all_overlaps

#' @export

dicer_overlaps <- function(dicer_dt, helix_df, chrom_name, reg_start){
   # convert helix positions to chromosome positions
   bin <- j <- X.End <- X.Start <- Y.Start <- Y.End <- paired_pos <- start <- width <- NULL

   helix_df <- helix_df %>% dplyr::mutate(i = i + reg_start, j = j + reg_start)
   ### paired bases are several nucleotides diff. than read starts... pad the ends of the helix
   if(nrow(helix_df) == 0){
     return(i_j_overlaps <- data.frame(r1_start = 0, r1_width = 0, r1_end = 0,
                                       r2_start = 0, r2_width = 0, r2_end = 0))
   }
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
   #check to make sure start and end values are not the same. If they are, remove

   new_startdf <- data.frame(matrix(ncol = 4, nrow = 0))
   new_enddf <- data.frame(matrix(ncol = 4, nrow = 0))
   for(i in 1:nrow(startdf)){
     if(startdf$i[i] != startdf$j[i])
       new_startdf <- rbind(new_startdf, startdf[i,])
   }
   for(i in 1:nrow(enddf)){
     if(enddf$i[i] != enddf$j[i])
       new_enddf <- rbind(new_enddf, enddf[i,])
   }

   startdf <- new_startdf
   enddf <- new_enddf


   tmp_df <- rbind(startdf, helix_df)

   final_helix_df <- rbind(tmp_df, enddf)

   #remove results where segment of paired bases is less than 15nt
   final_helix_df <- final_helix_df %>% dplyr::filter(j - i > 15)

   if(nrow(final_helix_df) == 0){
      i_j_overlaps <- data.frame(r1_start = 0, r1_width = 0, r1_end = 0,
                                 r2_start = 0, r2_width = 0, r2_end = 0)
      return(i_j_overlaps)
   }

   grouped_helix <- group_helix_res(final_helix_df$i, final_helix_df$j)
   filter_helix <- grouped_helix %>% dplyr::filter(X.End - X.Start > 17 | Y.Start - Y.End > 17)


   i_dat <- data.frame(start = numeric(0), end = numeric(0), rname = character(0))
   j_dat <- data.frame(start = numeric(0), end = numeric(0), rname = character(0))

   #print('dicer_dt')
   #print(head(dicer_dt))
   dicer_dt <- dicer_dt %>%
     data.frame() %>%
     dplyr::select(c(start, end, rname))

   if(nrow(filter_helix) > 0){
     for(i in 1:nrow(filter_helix)){

       # Make a range of position numbers from x.start to x.end for each row
       # If any of the dicer_dt start or end positions fall in that range, then store the index of which ones those are
       # Take the dicer_dt rows that have positions that fall in the range, and append them to a table called i_dat
       # Essentially filter dicer_dt for rows that have start or end positions withing the range of the current iteration of the helix range
       x_rng <- seq(filter_helix$X.Start[i], filter_helix$X.End[i])

       # Now do the exact same thing for y.start and y.end for each row and store those reads in j_dat
       #y_rng <- seq(filter_helix$Y.End[i], filter_helix$Y.Start[i]) # Seq doesn't need to be reversed here. if you're just searching for a match in the range, the range direction is irrelevent
       y_rng <- seq(filter_helix$Y.Start[i], filter_helix$Y.End[i])

       i_idx <- which(dicer_dt$start %in% x_rng)
       j_idx <- which(dicer_dt$start %in% y_rng)

       i_dat <- rbind(i_dat, dicer_dt[i_idx,])
       j_dat <- rbind(j_dat, dicer_dt[j_idx,])
     }

     # Now lets calculate the width and add that column.
     # This will be used to translate the read start and stop positions soon
     i_dat <- i_dat %>%
       dplyr::mutate(width = end - start + 1)

     #since we're converting to the paired position, we should take $j for i_dat
     i_dat$paired_pos <- final_helix_df$j[i_idx]

     i_dat <- i_dat %>%
       dplyr::mutate(paired_end = paired_pos + width - 1)
     i_dat <- i_dat %>%
       dplyr::select(rname, paired_pos, paired_end)
     i_dat <- i_dat %>%
       dplyr::rename(start = paired_pos, end = paired_end)

     j_dat <- j_dat %>%
       dplyr::relocate(rname)

      if(nrow(i_dat) > 2000 | nrow(j_dat) > 2000){
         #shuffle order of dts and randomly select sample
         #because find_overlaps can't handle large vector sizes....
         i_dat <- i_dat[sample(1:nrow(i_dat)),]
         j_dat <- j_dat[sample(1:nrow(j_dat)),]
         i_dat <- utils::head(i_dat, 2000)
         j_dat <- utils::head(j_dat, 2000)
      }
      i_dat <- i_dat %>% stats::na.omit(i_dat)
      j_dat <- j_dat %>% stats::na.omit(j_dat)

      if(!nrow(i_dat) == 0 && !nrow(j_dat) == 0){
        i_j_overlaps <- find_overlaps(i_dat, j_dat)
      } else {
        i_j_overlaps <- data.frame(r1_start = 0, r1_width = 0, r1_end = 0,
                                   r2_start = 0, r2_width = 0, r2_end = 0)
      }

   } else {
      i_j_overlaps <- data.frame(r1_start = 0, r1_width = 0, r1_end = 0,
                                 r2_start = 0, r2_width = 0, r2_end = 0)
   }

   #return table of overlapping read pairs
   return(i_j_overlaps)

}
