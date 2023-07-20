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
   #this may become a problem?
   new_startdf <- data.frame(matrix(ncol = 4, nrow = 0))
   new_enddf <- data.frame(matrix(ncol = 4, nrow = 0))
   for(i in 1:nrow(startdf)){
     print(paste0("start$i : ", startdf$i[i], " startdf$j : ", startdf$j[i]))
     if(startdf$i[i] != startdf$j[i])
       new_startdf <- rbind(new_startdf, startdf[i,])
   }
   for(i in 1:nrow(enddf)){
     print(paste0("enddf$i : ", enddf$i[i], " enddf$j : ", enddf$j[i]))
     if(enddf$i[i] != enddf$j[i])
       new_enddf <- rbind(new_enddf, enddf[i,])
   }

   startdf <- new_startdf
   enddf <- new_enddf

   print(startdf)
   print(enddf)

   tmp_df <- rbind(startdf, helix_df)
   print("binding helix")
   final_helix_df <- rbind(tmp_df, enddf)
   print("filtering helix j - i > 15")
   #remove results where segment of paired bases is less than 15nt
   final_helix_df <- final_helix_df %>% dplyr::filter(j - i > 15)

   if(nrow(final_helix_df) == 0){
      i_j_overlaps <- data.frame(r1_start = 0, r1_width = 0, r1_end = 0,
                                 r2_start = 0, r2_width = 0, r2_end = 0)
      return(i_j_overlaps)
   }
   print('making grouped_helix')
   grouped_helix <- group_helix_res(final_helix_df$i, final_helix_df$j)
   print('making filter_helix')
   filter_helix <- grouped_helix %>% dplyr::filter(X.End - X.Start > 17 | Y.Start - Y.End > 17)

   i_dat <- data.frame(rname = numeric(0), start = numeric(0), end = numeric(0), width = numeric(0))
   j_dat <- data.frame(rname = numeric(0), start = numeric(0), end = numeric(0), width = numeric(0))
   print('dicer_dt')
   dicer_dt <- dicer_dt %>%
     data.frame() %>%
     dplyr::select(c(start, end, rname))

   print('getting i,j idx')
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

       print('setting paired pos')
       i_dat$paired_pos <- final_helix_df$i[i_idx]

      #get reads that start or end at a paired position
       i_dat <- i_dat[,c("paired_pos", "end", "rname")]
       i_dat$width <- (i_dat$end - i_dat$paired_pos + 1)
       i_dat <- i_dat[,c("paired_pos", "width", "rname")]
       i_dat$end <- (i_dat$paired_pos + i_dat$width - 1)

      i_dat <- i_dat %>% dplyr::rename("start" = paired_pos)

      j_dat$paired_pos <- final_helix_df$j[j_idx]


      j_dat <- j_dat[c("paired_pos", "end", "rname")]
      j_dat$width <- (j_dat$end - j_dat$paired_pos + 1)
      j_dat <- j_dat[,c("paired_pos", "width", "rname")]
      j_dat$end <- (j_dat$paired_pos + j_dat$width - 1)


      j_dat <- j_dat %>% dplyr::rename("start" = paired_pos)

      if(nrow(i_dat) > 2000 | nrow(j_dat) > 2000){
         #shuffle order of dts and randomly select sample
         #because find_overlaps can't handle large vector sizes....
         i_dat <- i_dat[sample(1:nrow(i_dat)),]
         j_dat <- j_dat[sample(1:nrow(j_dat)),]
         i_dat <- utils::head(i_dat, 2000)
         j_dat <- utils::head(j_dat, 2000)
      }
      print('omitting NA')
      i_dat <- i_dat %>% stats::na.omit(i_dat)
      j_dat <- j_dat %>% stats::na.omit(j_dat)
      print('finding overlaps')
      print("Print this statement")
      print(head(i_dat))
      print(head(j_dat))
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
