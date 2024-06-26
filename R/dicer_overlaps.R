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

   helix_df <- helix_df %>% dplyr::mutate(i = i + (reg_start - 1), j = j + (reg_start - 1))
   ### paired bases are several nucleotides diff. than read starts... pad the ends of the helix
   if(nrow(helix_df) == 0){
     return(i_j_overlaps <- data.frame(r1_start = 0, r1_width = 0, r1_end = 0,
                                       r2_start = 0, r2_width = 0, r2_end = 0))
   }
   i_begin <- helix_df$i[1]
   i_end <- helix_df$i[nrow(helix_df)]
   j_begin <- helix_df$j[1]
   j_end <- helix_df$j[nrow(helix_df)]

   if(reg_start >= 8){
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

     #add the new intervals to helix_df
     tmp_df <- rbind(startdf, helix_df)

     final_helix_df <- rbind(tmp_df, enddf)

   } else {
     final_helix_df <- helix_df
   }


   #remove results where segment of paired bases is less than 15nt
   # for example the loop sequence
   final_helix_df <- final_helix_df %>% dplyr::filter(j - i > 15)

   if(nrow(final_helix_df) == 0){
      i_j_overlaps <- data.frame(r1_start = 0, r1_width = 0, r1_end = 0,
                                 r2_start = 0, r2_width = 0, r2_end = 0)
      return(i_j_overlaps)
   }

   # group together regions which are paired
   grouped_helix <- group_helix_res(final_helix_df$i, final_helix_df$j)
   filter_helix <- grouped_helix %>% dplyr::filter(X.End - X.Start > 17 | Y.Start - Y.End > 17)
   #write.table(grouped_helix, "grouped_helix.txt", sep = "\t", row.names = FALSE, quote = FALSE)

   i_dat <- data.frame(start = numeric(0), end = numeric(0), rname = character(0), paired_start = numeric(0), paired_end = numeric(0), num_shifted = numeric(0))
   j_dat <- data.frame(start = numeric(0), end = numeric(0), rname = character(0), paired_start = numeric(0), paired_end = numeric(0))

   #print('dicer_dt')
   #print(head(dicer_dt))
   dicer_dt <- dicer_dt %>%
     data.frame() %>%
     dplyr::select(c(start, end, rname))

  #i = 13
   ## Start
   if(nrow(filter_helix) > 0){
      for(i in 1:nrow(filter_helix)){
        print(paste0("i: ", i))
        #print(paste0("i: ", i))
        x_rng <- seq(filter_helix$X.Start[i], filter_helix$X.End[i])

        y_rng <- seq(filter_helix$Y.Start[i], filter_helix$Y.End[i])


        ############ testing #####################################
        #get i & j reads that start or end at pos in helix_df
        x_idx <- which(dicer_dt$start %in% x_rng | dicer_dt$end %in% x_rng)
        x_dat <- dicer_dt[x_idx,]

        #eliminate reads which extend to a highly unpaired region
        x_dat <- x_dat %>% dplyr::filter(start >= filter_helix$X.Start[i])

        y_idx <- which(dicer_dt$start %in% y_rng | dicer_dt$end %in% y_rng)

        #don't need to transform j reads, just extract
        y_dat <- dicer_dt[y_idx,]


        prior_end_idx <- integer()
        prior_start_idx <- integer()
        if(nrow(x_dat) > 0){
          for (j in 1:nrow(x_dat)) {

            #check if read start / read end in paired positions
            #this affects how the new paired position is determined
            #2/19 need to reverse this to get the ends of i that have paired pos

            start_idx <- which(final_helix_df$i == x_dat$start[j])
            end_idx <- which(final_helix_df$i == x_dat$end[j])

            #print(paste0("start_idx: ", start_idx))
            #print(paste0("end_idx: ", end_idx))

            #if read starts at paired pos but ends at non-paired pos
            if (identical(end_idx, integer(0)) & !identical(start_idx, integer(0))){
              found <- FALSE
              current_start <- x_dat$start[j]
              current_end <- x_dat$end[j]
              #print(paste0("current_end: ", current_end))
              new_end <- current_end - 1
              #print(paste0("new_end: ", new_end))
              while (!found) {
                #print("NOT FOUND")
                #print("!found")
                end_idx <- which(final_helix_df$i == new_end)
                #print("while end_idx: ", end_idx)
                if (identical(end_idx, integer(0))) {
                  new_end <- new_end - 1
                  shift <- new_end - current_end
                  #print(paste0("new_end: ", new_end))
                  #print(paste0("shift: ", shift))
                } else {
                  found <- TRUE
                  #print("FOUND")
                  prior_end_idx <- end_idx
                  shift <- new_end - current_end
                  #print(paste0("prior_end_idx: ", prior_end_idx))
                  #print(paste0("shift: ", shift))
                }

              }
              paired_start <- final_helix_df$j[end_idx] + shift
              paired_end <- final_helix_df$j[start_idx]

            } else if(!identical(end_idx, integer(0)) & identical(start_idx, integer(0))){ #if read ends at paired pos but starts at unpaired pos
              #end idx is found but start is not
              found <- FALSE

              # coordinates of the read
              current_start <- x_dat$start[j]
              current_end <- x_dat$end[j]
              #print(paste0("current_end: ", current_end))
              new_start <- current_start + 1
              #print(paste0("new_start: ", new_start))
              while (!found) {
                #print("NOT FOUND")
                #print("!found")
                start_idx <- which(final_helix_df$i == new_start)
                #print(paste0("start_idx: ", start_idx))
                #print("while end_idx: ", end_idx)
                if (identical(start_idx, integer(0))) {
                  new_start <- new_start + 1
                  shift <- new_start - current_start
                  #print(paste0("new_start: ", new_start))
                  #print(paste0("shift: ", shift))
                } else {
                  found <- TRUE
                  #print("FOUND")
                  prior_start_idx <- start_idx
                  shift <- new_start - current_start
                  #print(paste0("prior_start_idx: ", prior_start_idx))
                  #print(paste0("shift: ", shift))
                }

              }


              paired_end <- (final_helix_df$j[start_idx] + shift)
              paired_start <- (final_helix_df$j[end_idx])

            } else if(identical(start_idx, integer(0)) & identical(end_idx, integer(0))) { #if read begins and ends at unpaired pos
              current_start <- x_dat$start[j]
              current_end <- x_dat$end[j]
              rng <- seq(x_dat$start[j], x_dat$end[j])
              min_paired_idx <- min(which(rng %in% helix_df$i))
              max_paired_idx <- max(which(rng %in% helix_df$i))
              min_paired <- rng[min_paired_idx]
              max_paired <- rng[max_paired_idx]

              # determine whether the paired base before or after current_start is closer
              # avoids pairing starts or ends where there is a large section of unpaired bases between.

              alt_min_idx <- which(final_helix_df$i == min_paired)
              alt_min <- final_helix_df$i[alt_min_idx - 1]
              #alt max_idx?
              diff_alt <- abs(alt_min - current_start)
              diff_min <- abs(min_paired - current_start)

              if(diff_alt < diff_min){ #then use the previous paired pos
                left_bulge <- diff_alt
                right_bulge <- length(rng) - max_paired_idx
                min_paired <- alt_min
                min_paired_idx <- which(final_helix_df$i == min_paired)
                min_paired_j <- final_helix_df$j[min_paired_idx]
                min_paired_j_idx <- min_paired_idx
                max_paired_j_idx <- which(helix_df$i == max_paired)

              } else if(diff_min < diff_alt){ #then use the next paired pos
                left_bulge <- min_paired_idx - 1
                right_bulge <- length(rng) - max_paired_idx
                min_paired_j_idx <- which(helix_df$i == min_paired)
                max_paired_j_idx <- which(helix_df$i == max_paired)
              } else { #if the distance is equal use the min
                left_bulge <- min_paired_idx - 1
                right_bulge <- length(rng) - max_paired_idx
                min_paired_j_idx <- which(helix_df$i == min_paired)
                max_paired_j_idx <- which(helix_df$i == max_paired)
              }

              paired_end <- helix_df$j[min_paired_j_idx] + left_bulge
              paired_start <- helix_df$j[max_paired_j_idx] - right_bulge

            } else { #if read has paired start pos and paired end pos
              paired_start <- final_helix_df$j[end_idx]
              paired_end <- final_helix_df$j[start_idx]

            }

             x_dat$paired_start[j] <- paired_start
             x_dat$paired_end[j] <- paired_end

          } # End 'j' for loop

          i_dat <- rbind(i_dat, x_dat)
          j_dat <- rbind(j_dat, y_dat)
        } else { #end if(nrow(x_dat) > 0)
          next
        }
      } # End 'i' for loop

   }

     i_dat <- i_dat %>%
       dplyr::mutate(width = end - start + 1, paired_width = paired_end - paired_start + 1)

     i_dat <- i_dat %>%
       dplyr::select(rname, paired_start, paired_end, paired_width)


     i_dat <- i_dat %>%
       dplyr::rename(start = paired_start, end = paired_end, width = paired_width)

     # filter out transformed reads that are shorter than 18 and longer than 32
     # testing whether this improves dcr sig.
     i_dat <- i_dat %>% dplyr::filter(width >= 18 & width <= 32) %>% dplyr::select(-c(width))
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
        #i_j_overlaps <- find_overlaps(j_dat, i_dat)
        i_j_overlaps <- find_overlaps(i_dat, j_dat)
      } else {
        i_j_overlaps <- data.frame(r1_start = 0, r1_width = 0, r1_end = 0,
                                   r2_start = 0, r2_width = 0, r2_end = 0)
      }

   #return table of overlapping read pairs
   return(i_j_overlaps)

}
