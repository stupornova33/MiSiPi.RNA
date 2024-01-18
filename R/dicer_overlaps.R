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
   print(helix_df)
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

   #add the new intervals to helix_df
   tmp_df <- rbind(startdf, helix_df)

   final_helix_df <- rbind(tmp_df, enddf)

   #remove results where segment of paired bases is less than 15nt
   # for example the loop sequence
   final_helix_df <- final_helix_df %>% dplyr::filter(j - i > 15)

   if(nrow(final_helix_df) == 0){
      i_j_overlaps <- data.frame(r1_start = 0, r1_width = 0, r1_end = 0,
                                 r2_start = 0, r2_width = 0, r2_end = 0)
      return(i_j_overlaps)
   }

   grouped_helix <- group_helix_res(final_helix_df$i, final_helix_df$j)
   filter_helix <- grouped_helix %>% dplyr::filter(X.End - X.Start > 17 | Y.Start - Y.End > 17)
   #write.table(grouped_helix, "grouped_helix.txt", sep = "\t", row.names = FALSE, quote = FALSE)

   i_dat <- data.frame(start = numeric(0), end = numeric(0), rname = character(0), paired_pos = numeric(0))
   j_dat <- data.frame(start = numeric(0), end = numeric(0), rname = character(0), paired_pos = numeric(0))

   #print('dicer_dt')
   #print(head(dicer_dt))
   dicer_dt <- dicer_dt %>%
     data.frame() %>%
     dplyr::select(c(start, end, rname))



##################### OLD CODE #################
# leaving for now JIC
   #global_i_idx <- numeric()

   #if(nrow(filter_helix) > 0){
  #   for(i in 1:nrow(filter_helix)){

       # Make a range of position numbers from x.start to x.end for each row
       # If any of the dicer_dt start or end positions fall in that range, then store the index of which ones those are
       # Take the dicer_dt rows that have positions that fall in the range, and append them to a table called i_dat
       # Essentially filter dicer_dt for rows that have start or end positions withing the range of the current iteration of the helix range
   #    x_rng <- seq(filter_helix$X.Start[i], filter_helix$X.End[i])
  #     print("x_rng: ")
  #     print(x_rng)
       # Now do the exact same thing for y.start and y.end for each row and store those reads in j_dat
       #y_rng <- seq(filter_helix$Y.End[i], filter_helix$Y.Start[i]) # Seq doesn't need to be reversed here. if you're just searching for a match in the range, the range direction is irrelevent

  #     y_rng <- seq(filter_helix$Y.Start[i], filter_helix$Y.End[i])
  #     print("y_rng: ")
  #     print(y_rng)
  #     i_idx <- which(dicer_dt$start %in% x_rng)

       #keep track of which reads are found in the helix rngs over all iterations
       #for future use
  #     global_i_idx <- append(global_i_idx, i_idx)
  #     j_idx <- which(dicer_dt$start %in% y_rng)

  #     i_dat <- rbind(i_dat, dicer_dt[i_idx,])
  #     j_dat <- rbind(j_dat, dicer_dt[j_idx,])
  #   }

   ## Start
   if(nrow(filter_helix) > 0){
      for(i in 1:nrow(filter_helix)){
        #print(paste0("i: ", i))
        x_rng <- seq(filter_helix$X.Start[i], filter_helix$X.End[i])

        #y_rng <- seq(filter_helix$Y.End[i], filter_helix$Y.Start[i])
        y_rng <- seq(filter_helix$Y.Start[i], filter_helix$Y.End[i])

        #get i reads that start at a paired pos in helix_df
        x_idx <- which(dicer_dt$start %in% x_rng)
        x_dat <- dicer_dt[x_idx,] %>%
          dplyr::arrange(start)
        #get j reads that start at a paired pos in helix_df
        y_idx <- which(dicer_dt$start %in% y_rng)

        #don't need to transform j reads, just extract
        y_dat <- dicer_dt[y_idx,]

        paired_pos <- vector() # Possibly move above into i loop
        prior_i_idx <- integer()

        if(nrow(x_dat) > 0){
        for (j in 1:nrow(x_dat)) {
          #print(paste0("j: ", j))
          i_idx <- which(final_helix_df$i == x_dat$start[j])
          #if no results returned
          if(identical(i_idx, integer(0)) & !(j == 1)) {
             #get the index of the last paired pos
             i_idx <- prior_i_idx
            #else if no results and this is the first iteration
          } else if (identical(i_idx, integer(0)) & (j == 1)) {
              found <- FALSE
              current_start <- x_dat$start[j]
              new_start <- current_start - 1

              while (!found) {
                #print("!found")
                i_idx <- which(final_helix_df$i == new_start)
                #print("while i_idx: ", i_idx)
                if (identical(i_idx, integer(0))) {
                  new_start <- new_start - 1
                } else {
                  found <- TRUE
                  prior_i_idx <- i_idx
                }

              }
            #print(paste0("i_idx: ", i_idx))
          } else {

            # Set prior i_idx if we have a valid one
            prior_i_idx <- i_idx
          }
          paired_pos <- append(paired_pos, final_helix_df$j[i_idx])
        } # End 'j' for loop
        x_dat$paired_pos <- paired_pos
        i_dat <- rbind(i_dat, x_dat)
        j_dat <- rbind(j_dat, y_dat)
        } else { #end if(nrow(x_dat) > 0)
            next
      }
     } # End 'i' for loop

   } # End if statement
   ## END

     i_dat <- i_dat %>%
       dplyr::mutate(width = end - start + 1)

     #since we're converting to the paired position, we should take $j for i_dat

     i_dat <- i_dat %>%
       dplyr::mutate(paired_start = paired_pos - width + 1)
     i_dat <- i_dat %>%
       dplyr::select(rname, paired_pos, paired_start)
     i_dat <- i_dat %>%
       dplyr::rename(start = paired_start, end = paired_pos)

     #write.table(i_dat, "i_dat_after_convert.txt", sep = "\t", row.names = FALSE, quote = FALSE)
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

   #return table of overlapping read pairs
   return(i_j_overlaps)

}
