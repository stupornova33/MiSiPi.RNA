#' calculates overlaps for Dicer signature
#' @param dicer_dt a data frame
#' @param helix_df a data frame
#' @param chrom_name a string
#' @param reg_start an integer
#' @return all_overlaps

#' @export
#if this all works we can get rid of dicer_overlaps and its associated Rcpp functions

new_dcr_overlaps <- function(dicer_dt, helix_df, chrom_name, reg_start){
  r2_dt <- dicer_dt
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

  #i <- 1758
  final_reads_df <- data.frame()
  for(i in 1:nrow(r2_dt)){
    print(i)
    read <- r2_dt[i,] %>% dplyr::select(rname, start, end, width, first, count)
    #read_tbl <- data.frame(i = c(seq(read$start, read$end)))

    helix_sub <- final_helix_df %>% dplyr::filter(i >= read$start & i <= read$end) %>% dplyr::select(-c(length, value))
    if(nrow(helix_sub) == 0){
      print("helix_sub is empty")
      next
    }
    start <- read$start[1]
    end <- read$end[1]

    start_i_idx <- which(helix_sub$i == start)
    end_i_idx <- which(helix_sub$i == end)

    helix_i <- helix_sub %>% dplyr::select(c(i)) %>% dplyr::rename("fake_i" = "i")

    helix_i <- helix_i %>% dplyr::mutate(i = fake_i)
    rng <- data.frame(i = c(seq(helix_sub$i[1], helix_sub$i[nrow(helix_sub)])))

    ## fill in gaps in helix$i by merging it to a range from first value in helix$i to last value in helix$i
    testmerge_i <- merge(helix_i, rng, by = "i", all = TRUE)
    testmerge_i <- testmerge_i %>% dplyr::rename("pos_i" = "i", "i" = "fake_i")

    # check to see if start is paired and end is unpaired
    if(!identical(start_i_idx, integer(0)) & identical(end_i_idx, integer(0))) {
      #print("Read start is paired and read end is unpaired.")
      paired_start_idx <- which(final_helix_df$i == read$start)
      paired_start <- final_helix_df$j[paired_start_idx]

      #find the last paired position of read
      last_paired_i <- testmerge_i$i[nrow(testmerge_i)]
      j_idx <- which(final_helix_df$i == last_paired_i)

      last_paired_j <- final_helix_df$j[j_idx]
      shift <- read$end - last_paired_i

      paired_end_idx <- which(final_helix_df$j == last_paired_j)
      paired_end <- final_helix_df$j[paired_end_idx] - shift

      #read$paired_start <- paired_end
      #read$paired_end <- paired_start

    } else if(identical(start_i_idx, integer(0)) & !identical(end_i_idx, integer(0))) { # if end is paired and start is unpaired
      #print("Read start is unpaired and read end is paired.")
      paired_end_idx <- which(final_helix_df$i == read$end)
      paired_end <- final_helix_df$j[paired_end_idx]

      #find the first paired position of read
      first_paired_i <- testmerge_i$i[1]
      j_idx <- which(final_helix_df$i == first_paired_i)

      first_paired_j <- final_helix_df$j[j_idx]
      shift <- first_paired_i - read$start

      paired_start_idx <- which(final_helix_df$j == first_paired_j)
      paired_start <- final_helix_df$j[paired_start_idx] + shift

      #read$paired_start <- paired_end
      #read$paired_end <- paired_start

    } else if(!identical(start_i_idx, integer(0)) & !identical(end_i_idx, integer(0))) { # if both are paired
      #print("Read start and end are both paired.")
      paired_start <- helix_sub$j[end_i_idx]
      paired_end <- helix_sub$j[start_i_idx]

    } else { # both are unpaired
      #get the first and last paired positions
      # take difference of read_start vs first paired and read_end vs last paired
      print("Read start and end are both unpaired.")
      first_paired <- testmerge_i$i[1]
      last_paired <- testmerge_i$i[nrow(testmerge_i)]
      start_shift <- first_paired - read$start
      end_shift <- read$end - last_paired
      paired_start <- helix_sub$j[1] + start_shift
      paired_end <- helix_sub$j[nrow(helix_sub)] - end_shift
      #read$paired_start <- paired_end
      #read$paired_end <- paired_start

    }
    if(paired_start < paired_end){
      read$paired_start <- paired_start
      read$paired_end <- paired_end
    } else {
      read$paired_start <- paired_end
      read$paired_end <- paired_start
    }

#################### End dicer transform. Find overlapping reads

    final_reads_df <- rbind(read, final_reads_df)
  }

  if(nrow(final_reads_df) == 0){
    print("No ")
    return(i_j_overlaps <- data.frame(r1_start = 0, r1_width = 0, r1_end = 0,
                                      r2_start = 0, r2_width = 0, r2_end = 0))
  }

  i_dat <- final_reads_df %>%
    dplyr::mutate(width = end - start + 1, paired_width = paired_end - paired_start + 1) %>% dplyr::ungroup()

  i_dat <- i_dat %>%
    dplyr::select(c(rname, paired_start, paired_end, paired_width, first, count)) %>% dplyr::rename("start" = "paired_start", "end" = "paired_end", "width" = "paired_width")

  expand_i_dat <- rep_nonseq_reads(i_dat$count, i_dat$rname, i_dat$start, i_dat$end, i_dat$first)
  j_dat <- rep_nonseq_reads(r2_dt$count, r2_dt$rname, r2_dt$start, r2_dt$end, r2_dt$first)

  if(nrow(expand_i_dat) > 2000 | nrow(j_dat) > 2000){
    #shuffle order of dts and randomly select sample
    #because find_overlaps can't handle large vector sizes....
    i_dat <- expand_i_dat[sample(1:nrow(expand_i_dat)),]
    j_dat <- j_dat[sample(1:nrow(j_dat)),]
    i_dat <- utils::head(expand_i_dat, 2000)
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
