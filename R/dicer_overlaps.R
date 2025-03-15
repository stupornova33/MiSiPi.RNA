# calculates overlaps for Dicer signature
#
# dicer_dt - Data Frame
# helix_df - Data Frame
# chrom_name - Character
# original_start - Integer
# return - Data Frame

.dicer_overlaps <- function(dicer_dt, helix_df, chrom_name, original_start) {
  bin <- j <- X.End <- X.Start <- Y.Start <- Y.End <- paired_pos <- start <- width <- NULL

  if (nrow(helix_df) == 0) {
    return(.null_dicer_res())
  }

  #### Convert helix positions to chromosome positions ####
  # Positions returned by [RNAFold.exe or R4RNA::viennaToHelix] are indexed at 1 being the smallest
  # But viennaToHelix does keep track of the start and stop positions
  # We've passed in the start position here with the parameter original_start
  # So we're adding original_start - 1 to each helix position to get the actual positions
  helix_df <- helix_df %>%
    dplyr::mutate(
      i = i + original_start - 1,
      j = j + original_start - 1
    )

  #### Pad Data Frame ####

  i_begin <- helix_df$i[1]
  i_end <- helix_df$i[nrow(helix_df)]
  j_begin <- helix_df$j[1]
  j_end <- helix_df$j[nrow(helix_df)]

  if (original_start < 8) {
    final_helix_df <- helix_df
  } else {
    # This section is just adding up to 16 observations to the helix_df
    # Up to 8 before the first observation and up to 8 after the last one
    # They are just sequential from the observations, and any observation
    # Where i & j are the same get filtered out
    i_start_vals <- seq(i_begin - 8, i_begin - 1)
    i_end_vals <- seq(i_end + 1, i_end + 8)
    j_start_vals <- seq(j_begin + 8, j_begin + 1)
    j_end_vals <- seq(j_end - 1, j_end - 8)

    startdf <- data.frame(
      i = i_start_vals,
      j = j_start_vals,
      length = 1,
      value = NA
    )

    enddf <- data.frame(
      i = i_end_vals,
      j = j_end_vals,
      length = 1,
      value = NA
    )

    # check to make sure start and end values are not the same. If they are, remove
    startdf <- startdf %>%
      dplyr::filter(i != j)
    
    enddf <- enddf %>%
      dplyr::filter(i != j)

    # add the new intervals to helix_df
    final_helix_df <- dplyr::bind_rows(startdf, helix_df, enddf)
  }

  #### Filter ####
  # Originally we were removing results where the segment of paired
  # bases was less than 15nt (for example, the loop sequence)
  # Currently, the filtering is just to ensure the paired position is
  # greater than the original position
  final_helix_df <- final_helix_df %>%
    dplyr::filter(j - i > 0)

  if (nrow(final_helix_df) == 0) {
    return(.null_dicer_res())
  }

  # group together regions which are paired
  grouped_helix <- group_helix_res(final_helix_df$i, final_helix_df$j)
  
  # Filter 
  filter_helix <- grouped_helix %>%
    dplyr::filter(X.End - X.Start > 17 | Y.Start - Y.End > 17)

  i_dat <- data.frame(start = numeric(0), end = numeric(0), rname = character(0), n = numeric(0), paired_start = numeric(0), paired_end = numeric(0))
  j_dat <- data.frame(start = numeric(0), end = numeric(0), rname = character(0), n = numeric(0), paired_start = numeric(0), paired_end = numeric(0))

  ## This is not selecting columns as the data frame is currently grouped. might need to ungroup first
  dicer_dt <- dicer_dt %>%
    dplyr::ungroup() %>%
    dplyr::select(c(start, end, rname, n))

  ## Start
  if (nrow(filter_helix) > 0) {
    for (i in 1:nrow(filter_helix)) {
      x_rng <- seq(filter_helix$X.Start[i], filter_helix$X.End[i])

      y_rng <- seq(filter_helix$Y.Start[i], filter_helix$Y.End[i])


      ############ testing #####################################
      # get i & j reads that start or end at pos in helix_df
      x_idx <- which(dicer_dt$start %in% x_rng | dicer_dt$end %in% x_rng)
      x_dat <- dicer_dt[x_idx, ]
      
      # R gives warnings when trying to add to a new column at an indexed position
      # Initializing the columns now
      x_dat$paired_start <- NA
      x_dat$paired_end <- NA
      

      # eliminate reads which extend to a highly unpaired region
      x_dat <- x_dat %>%
        dplyr::filter(start >= filter_helix$X.Start[i])

      y_idx <- which(dicer_dt$start %in% y_rng | dicer_dt$end %in% y_rng)

      # don't need to transform j reads, just extract
      y_dat <- dicer_dt[y_idx, ]


      prior_end_idx <- integer()
      prior_start_idx <- integer()
      if (nrow(x_dat) > 0) {
        for (j in 1:nrow(x_dat)) {
          # check if read start / read end in paired positions
          # this affects how the new paired position is determined
          # 2/19 need to reverse this to get the ends of i that have paired pos

          start_idx <- which(final_helix_df$i == x_dat$start[j])
          end_idx <- which(final_helix_df$i == x_dat$end[j])

          # if read starts at paired pos but ends at non-paired pos
          if (identical(end_idx, integer(0)) & !identical(start_idx, integer(0))) {
            found <- FALSE
            current_start <- x_dat$start[j]
            current_end <- x_dat$end[j]
            new_end <- current_end - 1
            while (!found) {
              end_idx <- which(final_helix_df$i == new_end)
              if (identical(end_idx, integer(0))) {
                new_end <- new_end - 1
                shift <- new_end - current_end
              } else {
                found <- TRUE
                prior_end_idx <- end_idx
                shift <- new_end - current_end
              }
            }
            paired_start <- final_helix_df$j[end_idx] + shift
            paired_end <- final_helix_df$j[start_idx]
          } else if (!identical(end_idx, integer(0)) & identical(start_idx, integer(0))) { # if read ends at paired pos but starts at unpaired pos
            # end idx is found but start is not
            found <- FALSE

            # coordinates of the read
            current_start <- x_dat$start[j]
            current_end <- x_dat$end[j]
            new_start <- current_start + 1
            while (!found) {
              start_idx <- which(final_helix_df$i == new_start)
              if (identical(start_idx, integer(0))) {
                new_start <- new_start + 1
                shift <- new_start - current_start
              } else {
                found <- TRUE
                prior_start_idx <- start_idx
                shift <- new_start - current_start
              }
            }

            paired_end <- (final_helix_df$j[start_idx] + shift)
            paired_start <- (final_helix_df$j[end_idx])
          } else if (identical(start_idx, integer(0)) & identical(end_idx, integer(0))) { # if read begins and ends at unpaired pos
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
            # alt max_idx?
            diff_alt <- abs(alt_min - current_start)
            diff_min <- abs(min_paired - current_start)

            if (diff_alt < diff_min) { # then use the previous paired pos
              left_bulge <- diff_alt
              right_bulge <- length(rng) - max_paired_idx
              min_paired <- alt_min
              min_paired_idx <- which(final_helix_df$i == min_paired)
              min_paired_j <- final_helix_df$j[min_paired_idx]
              min_paired_j_idx <- min_paired_idx
              max_paired_j_idx <- which(helix_df$i == max_paired)
            } else if (diff_min < diff_alt) { # then use the next paired pos
              left_bulge <- min_paired_idx - 1
              right_bulge <- length(rng) - max_paired_idx
              min_paired_j_idx <- which(helix_df$i == min_paired)
              max_paired_j_idx <- which(helix_df$i == max_paired)
            } else { # if the distance is equal use the min
              left_bulge <- min_paired_idx - 1
              right_bulge <- length(rng) - max_paired_idx
              min_paired_j_idx <- which(helix_df$i == min_paired)
              max_paired_j_idx <- which(helix_df$i == max_paired)
            }

            paired_end <- helix_df$j[min_paired_j_idx] + left_bulge
            paired_start <- helix_df$j[max_paired_j_idx] - right_bulge
          } else { # if read has paired start pos and paired end pos
            paired_start <- final_helix_df$j[end_idx]
            paired_end <- final_helix_df$j[start_idx]
          }

          x_dat$paired_start[j] <- paired_start
          x_dat$paired_end[j] <- paired_end
        } # End 'j' for loop

        i_dat <- rbind(i_dat, x_dat)
        j_dat <- rbind(j_dat, y_dat)
      } else { # end if(nrow(x_dat) > 0)
        next
      }
    } # End 'i' for loop
  }

  i_dat <- i_dat %>%
    dplyr::mutate(paired_width = paired_end - paired_start + 1) %>%
    dplyr::select(rname, paired_start, paired_end, paired_width, n) %>%
    dplyr::rename(start = paired_start, end = paired_end, width = paired_width)

  # filter out transformed reads that are shorter than 18 and longer than 32
  # testing whether this improves dcr sig.
  i_dat <- i_dat %>%
    dplyr::filter(width >= 18 & width <= 32) %>%
    dplyr::select(-width)

  j_dat <- j_dat %>%
    dplyr::relocate(rname)

  
  if (nrow(i_dat) == 0 || nrow(j_dat) == 0) {
    return(.null_dicer_res())
  }
  
  i_j_overlaps <- .find_overlaps(i_dat, j_dat)

  # return table of overlapping read pairs
  return(i_j_overlaps)
}

.null_dicer_res <- function() {
  i_j_overlaps <- data.frame(
    r1_start = 0, r1_width = 0, r1_end = 0, r1_dupes = 0,
    r2_start = 0, r2_width = 0, r2_end = 0, r2_dupes = 0
  )
  return(i_j_overlaps)
}
