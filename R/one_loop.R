#only one loop script
#' run in process miRNA if one loop
#' @param vienna_vec a string
#' @param orig_loop_pos a string
#' @param converted a string
#' @return res_list a list
#' @export


one_loop <- function(vienna_vec, orig_loop_pos, converted){
   vienna <- collapse(vienna_vec)
   left_arm <- vienna_vec[1:orig_loop_pos[1]]
   right_arm <- vienna_vec[orig_loop_pos[2]:length(vienna_vec)]
   loop_pos <- orig_loop_pos
   num_loops <- length(orig_loop_pos)/2
   #reverse right arm so pairs sort of line up
   right_arm <- rev(right_arm)
   new_left_vienna <- left_arm
   new_right_vienna <- right_arm
   
   #get sequence 
   split_seq <- expand(converted)
   
      
   left_seq <- split_seq[1:loop_pos[1]]
   right_seq <- split_seq[loop_pos[2]:length(split_seq)]
   right_seq <- rev(right_seq)
   
   stem_list <- .align_stems(right_seq, left_seq, num_loops)
   
   left_arm <- stem_list$left_arm
   right_arm <- stem_list$right_arm
 
   
   df <- data.frame(left = expand(left_arm), right = expand(right_arm))
   match_pos <- which(df$left == df$right)
   
   
   right_arm <- Biostrings::RNAString(right_arm)
   right_arm <- Biostrings::complement(right_arm)
   right_arm <- as.character(unlist(unname(right_arm)))
      
   
   left_dashes <- unlist(stringi::stri_locate_all_regex(left_arm, "-")) %>% unique()
   right_dashes <-  unlist(stringi::stri_locate_all_regex(right_arm, "-")) %>% unique()
   
   left_seq <- expand(left_arm)
   right_seq <- expand(right_arm)
   
   #create new vienna with dashes
   tmp <- numeric(length(left_seq))
   tmp[match_pos] <- "A"
   tmp[left_dashes] <- "-"
   new_left_vienna <- gsub("0", "S", tmp)
   
   tmp <- numeric(length(right_seq))
   tmp[match_pos] <- "A"
   tmp[right_dashes] <- "-"
   new_right_vienna <- gsub("0", "S", tmp)
   
   
   #make df for loop
   loop_seq <- split_seq[(loop_pos[1] + 1):(loop_pos[length(loop_pos)] - 1)]
   
   left_loop_seq <- vector()
   right_loop_seq <- vector()
   
   div <- (length(loop_seq) - 2)/2
   if(length(loop_seq) %% 2 == 1){ #if loop length is odd
      loop_bool <- 'TRUE'
      half_len <- round(length(loop_seq)/2)
      div <- floor(div)
      
      left_loop_seq <- append(left_loop_seq, loop_seq[1:div])
      
      if(!div == 0){
         right_loop_seq <- append(right_loop_seq, loop_seq[((length(loop_seq) - div) + 1):length(loop_seq)])
         top_loop <- left_loop_seq
         bottom_loop <- right_loop_seq
      } else {
         right_loop_seq <- loop_seq[length(loop_seq)]
         top_loop <- left_loop_seq
         bottom_loop <- right_loop_seq
      }
      
      print(paste0('length left loop_seq :', length(left_loop_seq)))
      print(paste0('length right loop_seq : ', length(right_loop_seq)))
      
      #divide up loop sequence between arms
      
      mid_top_loop <- loop_seq[(div +1)]
      mid_loop <- loop_seq[(div + 2)]
      mid_bottom_loop <- loop_seq[(div + 3)]
      
   } else { #if length loop is even
      half_len <- length(loop_seq)/2
      loop_bool <- "FALSE"
      left_loop_seq <- append(left_loop_seq, loop_seq[1:(half_len - 1)])
      right_loop_seq <- append(right_loop_seq, loop_seq[((length(loop_seq) - half_len) + 2):length(loop_seq)])
      right_loop_seq <- rev(right_loop_seq)
      #loop_seq <- paste0(left_loop_seq, right_loop_seq)
      print(paste0('length left loop_seq :', length(left_loop_seq)))
      print(paste0('length right loop_seq : ', length(right_loop_seq)))
      
      #divide up loop sequence between arms
      top_loop <- left_loop_seq[1:div]
      mid_loop <- loop_seq[(div + 1):(div + 2)]
      bottom_loop <- right_loop_seq[1:length(right_loop_seq)]
      mid_top_loop <- mid_loop[1]
      mid_bottom_loop <- mid_loop[2]
   }
   
   
   if(length(left_seq) > length(right_seq)){
      l_bool <- "TRUE"
      longest_arm <- left_seq
      shortest_arm <- right_seq
   } else if (length(right_seq > left_seq)){
      l_bool <- "FALSE"
      longest_arm <- right_seq
      shortest_arm <- left_seq
   } else { #if equal lengths just pick one
      l_bool <- 'FALSE'
      longest_arm <- right_seq
      shortest_arm <- left_seq
   }
   
   print(paste0('length left seq: ', length(left_seq)))
   print(paste0('length right seq: ', length(right_seq)))
   other_char_vec <- rep(" ", times = length(longest_arm))
   seq_diff <- length(longest_arm) - length(shortest_arm)
   add_dashes <- rep("-", times = seq_diff)
   
   #tmp <- numeric(length(expand(left_arm)))
   #tmp[match_pos] <- 'A'
   #dash_pos <- unlist(stringi::stri_locate_all_regex(left_arm, "-")) %>% unique()
   #tmp[dash_pos] <- "-"
   #new_left_vienna <- gsub("0", "S", tmp)
   
   #tmp <- numeric(length(expand(right_arm)))
   #tmp[match_pos] <- 'A'
   #dash_pos <- unlist(stringi::stri_locate_all_regex(right_arm, "-")) %>% unique()
   #new_right_vienna <- gsub("0", "S", tmp)
   
   
   print('making top line')
   #for the top line
   dot_pos <- unlist(stringi::stri_locate_all_regex(collapse(new_left_vienna), "S")) %>% unique()
   dash_pos <- unlist(stringi::stri_locate_all_regex(collapse(new_left_vienna), "-")) %>% unique()
   #paired_pos <- unlist(stringi::stri_locate_all_regex(collapse(new_left_vienna), "A")) %>% unique()
   
   tmp <- numeric(length(new_left_vienna))
   tmp[dot_pos] <- left_seq[dot_pos]
   tmp[dash_pos] <- "-"
   tmp[is.na(tmp)] <- 0
   top_vec <- gsub("0", " ", tmp)
   
   print('making mid_top')
   #mid top line
   paired_pos <- unlist(stringi::stri_locate_all_regex(collapse(new_left_vienna), "A")) %>% unique()
   tmp <- numeric(length(new_left_vienna))
   tmp[paired_pos] <- left_seq[paired_pos]
   mid_top_vec <- gsub("0", " ", tmp)
   
   print('making middle')
   #middle line
   tmp <- numeric(length(longest_arm))
   tmp[paired_pos] <- "|"
   char_vec <- gsub("0", " ", tmp)
   
   #mid bottom
   print('making mid_bottom')
   tmp <- numeric(length(new_right_vienna))
   paired_pos <- unlist(stringi::stri_locate_all_regex(collapse(new_right_vienna), "A")) %>% unique()
   tmp[paired_pos] <- right_seq[paired_pos]
   mid_bottom_vec <- gsub("0", " ", tmp)
   
   print('making bottom')
   #bottom
   tmp <- numeric(length(longest_arm))
   dot_pos <- unlist(stringi::stri_locate_all_regex(collapse(new_right_vienna), "S")) %>% unique()
   dash_pos <- unlist(stringi::stri_locate_all_regex(collapse(new_right_vienna), "-")) %>% unique()
   tmp[dot_pos] <- right_seq[dot_pos]
   tmp[dash_pos] <- "-"
   bottom_vec <- gsub("0", " ", tmp)
   
   
   if(loop_bool == 'TRUE'){ #add the nt at the very end of the bubble
      print("loop_bool == TRUE")
      top_vec <- c(top_vec, top_loop)
      mid_top <- c(mid_top_vec, rep(" ", length(top_loop)), mid_top_loop)
      char_vec <- c(char_vec, rep(" ", length(top_loop)), mid_loop)
      mid_bottom <- c(mid_bottom_vec, rep(" ", length(bottom_loop)), mid_bottom_loop)
      bottom_vec <- c(bottom_vec, bottom_loop)
   } else {
      top_vec <- c(top_vec, top_loop)
      mid_top <- c(mid_top_vec, rep(" ", length(top_loop)), mid_top_loop)
      char_vec <- c(char_vec, rep(" ", length(top_loop)), rep(" ", length(mid_top_loop)))
      mid_bottom <- c(mid_bottom_vec, rep(" ", length(bottom_loop)), mid_bottom_loop)
      bottom_vec <- c(bottom_vec, bottom_loop)
   }

   
   print('making all dfs')
   
   top_df <- data.frame("x" = c(1:length(top_vec)), "y" = c(5.1), "text" = top_vec)
   mid_top_df <- data.frame("x" = c(1:length(mid_top)), "y" = c(4.9), "text" = mid_top)
   char_df <- data.frame("x" = c(1:length(char_vec)), "y" = c(4.7), "text" = char_vec)
   mid_bottom_df <- data.frame("x" = c(1:length(mid_bottom)), "y" = c(4.5), "text" = mid_bottom)
   bottom_df <- data.frame("x" = c(1:length(bottom_vec)), "y" = c(4.3), "text" = bottom_vec)
   #f_df <- rbind(top_df, mid_top_df, char_df, mid_bottom_df, bottom_df)
   res_list <- list(top = top_df, mid_top = mid_top_df, char_df = char_df, mid_bottom = mid_bottom_df, bottom = bottom_df)
   return(res_list)

}

