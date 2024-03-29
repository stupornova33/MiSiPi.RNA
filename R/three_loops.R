#' run in process miRNA if three loops
#' @param vienna_vec a string
#' @param orig_loop_pos a string
#' @param converted a string
#' @return res_list a list
#' @export

three_loops <- function(vienna_vec, orig_loop_pos, converted){
   num_loops <- (length(orig_loop_pos)/2)
      print(paste0('num loops: ', num_loops))
      vienna <- collapse(vienna_vec)
      mid_loop_pos <- (length(orig_loop_pos) /2) - 1 #the likely terminal loop
      mid_loop_1 <- orig_loop_pos[mid_loop_pos] #the left loop coord
      mid_loop_2 <- orig_loop_pos[mid_loop_pos + num_loops] #right loop coord
      left_loop <- orig_loop_pos[mid_loop_pos - 1] #the loop to the left of the terminal loop
      right_loop <- orig_loop_pos[mid_loop_pos + 1] #the loop to the right
      
      left_arm <- vienna_vec[1:mid_loop_1]
      right_arm <- vienna_vec[mid_loop_2:length(vienna_vec)]
      
      #reverse right arm so pairs sort of line up
      right_arm <- rev(right_arm)
      left_vienna <- collapse(left_arm)
      right_vienna <- collapse(right_arm)
      split_seq <- expand(converted)
   
   
      tmp_pos <- unlist(stringi::stri_locate_all_regex(vienna, "\\)\\.*\\("))
       
      tmp_pos_1 <- c(unlist(stringi::stri_locate_all_regex(vienna, "\\)\\.*\\("))[1], unlist(stringi::stri_locate_all_regex(vienna, "\\)\\.*\\("))[3])
      tmp_pos_2 <- c(unlist(stringi::stri_locate_all_regex(vienna, "\\)\\.*\\("))[2], unlist(stringi::stri_locate_all_regex(vienna, "\\)\\.*\\("))[4])

      
      split_seq <- expand(converted)
      loop_seq <- split_seq[(mid_loop_1 + 1):(mid_loop_2 - 1)]
     
      
      left_stem_seq <- split_seq[tmp_pos_1[2]:mid_loop_1]
      right_stem_seq <- split_seq[mid_loop_2:(mid_loop_2 + length(left_stem_seq) -1)]
      right_stem_seq <- rev(right_stem_seq)
      
      l_other_seq <- split_seq[1:(tmp_pos_1[2] - 1)]
      r_other_seq <- split_seq[(mid_loop_2 + length(left_stem_seq)):length(split_seq)]
      r_other_seq <- rev(r_other_seq)
      
      new_left_vienna <- vienna_vec[tmp_pos_1[2]:mid_loop_1]
      new_right_vienna <- vienna_vec[mid_loop_2:(mid_loop_2 + length(left_stem_seq) -1)]
   
      new_left_vienna <- collapse(new_left_vienna)
      new_right_vienna <- collapse(new_right_vienna)
 
      
      stem_list <- align_stems(left_stem_seq, right_stem_seq, num_loops)
      left_arm <- stem_list$left_arm
      right_arm <- stem_list$right_arm
      
     
      major_loop_seq <- split_seq[(mid_loop_1 + 1):(mid_loop_2 - 1)]
      
      left_dashes <- unlist(stringi::stri_locate_all_regex(left_arm, "-")) %>% unique()
      right_dashes <-  unlist(stringi::stri_locate_all_regex(right_arm, "-")) %>% unique()
      
      
      df <- data.frame(left = expand(left_arm), right = expand(right_arm))
      match_pos <- which(df$left == df$right)
      
      
      right_arm <- Biostrings::RNAString(right_arm)
      right_arm <- Biostrings::complement(right_arm)
      right_arm <- as.character(unlist(unname(right_arm)))
      
      tmp <- numeric(length(expand(left_arm)))
      tmp[match_pos] <- 'A'
      dash_pos <- unlist(stringi::stri_locate_all_regex(left_arm, "-")) %>% unique()
      tmp[dash_pos] <- "-"
      new_left_vienna <- gsub("0", "S", tmp)
      
      tmp <- numeric(length(expand(right_arm)))
      tmp[match_pos] <- 'A'
      dash_pos <- unlist(stringi::stri_locate_all_regex(right_arm, "-")) %>% unique()
      new_right_vienna <- gsub("0", "S", tmp)
      
      
      
      #insert dashes into sequence in correct place
      if('FALSE' %in% is.na(left_dashes)){
         what <- rep("-", length(left_dashes))
         l_stem_seq <- insert_vector(left_stem_seq, left_dashes, what)
      } else {
         l_stem_seq <- left_stem_seq
      }
      
      if('FALSE' %in% is.na(right_dashes)){
         what <- rep("-", length(right_dashes)) 
         r_stem_seq <- insert_vector(right_stem_seq, right_dashes, what)
         r_stem_seq <- rev(r_stem_seq)
      } else {
         r_stem_seq <- right_stem_seq
         r_stem_seq <- rev(r_stem_seq)
      }
      
      #divide the loop sequence between the top and bottom
      left_loop_seq <- vector()
      right_loop_seq <- vector()
      
      div <- (length(major_loop_seq) - 2)/2
      if(length(major_loop_seq) %% 2 == 1){ #if loop length is odd
         loop_bool <- 'TRUE'
         half_len <- round(length(major_loop_seq)/2)
         div <- floor(div)
         
         left_loop_seq <- append(left_loop_seq, major_loop_seq[1:div])
         
         if(!div == 0){
            right_loop_seq <- append(right_loop_seq, major_loop_seq[((length(major_loop_seq) - div) + 1):length(major_loop_seq)])
            top_loop <- left_loop_seq
            bottom_loop <- right_loop_seq
         } else {
            right_loop_seq <- major_loop_seq[length(major_loop_seq)]
            top_loop <- left_loop_seq
            bottom_loop <- right_loop_seq
         }
         
         print(paste0('length left loop_seq :', length(left_loop_seq)))
         print(paste0('length right loop_seq : ', length(right_loop_seq)))
         
         #divide up loop sequence between arms
         
         mid_top_loop <- major_loop_seq[(div +1)]
         mid_loop <- major_loop_seq[(div + 2)]
         mid_bottom_loop <- major_loop_seq[(div + 3)]
         
      } else { #if length loop is even
         half_len <- length(major_loop_seq)/2
         loop_bool <- "FALSE"
         left_loop_seq <- append(left_loop_seq, major_loop_seq[1:(half_len - 1)])
         right_loop_seq <- append(right_loop_seq, major_loop_seq[((length(major_loop_seq) - half_len) + 2):length(major_loop_seq)])
         right_loop_seq <- rev(right_loop_seq)
         #loop_seq <- paste0(left_loop_seq, right_loop_seq)
         print(paste0('length left loop_seq :', length(left_loop_seq)))
         print(paste0('length right loop_seq : ', length(right_loop_seq)))
         
         #divide up loop sequence between arms
         top_loop <- left_loop_seq[1:div]
         mid_loop <- major_loop_seq[(div + 1):(div + 2)]
         bottom_loop <- right_loop_seq[1:length(right_loop_seq)]   
         mid_top_loop <- mid_loop[1]
         mid_bottom_loop <- mid_loop[2]
      }
      
      #determine which arm is longest to fill the other minor stem with "-"
      if(length(l_other_seq) > length(r_other_seq)){  
         l_bool <- 'FALSE'
         diff <- length(l_other_seq) - length(r_other_seq)
         dash_vec <- rep("-", times = diff)
         other_char_vec <- rep(" ", times = length(l_other_seq))
      } else if(length(r_other_seq) > length(l_other_seq)){
         l_bool <- 'TRUE'
         diff <- length(r_other_seq) - length(l_other_seq)
         dash_vec <- rep("-", times = diff)
         other_char_vec <- rep(" ", times = length(r_other_seq))
      } else { #they're the same length
         l_bool <- 'FALSE'
         dash_vec <- vector()
         other_char_vec <- rep(" ", times = length(l_other_seq))
      }
      
      
   
      print('making top line')
      #for the top line (just for stem loop structure)
      dot_pos <- unlist(stringi::stri_locate_all_regex(collapse(new_left_vienna), "S")) %>% unique()
      dash_pos <- unlist(stringi::stri_locate_all_regex(collapse(new_left_vienna), "-")) %>% unique()
   
      tmp <- numeric(length(l_stem_seq))
      tmp[dot_pos] <- l_stem_seq[dot_pos]

      tmp[dash_pos] <- "-"
      tmp[is.na(tmp)] <- 0
      top_vec <- gsub("0", " ", tmp)
      
      if(l_bool == 'TRUE'){
         #add the dashes where the deleted sequence is
         top_vec <- c(dash_vec, l_other_seq, top_vec)
         other_char_vec <- rep(' ', times = length(r_other_seq))
         
      } else {
         top_vec <- c(l_other_seq, top_vec)
         other_char_vec <- rep(' ', times = length(l_other_seq))
      }
      
      print('making mid_top')
      #mid top line
      paired_pos <- unlist(stringi::stri_locate_all_regex(collapse(new_left_vienna), "A")) %>% unique()
   
      tmp <- numeric(length(l_stem_seq))
      tmp[paired_pos] <- l_stem_seq[paired_pos]
      mid_top_vec <- gsub("0", " ", tmp)
      mid_top_vec <- c(other_char_vec, mid_top_vec)
      
      print('making middle')
      #middle line
      tmp <- numeric(length(l_stem_seq))
      tmp[paired_pos] <- "|"
      char_vec <- gsub("0", " ", tmp)
      char_vec <- c(other_char_vec, char_vec)
      
      #mid bottom
      print('making mid_bottom')
      tmp <- numeric(length(l_stem_seq))
      paired_pos <- unlist(stringi::stri_locate_all_regex(collapse(new_right_vienna), "A")) %>% unique()
      tmp[paired_pos] <- r_stem_seq[paired_pos]
      mid_bottom_vec <- gsub("0", " ", tmp)
      mid_bottom_vec <- c(other_char_vec, mid_bottom_vec)
      
      print('making bottom')
      #bottom
      tmp <- numeric(length(l_stem_seq))
      dot_pos <- unlist(stringi::stri_locate_all_regex(collapse(new_right_vienna), "S")) %>% unique()
      dash_pos <- unlist(stringi::stri_locate_all_regex(collapse(new_right_vienna), "-")) %>% unique()
      tmp[dot_pos] <- r_stem_seq[dot_pos]
      tmp[dash_pos] <- "-"
      bottom_vec <- gsub("0", " ", tmp)
      
      if(l_bool == TRUE){
         bottom_vec <- c(r_other_seq, bottom_vec)
      } else {
         bottom_vec <- c(dash_vec, r_other_seq, bottom_vec)
         
      }
      
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
