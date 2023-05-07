#' function to find max_pi heatmap cell values
#' takes a matrix created from piRNA function
#' finds max row, column, and count of highest cell
#' returns highest_pi_row, highest_pi_col, second_pi_row, second_pi_col, highest_pi_count, second_pi_count as a data frame
#' 
#' @param si_res a data table
#' @return max_pi_df



#' @export 
#' 
#get_max_pi_heat<- function(pi_res){
#   max_results <- unlist(unname(which(pi_res == max(pi_res), arr.ind = TRUE)))
#   if(length(max_results) == 0){
#      highest_size <- NA
#      second_highest_size <- NA
#      highest_count <- NA
#      second_highest_count <- NA
      
#   } else if (length(max_results) == 2) {
#      highest_size <- max_results + 15
#      pos <- c(max_results[1,1], max_results[1,2])
#      highest_count <- pi_res[pos[1], pos[2]]
#      pi_res[pos[1], pos[2]] <- 0  #remove highest results
#      new_max_results <- unlist(unname(which(pi_res == max(pi_res), arr.ind = TRUE)))
      
#      if(length(new_max_results) == 4){
#         set.seed(901)
#         randnum <- round(runif(1, 1, nrow(new_max_results))) #generate a random number between 1 and nrows of results
#         second_highest_size <- new_max_results[randnum,] + 15
#         pos <- c(new_max_results[1,1], new_max_results[1,2])
#         second_highest_count <- pi_res[pos[1], pos[2]]
         
#      } else if(length(new_max_results == 2)) {
#         second_highest_size <- c(new_max_results) + 15
#         pos <- c(new_max_results[1,1], new_max_results[1,2])
#         second_highest_count <- pi_res[pos[1], pos[2]]
#      } else {
#         second_highest_size <- NA
#         second_highest_count <- NA
#      } 
      
#   } else if(length(max_results) == 4){ #if two max results, choose from max pi_results at random
#      set.seed(123)
#      randnum <- round(runif(1, 1, nrow(max_results)))
#      highest_size <- max_results[randnum,] + 15
#      pos <- c(max_results[1,1], max_results[1,2]) #remove from the pi_results
#      highest_count <- pi_res[pos[1], pos[2]]
#      pi_res[pos[1], pos[2]] <- 0
#      new_max_results <- unlist(unname(which(pi_res == max(pi_res), arr.ind = TRUE))) #again choose the highest result
      
#      if(length(new_max_results) == 2){
#         second_highest_size <- c(new_max_results) + 15
#         pos <- c(new_max_results[1,1], new_max_results[1,2])
#         second_highest_count <- pi_res[pos[1], pos[2]]
#      } else if(length(new_max_results == 4)){
#         set.seed(375)
#         randnum <- round(runif(1,1, nrow(new_max_results)))
#         second_highest_size <- max_results[randnum,] + 15
#         pos <- c(new_max_results[1,1], new_max_results[1,2])
#         second_highest_count <- pi_res[pos[1], pos[2]]         
#      } else {
#         second_highest_size <- NA
#         second_highest_count <- NA
#      }
      
#   } else { #if more than 2 max results
#      set.seed(923)
#      randnum <- round(runif(1,1,nrow(max_results)))
#      highest_size <- max_results[randnum,] + 15
#      pos <- c(max_results[1,1], max_results[1,2])
#      highest_count <- pi_res[pos[1], pos[2]]
#      pi_res[pos[1], pos[2]] <- 0
#      new_max_results <- unlist(unname(which(pi_res == max(pi_res), arr.ind = TRUE)))
#      set.seed(87)
#      randnum <- round(runif(1,1,nrow(new_max_results)))
#      second_highest_size <- new_max_results[randnum,] + 15
#      pos <- c(new_max_results[1,1], new_max_results[1,2])
#      second_highest_count <- pi_res[pos[1], pos[2]]
#   }
#   max_pi_df <- data.frame(highest_pi_row = highest_size[1], highest_pi_col = highest_size[2], highest_pi_count = highest_count, 
#                           second_pi_row = second_highest_size[1], second_pi_col = second_highest_size[2], second_pi_count = second_highest_count)
#   return(max_pi_df)
#}


get_max_pi_heat <- function(pi_res){
   pi_table <- as.matrix(pi_res[[1]])
   
   max_results <- as.matrix(unlist(unname(which(pi_res[[1]] == max(pi_res[[1]]), arr.ind = TRUE))))
   
   if(nrow(max_results) == 0){
      highest_pi_count <- 0
      highest_pi_row <- NA
      highest_pi_col <- NA
   } else if(nrow(max_results) == 1){
      highest_pi_row <- max_results[1,1]
      highest_pi_col <- max_results[1,2]
      highest_pi_count <- pi_table[highest_pi_row, highest_pi_col]
   } else {
      counts <- pi_table[max_results]
      highest_pi_count <- mean(counts)
      highest_pi_row <- mean(max_results[,1])
      highest_pi_col <- mean(max_results[,2])
      
   }
   
   highest_pi_row <- highest_pi_row + 15
   highest_pi_col <- highest_pi_col + 15
   
   max_pi_df <- data.frame(highest_pi_row = highest_pi_row, highest_pi_col = highest_pi_col, highest_pi_count = highest_pi_count)
   return(max_pi_df)
   
   
   
   
   
}


















