#' function to find max_si heatmap cell values
#' takes a matrix created from siRNA function
#' finds max row, column, and count of highest cell
#' returns highest_si_row, highest_si_col, second_si_row, second_si_col, highest_si_count, second_si_count as a data frame
#' 
#' @param si_res a data table
#' @return max_si_df



#' @export

#get_max_si_heat<- function(si_res){
#  max_results <- unlist(unname(which(si_res == max(si_res), arr.ind = TRUE)))
#  if(length(max_results) == 0){
#      highest_size <- NA
#      second_highest_size <- NA
#      highest_count <- NA
#      second_highest_count <- NA
      
#    } else if (length(max_results) == 2) {
#         highest_size <- max_results + 15
#         pos <- c(max_results[1,1], max_results[1,2])
#         highest_count <- si_res[pos[1], pos[2]]
#         si_res[pos[1], pos[2]] <- 0  #remove highest results
#         new_max_results <- unlist(unname(which(si_res == max(si_res), arr.ind = TRUE)))
         
#         if(length(new_max_results) == 4){
#            set.seed(901)
#            randnum <- round(runif(1, 1, nrow(new_max_results))) #generate a random number between 1 and nrows of results
#            second_highest_size <- new_max_results[randnum,] + 15
#            pos <- c(new_max_results[1,1], new_max_results[1,2])
#            second_highest_count <- si_res[pos[1], pos[2]]
            
#          } else if(length(new_max_results == 2)) {
#               second_highest_size <- c(new_max_results) + 15
#               pos <- c(new_max_results[1,1], new_max_results[1,2])
#               second_highest_count <- si_res[pos[1], pos[2]]
#          } else {
#               second_highest_size <- NA
#               second_highest_count <- NA
#          } 
         
#     } else if(length(max_results) == 4){ #if two max results, choose from max si_results at random
#          set.seed(123)
#          randnum <- round(runif(1, 1, nrow(max_results)))
#          highest_size <- max_results[randnum,] + 15
#          pos <- c(max_results[1,1], max_results[1,2]) #remove from the si_results
#          highest_count <- si_res[pos[1], pos[2]]
#          si_res[pos[1], pos[2]] <- 0
#          new_max_results <- unlist(unname(which(si_res == max(si_res), arr.ind = TRUE))) #again choose the highest result
          
#          if(length(new_max_results) == 2){
#             second_highest_size <- c(new_max_results) + 15
#             pos <- c(new_max_results[1,1], new_max_results[1,2])
#             second_highest_count <- si_res[pos[1], pos[2]]
#          } else if(length(new_max_results == 4)){
#             set.seed(375)
#             randnum <- round(runif(1,1, nrow(new_max_results)))
#             second_highest_size <- max_results[randnum,] + 15
#             pos <- c(new_max_results[1,1], new_max_results[1,2])
#             second_highest_count <- si_res[pos[1], pos[2]]         
#          } else {
#             second_highest_size <- NA
#             second_highest_count <- NA
#          }
          
#      } else if(length(max_results) == 6){ #if more than 2 max results
#         set.seed(923)
#         randnum <- round(runif(1,1,nrow(max_results)))
#         highest_size <- max_results[randnum,] + 15
#         pos <- c(max_results[1,1], max_results[1,2])
#         highest_count <- si_res[pos[1], pos[2]]
#         si_res[pos[1], pos[2]] <- 0
#         new_max_results <- unlist(unname(which(si_res == max(si_res), arr.ind = TRUE)))
#         set.seed(87)
#         randnum <- round(runif(1,1,nrow(new_max_results)))
#         second_highest_size <- new_max_results[randnum,] + 15
#         pos <- c(new_max_results[1,1], new_max_results[1,2])
#         second_highest_count <- si_res[pos[1], pos[2]]
#   } else{
#      highest_size <- c(NA, NA)
#      highest_count <- NA
#      second_highest_size <- c(NA, NA)
#      second_highest_count <- NA
      
#   }
#   max_si_df <- data.frame(highest_si_row = highest_size[1], highest_si_col = highest_size[2], highest_si_count = highest_count, 
#                           second_si_row = second_highest_size[1], second_si_col = second_highest_size[2], second_highest_si_count = second_highest_count)
#   return(max_si_df)
#}

get_max_si_heat<- function(si_res){
   
   max_results <- as.matrix(unlist(unname(which(si_res[[1]] == max(si_res[[1]]), arr.ind = TRUE))))
   si_table <- as.matrix(si_res[[1]])
   if(nrow(max_results) == 0){
      highest_si_count <- 0
      highest_si_row <- NA
      highest_si_col <- NA
   } else if(nrow(max_results) == 1){
      highest_si_row <- max_results[1,1]
      highest_si_col <- max_results[1,2]
      highest_si_count <- si_table[highest_si_row, highest_si_col]
   } else {
      counts <- si_table[max_results]
      highest_si_count <- mean(counts)
      #max_results <- max_results + 15
      highest_si_row <- mean(max_results[,1])
      highest_si_col <- mean(max_results[,2])
      
   }
   highest_si_row <- highest_si_row + 15
   highest_si_col <- highest_si_col + 15
   
   max_si_df <- data.frame(highest_si_row = highest_si_row, highest_si_col = highest_si_col, highest_si_count = highest_si_count)
   return(max_si_df)
      
}
   

