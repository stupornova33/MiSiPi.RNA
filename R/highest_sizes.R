#' function to find the read size of hte highest abundance
#' takes a table of read sizes and counts
#' returns highest and second highest read sizes
#' 
#' @param read_dist a table 
#' @return a list



#' @export

#highest_sizes <- function(read_dist){
#   max <- which(read_dist$count == max(read_dist$count))
#   if(length(max) > 2){
#      return(list(NA, NA))
#   } else if(length(max) == 2){
#     highest_size <- read_dist$width[max[1]]
#     second_highest <- read_dist$width[max[2]]
#   } else {
#   highest_size <- read_dist$width[max]
#   read_dist <- read_dist[-c(max),]
#   max2 <- which(read_dist$count == max(read_dist$count))
#   second_highest <- read_dist$width[max]
#   }
#   return(list(highest_size, second_highest))
#}


#try to take the average 

highest_sizes <- function(read_dist){
   max <- which(read_dist$count == max(read_dist$count))
   if(length(max) > 1){
      #take the average of the first two results
      highest_size <- mean(read_dist$width[max[1]], read_dist$width[max[2]])
   } else {
      highest_size <- read_dist$width[max]
   }
   
   return(highest_size)
}