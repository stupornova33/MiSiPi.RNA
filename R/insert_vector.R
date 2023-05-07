#' insert vector
#' inserts characters into string vector at desired position
#' @param input a string vector
#' @param position an integer
#' @param values a character
#' @return res a vector 
#' @export


#insert vector
insert_vector <- function(input, position, values) {
   options(warn=-1)
   #Create result vector
   res <- numeric(length(input) + length(values))
   res[position] <- values
   #Create an index of input vector
   inds1 <- seq_along(res)
   inds2 <- !(inds1 %in% position)
   
   #Insert input vector
   res[inds2] <- input
   return(res)
}