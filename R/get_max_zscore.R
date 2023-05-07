#' function to find max zscore
#' takes a table with one column of zscores and one column of overlaps
#' returns max_zscore
#' 
#' @param Z_score a column of a data table containing zscores
#' @param Overlap a column of a data table containing overlap values
#' @return max_z



#' @export


get_max_zscore <- function(Z_score, Overlap){
   if(!sum(Z_score) == "NaN"){
      max <- which(Z_score == max(Z_score))
      if(length(max) > 1){
         overlap_df <- Overlap[max]
         maxz_overlap <- mean(overlap_df)
         z_df <- Z_score[max]
         max_zscore <- mean(Z_score)
      } else {
         max_zscore <- Z_score[max]
         maxz_overlap <- Overlap[max]
      }
   } else {
      max_zscore <- NA
      maxz_overlap <- NA
   }
   
   return(list(maxz_overlap, max_zscore))
   
}
