#' function to find max z_score
#' takes a table with one column of z_scores and one column of overlaps
#' returns max_zscore
#' 
#' @param z_score a column of a data table containing z_scores
#' @param Overlap a column of a data table containing overlap values
#' @return max_z

#' @export

get_max_zscore <- function(z_score, Overlap){
   if(!sum(z_score) == "NaN"){
      max <- which(z_score == max(z_score))
      if(length(max) > 1){
         overlap_df <- Overlap[max]
         maxz_overlap <- mean(overlap_df)
         z_df <- z_score[max]
         max_zscore <- mean(z_score)
      } else {
         max_zscore <- z_score[max]
         maxz_overlap <- Overlap[max]
      }
   } else {
      max_zscore <- NA
      maxz_overlap <- NA
   }
   
   return(list(maxz_overlap, max_zscore))
   
}
