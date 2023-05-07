#' miRNA structure processing function
#' @param converted a string
#' @param vienna a string

#' @return f_df
#' @export


process_miRNA_struct <- function(vienna, converted){

   #returns left loop coords before right loop coords
   # if there is third loop?
   orig_loop_pos <- unlist(stringi::stri_locate_all_regex(vienna, '\\(\\.*\\)'))
   vienna_vec <- unlist(strsplit(vienna, ''))
   len <- length(vienna_vec)
   num_loops <- length(orig_loop_pos)/2
   if(num_loops == 1){
      #call one loop script
      res_list <- one_loop(vienna_vec, orig_loop_pos, converted)
   } else if(num_loops == 2){
      res_list <- two_loops(vienna_vec, orig_loop_pos, converted)
      
   } else if(num_loops == 3) { #if 3 loops
      #call three loop script
      res_list <- three_loops(vienna_vec, orig_loop_pos, converted)
    } else {  #num loops 4 or more
      res_list <- NA
    }
return(res_list)
}
