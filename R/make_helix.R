#' function to convert vienna to helix
#' @param fold_list
#' @return list

#' @export


make_helix <- function(fold_list){
   vienna <- unlist(unname(fold_list))
   helix <- R4RNA::viennaToHelix(unlist(unname(fold_list)))
   return(list(vienna, helix))
}
