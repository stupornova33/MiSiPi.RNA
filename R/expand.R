#' expand
#' for subsetting a string 
#' @param string a string
#' @return vienna
#' @export


expand <- function(string){
   vienna <- unlist(strsplit(string, ''))
   return(vienna)
}
