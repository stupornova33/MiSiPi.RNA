#' @title exists_not_empty
#'
#' @description Checks to see if a given file exists and is not empty
#'
#' @param file_name Filename to be tested. If path is not given, working directory will be the path
#'
#' @return TRUE or FALSE
#'
#' @examples
#' if (exists_not_empty(file_name)) { print("True") }
#' @export

exists_not_empty <- function(file_name) {
  if (file.exists(file_name) & file.size(file_name) > 0) {
    return(TRUE)
  } else {
    return(FALSE)
  }
}
