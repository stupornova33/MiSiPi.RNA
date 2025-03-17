#' Wrapper function that calls core functions
#' @param vars a list of parameters created by set_vars
#' @param method a string indicating which processing method should be used
#' @return plots

#' @export

misipi_rna <- function(vars, method = c("all", "miRNA", "piRNA", "siRNA")) {
  method <- match.arg(method)
  
  if (method == "all") {
    run_all(vars)
  } else if (method == "miRNA") {
    miRNA(vars)
  } else if (method == "piRNA") {
    piRNA(vars)
  } else if (method == "siRNA") {
    siRNA(vars)
  } else {
    stop(paste("`method`:", method, "is not valid."))
  }
}
