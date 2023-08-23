#' run the run_all function
#' @param vars a list
#' @importFrom Rcpp sourceCpp
#' @return plots

#' @export

run_all_function <- function(vars){

`%>%` <- magrittr::`%>%`

mapply(new_run_all, vars[[1]], vars[[2]], vars[[3]], vars[[5]], vars[[4]], vars[[10]], vars[[9]], vars[[8]],
          vars[[13]], vars[[12]], vars[[6]], vars[[7]], vars[[14]], vars[[15]], vars[[11]])

}

