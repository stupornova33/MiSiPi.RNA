#' run the run_all function
#' @param vars a list
#' @importFrom Rcpp sourceCpp
#' @return plots

#' @export

run_all_function <- function(vars){

                   #chrom name start       stop       chrom      length    input       bed file   genome    min read count
mapply(new_run_all, vars[[1]], vars[[2]], vars[[3]], vars[[5]], vars[[4]], vars[[10]], vars[[11]],vars[[9]], vars[[8]],
       # si pal       pi pal      plot output RNA fold   annotatebed weight      bed file   annot file
          vars[[13]], vars[[12]], vars[[6]], vars[[7]], vars[[14]], vars[[15]], vars[[16]])

}


#' run the run_all function
#' @param vars a list
#' @importFrom Rcpp sourceCpp
#' @return plots

#' @export

misipi_rna <- function(vars){

`%>%` <- magrittr::`%>%`
                   #chrom name start       stop       chrom      length    input       bed file   genome    min read count
mapply(new_run_all, vars[[1]], vars[[2]], vars[[3]], vars[[5]], vars[[4]], vars[[10]], vars[[11]],vars[[9]], vars[[8]],
       # si pal       pi pal      plot output RNA fold   annotatebed weight      bed file   annot file
          vars[[13]], vars[[12]], vars[[6]], vars[[7]], vars[[14]], vars[[15]], vars[[16]])

}

