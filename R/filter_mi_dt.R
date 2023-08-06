#' function to filter data table for short hairpin function
#' makes chrom obj into bam df and filters reads between 18-25 nt
#' summarizes count of grouped reads
#' returns filter dt
#'
#' @param chrom a chrom object
#' @param chrom_name a string
#' @return filter_dt

#' @export


filter_mi_dt <- function(chrom, chrom_name){
   width <- pos <- start <- end <- NULL
  filter_dt <- data.table::setDT(makeBamDF(chrom)) %>%
    base::subset(width <= 26 & width >= 18) %>%
    dplyr::mutate(start = pos, end = pos + width) %>%
    dplyr::group_by_at(dplyr::vars(start, end)) %>%
    dplyr::summarize(count = dplyr::n())
  filter_dt <- filter_dt %>% dplyr::mutate(rname = chrom_name)
  return(filter_dt)
}
