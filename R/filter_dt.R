#' function to filter dt based on passed in min and max widths
#' selects only reads beginning with T for piRNAs
#' takes two a df and a min and max width
#' returns a filtered dt
#' 
#' @param chrom a chrom obj
#' @return filtered_df

#' @export



filter_dt <- function(chrom){
   filter_dt <- data.table::setDT(makeBamDF(chrom)) %>%
      base::subset(width <= 32 & width >= 18) %>% 
      dplyr::filter(first == "T") %>%
      dplyr::mutate(start = pos, end = pos + width - 1) %>%
      dplyr::select(-c(pos, first, seq))
   return(filter_dt)
}
