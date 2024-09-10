#' Make read size distribution table
#' @param dat a data frame containing read sizes and counts
#' @param prefix a string
#' @param wkdir a string
#' @param strand a string, "+" or "-". Default is NULL.
#' @param type a string
#' @return nothing
#' @export
#'

output_readsize_dist <- function(dat, prefix, wkdir, strand, type){

  dat <- as.data.frame(dat)

  empty_tbl <- data.frame(width = c(seq(18,32,1)), count = 0)
  sizes <- merge(dat, empty_tbl, by = c("width"), all = TRUE) %>% dplyr::select(c(width, count.x)) %>% dplyr::rename("count" = "count.x") %>%
    dplyr::select(c(count))

  sizes[is.na(sizes)] <- 0
  sizes <- t(sizes)

  cols <- c(seq(18,32,1))

  colnames(sizes) <- cols

  sizes <- as.data.frame(sizes)

  sizes$locus <- prefix
  sizes <- sizes %>% dplyr::select(c(locus, "18", "19", "20", "21", "22", "23","24","25","26","27","28","29","30", "31", "32"))

  if(!is.null(strand)){
    name <- paste0(wkdir, type, "_", strand, "_read_size_distributions.txt")
    suppressWarnings(
      if(!file.exists(name)){
        utils::write.table(sizes, file = paste0(wkdir, type, "_", strand, "_read_size_distributions.txt"), sep = "\t", quote = FALSE, append = T, col.names = T, na = "NA", row.names = F)
      } else {
        utils::write.table(sizes, file = paste0(wkdir, type, "_", strand, "_read_size_distributions.txt"), quote = FALSE, sep = "\t", col.names = F, append = TRUE, na = "NA", row.names = F)
      })
  } else {
      name <- paste0(wkdir,type, "_read_size_distributions.txt")
      suppressWarnings(
        if(!file.exists(name)){
          utils::write.table(sizes, file = paste0(wkdir, type, "_read_size_distributions.txt"), sep = "\t", quote = FALSE, append = T, col.names = T, na = "NA", row.names = F)
        } else {
          utils::write.table(sizes, file = paste0(wkdir, type, "_read_size_distributions.txt"), quote = FALSE, sep = "\t", col.names = F, append = TRUE, na = "NA", row.names = F)
        })
    }

}
