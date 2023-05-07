#' extracts read data and plots counts by read size 
#' for multiple loci in a bed file
#' @param bam_file a a string
#' @param logfile a string
#' @param chrom_name a list of chromosome names
#' @param reg_start a list of integers
#' @param reg_stop a list of integers
#' @return a heatmap plot
#' @export


heatmap_dendrogram <- function(bam_file, logfile, chrom_name, reg_start, reg_stop){
  
#plot heatmap dendrogram for read sizes of all loci in input bed file
  bam_obj <- OpenBamFile(bam_file, logfile)
  bam_header <- Rsamtools::scanBamHeader(bam_obj)
  chr_name <- names(bam_header[['targets']])
  chr_length <- unname(bam_header[['targets']])
  bam_header <- NULL
  

  rng_df <- data.frame(chrom_name = c(unlist(vars[1])), start = c(unlist(vars[2])), end = c(unlist(vars[3])))
  
  get_reads <- function(i) {
    chromM <- getChrMinus(bam_obj, chrom_name = c(rng_df[i,1]), reg_start = c(rng_df[i,2]), reg_stop = c(rng_df[i,3]))
    chromP <- getChrPlus(bam_obj, chrom_name = c(rng_df[i,1]), reg_start = c(rng_df[i,2]), reg_stop = c(rng_df[i,3]))
    
    forward_dt <- data.table::setDT(makeBamDF(chromP)) %>%
      subset(width <= 32 & width >= 15) %>% 
      dplyr::mutate(start = pos, end = pos + width - 1) %>% dplyr::distinct() %>%
      dplyr::select(-c(pos, first, seq)) %>% dplyr::mutate(name = paste0(rname, ":", rng_df[i,2], "-", rng_df[i,3])) %>%
      dplyr::select(-c(rname))
    
    
    reverse_dt <- data.table::setDT(makeBamDF(chromM)) %>%
      subset(width <= 32 & width >= 15) %>%
      dplyr::mutate(start = pos, end = pos + width - 1) %>% dplyr::distinct() %>%
      dplyr::select(-c(pos, first,seq)) %>% dplyr::mutate(name = paste0(rname, ":", rng_df[i,2], "-", rng_df[i,3])) %>% 
      dplyr::select(-c(rname))
    
    all_dt <- rbind(forward_dt, reverse_dt)
    all_dt <- all_dt %>% dplyr::count(width, name)
    
    
    empty_tbl <- data.frame(width = c(seq(15, 32)), n = 0)
    merged <- merge(all_dt, empty_tbl, by = c('width'), all = TRUE) %>% dplyr::select(-c(n.y))
    
    idx <- which(is.na(merged$name) == "TRUE")
    merged$name[idx] <- paste0(rng_df[i,1], ":", rng_df[i,2], "-", rng_df[i,3])
    merged[is.na(merged)] <- 0
    return(merged)  
  }
  
  #dt <- mapply(get_reads, bam_obj, rng_df$chr_name, rng_df$start, rng_df$stop)
  dt <- lapply(seq(nrow(rng_df)), get_reads)
  combined_dt <- dplyr::bind_rows(dt) %>% dplyr::rename('locus' = name)

 #combined_dt <- combined_dt %>% dplyr::mutate(across(where(is.numeric), scale))

  mat <- matrix(combined_dt$n.x, ncol = 18)
  colnames(mat) <- as.character(unique(combined_dt$width))

    
  heatmap <- pheatmap::pheatmap(mat, cluster_cols = FALSE, show_rownames = F,  show_colnames = T)
  
  
  return(heatmap)
  
}
