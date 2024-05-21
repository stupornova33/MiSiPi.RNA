#' write_proper_overhangs function
#' if write_fastas == TRUE, outputs read pairs with 3' 2nt overhang as a fasta file
#' @param dt1 a data frame of reads
#' @param dt2 a data frame of reads
#' @param wkdir a string specifying the directory to write the file to
#' @param prefix a string specifying the prefix of the file name
#' @param overlaps Specifies whether file types for plots are png or pdf. Default is pdf.
#' @param suffix a string specifying the suffix of the file name
#' @return nothing

#' @export

write_proper_overhangs <- function(dt1, dt2, wkdir, prefix, overlaps, suffix){

  proper_overhangs <- overlaps %>% dplyr::filter(r2_start - r1_start == 2 & r2_end - r1_end == 2)
  rname <- dt1$rname[1]

  uniq_forward <- dt1 %>% dplyr::distinct()
  uniq_reverse <- dt2 %>% dplyr::distinct()
  tmp <- rbind(uniq_forward, uniq_reverse)
  uniq_overhangs <- proper_overhangs %>% dplyr::distinct()
  rreads <- data.frame()
  freads <- data.frame()
  for(i in 1:nrow(uniq_overhangs)){
    tmp_r <- tmp[which(tmp$start == uniq_overhangs$r1_start[i] & tmp$end == uniq_overhangs$r1_end[i]), ] %>%
      dplyr::distinct(start, end, .keep_all = TRUE)
    tmp_f <- tmp[which(tmp$start == uniq_overhangs$r2_start[i] & tmp$end == uniq_overhangs$r2_end[i]), ] %>%
      dplyr::distinct(start, end, .keep_all = TRUE)

    rreads <- rbind(rreads, tmp_r)
    freads <- rbind(freads, tmp_f)
  }


  proper_overhangs <- NULL

  #if(suffix != "_hairpin"){
  rreads <- rreads %>% dplyr::rename("r1_start" = "start", "r1_end" = "end", "r1_seq" = seq) %>% dplyr::select(-c(width, first,rname))
  freads <- freads %>% dplyr::rename("r2_start" = "start", "r2_end" = "end", "r2_seq" = seq) %>% dplyr::select(-c(width, first,rname))

  #}

  if(nrow(rreads) > 0 & nrow(freads) > 0){
    paired_seqs <- cbind(rreads, freads)
    paired_seqs <- paired_seqs %>%
    dplyr::mutate(read1_seq = paste0(rname, ":", paired_seqs$r1_start, "-", paired_seqs$r1_end, " ", paired_seqs$r1_seq),
                  read2_seq = paste0(rname, ":", paired_seqs$r1_start, "-", paired_seqs$r1_end, " " , paired_seqs$r2_seq))

   fastas <- paired_seqs %>% dplyr::select(c(read1_seq, read2_seq)) %>%
     dplyr::transmute(col1 = paste0(read1_seq, ",", read2_seq)) %>%
     tidyr::separate_rows(col1, sep = ",")

   fastas <- stringi::stri_split_regex(fastas$col1, " ")


  } else {
    if(nrow(rreads) > 0){
      paired_seqs <- rreads
      paired_seqs <- paired_seqs %>%
        dplyr::mutate(read1_seq = paste0(">", rname, ":", rreads$r1_start, "-", rreads$r1_end, " ", paired_seqs$r1_seq)) %>%
        dplyr::ungroup()

      fastas <- paired_seqs %>% dplyr::select(c(read1_seq)) %>%
        dplyr::transmute(col1 = paste0(read1_seq, " ")) #%>%
        #tidyr::separate_rows(col1, sep = " ")

      fastas <- stringi::stri_split_regex(fastas$col1, " ")

    } else {
      paired_seqs <- freads
      paired_seqs <- paired_seqs %>%
        dplyr::mutate(read2_seq = paste0(">",rname, ":", r2_start, "-", r2_end, " ", r2_seq)) %>%
        dplyr::ungroup()

      fastas <- paired_seqs %>% dplyr::select(c(read2_seq)) %>%
        dplyr::transmute(col1 = paste0(read2_seq, " ")) #%>%
        #tidyr::separate_rows(col1, sep = " ")

      fastas <- stringi::stri_split_regex(fastas$col1, " ")
    }
  }

  rreads <- NULL
  freads <- NULL



  write.table(unlist(fastas), file = paste0(prefix, "_dicer", suffix, ".fa"), sep = " ", quote = FALSE, row.names = FALSE, col.names = FALSE)
  paired_seqs <- NULL
  fastas <- NULL


}
