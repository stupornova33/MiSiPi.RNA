# write_proper_overhangs function
# if write_fastas == TRUE, outputs read pairs with 3' 2nt overhang as a fasta file
# @param df1 a data frame of reads
# @param df2 a data frame of reads
# @param wkdir a string specifying the directory to write the file to
# @param prefix a string specifying the prefix of the file name
# @param overlaps
# @param suffix a string specifying the suffix of the file name
# @return nothing

.write_proper_overhangs <- function(df1, df2, wkdir, prefix, overlaps, suffix) {
  proper_overhangs <- overlaps %>%
    dplyr::filter(p5_overhang == 2 & p3_overhang == 2) %>%
    dplyr::select(-c(r1_dupes, r2_dupes))

  rname <- df1$rname[1]

  # We might be able to get rid of this and pass in the summarized dts
  # In .siRNA, at least they are available
  # TODO check other call locations
  uniq_forward <- dplyr::distinct(df1)
  uniq_reverse <- dplyr::distinct(df2)

  r1_matches <- proper_overhangs %>%
    dplyr::select(r1_start, r1_end) %>%
    dplyr::distinct() %>%
    dplyr::rename(start = r1_start, end = r1_end)
  
  r2_matches <- proper_overhangs %>%
    dplyr::select(r2_start, r2_end) %>%
    dplyr::distinct() %>%
    dplyr::rename(start = r2_start, end = r2_end)
  
  # Inner joins to find all matching rows
  freads <- dplyr::inner_join(r1_matches, uniq_forward, by = c("start", "end")) %>%
    dplyr::distinct(start, end, .keep_all = TRUE) %>%
    dplyr::relocate(rname)
  
  rreads <- dplyr::inner_join(r2_matches, uniq_reverse, by = c("start", "end")) %>%
    dplyr::distinct(start, end, .keep_all = TRUE) %>%
    dplyr::relocate(rname)
  
  proper_overhangs <- NULL

  # if(suffix != "_hairpin"){
  rreads <- rreads %>%
    dplyr::rename("r1_start" = "start", "r1_end" = "end", "r1_seq" = seq) %>%
    dplyr::select(-c(width, first, rname))

  freads <- freads %>%
    dplyr::rename("r2_start" = "start", "r2_end" = "end", "r2_seq" = seq) %>%
    dplyr::select(-c(width, first, rname))

  # }

  if (nrow(freads) == 0) {
    freads <- NULL
  }
  if (nrow(rreads) == 0) {
    rreads <- NULL
  }
  
  fastas <- .compile_fastas(freads, rreads, rname)

  rreads <- NULL
  freads <- NULL

  filename <- paste0(prefix, "_dicer", suffix, ".fa")

  write.table(unlist(fastas),
    file = file.path(wkdir, filename),
    sep = " ",
    quote = FALSE,
    row.names = FALSE,
    col.names = FALSE
  )

  fastas <- NULL
}


.compile_fastas <- function(freads, rreads, rname) {
  if (is.null(freads) & is.null(rreads)) return(NULL)
  if (is.null(freads)) { # Only rreads present
    fastas <- rreads %>%
      dplyr::mutate(read1_seq = paste0(">", rname, ":", r1_start, "-", r1_end, " ", r1_seq)) %>%
      dplyr::select(read1_seq) %>%
      dplyr::transmute(col1 = paste0(read1_seq, " "))
    
  } else if (is.null(rreads)) { # Only freads present
    fastas <- freads %>%
      dplyr::mutate(read2_seq = paste0(">", rname, ":", r2_start, "-", r2_end, " ", r2_seq)) %>%
      dplyr::select(read2_seq) %>%
      dplyr::transmute(col1 = paste0(read2_seq, " "))

  } else { # Both rreads and freads present
    fastas <- dplyr::bind_cols(rreads, freads) %>%
      dplyr::mutate(
        read1_seq = paste0(">", rname, ":", r1_start, "-", r1_end, " ", r1_seq),
        read2_seq = paste0(">", rname, ":", r2_start, "-", r2_end, " ", r2_seq)
      ) %>%
      dplyr::select(c(read1_seq, read2_seq)) %>%
      dplyr::transmute(col1 = paste0(read1_seq, ",", read2_seq)) %>%
      tidyr::separate_rows(col1, sep = ",")
  }
  
  fastas <- stringi::stri_split_regex(fastas$col1, " ")
  
  return(fastas)
}
