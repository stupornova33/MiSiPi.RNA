# write_proper_pairs function
# if write_fastas == TRUE, outputs read pairs with 3' 2nt overhang for miRNA and siRNA
#     or 10nt overlap for piRNA as a fasta file
# @param df1 a data frame of reads
#   - When called from .siRNA() or .piRNA(), this will be the sense strand of reads
#   - When called from .miRNA(), this will be whatever strand of reads is currently being processed
# @param df2 a data frame of reads
#   - When called from .siRNA() or .piRNA(), this will be the antisense strand of reads
#   - When called from .miRNA(), this will be NULL as only one strand is being processed
# @param wkdir a string specifying the directory to write the file to
# @param prefix a string specifying the prefix of the file name
# @param overlaps
# @param suffix a string specifying the suffix of the file name
# @param calling_func The function that called this function
# @param strand When called from miRNA, this will specify the strand
# @return nothing

.write_proper_pairs <- function(df1, df2 = NULL, wkdir, prefix, overlaps, suffix, calling_func = c("mi", "si", "pi"), strand = NULL) {
  
  calling_func <- match.arg(calling_func)

  # In the case of miRNA and siRNA these would be overlapping reads with a proper dicer signature
  # In the case of piRNA this would be overlapping reads with the proper 10nt overlap
  proper_pairs <- .get_proper_pairs(overlaps, calling_func)
  
  if (nrow(proper_pairs) == 0) return()
  
  rname <- df1$rname[1]

  uniq_forward <- dplyr::distinct(df1) %>%
    dplyr::select(rname, start, end, seq)
  
  if (is.null(df2)) {
    uniq_reverse <- uniq_forward
  } else {
    uniq_reverse <- dplyr::distinct(df2) %>%
      dplyr::select(rname, start, end, seq)
  }
  
  # Define strands for fasta output
  # They will be the same for miRNA calls
  if (calling_func == "mi") {
    strand_1 <- strand
    strand_2 <- strand
  } else {
    strand_1 <- "+"
    strand_2 <- "-"
  }
  
  df1_key <- uniq_forward %>%
    dplyr::select(
      r1_start = start,
      r1_end = end,
      r1_seq = seq
    )
  
  uniq_forward <- NULL
  
  df2_key <- uniq_reverse %>%
    dplyr::select(
      r2_start = start,
      r2_end = end,
      r2_seq = seq
    )
  
  uniq_reverse <- NULL
  
  sequence_pairs <- proper_pairs %>%
    dplyr::mutate(rname = rname) %>%
    dplyr::left_join(df1_key, by = c("r1_start", "r1_end")) %>%
    dplyr::left_join(df2_key, by = c("r2_start", "r2_end"))
  
  if (calling_func == "mi") {
    sequence_pairs <- sequence_pairs %>%
      na.omit()
  }
  
  sequence_pairs <- sequence_pairs %>%
    dplyr::mutate(pair = dplyr::row_number()) %>%
    dplyr::group_by(r1_start, r1_end, r2_start, r2_end) %>%
    dplyr::mutate(group = dplyr::cur_group_id()) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(
      r1_header = glue::glue("{rname}:{r1_start}-{r1_end}|strand={strand_1}|group={group}|pair={pair}"),
      r2_header = glue::glue("{rname}:{r2_start}-{r2_end}|strand={strand_2}|group={group}|pair={pair}")
    ) %>%
    dplyr::select(r1_header, r1_seq, r2_header, r2_seq)
  
  proper_pairs <- NULL
  df1_key <- NULL
  df2_key <- NULL
  
  if (nrow(sequence_pairs) == 0) return()
  
  sequence_pairs <- sequence_pairs %>%
    dplyr::mutate(.row_id = dplyr::row_number())
  
  df1_seqs <- sequence_pairs %>%
    dplyr::select(header = r1_header, seq = r1_seq, .row_id) %>%
    dplyr::mutate(source = "df1")
  
  df2_seqs <- sequence_pairs %>%
    dplyr::select(header = r2_header, seq = r2_seq, .row_id) %>%
    dplyr::mutate(source = "df2")
  
  sequence_pairs <- NULL
  
  fasta_df <- dplyr::bind_rows(df1_seqs, df2_seqs) %>%
    dplyr::arrange(.row_id, source) %>%
    dplyr::select(-c(.row_id, source))
  
  df1_seqs <- NULL
  df2_seqs <- NULL
  
  
  if (calling_func == "mi" | calling_func == "si") {
    filename <- paste0(prefix, "_dicer", suffix, ".fa")
  } else {
    filename <- paste0(prefix, "_piRNA_pairs.fa")
  }
  
  write.table(
    fasta_df,
    file = file.path(wkdir, filename),
    sep = " ",
    quote = FALSE,
    row.names = FALSE,
    col.names = FALSE
  )

  fasta_df <- NULL
}

.get_proper_pairs <- function(overlaps, method = c("mi", "si", "pi")) {
  method <- match.arg(method)
  
  if (method == "mi" | method == "si") {
    if (method == "mi") {
      overlaps <- overlaps %>%
        dplyr::mutate(
          p5_overhang = r1_start - r2_start,
          p3_overhang = r1_end - r2_end
        )
    }
    
    proper_pairs <- overlaps %>%
      dplyr::filter(p5_overhang == 2 & p3_overhang == 2) %>%
      dplyr::select(-c(r1_dupes, r2_dupes, p5_overhang, p3_overhang))
    
  } else {
    # For piRNA pairs
    proper_pairs <- overlaps %>%
      dplyr::mutate(overlap = dplyr::case_when(
        r1_start > r2_start ~ (r2_end - r1_start), r1_start < r2_start ~ (r1_end - r2_start),
        (r1_start >= r2_start & r1_end <= r2_end) ~ (r1_end - r1_start + 1),
        (r2_start >= r1_start & r2_end <= r1_end) ~ (r2_end - r2_start + 1)
      )) %>%
      dplyr::filter(overlap == 10) %>%
      dplyr::select(-c(r1_dupes, r2_dupes, overlap))
  }
  
  return(proper_pairs)
}
