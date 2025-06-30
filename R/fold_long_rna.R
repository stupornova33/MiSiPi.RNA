# function to run RNAfold
# processes output of RNA fold to get MFE and vien struct
# returns list of values for each region
# @param chrom_name a string
# @param start a whole number
# @param stop a whole number
# @param converted a vector containing a sequence
# @param path_to_RNAfold a string
# @param wkdir a string
# @return list

.fold_long_rna <- function(chrom_name, start, stop, converted, path_to_RNAfold, wkdir) {
  
  converted_df <- data.frame(converted)
  
  region_string <- paste0(chrom_name, "-", start - 1, "_", stop - 1)
  fasta_header <- paste0(">", region_string)
  
  colnames(converted_df) <- fasta_header
  
  write.table(converted_df,
    file = file.path(wkdir, "converted.fasta"),
    sep = "\n",
    append = FALSE,
    row.names = FALSE,
    quote = FALSE
  )

  syscheck <- Sys.info()["sysname"]

  if (syscheck == "Windows") {
    fold <- system2(
      command = path_to_RNAfold,
      args = file.path(wkdir, "converted.fasta"),
      stdout = TRUE,
      wait = TRUE,
      invisible = TRUE
    )
  } else if (syscheck == "Linux" | syscheck == "Darwin") {
    fold <- system(
      command = paste(path_to_RNAfold, file.path(wkdir, "converted.fasta"), sep = " "),
      intern = TRUE
    )
  } else {
    warning("Operating system is not Windows or Linux/Darwin. Returning empty handed.")
    return(NULL)
  }

  # RNAFold - Version 2.7.0
  # The return object from RNAfold can have anywhere from 0-4 elements for a single input sequence
  #
  # Length 0:  If the return object is empty, that means no sequence was present in the input file
  #            This is true even if a fasta header is present but no sequence is present.
  #
  # Length 1:  The only time I've encountered this so far is when there was an error running RNAfold
  #            One example of such an error is if the input file passed to RNAfold wasn't present or
  #            Couldn't be read for some reason.
  #            FIELD 1: [ERROR] - Error Message
  #
  # Length 2:  This should only occur when a sequence is present in the input file but no fasta header
  #            is present.
  #            FIELD 1: Input Sequence
  #            FIELD 2: Dot Bracket and MFE
  #
  # Length 3:  This occurs when a single properly formatted fasta header and sequence are present in
  #            the input file. This could also occur when a warning is present but no fasta header
  #            was present.
  #            FIELD 1: Fasta Header or [WARNING] - Warning Message
  #            FIELD 2: Input Sequence
  #            FIELD 3: Dot Bracket and MFE
  #
  # Length 4+: This can occur in 2 known situations. The first is when multiple input sequences are
  #            passed into RNAfold. The second is when a warning is present along with the results.
  #            For the purposes of this test, the fields below assume only a single input sequence was
  #            provided.
  #            FIELD 1: [WARNING] - Warning Message
  #            FIELD 2: Fasta Header
  #            FIELD 3: Input Sequence
  #            FIELD 4: Dot Bracket and MFE
  # Length 4+ and length sequence > 8000? 
  #            FIELD 1: FASTA header 
  #            FIELD 2: Sequence part 1
  #            FIELD 3: Sequence part 2
  #            FIELD 4: Dot bracket part 1
  #            FIELD 5: Dot bracket part 2 and MFE
  # 
  
  
  
  if (length(fold) == 0) {                        # EMPTY RESULTS
    warning("Empty RNAfold results, check input")
    
  } else if (length(fold) == 1 &&
             length(grep("ERROR", fold)) > 0) {   # ERROR
    warning("Error running RNAfold, check input")
    
  } else if (length(fold) == 2) {                 # Missing Fasta Header
    output_file <- "rna.ps"
    vien_seq <- fold[1]
    vien_mfe <- stringi::stri_split_fixed(fold[2], pattern = " ", n = 2)[[1]][2]
    vien_struct <- stringi::stri_split_fixed(fold[2], pattern = " ", n = 2)[[1]][1]
    
  } else if (length(fold) == 3 &&
             length(grep("WARNING", fold)) == 0) { # Normal result
    output_file <- paste0(region_string, "_ss.ps")
    vien_seq <- fold[2]
    vien_mfe <- stringi::stri_split_fixed(fold[3], pattern = " ", n = 2)[[1]][2]
    vien_struct <- stringi::stri_split_fixed(fold[3], pattern = " ", n = 2)[[1]][1]
    
  } else if (length(fold) == 3 &&
             length(grep("WARNING", fold)) > 0) { # Missing Fasta Header with Warning
    output_file <- "rna.ps"
    vien_seq <- fold[2]
    vien_mfe <- stringi::stri_split_fixed(fold[3], pattern = " ", n = 2)[[1]][2]
    vien_struct <- stringi::stri_split_fixed(fold[3], pattern = " ", n = 2)[[1]][1]
    
  } else if (length(fold) == 4 &&
             length(grep("WARNING", fold)) > 0) { # Warning result 
    output_file <- paste0(region_string, "_ss.ps")
    vien_seq <- fold[3]
    vien_mfe <- stringi::stri_split_fixed(fold[4], pattern = " ", n = 2)[[1]][2]
    vien_struct <- stringi::stri_split_fixed(fold[4], pattern = " ", n = 2)[[1]][1]
    
  } else if(length(fold) > 4 &&
            identical(grep("WARNING", fold), integer(0)) ){
    output_file <- paste0(region_string, "_ss.ps")
    vien_seq <- paste0(fold[2], fold[3])
    vien_mfe <- stringi::stri_split_fixed(fold[5], pattern = " ", n = 2)[[1]][2]
    struct <- paste0(fold[4], fold[5])
    vien_struct <- stringi::stri_split_fixed(struct, pattern = " ", n = 2)[[1]][1]
  }
  
  ct <- RRNA::makeCt(vien_struct, vien_seq)

  # Using capture.output to silence the excessive console output from ct2coord
  utils::capture.output(coord <- RRNA::ct2coord(ct), file = nullfile())
  # RRNA::RNAPlot(coord, nt = TRUE)
  # split the string to get mfe
  # Split based on space and remove parentheses
  mfe <- gsub(" ", "", gsub("[)]", "", gsub("[(]", "", vien_mfe)))
  mfe <- as.numeric(mfe)
  start <- start
  stop <- stop

  # RNAfold creates files in the base directory in which MiSiPi was originally called
  # Move to results siRNA directory
  file.rename(from = output_file,
              to = file.path(wkdir, output_file))

  return(list(start, stop, mfe, coord, vien_struct))
}
