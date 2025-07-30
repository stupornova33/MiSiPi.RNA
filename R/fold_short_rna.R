# function to run RNAfold
# processes output of RNA fold to get MFE and vien struct
# returns list of values for each region
# @param start an integer
# @param stop an integer
# @param converted a vector containing a sequence
# @param path_to_RNAfold a string
# @param chrom_name a string
# @param wkdir a string
# @return list

.fold_short_rna <- function(start, stop, converted, path_to_RNAfold, chrom_name, wkdir) {
  write.table(converted,
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
    warning("Operating system is not Windows or Linux/Darwin. Returning empty handed")
    return(NULL)
  }

  # Delete unwanted .ps file
  ps_filename_to_remove <- paste0(chrom_name, "-", start - 1, "_", stop - 1, "_ss.ps")
  file.remove(ps_filename_to_remove)

  
  if (length(fold) == 0) {                        # EMPTY RESULTS
    warning("Empty RNAfold results, check input")
    
  } else if (length(fold) == 1 &&
             length(grep("ERROR", fold)) > 0) {   # ERROR
    warning("Error running RNAfold, check input")
    
  } else if (length(fold) == 2) {                 # Missing Fasta Header
    vien_seq <- fold[1]
    vien_mfe <- stringi::stri_split_fixed(fold[2], pattern = " ", n = 2)[[1]][2]
    vien_struct <- stringi::stri_split_fixed(fold[2], pattern = " ", n = 2)[[1]][1]
    
  } else if (length(fold) == 3 &&
             length(grep("WARNING", fold)) == 0) { # Normal result
    vien_seq <- fold[2]
    vien_mfe <- stringi::stri_split_fixed(fold[3], pattern = " ", n = 2)[[1]][2]
    vien_struct <- stringi::stri_split_fixed(fold[3], pattern = " ", n = 2)[[1]][1]
    
  } else if (length(fold) == 3 &&
             length(grep("WARNING", fold)) > 0) { # Missing Fasta Header with Warning
    vien_seq <- fold[2]
    vien_mfe <- stringi::stri_split_fixed(fold[3], pattern = " ", n = 2)[[1]][2]
    vien_struct <- stringi::stri_split_fixed(fold[3], pattern = " ", n = 2)[[1]][1]
    
  } else if (length(fold) == 4 &&
             length(grep("WARNING", fold)) > 0) { # Warning result 
    vien_seq <- fold[3]
    vien_mfe <- stringi::stri_split_fixed(fold[4], pattern = " ", n = 2)[[1]][2]
    vien_struct <- stringi::stri_split_fixed(fold[4], pattern = " ", n = 2)[[1]][1]
    
  }
  
  ct <- RRNA::makeCt(vien_struct, vien_seq)
  
  # Using capture.output to silence the excessive console output from ct2coord
  # check to see if no bases are paired first
  if (sum(ct$bound) > 0) {
    utils::capture.output(coord <- RRNA::ct2coord(ct), file = nullfile())
    mfe <- gsub(" ", "", gsub("[)]", "", gsub("[(]", "", vien_mfe)))
    mfe <- as.numeric(mfe)
  } else {
    coord <- NULL
    mfe <- 0
  }

  start <- start
  stop <- stop
  converted <- converted
  
  return(list(
    "start" = start,
    "stop" = stop,
    "mfe" = mfe,
    "vienna" = vien_struct,
    "converted" = converted,
    "extracted_df" = coord
  ))
}
