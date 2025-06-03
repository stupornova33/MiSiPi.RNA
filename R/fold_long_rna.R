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
    warning("Operating system is not Windows or Linux/Darwin. Returning empty handed.")
    return(NULL)
  }

  if (length(grep("WARNING", fold)) > 0) { # warning can be ignored
    fold <- fold[1:(grep("WARNING", fold)) - 1]
  }

  # RNAfold returns 3 objects in some cases and only 2 in others
  # Need to extract the vienna things differently based on this
  if (length(fold) > 2) {
    vien_split <- stringi::stri_split_fixed(fold[4], pattern = " ")[[1]][2]
    vien_struct <- paste0(fold[3], stringi::stri_split_fixed(fold[4], pattern = " ")[[1]][1])
  } else {
    # In some cases, RNAFold returns a string with a different amount of spaces in the result
    # This can throw of string splitting
    # Adding the parameter n = 2 ensures that only 1 space is considered when
    # Splitting the result string into its separate parts
    vien_split <- stringi::stri_split_fixed(fold[[2]], pattern = " ", n = 2)[[1]][2]
    vien_struct <- stringi::stri_split_fixed(fold[[2]], pattern = " ", n = 2)[[1]][1]
  }

  ct <- RRNA::makeCt(vien_struct, fold[1])

  # Using capture.output to silence the excessive console output from ct2coord
  utils::capture.output(coord <- RRNA::ct2coord(ct), file = nullfile())
  # RRNA::RNAPlot(coord, nt = TRUE)
  # split the string to get mfe
  # Split based on space and remove parentheses
  mfe <- gsub(" ", "", gsub("[)]", "", gsub("[(]", "", vien_split)))
  mfe <- as.numeric(mfe)
  start <- start
  stop <- stop

  ps_file_to_remove <- file.path(dirname(dirname(wkdir)), "rna.ps")
  file.remove(ps_file_to_remove)

  return(list(start, stop, mfe, coord, vien_struct))
}
