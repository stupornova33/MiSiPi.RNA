#' function to run RNAfold
#' processes output of RNA fold to get MFE and vien struct
#' returns list of values for each region
#' @param start an integer
#' @param stop an integer
#' @param converted a vector containing a sequence
#' @param path_to_RNAfold a string
#' @param chrom_name a string
#' @param wkdir a string
#' @return list
#' @export

fold_short_rna <- function(start, stop, converted, path_to_RNAfold, chrom_name, wkdir) {
  write.table(converted, file = file.path(wkdir, "converted.fasta"), sep = "\n", append = FALSE, row.names = FALSE, quote = FALSE)
  #fold <- system2(command = path_to_RNAfold, args = "converted.fasta", stdout= TRUE, wait = TRUE, invisible = TRUE)

  syscheck <- unlist(unname(Sys.info()[1]))

  if (syscheck == "Windows") {
    fold <- system2(command = path_to_RNAfold,
                    args = file.path(wkdir, "converted.fasta"),
                    stdout= TRUE,
                    wait = TRUE,
                    invisible = TRUE)
  } else if (syscheck == "Linux" | syscheck == "Darwin") {
    fold <- system(command = path_to_RNAfold, input = converted, intern = TRUE)
  } else {
    print("Operating system is not Windows or Linux. Halt.")
    return()
  }

  # Delete unwanted .ps file
  ps_filename_to_remove <- paste0(chrom_name, "-", start - 1, "_", stop - 1, "_ss.ps")
  file.remove(file.path(ps_filename_to_remove))
  
  pos_vec <- seq(2, length(fold), by = 2)
  sep_seqs <- function(x) {
    vien_split <- stringi::stri_split_fixed(fold[3], pattern = " ")[[1]][2]
    vien_struct <- stringi::stri_split_fixed(fold[3], pattern = " ")[[1]][1]
    return(list(vien_split, vien_struct))
  }

  test <- lapply(pos_vec, sep_seqs)

  make_fold <- function(x) {
    vien_struct <- test[[x]][[2]]
    print(vien_struct)
    vien_split <- test[[x]][[1]]
    print(vien_split)
    ct <- RRNA::makeCt(vien_struct, fold[[1]])
    # Using capture.output to silence the excessive console output from ct2coord
    utils::capture.output(coord <- RRNA::ct2coord(ct), file = nullfile())

    #split the string to get mfe
    mfe <- gsub(' ', '', gsub('[)]','', gsub('[(]', '', vien_split)))
    print(mfe)
    mfe <- as.numeric(mfe)
    start <- start
    stop <- stop
    converted <- converted
    return(list("start" = start,
                "stop" = stop,
                "mfe" = mfe,
                "vienna" = vien_struct,
                "converted" = converted,
                "extracted_df" = coord))
  }

  res <- lapply(1:length(test), make_fold)

  return(res)
}
