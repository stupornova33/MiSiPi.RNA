#' function to run RNAfold
#' processes output of RNA fold to get MFE and vien struct
#' returns list of values for each region
#' @param chrom_name a string
#' @param start a whole number
#' @param stop a whole number
#' @param converted a vector containing a sequence
#' @param path_to_RNAfold a string
#' @return list
#' @export


fold_long_rna <- function(chrom_name, start, stop, converted, path_to_RNAfold){
   prefix <- paste0(chrom_name,"-", start, "-", stop)
   cat(converted, file = "converted.fasta", sep = "\n", append = FALSE)

   syscheck <- unlist(unname(Sys.info()[1]))

   if(syscheck == "Windows") {
     fold <- system2(command = path_to_RNAfold, args = "converted.fasta", stdout= TRUE, wait = TRUE, invisible = TRUE)
   } else if(syscheck == "Linux"){
     fold <- system(command = path_to_RNAfold, input = converted)
   } else {
     print("Operating system is not Windows or Linux. Halt.")
     return()
   }
    if(length(grep("WARNING", fold)) > 0){   #warning can be ignored
         fold <- fold[1:(grep("WARNING", fold)) - 1]
      }

   if(length(fold) > 2){
      vien_split <- stringi::stri_split_fixed(fold[4], pattern = " ")[[1]][2]
      vien_struct <- paste0(fold[3],stringi::stri_split_fixed(fold[4], pattern = " ")[[1]][1])

   } else {
      vien_split <- stringi::stri_split_fixed(fold[[2]], pattern = " ")[[1]][2]
      vien_struct <- stringi::stri_split_fixed(fold[[2]], pattern = " ")[[1]][1]

   }

   ct <- RRNA::makeCt(vien_struct, fold[1])
   coord <- RRNA::ct2coord(ct)
   #RRNA::RNAPlot(coord, nt = TRUE)
   #split the string to get mfe
   mfe <- gsub(' ', '', gsub('[)]','', gsub('[(]', '', vien_split)))
   mfe <- as.numeric(mfe)
   start <- start
   stop <- stop


   return(list(start, stop, mfe, coord, vien_struct))
}
