
fold_short_rna <- function(start, stop, converted, path_to_RNAfold){
   cat(converted, file = "converted.fasta", sep = "\n", append = FALSE)
   #fold <- system2(command = path_to_RNAfold, args = "converted.fasta", stdout= TRUE, wait = TRUE, invisible = TRUE)

    syscheck <- unlist(unname(Sys.info()[1]))

   if(syscheck == "Windows") {
     fold <- system2(command = path_to_RNAfold, args = "converted.fasta", stdout= TRUE, wait = TRUE, invisible = TRUE)
   } else if(syscheck == "Linux"){
     fold <- system(command = path_to_RNAfold, input = converted, stdout= TRUE, wait = TRUE, invisible = TRUE)
   } else {
     print("Operating system is not Windows or Linux. Halt.")
     return()
   }

   pos_vec <- seq(2, length(fold), by = 2)
   sep_seqs <- function(x){
      vien_split <- stringi::stri_split_fixed(fold[2], pattern = " ")[[1]][2]
      vien_struct <- stringi::stri_split_fixed(fold[2], pattern = " ")[[1]][1]
      return(list(vien_split, vien_struct))
   }

   test <- lapply(pos_vec, sep_seqs)

   make_fold <- function(x){
      vien_struct <- test[[x]][[2]]
      print(vien_struct)
      vien_split <- test[[x]][[1]]
      print(vien_split)
      ct <- RRNA::makeCt(vien_struct, fold[[1]])
      coord <- RRNA::ct2coord(ct)

      #split the string to get mfe
      mfe <- gsub(' ', '', gsub('[)]','', gsub('[(]', '', vien_split)))
      print(mfe)
      mfe <- as.numeric(mfe)
      start <- start
      stop <- stop
      converted <- converted
      return(list(mfe, start, stop, coord, vien_struct, converted))
   }

   res <- lapply(1:length(test), make_fold)

   return(res)
}
