#' ps2png
#' @param path_to_magick_exe A string specifying the full path to the ImageMagick .exe file
#' @param file_dir A string specifying the path to the .ps files to be converted. This will be the location of the converted .png files
#' @return nothing
#' @export



ps2png <- function(path_to_magick_exe, file_dir){
  setwd(file_dir)
  ps_files <- list.files(file_dir, pattern = "_ss.ps")
  last <- tail(unlist(strsplit(file_dir, "")),1)
  

  for(i in 1:length(ps_files)){
    if(last == "\\/"){
      input_file <- paste0(file_dir, ps_files[i])
      output_file_pref <- unlist(strsplit(ps_files[i], ".ps"))
      out_file <- paste0(file_dir, output_file_pref, ".png")
    } else {
      input_file <- paste0(file_dir, "\\/", ps_files[i])
      output_file_pref <- unlist(strsplit(ps_files[i], ".ps"))
      out_file <- paste0(file_dir,"\\/", output_file_pref, ".png")
      
    }

    args <- paste0( "-density 300 ",input_file, " -quality 100 ", out_file)
    system2(path_to_magick_exe, args = args)
  }

}
