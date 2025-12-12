#' eps2png
#' @param path_to_magick_exe A string specifying the full path to the ImageMagick .exe file
#' @param file_dir A string specifying the path to the .eps files to be converted. This will be the location of the converted .png files
#' @return nothing
#' @export

eps2png <- function(path_to_magick_exe, file_dir) {
  old_dir <- setwd(file_dir)
  eps_files <- list.files(file_dir, pattern = "_ss.eps")

  for (i in 1:length(eps_files)) {
    input_file <- file.path(file_dir, eps_files[i])
    out_file <- paste0(tools::file_path_sans_ext(input_file), ".png")
    args <- paste0("-density 300 ", input_file, " -quality 100 ", out_file)
    #try(
      system2(path_to_magick_exe, args = args, stderr = NULL)#,
    #  silent = TRUE
    #)
    
  }
  
  ps_files <- list.files(file_dir, pattern = "_ss.ps")
  for (i in 1:length(ps_files)) {
    input_file <- file.path(file_dir, ps_files[i])
    out_file <- paste0(tools::file_path_sans_ext(input_file), ".png")
    args <- paste0("-density 300 ", input_file, " -quality 100 ", out_file)
    #try(
    system2(path_to_magick_exe, args = args, stderr = NULL)#,
    #  silent = TRUE
    #)
    
  }
  
  setwd(old_dir)
}
