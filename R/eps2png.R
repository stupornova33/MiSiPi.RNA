#' eps2png
#' @param path_to_magick_exe A string specifying the full path to the ImageMagick .exe file
#' @param file_dir A string specifying the path to the .eps files to be converted. This will be the location of the converted .png files
#' @return nothing
#' @export

eps2png <- function(path_to_magick_exe, file_dir) {
  old_dir <- getwd()
  setwd(file_dir)
  eps_files <- list.files(file_dir, pattern = "_ss.eps")

  for (i in 1:length(eps_files)) {
    input_file <- file.path(file_dir, eps_files[i])
    print(input_file)
    file.exists(input_file)
    output_file_pref <- unlist(strsplit(input_file, ".eps"))
    out_file <- paste0(file_dir, "/", output_file_pref, ".png")
    args <- paste0("-density 300 ", input_file, " -quality 100 ", out_file)
    system2(path_to_magick_exe, args = args)
  }
  setwd(old_dir)
}