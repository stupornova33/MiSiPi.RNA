.write.quiet <- function(x, file = "") {
  append <- FALSE
  col.names <- TRUE
  if (file.exists(file) & file.size(file) > 0) {
    append <- TRUE
    col.names <- FALSE
  }
  suppressWarnings(
    write.table(x = x, file = file, sep = "\t", quote = FALSE, append = append, col.names = col.names, row.names = FALSE, na = "NA")
  )
}
