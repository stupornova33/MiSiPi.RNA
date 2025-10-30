
.get_negative_bed_coords <- function(start, stop) {
  return(start < 0 | stop < 0)
}

.get_inverted_bed_coords <- function(start, stop) {
  return(start > stop)
}

.get_missing_bed_chr <- function(bed_chr, bam_chr) {
  return(!(bed_chr %in% bam_chr))
}
