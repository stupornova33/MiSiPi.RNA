# function to create a GR obj from a data table of reads
# needed for finding overlaps
# returns a GR obj
#
# @param data a data table
# @param name a string passed in by user or bed file
# @param reg_start a whole number
# @param reg_stop a whole number
# @return GRobj

.makeGRobj <- function(data, name, reg_start, reg_stop) {
  GRobj <- GenomicRanges::makeGRangesFromDataFrame(data,
    keep.extra.columns = TRUE,
    ignore.strand = TRUE,
    seqnames.field = name,
    start.field = reg_start,
    end.field = reg_stop
  )
  return(GRobj)
}
