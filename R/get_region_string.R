# Function to revert the bed file positions back to zero based half-open for writing results
#
# @param chrom Character: The region name
# @param start Integer: The one based region start position
# @param stop Integer: The region stop position
# @return Character: The region name and zero based positions

.get_region_string <- function(chrom, start, stop) {
  # While normal bed notation would list this as [chromosome:start-end],
  # This string is used to name some files, and colons are not allowed in file names
  # Underscore is being used in place of colon
  # Only the start position needs to be decremented as bed files use half-open coordinates
  return(paste0(chrom, "-", start - 1, "_", stop))
}
