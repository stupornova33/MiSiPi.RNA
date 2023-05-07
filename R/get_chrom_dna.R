#' function to subset one chromosome sequence from a genome file
#' takes genome object
#' returns vector with dna sequence
#' 
#' @param genome a genome fasta file read by Biostrings
#' @return dna

#' @export

get_chrom_dna <- function(genome){
   name_list <- names(genome)
   name_split <- vector()
   for(i in 1:length(name_list)){
      name_split[i] <- stringr::str_split(name_list[[i]], " ")[[1]][1]
   }
   names(genome) <- name_split
   
   idx <- which(chr_name == chrom_name)
   mygranges <- GenomicRanges::GRanges(
      seqnames = c(chrom_name),
      ranges = IRanges::IRanges(start=c(1), end=c(chr_length[idx])),
      strand = c("-"))
   
   chr_dna <- BSgenome::getSeq(genome, mygranges)
   size <- nrow(read_pileups)
   
   dna <- as.character(chr_dna)
   return(dna)
}