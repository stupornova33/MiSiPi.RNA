#' rna_plot
#' calls rnafold and rnaplot from viennaRNA and colors the most abundant reads on the structure.
#' @param path_to_RNAfold a string, the full path to the RNAfold executable.
#' @param path_to_RNAplot a string, the full path to the RNAplot executable.
#' @param wkdir a string, the path to the desired working directory.
#' @param pos_df a data.frame containing the starts and ends of the abundant reads.
#' @param colors a list of two colors for highlighting abundant reads on the structure.
#' @return nothing
#' @export

rna_plot <- function(path_to_RNAfold, path_to_RNAplot, wkdir, pos_df, colors){

  #need to add the arg in fold to "fold_short_rna"
  fold <- system2(command = path_to_RNAfold, args = paste0(wkdir, "converted.fasta --outfile=", "converted.txt"), stdout = TRUE, wait = TRUE, invisible = TRUE)

  #example arguments
  #a <- '--pre="1 15 8 RED omark" test.txt'

  pos1 <- pos_df$r1_start
  pos2 <- pos_df$r1_end
  pos3 <- pos_df$r2_start
  pos4 <- pos_df$r2_end

  final_arg <- paste0('--pre="', pos1, ' ', pos2, ' ', 8, ' ', colors[1], ' ', 'omark ', pos3, ' ', pos4, ' ', 9, ' ', colors[2], ' ', 'omark" ', 'converted.txt')

  system2(command = path_to_RNAplot,
          args = final_arg,
          stdout = TRUE, wait = TRUE, invisible = TRUE)

}













