# Align two strings
# right_stem_seq - Character
# left_stem_seq - Character
# num_loops - Integer

.align_stems <- function(right_stem_seq, left_stem_seq, num_loops) {
  # convert parentheses to match right arm

  a <- c(5, -4, -4, -4)
  t <- c(-4, 5, -4, -4)
  g <- c(-4, -4, 5, -4)
  c <- c(-4, -4, -4, 5)

  dnafull <- matrix(c(a, t, g, c), ncol = 4, nrow = 4)
  rownames(dnafull) <- c("A", "U", "G", "C")
  colnames(dnafull) <- c("A", "U", "G", "C")

  l <- collapse(left_stem_seq)
  r <- collapse(right_stem_seq)

  l2 <- Biostrings::BString(l)
  r2 <- Biostrings::RNAString(r)
  r2 <- Biostrings::complement(r2)
  r2 <- Biostrings::BString(r2)

  # this may be best parameters for pairwise alignment? Came from default Needle
  pa1 <- Biostrings::pairwiseAlignment(l2, r2,
    substitutionMatrix = dnafull, type = "local-global",
    gapOpening = 10, gapExtension = 0.5
  )

  S1 <- as.character(Biostrings::alignedSubject(pa1))[[1]]
  P1 <- as.character(Biostrings::alignedPattern(pa1))[[1]]

  left_arm <- P1
  right_arm <- S1

  return(list(left_arm = left_arm, right_arm = right_arm))
}
