#' miRNA function
#' @param converted a string
#' @param path_to_RNAfold a string

#' @return plots
#' @export


plot_miRNA_struct <- function(converted, path_to_RNAfold){

fold_list <- t(mapply(fold_rna, "chrom", 1, 107, converted, path_to_RNAfold))

vienna <- unlist(unname(fold_list[,5]))
len <- length(strsplit(vienna, '')[[1]])

#returns left loop coords before right loop coords
# if there is third loop?
orig_loop_pos <- unlist(stringi::stri_locate_all_regex(vienna, '\\(\\.*\\)'))
vienna_vec <- unlist(strsplit(vienna, ''))

if(length(orig_loop_pos) == 2){
   left_arm <- vienna_vec[1:orig_loop_pos[1]]
   right_arm <- vienna_vec[orig_loop_pos[2]:length(vienna_vec)]
   loop_pos <- orig_loop_pos
   #reverse right arm so pairs sort of line up
   right_arm <- rev(right_arm)
   
   #collapse for regex matching
   
   collapse <- function(string){
      vienna <- paste(string, collapse = "")
      return(vienna)
   }
   left_vienna <- collapse(left_arm)
   right_vienna <- collapse(right_arm)
   
   #convert parentheses to match right arm
   l <- gsub("[(]", ")", left_vienna)
   l <- gsub("[)]", "A", l)
   l <- gsub("[.]", "S", l)
   
   r <- gsub("[(]", ")", right_vienna)
   r <- gsub("[)]", "(", r)
   r <- gsub("[(]", "A", r)
   r <- gsub("[.]", "S", r)

} else { #if num loops is 2 or more
   loop_1_pos <- c(orig_loop_pos[1], orig_loop_pos[3])
   loop_2_pos <- c(orig_loop_pos[2], orig_loop_pos[4])
   half_len <- len/2
   dist1 <- abs(loop_1_pos[1] - half_len)
   dist2 <- abs(loop_2_pos[2] - half_len)
   if(dist1 < dist2){
      loop_pos <- loop_1_pos
   } else {
      loop_pos <- loop_2_pos
   }
   
   left_arm <- vienna_vec[1:loop_pos[1]]
   right_arm <- vienna_vec[loop_pos[2]:length(vienna_vec)]
   
   #reverse right arm so pairs sort of line up
   right_arm <- rev(right_arm)
   
   #collapse for regex matching
   
   collapse <- function(string){
      vienna <- paste(string, collapse = "")
      return(vienna)
   }
   
   left_vienna <- collapse(left_arm)
   right_vienna <- collapse(right_arm)
   
   
   minor_loop_end <- unlist(stringi::stri_locate_all_regex(vienna, '\\)\\.*\\('))
   if(minor_loop_end[1] > half_len) { #minor loop is on right arm
      minor_loop_pos <- unlist(stringi::stri_locate_all_regex(right_vienna, '\\(\\.*\\)'))
   } else {
      minor_loop_pos <- unlist(stringi::stri_locate_all_regex(left_vienna, '\\(\\.*\\)'))
   }
   
   #only want to align major stem
   r_major_stem <- right_arm[minor_loop_pos[2]:length(right_arm)]
   r_other_vienna <- right_arm[1:minor_loop_pos[2]]
   
   #may need to change this
   l_major_stem <- left_arm[(length(left_arm) - (length(r_major_stem))): length(left_arm)]
   l_other_vienna <- left_arm[1:(length(left_arm) - (length(r_major_stem)))]
   
   new_right_vienna <- collapse(r_major_stem)
   new_left_vienna <- collapse(l_major_stem)
   #convert parentheses to match right arm
   l <- gsub("[(]", ")", new_left_vienna)
   l <- gsub("[)]", "A", l)
   l <- gsub("[.]", "S", l)
   
   r <- gsub("[(]", ")", new_right_vienna)
   r <- gsub("[)]", "(", r)
   r <- gsub("[(]", "A", r)
   r <- gsub("[.]", "S", r)
   
}

## end processing by number of loops
l2 <- Biostrings::BString(l)
r2 <- Biostrings::BString(r)

pa1 <- Biostrings::pairwiseAlignment(l2,r2, substitutionMatrix=NULL, type= "local-global",
                        gapOpening=1, gapExtension=2)
pa1


Biostrings::writePairwiseAlignments(pa1, file = 'aligned_vienna.txt')

test <- list(readLines("aligned_vienna.txt"))

P1 <- test[[1]][23] 
S1 <- test[[1]][25]

P1 <- gsub('P1', '', P1)
S1 <- gsub('S1', '', S1)

P1 <- gsub('[[:digit:]]+', '', P1)
S1 <- gsub('[[:digit:]]+', '', S1)

left_arm <- gsub(' ', '', P1)
right_arm <- gsub(' ', '', S1)

#get sequence 
split_seq <- unlist(strsplit(seq, ''))

if(length(orig_loop_pos) == 2){
   left_seq <- split_seq[1:loop_pos[1]]
   right_seq <- split_seq[loop_pos[2]:length(vienna_vec)]
   right_seq <- rev(right_seq)
   
   left_collapsed <- paste(left_arm, collapse = "")
   right_collapsed <- paste(right_arm, collapse = "")
   
   left_dashes <- unlist(stringi::stri_locate_all_regex(left_collapsed, "-")) %>% unique()
   right_dashes <-  unlist(stringi::stri_locate_all_regex(right_collapsed, "-")) %>% unique()
   
   
} else {
   left_seq <- split_seq[1:loop_pos[1]]
   right_seq <- split_seq[loop_pos[2]:length(split_seq)]
   right_seq <- rev(right_seq)
   
   r_stem_seq <- right_seq[minor_loop_pos[2]:length(right_seq)]
   r_other_seq <- right_seq[1:minor_loop_pos[2]]
   
   #get major stem seqs
   #may need to change this
   l_stem_seq <- left_seq[(length(left_seq) - (length(r_major_stem))): length(left_seq)]
   l_other_seq <- left_seq[1:(length(left_seq) - (length(r_major_stem)))]   
   
   # ?
   #left_collapsed <- paste(left_arm, collapse = "")
   #right_collapsed <- paste(right_arm, collapse = "")
   
   left_dashes <- unlist(stringi::stri_locate_all_regex(left_arm, "-")) %>% unique()
   right_dashes <-  unlist(stringi::stri_locate_all_regex(right_arm, "-")) %>% unique()
   
}





add_dashes <- function(x){
   print(left_dashes[x])
   append(left_seq, "-", after = left_dashes[x] - 1)
}

test <- lapply(1:length(left_dashes), add_dashes)

#insert dashes into sequence in correct place
for(i in 1:length(left_dashes)){
   print(left_dashes[i] - 1 )
   left_seq <- append(left_seq, "-", after = left_dashes[i] - 1)
   
}

for(i in 1:length(right_dashes)){
   print(right_dashes[i] - 1 )
   right_seq <- append(right_seq, "-", after = right_dashes[i] - 1)
   
}


#convert A's and S's back to parentheses
l <- gsub("A", "\\)", left_arm)
l <- unlist(strsplit(gsub("S", ".", l), ''))

r <- gsub( "A", "\\)", right_arm)
r <- unlist(strsplit(gsub( "S", ".", r), ''))


#make df for loop
loop_seq <- split_seq[(loop_pos[1] + 1):(loop_pos[length(loop_pos)] - 1)]
loop_vienna <- vienna_vec[(loop_pos[1] + 1):(loop_pos[length(loop_pos)] - 1)]

left_loop_seq <- vector()
right_loop_seq <- vector()
left_loop_vienna <- vector()
right_loop_vienna <- vector()
if(length(loop_seq) %% 2 == 1){ #if loop length is odd
   half_len <- round(length(loop_vienna)/2)
   left_loop_seq <- append(left_loop_seq, loop_seq[1:half_len])
   right_loop_seq <- append(right_loop_seq, loop_seq[half_len:length(loop_seq)])
   left_loop_vienna <- append(left_loop_vienna, loop_vienna[1:half_len])
   right_loop_vienna <- append(right_loop_vienna, loop_vienna[half_len:length(loop_seq)])

} else { #if length loop is even
   half_len <- length(loop_vienna)/2
   left_loop_seq <- append(left_loop_seq, loop_seq[1:half_len])
   right_loop_seq <- append(right_loop_seq, loop_seq[half_len:length(loop_seq)])
   left_loop_vienna <- append(left_loop_vienna, loop_vienna[1:half_len])
   right_loop_vienna <- append(right_loop_vienna, loop_vienna[half_len:length(loop_seq)])

   
}

#create middle row of plot
char_vec <- vector()
for(i in 1:length(l)){
   if((l[i] == ")") || (l[i] == "(")){
      char_vec <- append(char_vec, "|")
   } else {
      char_vec <- append(char_vec, " ")
   }
}


#make the top row of plot (unpaired bases on left arm)
top_vec <- vector()
mid_top <- vector()
for(i in 1:length(l)){
   if(l[i] == "."|| l[i] == "-"){
      top_vec <- append(top_vec, left_seq[i])
      mid_top <- append(mid_top, " ")
   } else {
      top_vec <- append(top_vec, " ")
      mid_top <- append(mid_top, left_seq[i])
   }
}


#make bottom row of plot (unpaired bases on right arm)
bottom_vec <- vector()
mid_bottom <- vector()
for(i in 1:length(r)){
   if(r[i] == "."|| r[i] == "-"){
      bottom_vec <- append(bottom_vec, right_seq[i])
      mid_bottom <- append(mid_bottom, " ")
   } else {
      bottom_vec <- append(bottom_vec, " ")
      mid_bottom <- append(mid_bottom, right_seq[i])
   }
}

div <- (length(loop_seq) - 2)/2
if(length(loop_seq) %% 2 == 0){ #if length is even
   top_loop <- loop_seq[1:div]
   mid_loop <- loop_seq[(div + 1):(div + 2)]
   bottom_loop <- loop_seq[(div + 3):length(loop_seq)]
} else { #if odd
   top_loop <- loop_seq[1:div]
   mid_loop <- loop_seq[(div + 1):(div + 1)]
   bottom_loop <- loop_seq[(div + 2):length(loop_seq)]
}

mid_top_loop <- mid_loop[1]
mid_bottom_loop <- mid_loop[2]

if(length(orig_loop_pos) >2){
   if(length(l_other_seq) > length(r_other_seq)){
      l_bool <- "TRUE"
      longest_arm <- l_other_seq
      shortest_arm <- r_other_seq
   } else {
      l_bool <- "FALSE"
      longest_arm <- r_other_seq
      shortest_arm <- l_other_seq
   }
   
   other_char_vec <- rep(" ", times = length(longest_arm))
   seq_diff <- length(longest_arm) - length(shortest_arm)
   add_dashes <- rep("-", times = seq_diff)
   
   #concatenate dfs
   if(l_bool == 'FALSE'){
      new_top_vec <- c(add_dashes, l_other_seq, top_vec)
      mid_top_vec <- c(other_char_vec, mid_top)
      new_char_vec <- c(other_char_vec, char_vec)
      mid_bottom_vec <- c(other_char_vec, mid_bottom)
      new_bottom_vec <- c(r_other_seq, bottom_vec)
      bottom_vec <- new_bottom_vec
      mid_bottom <- mid_bottom_vec
      top_vec <- new_top_vec
      mid_top <- mid_top_vec
      char_vec <- new_char_vec
      
   } else {
      new_top_vec <- c(add_dashes, r_other_seq, top_vec)
      mid_top_vec <- c(other_char_vec, mid_top)
      new_char_vec <- c(other_char_vec, char_vec)
      mid_bottom_vec <- c(other_char_vec, mid_bottom)
      new_bottom_vec <- c(l_other_seq, bottom_vec)
      bottom_vec <- new_bottom_vec
      mid_bottom <- mid_bottom_vec
      top_vec <- new_top_vec
      mid_top <- mid_top_vec
      char_vec <- new_char_vec
      
   }
   
}



top_df <- data.frame("x" = c(1:length(top_vec)), "y" = c(5), "text" = top_vec)
top_loop_df <- data.frame("x" = seq((length(top_vec) + 1), (length(top_vec) + length(top_loop))), "y" = c(5), "text" = top_loop)
mid_top_loop_df <- data.frame("x" = seq((length(mid_top) + length(top_loop) + 1), (length(mid_top) + length(top_loop) + length(mid_top_loop))), "y" = c(4.9), "text" = mid_top_loop)
mid_top_df <- data.frame("x" = c(1:length(mid_top)), "y" = c(4.9), "text" = mid_top)

char_df <- data.frame("x" = c(1:length(char_vec)), "y" = c(4.8), "text" = char_vec)

mid_bottom_df <- data.frame("x" = c(1:length(mid_bottom)), "y" = c(4.7), "text" = mid_bottom)
mid_bottom_loop_df <- data.frame("x" = seq(((length(mid_bottom) + length(top_loop)) + 1), (length(mid_top) + length(top_loop) + length(mid_bottom_loop))), "y" = c(4.7), "text" = mid_bottom_loop)
bottom_df <- data.frame("x" = c(1:length(bottom_vec)), "y" = c(4.6), "text" = bottom_vec)
bottom_loop_df <- data.frame("x" = seq((length(bottom_vec) + 1), (length(bottom_vec) + length(bottom_loop))), "y" = c(4.6), "text" = bottom_loop)

f_df <- rbind(top_df, top_loop_df, mid_top_loop_df, mid_top_df, char_df, mid_bottom_df, mid_bottom_loop_df, bottom_df, bottom_loop_df)


plot_miRNA_text <- function(f_df){
   g <- ggplot2::ggplot(f_df, ggplot2::aes(x,y, label = text)) +
      ggplot2::geom_text() + 
      ggplot2::xlim(0, max(length(top_vec), (length(top_vec) + length(top_loop) + length(mid_top_loop)))) + 
      
      ggplot2::scale_y_continuous(breaks = seq(2,5, by = 1), limits = c(2,5))+ 
      ggplot2::theme(plot.margin = ggplot2::unit(c(0, 0, 0, 0), "cm")) +
      ggplot2::theme_void()
   
   return(g)
}


plot <- plot_miRNA_text(f_df)

}