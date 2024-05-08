#' miRNA_function
#' processes reads according to phased miRNA algorithm
#' plots output
#' returns plots
#' @param chrom_name a string
#' @param reg_start a whole number
#' @param reg_stop a whole number
#' @param chromosome chrom names extracted from bam file
#' @param length length of chromosome of interest
#' @param strand a character passed in, "+" or "-"
#' @param min_read_count a whole number
#' @param genome_file a fasta file of chrom sequences
#' @param bam_file a BAM file
#' @param logfile a string
#' @param wkdir a string
#' @param plot_output a string, default = TRUE
#' @param path_to_RNAfold a string
#' @param write_fastas a string, "T" or "F"
#' @param weight_reads Determines whether read counts will be weighted and with which method. Valid options are "weight_by_prop", "locus_norm", a user-defined value, or "none". See MiSiPi documentation for descriptions of the weighting methods.
#' @param out_type The type of file to write the plots to. Options are "png" or "pdf". Default is PDF.
#' @return plots
#' @export

run_miRNA_function <- function(chrom_name, reg_start, reg_stop, chromosome, length,
                           strand, min_read_count, genome_file, bam_file, logfile, wkdir, plot_output, path_to_RNAfold, write_fastas, weight_reads, out_type){
   print(paste0(chrom_name, "-",reg_start,"-", reg_stop))
   cat(file = paste0(wkdir,logfile), paste0("chrom_name: ", chrom_name, " reg_start: ", reg_start, " reg_stop: ", reg_stop, "\n"), append = TRUE)
   pos <- count <- count.x <- count.y <- end <- r1_end <- r1_start <- dist <- r2_end <- r2_start <- lstop <- lstart <- r1_seq <- loop_seq <- r2_seq <- start <- whole_seq <- width <- NULL

   # do not run locus if length is > 300 - not a miRNA. Also avoids issue where user provides coordinates of miRNA cluster.
   if(reg_stop - reg_start > 300){
      cat(file = paste0(wkdir, logfile), "length of region is greater than 300. \n", append = TRUE)
      return(null_mi_res())
   }

   bam_obj <- OpenBamFile(bam_file, logfile)
   bam_header <- Rsamtools::scanBamHeader(bam_obj)
   chr_name <- names(bam_header[['targets']])
   chr_length <- unname(bam_header[['targets']])

   bam_header <- NULL

   which <- GenomicRanges::GRanges(seqnames=chrom_name, IRanges::IRanges(reg_start, reg_stop))

   mygranges <- GenomicRanges::GRanges(
      seqnames = c(chrom_name),
      ranges = IRanges::IRanges(start=c(1), end=c(length)))
   print(paste0("mygranges: ", mygranges))

   geno_seq <- Rsamtools::scanFa(genome_file, mygranges)
   geno_seq <- as.character(unlist(Biostrings::subseq(geno_seq, start = 1, end = length)))

   #for the read size distribution plot
   chrom_m <- getChrMinus(bam_obj, chrom_name, reg_start, reg_stop)
   chrom_p <- getChrPlus(bam_obj, chrom_name, reg_start, reg_stop)

   read_dist <- get_read_dist(bam_obj, chrom_name, reg_start, reg_stop)

   # Moved this code block up so that the getChr functions don't have to be called again
   if(strand == "-"){
      chrom <- chrom_m
   } else {
      chrom <- chrom_p
   }
   chrom_m <- NULL
   chrom_p <- NULL

   if(strand == "-"){
      bam_scan <- Rsamtools::ScanBamParam(flag = Rsamtools::scanBamFlag(isMinusStrand = TRUE), what=c('rname', 'pos', 'qwidth'), which=which)
   } else if(strand == "+"){
      bam_scan <- Rsamtools::ScanBamParam(flag = Rsamtools::scanBamFlag(isMinusStrand = FALSE), what=c('rname', 'pos', 'qwidth'), which=which)
   } else{
      bam_scan <- Rsamtools::ScanBamParam(what=c('rname', 'pos', 'strand','qwidth'), which=which)
   }


   ########################################################## main logic ################################################################
   ## make the read data tables

   # Processing one strand, so make copy of df to transform
   # 1/9/24 added seq back in for writing miRNA pairs (arms)

   chrom_df <- make_si_BamDF(chrom) %>%
     subset(width <= 32 & width >= 18) %>%
     dplyr::mutate(start = pos, end = pos + width - 1) %>%
     dplyr::select(-c(pos)) %>%
     dplyr::group_by_all() %>%
     dplyr::summarize(count = dplyr::n())
   filter_r2_dt <- chrom_df %>% dplyr::mutate(rname = chrom_name)

   if(nrow(filter_r2_dt) == 0){
      return(null_mi_res())
   } else {
      if(weight_reads == "T"){
        r2_dt <- weight_by_prop(filter_r2_dt, chrom_name)
        r1_dt <- weight_by_prop(filter_r2_dt, chrom_name)
      } else if(weight_reads == "Locus_norm" | weight_reads == "locus_norm"){
        locus_length <- reg_stop - reg_start
        locus_read_count <- sum(filter_r2_dt$count)
        r2_dt <- locus_norm(filter_r2_dt, locus_read_count)
        r1_dt <- locus_norm(filter_r2_dt, locus_read_count)
      } else {
        r2_dt <- no_weight(filter_r2_dt, chrom_name)
        r1_dt <- no_weight(filter_r2_dt, chrom_name)
      }
      #transform ends of one set of reads
      r1_dt <- r1_dt %>% dplyr::mutate(end = end + 59)
      filter_r2_dt <- NULL
      if(nrow(r2_dt) == 0){
        return(null_mi_res())
      }
   }

   chrom <- NULL

   if(nrow(r1_dt) == 0 || nrow(r2_dt) == 0){
      cat(paste0(wkdir, logfile), "After filtering for width and strand, zero reads remain. Please check bam BAM file.\n", append = TRUE)
      return(null_mi_res())
   }

   pileup_start <- min(r1_dt$start)
   pileup_stop <- max(r2_dt$end)
   pileups <- get_read_pileups(pileup_start, pileup_stop, bam_scan, bam_file)
   #dt <- pileups %>% dplyr::group_by(pos) %>% dplyr::summarise(count = sum(count))

   empty_table <- data.frame(pos = c(seq(pileup_start, pileup_stop)), count = c(0))

   dt_table <- merge(empty_table, pileups, by = 'pos', all.x = TRUE) %>% dplyr::select(-c(count.x)) %>% dplyr::rename('count' = count.y)
   dt_table[is.na(dt_table)] = 0


   # returns overlaps
   # using find_overlaps instead of find_hp_overlaps
   # find_overlaps only requires one df passed in and doesn't transform end of reads to original
   # find_hp_overlaps requires two dfs and automatically transforms ends of reads
   overlaps_tmp <- find_overlaps(r1_dt, r2_dt)

   overlaps <- overlaps_tmp %>%
      dplyr::mutate(r1_end = r1_end - 59) %>%
      dplyr::mutate(r1_width = r1_end - r1_start + 1)

   overlaps_tmp <- NULL
   #need to filter to remove self matches, get dist between end of R1 and start of R2
   size <- nrow(overlaps)
   overlaps <- get_nearby(overlaps$r1_start, overlaps$r1_end, overlaps$r2_start, overlaps$r2_end, 60, size)

   # only want read pairs where start position of read 2 is after the end of read 1
   overlaps <- overlaps %>%
      dplyr::filter(dist > 2) %>%
      #dplyr::mutate(end_r1 = start_r1 + widthx - 1, end_r2 = start_r2 + widthy - 1) %>%
      dplyr::select(-c(dist))

   # write out pairs here

   if(nrow(overlaps) == 0){
      cat(paste0(wkdir, logfile), "No overlapping reads found.\n", append = TRUE)
      return(null_mi_res())
   }

   # Keep only distinct rows in overlaps, but mutate in a column to keep track of the number of duplicates that were present
   #overlaps_group <- overlaps %>%
   #    dplyr::group_by_all() %>%
   #    dplyr::count()

   read_pileups <- getPileupsMap(dt_table$pos, dt_table$count, overlaps$start_r1, overlaps$end_r1, overlaps$start_r2, overlaps$end_r2)# %>%

   overlaps <- NULL

   if(nrow(read_pileups) == 0){
      return(null_mi_res())
   }

   read_pileups <- read_pileups %>%
      dplyr::mutate(width_r1 = (r1_end - r1_start) + 1, width_r2 = (r2_end - r2_start) + 1)

   size <- nrow(read_pileups)
   r1_df <- data.frame(Start = read_pileups$r1_start, Stop = read_pileups$r1_end)
   r2_df <- data.frame(Start = read_pileups$r2_start, Stop = read_pileups$r2_end)


   ## define coordinates of loop sequence
   ## loop is the sequence between the read 1 and read 2
   loop_coord <- data.frame("r1_start" = r1_df$Start, "r1_stop" = r1_df$Stop, "lstart" = r1_df$Stop + 1,
                            "lstop" = r2_df$Start - 1, "r2_start" = r2_df$Start, "r2_stop" = r2_df$Stop)

   #get rid of rows with loop length less than 3
   loop_coord <- loop_coord %>% dplyr::filter((lstop - lstart + 1) > 2)
   if(nrow(loop_coord) == 0){
      return(null_mi_res())
   }

   #remove results where loop sequence has greater than 5% of total read count
   total_count <- sum(pileups$count)


   df <- loop_coord
   loop_seqs <- getFastas(geno_seq, df$lstart, df$lstop, nrow(df))
   r1_seqs <- getFastas(geno_seq, df$r1_start, df$r1_stop, nrow(df))
   r2_seqs <- getFastas(geno_seq, df$r2_start, df$r2_stop, nrow(df))

   final_coord <- data.frame(start = r1_seqs$start, stop = r2_seqs$stop, r1_seq = r1_seqs$Seq, loop_seq = loop_seqs$Seq, r2_seq = r2_seqs$Seq) %>%
      dplyr::mutate(whole_seq = stringr::str_c(r1_seq, loop_seq, r2_seq)) %>%
      dplyr::select(c(start, stop, whole_seq)) %>%
      dplyr::mutate(width = stop - start + 1) %>% dplyr::distinct()



   ## some results overlapping, get longest result for each read group
   longest_seqs <- as.data.frame(final_coord) %>%
                                    dplyr::mutate(width = (stop - start) + 1) %>%
                                    dplyr::group_by(start) %>%
                                    dplyr::filter(width == max(width)) %>%
                                    dplyr::group_by(stop) %>%
                                    dplyr::filter(width == max(width))
   #merge overlapping results based on position using IRanges
   #get the new sequence using getFastas


   #ir <- IRanges::IRanges(start = longest_seqs$start, end = longest_seqs$stop, width = longest_seqs$width)
   #reduce_ir <- IRanges::reduce(ir)

   longest_seqs$ave_count <- 0
   ## get the read abundance for each set of coordinates
   for(i in 1:nrow(longest_seqs)){
     df <- dt_table %>% dplyr::filter(pos >= longest_seqs$start[i] & pos <= longest_seqs$stop[i])
     longest_seqs$ave_count[i] <- round(mean(df$count))

   }


   if(nrow(longest_seqs > 1)){
     out <- longest_seqs %>% dplyr::mutate(chrom = chrom_name) %>%
        dplyr::select(c(chrom, start, stop, ave_count))
     write.table(out, file = paste0(wkdir, "alt_miRNAs_coord.bed"), sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
     cat(file = paste0(wkdir,logfile), "Writing potential alternative miRNA start and stop coordinates to alt_miRNAs_coord.bed.", append = TRUE)
   }

   # select the result with the highest read abundance
   most_abundant <- which(longest_seqs$ave_count == max(longest_seqs$ave_count))

   #if potential miRNAs have tied abundances, choose first
   if(length(most_abundant) > 1){
     most_abundant <- most_abundant[1]
   }

   # this nomenclature and it's related apply function should be changed in the future
   # used to iterate over multiple rows of df, but changed to select most abundant read pair
   reduced_df <- longest_seqs[most_abundant,]

   reduced_seqs <- getFastas(geno_seq, reduced_df$start, reduced_df$stop, nrow(reduced_df))

   #RNAfold wants U's not T's, so convert to U
   converted <- list(convertU(reduced_seqs$Seq, nrow(reduced_seqs)))

   converted <- data.frame('V1' = unname(unlist(converted)))
   reduced_df$converted <- converted$V1

   geno_seq <- NULL

   ## calculate the pileups for each read for the coverage plots
   #for the coverage plot
   apply_outer <- function(i) {
      cov_vec <- vector(length = 0L)
      print(paste0("start: ", reduced_df$start[i], " stop: ", reduced_df$stop[i]))
      seq_rng <- seq(reduced_df$start[i], reduced_df$stop[i])
      print(paste0("length seq_rng: ", length(seq_rng)))

      # Define the inner loop vector to iterate over
      cov_vec <- sapply(seq_rng, apply_inner)
      print(paste0("length cov_vec: ", length(cov_vec)))

      list(seq_rng, cov_vec)
   }

   apply_inner <- function(j) {
      cur_vec <- which(dt_table$pos == j)
      if (length(cur_vec) > 0){
         tmp_cov <- dt_table$count[cur_vec]
      } else {
         tmp_cov <- 0
      }
      return(tmp_cov)
   }

   # Define the outer loop vector to iterate over
   seq_pileup <- lapply(1:nrow(reduced_df), apply_outer)

   print('making fold_list')
   ################################################################################################################
   #fold_short_rna folds a list of sequences whereas fold_long_rna only folds one
   fold_list <- mapply(fold_short_rna, reduced_df$start, reduced_df$stop, reduced_df$converted, path_to_RNAfold)

   #extracts the needed fields from the fold list, store in more compact list
   make_reduced_list <- function(x){
      start <- fold_list[[x]][[2]]
      stop <- fold_list[[x]][[3]]
      mfe <- fold_list[[x]][[1]]
      vienna <- fold_list[[x]][[5]]
      helix <- R4RNA::viennaToHelix(vienna)
      converted <- fold_list[[x]][[6]]
      extracted_df <- fold_list[[x]][[4]]
      return(list(start = start, stop = stop, mfe = mfe, vienna = vienna, helix = helix, converted = converted, extracted_df = extracted_df))
   }

   reduced_list <- lapply(1:length(fold_list), make_reduced_list)


   #make the plots for all the sequences in the "reduced_list"
   plot_rna_struct <- function(x){
      pos <- count <- count.x <- NULL
      prefix <- paste0(wkdir, chrom_name, "-", reduced_list[[x]]$start, "-", reduced_list[[x]]$stop)
      mfe <- reduced_list[[x]]$mfe
      print(prefix)

      # put the helix table into a data.frame
      helix_df <- data.frame(reduced_list[[x]]$helix)

      perc_paired <- (length(reduced_list[[x]]$helix$i)*2)/(reduced_list[[x]]$stop - reduced_list[[x]]$start)

      # transforms reads from one arm of hairpin to their paired position
      # makes a table of reads which are overlapping
      dicer_overlaps <- dicer_overlaps(r2_dt, reduced_list[[x]]$helix, chrom_name, reduced_list[[x]]$start)

      # summarize the counts by the # overlapping nucleotides
      z_res <- make_count_table(r1_dt$start, r1_dt$end, r1_dt$width, r2_dt$start, r2_dt$end, r2_dt$width)

      # make_count_table was originally written for piRNAs. Need to subtract 3 from each overlap size.
      z_res <- z_res %>% dplyr::mutate(overlap = overlap - 3)
      # create empty z_df
      z_df <- data.frame("Overlap" = z_res[ ,1], "Z_score" = calc_zscore(z_res$count))

      # calculate the zscores, if there are results
      if(is.na(dicer_overlaps[1,1]) | dicer_overlaps[1,1] == 0) {
         overhangs <- data.frame(shift = c(-4,-3,-2,-1,0,1,2,3,4), proper_count = c(0,0,0,0,0,0,0,0,0), improper_count = c(0,0,0,0,0,0,0,0,0))
         overhangs$zscore <- calc_zscore(overhangs$proper_count)
      } else {

         print('making overhangs')
         overhangs <- data.frame(calc_overhangs(dicer_overlaps$r1_start, dicer_overlaps$r1_end,
                                     dicer_overlaps$r2_start, dicer_overlaps$r2_width))
         overhangs$zscore <- calc_zscore(overhangs$proper_count)
      }

      # transform data frame from table to row
      overhang_output <- data.frame(t(overhangs$zscore))
      colnames(overhang_output) <- overhangs$shift
      overhang_output$locus <- paste0(chrom_name, "_", reduced_list[[x]]$start, "_", reduced_list[[x]]$stop)
      overhang_output <- overhang_output[, c(10, 1:9)]

      dice_file <- paste0(wkdir, "miRNA_dicerz.txt")
      suppressWarnings(
         if(!file.exists(dice_file)){
            utils::write.table(overhang_output, file = dice_file, sep = "\t", quote = FALSE, append = T, col.names = T, na = "NA", row.names = F)
         } else {
            utils::write.table(overhang_output, file = dice_file, quote = FALSE, sep = "\t", col.names = F, append = TRUE, na = "NA", row.names = F)
         }
      )


      if(plot_output == "T"){
         #make the plots
         dicer_sig <- plot_overhangz(overhangs, "+")
         print(overhangs)
         #make new pileups dt for structure

         # get the per-base coverage
         # returns a two column df with pos and coverage
         new_pileups <- get_read_pileups(reduced_list[[x]]$start, reduced_list[[x]]$stop, bam_scan, bam_file)  %>%
            dplyr::group_by(pos) %>% dplyr::summarise(count = sum(count))

         # make a table with the same positions but empty cov column for combining with pileups
         # necessary because the pileups table doesn't have all positions from the first nt to the last because
         # coverages of zero aren't reported
         # set these to zero below
         empty_table <- data.frame(pos = c(seq(reduced_list[[x]]$start, reduced_list[[x]]$stop)), count = c(0))


         struct_table <- merge(empty_table, new_pileups, by = 'pos', all.x = TRUE) %>%
            dplyr::select(-c(count.x)) %>% dplyr::rename('count' = count.y)
         #set NA to zero
         struct_table[is.na(struct_table)] = 0

         # the pileup df is the whole region of interest, so subset the df to match
         # the coordinates of the reduced_list (potential miRNAs)
         seq_cov_dat <- subset(struct_table, pos >= reduced_list[[x]]$start & pos <= reduced_list[[x]]$stop)

         seq_cov_dat$count <- seq_cov_dat$count + 1
         seq_cov_dat$logcount <- log2(seq_cov_dat$count)
         # assign colors to the coverage values
         dat_with_colorvals <- data.table::as.data.table(colorNtByDepth(seq_cov_dat))

         #######################################################################################
         graphics::par(mar = c(1,1,1,1))
         #print(reduced_list)

         # assign the color values to each nucloeotide in the predicted structure
         # this function is obscenely complicated
         # creates a 5-line text output for plotting similar to miRBase plots
         miRNA_list <- process_miRNA_struct(reduced_list[[x]]$vienna, reduced_list[[x]]$converted)


         # function to place dashes in the correct place
         insert_vector <- function(input, position, values) {
            options(warn=-1)
            #Create result vector
            res <- numeric(length(input) + length(values))
            res[position] <- values
            #Create an index of input vector
            inds1 <- seq_along(res)
            inds2 <- !(inds1 %in% position)

            #Insert input vector
            res[inds2] <- input
            return(res)
         }

         ### replace paired bases with their respective nucleotide and unpaired bases with "-"
         ### get positions where there is only a space
         ### combine color values with this data

         top <- unname(unlist(miRNA_list$top[3]))
         mid_top <- unname(unlist(miRNA_list$mid_top[3]))
         mid_bottom <- unname(unlist(miRNA_list$mid_bottom[3]))
         bottom <- unname(unlist(miRNA_list$bottom[3]))
         char <- unname(unlist(miRNA_list$char_df[3]))

         #get the positions of paired bases on each line
         top_pos <- which(top == "U" | top == "C" | top == "G" | top == "A")

         mid_top_pos <- which(mid_top == "U" | mid_top == "C" | mid_top == "G" | mid_top == "A")

         char_pos <- which(char == "U" | char == "C" | char == "G" | char == "A")

         bottom_pos <- which(bottom == "U" | bottom == "C" | bottom == "G" | bottom == "A")

         mid_bottom_pos <- which(mid_bottom == "U" | mid_bottom == "C" | mid_bottom == "G" | mid_bottom == "A")

         # get the position of unpaired bases on each line that has them
         top_dash_pos <- which(top == "-")
         bottom_dash_pos <- which(bottom == "-")

         # reverse the bottom half of the data because they will be reversed!
         rev_color_dat <- rev(dat_with_colorvals$bin)
         # this is the top half
         color_half1 <- dat_with_colorvals$bin[1:(length(mid_top_pos) + length(top_pos))]
         # this is the bottom half
         color_half2 <- rev_color_dat[1:(length(mid_bottom_pos) + length(bottom_pos))]

         # put the dashes in their place
         color_half1 <- insert_vector(color_half1, top_dash_pos, rep("-", length(top_dash_pos)))
         color_half2 <- insert_vector(color_half2, bottom_dash_pos, rep("-", length(bottom_dash_pos)))

         # get the number of unique expression values so we know how many colors to use
         num_uniq <- length(unique(dat_with_colorvals$bin)) + 1

         # use the positions gotten earler to create a vector containing the color by position
         mid_color_vec <- dat_with_colorvals$bin[char_pos]

         top_vec <- numeric(length(top))
         top_vec[top_pos] <- color_half1[top_pos]
         top_vec[(top_vec == 0)] <- num_uniq

         mid_top_vec <- numeric(length(mid_top))
         mid_top_vec[mid_top_pos] <- color_half1[mid_top_pos]
         mid_top_vec[(mid_top_vec == 0)] <- num_uniq

         char_vec <- numeric(length(char))
         char_vec[char_pos] <- mid_color_vec
         char_vec[(char_vec == 0)] <- num_uniq

         mid_bottom_vec <- numeric(length(mid_bottom))
         mid_bottom_vec[mid_bottom_pos] <- color_half2[mid_bottom_pos]
         mid_bottom_vec[(mid_bottom_vec == 0)] <- num_uniq

         bottom_vec <- numeric(length(bottom))
         bottom_vec[bottom_pos] <- color_half2[bottom_pos]
         bottom_vec[(bottom_vec == 0)] <- num_uniq


         # now put them all together in the correct order into a data frame for plotting
         top_df <- data.frame(x = unname(unlist(miRNA_list$top[1])), y = unname(unlist(miRNA_list$top[2])), text = unname(unlist(miRNA_list$top[3])),
                              color_bins = top_vec)

         mid_top_df <- data.frame(x = unname(unlist(miRNA_list$mid_top[1])), y = unname(unlist(miRNA_list$mid_top[2])), text = unname(unlist(miRNA_list$mid_top[3])),
                              color_bins = mid_top_vec)
         char_df <- data.frame(x = unname(unlist(miRNA_list$char[1])), y = unname(unlist(miRNA_list$char[2])), text = unname(unlist(miRNA_list$char[3])),
                              color_bins = char_vec)
         mid_bottom_df <- data.frame(x = unname(unlist(miRNA_list$mid_bottom[1])), y = unname(unlist(miRNA_list$mid_bottom[2])), text = unname(unlist(miRNA_list$mid_bottom[3])),
                              color_bins = mid_bottom_vec)
         bottom_df <- data.frame(x = unname(unlist(miRNA_list$bottom[1])), y = unname(unlist(miRNA_list$bottom[2])), text = unname(unlist(miRNA_list$bottom[3])),
                              color_bins = bottom_vec)

         # rbind the results
         # this results in a 4-column df with x-coordinate, y-coordinate (chosen by me), the text character to use, and the color bin
         miRNA_df <- rbind(top_df, mid_top_df, char_df, mid_bottom_df, bottom_df)
         ########################################################################################
         struct_plot <- plot_miRNA_struct(miRNA_df)

         density <- read_densityBySize(bam_obj, chrom_name, reg_start, reg_stop, bam_file, wkdir)
         density_plot <- plot_density(density, reg_start, reg_stop)

         dist_plot <- plot_sizes(read_dist)

         zplot <- plot_overlapz(z_df)

         if(reg_stop - reg_start <= 150){
           #if miRNA is long, use landscape orientation so the struct plot isn't squished.
           left_top <- cowplot::plot_grid(dist_plot, dicer_sig, ncol = 1, rel_widths = c(1,1), rel_heights = c(1,1), align = "vh", axis = "lrtb")
           right_top <- cowplot::plot_grid(NULL, density_plot,zplot, ncol = 1, rel_widths = c(1,1,1), rel_heights = c(0.4,1,1), align = "vh", axis = "lrtb")
           bottom <- cowplot::plot_grid(struct_plot)
           top <- cowplot::plot_grid(left_top, right_top, rel_heights = c(1,1), rel_widths = c(1,1), align = "vh", axis = "lrtb")


           all_plot <- cowplot::plot_grid(top,NULL, bottom, rel_widths = c(1,1,1), ncol = 1, rel_heights = c(1,0.1,0.5), align = "vh", axis = "lrtb")

           if(out_type == "png" || out_type == "PNG"){
             grDevices::png(file = paste0(prefix,"_", strand, "_combined.png"), height = 11, width = 8, units = "in", res = 300)
             print(all_plot)
             grDevices::dev.off()
           } else {
             grDevices::pdf(file = paste0(prefix,"_", strand, "_combined.pdf"), height = 11, width = 8)
             print(all_plot)
           }

          } else {
            top <- cowplot::plot_grid(dist_plot, dicer_sig, density_plot, zplot, ncol = 4, rel_widths = c(1,1,1.4,1), align = "vh", axis = "lrtb")
            #left_top <- cowplot::plot_grid(dist_plot, dicer_sig, ncol = 1, rel_widths = c(1,1), rel_heights = c(1,1), align = "vh", axis = "lrtb")
             #right_top <- cowplot::plot_grid(NULL, density_plot,zplot, ncol = 1, rel_widths = c(1,1,1), rel_heights = c(0.4,1,1), align = "vh", axis = "lrtb")
             bottom <- cowplot::plot_grid(struct_plot)
             #top <- cowplot::plot_grid(left_top, right_top, rel_heights = c(1,1), rel_widths = c(1,1), align = "vh", axis = "lrtb")
             all_plot <- cowplot::plot_grid(top,NULL, bottom, rel_widths = c(1,1,0.6), ncol = 1, rel_heights = c(1,0.1,0.8), align = "vh", axis = "lrtb")


           if(out_type == "png" || out_type == "png"){
            grDevices::png(file = paste0(prefix,"_", strand, "_combined.png"), height = 7, width = 15, units = "in", res = 300)
            print(all_plot)
            grDevices::dev.off()
           } else {
            grDevices::pdf(file = paste0(prefix,"_", strand, "_combined.pdf"), height = 7, width = 15)
            print(all_plot)
            grDevices::dev.off()
           }

         }

      }
      return(list("mfe" = mfe, "perc_paired" = perc_paired, "overhangs" = c(overhangs,z_df)))
   }

   results <- lapply(1:length(reduced_list), plot_rna_struct)

   return(results)
}
