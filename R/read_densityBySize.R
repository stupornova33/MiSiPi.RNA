#' function to filter reads by size and plot pileup density
#' @param bam_obj a string
#' @param chrom_name a string
#' @param reg_start a string
#' @param reg_stop a string
#' @param input_file a string
#' @param dir a string
#' @param logfile a string
#' @return all_df

#' @export

read_densityBySize <- function(bam_obj, chrom_name, reg_start, reg_stop, input_file, dir, logfile){
   
   pos <- width <- rname <- NULL
   filter_bamfile <- function(input_file, size1, size2,  strand){
     seqnames <- NULL
      which <- GenomicRanges::GRanges(seqnames=chrom_name, IRanges::IRanges(reg_start, reg_stop))
      filters <- S4Vectors::FilterRules(list(MinWidth=function(x) (BiocGenerics::width(x$seq) >= size1 & BiocGenerics::width(x$seq) <= size2)))
      if(strand == "+"){
         filename <- paste0(dir, size1, "_", size2, "_pos.bam")
         Rsamtools::filterBam(input_file, destination = filename, filter = filters, 
                        param = Rsamtools::ScanBamParam(flag = Rsamtools::scanBamFlag(isMinusStrand = FALSE),what=c('rname', 'pos', 'qwidth', 'seq'),  which = which))
 
      } else{
         filename <- paste0(dir, size1, "_", size2, "_neg.bam")
         Rsamtools::filterBam(input_file, destination = filename, filter = filters, 
                        param = Rsamtools::ScanBamParam(flag = Rsamtools::scanBamFlag(isMinusStrand = TRUE),what=c('rname', 'pos', 'qwidth', 'seq'),  which = which))

   }

      
   new_bam_obj <- OpenBamFile(filename, logfile)
   print(filename)
   
   return(new_bam_obj)
   }
   
   filtered_pos_18_19_bam <- filter_bamfile(input_file, 18, 19, "+")
   filtered_pos_20_22_bam <- filter_bamfile(input_file, 20, 22, "+")
   filtered_pos_23_25_bam <- filter_bamfile(input_file, 23, 25, "+")
   filtered_pos_26_32_bam <- filter_bamfile(input_file, 26, 32, "+")
   #filtered_pos_all_bam <- filter_bamfile(input_file, 18, 32, "+")
   
   filtered_neg_18_19_bam <- filter_bamfile(input_file, 18,19, "-")
   filtered_neg_20_22_bam <- filter_bamfile(input_file, 20, 22, "-")
   filtered_neg_23_25_bam <- filter_bamfile(input_file, 23, 25, "-")
   filtered_neg_26_32_bam <- filter_bamfile(input_file, 26, 32, "-")
   #filtered_neg_all_bam <- filter_bamfile(input_file, 18, 32, "-")
   
  
   
   make_bam_pileup <- function(bam, strand){
      seqnames <- pos <- count <- NULL
      which <- GenomicRanges::GRanges(seqnames=chrom_name, IRanges::IRanges(reg_start, reg_stop))
      if(strand == "-"){
         bam_scan <- Rsamtools::ScanBamParam(flag = Rsamtools::scanBamFlag(isMinusStrand = TRUE), what=c('rname', 'pos', 'qwidth'), which=which)
      } else {
         bam_scan <- Rsamtools::ScanBamParam(flag = Rsamtools::scanBamFlag(isMinusStrand = FALSE), what=c('rname', 'pos', 'qwidth'), which=which)
      }
      
      params <- Rsamtools::PileupParam(max_depth=4000, min_base_quality=20, min_mapq=0, min_nucleotide_depth=1, distinguish_strands=TRUE,
                                    distinguish_nucleotides=TRUE, ignore_query_Ns=TRUE, include_deletions=TRUE, include_insertions=FALSE, left_bins=NULL,
                                    query_bins=NULL, cycle_bins=NULL)  
      pileups <- Rsamtools::pileup(bam, index=(stringr::str_c(input_file, '','.bai')), scanBamParam=bam_scan, pileupParam=params)  %>%
         dplyr::select(-c(seqnames, strand))
      
      #Summarize duplicate positions including minor alleles
      dt <- pileups %>% dplyr::group_by(pos) %>% dplyr::summarise(count = sum(count)) 
      return(dt)
   }
   
   ## 18-19 nt   
   chromP <- getChrPlus(filtered_pos_18_19_bam, chrom_name, reg_start, reg_stop)
  
   chromM <- getChrMinus(filtered_neg_18_19_bam, chrom_name, reg_start, reg_stop)
   
   pos_18_19_dt <- data.table::setDT(makeBamDF(chromP)) %>%
      dplyr::mutate(start = pos, end = pos + width - 1) %>%
      dplyr::select(-c(rname, pos))
   
   neg_18_19_dt <- data.table::setDT(makeBamDF(chromM)) %>%
      dplyr::mutate(start = pos, end = pos + width - 1) %>%
      dplyr::select(-c(rname, pos))
   
   pos_18_19_pileup <- make_bam_pileup(filtered_pos_18_19_bam, "+")
   neg_18_19_pileup <- make_bam_pileup(filtered_neg_18_19_bam, "-")
   
   ##20-22 nt
   chromP <- getChrPlus(filtered_pos_23_25_bam, chrom_name, reg_start, reg_stop)
   
   chromM <- getChrMinus(filtered_pos_23_25_bam, chrom_name, reg_start, reg_stop)
   
   pos_20_22_dt <- data.table::setDT(makeBamDF(chromP)) %>%
      dplyr::mutate(start = pos, end = pos + width - 1) %>%
      dplyr::select(-c(rname, pos))
   
   neg_20_22_dt <- data.table::setDT(makeBamDF(chromM)) %>%
      dplyr::mutate(start = pos, end = pos + width - 1) %>%
      dplyr::select(-c(rname, pos))
   
   pos_20_22_pileup <- make_bam_pileup(filtered_pos_20_22_bam, "+")
   neg_20_22_pileup <- make_bam_pileup(filtered_neg_20_22_bam, "-")
   
   ##23-25nt
   
   chromP <- getChrPlus(filtered_pos_23_25_bam, chrom_name, reg_start, reg_stop)
   
   chromM <- getChrMinus(filtered_neg_23_25_bam, chrom_name, reg_start, reg_stop)
   
   pos_23_25_dt <- data.table::setDT(makeBamDF(chromP)) %>%
      dplyr::mutate(start = pos, end = pos + width - 1) %>%
      dplyr::select(-c(rname, pos))
   
   neg_23_25_dt <- data.table::setDT(makeBamDF(chromM)) %>%
      dplyr::mutate(start = pos, end = pos + width - 1) %>%
      dplyr::select(-c(rname, pos))
   
   pos_23_25_pileup <- make_bam_pileup(filtered_pos_23_25_bam, "+")
   neg_23_25_pileup <- make_bam_pileup(filtered_neg_23_25_bam, "_")
   
   ## 26-32
   chromP <- getChrPlus(filtered_pos_26_32_bam, chrom_name, reg_start, reg_stop)
   
   chromM <- getChrMinus(filtered_neg_26_32_bam, chrom_name, reg_start, reg_stop)
   
   pos_26_32_dt <- data.table::setDT(makeBamDF(chromP)) %>%
      dplyr::mutate(start = pos, end = pos + width - 1) %>%
      dplyr::select(-c(rname, pos))
   
   neg_26_32_dt <- data.table::setDT(makeBamDF(chromM)) %>%
      dplyr::mutate(start = pos, end = pos + width - 1) %>%
      dplyr::select(-c(rname, pos))
   
   pos_26_32_pileup <- make_bam_pileup(filtered_pos_26_32_bam, "+")
   neg_26_32_pileup <- make_bam_pileup(filtered_neg_26_32_bam, "-")
   

   empty_dat <- data.frame(pos = c(seq(reg_start, reg_stop)))
   
   merge_neg_table <- function(empty_dat, res_dat){
      count <- NULL
      neg_res <- merge(empty_dat, res_dat, by = "pos", all.x = TRUE) %>% dplyr::mutate(count = count *-1) %>%
         dplyr::mutate(size = "all")
      neg_res["count"][is.na(neg_res["count"])] <- 0
      return(neg_res)
   }
   
   merge_pos_table <- function(empty_dat, res_dat){
      pos_res <- merge(empty_dat, res_dat, by = "pos", all.x = TRUE) %>%
         dplyr::mutate(size = "all")
      pos_res["count"][is.na(pos_res["count"])] <- 0
      
      return(pos_res)
   }

   pos_18_19_res <- merge_pos_table(empty_dat, pos_18_19_pileup)
   neg_18_19_res <- merge_neg_table(empty_dat, neg_18_19_pileup)
   
   pos_20_22_res <- merge_pos_table(empty_dat, pos_20_22_pileup)
   neg_20_22_res <- merge_neg_table(empty_dat, neg_20_22_pileup)
   
   pos_23_25_res <- merge_pos_table(empty_dat, pos_23_25_pileup)
   neg_23_25_res <- merge_neg_table(empty_dat, neg_23_25_pileup)
   
   pos_26_32_res <- merge_pos_table(empty_dat, pos_26_32_pileup)
   neg_26_32_res <- merge_neg_table(empty_dat, neg_26_32_pileup)
   
   all_df <- data.frame(position = pos_26_32_res$pos, pos_26_32 = pos_26_32_res$count, neg_26_32 = neg_26_32_res$count,
                        pos_23_25 = pos_23_25_res$count, neg_23_25 = neg_23_25_res$count,
                        pos_20_22 = pos_20_22_res$count, neg_20_22 = neg_20_22_res$count,
                        pos_18_19 = pos_18_19_res$count, neg_18_19 = neg_18_19_res$count)
                        
   return(all_df)
   
}
