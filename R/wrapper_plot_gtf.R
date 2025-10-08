# takes a gtf file
# outputs a plot
# @param gtf_file a string path leading to the gtf file
# @param chrom_name a string
# @param reg_start an integer
# @param reg_stop an integer
# @param logfile
#
# @return plot

.wrapper_plot_gtf <- function(gtf_file, chrom_name, reg_start, reg_stop, logfile) {
  # Load GTF
  gtf <- read.csv(gtf_file, header = FALSE, sep = "\t", fill = TRUE, quote="")
  colnames(gtf)[1:9] <- c("seqname","source","feature","start","end","score","strand","frame","attribute")
  
  # Filter for region
  gtf <- gtf %>% 
    # filter features that 
    dplyr::filter(seqname == chrom_name, start >= reg_start & start <= reg_stop | end >= reg_start & end <= reg_stop)  %>%
    dplyr::filter(feature != "start_codon" & feature != "CDS" & feature != "stop_codon")
  
  if (nrow(gtf) == 0) {
    message("No features found in region.")
    return(NULL)
  }
  
  # Remove quotation marks
  gtf$attribute <- gsub("\"", "", gtf$attribute)
  
  exons_present <- any(gtf$feature == "exon")
  transcripts_present <- any(gtf$feature == "transcript")
  genes_present <- any(gtf$feature == "gene")
  
  if (exons_present & !transcripts_present & !genes_present) { # exons only
    gtf_plot <- .plot_exons_only(gtf, reg_start, reg_stop)
  } else if (exons_present & transcripts_present & !genes_present) { # exons and transcripts
    gtf_plot <- .plot_transcripts_exons(gtf, reg_start, reg_stop)
  } else if (!exons_present & transcripts_present & genes_present) { # genes and transcripts
    gtf_plot <- .plot_genes_transcripts(gtf, reg_start, reg_stop)
  } else if (exons_present & transcripts_present & genes_present) { # exons and transcripts and genes
    gtf_plot <- .plot_genes_exon_transcripts(gtf, reg_start, reg_stop)
  } else if (!exons_present & transcripts_present & !genes_present) { # transcripts only
    gtf_plot <- .plot_transcripts_only(gtf, reg_start, reg_stop)
  } else if (!exons_present & !transcripts_present & genes_present) { # genes only
    gtf_plot <- .plot_genes_only(gtf, reg_start, reg_stop)
  } else { # exons and genes only
    gtf_plot <- .plot_genes_exons(gtf, reg_start, reg_stop)
  }
  
  return(gtf_plot)
}
  
  
  
  
  
  
  
  
  
  
  
  
  
  