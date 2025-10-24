# Various helper functions for working with and validating gtfs

# For the purposes of MiSiPi.RNA, certain aspects of gtf files need to be validated
#   in order to ensure proper annotation plotting. The files do not need to be
#   validated completely in order to work with this package. This is just a
#   stripped down validation to ensure proper functionality

# Assumes that the file exists as checked by set_vars prior to calling

# Input: gtf_file [character] - A path to a gtf file
# Return: [list] - A list containing 3 objects
#                - 1. is_valid [logical] - Is the gtf valid for the purposes of MiSiPi.RNA?
#                - 2. msg [character] - A message detailing why the gtf is not valid (if it is not valid)
#                - 3. gtf [data.frame] - A data frame containing the read in gtf and filtered for exons, genes, and transcript features

.validate_gtf <- function(gtf_file) {
  
  # Check that it has a consistent number of columns and return early if not
  gtf_columns_vector <- utils::count.fields(gtf_file, sep = "\t", quote = "")
  
  gtf_columns_consistent <- length(unique(gtf_columns_vector)) == 1
  if (!gtf_columns_consistent) {
    msg <- "The gtf file does not have the same number of columns throughout."
    return(list(is_valid = FALSE, msg = msg, gtf = NULL))
  }
  
  # Check that it has 9 columns
  # Just checking the first row here since they are all the same at this point
  gtf_num_columns <- gtf_columns_vector[1] == 9
  if (gtf_columns_consistent & !gtf_num_columns) {
    msg <- glue::glue("The gtf file appears to have {gtf_columns_vector[1]} columns but should have 9.")
    return(list(is_valid = FALSE, msg = msg, gtf = NULL))
  }
  
  # Now that initial validation is complete, read in the file to memory
  # and filter for exons, genes, and transcripts
  gtf <- readr::read_tsv(
    file = gtf_file,
    comment = "#",
    col_names = c("seqname","source","feature","start","end","score","strand","frame","attribute"),
    show_col_types = FALSE
  ) %>%
    dplyr::filter(feature == "exon" | feature == "gene" | feature == "transcript")
  
  # Remove quotation marks
  gtf$attribute <- stringr::str_remove_all(gtf$attribute, "\"")
  
  # Check if any data remains after filtering
  gtf_has_rows <- nrow(gtf) > 0
  if (!gtf_has_rows) {
    msg <- "After filtering for exons, genes, and transcripts, the gtf file has no observations."
    return(list(is_valid = FALSE, msg = msg, gtf = NULL))
  }
  
  # For exon, gene, and transcript feature observations:
  # Check that the attribute column contains gene_id's
  
  # The pattern (^|;\\s*)key\\s works like this (using matching groups):
  # Match for key plus a space (gene_id or transcript_id in this case) either at the beginning of the string,
  #   or match for key plus a space after a semicolon (option space allowed here)
  #   Since the attribute column is formatted Key1 Value1; Key2 Value2; etc
  #   This will quickly search for the key given in each Key Value pair.
  has_gene_ids <- all(stringi::stri_detect_regex(gtf$attribute, "(^|;\\s*)gene_id\\s"))
  
  # For exon and transcript feature observations:
  # Check that the attribute column contains transcript_id's
  sans_genes <- gtf %>%
    dplyr::filter(feature != "gene")
  
  has_transcript_ids <- all(stringi::stri_detect_regex(sans_genes$attribute, "(^|;\\s*)transcript_id\\s"))
  
  sans_genes <- NULL
  
  if (!has_gene_ids | !has_transcript_ids) {
    
    if (!has_gene_ids & !has_transcript_ids) {
      msg <- "Attribute field is missing gene_id and transcript_id in one or more observations."
    } else if (!has_gene_ids) {
      msg <- "Attribute field is missing gene_id in one or more observations."
    } else {
      msg <- "Attribute field is missing transcript_id in one or more observations."
    }
    
    return(list(is_valid = FALSE, msg = msg, gtf = NULL))
  }
  
  return(list(is_valid = TRUE, msg = NULL, gtf = gtf))
}

.get_gene_id <- Vectorize(function(x) {
  tmp <- unlist(strsplit(x, ";"))
  
  # Locate gene_id in the attributes column
  gene_idx <- grep("gene_id", tmp)
  gene_field <- tmp[gene_idx]
  
  # Strip out the text gene_id
  gene_id <- gsub("gene_id ", "", gene_field)
  
  # In case gene_id wasn't the first attribute, strip another space
  gene_id <- gsub(" ", "", gene_id)
  
  return(gene_id)
})

.get_transcript_id <- Vectorize(function(x) {
  tmp <- unlist(strsplit(x, ";"))
  
  # Locate transcript_id in the attributes column
  transcript_idx <- grep(" transcript_id", tmp)
  
  # If transcript_id isn't in attributes column, set to NA
  if (identical(transcript_idx, integer(0))) {
    return(NA)
  }
  
  transcript_field <- tmp[transcript_idx]
  
  # Strip out the text transcript_id
  transcript_id <- gsub(" transcript_id ", "", transcript_field)
  # In case transcript_id wasn't the first attribute, strip another space
  # transcript_id <- gsub(" ", "", transcript_id)
  
  # If no actual transcript_id exists, set to NA
  # This would occur when the attribute name "transcript_id" is present but no ID follows it
  if (transcript_id == "") {
    transcript_id <- NA
  }
  
  return(transcript_id)
})


