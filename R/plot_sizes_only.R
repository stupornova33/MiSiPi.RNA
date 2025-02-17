# plot_sizes_only
# A function which will process the reads at a locus but only plot the size distribution and will not run the miRNA, siRNA or piRNA modules.
# @param wkdir A string - the path of the output directory.
# @param geno A string specifying the genotype of the sample, for appending as a suffix to the output file.
# @param bam_file A string - the name of the BAM file.
# @param bed_file A string - the name of the 3-column BED file.
# @return results

.plot_sizes_only <- function(wkdir, geno, bam_file, bed_file) {
  roi <- read.csv(bed_file, header = FALSE, sep = "\t")

  for (i in 1:nrow(roi)) {
    chrom_name <- roi$V1[i]
    reg_start <- roi$V2[i]
    reg_stop <- roi$V3[i]
    logfile <- "plot_sizes.log.txt"
    prefix <- .get_region_string(chrom_name, reg_start, reg_stop)

    # use Rsamtools to process the bam file
    bam_obj <- .open_bam(bam_file, logfile)
    bam_header <- Rsamtools::scanBamHeader(bam_obj)
    chr_name <- names(bam_header[["targets"]])
    chr_length <- unname(bam_header[["targets"]])
    bam_header <- NULL

    # extract reads by strand
    # this creates a list object
    chromP <- .get_chr(bam_obj, chrom_name, reg_start, reg_stop, strand = "plus")
    chromM <- .get_chr(bam_obj, chrom_name, reg_start, reg_stop, strand = "minus")

    forward_dt <- data.table::setDT(.make_si_BamDF(chromP)) %>%
      subset(width <= 32 & width >= 18) %>%
      dplyr::rename(start = pos) %>%
      dplyr::mutate(end = start + width - 1) %>%
      dplyr::group_by_all() %>%
      # get the number of times a read occurs
      dplyr::summarize(count = dplyr::n())

    reverse_dt <- data.table::setDT(.make_si_BamDF(chromM)) %>%
      subset(width <= 32 & width >= 18) %>%
      dplyr::rename(start = pos) %>%
      dplyr::mutate(end = start + width - 1) %>%
      dplyr::group_by_all() %>%
      dplyr::summarize(count = dplyr::n())

    dist <- .get_weighted_read_dist(forward_dt, reverse_dt)

    # size_plot <- .plot_sizes(dist)

    # grDevices::pdf(file = paste0(wkdir, prefix, "_sizes.pdf"), height = 4, width = 4)
    # print(size_plot)
    # grDevices::dev.off()

    ###########################################################################################################################

    # plot histogram with AUC
    forward_dt <- .weight_reads(forward_dt, weight_reads = "none", 0L, 0L) %>%
      dplyr::mutate(width = end - start + 1)
    
    reverse_dt <- .weight_reads(reverse_dt, weight_reads = "none", 0L, 0L) %>%
      dplyr::mutate(width = end - start + 1)

    all_data <- rbind(forward_dt, reverse_dt)
    all_data <- na.omit(all_data)

    if (nrow(all_data) > 0) {
      png(paste0(wkdir, prefix, "_", geno, "_density.png"), height = 5, width = 5, units = "in", res = 300)
      d <- density(all_data$width, bw = 0.4, kernel = "biweight")
      # auc <- sum(diff(d$x) * (head(d$y,-1)+tail(d$y,-1)))/2
      auc <- MESS::auc(x = d$x, y = d$y, from = 18, to = 32, type = "spline")
      plot(d, main = "Density", sub = paste0("AUC = ", auc), xlim = c(18, 32), ylim = c(0, 1))
      dev.off()
    } else {
      cat(paste0(geno, " ", prefix, "had no 0 reads.\n"), file = logfile, append = TRUE)
      next
    }
  }
}
