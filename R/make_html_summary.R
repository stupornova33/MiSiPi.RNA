#' Make_html_summary
#' Takes a file path and creates a graphic summary of all loci from a MiSiPi run
#' @param path_to_tables A string specifying the full path to the folder containing the table outputs from MiSiPi.RNA
#' @param type A string specifying which module of MiSiPi.RNA should be plotted in the summary. Options are one of: "miRNA", "siRNA", or "piRNA"
#' @param ml_plots A string specifying whether probability plots are available in the directory and should be included in the outputs.
#' @return nothing
#' @export

make_html_summary <- function(path_to_tables, type, ml_plots = FALSE) {
  # Process the input data tables and make summary plots
  # need to make an option for pdf or png files
  # All summary files get the ML table and read size distribution plot
  if (!dir.exists(paste0(path_to_tables, "summary_inputs"))) {
    dir.create(paste0(path_to_tables, "summary_inputs"))
  }

  linked_exist <- list.files(path_to_tables, pattern = "linked_*")
  
  suppressMessages(
    if(!identical(linked_exist, character(0))){
      file.remove(paste0(path_to_tables, linked_exist))
    }
  )
  input_dir <- paste0(path_to_tables, "summary_inputs/")

  ml_file <- list.files(path_to_tables, pattern = "_ml.txt")
  ml_tab <- read.csv(paste0(path_to_tables, ml_file), sep = "\t", header = TRUE)
  size_dist_tab <- list.files(path_to_tables, pattern = "_size_distributions.txt")
  size_dist_tab <- read.table(paste0(path_to_tables, size_dist_tab), sep = "\t", header = TRUE)

  RowVar <- function(x, ...) {
    rowSums((x - rowMeans(x, ...))^2, ...) / (dim(x)[2] - 1)
  }

  size_dist_mat <- size_dist_tab[, 2:ncol(size_dist_tab)]
  var <- RowVar(size_dist_mat)
  idx <- which(var == 0)

  if (length(idx) > 0) {
    size_dist_mat <- size_dist_mat[-c(idx), ]
  }

  heatmap <- pheatmap::pheatmap(size_dist_mat, main = "Read size distribution", cluster_cols = FALSE, show_rownames = F, show_colnames = T, 
                                scale = "row", labels_col =  c("16", "", "18", "", "20", "", "22", "", "24", "", "26", "", "28", "", "30", "", "32"))

  ggplot2::ggsave(heatmap, file = paste0(input_dir, "size_dist_heatmap.png"))


  nplots <- vector()
  ################################################################# siRNA ###################################################################
  if (type == "siRNA" || type == "sirna") {
    dicerz <- read.table(paste0(path_to_tables, "siRNA_outputs/", "siRNA_dicerz.txt"), header = TRUE)
    dicerz[is.na(dicerz)] <- -33
    dicerz <- dicerz %>% dplyr::select(-c(locus))
    colnames(dicerz) <- c("-4", "-3", "-2", "-1", "0", "1", "2", "3", "4")

    var <- RowVar(dicerz)
    idx <- which(var == 0)

    if (length(idx) > 0) {
      dicerz <- dicerz[-c(idx), ]
    }

    # if row clustering is to be done, need to remove rows with no variance
    if (nrow(dicerz) > 2) {
      dicer_heat <- pheatmap::pheatmap(dicerz, main = "2nt Overhang Z-score", fontsize_col = 14, cluster_cols = FALSE, show_rownames = FALSE, cluster_rows = TRUE, scale = "row")

      # png(filename=paste0(path_to_tables, "si_dicerz_heatmap.png"))
      # print(dicer_heat)
      # dev.off()
      ggplot2::ggsave(dicer_heat, file = paste0(input_dir, "si_dicerz_heatmap.png"))
      nplots <- append(nplots, paste0(input_dir, "si_dicerz_heatmap.png"))
    } else {
      print("siRNA dicer heat table contained two few loci to make a clustered heatmap.")
    }

    minus_hp_phasedz <- read.table(paste0(path_to_tables, "siRNA_outputs/", "minus_hp_phasedz.txt"), header = TRUE)
    minus_hp_phasedz[is.na(minus_hp_phasedz)] <- -33

    minus_hp_phasedz <- minus_hp_phasedz[,-c(1)]
    minus_hp_phasedz <- minus_hp_phasedz[, 1:10]
    
    var <- RowVar(minus_hp_phasedz)
    idx <- which(var == 0)
    if (length(idx) > 0) {
      minus_hp_phasedz <- minus_hp_phasedz[-c(idx), ]
    }

    if (nrow(minus_hp_phasedz) > 2) {
      minus_hp_heat <- pheatmap::pheatmap(minus_hp_phasedz, main = "Minus hp phased z-scores", fontsize = 10, cluster_cols = FALSE, 
                                          show_rownames = FALSE, cluster_rows = TRUE, scale = "row", 
                                          labels_col = c("1","", "3", "", "5", "", "7", "", "9", ""))
      
      ggplot2::ggsave(minus_hp_heat, file = paste0(input_dir, "minus_hp_phasedz_heat.png"))
      # nplots <- append(nplots, paste0(path_to_tables, 'minus_hp_heatmap.png'))
      nplots <- append(nplots, paste0(input_dir, "minus_hp_phasedz_heat.png"))
    } else {
      print("Minus_hp_phasedz table contained two few loci to make a clustered heatmap.")
    }

    plus_hp_phasedz <- read.table(paste0(path_to_tables, "siRNA_outputs/", "plus_hp_phasedz.txt"), header = TRUE)

    plus_hp_phasedz <- plus_hp_phasedz %>% dplyr::select(-c(V1))
    plus_hp_phasedz <- na.omit(plus_hp_phasedz)
    #colnames(plus_hp_phasedz) <- c(seq(1, 10))
    plus_hp_phasedz <- plus_hp_phasedz[, 1:10]
    var <- RowVar(plus_hp_phasedz)
    idx <- which(var == 0)
    if (length(idx) > 0) {
      plus_hp_phasedz <- plus_hp_phasedz[-c(idx), ]
    }


    # plus_hp_phasedz <- plus_hp_phasedz[,1:10]

    if (nrow(plus_hp_phasedz) > 2) {
      plus_hp_heat <- pheatmap::pheatmap(plus_hp_phasedz, main = "Plus hp phased z-scores", fontsize = 10, cluster_cols = FALSE, show_rownames = FALSE, 
                                         cluster_rows = TRUE, scale = "row",  labels_col = c("1","", "3", "", "5", "", "7", "", "9", ""))
      ggplot2::ggsave(plus_hp_heat, file = paste0(input_dir, "plus_hp_phasedz_heat.png"))

      nplots <- append(nplots, paste0(input_dir, "plus_hp_phasedz_heat.png"))
    } else {
      print("Plus_hp_phasedz table contained two few loci to make a clustered heatmap.")
    }

    widths <- paste(rep("'300px'", times = length(nplots)), collapse = ",")
    heights <- paste(rep("'350px'", times = length(nplots)), collapse = ",")

    ### set fig out heights and widths
    cat_sizes <- paste0("```{r, echo = FALSE, out.width=c('350px',", widths, "),", "out.height = c('400px',", heights, "), fig.show='hold' }\n")

    nplots <- paste(nplots, collapse = '", "')
    ## set cat statements for files
    cat_names <- paste0('knitr::include_graphics(c("', input_dir, 'size_dist_heatmap.png", "', nplots, '"))\n')


    ############################################################### miRNA ##################################################################
  } else if (type == "miRNA" || type == "mirna") {
    mirna_dicerz <- read.table(paste0(path_to_tables, "miRNA_outputs/", "miRNA_dicerz.txt"), header = TRUE) %>% dplyr::select(-c(locus))
    colnames(mirna_dicerz) <- c(seq(-4, 4, by = 1))
    mirna_dicerz[is.na(mirna_dicerz)] <- -33

    var <- RowVar(mirna_dicerz)
    idx <- which(var == 0)
    if (length(idx) > 0) {
      mirna_dicerz <- mirna_dicerz[-c(idx), ]
    }
    if (nrow(mirna_dicerz) > 2) {
      mirna_dicerz_heat <- pheatmap::pheatmap(mirna_dicerz, main = "miRNA Dicer Z-scores", fontsize = 12, cluster_cols = FALSE, show_rownames = FALSE, cluster_rows = TRUE, scale = "row")
      ggplot2::ggsave(mirna_dicerz_heat, file = paste0(input_dir, "mirna_dicerz_heatmap.png"))

      nplots <- append(nplots, paste0(input_dir, "mirna_dicerz_heatmap.png"))
    } else {
      print("miRNA dicerz table contained two few loci to make a clustered heatmap.")
    }

    widths <- paste(rep("'300px'", times = length(nplots)), collapse = ",")
    heights <- paste(rep("'350px'", times = length(nplots)), collapse = ",")

    ### set fig out heights and widths
    cat_sizes <- paste0("```{r, echo = FALSE, out.width=c('350px',", widths, "),", "out.height = c('400px',", heights, "), fig.show='hold' }\n")

    nplots <- paste(nplots, collapse = '", "')
    ## set cat statements for files
    cat_names <- paste0('knitr::include_graphics(c("', input_dir, 'size_dist_heatmap.png", "', nplots, '"))\n')

    ################################################################# piRNA #####################################################################
  } else {
    all_pirna_phasedz <- read.table(paste0(path_to_tables, "piRNA_outputs/", "all_phased_piRNA_zscores.txt"), header = TRUE) %>% dplyr::select(-c(locus))

    all_pirna_phasedz[is.na(all_pirna_phasedz)] <- -33
    var <- RowVar(all_pirna_phasedz)
    idx <- which(var == 0)

    if (length(idx) > 0) {
      all_pirna_phasedz <- all_pirna_phasedz[-c(idx), ]
    }

    if (nrow(all_pirna_phasedz) > 2) {
      pi_phased_heat <- pheatmap::pheatmap(all_pirna_phasedz, main = "phased piRNA Z-scores", fontsize_col = 12, 
                                           cluster_cols = FALSE, show_rownames = FALSE, cluster_rows = TRUE, scale = "row",
                                           labels_col = c(
                                             "0", "", "", "3", "", "", "6", "", "", "9", "", "", "12", "", "", "15", "", "", "18", "", "", "21", "", "", "24", "", "",
                                             "27", "", "", "30", "", "", "33", "", "", "36", "", "", "39", "", "", "42", "", "", "45", "", "", "48", "", ""
                                           ))

      ggplot2::ggsave(pi_phased_heat, file = paste0(input_dir, "pi_phased_heat.png"))

      nplots <- append(nplots, paste0(input_dir, "pi_phased_heat.png"))
    } else {
      print("piRNA phasedz table contained two few loci to make a clustered heatmap.")
    }

    max_piz_overlap <- read.table(paste0(path_to_tables, "piRNA_outputs/", "piRNA_alloverlaps_counts.txt"), header = TRUE)

    colnames(max_piz_overlap) <- c(
      "4", "", "6", "", "8", "", "10", "", "12", "", "14", "", "16", "", "18", "", "20",
      "", "22", "", "24", "", "26", "", "28", "", "30"
    )

    max_piz_overlap[is.na(max_piz_overlap)] <- -33
    var <- RowVar(max_piz_overlap)
    idx <- which(var == 0)

    if (length(idx) > 0) {
      max_piz_overlap <- max_piz_overlap[-c(idx), ]
    }

    if (nrow(max_piz_overlap) > 0) {
      pi_overlap <- pheatmap::pheatmap(max_piz_overlap, main = "Overlap Counts", fontsize_col = 12, cluster_cols = FALSE, 
                                       show_rownames = FALSE, cluster_rows = TRUE, scale = "row",
                                       labels_col = c("4", "", "6", "", "8", "", "10", "", "12", "", "14", "", "16", "", "18", "", "20",
                                         "", "22", "", "24", "", "26", "", "28", "", "30"))
      ggplot2::ggsave(pi_overlap, file = paste0(input_dir, "pi_overlapz_heat.png"))

      # nplots <- append(nplots, paste0(path_to_tables, "max_piz_overlap_heat.png"))
      nplots <- append(nplots, paste0(input_dir, "pi_overlapz_heat.png"))
    } else {
      print("max_piz_overlap contained two few loci to make a clustered heatmap.")
    }

    minus_phasedz <- read.table(paste0(path_to_tables, "piRNA_outputs/", "phased_minus_piRNA_zscores.txt"), header = FALSE) %>% dplyr::select(c(V2, V3, V4, V5))

    minus_phasedz[is.na(minus_phasedz)] <- -33

    colnames(minus_phasedz) <- c("0nt", "1nt", "2nt", "3nt")

    var <- RowVar(minus_phasedz)
    idx <- which(var == 0)

    if (length(idx) > 0) {
      minus_phasedz <- minus_phasedz[-c(idx), ]
    }

    if (nrow(minus_phasedz) > 2) {
      minus_phasedz_heat <- pheatmap::pheatmap(minus_phasedz, main = "Phasing Z-scores (minus strand)", border_color = "NA", fontsize_col = 12, cluster_cols = FALSE, show_rownames = FALSE, cluster_rows = TRUE, scale = "row")

      ggplot2::ggsave(minus_phasedz_heat, file = paste0(input_dir, "minus_phasedz_heat.png"))

      nplots <- append(nplots, paste0(input_dir, "minus_phasedz_heat.png"))
    } else {
      print("piRNA minus_phasedz table contained two few loci to make a clustered heatmap.")
    }


    plus_phasedz <- read.table(paste0(path_to_tables, "piRNA_outputs/", "phased_plus_piRNA_zscores.txt"), header = TRUE, sep = "\t") %>%
      dplyr::select(V2, V3, V4, V5)
    plus_phasedz[is.na(plus_phasedz)] <- -33
    colnames(plus_phasedz) <- c("0nt", "1nt", "2nt", "3nt")

    var <- RowVar(plus_phasedz)
    idx <- which(var == 0)

    if (length(idx) > 0) {
      plus_phasedz <- plus_phasedz[-c(idx), ]
    }

    if (nrow(plus_phasedz) > 0) {
      plus_phasedz_heat <- pheatmap::pheatmap(plus_phasedz, main = "Phasing Z-scores (plus strand)", fontsize_col = 12, border_color = "NA", cluster_cols = FALSE, show_rownames = FALSE, cluster_rows = TRUE, scale = "row")
      ggplot2::ggsave(plus_phasedz_heat, file = paste0(input_dir, "plus_phasedz_heatmap.png"))

      nplots <- append(nplots, paste0(input_dir, "plus_phasedz_heatmap.png"))
    } else {
      print("piRNA plus_phasedz table contained two few loci to make a clustered heatmap.")
    }
    ### set fig out heights and widths

    widths <- paste(rep("'250px'", times = length(nplots)), collapse = ",")
    heights <- paste(rep("'300px'", times = length(nplots)), collapse = ",")

    ### set fig out heights and widths
    cat_sizes <- paste0("```{r, echo = FALSE, out.width=c('300px',", widths, "),", "out.height = c('350px',", heights, "), fig.show='hold' }\n")

    nplots <- paste(nplots, collapse = '", "')
    ## set cat statements for files
    cat_names <- paste0('knitr::include_graphics(c("', input_dir, 'size_dist_heatmap.png", "', nplots, '"))\n')
  }

  ########################################### Find and link the type-specific plots to the ML table ####################################################

  ml_tab$mi_plus_col <- ""
  ml_tab$mi_minus_col <- ""
  ml_tab$si_col <- ""
  ml_tab$pi_col <- ""
  ml_tab$prob_col <- ""

  ml_tab <- ml_tab %>% dplyr::select(c(locus, mi_plus_col, mi_minus_col, pi_col, si_col, prob_col), 7:ncol(ml_tab))

  si_files <- list.files(paste0(path_to_tables, "siRNA_outputs/"), pattern = "_si_plot")
  mi_plus_files <- list.files(paste0(path_to_tables, "miRNA_outputs/"), pattern = "_\\+_combined")
  mi_minus_files <- list.files(paste0(path_to_tables, "miRNA_outputs/"), pattern = "_-_combined")
  pi_files <- list.files(paste0(path_to_tables, "piRNA_outputs/"), pattern = "_pi-zscore")

  if (ml_plots == TRUE) {
    radar_files <- list.files(paste0(path_to_tables, "radar_plots/"), pattern = "_prob.png")
  }

  pref <- "<a href="
  suff <- ' style="color:blue; text-decoration:none;">'


  for (i in 1:nrow(ml_tab)) {
    loc <- ml_tab$locus[i]
    mi_plus_idx <- grep(loc, mi_plus_files)
    if (!identical(mi_plus_idx, integer(0))) {
      mi_plus_file <- paste0(pref, "'", "miRNA_outputs/", mi_plus_files[mi_plus_idx], "'", suff, "miRNA plus</a>")
    } else {
      mi_plus_file <- NA
    }
    ml_tab$mi_plus_col[i] <- mi_plus_file

    mi_minus_idx <- grep(loc, mi_minus_files)
    if (!identical(mi_minus_idx, integer(0))) {
      mi_minus_file <- paste0(pref, "'", "miRNA_outputs/", mi_minus_files[mi_minus_idx], "'", suff, "miRNA minus</a>")
    } else {
      mi_minus_file <- NA
    }
    ml_tab$mi_minus_col[i] <- mi_minus_file

    pi_idx <- grep(loc, pi_files)
    if (!identical(pi_idx, integer(0))) {
      pi_file <- paste0(pref, "'", "piRNA_outputs/", pi_files[pi_idx], "'", suff, "piRNA</a>")
    } else {
      pi_file <- NA
    }
    ml_tab$pi_col[i] <- pi_file

    si_idx <- grep(loc, si_files)
    if (!identical(si_idx, integer(0))) {
      si_file <- paste0(pref, "'", "siRNA_outputs/", si_files[si_idx], "'", suff, "siRNA</a>")
    } else {
      si_file <- NA
    }
    ml_tab$si_col[i] <- si_file

    if (ml_plots == TRUE) {
      radar_idx <- grep(loc, radar_files)
      if (!identical(radar_idx, integer(0))) {
        radar_file <- paste0(pref, "'", "radar_plots/", radar_files[radar_idx], "'", suff, "probability<a>")
      } else {
        radar_file <- NA
      }
      ml_tab$prob_col[i] <- radar_file
    }
  }



  write.table(ml_tab, file = paste0(path_to_tables, "linked_", ml_file), sep = "\t", row.names = FALSE, col.names = TRUE, na = "NA")

  ########################################################## write and render the Rmd ##################################################################
  ml_cat_stat <- paste0('tab <- read.csv("', "linked_", ml_file, '", sep = "\t", header = TRUE)\n')

  # test <- rvest::read_html(paste0(path_to_tables, "linked_", ml_file))

  sink(paste0(path_to_tables, type, "_", "misipi_summary_page.Rmd"))
  cat(c("---\n"), append = TRUE)
  cat(c("title: 'MiSiPi.RNA Summary'\n"), append = TRUE)
  cat(c("output:\n"), append = TRUE)
  cat(c("   html_document:\n"), append = TRUE)
  cat(c("toc: no\n"), append = TRUE)
  cat(c("---\n"), append = TRUE)
  cat(c("<style type='text/css'>\n"), append = TRUE)
  cat(c("   .main-container {\n"), append = TRUE)
  cat(c("max-width: 1500px;\n"), append = TRUE)
  # cat(c("height: 1000px;\n"), append = TRUE)
  cat(c("margin-left: 100px;\n"), append = TRUE)
  cat(c("margin-right: 50px;\n"), append = TRUE)
  # cat(c("margin-bottom: 40px;\n"),append = TRUE)
  cat(c("}\n"), append = TRUE)
  cat(c(".dataTables_scrollBody {\n"), append = TRUE)
  cat(c("height: auto !important;\n"), append = TRUE)
  cat(c("}\n"), append = TRUE)
  cat(c("</style>\n"), append = TRUE)
  cat("\\newpage\n", append = TRUE)
  # cat("```{r, echo = FALSE, out.width=c('350px', '300px', '300px', '300px','300px'), out.height = c('400px','350px', '350px', '350px','350px'), fig.show='hold' }\n",append = TRUE)
  # cat('knitr::include_graphics("pi_overlap_heatmap.png", "pi_heatmap.png", "si_dicerz_heatmap.png", "plus_pi_phased_heat.png", "minus_pi_phased_heat.png"))\n', append = TRUE)
  cat(cat_sizes, append = TRUE)
  cat(cat_names, append = TRUE)
  cat("```\n", append = TRUE)

  cat("```{r, echo = FALSE, warning = FALSE, message = FALSE}\n", append = TRUE)

  cat(ml_cat_stat, append = TRUE)
  cat('dat <- DT::datatable(tab, options = list(), class = "display",
                   callback = DT::JS("return table;"),
                   caption = NULL, filter = c("top"), escape = FALSE,
                   style = "auto", width = NULL, elementId = NULL,
                   fillContainer = TRUE,
                   autoHideNavigation = getOption("DT.autoHideNavigation", NULL),
                   selection = c("multiple", "single", "none"), extensions = list(),
                   plugins = NULL, editable = FALSE)\n', append = TRUE)
  cat("dat\n", append = TRUE)
  cat("```\n", append = TRUE)

  sink()

  rmarkdown::render(paste0(path_to_tables, type, "_", "misipi_summary_page.Rmd"))

  file.remove(paste0(path_to_tables, "linked_", ml_file))
  rmds <- list.files(path_to_tables, pattern = ".Rmd")
  file.remove(paste0(path_to_tables, rmds))
}
