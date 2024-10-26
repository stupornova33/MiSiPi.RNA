#' Make_html_summary
#' Takes a file path and creates a graphic summary of all loci from a MiSiPi run
#' @param path_to_tables A string specifying the full path to the folder containing the table outputs from MiSiPi.RNA
#' @param type A string specifying which module of MiSiPi.RNA should be plotted in the summary. Options are one of: "miRNA", "siRNA", or "piRNA"
#' @return nothing
#' @export


make_html_summary <- function(path_to_tables, type){

  # Expects: a read size dist table
  # output ML table from misipi_rna
  #
  if(type == "siRNA" || type == "sirna"){
    tables_needed <- c(siRNA_dicerz, plus_hp_phasedz, minus_hp_phasedz)
  } else if(type == "miRNA" || type == "mirna"){
    tables_needed <- c(miRNA_dicerz)
  } else {
    tables_needed <- c(all_phased_zscores, pi_heatmap, minus_phasedz, plus_phasedz)
  }


  ml_tab <- list.files(path_to_tables, pattern = "_ML.txt")

  size_dist_tab <- list.files(path_to_tables, pattern = "_size_distributions.txt")
  tab <- read.csv(paste0(path_to_tables, "piRNA_read_size_distributions.txt"), sep = "\t", header = FALSE)

  RowVar <- function(x, ...) {
    rowSums((x - rowMeans(x, ...))^2, ...)/(dim(x)[2] - 1)
  }

  tab_mat <- tab[,2:ncol(tab)]
  var <- RowVar(tab_mat)
  idx <- which(var == 0)

  if(length(idx) > 0){
    tab_mat <- tab_mat[-c(idx),]
  }

  names <- c("16","17","18","19","20","21","22","23","24","25","26","27","28","29","30", "31", "32")
  colnames(tab_mat) <- names
  heatmap <- pheatmap::pheatmap(tab_mat, main = "Read size distribution", cluster_cols = FALSE, show_rownames = F,  show_colnames = T, scale = "row")

  png(filename="pi_heatmap.png")
  heatmap
  dev.off()

  # read in other plots and tables
  si_dicerz <- read.table(paste0(path_to_tables, "siRNA_dicerz.txt"), header = TRUE, sep = "\t") %>% dplyr::select(-c(locus))

  colnames(si_dicerz) <- c("-4","-3","-2","-1","0","1","2","3","4")

  var <- RowVar(si_dicerz)
  idx <- which(var == 0)

  if(length(idx) > 0){
    si_dicerz <- si_dicerz[-c(idx),]
  }
  #if row clustering is to be done, need to remove rows with no variance
  dicer_heat <- pheatmap::pheatmap(si_dicerz, main = "2nt Overhang Z-score", fontsize_col = 15,cluster_cols = FALSE, show_rownames = FALSE, cluster_rows = FALSE, scale = "row")

  png(filename="si_dicerz_heatmap.png")
  dicer_heat
  dev.off()


  ####### piRNA phasing

  m_pi_phased_tab <- read.table(paste0(path_to_tables, "phased_minus_piRNA_zscores.txt"), header = TRUE, sep = '\t') %>%
    dplyr::select(V2,V3,V4,V5,V6)
  m_pi_phased_tab[is.na(m_pi_phased_tab)] <- -33
  colnames(m_pi_phased_tab) <- c("1nt", "2nt", "3nt","4nt", "5nt")

  var <- RowVar(m_pi_phased_tab)
  idx <- which(var == 0)

  if(length(idx) > 0){
    m_pi_phased_tab <- m_pi_phased_tab[-c(idx),]
  }

  m_pi_phased_heat <- pheatmap::pheatmap(m_pi_phased_tab, main = "Phased piRNA Z-scores (minus strand)", fontsize_col = 15, cluster_cols = FALSE, show_rownames = FALSE, cluster_rows = TRUE, scale = "row")

  png(filename="minus_pi_phased_heat.png")
  m_pi_phased_heat
  dev.off()



  p_pi_phased_tab <- read.table(paste0(path_to_tables, "phased_plus_piRNA_zscores.txt"), header = TRUE, sep = '\t') %>%
    dplyr::select(V2,V3,V4,V5,V6)
  p_pi_phased_tab[is.na(p_pi_phased_tab)] <- -33
  colnames(p_pi_phased_tab) <- c("1nt", "2nt", "3nt","4nt", "5nt")

  var <- RowVar(p_pi_phased_tab)
  idx <- which(var == 0)

  if(length(idx) > 0){
    p_pi_phased_tab <- p_pi_phased_tab[-c(idx),]
  }

  p_pi_phased_heat <- pheatmap::pheatmap(p_pi_phased_tab, main = "Phased piRNA Z-scores (plus strand)",fontsize_col = 15,  cluster_cols = FALSE, show_rownames = FALSE, cluster_rows = TRUE, scale = "row")

  png(filename="plus_pi_phased_heat.png")
  m_pi_phased_heat
  dev.off()


  ########################################################## write and render the Rmd #########################################################


  sink(paste0(path_to_tables, "misipi_summary_page.Rmd"))
  cat(c("---\n"), append = TRUE)
  cat(c("title: 'MiSiPi.RNA Summary'\n"), append = TRUE)
  cat(c("output:\n"), append = TRUE)
  cat(c("   html_document:\n"), append = TRUE)
  cat(c("toc: no\n"), append = TRUE)
  cat(c("---\n"), append = TRUE)
  cat(c("<style type='text/css'>\n"),append = TRUE)
  cat(c("   .main-container {\n"),append = TRUE)
  cat(c("max-width: 1800px;\n"),append = TRUE)
  cat(c("margin-left: auto;\n"),append = TRUE)
  cat(c("margin-right: auto;\n"),append = TRUE)
  cat(c("}\n"),append = TRUE)
  cat(c("</style>\n"),append = TRUE)
  cat("\\newpage\n",append = TRUE)
  cat("```{r, echo = FALSE, out.width=c('350px', '300px', '300px', '300px','300px'), out.height = c('400px','350px', '350px', '350px','350px'), fig.show='hold' }\n",append = TRUE)
  cat('knitr::include_graphics("pi_overlap_heatmap.png", "pi_heatmap.png", "si_dicerz_heatmap.png", "plus_pi_phased_heat.png", "minus_pi_phased_heat.png"))\n', append = TRUE)

  cat("```\n",append = TRUE)

  cat('```{r, echo = FALSE, warning = FALSE, message = FALSE}\n', append = TRUE)

  cat('tab <- read.table("all_dsim_ML.txt", header = TRUE)\n')
  cat('dat <- DT::datatable(tab, options = list(), class = "display",
                   callback = DT::JS("return table;"),
                   caption = NULL, filter = c("none", "bottom", "top"), escape = TRUE,
                   style = "auto", width = NULL, height = NULL, elementId = NULL,
                   fillContainer = getOption("DT.fillContainer", NULL),
                   autoHideNavigation = getOption("DT.autoHideNavigation", NULL),
                   selection = c("multiple", "single", "none"), extensions = list(),
                   plugins = NULL, editable = FALSE)\n', append = TRUE)
  cat('dat\n', append = TRUE)
  cat('```\n', append = TRUE)

  sink()



  rmarkdown::render(paste0(path_to_tables, "misipi_summary_page.Rmd"))



}

