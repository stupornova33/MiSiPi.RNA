#' ml_probability
#' @param path_to_table A string specifying the full path to the folder which contains the _ml table. This directory will be the location of the probability plots.
#' @param table A string specifying the full path to the _ml table produced by misipi_rna
#' @return nothing
#' @export

ml_probability <- function(path_to_table, table){
  wkdir <- setwd(path_to_table)
  
  ml_table <- read.table(paste0(path_to_table, table), header = TRUE)
  names <- ml_table$locus
  ml_table <- ml_table %>% dplyr::select(-c(locus, locus_length))
  
  ml_table$mirna_mfe <- abs(ml_table$mirna_mfe)
  
  ml_table$hp_mfe <- abs(ml_table$hp_mfe)
  
  
  model_path <- system.file("extdata", "xgb_all_classes.rds", package="MiSiPi.RNA")
  all_model <- xgboost::xgb.load(model_path)
  
  # columns are always in alphabetical order
  #hpRNA not_hpRNA
  
  options(scipen = 999)

  
  all_pred <- predict(all_model, as.matrix(ml_table))
  all_df <- as.data.frame(matrix(all_pred, ncol = 5, byrow = TRUE))
  #colnames are c("prob_cis", "prob_contam", "prob_hp", "prob_mi", "prob_pi")
  colnames(all_df) <- c("Prob. cisNAT", "Prob. contam", "Prob. hpRNA", "Prob. miRNA", "Prob. piRNA")
  
  all_df <- all_df*100
  
  max <- c(100,100, 100, 100,100)
  min <- c(0, 0,0, 0,0)
  
  
  if(!dir.exists(paste0(wkdir, "/radar_plots/"))){
    dir.create(paste0(wkdir, "/radar_plots/"))
    
  }
  
  radar_dir <- paste0(wkdir, "/radar_plots/")
  for(i in 1:nrow(all_df)){
    new_df <- rbind(max, min, all_df[i,])
    
    rownames(new_df) <- c("Max.", "Min.", "Values")
    
    png(paste0(radar_dir, names[i], "_prob.png"), width = 850, height = 800, units = "px")
    fmsb::radarchart(new_df[1:3,], axistype =2, seg = 5, pcol=rgb(0.2,0.5,0.5,0.9) , pfcol=rgb(0.2,0.5,0.5,0.5) , plwd=4,
                     cglcol="darkgrey", cglty=1, axislabcol="darkgrey", caxislabels=seq(0,100,10), cglwd=0.8,
                     #custom labels
                     vlcex=2, palcex = 2)
    dev.off()
    
    
    
    
  }
  
  print("ML probability plots have been made.")
}
