library(dplyr)

# read in output from misipi_rna() plus a column containing the class
# and a column containing the library size of the BAM file used to make the table
tbl <- read.table("C:/Users/tmjar/Desktop/MiSiPi_Stuff/machine_learning_species/new_all_species/table_versions/new_all_species_training_more_cisnats.txt", sep = "\t", header = TRUE)

# make sure rows are unique
tbl <- tbl %>% dplyr::distinct(log_shap_p,auc, strand_bias,perc_GC,
                               ave_size, perc_first_nucT, perc_A10, highest_si_col, si_dicerz, num_si_dicer_reads,
                               hp_perc_paired, hp_phasedz, hp_mfe, hp_dicerz, mi_perc_paired, mirna_dicerz, mirna_mfe,
                               mirna_overlapz, pingpong_col, max_pi_count, max_piz_overlap, pi_phasedz, pi_phased26z,
                               unique_read_bias, type, species, libsize, X, .keep_all = TRUE)

# Deal with missing max_piz_overlaps
idx <- which(is.na(tbl$max_piz_overlap))
tbl$max_piz_overlap[idx] <- -33

idx <- which(is.na(tbl$hp_mfe))
idx2 <- which(tbl$hp_mfe == 0)

#set missing hp MFE's to unlikely number
tbl$hp_mfe[idx] <- 10000
tbl$hp_mfe[idx2] <- 10000

idx <- which(tbl$mirna_mfe == 0)
tbl$mirna_mfe[idx] <- 1000


for(i in 1:nrow(tbl)){
  if(tbl$max_pi_count[i] == -33){
    tbl$max_pi_count[i] <- 0
  }
}

# Group all types of contamination into one class
for(i in 1:nrow(tbl)){
  if(tbl$type[i] == "snoRNA" || tbl$type[i] == "tRNA" || tbl$type[i] == "rRNA" || tbl$type[i] == "sno"){
    tbl$type[i] <- "contamination"
  }
}

#normalize counts by libsize

tbl <- tbl %>% dplyr::mutate(pi_norm = (max_pi_count/libsize/1000000), si_norm = (num_si_dicer_reads/libsize/1000000))


all_tbl <- tbl %>% dplyr::select(-c(locus_length, species, X, max_pi_count, num_si_dicer_reads, libsize))
all_tbl$pi_norm <- as.numeric(all_tbl$pi_norm)
all_tbl$si_norm <- as.numeric(all_tbl$si_norm)

#write.table(all_tbl, file = "new_ML_contam_more_libsize_norm_mod.txt", sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
# install required packages if running for the first time
library(dplyr)   # for data cleaning
library(CatEncoders)
library(xgboost)
library(e1071)

# The MFE doesn't need to be negative, so take absolute value
all_tbl$mirna_mfe <- abs(all_tbl$mirna_mfe)

all_tbl$hp_mfe <- abs(all_tbl$hp_mfe)


# remove the locus name and
# Hot encode the type class
tbl <- all_tbl %>% dplyr::select(-c(locus))
train_labs <- LabelEncoder.fit(tbl$type)
tbl$type <- transform(train_labs, tbl$type)
# tranform to be zero-based
tbl$type <- tbl$type - 1


# check whether there is a class imbalance
table(tbl$type)                     # Raw class counts
prop.table(table(tbl$type)) * 100  #

# If serious class imbalance exists, consider methods to either re-sample,
# add more data, or use class weighting.
# tbl$type <- as.factor(tbl$type)  # caret requires a factor

# or can use SMOTE:
# Better than simple upsampling — creates synthetic examples of minority classes.
# remotes::install_github("cran/DMwR")
#  library(DMwR)
# #
# tbl$type <- as.factor(tbl$type)
# balanced_data <- SMOTE(type ~ ., data = tbl, perc.over = 100, perc.under = 200)
# tbl <- balanced_data

# Create Training and Test data -
set.seed(100)  # setting seed to reproduce results of random sampling
# Segregate 80% of data for training and leave 20% out for cross validation
trainingRowIndex <- sample(1:nrow(tbl), 0.8*nrow(tbl))# row indices for training data
trainingData <- tbl[trainingRowIndex, ] #%>% dplyr::select(-c(locus, species))  # model training data
testData  <- tbl[-trainingRowIndex, ]  #%>% dplyr::select(-c(locus, species))  # test data


train_x <- model.matrix(type ~ ., trainingData)[, -1]

test_x <- model.matrix(type ~ ., testData)[, -1]


train_y <- trainingData$type


test_y <- testData$type


train_y <- as.integer(as.character(train_y))  # safely convert factor to numeric
test_y  <- as.integer(as.character(test_y))
# convert the train and test data into xgboost matrix type.
dtrain = xgb.DMatrix(data=train_x, label=train_y)
dval = xgb.DMatrix(data=test_x, label=test_y)

# Store X and Y for later use
#y = train_y
#x = train_x

options(scipen = 999)

#featurePlot(x = x,
#            y = as.factor(y),
#            plot = "density",
#            strip = strip.custom(par.strip.text = list(cex = .7)),
#            scales = list(x = list(relation="free"),
#                          y = list(relation="free")))

## read in the cross-validation table and transform the values same as the training data
# I made my own cross-validation with completely unseen data to get a better idea of real life accuracy,
# as the held-back training data seemed to have inflated accuracies.
# This table is formatted the same as the ml table.

dat <- read.table("C:/Users/tmjar/Desktop/MiSiPi_Stuff/machine_learning_species/new_all_species/table_versions/all_species_w_whitefly_CV.txt", header = TRUE)

# Group all types of contamination as one class
dat <- dat %>% dplyr::select(-c(locus, locus_length, X,species))
for(i in 1:nrow(dat)){
  if(dat$type[i] == "snoRNA" || dat$type[i] == "tRNA" || dat$type[i] == "rRNA" || dat$type[i] == "sno" || dat$type[i] == "contamination"){
    dat$type[i] <- "contamination"
  }
}


dat <- dat %>% dplyr::mutate(pi_norm = "", si_norm = "")

#normalize counts by libsize

for(i in 1:nrow(dat)){
  #idx <- which(df$dataset == dat$X[i])
  dat$pi_norm[i] <- dat$max_pi_count[i]/dat$libsize[i]/1000000
  dat$si_norm[i] <- dat$num_si_dicer_reads[i]/dat$libsize[i]/1000000

}
dat$pi_norm <- as.numeric(dat$pi_norm)
dat$si_norm <- as.numeric(dat$si_norm)
dat <- dat %>% dplyr::select(-c(libsize, max_pi_count, num_si_dicer_reads))


test_labs <- LabelEncoder.fit(dat$type)
dat$type <- transform(test_labs, dat$type)
dat$type <- dat$type - 1

idx <- which(is.na(dat$max_piz_overlap))
dat$max_piz_overlap[idx] <- -33

# Set missing MFE to unlikely number

idx1 <- which(is.na(dat$hp_mfe))
idx2 <- which(dat$hp_mfe == 0)
dat$hp_mfe[idx1] <- 10000
dat$hp_mfe[idx2] <- 10000

idx1 <- which(is.na(dat$mirna_mfe))
idx2 <- which(dat$mirna_mfe == 0)


dat$mirna_mfe[idx1] <- 1000
dat$mirna_mfe[idx2] <- 1000

dat$hp_mfe <- abs(dat$hp_mfe)


new_data_y <- dat$type
#new_data_y <- new_data_y - 1
new_data_x <- dat %>% dplyr::select(-type)



####### Do the training

rounds <- c(700,800,900,1000)
etas1 <- c(seq(0.01, 0.3, 0.05))
#etas <- c(etas1, etas2)
#etas <- c(seq(0.01, 0.4, 0.1))
subs <- c(0.7,0.75,0.8, 0.85,0.9, 0.95, 0.1)
depths <- c(5,6,7,8,9,10,11)
lambdas <- c(0.01,0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09,0.1,1,0)
alphas <- c(0,0.1,0.5,1.0)
gammas <- c(0.05,0.06,0.07,0.08,0.09,0.1,0.3,0.5,0.7,0.9,1.0)

len <- length(gammas)
n_iter <- len

pb <- txtProgressBar(min = 0,
                     max = n_iter,
                     style = 3,
                     width = 50,
                     char = "=")

results <- data.frame(par = numeric(length = length(len)), mean_accuracy = numeric(length = length(len)))





for( i in 1:length(gammas)){
  set.seed(1234)
  #run with softmax for quick accuracy check, pick best model, re-run with softprob and save
  cur_model <- xgboost(data = dtrain, max.depth = 6, subsample = 0.85, eta = 0.01, nthread = 2, lambda = 0.01, alpha = 0.5,
                       gamma = 0.05, eval_metric = "mlogloss", nrounds = 700, objective = "multi:softmax", num_class = 5, verbose = 0,
                       early_stopping_rounds = 10, print_every_n = 10)
  #prob_model <- xgboost(data = dtrain, max.depth = 6, eval_metric = "mlogloss", eta = 0.01, subsample = 0.85, nthread = 2,
  #                      alpha = 0.1, lambda = 0.01, gamma = 0.05,
  #                      nrounds = 700, objective = "multi:softprob", num_class = 5, verbose = FALSE)
  real_pred <- predict(cur_model, as.matrix(new_data_x))
  #df <- as.data.frame(matrix(real_pred, ncol = 2, byrow = TRUE)) #%>% dplyr::select(c(V1))

  prob_mat <- matrix(real_pred, ncol = 5, nrow = nrow(new_data_x))
  cur_err <- mean(real_pred != new_data_y)
  cur_accuracy <- 1 - cur_err
  print(cur_accuracy)

  results[i,] <- c(gammas[i], cur_accuracy)
  setTxtProgressBar(pb, i)
}



# To check whether model is overfit
# log <- model$evaluation_log
# plot(log$iter, log$train_mlogloss, type = "l", col = "blue")
# plot(log$iter, log$eval_mlogloss, type = "l", col = "red")
# legend("topright", legend = c("Train", "Validation"), col = c("blue", "red"), lty = 1)
#

# save the model
xgb.save(prob_model, "C:/Users/tmjar/Desktop/MiSiPi.RNA/inst/extdata/xgb_tuned.rds")

# load the model
all_model <- xgb.load("C:/Users/tmjar/Desktop/MiSiPi_Stuff/machine_learning_species/new_all_species/xgb_balanced.rds")

# Calculate Type I and Type II error rates
all_pred <- predict(all_model, as.matrix(new_data_x))
all_df <- as.data.frame(matrix(all_pred, ncol = 5, byrow = TRUE))
#all_pred <- all_pred[,1]
colnames(all_df) <- c("prob_cis", "prob_contam", "prob_hp", "prob_mi", "prob_pi")
all_df$real_type <- dat$type

for(i in 1:nrow(all_df)){
  real <- all_df$real_type[i]
  #true and contam
  #true and not contam
  #false and contam
  #false and not contam
  max <- which(all_df[i,1:5] == max(all_df[i,1:5]))
  if((max - 1) == real){
    all_df$cond[i] <- "TRUE"
  } else {
    all_df$cond[i] <- "FALSE"
  }

  if(all_df$real_type[i] == 1){
    all_df$contam[i] <- "contam"
  } else {
    all_df$contam[i] <- "not_contam"
  }
}


