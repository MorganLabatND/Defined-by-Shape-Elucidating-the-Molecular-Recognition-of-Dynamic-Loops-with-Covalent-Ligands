setwd('/lsi/home/cemarkey/RerunningML/ConMod')

# INSTALL AND LOAD ALL NECESSARY PACKAGES AND LIBRARIES -----------------------------------------------------------------
install.packages(c("stringi", 'glmnet', 'ModelMetrics', 'pls', 'randomForest', 'kernlab', 'earth', 'prospectr'), repos="http://cran.r-project.org", lib="/lsi/home/cemarkey/RerunningML/Anaconda_libs/")
install.packages('caret', dependencies = c("Depends", "Imports", "Suggests"), repos="http://cran.r-project.org", lib="/lsi/home/cemarkey/RerunningML/Anaconda_libs/")
install.packages('ggplot2', dependencies = TRUE, repos="http://cran.r-project.org", lib="/lsi/home/cemarkey/RerunningML/Anaconda_libs/")
install.packages('withr', dependencies = TRUE, repos="http://cran.r-project.org", lib="/lsi/home/cemarkey/RerunningML/Anaconda_libs/")
install.packages('Formula', dependencies = TRUE, repos="http://cran.r-project.org", lib="/lsi/home/cemarkey/RerunningML/Anaconda_libs/")
install.packages('plotmo', dependencies = TRUE, repos="http://cran.r-project.org", lib="/lsi/home/cemarkey/RerunningML/Anaconda_libs/")
install.packages('plotrix', dependencies = TRUE, repos="http://cran.r-project.org", lib="/lsi/home/cemarkey/RerunningML/Anaconda_libs/")
install.packages('TeachingDemos', dependencies = TRUE, repos="http://cran.r-project.org", lib="/lsi/home/cemarkey/RerunningML/Anaconda_libs/")
install.packages('recipes', dependencies = TRUE, repos="http://cran.r-project.org", lib="/lsi/home/cemarkey/RerunningML/Anaconda_libs/")
install.packages('dplyr', dependencies = TRUE, repos="http://cran.r-project.org", lib="/lsi/home/cemarkey/RerunningML/Anaconda_libs/")
install.packages('crayon', dependencies = TRUE, repos="http://cran.r-project.org", lib="/lsi/home/cemarkey/RerunningML/Anaconda_libs/")
install.packages('Matrix', dependencies = TRUE, repos="http://cran.r-project.org", lib="/lsi/home/cemarkey/RerunningML/Anaconda_libs/")
install.packages('lattice', dependencies = TRUE, repos="http://cran.r-project.org", lib="/lsi/home/cemarkey/RerunningML/Anaconda_libs/")
library(lattice, lib="/lsi/home/cemarkey/RerunningML/Anaconda_libs/")
library(Matrix, lib="/lsi/home/cemarkey/RerunningML/Anaconda_libs/")
library(crayon, lib="/lsi/home/cemarkey/RerunningML/Anaconda_libs/")
library(dplyr, lib="/lsi/home/cemarkey/RerunningML/Anaconda_libs/")
library(recipes, lib="/lsi/home/cemarkey/RerunningML/Anaconda_libs/")
library(Formula, lib="/lsi/home/cemarkey/RerunningML/Anaconda_libs/")
library(TeachingDemos, lib="/lsi/home/cemarkey/RerunningML/Anaconda_libs/")
library(plotrix, lib="/lsi/home/cemarkey/RerunningML/Anaconda_libs/")
library(plotmo, lib="/lsi/home/cemarkey/RerunningML/Anaconda_libs/")
library(withr, lib="/lsi/home/cemarkey/RerunningML/Anaconda_libs/")
library(ggplot2, lib="/lsi/home/cemarkey/RerunningML/Anaconda_libs/")
library(stringi, lib="/lsi/home/cemarkey/RerunningML/Anaconda_libs/")
library(caret, lib="/lsi/home/cemarkey/RerunningML/Anaconda_libs/")
library(glmnet, lib="/lsi/home/cemarkey/RerunningML/Anaconda_libs/")
library(ModelMetrics, lib="/lsi/home/cemarkey/RerunningML/Anaconda_libs/")
library(pls, lib="/lsi/home/cemarkey/RerunningML/Anaconda_libs/")
library(randomForest, lib="/lsi/home/cemarkey/RerunningML/Anaconda_libs/")
library(kernlab, lib="/lsi/home/cemarkey/RerunningML/Anaconda_libs/")
library(earth, lib="/lsi/home/cemarkey/RerunningML/Anaconda_libs/")
library(prospectr, lib="/lsi/home/cemarkey/RerunningML/Anaconda_libs/")

# Read in each bent training set --------------------------------------------------------------
myTrainFiles <- list.files(pattern="trainSet.*csv")
trainDat<-lapply(myTrainFiles,read.csv)
fileNumbers <- regmatches(myTrainFiles, regexpr( "\\d+", myTrainFiles))

# Loop through each training set and implement external scrambling workflow --------------------
for (i in 1:length(myTrainFiles)){
  
  # Load current training set
  trainingSet <- trainDat[[i]]
  trainingSet <- subset(trainingSet, select = -c(Group))
  colnames(trainingSet)[3] <- "y"
  
  ## IDENTIFY AND REMOVE INVARIANT DESCRIPTORS
  set.seed(1234)
  
  # Create a new training data set with invariant descriptors removed
  all_desc <- colnames(trainingSet)
  invar_desc_removed_trainSet <- trainingSet[vapply(trainingSet, function(z) length(unique(z))>1, logical(1L))]
  
  # Find column indices corresponding to invariant descriptors
  reduced_desc <- colnames(invar_desc_removed_trainSet)
  invariant_desc <- setdiff(all_desc, reduced_desc)
  invariant_inds <- integer()
  
  for (a in 1:length(invariant_desc)){
    index <- which(all_desc == invariant_desc[a])
    invariant_inds <- c(invariant_inds, index)
  }
  
  # Save remaining descriptors as a model matrix
  x_train = model.matrix(y ~.-ID-SMILES, data = invar_desc_removed_trainSet)
  x_train = x_train[,-1]
  
  # Save molecule identifier information and Tethering (y) values
  ID = invar_desc_removed_trainSet$ID
  SMILES = invar_desc_removed_trainSet$SMILES
  y = as.numeric(invar_desc_removed_trainSet$y)
  
  # IDENTIFY AND REMOVE PERFECTLY CORRELATED DESCRIPTORS
  # Generate matrix containing correlations between each descriptors and Tethering
  bent_desc_tethering_cormat = cor(x_train, y, method = "kendall")
  
  # Generate matrix containing correlations between all descriptors and each other
  bent_descriptor_cormat = cor(x_train, method = "kendall")
  bent_descriptor_cormat = round(bent_descriptor_cormat, 3)
  
  # Set diagonal and upper triangle of matrix to 0 to prevent repeat analysis
  bent_descriptor_cormat[upper.tri(bent_descriptor_cormat)] <- 0
  diag(bent_descriptor_cormat) <- 0

  # Initialize empty list to store indices of descriptors to be removed
  ind_to_remove <- list()
  d = NULL
  
  # Iterate through matrix to identify perfectly correlated descriptors pairs, removing descriptor 
  # of pair with lower correlation to Tethering
  for (k in 1:ncol(bent_descriptor_cormat)){
    for (j in 1:nrow(bent_descriptor_cormat)){
      if (abs(bent_descriptor_cormat[j,k]) == 1.000){
        if (bent_desc_tethering_cormat[k] > bent_desc_tethering_cormat[j]){
          ind_to_remove <- c(ind_to_remove, j)
        }
        if (bent_desc_tethering_cormat[j] > bent_desc_tethering_cormat[k]){
          ind_to_remove <- c(ind_to_remove, k)
        }
        if (bent_desc_tethering_cormat[j] == bent_desc_tethering_cormat[k]){
          ind_to_remove <- c(ind_to_remove,j)
        }
        print(paste(colnames(bent_descriptor_cormat)[k], "-" , rownames(bent_descriptor_cormat)[j], ": ", bent_descriptor_cormat[j,k]))
        d = rbind(d, data.frame(colnames(bent_descriptor_cormat)[k], bent_desc_tethering_cormat[k], rownames(bent_descriptor_cormat)[j], bent_desc_tethering_cormat[j]))
      }
    }
  }
  
  # Generate ordered list of descriptor indicies removed via analysis
  new_ind_to_remove <- unique(ind_to_remove)
  new_ind_to_remove <- unlist(new_ind_to_remove)
  new_ind_to_remove <- sort(new_ind_to_remove)
  
  ## IDENTIFY AND REMOVE DESCRIPTORS CORRELATED > |0.8|
  # Update correlation matrices with new descriptors removed
  x_train = x_train[,-c(new_ind_to_remove)]
  bent_desc_tethering_cormat = cor(x_train, y, method = "kendall")
  bent_descriptor_cormat = cor(x_train, method = "kendall")
  bent_descriptor_cormat = round(bent_descriptor_cormat, 3)
  bent_descriptor_cormat[upper.tri(bent_descriptor_cormat)] <- 0
  diag(bent_descriptor_cormat) <- 0

  # Generate a vector of descriptors to remove to reduce pair-wise correlations
  built_in_inds_to_remove_80 <- findCorrelation(bent_descriptor_cormat, cutoff = 0.800)
  
  # Initialize empty list to store indices of descriptors to be removed
  ind_to_remove_80 <- list()
  
  # Iterate through matrix to identify correlated descriptors pairs, removing descriptor 
  # of pair with lower correlation to Tethering
  for (K in 1:ncol(bent_descriptor_cormat)){
    for (j in 1:nrow(bent_descriptor_cormat)){
      if (abs(bent_descriptor_cormat[j,K]) >= 0.800 & !(j %in% ind_to_remove_80) & !(K %in% ind_to_remove_80)){
        if (abs(bent_desc_tethering_cormat[K]) > abs(bent_desc_tethering_cormat[j]) & !(j %in% ind_to_remove_80)){
          ind_to_remove_80 <- c(ind_to_remove_80, j)
        }
        if (abs(bent_desc_tethering_cormat[j]) > abs(bent_desc_tethering_cormat[K]) & !(K %in% ind_to_remove_80)){
          ind_to_remove_80 <- c(ind_to_remove_80, K)
        }
        if (abs(bent_desc_tethering_cormat[j]) == abs(bent_desc_tethering_cormat[K]) & !(K %in% ind_to_remove_80) & !(j %in% ind_to_remove_80)){
          if (j %in% built_in_inds_to_remove_80 & !(K %in% built_in_inds_to_remove_80)){
            ind_to_remove_80 <- c(ind_to_remove_80, j)
          }
          if (K %in% built_in_inds_to_remove_80 & !(j %in% built_in_inds_to_remove_80)){
            ind_to_remove_80 <- c(ind_to_remove_80, K)
          }
          else{
            ind_to_remove_80 <- c(ind_to_remove_80, K, j)
          }
        }
      }
    }
  }
  
  # Generate ordered list of descriptor indicies removed via analysis
  new_ind_to_remove_80 <- unique(ind_to_remove_80)
  new_ind_to_remove_80 <- unlist(new_ind_to_remove_80)
  new_ind_to_remove_80 <- sort(new_ind_to_remove_80)
  
  # Save list of descriptors remaining post analysis and generate updated data set
  remaining_descriptors_80 <- x_train[,-c(new_ind_to_remove_80)]
  updated_data <- trainingSet[,c("ID", "SMILES", "y", colnames(remaining_descriptors_80))]
  
  ## SCALE DESCRIPTORS
  # Generate model matrix of training set containing remaining descriptors
  x_new = model.matrix(y ~ .-ID-SMILES, data = updated_data)
  
  # Scale all descriptor values in training set model matrix
  x_means <- colMeans(x_new)
  x_sd <- apply(x_new, 2, sd)
  for (K in 1:(ncol(x_new))){
    for (j in 1:(nrow(x_new))){
      scaled_val <- (x_new[j,K]-x_means[K])/(x_sd[K])
      x_new[j,K] <- scaled_val
    }
  }
  x_new[,1] = 1
  
  ## LASSO FEATURE SELECTION
  # Fit a cross validated lasso model between scaled descriptors and Tethering
  fit = glmnet(x_new, y, alpha = 1, standardize = FALSE)
  cvfit = cv.glmnet(x_new, y, alpha = 1, keep = TRUE, standardize = FALSE, nfolds = length(y))

  # Identify descriptors associated with the lasso model with the minimum lambda value
  beta <- coef(cvfit, s = "lambda.min")
  ii <- which(beta!=0)
  lambda_min_descriptors <- beta[ii, 1]
  lambda_min_chosen_desc <- names(lambda_min_descriptors)
  lambda_min_chosen_desc <- lambda_min_chosen_desc[-1]
  
  # Generate an updated training set containing only descriptors selected by lasso lambda min model
  if (length(lambda_min_chosen_desc) >= 1){
    subset(x_new, select=c(lambda_min_chosen_desc))
    lassoMin_updatedData <- cbind.data.frame(ID, SMILES, y, subset(x_new, select=c(lambda_min_chosen_desc)))
  }

  # Identify descriptors associated with the lasso model with a lambda value 1 standard error from the minimum
  beta <- coef(cvfit, s = "lambda.1se")
  ii <- which(beta!=0)
  lambda_1se_descriptors <- beta[ii, 1]
  lambda_1se_chosen_desc <- names(lambda_1se_descriptors)
  lambda_1se_chosen_desc <- lambda_1se_chosen_desc[-1]
  
  # Generate an updated training set containing only descriptors selected by lasso lambda 1se model
  if (length(lambda_1se_chosen_desc) >= 1){
    subset(x_new, select=c(lambda_1se_chosen_desc))
    lasso1se_updatedData <- cbind.data.frame(ID, SMILES, y, subset(x_new, select=c(lambda_1se_chosen_desc)))
  }
  
  # Generate a training set containing all descriptors before feature selection applied
  noFS_updated_data <- cbind.data.frame(ID, SMILES, y, x_new[,-1])
  
  ## MARS FEATURE SELECTION
  # Fit a MARS model between scaled descriptors and Tethering
  earth.mod <- earth(x_new, y)
  
  # Function to identify names of descriptors used in MARS model
  get.used.pred.names <- function(obj) # obj is an earth object
  {
    any1 <- function(x) any(x != 0)    # like any but no warning if x is double
    names(which(apply(obj$dirs[obj$selected.terms, , drop=FALSE], 2, any1)))
  }
  
  # Identify descriptors used in MARS model
  mars_chosen_desc <- get.used.pred.names(earth.mod)

  # Generate a training set containing only descriptors selected by mars model
  subset(x_new, select=c(mars_chosen_desc))
  mars_updatedData <- cbind.data.frame(ID, SMILES, y, subset(x_new, select=c(mars_chosen_desc)))
  
  
  ## INITIALIZE EMPTY DATAFRAMES TO HOLD SCRAMBLED MACHINE LEARNING MODEL STATS AND LABEL COLUMNS PROPERLY
  lm_lr_stats <- cbind.data.frame("NA", "NA", "NA")
  l1_lr_stats <- cbind.data.frame("NA", "NA", "NA")
  mars_lr_stats <- cbind.data.frame("NA", "NA", "NA")
  lm_pls_stats <- cbind.data.frame("NA", "NA", "NA")
  l1_pls_stats <- cbind.data.frame("NA", "NA", "NA")
  mars_pls_stats <- cbind.data.frame("NA", "NA", "NA")
  l1_rf_stats <- cbind.data.frame("NA", "NA", "NA")
  lm_svm_stats <- cbind.data.frame("NA", "NA", "NA") 
  l1_svm_stats <- cbind.data.frame("NA", "NA", "NA")
  mars_svm_stats <- cbind.data.frame("NA", "NA", "NA")
  
  colnames(lm_lr_stats) <- rbind("R2", "MAE ", "RMSE")
  colnames(l1_lr_stats) <- rbind("R2", "MAE ", "RMSE")
  colnames(mars_lr_stats) <- rbind("R2", "MAE ", "RMSE")
  colnames(lm_pls_stats) <- rbind("R2", "MAE ", "RMSE")
  colnames(l1_pls_stats) <- rbind("R2", "MAE ", "RMSE")
  colnames(mars_pls_stats) <- rbind("R2", "MAE ", "RMSE")
  colnames(l1_rf_stats) <- rbind("R2", "MAE ", "RMSE")
  colnames(lm_svm_stats) <- rbind("R2", "MAE ", "RMSE")
  colnames(l1_svm_stats) <- rbind("R2", "MAE ", "RMSE")
  colnames(mars_svm_stats) <- rbind("R2", "MAE ", "RMSE")
  
  # Generate and save 45 lists containing the training set y (Tethering) values in randomly scrambled orders
  list_of_scrambs <- list()
  for (J in 1:45){
    scrambled_y <- sample(y, size = length(y), replace = FALSE)
    list_of_scrambs[[J]] <- c(scrambled_y)
  }
  
  ## LASSO MIN LINEAR REGRESSION
  # Check if lasso min feature selection returned any descriptors before proceeding
  if (length(lambda_min_chosen_desc) >= 1){
    # Initialize empty lists to store fit and error statistics for scrambled training set predictions made by model
    lm_lr_int_list_of_r2 <- list()
    lm_lr_int_list_of_rmse <- list()
    lm_lr_int_list_of_mae <- list()
    
    # Train a linear model using each of the scrambled training sets and calculate evaluation statistics
    for (Z in 1:length(list_of_scrambs)){
      # Generate a dataset containing current scrambled y values and properly ordered descriptors
      int_scramb_lr_lassoMin_updatedData <- cbind.data.frame(list_of_scrambs[[Z]], lassoMin_updatedData[,-3])
      colnames(int_scramb_lr_lassoMin_updatedData)[1] <- "y"
      # Specify cross validation method as leave-one-out cross validation
      train.control <- trainControl(method = "LOOCV")
      # Fit a cross validated linear model between descriptors and Tethering
      int_scramb_LassoMin_linearModel <- train(y~.-ID-SMILES, data = int_scramb_lr_lassoMin_updatedData, method = "lm", trControl = train.control)
      int_scramb_LassoMin_linearFittedVals <- int_scramb_LassoMin_linearModel$finalModel$fitted.values
      # Store fit and error statistics for training set predictions made by model
      int_scramb_LassoMin_linearR2 <- summary(int_scramb_LassoMin_linearModel$finalModel)$r.squared
      lm_lr_int_list_of_r2 <- c(lm_lr_int_list_of_r2, int_scramb_LassoMin_linearR2)
      int_scramb_LassoMin_linear_RMSE_train <- rmse(list_of_scrambs[[Z]], int_scramb_LassoMin_linearFittedVals)
      lm_lr_int_list_of_rmse <- c(lm_lr_int_list_of_rmse, int_scramb_LassoMin_linear_RMSE_train)
      int_scramb_LassoMin_linear_MAE_train <- mae(list_of_scrambs[[Z]], int_scramb_LassoMin_linearFittedVals)
      lm_lr_int_list_of_mae <- c(lm_lr_int_list_of_mae, int_scramb_LassoMin_linear_MAE_train)
    }
    
    # Store mean of fit and error statistics for each scrambled model predictions
    BentLassoMin_lm_R2 <- mean(unlist(lm_lr_int_list_of_r2), na.rm = TRUE)
    BentLassoMin_lm_RMSE <- mean(unlist(lm_lr_int_list_of_rmse), na.rm = TRUE)
    BentLassoMin_lm_MAE <- mean(unlist(lm_lr_int_list_of_mae), na.rm = TRUE)
    
    # Save and prep all fit and error statistics for output
    lm_lr_stats <- cbind.data.frame(BentLassoMin_lm_R2, BentLassoMin_lm_MAE, BentLassoMin_lm_RMSE)
    colnames(lm_lr_stats) <- rbind("R2", "MAE ", "RMSE")
  } 
  
  ## LASSO 1se LINEAR REGRESSION
  # Check if lasso 1se feature selection returned any descriptors before proceeding
  if (length(lambda_1se_chosen_desc) >= 1){
    # Initialize empty lists to store fit and error statistics for scrambled training set predictions made by model
    l1_lr_int_list_of_r2 <- list()
    l1_lr_int_list_of_rmse <- list()
    l1_lr_int_list_of_mae <- list()
    
    # Train a linear model using each of the scrambled training sets and calculate evaluation statistics
    for (Z in 1:length(list_of_scrambs)){
      # Generate a dataset containing current scrambled y values and properly ordered descriptors
      int_scramb_lr_lasso1se_updatedData <- cbind.data.frame(list_of_scrambs[[Z]], lasso1se_updatedData[,-3])
      colnames(int_scramb_lr_lasso1se_updatedData)[1] <- "y"
      # Specify cross validation method as leave-one-out cross validation
      train.control <- trainControl(method = "LOOCV")
      # Fit a cross validated linear model between descriptors and Tethering
      int_scramb_Lasso1se_linearModel <- train(y~.-ID-SMILES, data = int_scramb_lr_lasso1se_updatedData, method = "lm", trControl = train.control)
      int_scramb_Lasso1se_linearFittedVals <- int_scramb_Lasso1se_linearModel$finalModel$fitted.values
      # Store fit and error statistics for training set predictions made by model
      int_scramb_Lasso1se_linearR2 <- summary(int_scramb_Lasso1se_linearModel$finalModel)$r.squared
      l1_lr_int_list_of_r2 <- c(l1_lr_int_list_of_r2, int_scramb_Lasso1se_linearR2)
      int_scramb_Lasso1se_linear_RMSE_train <- rmse(list_of_scrambs[[Z]], int_scramb_Lasso1se_linearFittedVals)
      l1_lr_int_list_of_rmse <- c(l1_lr_int_list_of_rmse, int_scramb_Lasso1se_linear_RMSE_train)
      int_scramb_Lasso1se_linear_MAE_train <- mae(list_of_scrambs[[Z]], int_scramb_Lasso1se_linearFittedVals)
      l1_lr_int_list_of_mae <- c(l1_lr_int_list_of_mae, int_scramb_Lasso1se_linear_MAE_train)
    }
    
    # Store mean of fit and error statistics for each scrambled model predictions
    BentLasso1se_lm_R2 <- mean(unlist(l1_lr_int_list_of_r2), na.rm = TRUE)
    BentLasso1se_lm_RMSE <- mean(unlist(l1_lr_int_list_of_rmse), na.rm = TRUE)
    BentLasso1se_lm_MAE <- mean(unlist(l1_lr_int_list_of_mae), na.rm = TRUE)
    
    # Save and prep all fit and error statistics for output
    l1_lr_stats <- cbind.data.frame(BentLasso1se_lm_R2, BentLasso1se_lm_MAE, BentLasso1se_lm_RMSE)
    colnames(l1_lr_stats) <- rbind("R2", "MAE ", "RMSE")
  } 
  
  ## MARS LINEAR REGRESSION
  # Initialize empty lists to store fit and error statistics for scrambled training set predictions made by model
  mars_lr_int_list_of_r2 <- list()
  mars_lr_int_list_of_rmse <- list()
  mars_lr_int_list_of_mae <- list()
  
  # Train a linear model using each of the scrambled training sets and calculate evaluation statistics
  for (Z in 1:length(list_of_scrambs)){
    # Generate a dataset containing current scrambled y values and properly ordered descriptors
    int_scramb_lr_mars_updatedData <- cbind.data.frame(list_of_scrambs[[Z]], mars_updatedData[,-3])
    colnames(int_scramb_lr_mars_updatedData)[1] <- "y"
    # Specify cross validation method as leave-one-out cross validation
    train.control <- trainControl(method = "LOOCV")
    # Fit a cross validated linear model between descriptors and Tethering
    int_scramb_Mars_linearModel <- train(y~.-ID-SMILES, data = int_scramb_lr_mars_updatedData, method = "lm", trControl = train.control)
    int_scramb_Mars_linearFittedVals <- int_scramb_Mars_linearModel$finalModel$fitted.values
    # Store fit and error statistics for training set predictions made by model
    int_scramb_Mars_linearR2 <- summary(int_scramb_Mars_linearModel$finalModel)$r.squared
    mars_lr_int_list_of_r2 <- c(mars_lr_int_list_of_r2, int_scramb_Mars_linearR2)
    int_scramb_Mars_linear_RMSE_train <- rmse(list_of_scrambs[[Z]], int_scramb_Mars_linearFittedVals)
    mars_lr_int_list_of_rmse <- c(mars_lr_int_list_of_rmse, int_scramb_Mars_linear_RMSE_train)
    int_scramb_Mars_linear_MAE_train <- mae(list_of_scrambs[[Z]], int_scramb_Mars_linearFittedVals)
    mars_lr_int_list_of_mae <- c(mars_lr_int_list_of_mae, int_scramb_Mars_linear_MAE_train)
  }
  
  # Store mean of fit and error statistics for each scrambled model predictions
  BentMars_lm_R2 <- mean(unlist(mars_lr_int_list_of_r2), na.rm = TRUE)
  BentMars_lm_RMSE <- mean(unlist(mars_lr_int_list_of_rmse), na.rm = TRUE)
  BentMars_lm_MAE <- mean(unlist(mars_lr_int_list_of_mae), na.rm = TRUE)
  
  # Save and prep all fit and error statistics for output
  mars_lr_stats <- cbind.data.frame(BentMars_lm_R2, BentMars_lm_MAE, BentMars_lm_RMSE)
  colnames(mars_lr_stats) <- rbind("R2", "MAE ", "RMSE")

  ## LASSO MIN PLS
  # Check if lasso min feature selection returned multiple descriptors before proceeding
  if (length(lambda_min_chosen_desc) >= 2){
    # Initialize empty lists to store fit and error statistics for scrambled training set predictions made by model
    lm_pls_int_list_of_r2 <- list()
    lm_pls_int_list_of_rmse <- list()
    lm_pls_int_list_of_mae <- list()
    
    # Train a pls model using each of the scrambled training sets and calculate evaluation statistics
    for (Z in 1:length(list_of_scrambs)){
      # Generate a dataset containing current scrambled y values and properly ordered descriptors
      int_scramb_pls_lassoMin_updatedData <- cbind.data.frame(list_of_scrambs[[Z]], lassoMin_updatedData[,-3])
      colnames(int_scramb_pls_lassoMin_updatedData)[1] <- "y"
      # Specify cross validation method as leave-one-out cross validation
      train.control <- trainControl(method = "LOOCV")
      # Fit a cross validated pls model between descriptors and Tethering
      BentLassoMin_PLSModel <- train(y~.-ID-SMILES, data = int_scramb_pls_lassoMin_updatedData, method = "pls", trControl = train.control)
      int_scramb_LassoMin_plsModel <- BentLassoMin_PLSModel$finalModel
      int_scramb_LassoMin_plsFittedVals <- int_scramb_LassoMin_plsModel$fitted.values[,,BentLassoMin_PLSModel$bestTune$ncomp]
      # Store fit and error statistics for training set predictions made by model
      int_scramb_LassoMin_plsR2 <- (cor(list_of_scrambs[[Z]], int_scramb_LassoMin_plsFittedVals))^2
      lm_pls_int_list_of_r2 <- c(lm_pls_int_list_of_r2, int_scramb_LassoMin_plsR2)
      int_scramb_LassoMin_pls_RMSE_train <- rmse(list_of_scrambs[[Z]], int_scramb_LassoMin_plsFittedVals)
      lm_pls_int_list_of_rmse <- c(lm_pls_int_list_of_rmse, int_scramb_LassoMin_pls_RMSE_train)
      int_scramb_LassoMin_pls_MAE_train <- mae(list_of_scrambs[[Z]], int_scramb_LassoMin_plsFittedVals)
      lm_pls_int_list_of_mae <- c(lm_pls_int_list_of_mae, int_scramb_LassoMin_pls_MAE_train)
    }
    
    # Store mean of fit and error statistics for each scrambled model predictions
    BentLassoMin_pls_R2 <- mean(unlist(lm_pls_int_list_of_r2), na.rm = TRUE)
    BentLassoMin_pls_RMSE <- mean(unlist(lm_pls_int_list_of_rmse), na.rm = TRUE)
    BentLassoMin_pls_MAE <- mean(unlist(lm_pls_int_list_of_mae), na.rm = TRUE)
    
    # Save and prep all fit and error statistics for output
    lm_pls_stats <- cbind.data.frame(BentLassoMin_pls_R2, BentLassoMin_pls_MAE, BentLassoMin_pls_RMSE)
    colnames(lm_pls_stats) <- rbind("R2", "MAE ", "RMSE")
  } 
  
  ## LASSO 1SE PLS
  # Check if lasso 1se feature selection returned multiple descriptors before proceeding
  if (length(lambda_1se_chosen_desc) >= 2){
    # Initialize empty lists to store fit and error statistics for scrambled training set predictions made by model
    l1_pls_int_list_of_r2 <- list()
    l1_pls_int_list_of_rmse <- list()
    l1_pls_int_list_of_mae <- list()
    
    # Train a pls model using each of the scrambled training sets and calculate evaluation statistics
    for (Z in 1:length(list_of_scrambs)){
      # Generate a dataset containing current scrambled y values and properly ordered descriptors
      int_scramb_pls_lasso1se_updatedData <- cbind.data.frame(list_of_scrambs[[Z]], lasso1se_updatedData[,-3])
      colnames(int_scramb_pls_lasso1se_updatedData)[1] <- "y"
      # Specify cross validation method as leave-one-out cross validation
      train.control <- trainControl(method = "LOOCV")
      # Fit a cross validated pls model between descriptors and Tethering
      BentLasso1se_PLSModel <- train(y~.-ID-SMILES, data = int_scramb_pls_lasso1se_updatedData, method = "pls", trControl = train.control)
      int_scramb_Lasso1se_plsModel <- BentLasso1se_PLSModel$finalModel
      int_scramb_Lasso1se_plsFittedVals <- int_scramb_Lasso1se_plsModel$fitted.values[,,BentLasso1se_PLSModel$bestTune$ncomp]
      # Store fit and error statistics for training set predictions made by model
      int_scramb_Lasso1se_plsR2 <- (cor(list_of_scrambs[[Z]], int_scramb_Lasso1se_plsFittedVals))^2
      l1_pls_int_list_of_r2 <- c(l1_pls_int_list_of_r2, int_scramb_Lasso1se_plsR2)
      int_scramb_Lasso1se_pls_RMSE_train <- rmse(list_of_scrambs[[Z]], int_scramb_Lasso1se_plsFittedVals)
      l1_pls_int_list_of_rmse <- c(l1_pls_int_list_of_rmse, int_scramb_Lasso1se_pls_RMSE_train)
      int_scramb_Lasso1se_pls_MAE_train <- mae(list_of_scrambs[[Z]], int_scramb_Lasso1se_plsFittedVals)
      l1_pls_int_list_of_mae <- c(l1_pls_int_list_of_mae, int_scramb_Lasso1se_pls_MAE_train)
    }
    
    # Store mean of fit and error statistics for each scrambled model predictions
    BentLasso1se_pls_R2 <- mean(unlist(l1_pls_int_list_of_r2), na.rm = TRUE)
    BentLasso1se_pls_RMSE <- mean(unlist(l1_pls_int_list_of_rmse), na.rm = TRUE)
    BentLasso1se_pls_MAE <- mean(unlist(l1_pls_int_list_of_mae), na.rm = TRUE)
    
    # Save and prep all fit and error statistics for output
    l1_pls_stats <- cbind.data.frame(BentLasso1se_pls_R2, BentLasso1se_pls_MAE, BentLasso1se_pls_RMSE)
    colnames(l1_pls_stats) <- rbind("R2", "MAE ", "RMSE")
  } 

  ## MARS PLS
  # Initialize empty lists to store fit and error statistics for scrambled training set predictions made by model
  mars_pls_int_list_of_r2 <- list()
  mars_pls_int_list_of_rmse <- list()
  mars_pls_int_list_of_mae <- list()
  
  # Train a pls model using each of the scrambled training sets and calculate evaluation statistics
  for (Z in 1:length(list_of_scrambs)){
    # Generate a dataset containing current scrambled y values and properly ordered descriptors
    int_scramb_pls_masr_updatedData <- cbind.data.frame(list_of_scrambs[[Z]], mars_updatedData[,-3])
    colnames(int_scramb_pls_masr_updatedData)[1] <- "y"
    # Specify cross validation method as leave-one-out cross validation
    train.control <- trainControl(method = "LOOCV")
    # Fit a cross validated pls model between descriptors and Tethering
    BentMars_PLSModel <- train(y~.-ID-SMILES, data = int_scramb_pls_masr_updatedData, method = "pls", trControl = train.control)
    int_scramb_Mars_plsModel <- BentMars_PLSModel$finalModel
    int_scramb_Mars_plsFittedVals <- int_scramb_Mars_plsModel$fitted.values[,,BentMars_PLSModel$bestTune$ncomp]
    # Store fit and error statistics for training set predictions made by model
    int_scramb_Mars_plsR2 <- (cor(list_of_scrambs[[Z]], int_scramb_Mars_plsFittedVals))^2
    mars_pls_int_list_of_r2 <- c(mars_pls_int_list_of_r2, int_scramb_Mars_plsR2)
    int_scramb_Mars_pls_RMSE_train <- rmse(list_of_scrambs[[Z]], int_scramb_Mars_plsFittedVals)
    mars_pls_int_list_of_rmse <- c(mars_pls_int_list_of_rmse, int_scramb_Mars_pls_RMSE_train)
    int_scramb_Mars_pls_MAE_train <- mae(list_of_scrambs[[Z]], int_scramb_Mars_plsFittedVals)
    mars_pls_int_list_of_mae <- c(mars_pls_int_list_of_mae, int_scramb_Mars_pls_MAE_train)
  }
  
  # Store mean of fit and error statistics for each scrambled model predictions
  BentMars_pls_R2 <- mean(unlist(mars_pls_int_list_of_r2), na.rm = TRUE)
  BentMars_pls_RMSE <- mean(unlist(mars_pls_int_list_of_rmse), na.rm = TRUE)
  BentMars_pls_MAE <- mean(unlist(mars_pls_int_list_of_mae), na.rm = TRUE)
  
  # Save and prep all fit and error statistics for output
  mars_pls_stats <- cbind.data.frame(BentMars_pls_R2, BentMars_pls_MAE, BentMars_pls_RMSE)
  colnames(mars_pls_stats) <- rbind("R2", "MAE ", "RMSE")
  
  ## LASSO 1SE RF
  # Check if lasso 1se feature selection returned any descriptors before proceeding
  if (length(lambda_1se_chosen_desc) >= 1){
    # Initialize empty lists to store fit and error statistics for scrambled training set predictions made by model
    l1_rf_int_list_of_r2 <- list()
    l1_rf_int_list_of_rmse <- list()
    l1_rf_int_list_of_mae <- list()
    
    # Train a random forest model using each of the scrambled training sets and calculate evaluation statistics
    for (Z in 1:length(list_of_scrambs)){
      # Generate a dataset containing current scrambled y values and properly ordered descriptors
      int_scramb_rf_lasso1se_updatedData <- cbind.data.frame(list_of_scrambs[[Z]], lasso1se_updatedData[,-3])
      colnames(int_scramb_rf_lasso1se_updatedData)[1] <- "y"
      # Specify cross validation method as leave-one-out cross validation
      train.control <- trainControl(method = "LOOCV")
      # Fit a cross validated rf model between descriptors and Tethering
      BentLasso1se_RFModel <- train(y~.-ID-SMILES, data = int_scramb_rf_lasso1se_updatedData, method = 'rf', trControl = train.control, tuneLength = length(lambda_1se_chosen_desc))
      final_BentLasso1se_rfModel <- BentLasso1se_RFModel$finalModel
      int_scramb_Lasso1se_rfFittedVals <- final_BentLasso1se_rfModel$predicted
      # Store fit and error statistics for training set predictions made by model
      int_scramb_Lasso1se_rfR2 <- (cor(list_of_scrambs[[Z]], int_scramb_Lasso1se_rfFittedVals))^2
      l1_rf_int_list_of_r2 <- c(l1_rf_int_list_of_r2, int_scramb_Lasso1se_rfR2)
      int_scramb_Lasso1se_rfRMSE <- rmse(list_of_scrambs[[Z]], int_scramb_Lasso1se_rfFittedVals)
      l1_rf_int_list_of_rmse <- c(l1_rf_int_list_of_rmse, int_scramb_Lasso1se_rfRMSE)
      int_scramb_Lasso1se_rfMAE <- mae(list_of_scrambs[[Z]], int_scramb_Lasso1se_rfFittedVals)
      l1_rf_int_list_of_mae <- c(l1_rf_int_list_of_mae, int_scramb_Lasso1se_rfMAE)
    }
    
    # Store mean of fit and error statistics for each scrambled model predictions
    BentLasso1se_rf_R2 <- mean(unlist(l1_rf_int_list_of_r2), na.rm = TRUE)
    BentLasso1se_rf_RMSE <- mean(unlist(l1_rf_int_list_of_rmse), na.rm = TRUE)
    BentLasso1se_rf_MAE <- mean(unlist(l1_rf_int_list_of_mae), na.rm = TRUE)
    
    # Save and prep all fit and error statistics for output
    l1_rf_stats <- cbind.data.frame(BentLasso1se_rf_R2, BentLasso1se_rf_MAE, BentLasso1se_rf_RMSE)
    colnames(l1_rf_stats) <- rbind("R2", "MAE ", "RMSE")
  }
  
  ## LASSO MIN SVM
  # Check if lasso min feature selection returned any descriptors before proceeding
  if (length(lambda_min_chosen_desc) >= 1){
    # Initialize empty lists to store fit and error statistics for scrambled training set predictions made by model
    lm_svm_int_list_of_r2 <- list()
    lm_svm_int_list_of_rmse <- list()
    lm_svm_int_list_of_mae <- list()
    
    # Train a svm model using each of the scrambled training sets and calculate evaluation statistics
    for (Z in 1:length(list_of_scrambs)){
      # Generate a dataset containing current scrambled y values and properly ordered descriptors
      int_scramb_svm_lassoMin_updatedData <- cbind.data.frame(list_of_scrambs[[Z]], lassoMin_updatedData[,-3])
      colnames(int_scramb_svm_lassoMin_updatedData)[1] <- "y"
      # Specify cross validation method as leave-one-out cross validation
      train.control <- trainControl(method = "LOOCV")
      # Fit a cross validated svm model between descriptors and Tethering
      BentLassoMin_SVMModel <- train(y~.-ID-SMILES, data = int_scramb_svm_lassoMin_updatedData, method = 'svmRadial', trControl = train.control, scale = FALSE, tuneLength = 15)
      int_scramb_LassoMin_svmModel <- BentLassoMin_SVMModel$finalModel
      int_scramb_LassoMin_svmFittedVals <- int_scramb_LassoMin_svmModel@fitted
      # Store fit and error statistics for training set predictions made by model
      int_scramb_LassoMin_svmR2 <- (cor(list_of_scrambs[[Z]], int_scramb_LassoMin_svmFittedVals))^2
      lm_svm_int_list_of_r2 <- c(lm_svm_int_list_of_r2, int_scramb_LassoMin_svmR2)
      int_scramb_LassoMin_svm_RMSE_train <- rmse(list_of_scrambs[[Z]], int_scramb_LassoMin_svmFittedVals)
      lm_svm_int_list_of_rmse <- c(lm_svm_int_list_of_rmse, int_scramb_LassoMin_svm_RMSE_train)
      int_scramb_LassoMin_svm_MAE_train <- mae(list_of_scrambs[[Z]], int_scramb_LassoMin_svmFittedVals)
      lm_svm_int_list_of_mae <- c(lm_svm_int_list_of_mae, int_scramb_LassoMin_svm_MAE_train)
    }
    
    # Store mean of fit and error statistics for each scrambled model predictions
    BentLassoMin_svm_R2 <- mean(unlist(lm_svm_int_list_of_r2), na.rm = TRUE)
    BentLassoMin_svm_RMSE <- mean(unlist(lm_svm_int_list_of_rmse), na.rm = TRUE)
    BentLassoMin_svm_MAE <- mean(unlist(lm_svm_int_list_of_mae), na.rm = TRUE)
    
    # Save and prep all fit and error statistics for output
    lm_svm_stats <- cbind.data.frame(BentLassoMin_svm_R2, BentLassoMin_svm_MAE, BentLassoMin_svm_RMSE)
    colnames(lm_svm_stats) <- rbind("R2", "MAE ", "RMSE")
  }  
  
  ## LASSO 1SE SVM
  # Check if lasso 1se feature selection returned any descriptors before proceeding
  if (length(lambda_1se_chosen_desc) >= 1){
    # Initialize empty lists to store fit and error statistics for scrambled training set predictions made by model
    l1_svm_int_list_of_r2 <- list()
    l1_svm_int_list_of_rmse <- list()
    l1_svm_int_list_of_mae <- list()
    
    # Train a svm model using each of the scrambled training sets and calculate evaluation statistics
    for (Z in 1:length(list_of_scrambs)){
      # Generate a dataset containing current scrambled y values and properly ordered descriptors
      int_scramb_svm_lasso1se_updatedData <- cbind.data.frame(list_of_scrambs[[Z]], lasso1se_updatedData[,-3])
      colnames(int_scramb_svm_lasso1se_updatedData)[1] <- "y"
      # Specify cross validation method as leave-one-out cross validation
      train.control <- trainControl(method = "LOOCV")
      # Fit a cross validated svm model between descriptors and Tethering
      BentLasso1se_SVMModel <- train(y~.-ID-SMILES, data = int_scramb_svm_lasso1se_updatedData, method = 'svmRadial', trControl = train.control, scale = FALSE, tuneLength = 15)
      int_scramb_Lasso1se_svmModel <- BentLasso1se_SVMModel$finalModel
      int_scramb_Lasso1se_svmFittedVals <- int_scramb_Lasso1se_svmModel@fitted
      # Store fit and error statistics for training set predictions made by model
      int_scramb_Lasso1se_svmR2 <- (cor(list_of_scrambs[[Z]], int_scramb_Lasso1se_svmFittedVals))^2
      l1_svm_int_list_of_r2 <- c(l1_svm_int_list_of_r2, int_scramb_Lasso1se_svmR2)
      int_scramb_Lasso1se_svm_RMSE_train <- rmse(list_of_scrambs[[Z]], int_scramb_Lasso1se_svmFittedVals)
      l1_svm_int_list_of_rmse <- c(l1_svm_int_list_of_rmse, int_scramb_Lasso1se_svm_RMSE_train)
      int_scramb_Lasso1se_svm_MAE_train <- mae(list_of_scrambs[[Z]], int_scramb_Lasso1se_svmFittedVals)
      l1_svm_int_list_of_mae <- c(l1_svm_int_list_of_mae, int_scramb_Lasso1se_svm_MAE_train)
    }
    
    # Store mean of fit and error statistics for each scrambled model predictions
    BentLasso1se_svm_R2 <- mean(unlist(l1_svm_int_list_of_r2), na.rm = TRUE)
    BentLasso1se_svm_RMSE <- mean(unlist(l1_svm_int_list_of_rmse), na.rm = TRUE)
    BentLasso1se_svm_MAE <- mean(unlist(l1_svm_int_list_of_mae), na.rm = TRUE)
    
    # Save and prep all fit and error statistics for output
    l1_svm_stats <- cbind.data.frame(BentLasso1se_svm_R2, BentLasso1se_svm_MAE, BentLasso1se_svm_RMSE)
    colnames(l1_svm_stats) <- rbind("R2", "MAE ", "RMSE")
  }
  
  ## MARS SVM
  # Initialize empty lists to store fit and error statistics for scrambled training set predictions made by model
  mars_svm_int_list_of_r2 <- list()
  mars_svm_int_list_of_rmse <- list()
  mars_svm_int_list_of_mae <- list()
  
  # Train a svm model using each of the scrambled training sets and calculate evaluation statistics
  for (Z in 1:length(list_of_scrambs)){
    # Generate a dataset containing current scrambled y values and properly ordered descriptors
    int_scramb_svm_mars_updatedData <- cbind.data.frame(list_of_scrambs[[Z]], mars_updatedData[,-3])
    colnames(int_scramb_svm_mars_updatedData)[1] <- "y"
    # Specify cross validation method as leave-one-out cross validation
    train.control <- trainControl(method = "LOOCV")
    # Fit a cross validated svm model between descriptors and Tethering
    BentMars_SVMModel <- train(y~.-ID-SMILES, data = int_scramb_svm_mars_updatedData, method = 'svmRadial', trControl = train.control, scale = FALSE, tuneLength = 15)
    int_scramb_Mars_svmModel <- BentMars_SVMModel$finalModel
    int_scramb_Mars_svmFittedVals <- int_scramb_Mars_svmModel@fitted
    # Store fit and error statistics for training set predictions made by model
    int_scramb_Mars_svmR2 <- (cor(list_of_scrambs[[Z]], int_scramb_Mars_svmFittedVals))^2
    mars_svm_int_list_of_r2 <- c(mars_svm_int_list_of_r2, int_scramb_Mars_svmR2)
    int_scramb_Mars_svm_RMSE_train <- rmse(list_of_scrambs[[Z]], int_scramb_Mars_svmFittedVals)
    mars_svm_int_list_of_rmse <- c(mars_svm_int_list_of_rmse, int_scramb_Mars_svm_RMSE_train)
    int_scramb_Mars_svm_MAE_train <- mae(list_of_scrambs[[Z]], int_scramb_Mars_svmFittedVals)
    mars_svm_int_list_of_mae <- c(mars_svm_int_list_of_mae, int_scramb_Mars_svm_MAE_train)
  }
  
  # Store mean of fit and error statistics for each scrambled model predictions
  BentMars_svm_R2 <- mean(unlist(mars_svm_int_list_of_r2), na.rm = TRUE)
  BentMars_svm_RMSE <- mean(unlist(mars_svm_int_list_of_rmse), na.rm = TRUE)
  BentMars_svm_MAE <- mean(unlist(mars_svm_int_list_of_mae), na.rm = TRUE)
  
  # Save and prep all fit and error statistics for output
  mars_svm_stats <- cbind.data.frame(BentMars_svm_R2, BentMars_svm_MAE, BentMars_svm_RMSE)
  colnames(mars_svm_stats) <- rbind("R2", "MAE ", "RMSE")
  
  
  # Combine all calculated fit and error statistics for current training set into one dataframe
  statistics <- rbind.data.frame(lm_lr_stats, l1_lr_stats, mars_lr_stats, lm_pls_stats, l1_pls_stats, mars_pls_stats, l1_rf_stats, lm_svm_stats, l1_svm_stats, mars_svm_stats)
  new_stats <- t(statistics)
  colnames(new_stats) <- rbind("Lasso Min Linear Regression",  "Lasso 1se Linear Regression", "Mars Linear Regression", "Lasso Min PLS", "Lasso 1se PLS", "Mars PLS", "Lasso 1se RF", "Lasso Min SVM", "Lasso 1se SVM", "Mars SVM")
  # Generate output file summarizing scrambled model statistics for current training set
  fileName = paste0("BentScramblingStats", "_", fileNumbers[i], ".csv")
  write.csv(new_stats, file = fileName)
}

