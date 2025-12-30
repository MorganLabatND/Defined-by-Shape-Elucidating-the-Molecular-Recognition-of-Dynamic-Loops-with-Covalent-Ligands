setwd('/lsi/groups/amapplab/chloemarkey')

# INSTALL AND LOAD ALL NECESSARY PACKAGES AND LIBRARIES -----------------------------------------------------------------
install.packages(c('stringi', 'glmnet', 'ModelMetrics', 'pls', 'randomForest', 'kernlab', 'earth', 'prospectr'), repos="http://cran.r-project.org", lib="/lsi/groups/amapplab/chloemarkey/Anaconda_libs/")
install.packages('caret', dependencies = c("Depends", "Imports", "Suggests"), repos="http://cran.r-project.org", lib="/lsi/groups/amapplab/chloemarkey/Anaconda_libs/")
install.packages('ggplot2', dependencies = TRUE, repos="http://cran.r-project.org", lib="/lsi/groups/amapplab/chloemarkey/Anaconda_libs/")
install.packages('withr', dependencies = TRUE, repos="http://cran.r-project.org", lib="/lsi/groups/amapplab/chloemarkey/Anaconda_libs/")
install.packages('Formula', dependencies = TRUE, repos="http://cran.r-project.org", lib="/lsi/groups/amapplab/chloemarkey/Anaconda_libs/")
install.packages('plotmo', dependencies = TRUE, repos="http://cran.r-project.org", lib="/lsi/groups/amapplab/chloemarkey/Anaconda_libs/")
install.packages('plotrix', dependencies = TRUE, repos="http://cran.r-project.org", lib="/lsi/groups/amapplab/chloemarkey/Anaconda_libs/")
install.packages('TeachingDemos', dependencies = TRUE, repos="http://cran.r-project.org", lib="/lsi/groups/amapplab/chloemarkey/Anaconda_libs/")
install.packages('recipes', dependencies = TRUE, repos="http://cran.r-project.org", lib="/lsi/groups/amapplab/chloemarkey/Anaconda_libs/")
install.packages('dplyr', dependencies = TRUE, repos="http://cran.r-project.org", lib="/lsi/groups/amapplab/chloemarkey/Anaconda_libs/")
install.packages('crayon', dependencies = TRUE, repos="http://cran.r-project.org", lib="/lsi/groups/amapplab/chloemarkey/Anaconda_libs/")
install.packages('Matrix', dependencies = TRUE, repos="http://cran.r-project.org", lib="/lsi/groups/amapplab/chloemarkey/Anaconda_libs/")
install.packages('lattice', dependencies = TRUE, repos="http://cran.r-project.org", lib="/lsi/groups/amapplab/chloemarkey/Anaconda_libs/")
library(lattice, lib="/lsi/groups/amapplab/chloemarkey/Anaconda_libs/")
library(Matrix, lib="/lsi/groups/amapplab/chloemarkey/Anaconda_libs/")
library(crayon, lib="/lsi/groups/amapplab/chloemarkey/Anaconda_libs/")
library(dplyr, lib="/lsi/groups/amapplab/chloemarkey/Anaconda_libs/")
library(recipes, lib="/lsi/groups/amapplab/chloemarkey/Anaconda_libs/")
library(Formula, lib="/lsi/groups/amapplab/chloemarkey/Anaconda_libs/")
library(TeachingDemos, lib="/lsi/groups/amapplab/chloemarkey/Anaconda_libs/")
library(plotrix, lib="/lsi/groups/amapplab/chloemarkey/Anaconda_libs/")
library(plotmo, lib="/lsi/groups/amapplab/chloemarkey/Anaconda_libs/")
library(withr, lib="/lsi/groups/amapplab/chloemarkey/Anaconda_libs/")
library(ggplot2, lib="/lsi/groups/amapplab/chloemarkey/Anaconda_libs/")
library(glmnet, lib="/lsi/groups/amapplab/chloemarkey/Anaconda_libs/")
library(stringi, lib="/lsi/groups/amapplab/chloemarkey/Anaconda_libs/")
library(caret, lib="/lsi/groups/amapplab/chloemarkey/Anaconda_libs/")
library(ModelMetrics, lib="/lsi/groups/amapplab/chloemarkey/Anaconda_libs/")
library(pls, lib="/lsi/groups/amapplab/chloemarkey/Anaconda_libs/")
library(randomForest, lib="/lsi/groups/amapplab/chloemarkey/Anaconda_libs/")
library(kernlab, lib="/lsi/groups/amapplab/chloemarkey/Anaconda_libs/")
library(earth, lib="/lsi/groups/amapplab/chloemarkey/Anaconda_libs/")
library(prospectr, lib="/lsi/groups/amapplab/chloemarkey/Anaconda_libs/")

# Read in molecule data ---------------------------------------------------------------------------
# Open corresponding training/test set pairs
myTestFiles <- list.files(pattern="testSet.*csv")
myTrainFiles <- list.files(pattern="trainSet.*csv")
testDat<-lapply(myTestFiles,read.csv)
trainDat<-lapply(myTrainFiles,read.csv)
fileNumbers <- regmatches(myTrainFiles, regexpr( "\\d+", myTrainFiles))
# Open a set of ZINC molecules
zinc_chunk0s2_data <- read.csv('/lsi/groups/amapplab/chloemarkey/chunk0_set2_0330boltzmanned.csv')

# Loop through each training/test data set pair and implement machine learning workflow --------------------
for (i in 1:length(myTrainFiles)){
  # Load current training and test set
  trainingSet <- trainDat[[i]]
  trainingSet <- subset(trainingSet, select = -c(Group))
  colnames(trainingSet)[3] <- "y"
  
  testingSet <- testDat[[i]]
  testingSet <- subset(testingSet, select = -c(Group))
  colnames(testingSet)[3] <- "y"
  
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
    index <- which(colnames(testingSet) == invariant_desc[a])
    invariant_inds <- c(invariant_inds, index)
  }
  
  # Save remaining descriptors as a model matrix
  x_train = model.matrix(y ~.-ID-SMILES, data = invar_desc_removed_trainSet)
  x_train = x_train[,-1]
  
  # Save molecule identifier information and Tethering (y) values
  ID = invar_desc_removed_trainSet$ID
  SMILES = invar_desc_removed_trainSet$SMILES
  y = as.numeric(invar_desc_removed_trainSet$y)
  
  ## IDENTIFY AND REMOVE PERFECTLY CORRELATED DESCRIPTORS
  # Generate matrix containing correlations between each descriptors and Tethering
  bent_desc_tethering_cormat = cor(x_train, y, method = "kendall")
  
  # Generate matrix containing correlations between all descriptors and each other
  bent_descriptor_cormat = cor(x_train, method = "kendall")
  bent_descriptor_cormat = round(bent_descriptor_cormat, 3)
  
  # Set diagonal and upper triangle of matrix to 0 to prevent repeat analysis
  bent_descriptor_cormat[upper.tri(bent_descriptor_cormat)] <- 0
  diag(bent_descriptor_cormat) <- 0
  descriptor_list <- colnames(bent_descriptor_cormat)
  
  # Initialize empty list to store indices of descriptors to be removed
  ind_to_remove <- list()

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
      }
    }
  }
  
  # Generate ordered list of descriptor indicies removed via analysis
  new_ind_to_remove <- unique(ind_to_remove)
  new_ind_to_remove <- unlist(new_ind_to_remove)
  new_ind_to_remove <- sort(new_ind_to_remove)
  
  # Save list of descriptors removed via analysis
  removed_descriptors <- descriptor_list[new_ind_to_remove]
  
  ## IDENTIFY AND REMOVE DESCRIPTORS CORRELATED > |0.8|
  # Update correlation matrices with new descriptors removed
  x_train = x_train[,-c(new_ind_to_remove)]
  bent_desc_tethering_cormat = cor(x_train, y, method = "kendall")
  bent_descriptor_cormat = cor(x_train, method = "kendall")
  bent_descriptor_cormat = round(bent_descriptor_cormat, 3)
  bent_descriptor_cormat[upper.tri(bent_descriptor_cormat)] <- 0
  diag(bent_descriptor_cormat) <- 0
  descriptor_list <- colnames(bent_descriptor_cormat)
  
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
  
  # Save list of descriptors removed and descriptors remaining post analysis to generate updated data set
  removed_descriptors_80 <- descriptor_list[new_ind_to_remove_80]
  remaining_descriptors_80 <- x_train[,-c(new_ind_to_remove_80)]
  updated_data <- trainingSet[,c("ID", "SMILES", "y", colnames(remaining_descriptors_80))]
  
  ## SCALE DESCRIPTORS
  # Generate model matrix of training set containing remaining descriptors
  x_new = model.matrix(y ~ .-ID-SMILES, data = updated_data)
  
  # Scale all descriptor values in training set model matrix
  x_means <- colMeans(x_new)
  x_sd <- apply(x_new, 2, sd)
  ext_desc = model.matrix(Ind ~ .-SMILES, data = zinc_chunk0s2_data)
  cols_to_keep <- which(colnames(ext_desc) %in% colnames(x_new))
  cols_to_keep1 <- which(colnames(x_new) %in% colnames(ext_desc))
  x_new_ext <- ext_desc[,c(cols_to_keep)]
  x_train_dat_forExt <- x_new[,c(cols_to_keep1)]
  x_train_means_forExt <- colMeans(x_train_dat_forExt)
  x_train_sd_forExt <- apply(x_train_dat_forExt, 2, sd)
  for (K in 1:(ncol(x_new))){
    for (j in 1:(nrow(x_new))){
      scaled_val <- (x_new[j,K]-x_means[K])/(x_sd[K])
      x_new[j,K] <- scaled_val
    }
  }
  x_new[,1] = 1
  
  # Generate model matrix of test set containing remaining descriptors
  updated_test <- testingSet[,c("ID", "SMILES", "y", colnames(remaining_descriptors_80))]
  test_desc = model.matrix(y ~ .-ID-SMILES, data = updated_test)
  
  # Scale all descriptor values in test set model matrix
  for (k in 1:(ncol(test_desc))){
    for (j in 1:(nrow(test_desc))){
      scaled_val <- ((test_desc[j,k]-x_means[k]))/(x_sd[k])
      test_desc[j,k] <- scaled_val
    }
  }
  test_desc = test_desc[,-1]
  test_teth = updated_test$y
  bentTestSet <- cbind.data.frame(test_teth, test_desc)
  bentTestSet[is.na(bentTestSet)] = 0
  
  # Scale all descriptor values in ZINC set model matrix
  for (k in 1:(ncol(x_new_ext))){
    for (j in 1:(nrow(x_new_ext))){
      scaled_val <- ((x_new_ext[j,k]-x_train_means_forExt[k]))/(x_train_sd_forExt[k])
      x_new_ext[j,k] <- scaled_val
    }
  }
  x_new_ext = x_new_ext[,-1]
  bentExtSet <- x_new_ext
  bentExtSet[is.na(bentExtSet)] = 0
  
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
  
  # Generate an updated training, test, and ZINC set containing only descriptors selected by lasso lambda min model
  if (length(lambda_min_chosen_desc) >= 1){
    subset(x_new, select=c(lambda_min_chosen_desc))
    lassoMin_updatedData <- cbind.data.frame(ID, SMILES, y, subset(x_new, select=c(lambda_min_chosen_desc)))
    test_SubSet_lambdaMin <- bentTestSet[,c("test_teth", lambda_min_chosen_desc)]
    colnames(test_SubSet_lambdaMin)[1] <- "y"
    if (length(setdiff(lambda_min_chosen_desc, colnames(bentExtSet))) == 0){
      ext_SubSet_lambdaMin <- bentExtSet[,c(lambda_min_chosen_desc)]
      ext_SubSet_lambdaMin <- as.data.frame(ext_SubSet_lambdaMin)
      colnames(ext_SubSet_lambdaMin) <- c(lambda_min_chosen_desc)
    }
  }

  # Identify descriptors associated with the lasso model with a lambda value 1 standard error from the minimum
  beta <- coef(cvfit, s = "lambda.1se")
  ii <- which(beta!=0)
  lambda_1se_descriptors <- beta[ii, 1]
  lambda_1se_chosen_desc <- names(lambda_1se_descriptors)
  lambda_1se_chosen_desc <- lambda_1se_chosen_desc[-1]
  
  # Generate an updated training, test, and ZINC set containing only descriptors selected by lasso lambda 1se model
  if (length(lambda_1se_chosen_desc) >= 1){
    subset(x_new, select=c(lambda_1se_chosen_desc))
    lasso1se_updatedData <- cbind.data.frame(ID, SMILES, y, subset(x_new, select=c(lambda_1se_chosen_desc)))
    test_SubSet_lambda1se <- bentTestSet[,c("test_teth", lambda_1se_chosen_desc)]
    colnames(test_SubSet_lambda1se)[1] <- "y"
    if (length(setdiff(lambda_1se_chosen_desc, colnames(bentExtSet))) == 0){
      ext_SubSet_lambda1se <- bentExtSet[,c(lambda_1se_chosen_desc)]
      ext_SubSet_lambda1se <- as.data.frame(ext_SubSet_lambda1se)
      colnames(ext_SubSet_lambda1se) <- c(lambda_1se_chosen_desc)
    }
  }
  
  # Generate a training and test set containing all descriptors before feature selection applied
  noFS_updated_data <- cbind.data.frame(ID, SMILES, y, x_new[,-1])
  bentTestSet <- bentTestSet[,c("test_teth", colnames(remaining_descriptors_80))]
  colnames(bentTestSet)[1] <- "y"
  
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

  # Generate a training, test, and ZINC set containing only descriptors selected by mars model
  subset(x_new, select=c(mars_chosen_desc))
  mars_updatedData <- cbind.data.frame(ID, SMILES, y, subset(x_new, select=c(mars_chosen_desc)))
  test_SubSet_mars <- bentTestSet[,c("y", mars_chosen_desc)]
  if (length(setdiff(mars_chosen_desc, colnames(bentExtSet))) == 0){
    ext_SubSet_mars <- bentExtSet[,c(mars_chosen_desc)]
    ext_SubSet_mars <- as.data.frame(ext_SubSet_mars)
    colnames(ext_SubSet_mars) <- c(mars_chosen_desc)
  }
  
  ## INITIALIZE EMPTY DATAFRAMES TO HOLD MACHINE LEARNING MODEL STATS AND LABEL COLUMNS PROPERLY
  nf_lr_stats <- cbind.data.frame("NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA")
  lm_lr_stats <- cbind.data.frame("NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA")
  l1_lr_stats <- cbind.data.frame("NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA")
  mars_lr_stats <- cbind.data.frame("NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA")
  nf_pls_stats <- cbind.data.frame("NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA")
  lm_pls_stats <- cbind.data.frame("NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA")
  l1_pls_stats <- cbind.data.frame("NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA")
  mars_pls_stats <- cbind.data.frame("NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA")
  nf_rf_stats <- cbind.data.frame("NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA")
  lm_rf_stats <- cbind.data.frame("NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA")
  l1_rf_stats <- cbind.data.frame("NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA")
  mars_rf_stats <- cbind.data.frame("NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA")
  nf_svm_stats <- cbind.data.frame("NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA")
  lm_svm_stats <- cbind.data.frame("NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA") 
  l1_svm_stats <- cbind.data.frame("NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA")
  mars_svm_stats <- cbind.data.frame("NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA")
  
  colnames(nf_lr_stats) <- rbind("r2 train", "adj r2 train", "MAE train", "RMSE train", "r2 LOO", "MAE LOO", "RMSE LOO", "q2 test", "MAE test", "RMSE test", "r2 int scramb", "MAE int scramb", "RMSE int scramb")
  colnames(lm_lr_stats) <- rbind("r2 train", "adj r2 train", "MAE train", "RMSE train", "r2 LOO", "MAE LOO", "RMSE LOO", "q2 test", "MAE test", "RMSE test", "r2 int scramb", "MAE int scramb", "RMSE int scramb")
  colnames(l1_lr_stats) <- rbind("r2 train", "adj r2 train", "MAE train", "RMSE train", "r2 LOO", "MAE LOO", "RMSE LOO", "q2 test", "MAE test", "RMSE test", "r2 int scramb", "MAE int scramb", "RMSE int scramb")
  colnames(mars_lr_stats) <- rbind("r2 train", "adj r2 train", "MAE train", "RMSE train", "r2 LOO", "MAE LOO", "RMSE LOO", "q2 test", "MAE test", "RMSE test", "r2 int scramb", "MAE int scramb", "RMSE int scramb")
  colnames(nf_pls_stats) <- rbind("r2 train", "adj r2 train", "MAE train", "RMSE train", "r2 LOO", "MAE LOO", "RMSE LOO", "q2 test", "MAE test", "RMSE test", "r2 int scramb", "MAE int scramb", "RMSE int scramb")
  colnames(lm_pls_stats) <- rbind("r2 train", "adj r2 train", "MAE train", "RMSE train", "r2 LOO", "MAE LOO", "RMSE LOO", "q2 test", "MAE test", "RMSE test", "r2 int scramb", "MAE int scramb", "RMSE int scramb")
  colnames(l1_pls_stats) <- rbind("r2 train", "adj r2 train", "MAE train", "RMSE train", "r2 LOO", "MAE LOO", "RMSE LOO", "q2 test", "MAE test", "RMSE test", "r2 int scramb", "MAE int scramb", "RMSE int scramb")
  colnames(mars_pls_stats) <- rbind("r2 train", "adj r2 train", "MAE train", "RMSE train", "r2 LOO", "MAE LOO", "RMSE LOO", "q2 test", "MAE test", "RMSE test", "r2 int scramb", "MAE int scramb", "RMSE int scramb")
  colnames(nf_rf_stats) <- rbind("r2 train", "adj r2 train", "MAE train", "RMSE train", "r2 LOO", "MAE LOO", "RMSE LOO", "q2 test", "MAE test", "RMSE test", "r2 int scramb", "MAE int scramb", "RMSE int scramb")
  colnames(lm_rf_stats) <- rbind("r2 train", "adj r2 train", "MAE train", "RMSE train", "r2 LOO", "MAE LOO", "RMSE LOO", "q2 test", "MAE test", "RMSE test", "r2 int scramb", "MAE int scramb", "RMSE int scramb")
  colnames(l1_rf_stats) <- rbind("r2 train", "adj r2 train", "MAE train", "RMSE train", "r2 LOO", "MAE LOO", "RMSE LOO", "q2 test", "MAE test", "RMSE test", "r2 int scramb", "MAE int scramb", "RMSE int scramb")
  colnames(mars_rf_stats) <- rbind("r2 train", "adj r2 train", "MAE train", "RMSE train", "r2 LOO", "MAE LOO", "RMSE LOO", "q2 test", "MAE test", "RMSE test", "r2 int scramb", "MAE int scramb", "RMSE int scramb")
  colnames(nf_svm_stats) <- rbind("r2 train", "adj r2 train", "MAE train", "RMSE train", "r2 LOO", "MAE LOO", "RMSE LOO", "q2 test", "MAE test", "RMSE test", "r2 int scramb", "MAE int scramb", "RMSE int scramb")
  colnames(lm_svm_stats) <- rbind("r2 train", "adj r2 train", "MAE train", "RMSE train", "r2 LOO", "MAE LOO", "RMSE LOO", "q2 test", "MAE test", "RMSE test", "r2 int scramb", "MAE int scramb", "RMSE int scramb")
  colnames(l1_svm_stats) <- rbind("r2 train", "adj r2 train", "MAE train", "RMSE train", "r2 LOO", "MAE LOO", "RMSE LOO", "q2 test", "MAE test", "RMSE test", "r2 int scramb", "MAE int scramb", "RMSE int scramb")
  colnames(mars_svm_stats) <- rbind("r2 train", "adj r2 train", "MAE train", "RMSE train", "r2 LOO", "MAE LOO", "RMSE LOO", "q2 test", "MAE test", "RMSE test", "r2 int scramb", "MAE int scramb", "RMSE int scramb")
  
  ## NO FEATURE SELECTION LINEAR REGRESSION MODEL
  
  # Specify cross validation method as leave-one-out cross validation
  train.control <- trainControl(method = "LOOCV")
  # Fit a cross validated linear model between all scaled descriptors and Tethering
  BentNoFS_linearModel <- train(y~.-ID-SMILES, data = noFS_updated_data, method = "lm", trControl = train.control)
  
  # Store fit and error statistics for training set predictions made in the cross validation process
  BentNoFS_lm_r2_cv <- BentNoFS_linearModel$results$Rsquared
  BentNoFS_lm_MAE_cv <- BentNoFS_linearModel$results$MAE
  BentNoFS_lm_RMSE_cv <- BentNoFS_linearModel$results$RMSE
  
  # Store fit and error statistics for training set predictions made by model with optimized parameters
  BentNoFS_lm_RMSE_train <- rmse(y, BentNoFS_linearModel$finalModel$fitted.values)
  BentNoFS_lm_MAE_train <- mae(y, BentNoFS_linearModel$finalModel$fitted.values)
  BentNoFS_lm_r2_train <- summary(BentNoFS_linearModel$finalModel)$r.squared
  BentNoFS_lm_adj_r2_train <- summary(BentNoFS_linearModel$finalModel)$adj.r.squared
  
  # Use optimized model to make predictions of test set Tethering 
  BentNoFS_linearPreds <- predict(BentNoFS_linearModel$finalModel, newdata = bentTestSet)
  # Store fit and error statistics for test set predictions made by model with optimized parameters
  BentNoFS_lm_Q2_test <- (cor(bentTestSet$y, BentNoFS_linearPreds))^2
  BentNoFS_lm_RMSE_test <- rmse(bentTestSet$y, BentNoFS_linearPreds)
  BentNoFS_lm_MAE_test <- mae(bentTestSet$y, BentNoFS_linearPreds)

  # Generate a vector containing the y (Tethering) values in the training set in a randomized order
  BentNoFS_lm_scrambled_y <- sample(y, size = length(y), replace = FALSE)
  BentNoFS_lm_list_of_scrambs <- list()
  
  # Generate and save 45 scrambled list of training set y values
  for (j in 1:45){
    BentNoFS_lm_scrambled_y <- sample(y, size = length(y), replace = FALSE)
    BentNoFS_lm_list_of_scrambs[[j]] <- c(BentNoFS_lm_scrambled_y)
  }
  
  # Initialize empty lists to store fit and error statistics for scrambled training set predictions made by model with optimized parameters 
  BentNoFS_lm_list_of_r2 <- list()
  BentNoFS_lm_list_of_rmse <- list()
  BentNoFS_lm_list_of_mae <- list()
  
  # Use optimized model to make predictions of scrambled training set Tethering 
  for (j in 1:length(BentNoFS_lm_list_of_scrambs)){
    reshuffled_noFS_linear <- as.data.frame(cbind(BentNoFS_lm_list_of_scrambs[[j]], subset(noFS_updated_data, select = -c(ID, SMILES, y))))
    BentNoFS_shuffledPreds <- predict(BentNoFS_linearModel$finalModel, newdata = reshuffled_noFS_linear)
    RMSE_scramb <- rmse(BentNoFS_lm_list_of_scrambs[[j]], BentNoFS_shuffledPreds)
    BentNoFS_lm_list_of_rmse <- c(BentNoFS_lm_list_of_rmse, RMSE_scramb)
    MAE_scramb <- mae(BentNoFS_lm_list_of_scrambs[[j]], BentNoFS_shuffledPreds)
    BentNoFS_lm_list_of_mae <- c(BentNoFS_lm_list_of_mae, MAE_scramb)
    r2 <- (cor(BentNoFS_lm_list_of_scrambs[[j]], BentNoFS_shuffledPreds))^2
    BentNoFS_lm_list_of_r2 <- c(BentNoFS_lm_list_of_r2, r2)
  }
  
  # Store fit and error statistics for scrambled training set predictions made by model with optimized parameters
  BentNoFS_lm_scrambled_r2 <- mean(unlist(BentNoFS_lm_list_of_r2))
  BentNoFS_lm_scrambled_rmse <- mean(unlist(BentNoFS_lm_list_of_rmse))
  BentNoFS_lm_scrambled_mae <- mean(unlist(BentNoFS_lm_list_of_mae))
  # Store all fit and error statistics calculated for no feature selection/linear model
  nf_lr_stats <- cbind.data.frame(BentNoFS_lm_r2_train, BentNoFS_lm_adj_r2_train, BentNoFS_lm_MAE_train, BentNoFS_lm_RMSE_train, BentNoFS_lm_r2_cv, BentNoFS_lm_MAE_cv, BentNoFS_lm_RMSE_cv, BentNoFS_lm_Q2_test, BentNoFS_lm_MAE_test, BentNoFS_lm_RMSE_test, BentNoFS_lm_scrambled_r2, BentNoFS_lm_scrambled_mae, BentNoFS_lm_scrambled_rmse)
  colnames(nf_lr_stats) <- rbind("r2 train", "adj r2 train", "MAE train", "RMSE train", "r2 LOO", "MAE LOO", "RMSE LOO", "q2 test", "MAE test", "RMSE test", "r2 int scramb", "MAE int scramb", "RMSE int scramb")
  
  ## LASSO MIN LINEAR REGRESSION MODEL
  # Check if lasso min feature selection returned any descriptors before proceeding
  if (length(lambda_min_chosen_desc) >= 1){
    # Fit a cross validated linear model between lasso min selected descriptors and Tethering
    BentLassoMin_linearModel <- train(y~.-ID-SMILES, data = lassoMin_updatedData, method = "lm", trControl = train.control)
    
    # Store fit and error statistics for training set predictions made in the cross validation process
    BentLassoMin_lm_r2_cv <- BentLassoMin_linearModel$results$Rsquared
    BentLassoMin_lm_MAE_cv <- BentLassoMin_linearModel$results$MAE
    BentLassoMin_lm_RMSE_cv <- BentLassoMin_linearModel$results$RMSE
    
    # Store fit and error statistics for training set predictions made by model with optimized parameters
    BentLassoMin_lm_RMSE_train <- rmse(y, BentLassoMin_linearModel$finalModel$fitted.values)
    BentLassoMin_lm_MAE_train <- mae(y, BentLassoMin_linearModel$finalModel$fitted.values)
    BentLassoMin_lm_r2_train <- summary(BentLassoMin_linearModel$finalModel)$r.squared
    BentLassoMin_lm_adj_r2_train <- summary(BentLassoMin_linearModel$finalModel)$adj.r.squared
    BentLassoMin_trainPreds <- BentLassoMin_linearModel$finalModel$fitted.values
    
    # Use optimized model to make predictions of test set Tethering
    BentLassoMin_linearPreds <- predict(BentLassoMin_linearModel$finalModel, newdata = test_SubSet_lambdaMin)
    # Store fit and error statistics for test set predictions made by model with optimized parameters
    BentLassoMin_lm_Q2_test <- (cor(test_SubSet_lambdaMin$y, BentLassoMin_linearPreds))^2
    BentLassoMin_lm_RMSE_test <- rmse(test_SubSet_lambdaMin$y, BentLassoMin_linearPreds)
    BentLassoMin_lm_MAE_test <- mae(test_SubSet_lambdaMin$y, BentLassoMin_linearPreds)
    
    # Use optimized model to make predictions of ZINC set tethering and output predictions
    if (length(setdiff(lambda_min_chosen_desc, colnames(bentExtSet))) == 0){
      BentLassoMin_linearExtPreds <- predict(BentLassoMin_linearModel$finalModel, newdata = ext_SubSet_lambdaMin)
      extfFileName = paste0("Chunk0s2_BentExtPredsLassoMinLR", "_", fileNumbers[i], ".csv")
      write.csv(BentLassoMin_linearExtPreds, file = extfFileName)
    }
    
    # Generate a vector containing the y (Tethering) values in the training set in a randomized order
    BentLassoMin_lm_scrambled_y <- sample(y, size = length(y), replace = FALSE)
    BentLassoMin_lm_list_of_scrambs <- list()
  
    # Generate and save 45 scrambled list of training set y values
    for (J in 1:45){
      BentLassoMin_lm_scrambled_y <- sample(y, size = length(y), replace = FALSE)
      BentLassoMin_lm_list_of_scrambs[[J]] <- c(BentLassoMin_lm_scrambled_y)
    }
    
    # Initialize empty lists to store fit and error statistics for scrambled training set predictions made by model with optimized parameters 
    BentLassoMin_lm_list_of_r2 <- list()
    BentLassoMin_lm_list_of_rmse <- list()
    BentLassoMin_lm_list_of_mae <- list()
    
    # Use optimized model to make predictions of scrambled training set Tethering
    for (J in 1:length(BentLassoMin_lm_list_of_scrambs)){
      reshuffled_lassoMin_linear <- as.data.frame(cbind(BentLassoMin_lm_list_of_scrambs[[J]], subset(lassoMin_updatedData, select = -c(ID, SMILES, y))))
      BentLassoMin_lm_shuffledPreds <- predict(BentLassoMin_linearModel$finalModel, newdata = reshuffled_lassoMin_linear)
      RMSE_scramb <- rmse(BentLassoMin_lm_list_of_scrambs[[J]], BentLassoMin_lm_shuffledPreds)
      BentLassoMin_lm_list_of_rmse <- c(BentLassoMin_lm_list_of_rmse, RMSE_scramb)
      MAE_scramb <- mae(BentLassoMin_lm_list_of_scrambs[[J]], BentLassoMin_lm_shuffledPreds)
      BentLassoMin_lm_list_of_mae <- c(BentLassoMin_lm_list_of_mae, MAE_scramb)
      r2 <- (cor(BentLassoMin_lm_list_of_scrambs[[J]], BentLassoMin_lm_shuffledPreds))^2
      BentLassoMin_lm_list_of_r2 <- c(BentLassoMin_lm_list_of_r2, r2)
    }
    
    # Store fit and error statistics for scrambled training set predictions made by model with optimized parameters
    BentLassoMin_lm_scrambled_r2 <- mean(unlist(BentLassoMin_lm_list_of_r2))
    BentLassoMin_lm_scrambled_rmse <- mean(unlist(BentLassoMin_lm_list_of_rmse))
    BentLassoMin_lm_scrambled_mae <- mean(unlist(BentLassoMin_lm_list_of_mae))
    # Store all fit and error statistics calculated for lasso min/linear model
    lm_lr_stats <- cbind.data.frame(BentLassoMin_lm_r2_train, BentLassoMin_lm_adj_r2_train, BentLassoMin_lm_MAE_train, BentLassoMin_lm_RMSE_train, BentLassoMin_lm_r2_cv, BentLassoMin_lm_MAE_cv, BentLassoMin_lm_RMSE_cv, BentLassoMin_lm_Q2_test, BentLassoMin_lm_MAE_test, BentLassoMin_lm_RMSE_test, BentLassoMin_lm_scrambled_r2, BentLassoMin_lm_scrambled_mae, BentLassoMin_lm_scrambled_rmse)
    colnames(lm_lr_stats) <- rbind("r2 train", "adj r2 train", "MAE train", "RMSE train", "r2 LOO", "MAE LOO", "RMSE LOO", "q2 test", "MAE test", "RMSE test", "r2 int scramb", "MAE int scramb", "RMSE int scramb")
  }  
  
  ## LASSO 1SE LINEAR REGRESSION
  # Check if lasso 1se feature selection returned any descriptors before proceeding
  if (length(lambda_1se_chosen_desc) >= 1){
    # Fit a cross validated linear model between lasso 1se selected descriptors and Tethering
    BentLasso1se_linearModel <- train(y~.-ID-SMILES, data = lasso1se_updatedData, method = "lm", trControl = train.control)
    
    # Store fit and error statistics for training set predictions made in the cross validation process
    BentLasso1se_lm_r2_cv <- BentLasso1se_linearModel$results$Rsquared
    BentLasso1se_lm_MAE_cv <- BentLasso1se_linearModel$results$MAE
    BentLasso1se_lm_RMSE_cv <- BentLasso1se_linearModel$results$RMSE
    
    # Store fit and error statistics for training set predictions made by model with optimized parameters
    BentLasso1se_lm_RMSE_train <- rmse(y, BentLasso1se_linearModel$finalModel$fitted.values)
    BentLasso1se_lm_MAE_train <- mae(y, BentLasso1se_linearModel$finalModel$fitted.values)
    BentLasso1se_lm_r2_train <- summary(BentLasso1se_linearModel$finalModel)$r.squared
    BentLasso1se_lm_adj_r2_train <- summary(BentLasso1se_linearModel$finalModel)$adj.r.squared
    
    # Use optimized model to make predictions of test set Tethering
    BentLasso1se_linearPreds <- predict(BentLasso1se_linearModel$finalModel, newdata = test_SubSet_lambda1se)
    # Store fit and error statistics for test set predictions made by model with optimized parameters
    BentLasso1se_lm_Q2_test <- (cor(test_SubSet_lambda1se$y, BentLasso1se_linearPreds))^2
    BentLasso1se_lm_RMSE_test <- rmse(test_SubSet_lambda1se$y, BentLasso1se_linearPreds)
    BentLasso1se_lm_MAE_test <- mae(test_SubSet_lambda1se$y, BentLasso1se_linearPreds)

    # Use optimized model to make predictions of ZINC set tethering and output predictions
    if (length(setdiff(lambda_1se_chosen_desc, colnames(bentExtSet))) == 0){
      BentLasso1se_linearExtPreds <- predict(BentLasso1se_linearModel$finalModel, newdata = ext_SubSet_lambda1se)
      extfFileName = paste0("Chunk0s2_BentExtPredsLasso1seLR", "_", fileNumbers[i], ".csv")
      write.csv(BentLasso1se_linearExtPreds, file = extfFileName)
    }
    
    # Generate a vector containing the y (Tethering) values in the training set in a randomized order
    BentLasso1se_lm_scrambled_y <- sample(y, size = length(y), replace = FALSE)
    BentLasso1se_lm_list_of_scrambs <- list()
    
    # Generate and save 45 scrambled list of training set y values
    for (j in 1:45){
      BentLasso1se_lm_scrambled_y <- sample(y, size = length(y), replace = FALSE)
      BentLasso1se_lm_list_of_scrambs[[j]] <- c(BentLasso1se_lm_scrambled_y)
    }
    
    # Initialize empty lists to store fit and error statistics for scrambled training set predictions made by model with optimized parameters 
    BentLasso1se_lm_list_of_r2 <- list()
    BentLasso1se_lm_list_of_rmse <- list()
    BentLasso1se_lm_list_of_mae <- list()
    
    # Use optimized model to make predictions of scrambled training set Tethering
    for (j in 1:length(BentLasso1se_lm_list_of_scrambs)){
      BentLasso1se_lm_reshuffled <- as.data.frame(cbind(BentLasso1se_lm_list_of_scrambs[[j]], subset(lasso1se_updatedData, select = -c(ID, SMILES, y))))
      BentLasso1se_lm_shuffledPreds <- predict(BentLasso1se_linearModel$finalModel, newdata = BentLasso1se_lm_reshuffled)
      RMSE_scramb <- rmse(BentLasso1se_lm_list_of_scrambs[[j]], BentLasso1se_lm_shuffledPreds)
      BentLasso1se_lm_list_of_rmse <- c(BentLasso1se_lm_list_of_rmse, RMSE_scramb)
      MAE_scramb <- mae(BentLasso1se_lm_list_of_scrambs[[j]], BentLasso1se_lm_shuffledPreds)
      BentLasso1se_lm_list_of_mae <- c(BentLasso1se_lm_list_of_mae, MAE_scramb)
      r2 <- (cor(BentLasso1se_lm_list_of_scrambs[[j]], BentLasso1se_lm_shuffledPreds))^2
      BentLasso1se_lm_list_of_r2 <- c(BentLasso1se_lm_list_of_r2, r2)
    }
    
    # Store fit and error statistics for scrambled training set predictions made by model with optimized parameters
    BentLasso1se_lm_scrambled_r2 <- mean(unlist(BentLasso1se_lm_list_of_r2))
    BentLasso1se_lm_scrambled_rmse <- mean(unlist(BentLasso1se_lm_list_of_rmse))
    BentLasso1se_lm_scrambled_mae <- mean(unlist(BentLasso1se_lm_list_of_mae))
    # Store all fit and error statistics calculated for lasso 1se/linear model
    l1_lr_stats <- cbind.data.frame(BentLasso1se_lm_r2_train, BentLasso1se_lm_adj_r2_train, BentLasso1se_lm_MAE_train, BentLasso1se_lm_RMSE_train, BentLasso1se_lm_r2_cv, BentLasso1se_lm_MAE_cv, BentLasso1se_lm_RMSE_cv, BentLasso1se_lm_Q2_test, BentLasso1se_lm_MAE_test, BentLasso1se_lm_RMSE_test, BentLasso1se_lm_scrambled_r2, BentLasso1se_lm_scrambled_mae, BentLasso1se_lm_scrambled_rmse)
    colnames(l1_lr_stats) <- rbind("r2 train", "adj r2 train", "MAE train", "RMSE train", "r2 LOO", "MAE LOO", "RMSE LOO", "q2 test", "MAE test", "RMSE test", "r2 int scramb", "MAE int scramb", "RMSE int scramb")
  }
  
  ## MARS LINEAR REGRESSION
  # Fit a cross validated linear model between mars selected descriptors and Tethering
  BentMars_linearModel <- train(y~.-ID-SMILES, data = mars_updatedData, method = "lm", trControl = train.control)
  
  # Store fit and error statistics for training set predictions made in the cross validation process
  BentMars_lm_r2_cv <- BentMars_linearModel$results$Rsquared
  BentMars_lm_MAE_cv <- BentMars_linearModel$results$MAE
  BentMars_lm_RMSE_cv <- BentMars_linearModel$results$RMSE
    
  # Store fit and error statistics for training set predictions made by model with optimized parameters
  BentMars_lm_RMSE_train <- rmse(y, BentMars_linearModel$finalModel$fitted.values)
  BentMars_lm_MAE_train <- mae(y, BentMars_linearModel$finalModel$fitted.values)
  BentMars_lm_r2_train <- summary(BentMars_linearModel$finalModel)$r.squared
  BentMars_lm_adj_r2_train <- summary(BentMars_linearModel$finalModel)$adj.r.squared
  
  # Use optimized model to make predictions of test set Tethering
  BentMars_linearPreds <- predict(BentMars_linearModel$finalModel, newdata = test_SubSet_mars)
  # Store fit and error statistics for test set predictions made by model with optimized parameters
  BentMars_lm_Q2_test <- (cor(test_SubSet_mars$y, BentMars_linearPreds))^2
  BentMars_lm_RMSE_test <- rmse(test_SubSet_mars$y, BentMars_linearPreds)
  BentMars_lm_MAE_test <- mae(test_SubSet_mars$y, BentMars_linearPreds)

  # Generate a vector containing the y (Tethering) values in the training set in a randomized order
  BentMars_lm_scrambled_y <- sample(y, size = length(y), replace = FALSE)
  BentMars_lm_list_of_scrambs <- list()
    
  # Generate and save 45 scrambled list of training set y values
  for (j in 1:45){
    BentMars_lm_scrambled_y <- sample(y, size = length(y), replace = FALSE)
    BentMars_lm_list_of_scrambs[[j]] <- c(BentMars_lm_scrambled_y)
  }
  
  # Initialize empty lists to store fit and error statistics for scrambled training set predictions made by model with optimized parameters 
  BentMars_lm_list_of_r2 <- list()
  BentMars_lm_list_of_rmse <- list()
  BentMars_lm_list_of_mae <- list()
  
  # Use optimized model to make predictions of scrambled training set Tethering
  for (j in 1:length(BentMars_lm_list_of_scrambs)){
    BentMars_lm_reshuffled <- as.data.frame(cbind(BentMars_lm_list_of_scrambs[[j]], subset(mars_updatedData, select = -c(ID, SMILES, y))))
    BentMars_lm_shuffledPreds <- predict(BentMars_linearModel$finalModel, newdata = BentMars_lm_reshuffled)
    RMSE_scramb <- rmse(BentMars_lm_list_of_scrambs[[j]], BentMars_lm_shuffledPreds)
    BentMars_lm_list_of_rmse <- c(BentMars_lm_list_of_rmse, RMSE_scramb)
    MAE_scramb <- mae(BentMars_lm_list_of_scrambs[[j]], BentMars_lm_shuffledPreds)
    BentMars_lm_list_of_mae <- c(BentMars_lm_list_of_mae, MAE_scramb)
    r2 <- (cor(BentMars_lm_list_of_scrambs[[j]], BentMars_lm_shuffledPreds))^2
    BentMars_lm_list_of_r2 <- c(BentMars_lm_list_of_r2, r2)
  }
  
  # Store fit and error statistics for scrambled training set predictions made by model with optimized parameters
  BentMars_lm_scrambled_r2 <- mean(unlist(BentMars_lm_list_of_r2))
  BentMars_lm_scrambled_rmse <- mean(unlist(BentMars_lm_list_of_rmse))
  BentMars_lm_scrambled_mae <- mean(unlist(BentMars_lm_list_of_mae))
  # Store all fit and error statistics calculated for mars/linear model
  mars_lr_stats <- cbind.data.frame(BentMars_lm_r2_train, BentMars_lm_adj_r2_train, BentMars_lm_MAE_train, BentMars_lm_RMSE_train, BentMars_lm_r2_cv, BentMars_lm_MAE_cv, BentMars_lm_RMSE_cv, BentMars_lm_Q2_test, BentMars_lm_MAE_test, BentMars_lm_RMSE_test, BentMars_lm_scrambled_r2, BentMars_lm_scrambled_mae, BentMars_lm_scrambled_rmse)
  colnames(mars_lr_stats) <- rbind("r2 train", "adj r2 train", "MAE train", "RMSE train", "r2 LOO", "MAE LOO", "RMSE LOO", "q2 test", "MAE test", "RMSE test", "r2 int scramb", "MAE int scramb", "RMSE int scramb")
  
  ## NO FS PLS
  # Fit a cross validated PLS model between all scaled descriptors and Tethering
  BentNoFS_PLSModel <- train(y~.-ID-SMILES, data = noFS_updated_data, method = "pls", trControl = train.control)
  # Store fit and error statistics for training set predictions made in the cross validation process
  BentNoFS_pls_r2_cv <- BentNoFS_PLSModel$results$Rsquared[BentNoFS_PLSModel$bestTune$ncomp]
  BentNoFS_pls_MAE_cv <- BentNoFS_PLSModel$results$MAE[BentNoFS_PLSModel$bestTune$ncomp]
  BentNoFS_pls_RMSE_cv <- BentNoFS_PLSModel$results$RMSE[BentNoFS_PLSModel$bestTune$ncomp]
  
  # Store model with optimized parameters and training set predictions made by model
  final_BentNoFS_plsModel <- BentNoFS_PLSModel$finalModel
  BentNoFS_PLSfittedVals <- final_BentNoFS_plsModel$fitted.values[,,BentNoFS_PLSModel$bestTune$ncomp]
  
  # Store fit and error statistics for training set predictions made by model with optimized parameters
  BentNoFS_pls_r2_train <- (cor(y, BentNoFS_PLSfittedVals))^2
  n <- length(y)
  k <- ncol(x_new[,-1])
  BentNoFS_pls_adj_r2_train <- (1 - (((1-BentNoFS_pls_r2_train)*(n-1))/(n-k-1)))
  BentNoFS_pls_RMSE_train <- rmse(y, BentNoFS_PLSfittedVals)
  BentNoFS_pls_MAE_train <- mae(y, BentNoFS_PLSfittedVals)
  
  # Use optimized model to make predictions of test set Tethering
  BentNoFS_plsPreds <- predict(final_BentNoFS_plsModel, bentTestSet)
  BentNoFS_plsPreds <- BentNoFS_plsPreds[,,BentNoFS_PLSModel$bestTune$ncomp]
  # Store fit and error statistics for test set predictions made by model with optimized parameters
  BentNoFS_pls_Q2_test <- (cor(bentTestSet$y, BentNoFS_plsPreds))^2
  BentNoFS_pls_RMSE_test <- rmse(bentTestSet$y, BentNoFS_plsPreds)
  BentNoFS_pls_MAE_test <- mae(bentTestSet$y, BentNoFS_plsPreds)

  # Generate a vector containing the y (Tethering) values in the training set in a randomized order
  BentNoFS_pls_scrambled_y <- sample(y, size = length(y), replace = FALSE)
  BentNoFS_pls_list_of_scrambs <- list()
  
  # Generate and save 45 scrambled list of training set y values
  for (J in 1:45){
    BentNoFS_pls_scrambled_y <- sample(y, size = length(y), replace = FALSE)
    BentNoFS_pls_list_of_scrambs[[J]] <- c(BentNoFS_pls_scrambled_y)
  }
  
  # Initialize empty lists to store fit and error statistics for scrambled training set predictions made by model with optimized parameters 
  BentNoFS_pls_list_of_r2 <- list()
  BentNoFS_pls_list_of_adjr2 <- list()
  BentNoFS_pls_list_of_rmse <- list()
  BentNoFS_pls_list_of_mae <- list()
  
  # Use optimized model to make predictions of scrambled training set Tethering
  for (J in 1:length(BentNoFS_pls_list_of_scrambs)){
    BentNoFS_pls_reshuffled <- as.data.frame(cbind(BentNoFS_pls_list_of_scrambs[[J]], subset(noFS_updated_data, select = -c(ID, SMILES, y))))
    BentNoFS_pls_shuffledPreds <- predict(final_BentNoFS_plsModel, newdata = BentNoFS_pls_reshuffled)
    BentNoFS_pls_shuffledPreds <- BentNoFS_pls_shuffledPreds[,,BentNoFS_PLSModel$bestTune$ncomp]
    RMSE_scramb <- rmse(BentNoFS_pls_list_of_scrambs[[J]], BentNoFS_pls_shuffledPreds)
    BentNoFS_pls_list_of_rmse <- c(BentNoFS_pls_list_of_rmse, RMSE_scramb)
    MAE_scramb <- mae(BentNoFS_pls_list_of_scrambs[[J]], BentNoFS_pls_shuffledPreds)
    BentNoFS_pls_list_of_mae <- c(BentNoFS_pls_list_of_mae, MAE_scramb)
    r2 <- (cor(BentNoFS_pls_list_of_scrambs[[J]], BentNoFS_pls_shuffledPreds))^2
    BentNoFS_pls_list_of_r2 <- c(BentNoFS_pls_list_of_r2, r2)
  }
  
  # Store fit and error statistics for scrambled training set predictions made by model with optimized parameters
  BentNoFS_pls_scrambled_r2 <- mean(unlist(BentNoFS_pls_list_of_r2))
  BentNoFS_pls_scrambled_rmse <- mean(unlist(BentNoFS_pls_list_of_rmse))
  BentNoFS_pls_scrambled_mae <- mean(unlist(BentNoFS_pls_list_of_mae))

  # Store all fit and error statistics calculated for no feature selection/pls model
  nf_pls_stats <- cbind.data.frame(BentNoFS_pls_r2_train, BentNoFS_pls_adj_r2_train, BentNoFS_pls_MAE_train, BentNoFS_pls_RMSE_train, BentNoFS_pls_r2_cv, BentNoFS_pls_MAE_cv, BentNoFS_pls_RMSE_cv, BentNoFS_pls_Q2_test, BentNoFS_pls_MAE_test, BentNoFS_pls_RMSE_test, BentNoFS_pls_scrambled_r2, BentNoFS_pls_scrambled_mae, BentNoFS_pls_scrambled_rmse)
  colnames(nf_pls_stats) <- rbind("r2 train", "adj r2 train", "MAE train", "RMSE train", "r2 LOO", "MAE LOO", "RMSE LOO", "q2 test", "MAE test", "RMSE test", "r2 int scramb", "MAE int scramb", "RMSE int scramb")
  
  ## LASSO MIN PLS
  # Check if lasso min feature selection returned multiple descriptors before proceeding
  if (length(lambda_min_chosen_desc) >= 2){
    # Fit a cross validated linear model between lasso min selected descriptors and Tethering
    BentLassoMin_PLSModel <- train(y~.-ID-SMILES, data = lassoMin_updatedData, method = 'pls', trControl = train.control)
    
    # Store fit and error statistics for training set predictions made in the cross validation process
    BentLassoMin_pls_r2_cv <- BentLassoMin_PLSModel$results$Rsquared[BentLassoMin_PLSModel$bestTune$ncomp]
    BentLassoMin_pls_MAE_cv <- BentLassoMin_PLSModel$results$MAE[BentLassoMin_PLSModel$bestTune$ncomp]
    BentLassoMin_pls_RMSE_cv <- BentLassoMin_PLSModel$results$RMSE[BentLassoMin_PLSModel$bestTune$ncomp]
    
    # Store model with optimized parameters and training set predictions made by model
    final_BentLassoMin_plsModel <- BentLassoMin_PLSModel$finalModel
    BentLassoMin_PLSfittedVals <- final_BentLassoMin_plsModel$fitted.values[,,BentLassoMin_PLSModel$bestTune$ncomp]
    
    # Store fit and error statistics for training set predictions made by model with optimized parameters
    BentLassoMin_pls_r2_train <- (cor(y, BentLassoMin_PLSfittedVals))^2
    n <- length(y)
    k <- length(lambda_min_chosen_desc)
    BentLassoMin_pls_adj_r2_train <- (1 - (((1-BentLassoMin_pls_r2_train)*(n-1))/(n-k-1)))
    BentLassoMin_pls_RMSE_train <- rmse(y, BentLassoMin_PLSfittedVals)
    BentLassoMin_pls_MAE_train <- mae(y, BentLassoMin_PLSfittedVals)
    
    # Use optimized model to make predictions of test set Tethering
    BentLassoMin_plsPreds <- predict(final_BentLassoMin_plsModel, test_SubSet_lambdaMin)
    BentLassoMin_plsPreds <- BentLassoMin_plsPreds[,,BentLassoMin_PLSModel$bestTune$ncomp]
    # Store fit and error statistics for test set predictions made by model with optimized parameters
    BentLassoMin_pls_Q2_test <- (cor(test_SubSet_lambdaMin$y, BentLassoMin_plsPreds))^2
    BentLassoMin_pls_RMSE_test <- rmse(test_SubSet_lambdaMin$y, BentLassoMin_plsPreds)
    BentLassoMin_pls_MAE_test <- mae(test_SubSet_lambdaMin$y, BentLassoMin_plsPreds)
    
    # Use optimized model to make predictions of ZINC set tethering and output predictions
    if (length(setdiff(lambda_min_chosen_desc, colnames(bentExtSet))) == 0){
      BentLassoMin_plsExtPreds <- predict(final_BentLassoMin_plsModel, ext_SubSet_lambdaMin)
      BentLassoMin_plsExtPreds <- BentLassoMin_plsExtPreds[,,BentLassoMin_PLSModel$bestTune$ncomp]
      extfFileName = paste0("Chunk0s2_BentExtPredsLassoMinPLS", "_", fileNumbers[i], ".csv")
      write.csv(BentLassoMin_plsExtPreds, file = extfFileName)
    }
    
    # Generate a vector containing the y (Tethering) values in the training set in a randomized order
    BentLassoMin_pls_scrambled_y <- sample(y, size = length(y), replace = FALSE)
    BentLassoMin_pls_list_of_scrambs <- list()
    
    # Generate and save 45 scrambled list of training set y values
    for (J in 1:45){
      BentLassoMin_pls_scrambled_y <- sample(y, size = length(y), replace = FALSE)
      BentLassoMin_pls_list_of_scrambs[[J]] <- c(BentLassoMin_pls_scrambled_y)
    }
    
    # Initialize empty lists to store fit and error statistics for scrambled training set predictions made by model with optimized parameters 
    BentLassoMin_pls_list_of_r2 <- list()
    BentLassoMin_pls_list_of_rmse <- list()
    BentLassoMin_pls_list_of_mae <- list()
    
    # Use optimized model to make predictions of scrambled training set Tethering
    for (J in 1:length(BentLassoMin_pls_list_of_scrambs)){
      BentLassoMin_pls_reshuffled <- as.data.frame(cbind(BentLassoMin_pls_list_of_scrambs[[J]], subset(lassoMin_updatedData, select = -c(ID, SMILES, y))))
      BentLassoMin_pls_shuffledPreds <- predict(final_BentLassoMin_plsModel, newdata = BentLassoMin_pls_reshuffled)
      BentLassoMin_pls_shuffledPreds <- BentLassoMin_pls_shuffledPreds[,,BentLassoMin_PLSModel$bestTune$ncomp]
      RMSE_scramb <- rmse(BentLassoMin_pls_list_of_scrambs[[J]], BentLassoMin_pls_shuffledPreds)
      BentLassoMin_pls_list_of_rmse <- c(BentLassoMin_pls_list_of_rmse, RMSE_scramb)
      MAE_scramb <- mae(BentLassoMin_pls_list_of_scrambs[[J]], BentLassoMin_pls_shuffledPreds)
      BentLassoMin_pls_list_of_mae <- c(BentLassoMin_pls_list_of_mae, MAE_scramb)
      r2 <- (cor(BentLassoMin_pls_list_of_scrambs[[J]], BentLassoMin_pls_shuffledPreds))^2
      BentLassoMin_pls_list_of_r2 <- c(BentLassoMin_pls_list_of_r2, r2)
    }
    
    # Store fit and error statistics for scrambled training set predictions made by model with optimized parameters
    BentLassoMin_pls_scrambled_r2 <- mean(unlist(BentLassoMin_pls_list_of_r2))
    BentLassoMin_pls_scrambled_rmse <- mean(unlist(BentLassoMin_pls_list_of_rmse))
    BentLassoMin_pls_scrambled_mae <- mean(unlist(BentLassoMin_pls_list_of_mae))
    # Store all fit and error statistics calculated for lasso min/pls model
    lm_pls_stats <- cbind.data.frame(BentLassoMin_pls_r2_train, BentLassoMin_pls_adj_r2_train, BentLassoMin_pls_MAE_train, BentLassoMin_pls_RMSE_train, BentLassoMin_pls_r2_cv, BentLassoMin_pls_MAE_cv, BentLassoMin_pls_RMSE_cv, BentLassoMin_pls_Q2_test, BentLassoMin_pls_MAE_test, BentLassoMin_pls_RMSE_test, BentLassoMin_pls_scrambled_r2, BentLassoMin_pls_scrambled_mae, BentLassoMin_pls_scrambled_rmse)
    colnames(lm_pls_stats) <- rbind("r2 train", "adj r2 train", "MAE train", "RMSE train", "r2 LOO", "MAE LOO", "RMSE LOO", "q2 test", "MAE test", "RMSE test", "r2 int scramb", "MAE int scramb", "RMSE int scramb")
  }  
  
  ## LASSO 1SE PLS MODEL
  # Check if lasso 1se feature selection returned multiple descriptors before proceeding
  if (length(lambda_1se_chosen_desc) >= 2){
    # Fit a cross validated pls model between lasso 1se selected descriptors and Tethering
    BentLasso1se_PLSModel <- train(y~.-ID-SMILES, data = lasso1se_updatedData, method = 'pls', trControl = train.control)
    
    # Store fit and error statistics for training set predictions made in the cross validation process
    BentLasso1se_pls_r2_cv <- BentLasso1se_PLSModel$results$Rsquared[BentLasso1se_PLSModel$bestTune$ncomp]
    BentLasso1se_pls_MAE_cv <- BentLasso1se_PLSModel$results$MAE[BentLasso1se_PLSModel$bestTune$ncomp]
    BentLasso1se_pls_RMSE_cv <- BentLasso1se_PLSModel$results$RMSE[BentLasso1se_PLSModel$bestTune$ncomp]
    
    # Store model with optimized parameters and training set predictions made by model
    final_BentLasso1se_plsModel <- BentLasso1se_PLSModel$finalModel
    BentLasso1se_PLSfittedVals <- final_BentLasso1se_plsModel$fitted.values[,,BentLasso1se_PLSModel$bestTune$ncomp]
    
    # Store fit and error statistics for training set predictions made by model with optimized parameters
    BentLasso1se_pls_r2_train <- (cor(y, BentLasso1se_PLSfittedVals))^2
    n <- length(y)
    k <- length(lambda_1se_chosen_desc)
    BentLasso1se_pls_adj_r2_train <- (1 - (((1-BentLasso1se_pls_r2_train)*(n-1))/(n-k-1)))
    BentLasso1se_pls_RMSE_train <- rmse(y, BentLasso1se_PLSfittedVals)
    BentLasso1se_pls_MAE_train <- mae(y, BentLasso1se_PLSfittedVals)

    # Use optimized model to make predictions of test set Tethering
    BentLasso1se_plsPreds <- predict(final_BentLasso1se_plsModel, test_SubSet_lambda1se)
    BentLasso1se_plsPreds <- BentLasso1se_plsPreds[,,BentLasso1se_PLSModel$bestTune$ncomp]
    # Store fit and error statistics for test set predictions made by model with optimized parameters
    BentLasso1se_pls_Q2_test <- (cor(test_SubSet_lambda1se$y, BentLasso1se_plsPreds))^2
    BentLasso1se_pls_RMSE_test <- rmse(test_SubSet_lambda1se$y, BentLasso1se_plsPreds)
    BentLasso1se_pls_MAE_test <- mae(test_SubSet_lambda1se$y, BentLasso1se_plsPreds)
    
    # Use optimized model to make predictions of ZINC set Tethering
    if (length(setdiff(lambda_1se_chosen_desc, colnames(bentExtSet))) == 0){
      BentLasso1se_plsExtPreds <- predict(final_BentLasso1se_plsModel, ext_SubSet_lambda1se)
      BentLasso1se_plsExtPreds <- BentLasso1se_plsExtPreds[,,BentLasso1se_PLSModel$bestTune$ncomp]
      extfFileName = paste0("Chunk0s2_BentExtPredsLasso1sePLS", "_", fileNumbers[i], ".csv")
      write.csv(BentLasso1se_plsExtPreds, file = extfFileName)
    }
    
    # Generate a vector containing the y (Tethering) values in the training set in a randomized order
    BentLasso1se_pls_scrambled_y <- sample(y, size = length(y), replace = FALSE)
    BentLasso1se_pls_list_of_scrambs <- list()
  
    # Generate and save 45 scrambled list of training set y values
    for (J in 1:45){
      BentLasso1se_pls_scrambled_y <- sample(y, size = length(y), replace = FALSE)
      BentLasso1se_pls_list_of_scrambs[[J]] <- c(BentLasso1se_pls_scrambled_y)
    }
    
    # Initialize empty lists to store fit and error statistics for scrambled training set predictions made by model with optimized parameters 
    BentLasso1se_pls_list_of_r2 <- list()
    BentLasso1se_pls_list_of_rmse <- list()
    BentLasso1se_pls_list_of_mae <- list()
    
    # Use optimized model to make predictions of scrambled training set Tethering
    for (J in 1:length(BentLasso1se_pls_list_of_scrambs)){
      BentLasso1se_pls_reshuffled <- as.data.frame(cbind(BentLasso1se_pls_list_of_scrambs[[J]], subset(lasso1se_updatedData, select = -c(ID, SMILES, y))))
      BentLasso1se_pls_shuffledPreds <- predict(final_BentLasso1se_plsModel, newdata = BentLasso1se_pls_reshuffled)
      BentLasso1se_pls_shuffledPreds <- BentLasso1se_pls_shuffledPreds[,,BentLasso1se_PLSModel$bestTune$ncomp]
      RMSE_scramb <- rmse(BentLasso1se_pls_list_of_scrambs[[J]], BentLasso1se_pls_shuffledPreds)
      BentLasso1se_pls_list_of_rmse <- c(BentLasso1se_pls_list_of_rmse, RMSE_scramb)
      MAE_scramb <- mae(BentLasso1se_pls_list_of_scrambs[[J]], BentLasso1se_pls_shuffledPreds)
      BentLasso1se_pls_list_of_mae <- c(BentLasso1se_pls_list_of_mae, MAE_scramb)
      r2 <- (cor(BentLasso1se_pls_list_of_scrambs[[J]], BentLasso1se_pls_shuffledPreds))^2
      BentLasso1se_pls_list_of_r2 <- c(BentLasso1se_pls_list_of_r2, r2)
    }
    
    # Store fit and error statistics for scrambled training set predictions made by model with optimized parameters
    BentLasso1se_pls_scrambled_r2 <- mean(unlist(BentLasso1se_pls_list_of_r2))
    BentLasso1se_pls_scrambled_rmse <- mean(unlist(BentLasso1se_pls_list_of_rmse))
    BentLasso1se_pls_scrambled_mae <- mean(unlist(BentLasso1se_pls_list_of_mae))
    # Store all fit and error statistics calculated for lasso 1se/pls model
    l1_pls_stats <- cbind.data.frame(BentLasso1se_pls_r2_train, BentLasso1se_pls_adj_r2_train, BentLasso1se_pls_MAE_train, BentLasso1se_pls_RMSE_train, BentLasso1se_pls_r2_cv, BentLasso1se_pls_MAE_cv, BentLasso1se_pls_RMSE_cv, BentLasso1se_pls_Q2_test, BentLasso1se_pls_MAE_test, BentLasso1se_pls_RMSE_test, BentLasso1se_pls_scrambled_r2, BentLasso1se_pls_scrambled_mae, BentLasso1se_pls_scrambled_rmse)
    colnames(l1_pls_stats) <- rbind("r2 train", "adj r2 train", "MAE train", "RMSE train", "r2 LOO", "MAE LOO", "RMSE LOO", "q2 test", "MAE test", "RMSE test", "r2 int scramb", "MAE int scramb", "RMSE int scramb")
  }
  
  ## MARS PLS MODEL
  # Check if mars feature selection returned any descriptors before proceeding
  if (length(mars_chosen_desc) > 1){
    # Fit a cross validated pls model between mars selected descriptors and Tethering
    BentMars_PLSModel <- train(y~.-ID-SMILES, data = mars_updatedData, method = 'pls', trControl = train.control)
    
    # Store fit and error statistics for training set predictions made in the cross validation process
    BentMars_pls_r2_cv <- BentMars_PLSModel$results$Rsquared[BentMars_PLSModel$bestTune$ncomp]
    BentMars_pls_MAE_cv <- BentMars_PLSModel$results$MAE[BentMars_PLSModel$bestTune$ncomp]
    BentMars_pls_RMSE_cv <- BentMars_PLSModel$results$RMSE[BentMars_PLSModel$bestTune$ncomp]
    
    # Store model with optimized parameters and training set predictions made by model
    final_BentMars_plsModel <- BentMars_PLSModel$finalModel
    BentMars_PLSfittedVals <- final_BentMars_plsModel$fitted.values[,,BentMars_PLSModel$bestTune$ncomp]
    
    # Store fit and error statistics for training set predictions made by model with optimized parameters
    BentMars_pls_r2_train <- (cor(y, BentMars_PLSfittedVals))^2
    n <- length(y)
    k <- length(mars_chosen_desc)
    BentMars_pls_adj_r2_train <- (1 - (((1-BentMars_pls_r2_train)*(n-1))/(n-k-1)))
    BentMars_pls_RMSE_train <- rmse(y, BentMars_PLSfittedVals)
    BentMars_pls_MAE_train <- mae(y, BentMars_PLSfittedVals)

    # Use optimized model to make predictions of test set Tethering
    BentMars_plsPreds <- predict(final_BentMars_plsModel, test_SubSet_mars)
    BentMars_plsPreds <- BentMars_plsPreds[,,BentMars_PLSModel$bestTune$ncomp]
    # Store fit and error statistics for test set predictions made by model with optimized parameters
    BentMars_pls_Q2_test <- (cor(test_SubSet_mars$y, BentMars_plsPreds))^2
    BentMars_pls_RMSE_test <- rmse(test_SubSet_mars$y, BentMars_plsPreds)
    BentMars_pls_MAE_test <- mae(test_SubSet_mars$y, BentMars_plsPreds)

    # Generate a vector containing the y (Tethering) values in the training set in a randomized order
    BentMars_pls_scrambled_y <- sample(y, size = length(y), replace = FALSE)
    BentMars_pls_list_of_scrambs <- list()
    
    # Generate and save 45 scrambled list of training set y values
    for (J in 1:45){
      BentMars_pls_scrambled_y <- sample(y, size = length(y), replace = FALSE)
      BentMars_pls_list_of_scrambs[[J]] <- c(BentMars_pls_scrambled_y)
    }
    
    # Initialize empty lists to store fit and error statistics for scrambled training set predictions made by model with optimized parameters 
    BentMars_pls_list_of_r2 <- list()
    BentMars_pls_list_of_rmse <- list()
    BentMars_pls_list_of_mae <- list()
    
    # Use optimized model to make predictions of scrambled training set Tethering
    for (J in 1:length(BentMars_pls_list_of_scrambs)){
      BentMars_pls_reshuffled <- as.data.frame(cbind(BentMars_pls_list_of_scrambs[[J]], subset(mars_updatedData, select = -c(ID, SMILES, y))))
      BentMars_pls_shuffledPreds <- predict(final_BentMars_plsModel, newdata = BentMars_pls_reshuffled)
      BentMars_pls_shuffledPreds <- BentMars_pls_shuffledPreds[,,BentMars_PLSModel$bestTune$ncomp]
      RMSE_scramb <- rmse(BentMars_pls_list_of_scrambs[[J]], BentMars_pls_shuffledPreds)
      BentMars_pls_list_of_rmse <- c(BentMars_pls_list_of_rmse, RMSE_scramb)
      MAE_scramb <- mae(BentMars_pls_list_of_scrambs[[J]], BentMars_pls_shuffledPreds)
      BentMars_pls_list_of_mae <- c(BentMars_pls_list_of_mae, MAE_scramb)
      r2 <- (cor(BentMars_pls_list_of_scrambs[[J]], BentMars_pls_shuffledPreds))^2
      BentMars_pls_list_of_r2 <- c(BentMars_pls_list_of_r2, r2)
    }
    
    # Store fit and error statistics for scrambled training set predictions made by model with optimized parameters
    BentMars_pls_scrambled_r2 <- mean(unlist(BentMars_pls_list_of_r2))
    BentMars_pls_scrambled_rmse <- mean(unlist(BentMars_pls_list_of_rmse))
    BentMars_pls_scrambled_mae <- mean(unlist(BentMars_pls_list_of_mae))
    # Store all fit and error statistics calculated for mars/pls model
    mars_pls_stats <- cbind.data.frame(BentMars_pls_r2_train, BentMars_pls_adj_r2_train, BentMars_pls_MAE_train, BentMars_pls_RMSE_train, BentMars_pls_r2_cv, BentMars_pls_MAE_cv, BentMars_pls_RMSE_cv, BentMars_pls_Q2_test, BentMars_pls_MAE_test, BentMars_pls_RMSE_test, BentMars_pls_scrambled_r2, BentMars_pls_scrambled_mae, BentMars_pls_scrambled_rmse)
    colnames(mars_pls_stats) <- rbind("r2 train", "adj r2 train", "MAE train", "RMSE train", "r2 LOO", "MAE LOO", "RMSE LOO", "q2 test", "MAE test", "RMSE test", "r2 int scramb", "MAE int scramb", "RMSE int scramb")
  }
  
  ## NO FS RF MODEL
  # Fit a cross validated random forest model between all scaled descriptors and Tethering
  BentNoFS_RFModel <- train(y~.-ID-SMILES, data = noFS_updated_data, method = 'rf', trControl = train.control, tuneLength = ncol(x_new[,-1]))
  final_BentNoFS_rfModel <- BentNoFS_RFModel$finalModel

  # Store fit and error statistics for training set predictions made in the cross validation process
  BentNoFS_rf_r2_cv <- BentNoFS_RFModel$results$Rsquared[BentNoFS_RFModel$bestTune$mtry]
  BentNoFS_rf_MAE_cv <- BentNoFS_RFModel$results$MAE[BentNoFS_RFModel$bestTune$mtry]
  BentNoFS_rf_RMSE_cv <- BentNoFS_RFModel$results$RMSE[BentNoFS_RFModel$bestTune$mtry]
  
  # Store fit and error statistics for training set predictions made by model with optimized parameters
  BentNoFS_rf_r2_train <- (cor(final_BentNoFS_rfModel$predicted, y))^2
  n <- length(y)
  k <- ncol(x_new[,-1])
  BentNoFS_rf_adj_r2_train <- (1 - (((1-BentNoFS_rf_r2_train)*(n-1))/(n-k-1)))
  BentNoFS_rf_MAE_train <- mae(y, final_BentNoFS_rfModel$predicted)
  BentNoFS_rf_RMSE_train <- rmse(y, final_BentNoFS_rfModel$predicted)
  
  # Use optimized model to make predictions of test set Tethering
  BentNoFS_rfExtPreds <- predict(final_BentNoFS_rfModel, newdata = bentTestSet)
  # Store fit and error statistics for test set predictions made by model with optimized parameters
  BentNoFS_rf_Q2_test <- (cor(bentTestSet$y, BentNoFS_rfExtPreds))^2
  BentNoFS_rf_RMSE_test <- rmse(bentTestSet$y, BentNoFS_rfExtPreds)
  BentNoFS_rf_MAE_test <- mae(bentTestSet$y, BentNoFS_rfExtPreds)

  # Generate a vector containing the y (Tethering) values in the training set in a randomized order
  BentNoFS_rf_scrambled_y <- sample(y, size = length(y), replace = FALSE)
  BentNoFS_rf_list_of_scrambs <- list()
  
  # Generate and save 45 scrambled list of training set y values
  for (j in 1:45){
    BentNoFS_rf_scrambled_y <- sample(y, size = length(y), replace = FALSE)
    BentNoFS_rf_list_of_scrambs[[j]] <- c(BentNoFS_rf_scrambled_y)
  }
  
  # Initialize empty lists to store fit and error statistics for scrambled training set predictions made by model with optimized parameters 
  BentNoFS_rf_list_of_r2 <- list()
  BentNoFS_rf_list_of_rmse <- list()
  BentNoFS_rf_list_of_mae <- list()
  
  # Use optimized model to make predictions of scrambled training set Tethering
  for (j in 1:length(BentNoFS_rf_list_of_scrambs)){
    BentNoFS_rf_reshuffled <- as.data.frame(cbind(BentNoFS_rf_list_of_scrambs[[j]], subset(noFS_updated_data, select = -c(ID, SMILES, y))))
    BentNoFS_rf_shuffledPreds <- predict(final_BentNoFS_rfModel, newdata = BentNoFS_rf_reshuffled)
    RMSE_scramb <- rmse(BentNoFS_rf_list_of_scrambs[[j]], BentNoFS_rf_shuffledPreds)
    BentNoFS_rf_list_of_rmse <- c(BentNoFS_rf_list_of_rmse, RMSE_scramb)
    MAE_scramb <- mae(BentNoFS_rf_list_of_scrambs[[j]], BentNoFS_rf_shuffledPreds)
    BentNoFS_rf_list_of_mae <- c(BentNoFS_rf_list_of_mae, MAE_scramb)
    r2 <- (cor(BentNoFS_rf_list_of_scrambs[[j]], BentNoFS_rf_shuffledPreds))^2
    BentNoFS_rf_list_of_r2 <- c(BentNoFS_rf_list_of_r2, r2)
  }
  
  # Store fit and error statistics for scrambled training set predictions made by model with optimized parameters
  BentNoFS_rf_scrambled_r2 <- mean(unlist(BentNoFS_rf_list_of_r2))
  BentNoFS_rf_scrambled_rmse <- mean(unlist(BentNoFS_rf_list_of_rmse))
  BentNoFS_rf_scrambled_mae <- mean(unlist(BentNoFS_rf_list_of_mae))
  # Store all fit and error statistics calculated for no feature selection/random forest model
  nf_rf_stats <- cbind.data.frame(BentNoFS_rf_r2_train, BentNoFS_rf_adj_r2_train, BentNoFS_rf_MAE_train, BentNoFS_rf_RMSE_train, BentNoFS_rf_r2_cv, BentNoFS_rf_MAE_cv, BentNoFS_rf_RMSE_cv, BentNoFS_rf_Q2_test, BentNoFS_rf_MAE_test, BentNoFS_rf_RMSE_test, BentNoFS_rf_scrambled_r2, BentNoFS_rf_scrambled_mae, BentNoFS_rf_scrambled_rmse)
  colnames(nf_rf_stats) <- rbind("r2 train", "adj r2 train", "MAE train", "RMSE train", "r2 LOO", "MAE LOO", "RMSE LOO", "q2 test", "MAE test", "RMSE test", "r2 int scramb", "MAE int scramb", "RMSE int scramb")
  
  ## LASSO MIN RF MODEL
  # Check if lasso min feature selection returned any descriptors before proceeding
  if (length(lambda_min_chosen_desc) >= 1){
    # Fit a cross validated random forest model between lasso min selected descriptors and Tethering
    BentLassoMin_RFModel <- train(y~.-ID-SMILES, data = lassoMin_updatedData, method = 'rf', trControl = train.control, tuneLength = length(lambda_min_chosen_desc))
    final_BentLassoMin_rfModel <- BentLassoMin_RFModel$finalModel
  
    # Store fit and error statistics for training set predictions made in the cross validation process
    BentLassoMin_rf_r2_cv <- BentLassoMin_RFModel$results$Rsquared[BentLassoMin_RFModel$bestTune$mtry]
    BentLassoMin_rf_MAE_cv <- BentLassoMin_RFModel$results$MAE[BentLassoMin_RFModel$bestTune$mtry]
    BentLassoMin_rf_RMSE_cv <- BentLassoMin_RFModel$results$RMSE[BentLassoMin_RFModel$bestTune$mtry]
    
    # Store fit and error statistics for training set predictions made by model with optimized parameters
    BentLassoMin_rf_r2_train <- (cor(final_BentLassoMin_rfModel$predicted, y))^2
    n <- length(y)
    k <- length(lambda_min_chosen_desc)
    BentLassoMin_rf_adj_r2_train <- (1 - (((1-BentLassoMin_rf_r2_train)*(n-1))/(n-k-1)))
    BentLassoMin_rf_MAE_train <- mae(y, final_BentLassoMin_rfModel$predicted)
    BentLassoMin_rf_RMSE_train <- rmse(y, final_BentLassoMin_rfModel$predicted)
    
    # Use optimized model to make predictions of test set Tethering
    BentLambdaMin_rfExtPreds <- predict(final_BentLassoMin_rfModel, newdata = test_SubSet_lambdaMin)
    # Store fit and error statistics for test set predictions made by model with optimized parameters
    BentLassoMin_rf_Q2_test <- (cor(test_SubSet_lambdaMin$y, BentLambdaMin_rfExtPreds))^2
    BentLassoMin_rf_RMSE_test <- rmse(test_SubSet_lambdaMin$y, BentLambdaMin_rfExtPreds)
    BentLassoMin_rf_MAE_test <- mae(test_SubSet_lambdaMin$y, BentLambdaMin_rfExtPreds)
  
    # Generate a vector containing the y (Tethering) values in the training set in a randomized order
    BentLassoMin_rf_scrambled_y <- sample(y, size = length(y), replace = FALSE)
    BentLassoMin_rf_list_of_scrambs <- list()
    
    # Generate and save 45 scrambled list of training set y values
    for (J in 1:45){
      BentLassoMin_rf_scrambled_y <- sample(y, size = length(y), replace = FALSE)
      BentLassoMin_rf_list_of_scrambs[[J]] <- c(BentLassoMin_rf_scrambled_y)
    }
    
    # Initialize empty lists to store fit and error statistics for scrambled training set predictions made by model with optimized parameters 
    BentLassoMin_rf_list_of_r2 <- list()
    BentLassoMin_rf_list_of_rmse <- list()
    BentLassoMin_rf_list_of_mae <- list()
    
    # Use optimized model to make predictions of scrambled training set Tethering
    for (J in 1:length(BentLassoMin_rf_list_of_scrambs)){
      BentLassoMin_rf_reshuffled <- as.data.frame(cbind(BentLassoMin_rf_list_of_scrambs[[J]], subset(lassoMin_updatedData, select = -c(ID, SMILES, y))))
      BentLassoMin_rf_shuffledPreds <- predict(final_BentLassoMin_rfModel, newdata = BentLassoMin_rf_reshuffled)
      RMSE_scramb <- rmse(BentLassoMin_rf_list_of_scrambs[[J]], BentLassoMin_rf_shuffledPreds)
      BentLassoMin_rf_list_of_rmse <- c(BentLassoMin_rf_list_of_rmse, RMSE_scramb)
      MAE_scramb <- mae(BentLassoMin_rf_list_of_scrambs[[J]], BentLassoMin_rf_shuffledPreds)
      BentLassoMin_rf_list_of_mae <- c(BentLassoMin_rf_list_of_mae, MAE_scramb)
      r2 <- (cor(BentLassoMin_rf_list_of_scrambs[[J]], BentLassoMin_rf_shuffledPreds))^2
      BentLassoMin_rf_list_of_r2 <- c(BentLassoMin_rf_list_of_r2, r2)
    }
    
    # Store fit and error statistics for scrambled training set predictions made by model with optimized parameters
    BentLassoMin_rf_scrambled_r2 <- mean(unlist(BentLassoMin_rf_list_of_r2))
    BentLassoMin_rf_scrambled_rmse <- mean(unlist(BentLassoMin_rf_list_of_rmse))
    BentLassoMin_rf_scrambled_mae <- mean(unlist(BentLassoMin_rf_list_of_mae))
    # Store all fit and error statistics calculated for lasso min/random forest model
    lm_rf_stats <- cbind.data.frame(BentLassoMin_rf_r2_train, BentLassoMin_rf_adj_r2_train, BentLassoMin_rf_MAE_train, BentLassoMin_rf_RMSE_train, BentLassoMin_rf_r2_cv, BentLassoMin_rf_MAE_cv, BentLassoMin_rf_RMSE_cv, BentLassoMin_rf_Q2_test, BentLassoMin_rf_MAE_test, BentLassoMin_rf_RMSE_test, BentLassoMin_rf_scrambled_r2, BentLassoMin_rf_scrambled_mae, BentLassoMin_rf_scrambled_rmse)
    colnames(lm_rf_stats) <- rbind("r2 train", "adj r2 train", "MAE train", "RMSE train", "r2 LOO", "MAE LOO", "RMSE LOO", "q2 test", "MAE test", "RMSE test", "r2 int scramb", "MAE int scramb", "RMSE int scramb")
  }  
  
  ## LASSO 1SE RF
  # Check if lasso 1se feature selection returned any descriptors before proceeding
  if (length(lambda_1se_chosen_desc) >= 1){
    # Fit a cross validated random forest model between lasso 1se selected descriptors and Tethering
    BentLasso1se_RFModel <- train(y~.-ID-SMILES, data = lasso1se_updatedData, method = 'rf', trControl = train.control, tuneLength = length(lambda_1se_chosen_desc))
    final_BentLasso1se_rfModel <- BentLasso1se_RFModel$finalModel
    
    # Store fit and error statistics for training set predictions made in the cross validation process
    BentLasso1se_rf_r2_cv <- BentLasso1se_RFModel$results$Rsquared[BentLasso1se_RFModel$bestTune$mtry]
    BentLasso1se_rf_MAE_cv <- BentLasso1se_RFModel$results$MAE[BentLasso1se_RFModel$bestTune$mtry]
    BentLasso1se_rf_RMSE_cv <- BentLasso1se_RFModel$results$RMSE[BentLasso1se_RFModel$bestTune$mtry]
    
    # Store fit and error statistics for training set predictions made by model with optimized parameters
    BentLasso1se_rf_r2_train <- (cor(final_BentLasso1se_rfModel$predicted, y))^2
    n <- length(y)
    k <- length(lambda_1se_chosen_desc)
    BentLasso1se_rf_adj_r2_train <- (1 - (((1-BentLasso1se_rf_r2_train)*(n-1))/(n-k-1)))
    BentLasso1se_rf_MAE_train <- mae(y, final_BentLasso1se_rfModel$predicted)
    BentLasso1se_rf_RMSE_train <- rmse(y, final_BentLasso1se_rfModel$predicted)
    
    # Use optimized model to make predictions of test set Tethering
    BentLambda1se_rfExtPreds <- predict(final_BentLasso1se_rfModel, newdata = test_SubSet_lambda1se)
    # Store fit and error statistics for test set predictions made by model with optimized parameters
    BentLasso1se_rf_Q2_test <- (cor(test_SubSet_lambda1se$y, BentLambda1se_rfExtPreds))^2
    BentLasso1se_rf_RMSE_test <- rmse(test_SubSet_lambda1se$y, BentLambda1se_rfExtPreds)
    BentLasso1se_rf_MAE_test <- mae(test_SubSet_lambda1se$y, BentLambda1se_rfExtPreds)

    # Use optimized model to make predictions of ZINC set Tethering
    if (length(setdiff(lambda_1se_chosen_desc, colnames(bentExtSet))) == 0){
      print("making Lasso 1se ext pred RF")
      BentLasso1se_rfExtPreds <- predict(final_BentLasso1se_rfModel, newdata = ext_SubSet_lambda1se)
      extfFileName = paste0("Chunk0s2_BentExtPredsLasso1seRF", "_", fileNumbers[i], ".csv")
      write.csv(BentLasso1se_rfExtPreds, file = extfFileName)
    }
    
    # Generate a vector containing the y (Tethering) values in the training set in a randomized order
    BentLasso1se_rf_scrambled_y <- sample(y, size = length(y), replace = FALSE)
    BentLasso1se_rf_list_of_scrambs <- list()
    
    # Generate and save 45 scrambled list of training set y values
    for (J in 1:45){
      BentLasso1se_rf_scrambled_y <- sample(y, size = length(y), replace = FALSE)
      BentLasso1se_rf_list_of_scrambs[[J]] <- c(BentLasso1se_rf_scrambled_y)
    }
    
    # Initialize empty lists to store fit and error statistics for scrambled training set predictions made by model with optimized parameters 
    BentLasso1se_rf_list_of_r2 <- list()
    BentLasso1se_rf_list_of_rmse <- list()
    BentLasso1se_rf_list_of_mae <- list()
    
    # Use optimized model to make predictions of scrambled training set Tethering
    for (J in 1:length(BentLasso1se_rf_list_of_scrambs)){
      BentLasso1se_rf_reshuffled <- as.data.frame(cbind(BentLasso1se_rf_list_of_scrambs[[J]], subset(lasso1se_updatedData, select = -c(ID, SMILES, y))))
      BentLasso1se_rf_shuffledPreds <- predict(final_BentLasso1se_rfModel, newdata = BentLasso1se_rf_reshuffled)
      RMSE_scramb <- rmse(BentLasso1se_rf_list_of_scrambs[[J]], BentLasso1se_rf_shuffledPreds)
      BentLasso1se_rf_list_of_rmse <- c(BentLasso1se_rf_list_of_rmse, RMSE_scramb)
      MAE_scramb <- mae(BentLasso1se_rf_list_of_scrambs[[J]], BentLasso1se_rf_shuffledPreds)
      BentLasso1se_rf_list_of_mae <- c(BentLasso1se_rf_list_of_mae, MAE_scramb)
      r2 <- (cor(BentLasso1se_rf_list_of_scrambs[[J]], BentLasso1se_rf_shuffledPreds))^2
      BentLasso1se_rf_list_of_r2 <- c(BentLasso1se_rf_list_of_r2, r2)
    }
    
    # Store fit and error statistics for scrambled training set predictions made by model with optimized parameters
    BentLasso1se_rf_scrambled_r2 <- mean(unlist(BentLasso1se_rf_list_of_r2))
    BentLasso1se_rf_scrambled_rmse <- mean(unlist(BentLasso1se_rf_list_of_rmse))
    BentLasso1se_rf_scrambled_mae <- mean(unlist(BentLasso1se_rf_list_of_mae))
    # Store all fit and error statistics calculated for lasso 1se/random forest model
    l1_rf_stats <- cbind.data.frame(BentLasso1se_rf_r2_train, BentLasso1se_rf_adj_r2_train, BentLasso1se_rf_MAE_train, BentLasso1se_rf_RMSE_train, BentLasso1se_rf_r2_cv, BentLasso1se_rf_MAE_cv, BentLasso1se_rf_RMSE_cv, BentLasso1se_rf_Q2_test, BentLasso1se_rf_MAE_test, BentLasso1se_rf_RMSE_test, BentLasso1se_rf_scrambled_r2, BentLasso1se_rf_scrambled_mae, BentLasso1se_rf_scrambled_rmse)
    colnames(l1_rf_stats) <- rbind("r2 train", "adj r2 train", "MAE train", "RMSE train", "r2 LOO", "MAE LOO", "RMSE LOO", "q2 test", "MAE test", "RMSE test", "r2 int scramb", "MAE int scramb", "RMSE int scramb")
  }
  
  ## MARS RF
  # Fit a cross validated random forest model between mars selected descriptors and Tethering
  BentMars_RFModel <- train(y~.-ID-SMILES, data = mars_updatedData, method = 'rf', trControl = train.control, tuneLength = length(mars_chosen_desc))
  final_BentMars_rfModel <- BentMars_RFModel$finalModel
  
  # Store fit and error statistics for training set predictions made in the cross validation process
  BentMars_rf_r2_cv <- BentMars_RFModel$results$Rsquared[BentMars_RFModel$bestTune$mtry]
  BentMars_rf_MAE_cv <- BentMars_RFModel$results$MAE[BentMars_RFModel$bestTune$mtry]
  BentMars_rf_RMSE_cv <- BentMars_RFModel$results$RMSE[BentMars_RFModel$bestTune$mtry]
  
  # Store fit and error statistics for training set predictions made by model with optimized parameters
  BentMars_rf_r2_train <- (cor(final_BentMars_rfModel$predicted, y))^2
  n <- length(y)
  k <- length(mars_chosen_desc)
  BentMars_rf_adj_r2_train <- (1 - (((1-BentMars_rf_r2_train)*(n-1))/(n-k-1)))
  BentMars_rf_MAE_train <- mae(y, final_BentMars_rfModel$predicted)
  BentMars_rf_RMSE_train <- rmse(y, final_BentMars_rfModel$predicted)
  
  # Use optimized model to make predictions of test set Tethering
  BentMars_rfExtPreds <- predict(final_BentMars_rfModel, newdata = test_SubSet_mars)
  # Store fit and error statistics for test set predictions made by model with optimized parameters
  BentMars_rf_Q2_test <- (cor(test_SubSet_mars$y, BentMars_rfExtPreds))^2
  BentMars_rf_RMSE_test <- rmse(test_SubSet_mars$y, BentMars_rfExtPreds)
  BentMars_rf_MAE_test <- mae(test_SubSet_mars$y, BentMars_rfExtPreds)

  # Generate a vector containing the y (Tethering) values in the training set in a randomized order
  BentMars_rf_scrambled_y <- sample(y, size = length(y), replace = FALSE)
  BentMars_rf_list_of_scrambs <- list()
  
  # Generate and save 45 scrambled list of training set y values
  for (J in 1:45){
    BentMars_rf_scrambled_y <- sample(y, size = length(y), replace = FALSE)
    BentMars_rf_list_of_scrambs[[J]] <- c(BentMars_rf_scrambled_y)
  }
  
  # Initialize empty lists to store fit and error statistics for scrambled training set predictions made by model with optimized parameters 
  BentMars_rf_list_of_r2 <- list()
  BentMars_rf_list_of_rmse <- list()
  BentMars_rf_list_of_mae <- list()
  
  # Use optimized model to make predictions of scrambled training set Tethering
  for (J in 1:length(BentMars_rf_list_of_scrambs)){
    BentMars_rf_reshuffled <- as.data.frame(cbind(BentMars_rf_list_of_scrambs[[J]], subset(mars_updatedData, select = -c(ID, SMILES, y))))
    BentMars_rf_shuffledPreds <- predict(final_BentMars_rfModel, newdata = BentMars_rf_reshuffled)
    RMSE_scramb <- rmse(BentMars_rf_list_of_scrambs[[J]], BentMars_rf_shuffledPreds)
    BentMars_rf_list_of_rmse <- c(BentMars_rf_list_of_rmse, RMSE_scramb)
    MAE_scramb <- mae(BentMars_rf_list_of_scrambs[[J]], BentMars_rf_shuffledPreds)
    BentMars_rf_list_of_mae <- c(BentMars_rf_list_of_mae, MAE_scramb)
    r2 <- (cor(BentMars_rf_list_of_scrambs[[J]], BentMars_rf_shuffledPreds))^2
    BentMars_rf_list_of_r2 <- c(BentMars_rf_list_of_r2, r2)
  }
  
  # Store fit and error statistics for scrambled training set predictions made by model with optimized parameters
  BentMars_rf_scrambled_r2 <- mean(unlist(BentMars_rf_list_of_r2))
  BentMars_rf_scrambled_rmse <- mean(unlist(BentMars_rf_list_of_rmse))
  BentMars_rf_scrambled_mae <- mean(unlist(BentMars_rf_list_of_mae))
  # Store all fit and error statistics calculated for mars/random forest model
  mars_rf_stats <- cbind.data.frame(BentMars_rf_r2_train, BentMars_rf_adj_r2_train, BentMars_rf_MAE_train, BentMars_rf_RMSE_train, BentMars_rf_r2_cv, BentMars_rf_MAE_cv, BentMars_rf_RMSE_cv, BentMars_rf_Q2_test, BentMars_rf_MAE_test, BentMars_rf_RMSE_test, BentMars_rf_scrambled_r2, BentMars_rf_scrambled_mae, BentMars_rf_scrambled_rmse)
  colnames(mars_rf_stats) <- rbind("r2 train", "adj r2 train", "MAE train", "RMSE train", "r2 LOO", "MAE LOO", "RMSE LOO", "q2 test", "MAE test", "RMSE test", "r2 int scramb", "MAE int scramb", "RMSE int scramb")
  
  ## NO FS SVM
  # Fit a cross validated svm model between all scaled descriptors descriptors and Tethering
  BentNoFS_SVMModel <- train(y~.-ID-SMILES, data = noFS_updated_data, method = 'svmRadial', trControl = train.control, scale = FALSE, tuneLength = 15)
  
  # Store fit and error statistics for training set predictions made in the cross validation process
  BentNoFS_svm_r2_cv <- BentNoFS_SVMModel$results$Rsquared[which(BentNoFS_SVMModel$results$C == BentNoFS_SVMModel$bestTune$C)]
  BentNoFS_svm_MAE_cv <- BentNoFS_SVMModel$results$MAE[which(BentNoFS_SVMModel$results$C == BentNoFS_SVMModel$bestTune$C)]
  BentNoFS_svm_RMSE_cv <- BentNoFS_SVMModel$results$RMSE[which(BentNoFS_SVMModel$results$C == BentNoFS_SVMModel$bestTune$C)]
  
  # Store model with optimized parameters and training set predictions made by model
  final_BentNoFS_svmModel <- BentNoFS_SVMModel$finalModel
  BentNoFS_svmFittedVals <- final_BentNoFS_svmModel@fitted
  
  # Store fit and error statistics for training set predictions made by model with optimized parameters
  BentNoFS_svm_r2_train <- (cor(y, BentNoFS_svmFittedVals))^2
  n <- length(y)
  k <- ncol(x_new[,-1])
  BentNoFS_svm_adj_r2_train <- (1 - (((1-BentNoFS_svm_r2_train)*(n-1))/(abs(n-k-1))))
  BentNoFS_svm_RMSE_train <- rmse(y, BentNoFS_svmFittedVals)
  BentNoFS_svm_MAE_train <- mae(y, BentNoFS_svmFittedVals)

  # Use optimized model to make predictions of test set Tethering
  BentNoFS_svmExtPreds <- predict(final_BentNoFS_svmModel, newdata = bentTestSet[,-1])
  # Store fit and error statistics for test set predictions made by model with optimized parameters
  BentNoFS_svm_Q2_test <- (cor(bentTestSet$y, BentNoFS_svmExtPreds))^2
  BentNoFS_svm_RMSE_test <- rmse(bentTestSet$y, BentNoFS_svmExtPreds)
  BentNoFS_svm_MAE_test <- mae(bentTestSet$y, BentNoFS_svmExtPreds)

  # Generate a vector containing the y (Tethering) values in the training set in a randomized order
  BentNoFS_svm_scrambled_y <- sample(y, size = length(y), replace = FALSE)
  BentNoFS_svm_list_of_scrambs <- list()
  
  # Generate and save 45 scrambled list of training set y values
  for (J in 1:45){
    BentNoFS_svm_scrambled_y <- sample(y, size = length(y), replace = FALSE)
    BentNoFS_svm_list_of_scrambs[[J]] <- c(BentNoFS_svm_scrambled_y)
  }
  
  # Initialize empty lists to store fit and error statistics for scrambled training set predictions made by model with optimized parameters 
  BentNoFS_svm_list_of_r2 <- list()
  BentNoFS_svm_list_of_rmse <- list()
  BentNoFS_svm_list_of_mae <- list()
  
  # Use optimized model to make predictions of scrambled training set Tethering
  for (J in 1:length(BentNoFS_svm_list_of_scrambs)){
    BentNoFS_svm_reshuffled <- as.data.frame(cbind(BentNoFS_svm_list_of_scrambs[[J]], subset(noFS_updated_data, select = -c(ID, SMILES, y))))
    BentNoFS_svm_shuffledPreds <- predict(final_BentNoFS_svmModel, newdata = BentNoFS_svm_reshuffled[,-1])
    RMSE_scramb <- rmse(BentNoFS_svm_list_of_scrambs[[J]], BentNoFS_svm_shuffledPreds)
    BentNoFS_svm_list_of_rmse <- c(BentNoFS_svm_list_of_rmse, RMSE_scramb)
    MAE_scramb <- mae(BentNoFS_svm_list_of_scrambs[[J]], BentNoFS_svm_shuffledPreds)
    BentNoFS_svm_list_of_mae <- c(BentNoFS_svm_list_of_mae, MAE_scramb)
    r2 <- (cor(BentNoFS_svm_list_of_scrambs[[J]], BentNoFS_svm_shuffledPreds))^2
    BentNoFS_svm_list_of_r2 <- c(BentNoFS_svm_list_of_r2, r2)
  }
  
  # Store fit and error statistics for scrambled training set predictions made by model with optimized parameters
  BentNoFS_svm_scrambled_r2 <- mean(unlist(BentNoFS_svm_list_of_r2))
  BentNoFS_svm_scrambled_rmse <- mean(unlist(BentNoFS_svm_list_of_rmse))
  BentNoFS_svm_scrambled_mae <- mean(unlist(BentNoFS_svm_list_of_mae))
  # Store all fit and error statistics calculated for no feature selection/svm model
  nf_svm_stats <- cbind.data.frame(BentNoFS_svm_r2_train, BentNoFS_svm_adj_r2_train, BentNoFS_svm_MAE_train, BentNoFS_svm_RMSE_train, BentNoFS_svm_r2_cv, BentNoFS_svm_MAE_cv, BentNoFS_svm_RMSE_cv, BentNoFS_svm_Q2_test, BentNoFS_svm_MAE_test, BentNoFS_svm_RMSE_test, BentNoFS_svm_scrambled_r2, BentNoFS_svm_scrambled_mae, BentNoFS_svm_scrambled_rmse)
  colnames(nf_svm_stats) <- rbind("r2 train", "adj r2 train", "MAE train", "RMSE train", "r2 LOO", "MAE LOO", "RMSE LOO", "q2 test", "MAE test", "RMSE test", "r2 int scramb", "MAE int scramb", "RMSE int scramb")
  
  ## LASSO MIN SVM
  # Check if lasso min feature selection returned any descriptors before proceeding
  if (length(lambda_min_chosen_desc) >= 1){
    # Fit a cross validated svm model between lasso min selected descriptors and Tethering
    BentLassoMin_SVMModel <- train(y~.-ID-SMILES, data = lassoMin_updatedData, method = 'svmRadial', trControl = train.control, scale = FALSE, tuneLength = 15)
    
    # Store fit and error statistics for training set predictions made in the cross validation process
    BentLassoMin_svm_r2_cv <- BentLassoMin_SVMModel$results$Rsquared[which(BentLassoMin_SVMModel$results$C == BentLassoMin_SVMModel$bestTune$C)]
    BentLassoMin_svm_MAE_cv <- BentLassoMin_SVMModel$results$MAE[which(BentLassoMin_SVMModel$results$C == BentLassoMin_SVMModel$bestTune$C)]
    BentLassoMin_svm_RMSE_cv <- BentLassoMin_SVMModel$results$RMSE[which(BentLassoMin_SVMModel$results$C == BentLassoMin_SVMModel$bestTune$C)]
    
    # Store model with optimized parameters and training set predictions made by model
    final_BentLassoMin_svmModel <- BentLassoMin_SVMModel$finalModel
    BentLassoMin_svmFittedVals <- final_BentLassoMin_svmModel@fitted
    
    # Store fit and error statistics for training set predictions made by model with optimized parameters
    BentLassoMin_svm_r2_train <- (cor(y, BentLassoMin_svmFittedVals))^2
    n <- length(y)
    k <- length(lambda_min_chosen_desc)
    BentLassoMin_svm_adj_r2_train <- (1 - (((1-BentLassoMin_svm_r2_train)*(n-1))/(abs(n-k-1))))
    BentLassoMin_svm_RMSE_train <- rmse(y, BentLassoMin_svmFittedVals)
    BentLassoMin_svm_MAE_train <- mae(y, BentLassoMin_svmFittedVals)
  
    # Use optimized model to make predictions of test set Tethering
    BentLassoMin_svmExtPreds <- predict(final_BentLassoMin_svmModel, newdata = test_SubSet_lambdaMin[,-1])
    # Store fit and error statistics for test set predictions made by model with optimized parameters
    BentLassoMin_svm_Q2_test <- (cor(test_SubSet_lambdaMin$y, BentLassoMin_svmExtPreds))^2
    BentLassoMin_svm_RMSE_test <- rmse(test_SubSet_lambdaMin$y, BentLassoMin_svmExtPreds)
    BentLassoMin_svm_MAE_test <- mae(test_SubSet_lambdaMin$y, BentLassoMin_svmExtPreds)
  
    # Generate a vector containing the y (Tethering) values in the training set in a randomized order
    BentLassoMin_svm_scrambled_y <- sample(y, size = length(y), replace = FALSE)
    BentLassoMin_svm_list_of_scrambs <- list()
    
    # Generate and save 45 scrambled list of training set y values
    for (J in 1:45){
      BentLassoMin_svm_scrambled_y <- sample(y, size = length(y), replace = FALSE)
      BentLassoMin_svm_list_of_scrambs[[J]] <- c(BentLassoMin_svm_scrambled_y)
    }
    
    # Initialize empty lists to store fit and error statistics for scrambled training set predictions made by model with optimized parameters 
    BentLassoMin_svm_list_of_r2 <- list()
    BentLassoMin_svm_list_of_rmse <- list()
    BentLassoMin_svm_list_of_mae <- list()
    
    # Use optimized model to make predictions of scrambled training set Tethering
    for (J in 1:length(BentLassoMin_svm_list_of_scrambs)){
      BentLassoMin_svm_reshuffled <- as.data.frame(cbind(BentLassoMin_svm_list_of_scrambs[[J]], subset(lassoMin_updatedData, select = -c(ID, SMILES, y))))
      BentLassoMin_svm_shuffledPreds <- predict(final_BentLassoMin_svmModel, newdata = BentLassoMin_svm_reshuffled[,-1])
      RMSE_scramb <- rmse(BentLassoMin_svm_list_of_scrambs[[J]], BentLassoMin_svm_shuffledPreds)
      BentLassoMin_svm_list_of_rmse <- c(BentLassoMin_svm_list_of_rmse, RMSE_scramb)
      MAE_scramb <- mae(BentLassoMin_svm_list_of_scrambs[[J]], BentLassoMin_svm_shuffledPreds)
      BentLassoMin_svm_list_of_mae <- c(BentLassoMin_svm_list_of_mae, MAE_scramb)
      r2 <- (cor(BentLassoMin_svm_list_of_scrambs[[J]], BentLassoMin_svm_shuffledPreds))^2
      BentLassoMin_svm_list_of_r2 <- c(BentLassoMin_svm_list_of_r2, r2)
    }
    
    # Store fit and error statistics for scrambled training set predictions made by model with optimized parameters
    BentLassoMin_svm_scrambled_r2 <- mean(unlist(BentLassoMin_svm_list_of_r2))
    BentLassoMin_svm_scrambled_rmse <- mean(unlist(BentLassoMin_svm_list_of_rmse))
    BentLassoMin_svm_scrambled_mae <- mean(unlist(BentLassoMin_svm_list_of_mae))
    # Store all fit and error statistics calculated for lasso min/svm model
    lm_svm_stats <- cbind.data.frame(BentLassoMin_svm_r2_train, BentLassoMin_svm_adj_r2_train, BentLassoMin_svm_MAE_train, BentLassoMin_svm_RMSE_train, BentLassoMin_svm_r2_cv, BentLassoMin_svm_MAE_cv, BentLassoMin_svm_RMSE_cv, BentLassoMin_svm_Q2_test, BentLassoMin_svm_MAE_test, BentLassoMin_svm_RMSE_test, BentLassoMin_svm_scrambled_r2, BentLassoMin_svm_scrambled_mae, BentLassoMin_svm_scrambled_rmse)
    colnames(lm_svm_stats) <- rbind("r2 train", "adj r2 train", "MAE train", "RMSE train", "r2 LOO", "MAE LOO", "RMSE LOO", "q2 test", "MAE test", "RMSE test", "r2 int scramb", "MAE int scramb", "RMSE int scramb")
  }  
  
  ## LASSO 1SE SVM
  # Check if lasso 1se feature selection returned any descriptors before proceeding
  if (length(lambda_1se_chosen_desc) >= 1){
    # Fit a cross validated svm model between lasso 1se selected descriptors and Tethering
    BentLasso1se_SVMModel <- train(y~.-ID-SMILES, data = lasso1se_updatedData, method = 'svmRadial', trControl = train.control, scale = FALSE, tuneLength = 10)
    
    # Store fit and error statistics for training set predictions made in the cross validation process
    BentLasso1se_svm_r2_cv <- BentLasso1se_SVMModel$results$Rsquared[which(BentLasso1se_SVMModel$results$C == BentLasso1se_SVMModel$bestTune$C)]
    BentLasso1se_svm_MAE_cv <- BentLasso1se_SVMModel$results$MAE[which(BentLasso1se_SVMModel$results$C == BentLasso1se_SVMModel$bestTune$C)]
    BentLasso1se_svm_RMSE_cv <- BentLasso1se_SVMModel$results$RMSE[which(BentLasso1se_SVMModel$results$C == BentLasso1se_SVMModel$bestTune$C)]
    
    # Store model with optimized parameters and training set predictions made by model
    final_BentLasso1se_svmModel <- BentLasso1se_SVMModel$finalModel
    BentLasso1se_svmFittedVals <- final_BentLasso1se_svmModel@fitted
    
    # Store fit and error statistics for training set predictions made by model with optimized parameters
    BentLasso1se_svm_r2_train <- (cor(y, BentLasso1se_svmFittedVals))^2
    n <- length(y)
    k <- length(lambda_1se_chosen_desc)
    BentLasso1se_svm_adj_r2_train <- (1 - (((1-BentLasso1se_svm_r2_train)*(n-1))/(abs(n-k-1))))
    BentLasso1se_svm_RMSE_train <- rmse(y, BentLasso1se_svmFittedVals)
    BentLasso1se_svm_MAE_train <- mae(y, BentLasso1se_svmFittedVals)

    # Use optimized model to make predictions of test set Tethering
    BentLasso1se_svmExtPreds <- predict(final_BentLasso1se_svmModel, newdata = test_SubSet_lambda1se[,-1])
    # Store fit and error statistics for test set predictions made by model with optimized parameters
    BentLasso1se_svm_Q2_test <- (cor(test_SubSet_lambda1se$y, BentLasso1se_svmExtPreds))^2
    predictiveSD <- sqrt(mean((BentLasso1se_svmExtPreds-test_SubSet_lambda1se$y)^2))
    BentLasso1se_svm_RMSE_test <- rmse(test_SubSet_lambda1se$y, BentLasso1se_svmExtPreds)
    BentLasso1se_svm_MAE_test <- mae(test_SubSet_lambda1se$y, BentLasso1se_svmExtPreds)
    
    # Use optimized model to make predictions of ZINC set Tethering
    if (length(setdiff(lambda_1se_chosen_desc, colnames(bentExtSet))) == 0){
      BentLasso1se_svmExtPreds <- predict(final_BentLasso1se_svmModel, newdata = ext_SubSet_lambda1se)
      extfFileName = paste0("Chunk0s2_BentExtPredsLasso1seSVM", "_", fileNumbers[i], ".csv")
      write.csv(BentLasso1se_svmExtPreds, file = extfFileName)
    }
    
    # Generate a vector containing the y (Tethering) values in the training set in a randomized order
    BentLasso1se_svm_scrambled_y <- sample(y, size = length(y), replace = FALSE)
    BentLasso1se_svm_list_of_scrambs <- list()
    
    # Generate and save 45 scrambled list of training set y values
    for (J in 1:45){
      BentLasso1se_svm_scrambled_y <- sample(y, size = length(y), replace = FALSE)
      BentLasso1se_svm_list_of_scrambs[[J]] <- c(BentLasso1se_svm_scrambled_y)
    }
    
    # Initialize empty lists to store fit and error statistics for scrambled training set predictions made by model with optimized parameters 
    BentLasso1se_svm_list_of_r2 <- list()
    BentLasso1se_svm_list_of_rmse <- list()
    BentLasso1se_svm_list_of_mae <- list()
    
    # Use optimized model to make predictions of scrambled training set Tethering
    for (J in 1:length(BentLasso1se_svm_list_of_scrambs)){
      BentLasso1se_svm_reshuffled <- as.data.frame(cbind(BentLasso1se_svm_list_of_scrambs[[J]], subset(lasso1se_updatedData, select = -c(ID, SMILES, y))))
      BentLasso1se_svm_shuffledPreds <- predict(final_BentLasso1se_svmModel, newdata = BentLasso1se_svm_reshuffled[,-1])
      RMSE_scramb <- rmse(BentLasso1se_svm_list_of_scrambs[[J]], BentLasso1se_svm_shuffledPreds)
      BentLasso1se_svm_list_of_rmse <- c(BentLasso1se_svm_list_of_rmse, RMSE_scramb)
      MAE_scramb <- mae(BentLasso1se_svm_list_of_scrambs[[J]], BentLasso1se_svm_shuffledPreds)
      BentLasso1se_svm_list_of_mae <- c(BentLasso1se_svm_list_of_mae, MAE_scramb)
      r2 <- (cor(BentLasso1se_svm_list_of_scrambs[[J]], BentLasso1se_svm_shuffledPreds))^2
      BentLasso1se_svm_list_of_r2 <- c(BentLasso1se_svm_list_of_r2, r2)
    }
    
    # Store fit and error statistics for scrambled training set predictions made by model with optimized parameters
    BentLasso1se_svm_scrambled_r2 <- mean(unlist(BentLasso1se_svm_list_of_r2))
    BentLasso1se_svm_scrambled_rmse <- mean(unlist(BentLasso1se_svm_list_of_rmse))
    BentLasso1se_svm_scrambled_mae <- mean(unlist(BentLasso1se_svm_list_of_mae))
    # Store all fit and error statistics calculated for lasso 1se/svm model
    l1_svm_stats <- cbind.data.frame(BentLasso1se_svm_r2_train, BentLasso1se_svm_adj_r2_train, BentLasso1se_svm_MAE_train, BentLasso1se_svm_RMSE_train, BentLasso1se_svm_r2_cv, BentLasso1se_svm_MAE_cv, BentLasso1se_svm_RMSE_cv, BentLasso1se_svm_Q2_test, BentLasso1se_svm_MAE_test, BentLasso1se_svm_RMSE_test, BentLasso1se_svm_scrambled_r2, BentLasso1se_svm_scrambled_mae, BentLasso1se_svm_scrambled_rmse)
    colnames(l1_svm_stats) <- rbind("r2 train", "adj r2 train", "MAE train", "RMSE train", "r2 LOO", "MAE LOO", "RMSE LOO", "q2 test", "MAE test", "RMSE test", "r2 int scramb", "MAE int scramb", "RMSE int scramb")
  }
  
  ## MARS SVM
  # Fit a cross validated svm model between mars selected descriptors and Tethering
  BentMars_SVMModel <- train(y~.-ID-SMILES, data = mars_updatedData, method = 'svmRadial', trControl = train.control, scale = FALSE, tuneLength = 10)
  
  # Store fit and error statistics for training set predictions made in the cross validation process
  BentMars_svm_r2_cv <- BentMars_SVMModel$results$Rsquared[which(BentMars_SVMModel$results$C == BentMars_SVMModel$bestTune$C)]
  BentMars_svm_MAE_cv <- BentMars_SVMModel$results$MAE[which(BentMars_SVMModel$results$C == BentMars_SVMModel$bestTune$C)]
  BentMars_svm_RMSE_cv <- BentMars_SVMModel$results$RMSE[which(BentMars_SVMModel$results$C == BentMars_SVMModel$bestTune$C)]
  
  # Store model with optimized parameters and training set predictions made by model
  final_BentMars_svmModel <- BentMars_SVMModel$finalModel
  BentMars_svmFittedVals <- final_BentMars_svmModel@fitted
  
  # Store fit and error statistics for training set predictions made by model with optimized parameters
  BentMars_svm_r2_train <- (cor(y, BentMars_svmFittedVals))^2
  n <- length(y)
  k <- length(mars_chosen_desc)
  BentMars_svm_adj_r2_train <- (1 - (((1-BentMars_svm_r2_train)*(n-1))/(abs(n-k-1))))
  BentMars_svm_RMSE_train <- rmse(y, BentMars_svmFittedVals)
  BentMars_svm_MAE_train <- mae(y, BentMars_svmFittedVals)
  
  # Use optimized model to make predictions of test set Tethering
  BentMars_svmExtPreds <- predict(final_BentMars_svmModel, newdata = test_SubSet_mars[,-1])
  # Store fit and error statistics for test set predictions made by model with optimized parameters
  BentMars_svm_Q2_test <- (cor(test_SubSet_mars$y, BentMars_svmExtPreds))^2
  BentMars_svm_RMSE_test <- rmse(test_SubSet_mars$y, BentMars_svmExtPreds)
  BentMars_svm_MAE_test <- mae(test_SubSet_mars$y, BentMars_svmExtPreds)

  # Use optimized model to make predictions of ZINC set Tethering
  if (length(setdiff(mars_chosen_desc, colnames(bentExtSet))) == 0){
    BentMars_svmExtPreds <- predict(final_BentMars_svmModel, newdata = ext_SubSet_mars)
    extfFileName = paste0("Chunk0s2_BentExtPredsMarsSVM", "_", fileNumbers[i], ".csv")
    write.csv(BentMars_svmExtPreds, file = extfFileName)
  }
  
  # Generate a vector containing the y (Tethering) values in the training set in a randomized order
  BentMars_svm_scrambled_y <- sample(y, size = length(y), replace = FALSE)
  BentMars_svm_list_of_scrambs <- list()
  
  # Generate and save 45 scrambled list of training set y values
  for (J in 1:45){
    BentMars_svm_scrambled_y <- sample(y, size = length(y), replace = FALSE)
    BentMars_svm_list_of_scrambs[[J]] <- c(BentMars_svm_scrambled_y)
  }
  
  # Initialize empty lists to store fit and error statistics for scrambled training set predictions made by model with optimized parameters 
  BentMars_svm_list_of_r2 <- list()
  BentMars_svm_list_of_rmse <- list()
  BentMars_svm_list_of_mae <- list()
  
  # Use optimized model to make predictions of scrambled training set Tethering
  for (J in 1:length(BentMars_svm_list_of_scrambs)){
    BentMars_svm_reshuffled <- as.data.frame(cbind(BentMars_svm_list_of_scrambs[[J]], subset(mars_updatedData, select = -c(ID, SMILES, y))))
    BentMars_svm_shuffledPreds <- predict(final_BentMars_svmModel, newdata = BentMars_svm_reshuffled[,-1])
    RMSE_scramb <- rmse(BentMars_svm_list_of_scrambs[[J]], BentMars_svm_shuffledPreds)
    BentMars_svm_list_of_rmse <- c(BentMars_svm_list_of_rmse, RMSE_scramb)
    MAE_scramb <- mae(BentMars_svm_list_of_scrambs[[J]], BentMars_svm_shuffledPreds)
    BentMars_svm_list_of_mae <- c(BentMars_svm_list_of_mae, MAE_scramb)
    r2 <- (cor(BentMars_svm_list_of_scrambs[[J]], BentMars_svm_shuffledPreds))^2
    BentMars_svm_list_of_r2 <- c(BentMars_svm_list_of_r2, r2)
  }
  
  # Store fit and error statistics for scrambled training set predictions made by model with optimized parameters
  BentMars_svm_scrambled_r2 <- mean(unlist(BentMars_svm_list_of_r2))
  BentMars_svm_scrambled_rmse <- mean(unlist(BentMars_svm_list_of_rmse))
  BentMars_svm_scrambled_mae <- mean(unlist(BentMars_svm_list_of_mae))
  # Store all fit and error statistics calculated for mars/svm model
  mars_svm_stats <- cbind.data.frame(BentMars_svm_r2_train, BentMars_svm_adj_r2_train, BentMars_svm_MAE_train, BentMars_svm_RMSE_train, BentMars_svm_r2_cv, BentMars_svm_MAE_cv, BentMars_svm_RMSE_cv, BentMars_svm_Q2_test, BentMars_svm_MAE_test, BentMars_svm_RMSE_test, BentMars_svm_scrambled_r2, BentMars_svm_scrambled_mae, BentMars_svm_scrambled_rmse)
  colnames(mars_svm_stats) <- rbind("r2 train", "adj r2 train", "MAE train", "RMSE train", "r2 LOO", "MAE LOO", "RMSE LOO", "q2 test", "MAE test", "RMSE test", "r2 int scramb", "MAE int scramb", "RMSE int scramb")
 
  ## Generate dataframe summarizing model statistics for current train/test split
  statistics <- rbind.data.frame(nf_lr_stats, lm_lr_stats, l1_lr_stats, mars_lr_stats, nf_pls_stats, lm_pls_stats, l1_pls_stats, mars_pls_stats, nf_rf_stats, lm_rf_stats, l1_rf_stats, mars_rf_stats, nf_svm_stats, lm_svm_stats, l1_svm_stats, mars_svm_stats)
  new_stats <- t(statistics)
  colnames(new_stats) <- rbind("No FS Linear Regression", "Lasso Min Linear Regression",  "Lasso 1se Linear Regression", "Mars Linear Regression", "No FS PLS", "Lasso Min PLS", "Lasso 1se PLS", "Mars PLS", "No FS RF", "Lasso Min RF", "Lasso 1se RF", "Mars RF", "No FS SVM",  "Lasso Min SVM", "Lasso 1se SVM", "Mars SVM")
}

