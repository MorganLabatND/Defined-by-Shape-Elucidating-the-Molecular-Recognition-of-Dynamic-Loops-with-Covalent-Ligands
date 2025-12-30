setwd('/lsi/groups/amapplab/chloemarkey')

# INSTALL AND LOAD ALL NECESSARY PACKAGES AND LIBRARIES -----------------------------------------------------------------
install.packages(c("stringi", 'glmnet', 'ModelMetrics', 'pls', 'randomForest', 'kernlab', 'earth', 'prospectr'), repos="http://cran.r-project.org", lib="/lsi/groups/amapplab/chloemarkey/Anaconda_libs/")
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
library(stringi, lib="/lsi/groups/amapplab/chloemarkey/Anaconda_libs/")
library(caret, lib="/lsi/groups/amapplab/chloemarkey/Anaconda_libs/")
library(glmnet, lib="/lsi/groups/amapplab/chloemarkey/Anaconda_libs/")
library(ModelMetrics, lib="/lsi/groups/amapplab/chloemarkey/Anaconda_libs/")
library(pls, lib="/lsi/groups/amapplab/chloemarkey/Anaconda_libs/")
library(randomForest, lib="/lsi/groups/amapplab/chloemarkey/Anaconda_libs/")
library(kernlab, lib="/lsi/groups/amapplab/chloemarkey/Anaconda_libs/")
library(earth, lib="/lsi/groups/amapplab/chloemarkey/Anaconda_libs/")
library(prospectr, lib="/lsi/groups/amapplab/chloemarkey/Anaconda_libs/")

# DATA SET 22 ---------------------------------------------------------------------------
# Read in training set
trainingSet22 <- read.csv('/lsi/groups/amapplab/chloemarkey/trainSet22.csv')
trainingSet22 <- subset(trainingSet22, select = -c(Group))
colnames(trainingSet22)[3] <- "y"

# Read in test set
testingSet <- read.csv('/lsi/groups/amapplab/chloemarkey/testSet22.csv')
testingSet <- subset(testingSet, select = -c(Group))
colnames(testingSet)[3] <- "y"

# Read in ZINC set
extData <- read.csv('/lsi/groups/amapplab/chloemarkey/chunk0_secondHalf_0331boltzmanned.csv')

## IDENTIFY AND REMOVE INVARIANT DESCRIPTORS
set.seed(1234)

# Create a new training data set with invariant descriptors removed
all_desc <- colnames(trainingSet22)
invar_desc_removed_trainSet <- trainingSet22[vapply(trainingSet22, function(z) length(unique(z))>1, logical(1L))]
# Find column indices corresponding to invariant descriptors
reduced_desc <- colnames(invar_desc_removed_trainSet)
invariant_desc <- setdiff(all_desc, reduced_desc)
invariant_inds <- integer()

for (a in 1:length(invariant_desc)){
  index <- which(colnames(testingSet) == invariant_desc[a])
  invariant_inds <- c(invariant_inds, index)
}

# Save remaining descriptors as a model matrix
x_train22 = model.matrix(y ~.-ID-SMILES, data = invar_desc_removed_trainSet)
x_train22 = x_train22[,-1]

# Save molecule identifier information and Tethering (y) values
ID = invar_desc_removed_trainSet$ID
SMILES = invar_desc_removed_trainSet$SMILES
y = as.numeric(invar_desc_removed_trainSet$y)

## IDENTIFY AND REMOVE PERFECTLY CORRELATED DESCRIPTORS
# Generate matrix containing correlations between each descriptors and Tethering
bent_desc_tethering_cormat = cor(x_train22, y, method = "kendall")

# Generate matrix containing correlations between all descriptors and each other
bent_descriptor_cormat = cor(x_train22, method = "kendall")
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
x_train22 = x_train22[,-c(new_ind_to_remove)]
bent_desc_tethering_cormat = cor(x_train22, y, method = "kendall")
bent_descriptor_cormat = cor(x_train22, method = "kendall")
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
remaining_descriptors_80 <- x_train22[,-c(new_ind_to_remove_80)]
updated_data22 <- trainingSet22[,c("ID", "SMILES", "y", colnames(remaining_descriptors_80))]

## SCALE DESCRIPTORS
# Generate model matrix of training set containing remaining descriptors
x_new22 = model.matrix(y ~ .-ID-SMILES, data = updated_data22)

# Scale all descriptor values in training set model matrix
x_means22 <- colMeans(x_new22)
x_sd22 <- apply(x_new22, 2, sd)
ext_desc = model.matrix(Ind ~ .-SMILES, data = extData)
cols_to_keep <- which(colnames(ext_desc) %in% colnames(x_new22))
cols_to_keep1 <- which(colnames(x_new22) %in% colnames(ext_desc))
x_new_ext22 <- ext_desc[,c(cols_to_keep)]
x_train_dat_forExt22 <- x_new22[,c(cols_to_keep1)]
x_train_means_forExt22 <- colMeans(x_train_dat_forExt22)
x_train_sd_forExt22 <- apply(x_train_dat_forExt22, 2, sd)

for (K in 1:(ncol(x_new22))){
  for (j in 1:(nrow(x_new22))){
    scaled_val <- (x_new22[j,K]-x_means22[K])/(x_sd22[K])
    x_new22[j,K] <- scaled_val
  }
}
x_new22[,1] = 1

# Scale all descriptor values in ZINC set model matrix
for (k in 1:(ncol(x_new_ext22))){
  for (j in 1:(nrow(x_new_ext22))){
    scaled_val <- ((x_new_ext22[j,k]-x_train_means_forExt22[k]))/(x_train_sd_forExt22[k])
    x_new_ext22[j,k] <- scaled_val
  }
}
x_new_ext22 = x_new_ext22[,-1]
bentExtSet22 <- x_new_ext22
bentExtSet22[is.na(bentExtSet22)] = 0

# Generate a model matrix containing only training set descriptors 
trainingDesc22 <- subset(x_new22, select = c(vsurf_DD23, vsurf_DW23, vsurf_IW7, a_nS, chiral, GCUT_SLOGP_0, GCUT_SMR_0, PEOE_VSAm6, SlogP_VSA5))
# Calculate hats of training set descriptors
hats_train <- hat(trainingDesc22)

# Generate a model matrix containing only ZINC set descriptors 
extDesc22 <- subset(bentExtSet22, select = c(vsurf_DD23, vsurf_DW23, vsurf_IW7, a_nS, chiral, GCUT_SLOGP_0, GCUT_SMR_0, PEOE_VSAm6, SlogP_VSA5))

# Calculate hats for ZINC molecules by comparing descriptors to training set descriptors
ext22hats <- list()
for (i in 1:nrow(extDesc22)){
  newDat <- rbind(trainingDesc22, extDesc22[i,])
  hats <- hat(newDat)
  element <- tail(hats, n=1)
  ext22hats <- c(ext22hats, element)
}

# Generate output file summarizing hats for ZINC molecules/current data split
new_ext22hats <- t(ext22hats)
write.csv(new_ext22hats, file = "chunk0secondHalf_BentExtHats_22.csv")

# DATA SET 4 ---------------------------------------------------------------------------
# Read in training set
trainingSet4 <- read.csv('/lsi/groups/amapplab/chloemarkey/trainSet4.csv')
trainingSet4 <- subset(trainingSet4, select = -c(Group))
colnames(trainingSet4)[3] <- "y"

# Read in test set
testingSet <- read.csv('/lsi/groups/amapplab/chloemarkey/testSet4.csv')
testingSet <- subset(testingSet, select = -c(Group))
colnames(testingSet)[3] <- "y"

# Read in ZINC set
extData <- read.csv('/lsi/groups/amapplab/chloemarkey/chunk0_secondHalf_0331boltzmanned.csv')

## IDENTIFY AND REMOVE INVARIANT DESCRIPTORS
set.seed(1234)

# Create a new training data set with invariant descriptors removed
all_desc <- colnames(trainingSet4)
invar_desc_removed_trainSet <- trainingSet4[vapply(trainingSet4, function(z) length(unique(z))>1, logical(1L))]

# Find column indices corresponding to invariant descriptors
reduced_desc <- colnames(invar_desc_removed_trainSet)
invariant_desc <- setdiff(all_desc, reduced_desc)
invariant_inds <- integer()

for (a in 1:length(invariant_desc)){
  index <- which(colnames(testingSet) == invariant_desc[a])
  invariant_inds <- c(invariant_inds, index)
}

# Save remaining descriptors as a model matrix
x_train4 = model.matrix(y ~.-ID-SMILES, data = invar_desc_removed_trainSet)
x_train4 = x_train4[,-1]

# Save molecule identifier information and Tethering (y) values
ID = invar_desc_removed_trainSet$ID
SMILES = invar_desc_removed_trainSet$SMILES
y = as.numeric(invar_desc_removed_trainSet$y)

## IDENTIFY AND REMOVE PERFECTLY CORRELATED DESCRIPTORS
# Generate matrix containing correlations between each descriptors and Tethering
bent_desc_tethering_cormat = cor(x_train4, y, method = "kendall")

# Generate matrix containing correlations between all descriptors and each other
bent_descriptor_cormat = cor(x_train4, method = "kendall")
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
      #print(paste(colnames(bent_descriptor_cormat)[k], "-" , rownames(bent_descriptor_cormat)[j], ": ", bent_descriptor_cormat[j,k]))
      d = rbind(d, data.frame(colnames(bent_descriptor_cormat)[k], bent_desc_tethering_cormat[k], rownames(bent_descriptor_cormat)[j], bent_desc_tethering_cormat[j]))
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
x_train4 = x_train4[,-c(new_ind_to_remove)]
bent_desc_tethering_cormat = cor(x_train4, y, method = "kendall")
bent_descriptor_cormat = cor(x_train4, method = "kendall")
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
remaining_descriptors_80 <- x_train4[,-c(new_ind_to_remove_80)]
updated_data4 <- trainingSet4[,c("ID", "SMILES", "y", colnames(remaining_descriptors_80))]

## SCALE DESCRIPTORS
# Generate model matrix of training set containing remaining descriptors
x_new4 = model.matrix(y ~ .-ID-SMILES, data = updated_data4)

# Scale all descriptor values in training set model matrix
x_means4 <- colMeans(x_new4)
x_sd4 <- apply(x_new4, 2, sd)
ext_desc = model.matrix(Ind ~ .-SMILES, data = extData)
cols_to_keep <- which(colnames(ext_desc) %in% colnames(x_new4))
cols_to_keep1 <- which(colnames(x_new4) %in% colnames(ext_desc))
x_new_ext4 <- ext_desc[,c(cols_to_keep)]
x_train_dat_forExt4 <- x_new4[,c(cols_to_keep1)]
x_train_means_forExt4 <- colMeans(x_train_dat_forExt4)
x_train_sd_forExt4 <- apply(x_train_dat_forExt4, 2, sd)

for (K in 1:(ncol(x_new4))){
  for (j in 1:(nrow(x_new4))){
    scaled_val <- (x_new4[j,K]-x_means4[K])/(x_sd4[K])
    x_new4[j,K] <- scaled_val
  }
}
x_new4[,1] = 1

# Scale all descriptor values in ZINC set model matrix
for (k in 1:(ncol(x_new_ext4))){
  for (j in 1:(nrow(x_new_ext4))){
    scaled_val <- ((x_new_ext4[j,k]-x_train_means_forExt4[k]))/(x_train_sd_forExt4[k])
    x_new_ext4[j,k] <- scaled_val
  }
}
x_new_ext4 = x_new_ext4[,-1]
bentExtSet4 <- x_new_ext4
bentExtSet4[is.na(bentExtSet4)] = 0

# Generate a model matrix containing only training set descriptors 
trainingDesc4 <- subset(x_new4, select = c(vsurf_IW5, vsurf_IW6, vsurf_IW7, logP.o.w.))
# Calculate hats of training set descriptors
hats_train <- hat(trainingDesc4)

# Generate a model matrix containing only ZINC set descriptors 
extDesc4 <- subset(bentExtSet4, select = c(vsurf_IW5, vsurf_IW6, vsurf_IW7, logP.o.w.))

# Calculate hats for ZINC molecules by comparing descriptors to training set descriptors
ext4hats <- list()
for (i in 1:nrow(extDesc4)){
  newDat <- rbind(trainingDesc4, extDesc4[i,])
  hats <- hat(newDat)
  element <- tail(hats, n=1)
  ext4hats <- c(ext4hats, element)
}

# Generate output file summarizing hats for ZINC molecules/current data split
new_ext4hats <- t(ext4hats)
write.csv(new_ext4hats, file = "chunk0secondHalf_BentExtHats_4.csv")

# DATA SET 59 ---------------------------------------------------------------------------
# Read in training set
trainingSet59 <- read.csv('/lsi/groups/amapplab/chloemarkey/trainSet59.csv')
trainingSet59 <- subset(trainingSet59, select = -c(Group))
colnames(trainingSet59)[3] <- "y"

# Read in test set
testingSet <- read.csv('/lsi/groups/amapplab/chloemarkey/testSet59.csv')
testingSet <- subset(testingSet, select = -c(Group))
colnames(testingSet)[3] <- "y"

# Read in ZINC set
extData <- read.csv('/lsi/groups/amapplab/chloemarkey/chunk0_secondHalf_0331boltzmanned.csv')

## IDENTIFY AND REMOVE INVARIANT DESCRIPTORS
set.seed(1234)

# Create a new training data set with invariant descriptors removed
all_desc <- colnames(trainingSet59)
invar_desc_removed_trainSet <- trainingSet59[vapply(trainingSet59, function(z) length(unique(z))>1, logical(1L))]
# Find column indices corresponding to invariant descriptors
reduced_desc <- colnames(invar_desc_removed_trainSet)
invariant_desc <- setdiff(all_desc, reduced_desc)
invariant_inds <- integer()

for (a in 1:length(invariant_desc)){
  index <- which(colnames(testingSet) == invariant_desc[a])
  invariant_inds <- c(invariant_inds, index)
}

# Save remaining descriptors as a model matrix
x_train59 = model.matrix(y ~.-ID-SMILES, data = invar_desc_removed_trainSet)
x_train59 = x_train59[,-1]

# Save molecule identifier information and Tethering (y) values
ID = invar_desc_removed_trainSet$ID
SMILES = invar_desc_removed_trainSet$SMILES
y = as.numeric(invar_desc_removed_trainSet$y)

## IDENTIFY AND REMOVE PERFECTLY CORRELATED DESCRIPTORS
# Generate matrix containing correlations between each descriptors and Tethering
bent_desc_tethering_cormat = cor(x_train59, y, method = "kendall")

# Generate matrix containing correlations between all descriptors and each other
bent_descriptor_cormat = cor(x_train59, method = "kendall")
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
      #print(paste(colnames(bent_descriptor_cormat)[k], "-" , rownames(bent_descriptor_cormat)[j], ": ", bent_descriptor_cormat[j,k]))
      d = rbind(d, data.frame(colnames(bent_descriptor_cormat)[k], bent_desc_tethering_cormat[k], rownames(bent_descriptor_cormat)[j], bent_desc_tethering_cormat[j]))
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
x_train59 = x_train59[,-c(new_ind_to_remove)]
bent_desc_tethering_cormat = cor(x_train59, y, method = "kendall")
bent_descriptor_cormat = cor(x_train59, method = "kendall")
bent_descriptor_cormat = round(bent_descriptor_cormat, 3)
bent_descriptor_cormat[upper.tri(bent_descriptor_cormat)] <- 0
diag(bent_descriptor_cormat) <- 0
descriptor_list <- colnames(bent_descriptor_cormat)

# Generate a vector of descriptors to remove to reduce pair-wise correlations
built_in_inds_to_remove_80 <- findCorrelation(bent_descriptor_cormat, cutoff = 0.800)

# Initialize empty list to store indices of descriptors to be removed
ind_to_remove_80 <- list()

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
remaining_descriptors_80 <- x_train59[,-c(new_ind_to_remove_80)]
updated_data59 <- trainingSet59[,c("ID", "SMILES", "y", colnames(remaining_descriptors_80))]

## SCALE DESCRIPTORS
# Generate model matrix of training set containing remaining descriptors
x_new59 = model.matrix(y ~ .-ID-SMILES, data = updated_data59)

# Scale all descriptor values in training set model matrix
x_means59 <- colMeans(x_new59)
x_sd59 <- apply(x_new59, 2, sd)
ext_desc = model.matrix(Ind ~ .-SMILES, data = extData)
cols_to_keep <- which(colnames(ext_desc) %in% colnames(x_new59))
cols_to_keep1 <- which(colnames(x_new59) %in% colnames(ext_desc))
x_new_ext59 <- ext_desc[,c(cols_to_keep)]
x_train_dat_forExt59 <- x_new59[,c(cols_to_keep1)]
x_train_means_forExt59 <- colMeans(x_train_dat_forExt59)
x_train_sd_forExt59 <- apply(x_train_dat_forExt59, 2, sd)

for (K in 1:(ncol(x_new59))){
  for (j in 1:(nrow(x_new59))){
    scaled_val <- (x_new59[j,K]-x_means59[K])/(x_sd59[K])
    x_new59[j,K] <- scaled_val
  }
}
x_new59[,1] = 1

# Scale all descriptor values in ZINC set model matrix
for (k in 1:(ncol(x_new_ext59))){
  for (j in 1:(nrow(x_new_ext59))){
    scaled_val <- ((x_new_ext59[j,k]-x_train_means_forExt59[k]))/(x_train_sd_forExt59[k])
    x_new_ext59[j,k] <- scaled_val
  }
}
x_new_ext59 = x_new_ext59[,-1]
bentExtSet59 <- x_new_ext59
bentExtSet59[is.na(bentExtSet59)] = 0

# Generate a model matrix containing only training set descriptors 
trainingDesc59 <- subset(x_new59, select = c(GCUT_PEOE_1, GCUT_SMR_0, PEOE_VSA_POS))
# Calculate hats of training set descriptors
hats_train <- hat(trainingDesc59)

# Generate a model matrix containing only ZINC set descriptors 
extDesc59 <- subset(bentExtSet59, select = c(GCUT_PEOE_1, GCUT_SMR_0, PEOE_VSA_POS))

# Calculate hats for ZINC molecules by comparing descriptors to training set descriptors
ext59hats <- list()
for (i in 1:nrow(extDesc59)){
  newDat <- rbind(trainingDesc59, extDesc59[i,])
  hats <- hat(newDat)
  element <- tail(hats, n=1)
  ext59hats <- c(ext59hats, element)
}

# Generate output file summarizing hats for ZINC molecules/current data split
new_ext59hats <- t(ext59hats)
write.csv(new_ext59hats, file = "chunk0secondHalf_BentExtHats_59.csv")

# DATA SET 9 ---------------------------------------------------------------------------
# Read in training set
trainingSet9 <- read.csv('/lsi/groups/amapplab/chloemarkey/trainSet9.csv')
trainingSet9 <- subset(trainingSet9, select = -c(Group))
colnames(trainingSet9)[3] <- "y"

# Read in test set
testingSet <- read.csv('/lsi/groups/amapplab/chloemarkey/testSet9.csv')
testingSet <- subset(testingSet, select = -c(Group))
colnames(testingSet)[3] <- "y"

# Read in ZINC set
extData <- read.csv('/lsi/groups/amapplab/chloemarkey/chunk0_secondHalf_0331boltzmanned.csv')

## IDENTIFY AND REMOVE INVARIANT DESCRIPTORS
set.seed(1234)

# Create a new training data set with invariant descriptors removed
all_desc <- colnames(trainingSet9)
invar_desc_removed_trainSet <- trainingSet9[vapply(trainingSet9, function(z) length(unique(z))>1, logical(1L))]
# Find column indices corresponding to invariant descriptors
reduced_desc <- colnames(invar_desc_removed_trainSet)
invariant_desc <- setdiff(all_desc, reduced_desc)
invariant_inds <- integer()

for (a in 1:length(invariant_desc)){
  index <- which(colnames(testingSet) == invariant_desc[a])
  invariant_inds <- c(invariant_inds, index)
}

# Save remaining descriptors as a model matrix
x_train9 = model.matrix(y ~.-ID-SMILES, data = invar_desc_removed_trainSet)
x_train9 = x_train9[,-1]

# Save molecule identifier information and Tethering (y) values
ID = invar_desc_removed_trainSet$ID
SMILES = invar_desc_removed_trainSet$SMILES
y = as.numeric(invar_desc_removed_trainSet$y)

## IDENTIFY AND REMOVE PERFECTLY CORRELATED DESCRIPTORS
# Generate matrix containing correlations between each descriptors and Tethering
bent_desc_tethering_cormat = cor(x_train9, y, method = "kendall")

# Generate matrix containing correlations between all descriptors and each other
bent_descriptor_cormat = cor(x_train9, method = "kendall")
bent_descriptor_cormat = round(bent_descriptor_cormat, 3)
# Set diagonal and upper triangle of matrix to 0 to prevent repeat analysis
bent_descriptor_cormat[upper.tri(bent_descriptor_cormat)] <- 0
diag(bent_descriptor_cormat) <- 0
descriptor_list <- colnames(bent_descriptor_cormat)

# Initialize empty list to store indices of descriptors to be removed
ind_to_remove <- list()

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
x_train9 = x_train9[,-c(new_ind_to_remove)]
bent_desc_tethering_cormat = cor(x_train9, y, method = "kendall")
bent_descriptor_cormat = cor(x_train9, method = "kendall")
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
remaining_descriptors_80 <- x_train9[,-c(new_ind_to_remove_80)]
updated_data9 <- trainingSet9[,c("ID", "SMILES", "y", colnames(remaining_descriptors_80))]

## SCALE DESCRIPTORS
# Generate model matrix of training set containing remaining descriptors
x_new9 = model.matrix(y ~ .-ID-SMILES, data = updated_data9)

# Scale all descriptor values in training set model matrix
x_means9 <- colMeans(x_new9)
x_sd9 <- apply(x_new9, 2, sd)
ext_desc = model.matrix(Ind ~ .-SMILES, data = extData)
cols_to_keep <- which(colnames(ext_desc) %in% colnames(x_new9))
cols_to_keep1 <- which(colnames(x_new9) %in% colnames(ext_desc))
x_new_ext9 <- ext_desc[,c(cols_to_keep)]
x_train_dat_forExt9 <- x_new9[,c(cols_to_keep1)]
x_train_means_forExt9 <- colMeans(x_train_dat_forExt9)
x_train_sd_forExt9 <- apply(x_train_dat_forExt9, 2, sd)

for (K in 1:(ncol(x_new9))){
  for (j in 1:(nrow(x_new9))){
    scaled_val <- (x_new9[j,K]-x_means9[K])/(x_sd9[K])
    x_new9[j,K] <- scaled_val
  }
}
x_new9[,1] = 1

# Scale all descriptor values in ZINC set model matrix
for (k in 1:(ncol(x_new_ext9))){
  for (j in 1:(nrow(x_new_ext9))){
    scaled_val <- ((x_new_ext9[j,k]-x_train_means_forExt9[k]))/(x_train_sd_forExt9[k])
    x_new_ext9[j,k] <- scaled_val
  }
}
x_new_ext9 = x_new_ext9[,-1]
bentExtSet9 <- x_new_ext9
bentExtSet9[is.na(bentExtSet9)] = 0

# Generate a model matrix containing only training set descriptors 
trainingDesc9 <- subset(x_new9, select = c(vsurf_IW4, BCUT_SLOGP_0, BCUT_SMR_1, GCUT_PEOE_0, logP.o.w.))
# Calculate hats of training set descriptors
hats_train <- hat(trainingDesc9)

# Generate a model matrix containing only ZINC set descriptors 
extDesc9 <- subset(bentExtSet9, select = c(vsurf_IW4, BCUT_SLOGP_0, BCUT_SMR_1, GCUT_PEOE_0, logP.o.w.))

# Calculate hats for ZINC molecules by comparing descriptors to training set descriptors
ext9hats <- list()
for (i in 1:nrow(extDesc9)){
  newDat <- rbind(trainingDesc9, extDesc9[i,])
  hats <- hat(newDat)
  element <- tail(hats, n=1)
  ext9hats <- c(ext9hats, element)
}

# Generate output file summarizing hats for ZINC molecules/current data split
new_ext9hats <- t(ext9hats)
write.csv(new_ext9hats, file = "chunk0secondHalf_BentExtHats_9.csv")

# DATA SET 90 ---------------------------------------------------------------------------
# Read in training set
trainingSet90 <- read.csv('/lsi/groups/amapplab/chloemarkey/trainSet90.csv')
trainingSet90 <- subset(trainingSet90, select = -c(Group))
colnames(trainingSet90)[3] <- "y"

# Read in test set
testingSet <- read.csv('/lsi/groups/amapplab/chloemarkey/testSet90.csv')
testingSet <- subset(testingSet, select = -c(Group))
colnames(testingSet)[3] <- "y"

# Read in ZINC set
extData <- read.csv('/lsi/groups/amapplab/chloemarkey/chunk0_secondHalf_0331boltzmanned.csv')

## IDENTIFY AND REMOVE INVARIANT DESCRIPTORS
set.seed(1234)

# Create a new training data set with invariant descriptors removed
all_desc <- colnames(trainingSet90)
invar_desc_removed_trainSet <- trainingSet90[vapply(trainingSet90, function(z) length(unique(z))>1, logical(1L))]
# Find column indices corresponding to invariant descriptors
reduced_desc <- colnames(invar_desc_removed_trainSet)
invariant_desc <- setdiff(all_desc, reduced_desc)
invariant_inds <- integer()

for (a in 1:length(invariant_desc)){
  index <- which(colnames(testingSet) == invariant_desc[a])
  invariant_inds <- c(invariant_inds, index)
}

# Save remaining descriptors as a model matrix
x_train90 = model.matrix(y ~.-ID-SMILES, data = invar_desc_removed_trainSet)
x_train90 = x_train90[,-1]

# Save molecule identifier information and Tethering (y) values
ID = invar_desc_removed_trainSet$ID
SMILES = invar_desc_removed_trainSet$SMILES
y = as.numeric(invar_desc_removed_trainSet$y)

## IDENTIFY AND REMOVE PERFECTLY CORRELATED DESCRIPTORS
# Generate matrix containing correlations between each descriptors and Tethering
bent_desc_tethering_cormat = cor(x_train90, y, method = "kendall")

# Generate matrix containing correlations between all descriptors and each other
bent_descriptor_cormat = cor(x_train90, method = "kendall")
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
x_train90 = x_train90[,-c(new_ind_to_remove)]
bent_desc_tethering_cormat = cor(x_train90, y, method = "kendall")
bent_descriptor_cormat = cor(x_train90, method = "kendall")
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
remaining_descriptors_80 <- x_train90[,-c(new_ind_to_remove_80)]
updated_data90 <- trainingSet90[,c("ID", "SMILES", "y", colnames(remaining_descriptors_80))]

## SCALE DESCRIPTORS
# Generate model matrix of training set containing remaining descriptors
x_new90 = model.matrix(y ~ .-ID-SMILES, data = updated_data90)

# Scale all descriptor values in training set model matrix
x_means90 <- colMeans(x_new90)
x_sd90 <- apply(x_new90, 2, sd)
ext_desc = model.matrix(Ind ~ .-SMILES, data = extData)
cols_to_keep <- which(colnames(ext_desc) %in% colnames(x_new90))
cols_to_keep1 <- which(colnames(x_new90) %in% colnames(ext_desc))
x_new_ext90 <- ext_desc[,c(cols_to_keep)]
x_train_dat_forExt90 <- x_new90[,c(cols_to_keep1)]
x_train_means_forExt90 <- colMeans(x_train_dat_forExt90)
x_train_sd_forExt90 <- apply(x_train_dat_forExt90, 2, sd)

for (K in 1:(ncol(x_new90))){
  for (j in 1:(nrow(x_new90))){
    scaled_val <- (x_new90[j,K]-x_means90[K])/(x_sd90[K])
    x_new90[j,K] <- scaled_val
  }
}
x_new90[,1] = 1

# Scale all descriptor values in ZINC set model matrix
for (k in 1:(ncol(x_new_ext90))){
  for (j in 1:(nrow(x_new_ext90))){
    scaled_val <- ((x_new_ext90[j,k]-x_train_means_forExt90[k]))/(x_train_sd_forExt90[k])
    x_new_ext90[j,k] <- scaled_val
  }
}
x_new_ext90 = x_new_ext90[,-1]
bentExtSet90 <- x_new_ext90
bentExtSet90[is.na(bentExtSet90)] = 0

# Generate a model matrix containing only training set descriptors 
trainingDesc90 <- subset(x_new90, select = c(vsurf_DD23, vsurf_DW23, a_nS, chiral, GCUT_SLOGP_0, GCUT_SMR_0, SlogP_VSA0))
# Calculate hats of training set descriptors
hats_train <- hat(trainingDesc90)

# Generate a model matrix containing only ZINC set descriptors 
extDesc90 <- subset(bentExtSet90, select = c(vsurf_DD23, vsurf_DW23, a_nS, chiral, GCUT_SLOGP_0, GCUT_SMR_0, SlogP_VSA0))

# Calculate hats for ZINC molecules by comparing descriptors to training set descriptors
ext90hats <- list()
for (i in 1:nrow(extDesc90)){
  newDat <- rbind(trainingDesc90, extDesc90[i,])
  hats <- hat(newDat)
  element <- tail(hats, n=1)
  ext90hats <- c(ext90hats, element)
}

# Generate output file summarizing hats for ZINC molecules/current data split
new_ext90hats <- t(ext90hats)
write.csv(new_ext90hats, file = "chunk0secondHalf_BentExtHats_90.csv")



