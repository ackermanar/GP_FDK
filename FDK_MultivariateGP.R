# Description: multivariate genomic prediction
# rrBLUP package, caret for stratified grouping of folds
# by AJ Ackerman
# December 2021

###############################################################################

#Load packages
library(rrBLUP)
library(caret)
library(dplyr)
library(data.table)

#Load in data
pheno <- read.csv("FDK_ALLDATA_hmpEntries_NoLA.csv", header = TRUE)
kin <- as.matrix(read.table("LD08_numeric.kin.txt", header = TRUE, sep = "\t", row.names = 1, as.is=TRUE))

#Create empty dataframe for saving results for all phenotypes
GS_results <- data.frame(fold=1:5)

#Set variable to use phenotypes in columns 4:n (FDK traits) as covarites to predict phenotype in column 3 (DON)
for (i in 4:ncol(pheno)) {
  GS <- data.frame() #Create empty dataframe to save results for each individual phenotype

  #create pheno file for analysis, and append variable for additional covariate traits
  y <- pheno %>%
    select(taxa, Env, DON_ppm, i) %>%
    rename(trait = 4)
  y <- na.omit(y) #Omit rows with missing phenotypes

  #Use the training.kfold function in caret to create 5 stratified partitions for test & train
  training.kfold <- groupKFold(y$taxa, k = 5)

    #Mask the testing set within the trianing set with NA and create testing set with all available data
    for (j in training.kfold) {
      DON <- y$DON_ppm %>%
        replace(-j, "NA") #index every position not within current training fold (test set) and mask DON value with NA
      train <- data.frame(y$taxa, y$Env, DON, y$trait) #create dataframe with taxa, env, DON, and trait covariate
      train <- rename(train, trait = y.trait)
      train[,"DON"] <- as.numeric(train$DON) #Convert DON to numeric as NAs are included - will cause NA coercion error (Nonfatal)

      #Create test set by subtracting training set from y (no masking of DON values)
      test <- y %>%
        slice(-j)

      #Run kin.blup for prediction of genomic values
      trainset <- kin.blup(train,geno="y.taxa",pheno="DON",K=kin,fixed="y.Env", covariate = "trait", n.core=12)
      testset <- kin.blup(test,geno="taxa",pheno="DON_ppm",K=kin,fixed="Env", covariate = "trait", n.core=12)

      #Data wrangle results for analysis - I've tried nesting this with pipes but it always seems to cause errors
      trainset_res <- as.data.frame(trainset$pred) #convert results to data frame
      trainset_res <- setDT(trainset_res, keep.rownames = TRUE) #convert observations (taxa) to column one
      trainset_res <- rename(trainset_res, taxa = rn) #rename column one "taxa"
      trainset_res <- filter(trainset_res, trainset_res$taxa %in% test$taxa) #filter training set to only show predicted values (masked in training set with NA)

      #Repeat process for trestset
      testset_res <- as.data.frame(testset$pred)
      testset_res <- setDT(testset_res, keep.rownames = TRUE)
      testset_res <- rename(testset_res, taxa = rn)
      testset_res <- filter(testset_res, testset_res$taxa %in% test$taxa) #filter testset for only corresponding values to masked values in  training set (predicted values)

      #Join train set and test set for organized dataframe
      df <- left_join(trainset_res, testset_res, by ="taxa")

      #Find correlation of observed and predicted values using pearson correlations
      cor <- as.data.frame(cor(df[[2]], df[[3]], method = "pearson"))

      #Add correlation value to dataframe
      GS <- bind_rows(GS, cor)
    }

    #Coerce all GS dataframes into single df with all results for each phenotype
    newname <- paste0("Covariate_", i)
    GS <- rename(GS, !! newname := 1)
    GS_results <- bind_cols(GS_results, GS)
  }
