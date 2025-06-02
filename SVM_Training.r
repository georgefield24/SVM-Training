#load packages
library(e1071) # SVM, KNN
library(caTools) 
library(caret) #caret termed as Classification and Regression Training
library(randomForest) 
library(xgboost)
library(data.table)
library(ggplot2)
library(dplyr)
library(class) 
library(doMC)
library(parallel)
library(foreach)
library(kernlab)
library(LiblineaR)
library(ggalluvial)
library(doParallel)
library(caretEnsemble)
library(rlist)
library(plotly)
library(readxl)
library(ggpubr)
library(devtools)



#import metadata and DEqMS proteins
ImpGenesForSVM=read.table("/Users/george.field/Desktop/Data Analysis/Project 1/timsTOF_KTSP/Current1/GenesForMLNewComplete.txt")
MetaData=read_xls("/Users/george.field/Desktop/Data Analysis/Project 1/ImportData/MetaData.xls")
MetaData=MetaData[(MetaData$number %in% c(colnames(ImpGenesForSVM))), ]


#SVM of our proteins
ImpGenesForSVM=as.data.frame(t(ImpGenesForSVM))
ImpGenesForSVM$Subtype=MetaData$Proteomics
data=ImpGenesForSVM


#Extract subtype and column median centre data
Subtype=data$Subtype
log.data=log2(data[1:ncol(data)-1])
medians = apply(log.data, 2, median, na.rm = TRUE)
df_modified = sweep(log.data, 2, medians, FUN = "-")


# Combine with the last column which is unchanged
df_modified = cbind(df_modified, Subtype)
data=df_modified


# data matrix with grouping varaible 
data$Subtype = factor(data$Subtype, levels = 1:6, 
                      labels = c("subtype1", "subtype2", "subtype3", "subtype4", "subtype5", "subtype6"))


# automate the SVM , RFE and MCCV
automate_SVM_RFE_MCCV_Process_V1 <- function(data,
                                             NoItr, # run 100 times 
                                             cores, # number of cores to use 
                                             kfold, # Either the number of folds or number of resampling iterations
                                             Numrepeats, # For repeated k-fold cross-validation only
                                             #SubSetSize = 75, # on how many subset to test
                                             minFeaturesize, # minim feature size a subset to have
                                             topfeaturesize) {
  
  # step 1a: save or store the following object as list in each iteration
  vimp_list <- list()
  SankeyPlot_top200_list <- list()
  confusionMatrix_top200_list <- list()
  Sankey_ggplot_top200_list <- list()
  accuracy_train_top200_list <- list()
  accuracy_test_top200_list <- list()
  accuracy_train_bestSet_list <- list()
  accuracy_test_bestSet_list <- list()
  SankeyPlot_bestSet_list <- list()
  ggplot_REF_result_bestSet_list <- list()
  PredictedSamples <-list()
  accuracy_test_top200_list2 <- list()
  confusionMatrix_top200_list2 <- list()
  accuracy_train_top200_list2 <- list()
  
  
  # step 1b: run the iteration
  for (i in 1:NoItr) {
    set.seed(i)
    indexes <-
      createDataPartition(y = data$Subtype,
                          p = 0.8,
                          list = FALSE)
    trainData <- data[indexes,]
    testData <- data[-indexes,]
    
    
    trainData$Subtype <- as.factor(trainData$Subtype)
    testData$Subtype <- as.factor(testData$Subtype)
    
    
    trainVar <- trainData[,-ncol(trainData)]
    responseVar_train <- trainData$Subtype
    
    
    testVar <- testData[,-ncol(testData)]
    responseVar_test <- testData$Subtype
    
    
    registerDoParallel(cores = cores)
    numFeatures <- ncol(trainData) - 1
    selectedSizes = c(seq(numFeatures, minFeaturesize, by = -10), minFeaturesize)
    
    
    # step 1c: Define the control parameters for RFE
    modelType = c("svmLinear", "svmPoly", "svmRadial")
    rfeResults_svmLinear <- rfe(
      x = trainVar,
      y = responseVar_train,
      sizes = selectedSizes,
      rfeControl = rfeControl(
        functions = caretFuncs,
        method = "cv",
        number = kfold,
        repeats = Numrepeats,
        verbose = FALSE,
        allowParallel = TRUE
      ),
      method = modelType[1]
    )
    

    # step 2b : Identify the important features in the best subset
    vimp_list[[i]] <- data.frame(varImp(rfeResults_svmLinear))
    vimp_list[[i]] <- vimp_list[[i]][order(vimp_list[[i]][["Overall"]], decreasing = TRUE), , drop = FALSE]
    vimp_list[[i]][["var"]] <- rownames(vimp_list[[i]])
    bestVar = predictors(rfeResults_svmLinear) #rfeResults_svmLinear$optVariables
    
    
    if (all(bestVar %in% rownames(vimp_list[[i]]))) {
      vimp_list[[i]]  = cbind(vimp_list[[i]][bestVar, , drop = F], iteration = i)
      vimp_list[[i]]  = na.omit(vimp_list[[i]][seq(1:length(bestVar)), ])
      rownames(vimp_list[[i]]) = 1:nrow(vimp_list[[i]])
    } else {
      bestVar = bestVar[1:length(bestVar)]
    }
    
    
    vimp_list=list.remove(vimp_list, 4)
    
    
    # step 2c: plot the accuracy vs the different variable size using ggplot2
    ggplot_REF_result_bestSet_list[[i]] <-
      ggplot(rfeResults_svmLinear) +
      theme(
        plot.title = element_text(size = 16, hjust = 0.5),
        axis.title.x = element_text(size = 14, hjust = 0.5),
        axis.title.y = element_text(size = 14, hjust = 0.5)
      ) +
      geom_vline(xintercept = length(predictors(rfeResults_svmLinear))) +
      labs(y = "Accuracy",
           x = "Variables",
           title = "SVM linear recursive feature selection") +
      annotate(
        "text",
        x = length(predictors(rfeResults_svmLinear)),
        y = 0.95,
        label = paste(length(predictors(
          rfeResults_svmLinear
        )), "best variables"),
        color = "red"
      )
    
    
    # step 2d: Extract the accuracy score of the best training subset
    bestSubset = rfeResults_svmLinear[["bestSubset"]]
    accuracy_train_bestSet_list[[i]] <-
      cbind(
        subsetSize = rfeResults_svmLinear[["results"]][["Variables"]],
        accuracy = rfeResults_svmLinear[["results"]][["Accuracy"]],
        iteration = i
      ) %>% data.frame() %>%
      dplyr::filter(subsetSize == bestSubset)
    
    
    # step 2f: Extract the accuracy score of the best test subset size
    testData_bestSet = testData[, c("Subtype", predictors(rfeResults_svmLinear))]
    accuracy_test_bestSet_list[[i]] <- cbind(accuracy = confusionMatrix(
      predict(rfeResults_svmLinear, testData_bestSet),
      responseVar_test)[["overall"]][["Accuracy"]], iteration = i)
    
    
    # step 3a: predict the test data(Subset the selected features in the test datasets)
    testData_bestSet$Predicted <-
      predict(rfeResults_svmLinear, testData_bestSet)
    SankeyPlot_bestSet_list[[i]] <- cbind(
      data.frame(testData_bestSet[, c("Predicted", "Subtype")]), iteration = i)
    colnames(SankeyPlot_bestSet_list[[i]]) = c("Predicted", "Actual", "iteration")
    
    
    RetrainData = trainData[, c("Subtype", predictors(rfeResults_svmLinear))]
    RetestData =  testData[, c("Subtype", predictors(rfeResults_svmLinear))]
    
    
    trainControl <- caret::trainControl(method = "cv",
                                        number =  5,
                                        allowParallel = TRUE)
    
    
    svm_model <- caret::train(Subtype ~ .,
                              data = RetrainData,
                              trControl = trainControl,
                              method = modelType[1])
    
    
    #try alternative svm model
    svm_model2 <- svm(Subtype ~ ., 
                      data = RetrainData,
                      kernel = "linear",
                      probability = TRUE, 
                      decision.values = TRUE)
    
    
    #Make class predictions and get decision values on PreddNorm
    predictions2 <- predict(svm_model2, 
                            RetestData, 
                            decision.values = TRUE, 
                            probability = T)
    
    
    # step 3c: do the second prediction for top 200 features
    predictions <- predict(svm_model, RetestData)
    RetestData$Predicted <- predict(svm_model, RetestData)
    
    #add to retest data
    RetestData$predictions2 = predictions2
    
    probability=as.data.frame(attr(predictions2, "probabilities"))
    RetestData = cbind(RetestData, probability)
    
    decision_values = as.data.frame(attr(predictions2, "decision.values"))
    RetestData = cbind(RetestData, decision_values)
    
    
    # Function to calculate the ratio
    calculate_ratio <- function(row) {
      # Get scores as a vector from columns 1 to 6
      scores <- as.numeric(row[(ncol(RetestData)-20):(ncol(RetestData)-15)])
      
      # Sort scores to find the highest and second-highest
      sorted_scores <- sort(scores, decreasing = TRUE)
      highest_score <- sorted_scores[1]
      second_highest_score <- sorted_scores[2]
      
      # Calculate the ratio
      result <- (highest_score - second_highest_score) 
      return(result)
    }
    
    RetestData$SVMConfidence <- apply(RetestData, 1, calculate_ratio)
    RetestData$max_value = apply(RetestData[c((ncol(RetestData)-21):(ncol(RetestData)-16))], 1, max)
    
    # step 3d: return the following object as list in each iteration
    SankeyPlot_top200_list[[i]] <- cbind(data.frame(RetestData[, c("Predicted", "predictions2", "Subtype", "SVMConfidence","max_value")]),iteration = i)
    PredictedSamples[[i]] <- cbind(data.frame(rownames(RetestData), iteration = i))
    colnames(SankeyPlot_top200_list[[i]]) = c("Predicted", "predictions2", "Actual","iteration")
    confusionMatrix_top200_list[[i]] = confusionMatrix(RetestData$Predicted , reference = RetestData$Subtype)
    accuracy_train_top200_list[[i]] <- cbind(Accuracy = as.numeric(svm_model[["results"]][["Accuracy"]]),iteration = i)
    accuracy_test_top200_list[[i]] <- as.data.frame(cbind(Accuracy = as.numeric(confusionMatrix(predictions, RetestData$Subtype)[["overall"]][["Accuracy"]],iteration = i)))
    accuracy_test_top200_list2[[i]] <- as.data.frame(cbind(Accuracy = as.numeric(confusionMatrix(predictions2, RetestData$Subtype)[["overall"]][["Accuracy"]],iteration = i)))
    confusionMatrix_top200_list2[[i]] = confusionMatrix(RetestData$predictions2 , reference = RetestData$Subtype)
    accuracy_train_top200_list2[[i]] <- cbind(Accuracy = as.numeric(svm_model2[["results"]][["Accuracy"]]),iteration = i)
    
    
    # step 3f: plot sankey plot
    Sankey_ggplot_top200_list[[i]] <-
      ggplot(SankeyPlot_top200_list[[i]], aes(axis1 = Actual, axis2 = Predicted)) +
      geom_alluvium(aes(fill = Predicted)) +
      geom_stratum() +
      geom_text(stat = "stratum", aes(label = after_stat(stratum))) +
      theme_minimal() +
      ggtitle("Sankey Plot of Actual vs. Predicted Classes") +
      xlab("Classes") +
      ylab("") +
      theme(axis.text.x = element_blank())
    
    print(i)
    
  }
  
  
  # step 3g: return the data
  list(
    "vimp_list" = do.call(rbind, vimp_list),
    "SankeyPlot_top200_list" = do.call(rbind, SankeyPlot_top200_list),
    "confusionMatrix_top200_list" = confusionMatrix_top200_list,
    "Sankey_ggplot_top200_list" = Sankey_ggplot_top200_list,
    "accuracy_train_top200_list" = do.call(rbind, accuracy_train_top200_list),
    "accuracy_test_top200_list" = do.call(rbind, accuracy_test_top200_list),
    "accuracy_train_bestSet_list" = do.call(rbind, accuracy_train_bestSet_list),
    "accuracy_test_bestSet_list" = do.call(rbind, accuracy_test_bestSet_list),
    "SankeyPlot_bestSet_list" = do.call(rbind,SankeyPlot_bestSet_list),
    "ggplot_REF_result_bestSet_list" = ggplot_REF_result_bestSet_list,
    "PredictedSamples" = do.call(rbind, PredictedSamples),
    "accuracy_test_top200_list2" = do.call(rbind, accuracy_test_top200_list2),
    "confusionMatrix_top200_list2" = do.call(rbind, confusionMatrix_top200_list2),
    "accuracy_train_top200_list2" = do.call(rbind , accuracy_train_top200_list2)
    
  )
}


#Run function
ImpsvmREF_result3 <- 
  automate_SVM_RFE_MCCV_Process_V1(data = data, 
                                   kfold = 5,
                                   cores = 3,
                                   Numrepeats = 1,
                                   minFeaturesize = 75,
                                   #topfeaturesize = 250,
                                   NoItr = 100)




Impvimp=as.data.frame(ImpsvmREF_result3$vimp_list)
ImpvimpTable=as.data.frame(table(Impvimp$var))
ImpvimpTable=ImpvimpTable[order(ImpvimpTable$Freq, decreasing = T),]
rownames(ImpvimpTable)=NULL
ImpVimpFinal=ImpvimpTable[1:350,]
write.table(ImpVimpFinal, "/Users/george.field/Desktop/Data Analysis/Project 1/timsTOF_KTSP/Current1/ImpSVM_FeaturesComplete.txt", row.names = F)
