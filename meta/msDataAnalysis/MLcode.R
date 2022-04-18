




#########################building ML model for tumor data prediction

library(caret)
library(e1071)
protExp <- read.table(file = "meta/SE-MS/MSexpressionExpTableForML.txt", header = T, sep = "\t", stringsAsFactors = T)

#remove patient with missing CCND1 intensity value
protExp <- protExp[-162,]
#remove HLA-genes
protExp <- protExp[,-c(4,7)]

set.seed(5200)
inTrain <- createDataPartition(
  y = protExp$Type,
  p = .7,
  list = FALSE
)

training <- protExp[ inTrain,]
testing  <- protExp[-inTrain,]

set.seed(5200)
modelControl <- trainControl(method = "repeatedcv",
                             number = 10,
                             repeats = 5,
                             classProbs = TRUE
                             
)

set.seed(5200)
rfeCont <- rfeControl(functions = caretFuncs,
                      method = "repeatedcv",
                      number = 10,
                      repeats = 5,
                      verbose = FALSE)
######training
#svm
set.seed(5200)#5200
svmProfFit <- rfe(training[,2:ncol(training)], training$Type,
                  method="svmLinear",
                  sizes = c(1:20),#simulation will fit models with subset sizes of 1:20
                  rfeControl = rfeCont,
                  trControl = modelControl)

svmProfFit
svmProfFit$fit
head(svmProfFit$resample)


predictors(svmProfFit)

testerSvm <- testing[,colnames(testing) %in% c("Type" ,predictors(svmProfFit))] 
testClassesProSvm <- predict(svmProfFit, newdata = testing)
testSumProSvm <- confusionMatrix(data = testClassesProSvm$pred, reference = testerSvm$Type, positive = "Tumor")
testSumProSvm

#knn
set.seed(5200)#5200
knnProfFit <- rfe(training[,2:ncol(training)], training$Type,
                  method="knn",
                  sizes = c(1:20),
                  rfeControl = rfeCont,
                  trControl = modelControl)

knnProfFit
knnProfFit$fit
head(knnProfFit$resample)


predictors(knnProfFit)

testerKnn <- testing[,colnames(testing) %in% c("Type" ,predictors(knnProfFit))] 
testClassesProKnn <- predict(knnProfFit$fit, newdata = testing)
testSumProKnn <- confusionMatrix(data = testClassesProKnn, testing$Type, positive = "Tumor")
testSumProKnn


#rf
set.seed(5200)#5200
rfProfFit <- rfe(training[,2:ncol(training)], training$Type,
                 method="rf",
                 sizes = c(1:20),
                 rfeControl = rfeCont,
                 trControl = modelControl)

rfProfFit
rfProfFit$fit
head(rfProfFit$resample)

predictors(rfProfFit)

testClassesProRf <- predict(rfProfFit$fit, newdata = testing)
testSumProRf <- confusionMatrix(data = testClassesProRf, testing$Type, positive = "Tumor")
testSumProRf


#Logit
set.seed(5200)#5200
logitProfFit <- rfe(training[,2:ncol(training)], training$Type,
                    method="LogitBoost",
                    sizes = c(1:20),
                    rfeControl = rfeCont,
                    trControl = modelControl)

logitProfFit
logitProfFit$fit
head(logitProfFit$resample)

predictors(logitProfFit)


testClassesPrologit <- predict(logitProfFit, newdata = testing)
testSumPrologit <- confusionMatrix(data = testClassesPrologit$pred, testing$Type, positive = "Tumor")
testSumPrologit



#compare models
resamps <- resamples(list(RF = rfProfFit,
                          SVM = svmProfFit,
                          KNN = knnProfFit,
                          LOGIT = logitProfFit))

summary(resamps)
testRunsSum <- data.frame("RF" = c(testSumProRf$overall[1],testSumProRf$overall[2], testSumProRf$byClass[1], testSumProRf$byClass[2]),
                          "SVM" = c(testSumProSvm$overall[1],testSumProSvm$overall[2], testSumProSvm$byClass[1], testSumProSvm$byClass[2]),
                          "KNN" = c(testSumProKnn$overall[1],testSumProKnn$overall[2], testSumProKnn$byClass[1], testSumProKnn$byClass[2]),
                          "LOGIT" = c(testSumPrologit$overall[1],testSumPrologit$overall[2], testSumPrologit$byClass[1], testSumPrologit$byClass[2])
)

testRunsSum

dotplot(resamps)
dotplot(resamps, metric = "Accuracy")

save(rfProfFit, file = "rf_2ModelForTNdiff.Rdata")
save(svmProfFit, file = "svmModelForTNdiff.Rdata")
save(knnProfFit, file = "knnModelForTNdiff.Rdata")
save(logitProfFit, file = "logitModelForTNdiff.Rdata")

plotSvmTrainAcc <- plot(svmProfFit, type = c("g", "o"), main = "SVM", ylim = c(0.92,1))
plotKnnTrainAcc <- plot(knnProfFit, type = c("g", "o"), main = "KNN", ylim = c(0.92,1))
plotRfTrainAcc <- plot(rfProfFit, type = c("g", "o"), main = "RF", ylim = c(0.92,1))
plotLogitTrainAcc <- plot(logitProfFit, type = c("g", "o"), main = "LOGIT", ylim = c(0.92,1))
