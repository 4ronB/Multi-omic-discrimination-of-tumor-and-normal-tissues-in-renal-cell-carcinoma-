#load SVM model for ccRCC tissue prediction
load("svmModelForTNdiff.Rdata")
#load example dataset
#the matrix should have at least 9 columns with the 9 previously identified gene(PLOD2, PFKP,  PLIN2, IGFBP3, VEGFA, P4HA1, CCND1, VIM, ANXA1)
#the algorithm can accept more variables if all the 9 genes are inlcuded, rows are the patient identifiers, minimum accepted patient number is 1
d <- read.table(file = "exampleTable.txt", header = T, sep = "\t", stringsAsFactors = T)
#using the predict funtction the algorithm can differentiate between normal and malignant tissue
result <- predict(svmProfFit, newdata = d)
#the reuslt of the prediction is a dataframe, where each row represents a patient, the first column describes the most probable tissue type
#the second and third column gives the probability values of normal and tumor occurances
result

