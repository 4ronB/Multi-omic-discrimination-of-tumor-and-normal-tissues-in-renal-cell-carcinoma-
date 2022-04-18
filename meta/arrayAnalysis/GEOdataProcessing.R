library(readxl)
#open clinical data
clinTable <- read_xlsx(path = "GEOclinical.xlsx", sheet = 1)
clinTable <- as.data.frame(clinTable)
#separate normal and malignant patients
Cc <- clinTable[clinTable$Pair != "0",1]
Norm <- clinTable[clinTable$Pair != "0",7]
#load normalized expression table
load("kidney-clean-ann.RData")
expTable <- txtFile
rm(txtFile)
geneNames <- expTable$gene_sym
row.names(expTable) <- expTable$gene_sym
#make tumor and normal specific expression tables 
tumExp <- expTable[,c(colnames(expTable) %in% Cc)]
normExp <- expTable[,c(colnames(expTable) %in% Norm)]

data <- data.frame()
#prepare a summary table where all the deferentially expressed genes are identified
for (gene in 1:length(geneNames)){
  norm <- t(normExp[rownames(normExp) == geneNames[gene],])
  tum <- t(tumExp[rownames(tumExp) == geneNames[gene],])
  
  if (sum(norm) | sum(tum) > 0){
    mwt <- wilcox.test(x = norm, y = tum, paired = T, exact = F, conf.int = T)
  } else if (sum(norm) | sum(tum) == 0){
    mwt <- data.frame("p.value" = NA, "conf.int" = NA)
  }
  
  fivesNorm <- fivenum(norm)
  fivesTum <- fivenum(tum)
  table <- data.frame("Gene" = geneNames[gene],
                      "Wp" = mwt$p.value, 
                      "CImin" = mwt$conf.int[1], 
                      "CImax" = mwt$conf.int[2], 
                      "Mean-Norm" = mean(norm),
                      "Mean-Tum" = mean(tum),
                      "Med-Norm" = median(norm),
                      "Med-Tum" = median(tum),
                      "FC" = mean(tum)/mean(norm),
                      "Min-Norm" = fivesNorm[1],
                      "Q1-Norm" = fivesNorm[2],
                      "Med-Norm" = fivesNorm[3],
                      "Q3-Norm" = fivesNorm[4],
                      "Max-Norm" = fivesNorm[5],
                      "Min-Tum" = fivesTum[1],
                      "Q1-Tum" = fivesTum[2],
                      "Med-Tum" = fivesTum[3],
                      "Q3-Tum" = fivesTum[4],
                      "Max-Tum" = fivesTum[5])
  
  data <- rbind(data, table)
  
}
#save table
write.table(x = data, file = "kidney_GEO_sumMW_data.txt")

#median expression>1000, FC > 2, p < 0.01
#identification of top  genes

geneChipSum <- read.table(file = "kidney_GEO_sumMW_data.txt", header = T)
geneChipSum$pAdj <- p.adjust(p = geneChipSum$Wp, method = "fdr")
geneChipSumFiltered <- geneChipSum[geneChipSum$pAdj < 0.01 & geneChipSum$Med.Norm >= 1000 & geneChipSum$Med.Tum >= 1000 & geneChipSum$FC >= 2,]
geneChipSumFiltered <- geneChipSumFiltered[order(geneChipSumFiltered$FC,decreasing = T),]

write.table(x = geneChipSumFiltered, file = "topDiffGene-chip.txt", sep = "\t", row.names = F, col.names = T)
