####################################creating protein expression boxplots####################
library(readxl)
library(ggplot2)

RunlevelData <- read.table(file = "kidneyProteinSum.txt", header = T, sep = "\t")
ann <- read.table(file = "proteinAnn.txt", header = T, sep = "\t")

sumTab <- data.frame()
for (gene in 1:length(ann$GeneName)) {
  selectedGene <- ann$GeneName[gene]
  protName <- ann[ann$GeneName == selectedGene, 2]
  
  
  #prepare table for ggplot 
  table <- rbind(data.frame("ExpValues" = RunlevelData[RunlevelData$Protein == protName & RunlevelData$GROUP_ORIGINAL == "Normal",3], "Type" = "Normal"),
                 data.frame("ExpValues" = RunlevelData[RunlevelData$Protein == protName & RunlevelData$GROUP_ORIGINAL == "Tumor",3], "Type" = "Tumor"))
  colnames(table) <- c("ExpValues", "Type")
  
  #statistical analysis
  if (length(table[table$Type == "Normal",1]) == length(length(table[table$Type == "Tumor",1]))){
    sTest <- t.test(x = table[table$Type == "Tumor",1], y = table[table$Type == "Normal",1], paired = T)
  } else if (length(table[table$Type == "Normal",1]) != length(length(table[table$Type == "Tumor",1]))){
    sTest <- t.test(x = table[table$Type == "Tumor",1], y = table[table$Type == "Normal",1], paired = F)
  }
  
  fivesNorm <- fivenum(table[table$Type == "Normal",1])
  fivesTum <- fivenum(table[table$Type == "Tumor",1])
  
  statTab <- data.frame("Gene" = selectedGene,
                        "Protein" = protName,
                        "tValue" = sTest$statistic,
                        "p-Value" = sTest$p.value, 
                        "CImin" = sTest$conf.int[1], 
                        "CImax" = sTest$conf.int[2], 
                        "Mean-Norm" = mean(table[table$Type == "Normal",1]),
                        "Mean-Tum" = mean(table[table$Type == "Tumor",1]),
                        "log2FC" = mean(table[table$Type == "Tumor",1]) - mean(table[table$Type == "Normal",1]),
                        "FC" = 2^mean(table[table$Type == "Tumor",1])/2^mean(table[table$Type == "Normal",1]),
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
  
  sumTab <- rbind(sumTab, statTab)
  
}

sumTab$pAdj <- p.adjust(p = sumTab$p.Value, method = "fdr")
write.table(x = sumTab, file = "proteinSUMMARY.txt", sep = "\t", row.names = F, col.names = T)

#visualization
RunlevelData <- read.table(file = "kidneyProteinSum.txt", header = T, sep = "\t")
protSum <- read.table(file = "proteinSUMMARY.txt", header = T, sep = "\t")

for (gene in 1:length(protSum$Gene)) {
  selectedGene <- protSum$Gene[gene]
  protName <- protSum[protSum$Gene == selectedGene, 2]
  adjP <- protSum[protSum$Gene == selectedGene, 21]
  
  #prepare table for ggplot 
  table <- rbind(data.frame("ExpValues" = RunlevelData[RunlevelData$Protein == protName & RunlevelData$GROUP_ORIGINAL == "Normal",3], "Type" = "Normal"),
                 data.frame("ExpValues" = RunlevelData[RunlevelData$Protein == protName & RunlevelData$GROUP_ORIGINAL == "Tumor",3], "Type" = "Tumor"))
  colnames(table) <- c("ExpValues", "Type")
  
  #visualization
  boxPlot <- ggplot(data = table, aes(x = Type, y = ExpValues, size = Type, color = Type, fill = Type)) +
    geom_boxplot(outlier.alpha = 0) +
    geom_jitter(size = 3) +
    scale_color_manual( values =  c("goldenrod4", "lightsalmon4")) + 
    scale_size_manual(values = c(1,1)) +
    scale_fill_manual(values =  c("goldenrod1", "lightsalmon1")) +
    labs(x = NULL, y = paste0("Log2 intensity of ", selectedGene)) +
    geom_label(label = paste0("P = ", formatC(adjP, format = "e", digits = 2)), 
               fill = "white", color = "black",x =2.3, y = max(table$ExpValues), label.size = 0, size = 6) +
    theme(legend.position = "none", 
          panel.background = element_rect(fill = NA, colour = "grey50"),  
          panel.grid = element_blank(),
          axis.text.x = element_text(size = 26, colour = "black"),
          axis.text.y = element_text(size = 26, colour = "black"),
          axis.title.y = element_text(size = 26, face = "bold",  vjust = 2))
  #boxPlot
  #save the plot
  # ggsave(filename = paste0("plots/", selectedGene, "_log2Intensities_own", ".jpg"), 
  #         plot = boxPlot, width = 7.3 , height = 5.2, dpi = 1000)
  
  
  
}

