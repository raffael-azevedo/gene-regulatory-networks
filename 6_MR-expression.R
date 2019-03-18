library(affy)
library(limma)
library(dplyr)
library(tidyr)
library(ggplot2)
library(EnvStats)
load("/home/raffael/GDrive/Doutorado/Sepsis R/Paper/SIGNATURE_GSE26378/validation/commonRegulons.RData")
regulons <- commonGenes[intersection]
setwd("/home/raffael/Downloads/GSE4607_RAW/")
data <- ReadAffy()
exp <- rma(data)
GSE4607 <- exprs(exp)
#names_index <- setNames(c(rep("control", 15), rep("shock_ns", 9), rep("shock_s", 33)),colnames(GSE4607))
names_index <- setNames(c(rep("control", 15), rep("shock", 42)),colnames(GSE4607))
#colnames(GSE4607) <- c(rep("control", 15), rep("shock_ns", 9), rep("shock_s", 33))
a <- as.data.frame(rowMeans(GSE4607[,1:15]))
b <- as.data.frame(rowMeans(GSE4607[,16:57]))
# b <- as.data.frame(rowMeans(GSE4607[,16:24]))
# c <- as.data.frame(rowMeans(GSE4607[,25:57]))
dataset <- cbind(a,b)
colnames(dataset) <- c("control", "shock")
load("/home/raffael/GDrive/Doutorado/Sepsis R/Paper/SIGNATURE_GSE26378/GSE4607/RData/chunk1_nrData_GSE4607.RData")
PROBEID <- rownames(dataset)
rownames(dataset) <- NULL
dataset <- cbind(PROBEID, dataset)
dataset <- merge(dataset,gpl_new,by="PROBEID")
dataset <- dataset[,c(1,5,6,2,3,4)]
dataset <- dataset %>% mutate_if(is.factor, as.character)
rm(gpl_new, exp, HSexp, PROBEID,a,b,c,data)

temp <- as.data.frame(GSE4607)
temp$PROBEID <- rownames(temp)
#temp <- merge(temp, dataset[,1:3], by="PROBEID")
temp <- merge(temp, dataset[,c(1,4,5)], by="PROBEID")
#colnames(temp) <- c("PROBEID",rep("control", 15), rep("shock_ns", 9), rep("shock_s", 33),"ENTREZ", "SYMBOL")
comparison <- temp %>% gather(key = "key", value = "value", -PROBEID,-SYMBOL,-ENTREZ)
#comparison$key <- sapply(strsplit(comparison[,"key"],"\\."),"[",1)
comparison[,"tipo"] <- names_index[comparison[,"key"]]

best_probes <- comparison %>%
  group_by(PROBEID) %>%
  mutate(score = geoSD(value)^(1/geoMean(value))) %>%
  group_by(SYMBOL) %>%
  mutate(best = score==max(score))
#ou
#filter(score==max(score))
best_probes <- subset(best_probes, best == TRUE)

for(gene in intersection){
  print(gene)
  gene_subset <- subset(best_probes, SYMBOL == gene)
  pairwise.wilcox.test(gene_subset$value, gene_subset$tipo, p.adjust.method = "bonferroni") %>% print
}
ggplot(subset(best_probes, SYMBOL%in%intersection)) + geom_boxplot(aes(x = tipo, y = value)) + facet_wrap(SYMBOL~., nrow = 3,ncol = 5)