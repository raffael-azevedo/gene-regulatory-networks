# This script is not part of the gene regulatory network analysis (GRNA); Instead, this is an additional step.
# It is only useful when working with two or more network. This is interesting when one knows which are the master regulators
# common between two networks, and wants to identify regulated genes present in both regulons, in both networks.
# As this 

# chunk 0.1.1: get regulated genes for each network. ####
# GSE13904

load("~/GDrive/Doutorado/Sepsis R/Paper/SIGNATURE_GSE26378/GSE13904/RData/chunk3_TNI-TNA_analysis.RData")
sub.tnet <- tnet[,colnames(tnet)%in%mra$PROBEID]
symbols.subtnet <- c("ZNF331","KLF12","RORA","ZNF544","GATA3","ZNF551","IRF4","GCFC2",
         "ZNF253","ZNF529","ZNF235","SMAD3","TULP4","NKRF","ZNF234",
         "ZNF706","NR3C2","ZNF141","ZNF10","NFE2L3","REXO4","DBP","ZNF134",
         "ZKSCAN8","ZNF202","ZNF329","HOXB2","ZBTB25","CIITA","ZNF510",
         "SOX12","ZNF16","JUNB","TRIM25","PHTF1","MAFG","EPAS1","KLF5",
         "KLF7","MEF2A","MTF1","NFIL3","GAS7","CEBPB","HES1","FOSL2",
         "ZNF467","BCL6","RFX2")
idx.tnet <- subset(tna, tna$Regulon%in%symbols.subtnet)
PROBEID <- rownames(idx.tnet)
rownames(idx.tnet) <- NULL
idx.tnet <- cbind(PROBEID,idx.tnet)
colnames(idx.tnet)[c(1,2,9)] <- c("PROBEID", "SYMBOL","Adjusted.Pvalue")
sub.tnet <- tnet[,colnames(tnet)%in%idx.tnet$PROBEID]


regulated.genes <- list()
for(i in 1:ncol(sub.tnet))
{
  diff <- which(sub.tnet[,i] != 0)
  regulated.genes[[i]] <- sub.tnet[diff,i]
  temp <- as.data.frame(regulated.genes[[i]])
  colnames(temp) <- colnames(sub.tnet)[i]
  regulated.genes[[i]] <- temp
  rm(temp)
}

# chunk 0.1.2: load data and check for inconsistencies ####
GSE13904 <- regulated.genes
GSE13904_subtnet <- symbols.subtnet

# chunk 0.2.1: get regulated genes for each network. ####
# GSE4607

load("~/GDrive/Doutorado/Sepsis R/Paper/SIGNATURE_GSE26378/GSE4607/RData/chunk3_TNI-TNA_analysis.RData")
sub.tnet <- tnet[,colnames(tnet)%in%mra$PROBEID]
symbols.subtnet <- c("ZNF331","KLF12","RORA","ZNF544","GATA3","ZNF551","IRF4","GCFC2",
                     "ZNF253","ZNF529","ZNF235","SMAD3","TULP4","NKRF","ZNF234",
                     "ZNF706","NR3C2","ZNF141","ZNF10","NFE2L3","REXO4","DBP","ZNF134",
                     "ZKSCAN8","ZNF202","ZNF329","HOXB2","ZBTB25","CIITA","ZNF510",
                     "SOX12","ZNF16","JUNB","TRIM25","PHTF1","MAFG","EPAS1","KLF5",
                     "KLF7","MEF2A","MTF1","NFIL3","GAS7","CEBPB","HES1","FOSL2",
                     "ZNF467","BCL6","RFX2")
idx.tnet <- subset(tna, tna$Regulon%in%symbols.subtnet)
PROBEID <- rownames(idx.tnet)
rownames(idx.tnet) <- NULL
idx.tnet <- cbind(PROBEID,idx.tnet)
colnames(idx.tnet)[c(1,2,9)] <- c("PROBEID", "SYMBOL","Adjusted.Pvalue")
sub.tnet <- tnet[,colnames(tnet)%in%idx.tnet$PROBEID]


regulated.genes <- list()
for(i in 1:ncol(sub.tnet))
{
  diff <- which(sub.tnet[,i] != 0)
  regulated.genes[[i]] <- sub.tnet[diff,i]
  temp <- as.data.frame(regulated.genes[[i]])
  colnames(temp) <- colnames(sub.tnet)[i]
  regulated.genes[[i]] <- temp
  rm(temp)
}

GSE4607 <- regulated.genes
GSE4607_subtnet <- symbols.subtnet
rm(regulated.genes)
rm(symbols.subtnet)

intersection <- intersect(GSE13904_subtnet, GSE4607_subtnet)
load("~/GDrive/Doutorado/Sepsis R/Paper/SIGNATURE_GSE26378/GSE13904/RData/chunk1_nrData.RData")
rm(HSexp)
dictionary <- subset(gpl_new, gpl_new$SYMBOL%in%intersection)
dictionary$PROBEID <- as.character(dictionary$PROBEID)
dictionary$ENTREZ <- as.numeric(dictionary$ENTREZ)
dictionary$SYMBOL <- as.character(dictionary$SYMBOL)

#check
length(unique(dictionary$ENTREZ))==length(intersection)

# lasso 1: creates a dataframe for every master regulator present in the intersection of the two network and name the column after the MR ####
GSE13904_intesect <- list()
for (i in 1:length(GSE13904))
{
  temp <- colnames(GSE13904[[i]])[1]
  colnames(GSE13904[[i]])[1] <- dictionary[temp,3]
  rm(temp)
  if (!is.na(colnames(GSE13904[[i]][1])))
  {
    GSE13904_intesect[[i]] <- GSE13904[[i]]
  }
}

#save(GSE13904_intesect,dictionary, GSE13904_subtnet, file="~/GDrive/Doutorado/Sepsis R/Paper/SIGNATURE_GSE26378/GSE13904/RData/chunk2_regulatedGenes-subset.RData")

# lasso 2: the same of lasso 1, but for GSE4607 network ####
GSE4607_intesect <- list()
for (i in 1:length(GSE4607))
{
  temp <- colnames(GSE4607[[i]])[1]
  colnames(GSE4607[[i]])[1] <- dictionary[temp,3]
  rm(temp)
  if (!is.na(colnames(GSE4607[[i]][1])))
  {
    GSE4607_intesect[[i]] <- GSE4607[[i]]
  }
}

a <- c()
for(i in 1:length(GSE13904_intesect)){a[i] <- colnames(GSE13904_intesect[[i]])}
b <- c()
for(i in 1:length(GSE4607_intesect)){b[i] <- colnames(GSE4607_intesect[[i]])}
length(a)==length(b)

#save(GSE4607_intesect,GSE4607_subtnet,dictionary,file="~/GDrive/Doutorado/Sepsis R/Paper/SIGNATURE_GSE26378/GSE4607/RData/chunk2_regulatedGenes-subset_GSE4607.RData")

# lasso 3: this lasso formats the data into several forms, described below ####
# percentages: a named vector with the percentage of shared genes between regulons in the two networks;
# commonGenes: a list of dataframes, each one containing only the common up (or down) regulated, shared among the two networks;
# check regulated genes of inteserction
checkList <- list()
percentages <- c()
TFs <- c()
commonGenes <- list()
allGenes <- data.frame()
for (i in 1:length(GSE13904_intesect))
{
  temp1 <- GSE13904_intesect[[i]]
  temp2 <- GSE4607_intesect[[i]]
  allGenes <- as.data.frame(unique(c(rownames(temp1), rownames(temp2))))
  colnames(allGenes)[1] <- "allGenes"
  allGenes$GSE13904 <- 0
  allGenes$GSE4607 <- 0
  allGenes$allGenes <- as.character(allGenes$allGenes)
  for (j in 1:nrow(allGenes))
  {
    if (rownames(temp1)%in%allGenes$allGenes)
    {
      allGenes$GSE13904[match(rownames(temp1)[j],allGenes$allGenes)] <- 1
    }
    if (rownames(temp2)%in%allGenes$allGenes)
    {
      allGenes$GSE4607[match(rownames(temp2)[j],allGenes$allGenes)] <- 1
    }
  }
  colnames(allGenes)[1] <- colnames(GSE13904_intesect[[i]])
  checkList[[i]] <- allGenes
  TFs[i] <- colnames(GSE13904_intesect[[i]])
  percentages[i] <- (sum(rowSums(checkList[[i]][,c(2,3)])==2)/nrow(allGenes))*100
  names(percentages)[i] <- TFs[i]
  commonGenes[[i]] <- as.character(allGenes[which(rowSums((checkList[[i]][,c(2,3)]))==2),1])
  names(commonGenes) <- TFs
  a <- as.data.frame(commonGenes[[i]])
  a[,1] <- as.character(a[,1])
  colnames(a)[1] <- names(commonGenes[i])
  a <- merge(a, gpl_new, by.x = colnames(a)[1], by.y = "PROBEID")
  a <- merge(a, GSE4607_intesect[[i]], by.x = colnames(a)[1], by.y = "row.names")
  colnames(a)[4] <- "Regulation"
  commonGenes[[i]] <- a
}

library(dplyr)
commonGenes <- lapply(commonGenes, function(df) mutate_at(df, vars("ENTREZ"), as.numeric))
commonGenes <- lapply(commonGenes, function(df) mutate_at(df, vars("SYMBOL"), as.character))
commonGenes <- lapply(commonGenes, function(df) mutate_at(df, vars("Regulation"), as.numeric))

#save(percentages, TFs, intersection, GSE13904_intesect, GSE4607_intesect, commonGenes, dictionary, file="~/GDrive/Doutorado/Sepsis R/Paper/SIGNATURE_GSE26378/GSE13904/RData/chunk3_commonGenes_allDatasets.RData")
#barplot(sort(percentages, decreasing = T), ylim = c(0,40), las=2)
#duplicatedGenes <- allGenes[duplicated(allGenes)]

# chunk 2: MRs that separates groups within regulon activity ####

symbols.subtnet <- c("ZNF331","KLF12","RORA","ZNF544","GATA3","ZNF551","IRF4","GCFC2",
                     "ZNF253","ZNF529","ZNF235","SMAD3","TULP4","NKRF","ZNF234",
                     "ZNF706","NR3C2","ZNF141","ZNF10","NFE2L3","REXO4","DBP","ZNF134",
                     "ZKSCAN8","ZNF202","ZNF329","HOXB2","ZBTB25","CIITA","ZNF510",
                     "SOX12","ZNF16","JUNB","TRIM25","PHTF1","MAFG","EPAS1","KLF5",
                     "KLF7","MEF2A","MTF1","NFIL3","GAS7","CEBPB","HES1","FOSL2",
                     "ZNF467","BCL6","RFX2")
percentages_intersect <- percentages[names(percentages)%in%symbols.subtnet]
barplot(sort(percentages_intersect, decreasing = T), ylim = c(0,40), las=2)

# chunk 3: check within other diseases signatures for similarities in regulon activity ####
library(data.table)
outMRs <- list()
# multiple sclerosis - whole blood
gse4607_sig_gse21942 <- fread("/home/raffael/GDrive/Doutorado/Sepsis R/Paper/SIGNATURE_GSE21942/GSE4607_SIG_GSE21942.txt")
outMRs[["gse4607_sig_gse21942"]] <- symbols.subtnet[which(!symbols.subtnet%in%gse4607_sig_gse21942$Regulon)]
# rheumatoid arthritis - synovial tissue
gse4607_sig_gse48780 <- fread("/home/raffael/GDrive/Doutorado/Sepsis R/verify_step2/mra_candidates_GSE48780.txt",drop = 1)
outMRs[["gse4607_sig_gse48780"]] <- symbols.subtnet[which(!symbols.subtnet%in%gse4607_sig_gse48780$Regulon)]
# sepsis (network 2) - sepsis signature (skeletal muscle)
gse13904_sig_gse13205 <- fread("~/GDrive/Doutorado/Sepsis R/Paper/SIGNATURE_GSE13205/GSE13904/GSE13904_SIG_GSE13205.txt")
outMRs[["gse13904_sig_gse13205"]] <- symbols.subtnet[which(!symbols.subtnet%in%gse13904_sig_gse13205$Regulon)]
# sepsis (network 1) - sepsis signature (skeletal muscle)
gse4607_sig_gse13205 <- fread("~/GDrive/Doutorado/Sepsis R/Paper/SIGNATURE_GSE13205/GSE4607/GSE4607_SIG_GSE13205.txt")
outMRs[["gse4607_sig_gse13205"]] <- symbols.subtnet[which(!symbols.subtnet%in%gse4607_sig_gse13205$Regulon)]

intersection <- Reduce(intersect, outMRs)
barplot(sort(percentages[intersection], decreasing = T), ylim = c(0,40), las=2)
save(commonGenes, intersection, file = "/home/raffael/GDrive/Doutorado/Sepsis R/Paper/SIGNATURE_GSE26378/validation/commonRegulons.RData")

# chunk 4: compare genes within the signatures ####
setwd("/home/raffael/GDrive/Doutorado/Sepsis R/Paper/VALIDAÇÃO/")
load("/home/raffael/GDrive/Doutorado/Sepsis R/GSE26378/SIGNATURE/GSE26378_signature.RData")
GSE26378 <- subset(GSE26378_phenoIDs, GSE26378_phenoIDs$PROBEID%in%GSE26378_hits)
write.table(GSE26378,file = "GSE26378.txt", quote = F, sep = "\t", row.names = F, col.names = T)
rm(gpl_new, GSE26378_hits, GSE26378_pheno)
load("/home/raffael/GDrive/Doutorado/Sepsis R/verify_step2/GSE21942_signature.RData")
GSE21942 <- subset(GSE21942_pheno, GSE21942_pheno$PROBEID%in%GSE21942_hits)
write.table(GSE21942,file = "GSE21942.txt", quote = F, sep = "\t", row.names = F, col.names = T)
rm(gpl_new, GSE21942_hits, GSE21942_pheno)
# load("/home/raffael/GDrive/Doutorado/Sepsis R/verify_step2/GSE51808_signature.RData")
# GSE51808 <- subset(GSE51808_pheno, GSE51808_pheno$PROBEID%in%GSE51808_hits)
# write.table(GSE51808,file = "GSE51808.txt", quote = F, sep = "\t", row.names = F, col.names = T)
# rm(gpl_new, GSE51808_hits, GSE51808_pheno)
load("/home/raffael/GDrive/Doutorado/Sepsis R/verify_step2/chunk2_RTN_signatureGSE48780.RData")
GSE56649 <- subset(GSE56649_pheno, GSE56649_pheno$PROBEID%in%GSE56649_hits)
write.table(GSE56649,file = "GSE56649.txt", quote = F, sep = "\t", row.names = F, col.names = T)
rm(dataRTN, GSE56649_hits, GSE56649_pheno)


intersection2 <- Reduce(intersect, list(GSE21942$SYMBOL, GSE26378$SYMBOL, GSE56649$SYMBOL))


# chunk 5: differential functional enrichment ####
# this chunk gets the differentially expressed genes for (at least) the 15 MRs and identify the metabolic pathways and other
# functional enrichment gibberish

library(enrichR)
setwd("~/GDrive/Doutorado/Sepsis R/Paper/SIGNATURE_GSE26378/Funcional_Enrichment_GO-BP2018/subset")
dbs <- listEnrichrDbs()
dbs <- "GO_Biological_Process_2018"
temp_subset <- commonGenes[intersection]
for(i in 1:length(temp_subset)){
  enrichment <- enrichr(c(temp_subset[[i]][[3]]), dbs)
  gobp <- enrichment[[1]] %>% filter(Adjusted.P.value < 0.05)
  write.table(gobp, names(temp_subset[i]),quote = F, sep = "\t", row.names = F, col.names = T)
}
