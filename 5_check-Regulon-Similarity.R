# This script is not part of the gene regulatory network analysis (GRNA); Instead, this is an additional step.
# It is only useful when working with two or more network. This is interesting when one knows which are the master regulators
# common between two networks, and wants to identify regulated genes present in both regulons, in both networks.
# As this recquires at least two datasets (already treated by previous scripts) and (if possible) analyzed by the
# same platform. 
# DISCLAIMER: I will work with two datasets originated by the same technology (GPL570). Using two different
# technologies will add steps to the analysis, in order to equalize the gene identifiers. The other dataset for contrast
# will be GSE13904, which also is related to sepsis patients. Assuming that you came this far, you have the last RData
# generated at the end of script 4_Regulatory-networks-generation.R, we will just load them to start the work.
load("chunk3_TNI-TNA_analysis_GSE4607.RData")
rm(dataRTN,rtni,GSE4607_tree,rtna,RegScores)
# Chunk 1: Load regulated genes. ####
# For dataset: GSE4607 
# 1.1) Loads the regulated genes of GSE4607
load("~/regulatedGenes.RData")
# 1.2) Assign it to a new variable, since the variables names are the same for both datasets (by default)
GSE4607 <- regulated.genes
# 1.3) Loads the regulated genes of GSE13904
load("~/regulatedGenes.RData")
# 1.4) Assign it to a new variable, since the variables names are the same for both datasets (by default)
GSE13904 <- regulated.genes
# 1.5) This step is pretty much the same as stated in sections 1.1 to 1.3.5. Refer to script 4_Regulatory-networks-generation.R]
# for further questions.
master_regulators <- c("ZNF331","KLF12","RORA","ZNF544","GATA3","ZNF551","IRF4","GCFC2",
                                  "ZNF253","ZNF529","ZNF235","SMAD3","TULP4","NKRF","ZNF234",
                                  "ZNF706","NR3C2","ZNF141","ZNF10","NFE2L3","REXO4","DBP","ZNF134",
                                  "ZKSCAN8","ZNF202","ZNF329","HOXB2","ZBTB25","CIITA","ZNF510",
                                  "SOX12","ZNF16","JUNB","TRIM25","PHTF1","MAFG","EPAS1","KLF5",
                                  "KLF7","MEF2A","MTF1","NFIL3","GAS7","CEBPB","HES1","FOSL2",
                                  "ZNF467","BCL6","RFX2")
mra <- subset(tna, tna$Regulon%in%master_regulators)
PROBEID <- rownames(mra)
rownames(mra) <- NULL
mra <- cbind(PROBEID,mra)
colnames(mra)[c(1,2,9)] <- c("PROBEID", "SYMBOL","Adjusted.Pvalue")
# 1.6) Remember when you were told not to worry about names in regulated.genes object? Don't. We will begin to 
# add some gene identifiers (a.k.a., gene symbols). First for GSE4607
GSE4607_sub <- subset(GSE4607, unlist(lapply(seq_along(regulated.genes), 
                                      function(x){names(regulated.genes[[x]])}))%in%mra$PROBEID)
# 1.6.1) And further for GSE13904
GSE13904_sub <- subset(GSE13904, unlist(lapply(seq_along(regulated.genes),
                                      function(x){names(regulated.genes[[x]])}))%in%mra$PROBEID)
# 1.7) This step is messy, because we load a RData with expression data, but only needs the gpl_new file.
# If you are a better organized person than I, you should save the dictionary object in a different RData
# (or even a scrip, as I did for the TFs).
load("~/GDrive/Doutorado/Sepsis R/Paper/SIGNATURE_GSE26378/GSE13904/RData/chunk1_nrData.RData")
# 1.7) This will not be used, so, to the bin.
rm(HSexp)
# 1.7.1) Here the dictionary *ONLY* for the master regulators is made (which not means to get rid of the full dictionary).
dictionary <- subset(gpl_new, gpl_new$SYMBOL%in%master_regulators)
# 1.7.2) As the columns comes up as factors, is time for tidy it up, turning into: characters,
dictionary$PROBEID <- as.character(dictionary$PROBEID)
# 1.7.3) numerics,
dictionary$ENTREZ <- as.numeric(dictionary$ENTREZ)
# 1.7.4) and characters. Just to work with the proper data type.
dictionary$SYMBOL <- as.character(dictionary$SYMBOL)
# 1.8) Note that the dictionary is bigger than the master_regulators. This is due to more than one probe covering one gene;
# This step does a little sanity check, to see if the unique ENTREZ entries are an exact match of our master_regulators vector.
length(unique(dictionary$ENTREZ))==length(master_regulators)

# Chunk 2: Add the gene names to each data.frame inside the proper regulated.genes list. ####
# 2) This is made in a for loop, which runs for all elements of GSE4607_sub (could be the entire dataset), which gets
# the probe ID name of the master regulator, go to the dictionary and gets the gene symbol corresponding to that probe ID, in a specific position.
# Then, adds the gene symbol to the specific column name.
GSE4607_intersect <- list()
for (i in 1:length(GSE4607_sub))
{
  temp <- colnames(GSE4607_sub[[i]])[1]
  colnames(GSE4607_sub[[i]])[1] <- dictionary[temp,3]
  rm(temp)
  if (!is.na(colnames(GSE4607_sub[[i]][1])))
  {
    GSE4607_intersect[[i]] <- GSE4607_sub[[i]]
  }
}
# 2.2) The steps below (save and load) are checkpoints. If is necessary to re-do the analysis, parts of the
# data is stored in RData and can be loaded to speed up the re-analysis.
save(GSE4607_sub,GSE4607_intersect, dictionary, GSE13904, file="~/regulatedGenes-subsetGSE4607.RData")

# 2.3) This lasso does exactly the same thing for another dataset. It must be clear that each dataset have to pass this step.
GSE13904_intersect <- list()
for (i in 1:length(GSE13904_sub))
{
  temp <- colnames(GSE13904_sub[[i]])[1]
  colnames(GSE13904_sub[[i]])[1] <- dictionary[temp,3]
  rm(temp)
  if (!is.na(colnames(GSE13904_sub[[i]][1])))
  {
    GSE13904_intersect[[i]] <- GSE13904_sub[[i]]
  }
}
# 2.4) The steps below (save and load) are checkpoints. If is necessary to re-do the analysis, parts of the
# data is stored in RData and can be loaded to speed up the re-analysis.
save(GSE13904_intersect,dictionary, GSE13904_sub, GSE13904, file="~/regulatedGenes-subsetGSE13904.RData")
# 2.5) This step is just a sanity check, to see if there is any artifact left.
a <- c()
for(i in 1:length(GSE13904_intersect)){a[i] <- colnames(GSE13904_intersect[[i]])}
b <- c()
for(i in 1:length(GSE4607_intersect)){b[i] <- colnames(GSE4607_intersect[[i]])}
length(a)==length(b)

# Chunk 3: Dataset preparation. ####
# This lasso formats the data into several forms, described below 
# percentages: a named vector with the percentage of shared genes between regulons in the two networks;
# commonGenes: a list of dataframes, each one containing only the common up (or down) regulated, shared among the two networks;
# This chunk of nested loops does several things, that may look complex, but it is not: it runs for all TFs (i.e., length of the sub.tnet)
# and gets a vector of the unique genes between a given master regulator in two networks. Aside of the gene list,
# the script will create two columns (refering to the datasets). This will help to check wether a gene is regulated (1)
# or no (0). If in both rows they are regulated by the same master regulator (i.e., sum of rows = 2), it takes
# that gene as regulated by both master regulators. This goes for every master regulator within the sub.tnet
# and the results are added to a list (checkList). Simultaneously, the script also calculates the percentage
# of regulated genes between a given master regulator. After that, the script will discard the not commonly
# regulated genes (i.e., row sum =/= 2) and save the rest in a data.frame, called by the master regulator symbol, 
# and further stored in a list, named commonGenes. At this point, we will have a data.frame for every master regulator,
# with its regulated genes (all named after its gene symbol) with the regulation value associated to Probe ID and ENTREZ ID
# (along with gene symbol). As you may ask, the regulation is stored as the regulation value of GSE4607. 
# There is no easy way to get out of it. You can do a harmonic mean of this value. There is no consensus yet.
checkList <- list()
percentages <- c()
TFs <- c()
commonGenes <- list()
allGenes <- data.frame()
for (i in 1:length(GSE13904_intersect))
{
  temp1 <- GSE13904_intersect[[i]]
  temp2 <- GSE4607_intersect[[i]]
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
  colnames(allGenes)[1] <- colnames(GSE13904_intersect[[i]])
  checkList[[i]] <- allGenes
  TFs[i] <- colnames(GSE13904_intersect[[i]])
  percentages[i] <- (sum(rowSums(checkList[[i]][,c(2,3)])==2)/nrow(allGenes))*100
  names(percentages)[i] <- TFs[i]
  commonGenes[[i]] <- as.character(allGenes[which(rowSums((checkList[[i]][,c(2,3)]))==2),1])
  names(commonGenes) <- TFs
  a <- as.data.frame(commonGenes[[i]])
  a[,1] <- as.character(a[,1])
  colnames(a)[1] <- names(commonGenes[i])
  a <- merge(a, gpl_new, by.x = colnames(a)[1], by.y = "PROBEID")
  a <- merge(a, GSE4607_intersect[[i]], by.x = colnames(a)[1], by.y = "row.names")
  colnames(a)[4] <- "Regulation"
  commonGenes[[i]] <- a
}

# 3.1) In order to make all the dataset tidy and reproducible, we are going to treat it
# using dplyr. The commands below just change column classes of all data.frames inside the commonGenes list.
library(dplyr)
commonGenes <- lapply(commonGenes, function(df) mutate_at(df, vars("ENTREZ"), as.numeric))
commonGenes <- lapply(commonGenes, function(df) mutate_at(df, vars("SYMBOL"), as.character))
commonGenes <- lapply(commonGenes, function(df) mutate_at(df, vars("Regulation"), as.numeric))
# 3.2) This command plots the percentages in a bar plot, in a descending order. It is an interesting data,
# but not all interesting findings are useful. Think about it.
barplot(sort(percentages, decreasing = T), ylim = c(0,40), las=2)
# 3.3) The steps below (save and load) are checkpoints. If is necessary to re-do the analysis, parts of the
# data is stored in RData and can be loaded to speed up the re-analysis.
save(percentages, TFs, master_regulators, GSE13904_intersect, GSE4607_intersect, commonGenes, dictionary, file="~/chunk3_commonGenes_allDatasets.RData")


# Chunk 4: Check within other diseases signatures for similarities in regulon activity. ####
# In this chunk, the master regulators obtained with the original signature and datasets are compared with 
# master regulators obtained when interrogating the sepsis network with other genetic signatures, either from
# other inflammatory diseases or even from sepsis (but different tissues). This have an impact on the final number of master regulators.
# If you take into account that sepsis is an inflammatory disease, interrogate the sepsis network with generic inflammatory signatures
# would result in generic inflammatory master regulators. So, if they are present in the original set of master regulators, 
# ***THEORETICALLY*** they cannot be assigned to be specific SEPSIS MASTER REGULATORS.
# 4.1) Load the data.table package just because it is more efficient to read tabular files.
library(data.table)
# 4.2) Creates a list to store genes there are NOT in the below datasets.
outMRs <- list()
# 4.2.1) multiple sclerosis - whole blood
gse4607_sig_gse21942 <- fread("/home/raffael/GDrive/Doutorado/Sepsis R/Paper/SIGNATURE_GSE21942/GSE4607_SIG_GSE21942.txt")
outMRs[["gse4607_sig_gse21942"]] <- master_regulators[which(!master_regulators%in%gse4607_sig_gse21942$Regulon)]
# 4.2.2) rheumatoid arthritis - synovial tissue
gse4607_sig_gse48780 <- fread("/home/raffael/GDrive/Doutorado/Sepsis R/verify_step2/mra_candidates_GSE48780.txt",drop = 1)
outMRs[["gse4607_sig_gse48780"]] <- master_regulators[which(!master_regulators%in%gse4607_sig_gse48780$Regulon)]
# 4.2.3) sepsis (network 2) - sepsis signature (skeletal muscle)
gse13904_sig_gse13205 <- fread("~/GDrive/Doutorado/Sepsis R/Paper/SIGNATURE_GSE13205/GSE13904/GSE13904_SIG_GSE13205.txt")
outMRs[["gse13904_sig_gse13205"]] <- master_regulators[which(!master_regulators%in%gse13904_sig_gse13205$Regulon)]
# 4.2.4) sepsis (network 1) - sepsis signature (skeletal muscle)
gse4607_sig_gse13205 <- fread("~/GDrive/Doutorado/Sepsis R/Paper/SIGNATURE_GSE13205/GSE4607/GSE4607_SIG_GSE13205.txt")
outMRs[["gse4607_sig_gse13205"]] <- master_regulators[which(!master_regulators%in%gse4607_sig_gse13205$Regulon)]
# 4.3) Subset the entire list of master regulators, obtaining the final list. Notice that we came from a list > 200 of MRs
# and finished with (theoretically) only a few specific sepsis master regulators.
master_regulators_subset <- Reduce(intersect, outMRs)
# 4.4) This is the same code that generates the barplot as in 3.2, but this time only with the MRs in the intersection. Cool, right?
barplot(sort(percentages[master_regulators_subset], decreasing = T), ylim = c(0,40), las=2)
# 4.5) The steps below (save and load) are checkpoints. If is necessary to re-do the analysis, parts of the
# data is stored in RData and can be loaded to speed up the re-analysis. And if you came this far, YES, I COPY-PASTE THIS. Sue me.
save(commonGenes, master_regulators_subset, file = "~/commonRegulons.RData")

# Chunk 5: Functional Enrichment ####
# This chunk gets the differentially expressed genes for MRs and identify the metabolic pathways and other
# functional enrichment gibberish.

# 5.1) Starting with loading the package. This tool is also available online at http://amp.pharm.mssm.edu/Enrichr/.
library(enrichR)
# 5.2) Sets the working directory to store all output files separately. This is clearly optional, but is good for organization.
# This directory must be created BEFORE run this command, otherwise won't work.
setwd("~/Funcional_Enrichment_GO-BP2018/")
# 5.3) List the Dbs in Enricher is useful only if you don't know which database will use, if already know...
dbs <- listEnrichrDbs()
# 5.3.1) ...just parse it to a variable (or use it as character in the main function) 
dbs <- "GO_Biological_Process_2018"
# 5.4) Creates a temporary subset of the commonGenes list, which will be applied to the for loop below. 
temp_subset <- commonGenes[master_regulators_subset]
# 5.5) This loop simply gets all genes from each data.frame, submits to the enrichR server (which does de functional enrichment),
# then stores the result to a flat tabular file, for each regulon.
for(i in 1:length(temp_subset)){
  enrichment <- enrichr(c(temp_subset[[i]][[3]]), dbs)
  gobp <- enrichment[[1]] %>% filter(Adjusted.P.value < 0.05)
  write.table(gobp, names(temp_subset[i]),quote = F, sep = "\t", row.names = F, col.names = T)
}
