# Short script to clean and assemble in one data.frame all count tables for N samples in a RNA-seq project.
# Requirements: This script was adapted for featureCounts outputs. According to the manual, the read count 
# table includes annotation columns (‘Geneid’, ‘Chr’, ‘Start’, ‘End’, ‘Strand’and ‘Length’) 
# and data columns (eg. read counts for genes for each library). 

# Load required packages for file read and string manipulation
library(data.table)
library(stringr)

# First we will load metadata about the experiment.
metadata <- read.csv("~/GDrive/Doutorado/Sepsis R/RNA-seq/SraRunTable.txt")
# Then list all files in the directory. The loop (below in code) will iterate through this
files <- list.files("~/gene_counts/",pattern = "txt")
# In this case we will remove the first entry, which will be read manually for object start (counts).
files <- files[-1]
# First reads the object as a data.frame, because f*ck data.table format.
temp <- as.data.frame(fread("/home/raffael/gene_counts/SRR9642740.fastqAligned.sortedByCoord.out.bam_featureCounts.txt"))
# We only want the count column, which in this case, is the 7th column.
counts <- as.data.frame(temp[,7])
# Set the gene IDs as rownames
rownames(counts) <- temp$Geneid
# Change the column name to match the object we just loaded
colnames(counts)[1] <- "SRR9642740"

# Now is time to enter the directory containing the count tables and start the loop.
setwd("~/gene_counts")
# The loop does the following:
# 1) Extract the numbers that comes after the "SRR" in the file name;
# 2) Paste the numbers together;
# 3) Concatenates "SRR" with the number in the temp_name variable then
# 4) prints it to see if everything is ok;
# 5) Loads the file as the iterating variable "i" (as a data.frame, because f*ck data.tables);
# 6) Checks if the row names are the same between the first object (counts) and the current object (temp)
# 7) Gets the count columns (as data.frame... meh, don't bother to curse data.tables again);
# 8) Appends to the "count" object the "temp" object, containing the count table for the current object
# 9) Assigns the column name to the temp_name of the current object.
for(i in files)
{
  temp_name <- unlist(str_extract_all(i, "[^SRR]"))
  temp_name <- paste(temp_name[1:7], collapse = "")
  temp_name <- str_c("SRR", temp_name)
  print(temp_name)
  temp <- as.data.frame(fread(i))
  print(setequal(rownames(counts),temp$Geneid))
  temp <- as.data.frame(temp[,7])
  counts <- data.frame(temp_name = temp, counts)
  colnames(counts)[1] <- temp_name
  
}
# Invert the column order of counts data.frame. Somehow I managed to append  
# the columns backwards. 
counts <- counts[,order(ncol(counts):1)]
# Then we will make some sanity checks: 
# Test wether the run name order is the same in our 
# inverted data.frame and the original order from metadata; then, test if when 
# we add the names, we still get the same objects. May sound stupid, but is always
# good to be sure.
colnames(counts)==metadata$Run
x <- as.character(metadata$Run)
names(x) <- metadata$patient_group
y <- setNames(colnames(counts), metadata$patient_group)
x==y
names(x)==names(y)
# This step does the dictionary for all ENSEMBL (blergh) IDs to Entrez and gene symbol
# by using biomaRt default pipeline.
library(biomaRt)
ensembl <- useMart("ENSEMBL_MART_ENSEMBL")
ensembl <-  useDataset("hsapiens_gene_ensembl", mart=ensembl)
annot <- getBM(attributes=c("ensembl_gene_id", "entrezgene_id", "hgnc_symbol"), 
               filters = "ensembl_gene_id", 
               values = rownames(counts), 
               mart = ensembl)

# Then we save important objects to move on with edgeR pipeline.
save(counts, annot, metadata,file = "~/gene_counts/counts.RData")

# First, load the edgeR package.
library(edgeR)
# Now make the counts object as a DGEList object, which is a specific list that edgeR works with. 
# We also pass the patient groups as the "group" parameter.
dge_counts <- DGEList(counts = counts, group = metadata$patient_group)
# Then we will filter genes with low counts across all libraries.
keep <- filterByExpr(dge_counts)
# And further keep only the filtered genes, discarding low count ones. This step is necessary
# in order to avoid interference.
dge_counts <- dge_counts[keep, ,keep.lib.sizes=F]
# Next, we normalize the RNA counts in order to remove the RNA composition effect
# in order to avoid false positive down-regulated genes.
dge_counts <- calcNormFactors(dge_counts)
# For estimating differentially expressed genes, we must first provide an
# design matrix for edgeR fit the general linear model (GLM) to the counts.
Groups <- metadata$patient_group
design <- model.matrix(~0+Groups)
# Transform from raw counts to log2-counts per million (CPM).
dge_logcpm <- cpm(dge_counts, log = T)
# Update the .RData object with the dge_counts and the log CPM transformed counts.
save(counts, annot, metadata, dge_counts, dge_logcpm, design, file = "~/gene_counts/counts.RData")
# This step estimates the dispersion after data normalization, based on the design matrix.
dge_counts <- estimateDisp(dge_counts, design)
# By using the glmQLFit function, given raw counts, dispersions and the design matrix,
# glmQLFit() fits the negative binomial GLM for each tag and produces an object of class DGEGLM with some new components.
fit <- glmQLFit(dge_counts, design)
# Get the differentially expressed genes between conditions.
coef_test <- glmQLFTest(fit, coef = 2)
# Generate a table with differentially expressed genes and its associated statistics.
topper <- topTags(coef_test, n = 5000)
# Gets only the DEG table. 5000 is an arbitrary unit. 
DEGs <- topper$table
# Let's make an dictionary only for the differentially expressed genes.
ensembl <- useMart("ENSEMBL_MART_ENSEMBL")
ensembl <-  useDataset("hsapiens_gene_ensembl", mart=ensembl)
annot_DEGs <- getBM(attributes=c("ensembl_gene_id", "entrezgene_id", "hgnc_symbol"), 
               filters = "ensembl_gene_id", 
               values = rownames(DEGs), 
               mart = ensembl)

# This is a custom function which uses the package FactoMineR to generate principal
# component analysis. Can be modified in the future. It takes as argument the 
# logCPM (or any other expression unit matrix) and a factor list to pass as colors and
# labels.

my_pca <- function(logcpm, contrasts)
{
  library(FactoMineR)
  texp <- as.data.frame(t(logcpm))
  tPca <- PCA(texp, graph = F)
  scores <- tPca$ind$coord
  lev<-levels(contrasts)
  plot(scores[,1], scores[,2], 
       xlab="PCA 1", ylab="PCA 2",type="p", pch=19,
       col=1:length(lev), cex=1.0, 
       xlim=c(min(scores[,1])*1.1, max(scores[,1])*1.1),
       ylim=c(min(scores[,2]*1.1), max(scores[,2])*1.1))
  text(scores[,1]-7,scores[,2]-7, contrasts,cex=0.7)
}
