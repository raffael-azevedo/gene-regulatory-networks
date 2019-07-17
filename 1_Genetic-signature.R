# Script for generating genetic signature from a microarray dataset. This script must be edited to suit specific needs.
# Current version does not support fully automatic process. 
# In a few words, the script below is the standard for 1) reading, 2) normalizing, 3) fitting to linear model and 4) getting
# differentially expressed genes (DEGs).
# This is a case study, using a Affymetrix HGU133 Plus 2.0 platform on a sepsis dataset, obtained from Gene Expression Omnibus (GEO).

# Chunk 1: Read and treat data ####
# This chunk assumes one have already downloaded desired raw data from GEO and stored in a 
# particular directory, decompressed. First, load some required libraries. 

library(affy) # from Bioconductor
library(limma) # from Bioconductor
library(data.table) # from CRAN (part of tidyverse)

# 1.1) Set working directory to dataset directory
setwd("~/GSE26378_RAW/")
# 1.2) Load in .CEL files
data <- ReadAffy()
# 1.3) Performs data normalization using Robust Multiarray Analysis (RMA). See ?rma for further
# details.
exp <- rma(data)
# 1.4) Makes an expression matrix out from the normalized data (i.e., transforms Expression Data
# object in a standard expression matrix, which have genes in rows (observations) and samples
# in columns (variables). Expression is measured as Arbitrary Expression Units).
HSexp <- exprs(exp)
# 1.5) This step does a MDS plot. It is a necessary step to assess sample quality. The full
# quality control is defined in script 0-Quality-control. See more details in ?plotMDS.
plotMDS(HSexp)
# 1.6) This step loads the last GPL file provided by Affymetrix. The GPL contains information
# regarding each probe linking to an gene symbol, entrez, ensembl, and other identifiers.
# This step is necessary to build the dictionary for the platform. This dictionary will be
# crucial to perform future analysis. Note that 16 lines were skipped. This is due to extra info
# present in the file. Please be sure and check the original GPL file to ensure the correct number
# of lines to skip.
gpl570 <- fread("~/GPL570-55999.txt",skip = 16)
# 1.7) Further steps makes the GPL data tidy (but not using any Tidyverse package. Further improve?)
# 1.7.1) Lines below sets the columns of interest (Gene Symbol and ENTREZ identifiers) to temporary
# variables.
b <- gpl570$`Gene Symbol`
c <- gpl570$ENTREZ_GENE_ID
# 1.7.2) As seen in the original file, some lines have more than 1 observation, separated by '//'. 
# The data is reordered to a more clean format, tab delimited.
b <- as.data.frame(gsub("\\/.*","",b))
c <- as.data.frame(gsub("\\/.*","",c))
# 1.7.3) This line creates a data.frame called "dictionary", which holds the affy probe IDs and its
# respectives gene symbols and entrez identifiers. This object will be the 'phenoIDs' object of the RTN pipeline.
dictionary <- cbind.data.frame(gpl570$ID, c,b)
# 1.7.4) Lines below changes the column names to match its data, and assigns the first column (PROBEID)
# to rownames. This step is only necessary because the RTN requirements.
colnames(dictionary) <- c("PROBEID", "ENTREZ","SYMBOL")
rownames(dictionary) <- dictionary[,1]
# 1.7.5) Removes unecessary temporary variables.
rm(b,c,gpl570)

# Chunk 2: Differential Expression Analysis ####
# Herein begins the differential analysis. It uses the limma package to do measurements and statistical
# analysis. There are a few steps to do so, which are listed below.
# 2.0) The "treats" object is either a # text file setting each sample to its condition or a construct
# variable assigning factors to each sample name. In this case, for GSE26378 we have 21 control samples
# and 82 cases. Therefore, the "treats" object will be a factor of each sample, assigning it either to 
# case or control. It is easy when samples are ordered as such (first control samples, then case samples) 
# but it is often usual to encounter samples in a disarranged list (which is this case, by the way).
treats <- as.factor(c(rep("shock", 6), rep("control", 2), "shock", rep("control", 4),
                      rep("shock", 26), "control", rep("shock", 8), "control", rep("shock", 13),
                      "control", "control", rep("shock", 13), rep("control", 5), rep("shock", 12),
                      rep("control", 6), rep("shock", 3)))
# 2.1) This line does the model matrix, which establishes causes versus control. Based on this
# matrix, each comparison can be made and summarized correctly. The base function model.matrix creates a design (or model)
# matrix, e.g., by expanding factors to a set of dummy variables (depending on the contrasts) 
# and expanding interactions similarly. Also we have to modify the column names because this function
# inserts a fuzzy name. 
design <- model.matrix(~0+treats)
colnames(design) <- c("shock", "control")
# 2.2) Below line uses the lmFit() function to fit linear model for each gene given a series of arrays.
# It takes as arguments the expression matrix and the design object.
fit <- lmFit(HSexp,design)
# 2.3) This command establishes the contrasts. In this case, just a case-control was compared. 
# The "levels" arguments is set to the amount of factors defined in the design variable.
contrasts <- makeContrasts(shock-control, levels=design)
# 2.4) This command apply empirical Bayes smoothing to the standard errors, with a cutoff of 0.01.
# According to the function documentation, given a microarray linear model fit, 
# compute moderated t-statistics, moderated F-statistic, and log-odds of differential expression
# by empirical Bayes moderation of the standard errors towards a common value. It also deppends
# on the resulting fit of contrasts, in which Given a linear model fit to microarray data,
# compute estimated coefficients and standard errors for a given set of contrasts.
ct.fit <- eBayes(contrasts.fit(fit, contrasts),0.01)
# 2.5) The decideTests() function which summarize the results of the linear model,
# perform hypothesis tests and adjust thep-values for multiple testing. The choosen method is suitable
# for each analysis, as the p adjustment method.
res.fit<-decideTests(ct.fit,method="separate", adjust.method="BH", p.value=0.01)
# 2.5.1) This lines save as a data.frame an object necessary to originate three new objects:
# "pheno", "phenoIDs", and "hits". These three objects are mandatory to the RTN pipeline.
# This object stores the gene names, their logFC, their pvalue, adjusted pvalue,
# and differential condition (up, down or neutral regulation).
# The GSE26378 was used only as example purpose. Also, the 'pheno' is an object wich is used
# further for the RTN analysis. 
GSE26378 <- data.frame(logFC = ct.fit$coef, 
                             p.value = ct.fit$p.value, degenes = unclass(res.fit), 
                             stringsAsFactors = FALSE)
# 2.6) The 'hits' object is also an object which is used for the RTN analysis. It consists of 
# all genes which had a differential condition diverging from 0. 'hits' object just a character vector
# containing all differentially expressed probes.
GSE26378_hits <- subset(GSE26378, GSE26378$shock...control.2!=0
 & GSE26378$shock...control < -1 | GSE26378$shock...control > 1)
# 2.7) The 'pheno' object contains the expression value (in logFC) of all probes. It is a named vector,
# in which each expression value is named after its probe.
GSE26378_pheno <- GSE26378$shock...control
names(GSE26378_pheno) <- rownames(GSE26378)
# 2.7) The DEexp object is the final product of this analysis. It contains the GSE26378_pheno
# but with gene symbol and entrez IDs. It is not used by RTN, but is generated to save the DEGs.
# 2.7.1) Take the probe ids and assign to a temporary variable.
affy_hg_u133_plus_2 <- rownames(GSE26378)
# 2.7.2) Set the row names to NULL, cleaning the data.
rownames(GSE26378) <- NULL
# 2.7.3) Creates the DEexp object by binding the temp variable and the GSE26378 object,
# which holds the DEGs and its information (p-value, etc).
DEexp <- cbind(affy_hg_u133_plus_2,GSE26378)
# 2.7.4) Adds the information of ENTREZ and HGNC symbol to the dataset.
GSE26378 <- merge(dictionary,DEexp,by.x = "PROBEID",by.y = "affy_hg_u133_plus_2")
# 2.8) This line saves as a RData the dictionary file, pheno, phenoIDs and hit objects, which are needed to 
# RTN analysis, along with the resulting DEG analysis object.
save(GSE26378_hits,GSE26378_pheno, GSE26378,dictionary,file="GSE26378_signature.RData")