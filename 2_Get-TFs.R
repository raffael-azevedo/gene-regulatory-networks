# This part of the RTN analysis is needed to get information regarding transcription factors.
# This is an intermediate step, but needed to construct the RTN object. 

# Chunk 1: Load recquired packages ####
library(biomaRt) # from Bioconductor
library(Fletcher2013b) # from Bioconductor

# The main point is: Fletcher2013b has an object containing a list of TFs, linked to an Illumina probe ID. 
# The objective is to transform the probe IDs from the source to the probe ID of interest (in this case, Affymetrix HG U133 Plus 2.0)

# 1) GET TFs
# The Transcription Factors used for the network inference/construction were obtained from Fletcher2013b package.
# 1.1. Load "miscellaneous" data, containing (among other) the tfs object. The rest may be discarded. 
data("miscellaneous")
# 1.2. The TFs in the example dataset are in concordance with Illumina Human HT 12v4 Technology. In this step
# is used the biomaRt package to translate IDs from Illumina to Affymetrix (or any other platform used), also obtaining HGNC symbol.
tfs <- getBM(filters = "illumina_humanht_12_v4", attributes = c("affy_hg_u133_plus_2", "hgnc_symbol"), values = tfs, mart = ensembl)
# 1.3.1. Create temporary index to check if are duplicated probe IDs
index <- which(duplicated(tfs$affy_hg_u133_plus_2))
# 1.3.2. Create temporaty index to check if there are empty values in symbols column
index2 <- which(tfs$hgnc_symbol=="")
# 1.3.3. Create temporaty index to check if there are empty values in Affymetrix Probe ID column
index3 <- which(tfs$affy_hg_u133_plus_2=="")
# 1.4.1. Remove values based in index
tfs <- tfs[-index,]
# 1.4.2. Remove values based in index2
tfs <- tfs[-index2,]
# 1.4.3. Remove values based in index3
tfs <- tfs[-index3,]
# 1.5. Creates a temporary object to store Affy IDs related to the TFs
temp <- tfs$affy_hg_u133_plus_2
# 1.6. Assign the HGNC Symbols as names to the temp vector
names(temp) <- tfs$hgnc_symbol
# 1.7. Assigns the temp named vector to the tfs vector
tfs <- temp
# 1.8. Removes empty values that couldn't be mapped
tfs <- tfs[-which(names(tfs)=="")]
# 1.9. Removes temporary indexes and variables
rm(temp,index,index2,index3,risksites,randsites,chromlen,ESR1bdsites,FOXA1bdsites,GATA3bdsites,SPDEFbdsites,fimoESR1,fimoFOXA1,fimoGATA3,metaPCNA,consensus)
# 1.10) The steps below (save and load) are checkpoints. If is necessary to re-do the analysis, parts of the
# data is stored in RData and can be loaded to speed up the re-analysis.
save(tfs, "tfs.RData")
load("~/tfs.RData")