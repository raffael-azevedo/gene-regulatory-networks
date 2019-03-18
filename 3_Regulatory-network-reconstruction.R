# This script is set to work with the RTN bioconductor package, in order to conduct
# a regulatory transcription network reconstruction. For more information, please
# read the package vignette at Bioconductor and the README file. 
# This case study recquires the objects obtained from previous analysis (i.e., signature) and
# a dataset for network inference. For this last, it was used septic samples of patients from GSE4607.


# Chunk 1: LOAD PACKAGES AND RDATA #### 
# 1.1) The first step is to load the RData object which contains the objects created in the script 1.
load("~/GSE26378_signature.RData")
# 1.2) Load recquired libraries
library(affy) # from Bioconductor
library(RTN) # from Bioconductor
library(data.table) # from CRAN (part of tydiverse)

# 1.3) This command loads the raw data from septic patients of GSE4607 study and normalizes it using RMA.
# For more information refer to 1_Genetic-signature.R script on GitLab.
data <- ReadAffy()
exp <- rma(data)
GSE4607 <- exprs(exp)
# 1.3.1) The steps below (save and load) are checkpoints. If is necessary to re-do the analysis, parts of the
# data is stored in RData and can be loaded to speed up the re-analysis.
# save(GSE4607,file="chunk1_nrData_GSE4607.RData")
# load("chunk1_nrData_GSE4607.RData")

# Chunk 2: RTN ####
# For running RTN, a list containing 6 objects must be created: gexp, gexpIDs, pheno, phenoIDs, hits and tfs. 
# gexp: contains a NAMED gene expression MATRIX (MANDATORY); 
# gexpIDs: a DATA.FRAME with gexp annotation;
# pheno: a NAMED NUMERIC VECTOR with differential gene expression data
# phenoIDs: a DATA.FRAME with pheno annotation; 
# hits: a CHARACTER VECTOR with genes differentially expressed;
# tfs: a NAMED VECTOR with transcription factors (MANDATORY).

# This pipeline assumes that you already ran a Microarray differential expression pipeline (or a RNA-seq and have an count matrix)
# and have a differential expression matrix, a normalized expression set, contrasts (i.e., disease x control),
# and a dictionary for the probe IDs (platform.db package or GPL file from GEO). 

# In order to check requirements of RTN objects, is advisable to load toy dataset from RTN package
# and compare each object inside dataRTN before proceeding the analysis.

# 2.1) GET GEXPIDs
# 2.1.1 Gets the probe IDs that are differentially expressed, according to the pipeline 
diff_ids <- as.data.frame(rownames(GSE4607))
# 2.1.2 Set the column name to match the one present in the IDs object, in order to make easier to merge
colnames(diff_ids) <- "affy_hg_u133_plus_2"
# 2.1.3 Merge the differentially expressed probeIDs with the dictionary created before
gexpIDs <- merge(x=diff_ids, y=gpl_new, by.x="affy_hg_u133_plus_2",by.y="PROBEID")
# 2.1.4 Set row names for gexpIDs 
rownames(gexpIDs) <- gexpIDs$affy_hg_u133_plus_2
# 2.1.5.1 Gets all coordinates where within this column are blanks
index <- which((gexpIDs$ENTREZID)=="")
# 2.1.5.2 Gets all coordinates where within this column are blanks
index2 <- which((gexpIDs$SYMBOL)=="")
# 2.1.5.3 Changes all blank spaces to NAs in the second column
gexpIDs[index,2] <- NA
# 2.1.5.4 Changes all blank spaces to NAs in the third column
gexpIDs[index2,3] <- NA
# 2.1.5.5 Change column names
colnames(gexpIDs) <- c("PROBEID", "ENTREZ", "SYMBOL")

# 2.2) GET GEXP
# 2.2.1 Transforms DEexp matrix into data.frame for data handling
gexp <- as.data.frame(GSE4607)
# 2.2.2 Makes consistency test to get only entries present in gexpIDs
gexp <- subset(gexp, rownames(gexp)%in%rownames(gexpIDs))
# 2.2.3 Return the format to matrix, because of reasons
gexp <- as.matrix(gexp)

# 2.3) GET PHENO
# The 'pheno' object was already created in the previous step of the analysis, as an data.frame. For more information
# refer to 1_Genetic-signature.R script on GitLab. For this proceeding, it is needed only the probe ids
# and its respective logFC expression value. The final object should be an named vector.
pheno <- GSE26378_pheno

# 2.4) GET HITs
# The 'hits' object was already created in the previous step of the analysis. For more information
# refer to 1_Genetic-signature.R script on GitLab.
hits <- GSE26378_hits

# 2.5) GET PHENO IDs
# The 'phenoIDs' is the dictionary object that was created in the previous step of this analysis.
# For more information refer to 1_Genetic-signature.R script on GitLab.
phenoIDs <- gpl_new

# 2.6) GET TFs
# The 'tfs' object was already created in a previous step of the analysis. For more information
# # refer to 2_Get-TFs.R script on GitLab.
load("~/tfs.RData")
# 
# 2.7. MAKE THE LIST FOR USE IN RTN
# The dataRTN object is a list of the above created objects. As aforementioned, it is important
# to check data structure with the example data provided by the RTN package.
dataRTN <- list(gexp = gexp, gexpIDs = gexpIDs, pheno = pheno, phenoIDs = phenoIDs, hits = hits, tfs = tfs)

# 2.7.1) The steps below (save and load) are checkpoints. If is necessary to re-do the analysis, parts of the
# data is stored in RData and can be loaded to speed up the re-analysis.
# save(dataRTN,file = "chunk2_dataRTN_signatureGSE26378.RData")
# load("chunk2_dataRTN_signatureGSE26378.RData")

# Chunk 3: TNI ANALYSIS ####
# All information regarding commands were extracted from the RTN vignette.

# 3.1. Objects of class TNI provide a series of methods to do transcriptional network inference from high-throughput gene expression data. 
# In this 1st step, the generic function tni.preprocess is used to run several checks on the input data.
rtni <- tni.constructor(expData=dataRTN$gexp, regulatoryElements=dataRTN$tfs, rowAnnotation=dataRTN$gexpIDs)
# 3.2. The tni.permutation function takes the pre-processed TNI object and returns a transcriptional network
# inferred by mutual information (with multiple hypothesis testing corrections).
rtni<-tni.permutation(rtni)
# 3.3. In an additional step, unstable interactions can be removed by bootstrap analysis using the tni.bootstrap function, 
# which creates a consensus bootstrap network (referred here as refnet).
rtni<-tni.bootstrap(rtni)
# 3.4. In the TN each target can be linked to multiple TFs and regulation can occur as a result of both direct (TF-target) 
# and indirect interactions (TF-TF-target). The Data Processing Inequality (DPI) algorithm (Meyer, Lafitte, and Bontempi 2008) 
# is used to remove the weakest interaction in any triangle of two TFs and a common target gene, thus preserving the dominant TF-target pairs, 
# resulting in the filtered transcriptional network (referred here as tnet). The filtered TN has less complexity and highlights the most significant interactions.
rtni<-tni.dpi.filter(rtni)
# 3.5. The RTN package have several visualization formats (rmap, amap, amapDend) and each one is useful to see different characteristics
# of network. In this particular example, the chosen method is the 'amapDend', which returns an igraph object
# that is plotted as an 'tree-and-leaf' graph, an hierarchichal visualization of the transcription regulatory network.
GSE4607_tree<-tni.graph(rtni, tnet="dpi", gtype="amapDend", tfs=dataRTN$tfs)

# Chunk 4: TNA ANALYSIS ####
# TRANSCRIPTIONAL NETWORK ANALYSIS
# 4.1. Create a new TNA object. Objects of class TNA provide a series of methods to do enrichment analysis on transcriptional networks.
# 4.1.2. In this 1st step, the generic function tni2tna.preprocess is used to convert the preprocessed TNI object to TNA, 
# also running several checks on the input data.
rtna <- tni2tna.preprocess(object = rtni, phenotype = dataRTN$pheno, hits = dataRTN$hits, phenoIDs = dataRTN$phenoIDs)
# 4.2. Run MRA analysis pipeline. The tna.mra function takes the TNA object and returns the results of the Master Regulator Analysis (RMA) 
# over a list of regulons from a transcriptional network (with multiple hypothesis testing corrections).
# The MRA computes the overlap between the transcriptional regulatory unities (regulons) 
# and the input signature genes using the hypergeometric distribution (with multiple hypothesis testing corrections). 
rtna <- tna.mra(rtna)
# 4.3. Run overlap analysis pipeline. A simple overlap among all regulons can also be tested using the tna.overlap function.
rtna <- tna.overlap(rtna)
# 4.4. Run GSEA analysis pipeline (THIS STEP IS NOT MANDATORY TO GET THE MASTER REGULATORS). Alternatively, the gene set enrichment analysis (GSEA) can be used to assess if a given transcriptional regulatory unit
# is enriched for genes that are differentially expressed among 2 classes of microarrays (i.e. a differentially expressed phenotype).
# The GSEA uses a rank-based scoring metric in order to test the association between gene sets and the ranked phenotypic difference.
# Here regulons are treated as gene sets, an extension of the GSEA statistics as previously described.
rtna <- tna.gsea1(rtna, stepFilter=FALSE, nPermutations=10000)
# 4.5. Run two-tailed GSEA analysis pipeline (THIS STEP IS NOT MANDATORY TO GET THE MASTER REGULATORS). The two-tailed GSEA tests whether positive or negative targets for a TF are enriched
# at each extreme of a particular response (e.g. differentially expressed genes). The pipeline splits the regulon into a group of activated and a group of repressed genes, based on the Pearson’s correlation,
# and then asks how the two sets are distributed in the ranked list of genes.
rtna <- tna.gsea2(rtna,tfs=dataRTN$tfs, nPermutations=10000)
# 4.6. Get results. All results available in the TNA object can be retrieved using the tna.get function.
# In this particular event, the result of the Master Regulator Analysis is stored in the tna object, a data.frame.
tna <- tna.get(rtna,what="mra")

# Chunk 5: REGULON ACTIVITY ####
# 5.1. The new version of RTN can also check regulon activity within samples. The tni.gsea2 function
# gets the 'activity scores' for each sample.
RegScores <- tni.gsea2(rtni)
# 5.2. The heatmap plots all regulons of interest. In this case, it will make a heatmap for all master regulators
# identified in the TNA analysis.
heatmap(RegScores$dif[,c(tna$Regulon)])
# 5.2.1. The steps below (save and load) are checkpoints. If is necessary to re-do the analysis, parts of the
# data is stored in RData and can be loaded to speed up the re-analysis.
save(dataRTN,rtni,GSE4607_tree,rtna,tna,RegScores, file = "chunk3_TNI-TNA_analysis_GSE4607.RData")
