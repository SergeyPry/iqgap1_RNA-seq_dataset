# load libraries
library(BiocManager)

#BiocManager::install("DESeq2")
library(edgeR)
library(DESeq2)
library(ggplot2)
library(stringr)

# set the working folder to the current folder where this script is stored
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

 
# ###############
# # Functions
# ###############

## function to read and do quick processing of count data

readCountData <- function(filename, sampTab){
  count_data <- read.csv(filename, sep = ",", header = TRUE)
  rownames(count_data) <- count_data$gene
  count_data <- count_data[,2:ncol(count_data)]
  colnames(count_data) <- rownames(sampleTable)
  
  return(count_data)
}
# 
## a function to generate a table of gene hits with statistics and normalized counts
getDiffGenesTable_edgeR <- function(results, dataset, alpha = 0.05, lfcThreshold = 1) {
  
  # perform filtering of the data
  res <- as.data.frame(results[results$table$FDR < alpha & 
                                 abs(results$table$logFC) >= lfcThreshold, ])
  
  # obtain normalized counts
  norm_counts <- dataset$counts/dataset$samples$norm.factors
  
  # generating a final table of regulated genes
  resSig_full <- res[ order(-res$logFC), ]
  resSig_full$Fold_change <- 2^resSig_full$logFC
  
  norm_counts_sig <- norm_counts[rownames(resSig_full),]
  
  # combine significant hits table and normalized counts matrix
  output <- cbind(resSig_full, norm_counts_sig)
  
  # output the results
  return(output)
  
}

# ##########################
# # End Functions definition
# ##########################

#########################################################################
#                   STAR-RSEM data                                           #
#########################################################################

# # load and process the sample table - custom for each table
sampleTable <- read.csv("sample_table.txt", sep = "\t")
sampleTable$genotype <- factor(sampleTable$genotype, levels = c("wildtype", "iqgap1 mutant"))
rownames(sampleTable) <- sampleTable$sampleName

# # load the raw count data
count_data <- readCountData("iqgap1_53hpf_star_rsem_counts.csv", sampleTable)
head(count_data)
nrow(count_data)
# [1] 36351

# do basic quality controls
count_data <- count_data[rowSums(count_data == 0) < 3, ]
count_data <- count_data[rowSums(count_data) > 16, ]
nrow(count_data)
# [1] 26889

# genes
genetable <- data.frame(gene = rownames(count_data))

# edgeR differential analysis
edge_STAR <- DGEList(counts= count_data, 
                     samples=sampleTable, 
                     genes=genetable)

names(edge_STAR)

design_STAR <- model.matrix(~ genotype, edge_STAR$samples)

edge_STAR <- calcNormFactors(edge_STAR)
edge_STAR<- estimateDisp(edge_STAR, design_STAR)

# getting results
fit <- glmFit(edge_STAR, design_STAR)
lrt <- glmLRT(fit, coef=ncol(design_STAR))
tt.all <- topTags(lrt, n=nrow(edge_STAR), sort.by="none")

# getDiffGenesTable_edgeR <- function(results, dataset, alpha = 0.05, lfcThreshold = 1)
results_STAR <- getDiffGenesTable_edgeR(tt.all, edge_STAR, alpha = 0.05, lfcThreshold = 0.585)

dim(results_STAR)
head(results_STAR)

############# get results for all genes ###############################

results_full <- getDiffGenesTable_edgeR(tt.all, edge_STAR, alpha = 0.999, lfcThreshold = 0)

nrow(results_full)


################################### add gene information ############################

############### set up IDs tables ################

# map of the IDs to other parameters
gtf_df <- read.csv("v4_3_2geneinfo.txt", sep = '\t')

# add an Ensembl ID column without the version
gtf_df$EnsemblID <- str_replace(gtf_df$Ens99geneIDversion, '\\.\\d+', '')

# select the relevant ones
small_gtf <- gtf_df[,c('LLgeneID', 'EnsemblID', 'LLgeneSymbol', 'LLchr',	'LLstart',	'LLend')]



# Add LLgeneID column to enable merging
results_all <- results_STAR
results_all$LLgeneID <- rownames(results_all)
results_all <- merge(results_all, small_gtf, by = 'LLgeneID')

results_all <- results_all[order(results_all$logFC,  decreasing = TRUE),]

head(results_all)

colnames(results_all)

results_all <- results_all[, c("LLgeneID", "EnsemblID", "LLgeneSymbol",  "logFC", "Fold_change", "PValue", "FDR", "wt_1", "wt_2", "wt_3", "wt_4", "iqgap1_mut_1", "iqgap1_mut_2", "iqgap1_mut_3", "iqgap1_mut_4")]
rownames(results_all) <- results_all$LLgeneID

# filter out the LOC and XLOC genes

results_all <- results_all[!grepl("^XLOC", results_all$LLgeneSymbol),]
results_all <- results_all[!grepl("^LOC", results_all$LLgeneSymbol),]
nrow(results_all)

###################################### Output ########################################

# # write the table to an csv file

write.csv(results_all, "iqgap1_53hpf_star_rsem_DEG_results.csv")


########## add gene information to full dataset #######################

results_full$LLgeneID <- rownames(results_full)
results_full <- merge(results_full, small_gtf, by = 'LLgeneID')

results_full <- results_full[order(results_full$logFC,  decreasing = TRUE),]

head(results_full)

colnames(results_full)

results_full <- results_full[, c("LLgeneID", "EnsemblID", "LLgeneSymbol",  "logFC", "Fold_change", "PValue", "FDR", "wt_1", "wt_2", "wt_3", "wt_4", "iqgap1_mut_1", "iqgap1_mut_2", "iqgap1_mut_3", "iqgap1_mut_4")]
rownames(results_full) <- results_full$LLgeneID

# filter out the LOC and XLOC genes

results_full <- results_full[!grepl("^XLOC", results_full$LLgeneSymbol),]
results_full <- results_full[!grepl("^LOC", results_full$LLgeneSymbol),]
nrow(results_full)

###################################### Output ########################################

# # write the table to an csv file
write.csv(results_full, "iqgap1_53hpf_star_rsem_all-genes.csv")







