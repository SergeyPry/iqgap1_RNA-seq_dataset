# install packages if needed

# BiocManager::install("DESeq2", force = TRUE)
# BiocManager::install("AnnotationDbi")
# install.packages("htmlwidgets")
# install.packages("XML")
# install.packages("DBI")
# install.packages("bit64")
# install.packages("blob")
# install.packages("AnnotationDbi")


## load libraries
library(DESeq2)
library(ggplot2)
library(stringr)
library("IHW")

 
# set the working folder to the current folder where this script is stored
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

 
# ###############
# # Functions
# ###############
 
# # function to read and do quick processing of count data

readCountData <- function(filename, sampTab){
  count_data <- read.csv(filename, sep = ",", header = TRUE)
  count_data <- count_data[, 1:ncol(count_data)]
  rownames(count_data) <- count_data$gene
  count_data <- count_data[,-1]
  colnames(count_data) <- rownames(sampTab)

  return(count_data)
}



########################## Data preparation #################################

# # load and process the sample table - custom for each table
sampleTable <- read.csv("sample_table.txt", sep = "\t")
sampleTable$genotype <- factor(sampleTable$genotype, levels = c("wildtype", "iqgap1 mutant"))
rownames(sampleTable) <- sampleTable$sampleName
sampleTable

###################################
# iqgap1 mutants vs wildtype RNA-seq
###################################

# # load the raw count data
count_data <- readCountData("iqgap1_53hpf_star_rsem_counts.csv", sampleTable)
str(count_data)

head(count_data)

saved_rn <- rownames(count_data)

count_data <- as.data.frame(sapply(count_data, ceiling))
rownames(count_data) <- saved_rn

# # 1. Read in the count data into a DESeq2 object
dds <- DESeqDataSetFromMatrix(countData = count_data, colData = sampleTable, design = ~ genotype)

# # Pre-filtering the data
# # need to consider if the better way to count the number  of 0s or sum
# 
# # check the number of rows before
nrow(dds)
# 36351

dds <- dds[rowSums(counts(dds) == 0) < 3, ]
dds <- dds[rowSums(counts(dds)) > 16, ]

nrow(dds)
# 26967

# # 2. Perform vst transform to scale and stabilize variance 
vsd <- vst(dds)
 
# # 3.1 Plot PCA
plotPCA(vsd, "genotype") +  theme(axis.title=element_text(size=14,face="bold"),
                                  axis.text = element_text(size=14),
                                  legend.title = element_text(size=13, face = "bold"),
                                  legend.text = element_text(size=12))


ggsave("iqgap1_STAR-RSEM_PCA.png", dpi = 600)

################## Heatmap of all genes #######################################

complete_results_table <- read.csv("iqgap1_53hpf_star_rsem_DEG_results.csv")
row.names(complete_results_table) <- complete_results_table$LLgeneID

nrow(complete_results_table[complete_results_table$logFC > 0,])
# [1] 761

nrow(complete_results_table[complete_results_table$logFC < 0,])
# [1] 532

# Heatmaps
library(genefilter)
library(pheatmap)

rownames(vsd) <- complete_results_table$LLgeneID

# A heatmap of ALL regulated genes
vsd_all <- vsd[complete_results_table$LLgeneID,]
anno <- as.data.frame(colData(vsd_all)[, c("sampleName","genotype")])
colnames(vsd_all)

count_mat <- complete_results_table[, c("wt_1", "wt_2", "wt_3", "wt_4", "iqgap1_mut_1", "iqgap1_mut_2", "iqgap1_mut_3","iqgap1_mut_4")]
rownames(count_mat) <- complete_results_table$LLgeneID
colnames(count_mat) <- c("wt_1", "wt_2", "wt_3", "wt_4", "iqgap1_mut_1", "iqgap1_mut_2", "iqgap1_mut_3","iqgap1_mut_4")

mat_all2 <- as.matrix(count_mat)

pheatmap(mat_all2,
         method="complete", 
         annotation_col = anno, show_rownames = F, 
         annotation_legend = TRUE, legend=T, 
         cluster_cols=FALSE, cluster_rows = TRUE, scale="row", 
         treeheight_row = 0, treeheight_col = 0,
         color=colorRampPalette(c("navy", "white", "red"))(50))



