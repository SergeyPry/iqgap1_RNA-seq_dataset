if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
# orthogene is only available on Bioconductor>=3.14
if(BiocManager::version()<"3.14") BiocManager::install(version = "3.14")

#BiocManager::install("orthogene")

library(orthogene)

setwd("c:/00_Research_projects/11_RNA-seq/iqgap1_RNAseq_dataset/STAR-RSEM_data/GSEA")

method <- "gprofiler"

# full: iqgap1_53hpf_star_rsem_all-genes.csv
genes_df <- read.csv("iqgap1_53hpf_star_rsem_all-genes.csv", header = TRUE)

# IDs mapping - symbols
symbols <- genes_df$LLgeneSymbol

# remove duplicated symbols
non_duplicated_symb <- which(duplicated(symbols) == FALSE)
unique_symbols <- symbols[non_duplicated_symb]

# IDs mapping - EnsemblIDs

ensembl <- genes_df$EnsemblID

# remove duplicated symbols
non_duplicated_ens <- which(duplicated(ensembl) == FALSE)
unique_ensembl <- ensembl[non_duplicated_ens]


options(scipen=1)

mapped_ensembl <- orthogene::convert_orthologs(gene_df = unique_ensembl,
                                               gene_output = "columns",
                                        input_species = "zebrafish",
                                        output_species = "human",
                                        non121_strategy = "keep_both_species",
                                        method = method)


genes_df <- merge(genes_df, mapped_ensembl, by.x = "EnsemblID", by.y = "input_gene" )

new_colnames <- c("ortholog_gene", "EnsemblID", "LLgeneSymbol", "logFC", "Fold_change", "PValue", "FDR", "wt_1", "wt_2", "wt_3", "wt_4", "iqgap1_mut_1", "iqgap1_mut_2", "iqgap1_mut_3", "iqgap1_mut_4")

genes_df <-  genes_df[, new_colnames]

genes_df <- genes_df |> 
  arrange(desc(logFC))

# 2-fold: MZ-osgep-G177A_genes_edgeR_4_GSEA-v2.csv
# full: MZ-osgep-G177A_genes_edgeR_GSEA_full.csv

write.csv(genes_df, "iqgap1_mut_vs_wt_53hpf_genes_edgeR_GSEA_full.csv", quote = FALSE)
  
########################################################
