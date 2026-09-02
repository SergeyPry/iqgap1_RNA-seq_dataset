# https://bioinformatics.ccr.cancer.gov/docs/btep-coding-club/CC2025/GSEA/GSEA/


# BiocManager::install("clusterProfiler", version = "3.8")
# BiocManager::install("pathview")
# BiocManager::install("enrichplot")


library(clusterProfiler) #Bioconductor
library(tidyverse) #CRAN
library(org.Hs.eg.db) #Bioconductor
library(msigdbr)


# STEPS:

# 1. do orthogene for the full dataset and save the new file
# 2. read in this new file
# 3. run the script below.


setwd("c:/00_Research_projects/11_RNA-seq/iqgap1_RNAseq_dataset/STAR-RSEM_data/GSEA")

# full: iqgap1_mut_vs_wt_53hpf_genes_edgeR_GSEA_full.csv

genes_df <- read.csv("iqgap1_mut_vs_wt_53hpf_genes_edgeR_GSEA_full.csv", header = TRUE)

hallmark<- msigdbr(collection = "H")  %>% 
dplyr::select(gs_name, gene_symbol)

head(hallmark)

h2 <- msigdbr(collection = "H") %>% 
  dplyr::select(gs_name, gs_description)

ranked <- genes_df$logFC
names(ranked) <- genes_df$ortholog_gene


ranked2 <- sort(tapply(ranked, names(ranked), mean), decreasing = TRUE)
new_names <- names(ranked2)
new_ranked <- as.vector(ranked2)
names(new_ranked) <-  new_names



#set seed
set.seed(123)
#run GSEA
eh <- GSEA(new_ranked, TERM2GENE = hallmark, TERM2NAME = h2, seed=TRUE)

eh_df <-  as.data.frame(eh)

eh_df$negLogPvalue <- -log(eh_df$p.adjust, 10)

write.csv(eh_df, "GSEA_clusterProfiler_results_iqgap1-mut_53hpf.csv")


eh@result$Description <- eh_df$Description

enrichplot::gseaplot2(eh, geneSetID = c(1,2,3, 4), base_size = 14)
ggsave("GSEA_1to4.png", dpi = 300)


enrichplot::gseaplot2(eh, geneSetID = c(5,6,7,8), base_size = 14)
ggsave("GSEA_5to8.png", dpi = 300)


######################################## GSEA plots ############################

# simple bar plot

eh_df$ID <- str_replace(eh_df$ID, "HALLMARK_", "")

ggplot(eh_df, aes(negLogPvalue, fct_reorder(ID, negLogPvalue),
                                   fill= NES)) +
  geom_col() +
  geom_vline(xintercept=1.3, linetype="dashed", color="blue",size=1)+
  scale_fill_gradient(high = "darkorange", low = "lightblue") +
  scale_x_continuous(expand =c(0,0), breaks = c(0,1,2,3,4,5,6))+
  theme_bw() + 
  xlab("-log(P-value)") +
  ylab("Gene Set") +
  ggtitle("GSEA (Hallmark)")+
  theme(axis.text.x = element_text( size = 16),
        axis.title.x = element_text(size = 18),
        axis.text.y = element_text(size = 18),
        axis.title.y = element_blank(),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 14),
        panel.grid.major = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.major.x = element_line(colour = "grey70", linewidth = 0.1)
  ) 


ggsave("HALLMARK_sets-iqgap1-mut.png", dpi = 300)
  




