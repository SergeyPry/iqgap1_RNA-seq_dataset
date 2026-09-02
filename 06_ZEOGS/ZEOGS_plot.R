library(readr)
library(dplyr)
library("ggsci")
library("ggplot2")
library("gridExtra")

library(ggrepel)
library(forcats)
library(stringr)

setwd("c:/00_Research_projects/11_RNA-seq/iqgap1_RNAseq_dataset/STAR-RSEM_data/ZEOGS")

ZEOGS_df <- read.csv("ZEOGS_selected_df.tsv", sep = '\t')
head(ZEOGS_df)

ZEOGS_df$negLogPValue <- -log(ZEOGS_df$Pvalue, base = 10)


ZEOGS_df <-ZEOGS_df %>% arrange(negLogPValue)
values <- as.vector(ZEOGS_df$Term)
ZEOGS_df$Term <- factor(ZEOGS_df$Term, levels = values)


ggplot(data = ZEOGS_df, aes(x = Term, y = negLogPValue, fill = Term)) +
  geom_col() +
  ylab("Negative Log10 of P-value") +
  coord_flip()+
  theme_classic() +
  theme(axis.text.x = element_text( size = 16),
        axis.title.x = element_text(size = 16),
        axis.text.y = element_text(size = 17),
        axis.title.y = element_blank(),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 14),
        panel.grid.major = element_blank(),
        panel.grid.major.x = element_blank(),
        legend.position = "none"
  )

ggsave("ZEOGS_pathways.png", dpi = 300)



