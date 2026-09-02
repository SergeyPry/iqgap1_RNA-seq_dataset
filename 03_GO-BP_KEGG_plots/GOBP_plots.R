library(readr)
library(dplyr)
library("ggsci")
library("ggplot2")
library("gridExtra")

library(ggrepel)
library(forcats)
library(stringr)

setwd("c:/00_Research_projects/11_RNA-seq/iqgap1_RNAseq_dataset/STAR-RSEM_data/GO-BP_KEGG_plots")

GOBP_up <- read.csv("GOBP-UP.tsv", sep = '\t')

GOBP_up$negLogPValue <- -log(GOBP_up$PValue, base = 10)
head(GOBP_up)

GOBP_up <- GOBP_up %>% arrange(negLogPValue)
values <- as.vector(GOBP_up$Term)
GOBP_up$Term <- factor(GOBP_up$Term, levels = values)
##############
GOBP_down <- read.csv("GOBP-DOWN.tsv", sep = '\t')

GOBP_down$negLogPValue <- -log(GOBP_down$PValue, base = 10)
head(GOBP_down)

GOBP_down <- GOBP_down %>% arrange(negLogPValue)
values <- as.vector(GOBP_down$Term)
GOBP_down$Term <- factor(GOBP_down$Term, levels = values)


GOBP_full <- rbind(GOBP_up, GOBP_down)
GOBP_full$Direction <- c(rep("UP", nrow(GOBP_up)), rep("DOWN", nrow(GOBP_down)) )

GOBP_full$Direction <- factor(GOBP_full$Direction, levels = c("DOWN", "UP"))
GOBP_full <-GOBP_full %>% arrange(Direction, negLogPValue)
values <- as.vector(GOBP_full$Term)

GOBP_full$Term <- factor(GOBP_full$Term, levels = values)



############### Full dataset ###############################

ggplot(data = GOBP_full, aes(x = Term, y = negLogPValue, fill = Direction)) +
  geom_col() +
  scale_fill_manual(values=c(  "#56B4E9", "red")) +
  ylab("Negative Log10 of P-value") +
  coord_flip()+
  theme_classic() +
  theme(axis.text.x = element_text( size = 14),
        axis.title.x = element_text(size = 14),
        axis.text.y = element_text(size = 15),
        axis.title.y = element_blank(),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 14),
        panel.grid.major = element_blank(),
        panel.grid.major.x = element_blank()
  )


ggsave("GOBP_pathways_full.png", dpi = 300)
