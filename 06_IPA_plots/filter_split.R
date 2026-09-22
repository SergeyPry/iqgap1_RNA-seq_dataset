library(tidyverse)

# set the working folder to the current folder where this script is stored
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

all_data <- read.csv("iqgap1_53hpf_IPA_pathways.tsv")
head(all_data)

all_data <-  all_data[all_data$Zscore !="#NUM!",]
nrow(all_data)

all_data$Zscore <- as.numeric(all_data$Zscore)

above2_data <- all_data[all_data$negLogPvalue >= 1.3 & abs(all_data$Zscore) >= 1,]

above2_data$Abs_Zscore <- abs(above2_data$Zscore)

nrow(above2_data)

above2_data <-  above2_data |> 
                arrange(desc(negLogPvalue))


write.csv(above2_data, "IPA_results_53hpf-iqgap1-mutant.csv")

above2_data <- read.csv("IPA_results_53hpf-iqgap1-mutant.csv")

ggplot(data = above2_data, aes(x = reorder(Pathway, negLogPvalue), fill = Zscore )) +
  geom_col( aes(y=negLogPvalue), size=.1, color="black", alpha = 1, width = 0.7) + 
  scale_y_continuous(
    name = "-log(P-value)",
    breaks = c(0,2,4,6,8)) +
  scale_fill_gradient2(
    low = "blue", 
    mid = "white", 
    high = "darkorange", 
    midpoint = 0
  ) +
  coord_flip()+
  theme_classic() +
  theme(axis.text.x = element_text( size = 13),
        axis.title.x = element_text(size = 14),
        axis.text.y = element_text(size = 13),
        axis.title.y = element_blank(),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 14),
        panel.grid.major = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.major.x = element_line(colour = "grey70", linewidth = 0.1)
  )  

ggsave("IPA_pathways-iqgap1-mutant.png", dpi = 300)

