library(readr)
library(stringr)

# set the working folder to the current folder where this script is stored
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# get all files
base_dir <- getwd()
results_genes <- list()
dir_files <- list.files(base_dir)

# get gene files
gene_files <- dir_files[endsWith(dir_files, "2nd_pass_.genes.results")]

for(fname in gene_files){
  # read data from a file
  data <- read_tsv(fname, col_names = TRUE)
  
  # make sample name
  name_parts <- unlist(strsplit(fname,  '_2nd_'))
  sample_name <- name_parts[1]
  
  # store the data
  results_genes[[ sample_name ]] = data$expected_count
}


results_genes_df <- as.data.frame(results_genes)
rownames(results_genes_df) <- data$gene_id

head(results_genes_df)

colnames(results_genes_df) 

results_genes_df <- results_genes_df[, c("wt_1", "wt_2", "wt_3", "wt_4", "iqgap1_mut_1", "iqgap1_mut_2", "iqgap1_mut_3", "iqgap1_mut_4")]

head(results_genes_df)


write.csv(results_genes_df, "iqgap1_53hpf_star_rsem_counts.csv", quote = FALSE)
  
