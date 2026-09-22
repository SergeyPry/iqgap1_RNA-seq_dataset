## Repository for "Signaling scaffold protein Iqgap1 is required for vascular integrity and normal macrophage levels in zebrafish embryos but not for retinal vasculature structure"

A repository with the code, data and figures for this manuscript/paper. Data processing and the analysis contains the following parts: 


0. **Obtaining RNA-seq measurements.**
1. **Assembling RNA-seq data into a count matrix.**
2. **R script for the differential expression quantification by edgeR.**
3. **R script for generating the Principal Component Analysis plot and the global clustering heatmap.**
4. **R scripts for plotting barplots of Gene Ontology Biological Process and KEGG pathway terms.**
5. **R scripts for performing Gene Set Enrichment Analysis (GSEA) and plotting the results.**
6. **R script for plotting enrichment of Ingenuity Pathway Analysis (IPA) terms.**
7. **R script for plotting enrichment of Zebrafish Expression Ontology of Gene Sets (ZEOGS) terms.** 



### 0. Obtaining RNA-seq measurements.
Reads were aligned to the zebrafish GRCz11 reference genome using the 2-step STAR algorithm with the zebrafish transcriptome annotation published by [Lawson et al.](https://pubmed.ncbi.nlm.nih.gov/32831172/) followed by gene and transcript quantification using the RSEM algorithm as we described [before](https://github.com/SergeyPry/CPT_RNA-seq_zebrafish_paper).

### 1. Assembling RNA-seq data into a count matrix.
The code for this step starts with the gene expression measurements obtained in the previous step. The script [process_data.R](https://github.com/SergeyPry/iqgap1_RNA-seq_dataset/blob/main/01_process_data/process_data.R) collects all the files for gene abundance measurements and generates a data frame containing them. We provide the resulting file for download: [iqgap1_53hpf_star_rsem_counts.csv](https://github.com/SergeyPry/iqgap1_RNA-seq_dataset/blob/main/01_process_data/iqgap1_53hpf_star_rsem_counts.csv). 

### 2. R script for the differential expression quantification by edgeR.
The code for this step is included in the [edgeR_script.R](https://github.com/SergeyPry/iqgap1_RNA-seq_dataset/blob/main/02_analyze_data/edgeR_script.R). All the necessary input files for the script as well as the output files are stored in the same folder.

### 3. R script for generating the Principal Component Analysis (PCA) plot and the global clustering heatmap.
The [DESeq2-based script](https://github.com/SergeyPry/iqgap1_RNA-seq_dataset/blob/main/03_plot_PCA-heatmap/deseq2_script_genes.R) in this section takes in the raw count data to produce the [PCA plot](https://github.com/SergeyPry/iqgap1_RNA-seq_dataset/blob/main/03_plot_PCA-heatmap/iqgap1_STAR-RSEM_PCA.png). This script also produces a [global heatmap](https://github.com/SergeyPry/iqgap1_RNA-seq_dataset/blob/main/03_plot_PCA-heatmap/Heatmap_all.pdf) showing the overall pattern of differential gene expression.

### 4. R scripts for plotting barplots of Gene Ontology Biological Process and KEGG pathway terms.
The code for this step is included in the [GOBP_plots.R](https://github.com/SergeyPry/iqgap1_RNA-seq_dataset/blob/main/02_analyze_data/GOBP_plots.R) and [KEGG_plots.R](https://github.com/SergeyPry/iqgap1_RNA-seq_dataset/blob/main/02_analyze_data/KEGG_plots.R). All the necessary input files for the script and its output files are stored in the same folder. These scripts generate plots for Figure 6.

### 5. R scripts for performing Gene Set Enrichment Analysis (GSEA) and plotting the results.
The GSEA was done in two stages. First, the file [iqgap1_53hpf_star_rsem_all-genes.csv](https://github.com/SergeyPry/iqgap1_RNA-seq_dataset/blob/main/05_GSEA_plots/iqgap1_53hpf_star_rsem_all-genes.csv) containing fold changes and other information for all genes in the dataset was provided as the input to the [orthogene_converter.R](https://github.com/SergeyPry/iqgap1_RNA-seq_dataset/blob/main/05_GSEA_plots/orthogene_converter.R), which then produced the file [iqgap1_mut_vs_wt_53hpf_genes_edgeR_GSEA_full.csv])(https://github.com/SergeyPry/iqgap1_RNA-seq_dataset/blob/main/05_GSEA_plots/iqgap1_mut_vs_wt_53hpf_genes_edgeR_GSEA_full.csv) with a column containing human gene names. In the second part of GSEA, the script [gsea_hallmark_full.R](https://github.com/SergeyPry/iqgap1_RNA-seq_dataset/blob/main/05_GSEA_plots/gsea_hallmark_full.R) takes the output of the previous script, Hallmark pathways dataset from the **msigdbr** package to produce the main [GSEA output file](https://github.com/SergeyPry/iqgap1_RNA-seq_dataset/blob/main/05_GSEA_plots/GSEA_clusterProfiler_results_iqgap1-mut_53hpf.csv) as well as a [summary bar graph](https://github.com/SergeyPry/iqgap1_RNA-seq_dataset/blob/main/05_GSEA_plots/HALLMARK_sets-iqgap1-mut.png) and enrichment score plots ([GSEA_1to4.png](https://github.com/SergeyPry/iqgap1_RNA-seq_dataset/blob/main/05_GSEA_plots/GSEA_1to4.png) and [GSEA_5to8.png](https://github.com/SergeyPry/iqgap1_RNA-seq_dataset/blob/main/05_GSEA_plots/GSEA_5to8.png)).

### 6. R script for plotting enrichment of Ingenuity Pathway Analysis (IPA) terms.
The script [filter_split.R](https://github.com/SergeyPry/iqgap1_RNA-seq_dataset/blob/main/06_IPA_plots/filter_split.R) takes the raw IPA input file ([iqgap1_53hpf_IPA_pathways.tsv](https://github.com/SergeyPry/iqgap1_RNA-seq_dataset/blob/main/06_IPA_plots/iqgap1_53hpf_IPA_pathways.tsv)) and filters it based on Z-score being > 1 and negative Log10 of P-value > 1.3. The terms passing the filter were then stored in the output file ([IPA_results_53hpf-iqgap1-mutant.csv](https://github.com/SergeyPry/iqgap1_RNA-seq_dataset/blob/main/06_IPA_plots/IPA_results_53hpf-iqgap1-mutant.csv)). The script also produced an [enrichment plot](https://github.com/SergeyPry/iqgap1_RNA-seq_dataset/blob/main/06_IPA_plots/IPA_pathways-iqgap1-mutant.png) for the selected IPA pathways.

### 7. R script for plotting enrichment of Zebrafish Expression Ontology of Gene Sets (ZEOGS) terms.
The script [ZEOGS_plot.R](https://github.com/SergeyPry/iqgap1_RNA-seq_dataset/blob/main/07_ZEOGS/ZEOGS_plot.R) uses the most significant data from the ZEOGS output ([ZEOGS_selected_df.tsv](https://github.com/SergeyPry/iqgap1_RNA-seq_dataset/blob/main/07_ZEOGS/ZEOGS_selected_df.tsv)) and produces their [enrichment plot file](https://github.com/SergeyPry/iqgap1_RNA-seq_dataset/blob/main/07_ZEOGS/ZEOGS_pathways.png).