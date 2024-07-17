### Libraries
library(tidytree)
library(DESeq2)
library(circlize)
library(ComplexHeatmap)
library(clusterProfiler)
library(DOSE)
library(enrichplot)
library(org.Mm.eg.db)
library(magrittr)
library(pathview)
library(stats)
library(ggplot2)
library(data.table)

### Read Files ####

setwd("C:/Users/ammas/Documents/Bulk_RNAseq_UChicago_NLRP3")
load.path <- "C:/Users/ammas/Documents/Bulk_RNAseq_UChicago_NLRP3/Raw_Data/"

# Replace "path/to/file" with the actual path to your file
file_path <- paste0(load.path,"raw_combined_featureCounts.txt")

# Read the data into a data frame
feature_counts <- read.table(file_path, header = TRUE, sep = "\t", comment.char = "#")

#### Preparing Data and Metadata ####

# Extracting count data
count_data <- feature_counts[, 7:ncol(feature_counts)]

# Setting gene IDs as row names for count data
rownames(count_data) <- feature_counts$Geneid

# Creating metadata dataframe

metadata <- data.frame(
  sample = c("M1_U_1", "M1_U_2", "M1_ADT_3", "M1_ADT_1", "M1_ADT_2", "M1_U_3",
             "M1_NLRP3_1", "M1_NLRP3_2", "M1_Combo_3", "M1_Combo_1", "M1_Combo_2", "M1_NLRP3_3",
             "M2_U_1", "M2_U_2", "M2_U_3", "M2_ADT_1", "M2_ADT_2", "M2_ADT_3",
             "M2_NLRP3_1", "M2_NLRP3_2", "M2_NLRP3_3", "M2_Combo_1", "M2_Combo_2", "M2_Combo_3"),
  condition = c("M1_U", "M1_U", "M1_ADT", "M1_ADT", "M1_ADT", "M1_U",
                "M1_NLRP3", "M1_NLRP3", "M1_Combo", "M1_Combo", "M1_Combo", "M1_NLRP3",
                "M2_U", "M2_U", "M2_U", "M2_ADT", "M2_ADT", "M2_ADT",
                "M2_NLRP3", "M2_NLRP3", "M2_NLRP3", "M2_Combo", "M2_Combo", "M2_Combo"),
  file.name = c("AP.RS24.23013_S31_L006_R1_001.fastq_subread_results.bam",
                "AP.RS24.23014_S32_L006_R1_001.fastq_subread_results.bam",
                "AP.RS24.23015_S33_L006_R1_001.fastq_subread_results.bam",
                "AP.RS24.23016_S34_L006_R1_001.fastq_subread_results.bam",
                "AP.RS24.23017_S35_L006_R1_001.fastq_subread_results.bam",
                "AP.RS24.23018_S36_L006_R1_001.fastq_subread_results.bam",
                "AP.RS24.23019_S37_L006_R1_001.fastq_subread_results.bam",
                "AP.RS24.23020_S38_L006_R1_001.fastq_subread_results.bam",
                "AP.RS24.23021_S39_L006_R1_001.fastq_subread_results.bam",
                "AP.RS24.23022_S40_L006_R1_001.fastq_subread_results.bam",
                "AP.RS24.23023_S41_L006_R1_001.fastq_subread_results.bam",
                "AP.RS24.23024_S42_L006_R1_001.fastq_subread_results.bam",
                "AP.RS24.23025_S43_L006_R1_001.fastq_subread_results.bam",
                "AP.RS24.23026_S44_L006_R1_001.fastq_subread_results.bam",
                "AP.RS24.23027_S45_L006_R1_001.fastq_subread_results.bam",
                "AP.RS24.23028_S46_L006_R1_001.fastq_subread_results.bam",
                "AP.RS24.23029_S47_L006_R1_001.fastq_subread_results.bam",
                "AP.RS24.23030_S48_L006_R1_001.fastq_subread_results.bam",
                "AP.RS24.23031_S49_L006_R1_001.fastq_subread_results.bam",
                "AP.RS24.23032_S50_L006_R1_001.fastq_subread_results.bam",
                "AP.RS24.23033_S51_L006_R1_001.fastq_subread_results.bam",
                "AP.RS24.23034_S52_L006_R1_001.fastq_subread_results.bam",
                "AP.RS24.23035_S53_L006_R1_001.fastq_subread_results.bam",
                "AP.RS24.23036_S54_L006_R1_001.fastq_subread_results.bam")
)

metadata$cell_type <- ifelse(grepl("M1", metadata$condition), "M1", 
                             ifelse(grepl("M2", metadata$condition), "M2", NA))

# Remove the "M1_" and "M2_" prefixes from the condition column
metadata$condition <- gsub("M1_", "", metadata$condition)
metadata$condition <- gsub("M2_", "", metadata$condition)

# Align the sample names
rownames(metadata) <- metadata$file.name
metadata <- metadata[colnames(count_data), ]  # This reorders the metadata to match count data

#### DeSeqDataSet ####

# Create the DESeq dataset
# Setup DESeqDataSet
dds <- DESeqDataSetFromMatrix(countData = count_data,
                              colData = metadata,
                              design = ~ cell_type + condition + cell_type:condition)

# Running DESeq ####

dds <- DESeq(dds)

# Subset for M1
dds_M1 <- dds[metadata$cell_type == "M1",]
dds_M1 <- DESeqDataSetFromMatrix(countData = count_data[,metadata$cell_type == "M1"],
                                 colData = metadata[metadata$cell_type == "M1",],
                                 design = ~ condition)
dds_M1 <- DESeq(dds_M1)

# Subset for M2
dds_M2 <- dds[metadata$cell_type == "M2",]
dds_M2 <- DESeqDataSetFromMatrix(countData = count_data[,metadata$cell_type == "M2"],
                                 colData = metadata[metadata$cell_type == "M2",],
                                 design = ~ condition)
dds_M2 <- DESeq(dds_M2)

#### Heatmaps for genes of interest ##############
# Load necessary libraries
library(DESeq2)
library(pheatmap)

# Vector of mouse genes
mouse_genes <- c("Ikbkb", "Traf6", "Ifnl2", "Mavs", "Eftud2", "Ddx58", 
                 "Rela", "Irf3", "Il1a", "Ifih1", "Prkra", "Gsdmd", "Id2", 
                 "Tlr2", "Ifnb1", "Ap3b1", "Nlrp3", "Nod2", "Il10rb", 
                 "Tlr9", "Il17ra", "Sod1", "Hrh4", "Ifnar2", "Sdf2l1", 
                 "Vpreb1", "Pnp", "Igll1", "Rnase7", "Siglec15", "Psmb5", 
                 "Cebpe", "Il25", "Cmtm5", "Psme1", "Psme2", "Nfatc4", 
                 "Mx1", "Lif", "Ctsg", "Gzmb", "Tcn2", "Tnfrsf11a", 
                 "Nfkbia", "Tmem173", "Stat1", "Il22ra1", "Bcl2l1", 
                 "Jak1", "Jak2", "Jak3", "Tyk2")
# Subset the metadata for M1 samples
metadata_M1 <- metadata %>% filter(grepl("^M1", sample))

# Extract normalized counts for M1
normalized_counts_M1 <- counts(dds_M1, normalized=TRUE)

# Filter the count data to include only the mouse genes
selected_genes_M1 <- normalized_counts_M1[rownames(normalized_counts_M1) %in% mouse_genes,]

# Create a mapping of file names to conditions
file_to_condition_M1 <- setNames(metadata_M1$condition, metadata_M1$file.name)

# Split file names by condition
file_groups_M1 <- split(names(file_to_condition_M1), file_to_condition_M1)

# Aggregate counts by condition for M1
avg_counts_M1 <- sapply(file_groups_M1, function(files) {
  rowMeans(selected_genes_M1[, files, drop=FALSE])
})

# Subset the metadata for M2 samples
metadata_M2 <- metadata %>% filter(grepl("^M2", sample))

# Extract normalized counts for M2
normalized_counts_M2 <- counts(dds_M2, normalized=TRUE)

# Filter the count data to include only the mouse genes
selected_genes_M2 <- normalized_counts_M2[rownames(normalized_counts_M2) %in% mouse_genes,]

# Create a mapping of file names to conditions
file_to_condition_M2 <- setNames(metadata_M2$condition, metadata_M2$file.name)

# Split file names by condition
file_groups_M2 <- split(names(file_to_condition_M2), file_to_condition_M2)

# Aggregate counts by condition for M2
avg_counts_M2 <- sapply(file_groups_M2, function(files) {
  rowMeans(selected_genes_M2[, files, drop=FALSE])
})

### Heatmaps

# Create and save a high-resolution heatmap for M1
png("Heatmap_M1.png", width = 2500, height = 3200, res = 300)
pheatmap(avg_counts_M1, 
         cluster_rows = FALSE, 
         cluster_cols = T, 
         scale = "row", 
         main = "Heatmap of Selected Genes in M1 Cells Grouped by Condition")
dev.off()

# Create and save a high-resolution heatmap for M2
png("Heatmap_M2.png", width = 2500, height = 3200, res = 300)
pheatmap(avg_counts_M2, 
         cluster_rows = FALSE, 
         cluster_cols = T, 
         scale = "row", 
         main = "Heatmap of Selected Genes in M2 Cells Grouped by Condition")
dev.off()

### Do the same for M1 vs M2 U
# Subset the metadata for M1 samples
# Create a new column combining cell_type and condition
metadata$cell_condition <- paste(metadata$cell_type, metadata$condition, sep = "_")


metadata_U <- metadata

metadata_U$cell_condition <- ifelse(metadata_U$cell_condition %in% c("M1_U", "M2_U"), 
                                  metadata_U$cell_condition, "Other")
# Extract normalized counts
normalized_counts <- counts(dds, normalized=TRUE)

# Filter the count data to include only the mouse genes
selected_genes <- normalized_counts[rownames(normalized_counts) %in% mouse_genes,]

# Create a mapping of file names to conditions
file_to_condition_U <- setNames(metadata_U$cell_condition, metadata_U$file.name)

# Split file names by condition
file_groups_U <- split(names(file_to_condition_U), file_to_condition_U)

# Aggregate counts by condition for M1
avg_counts_U <- sapply(file_groups_U, function(files) {
  rowMeans(selected_genes[, files, drop=FALSE])
})
avg_counts_U <- avg_counts_U[, !colnames(avg_counts_U) %in% "Other"]


# Create and save a high-resolution heatmap for M2
png("Heatmap_U.png", width = 2500, height = 3200, res = 300)
pheatmap(avg_counts_U, 
         cluster_rows = FALSE, 
         cluster_cols = T, 
         scale = "row", 
         main = "Heatmap of Selected Genes in M2 Cells Grouped by Condition")
dev.off()

##########################################################
#### Differential Expression M1 ####

# Define the list of contrasts
contrasts_list <- list(
  c("Combo", "NLRP3"),
  c("Combo", "ADT"),
  c("Combo", "U"),
  c("NLRP3", "ADT"),
  c("NLRP3", "U"),
  c("ADT", "U")
)


# Specify the directory to store results
output_dir <- "C:/Users/axi313/Documents/Bulk_RNAseq_UChicago_NLRP3/DGE/M1"

# Loop through the contrasts and save results without NA
for (contrast in contrasts_list) {
  # Construct the name for the result object
  result_name <- paste("M1", contrast[1], "vs", contrast[2], "DGE", sep=".")
  
  # Perform the differential expression analysis
  result <- results(dds_M1, contrast=c("condition", contrast[1], contrast[2]))
  
  # Remove rows with any NA values
  result_clean <- na.omit(as.data.frame(result))
  
  # Assign the cleaned results to a dynamically named variable
  assign(result_name, result_clean)
  
  # Save the cleaned results to a CSV file
  write.csv(result_clean, 
            file.path(output_dir, paste0(result_name, ".csv")))
}

#### Differential Expression M2 ####

# Define the list of contrasts
contrasts_list <- list(
  c("Combo", "NLRP3"),
  c("Combo", "ADT"),
  c("Combo", "U"),
  c("NLRP3", "ADT"),
  c("NLRP3", "U"),
  c("ADT", "U")
)


# Specify the directory to store results
output_dir <- "C:/Users/axi313/Documents/Bulk_RNAseq_UChicago_NLRP3/DGE/M2"

# Loop through the contrasts and save results without NA
for (contrast in contrasts_list) {
  # Construct the name for the result object
  result_name <- paste("M2", contrast[1], "vs", contrast[2], "DGE", sep=".")
  
  # Perform the differential expression analysis
  result <- results(dds_M2, contrast=c("condition", contrast[1], contrast[2]))
  
  # Remove rows with any NA values
  result_clean <- na.omit(as.data.frame(result))
  
  # Assign the cleaned results to a dynamically named variable
  assign(result_name, result_clean)
  
  # Save the cleaned results to a CSV file
  write.csv(result_clean, 
            file.path(output_dir, paste0(result_name, ".csv")))
}

#### Differential Expression M2 Vs M1####

# Define the directory to store results
output_dir <- "C:/Users/axi313/Documents/Bulk_RNAseq_UChicago_NLRP3/DGE/M2vsM1"

# Get a list of unique conditions available
conditions <- unique(metadata$condition)

# Loop through each condition and save results
for (condition in conditions) {
  # Construct the name for the result object
  result_name <- paste("M2vsM1", condition, "DGE", sep=".")
  
  # Subset the dataset for the current condition
  dds_subset <- dds[metadata$condition == condition,]
  
  # Perform differential expression analysis between M2 and M1 for the current condition
  result <- results(dds_subset, contrast=c("cell_type", "M2", "M1"))
  
  # Remove rows with any NA values
  result_clean <- na.omit(as.data.frame(result))
  
  # Assign the cleaned results to a dynamically named variable
  assign(result_name, result_clean)
  
  # Save the cleaned results to a CSV file
  write.csv(result_clean, 
            file.path(output_dir, paste0(result_name, ".csv")))
}

### Pathway Analysis ####

# Define the function to filter and create a named vector for pathway analysis
create_genelist <- function(dif_genes, log2fc_threshold = 1.5, adj_pvalue_threshold = 0.05) {
  # Filter for significant genes based on adjusted p-value and fold change thresholds
  filtered_genes <- dif_genes[abs(dif_genes$log2FoldChange) > log2fc_threshold & 
                                dif_genes$padj < adj_pvalue_threshold, ]
  
  # Create a named vector of log2 fold changes with gene names as vector namesn
  
  genelist <- setNames(filtered_genes$log2FoldChange, rownames(filtered_genes))
  
  # Remove any NA values
  genelist <- na.omit(genelist)
  
  # Sort the list in decreasing order (required for clusterProfiler)
  genelist <- sort(genelist, decreasing = TRUE)
  
  # Return the final named vector
  return(genelist)
}

create_genelist_Entrez <-function(dif_genes, log2fc_threshold = 1.5, adj_pvalue_threshold = 0.05)  {
  # Filter for significant genes based on adjusted p-value and fold change thresholds
  filtered_genes <- dif_genes[abs(dif_genes$log2FoldChange) > log2fc_threshold & 
                                dif_genes$padj < adj_pvalue_threshold, ]
  
  ## Ad columns with EntrezID
  x <- bitr(rownames(filtered_genes), fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Mm.eg.db",drop = F)
  filtered_genes$SYMBOL <- rownames(filtered_genes)
  filtered_genes <- merge(filtered_genes,x,by="SYMBOL")
  filtered_genes <- na.omit(filtered_genes) #Remove Missing Values
  genelist <- extract2(filtered_genes, 'log2FoldChange') %>% set_names(filtered_genes$ENTREZID)
  # sort the list in decreasing order (required for clusterProfiler)
  genelist = sort(genelist, decreasing = TRUE)
  return(genelist)
}

dge_result_names <- c("M1.Combo.vs.NLRP3.DGE", "M1.Combo.vs.ADT.DGE", 
                      "M1.Combo.vs.U.DGE", "M1.NLRP3.vs.ADT.DGE",
                      "M1.NLRP3.vs.U.DGE", "M1.ADT.vs.U.DGE", 
                      "M2.Combo.vs.NLRP3.DGE", "M2.Combo.vs.ADT.DGE", 
                      "M2.Combo.vs.U.DGE", "M2.NLRP3.vs.ADT.DGE",
                      "M2.NLRP3.vs.U.DGE", "M2.ADT.vs.U.DGE",
                      "M2vsM1.Combo.DGE", "M2vsM1.NLRP3.DGE", 
                      "M2vsM1.ADT.DGE", "M2vsM1.U.DGE")


for (names in dge_result_names) {
  dge <- get(names)
  condition <- gsub('.DGE','',names)
  g1 <- create_genelist(dge)
  g2 <- create_genelist_Entrez(dge)
  assign(paste0(condition,'.genelist'), g1)
  assign(paste0(condition,'.egenelist'), g2)
  
}

#### Create folders to save output ####
base_dir <- "C:/Users/axi313/Documents/Bulk_RNAseq_UChicago_NLRP3/Pathway_Analysis"

# Create main folders
main_folders <- c("M1", "M2", "M2vsM1")
for (folder in main_folders) {
  dir.create(file.path(base_dir, folder), recursive = TRUE, showWarnings = FALSE)
}

# Function to determine the main folder
get_main_folder <- function(name) {
  if (grepl("^M2vsM1", name)) {
    return("M2vsM1")
  } else if (grepl("^M2", name)) {
    return("M2")
  } else {
    return("M1")
  }
}

# Loop through the result names to create subdirectories
for (name in dge_result_names) {
  clean_name <- gsub(".DGE", "", name) # Remove the '.DGE' part
  main_folder <- get_main_folder(clean_name) # Determine the main folder
  # Create the subdirectory within the appropriate main folder
  dir.create(file.path(base_dir, main_folder, clean_name), recursive = TRUE, showWarnings = FALSE)
}

####################### Working Dirs To Save Output ##########################
### M2.Combo.vs.NLRP3, "M1.NLRP3.vs.U", "M2.ADT.vs.U", "M1.NLRP3.vs.ADT", "M1.Combo.vs.ADT
#have too few DEGs - M2.NLRP3.vs.U no sig enrichment, "M2vsM1.Combo"

all.sample.names <-c("M1.Combo.vs.NLRP3", "M1.Combo.vs.U",
                     "M1.NLRP3.vs.ADT","M1.ADT.vs.U", "M2.Combo.vs.ADT", 
                     "M2.Combo.vs.U", "M2.NLRP3.vs.ADT",
                     "M2.NLRP3.vs.U", "M2vsM1.Combo", 
                     "M2vsM1.NLRP3", "M2vsM1.ADT", 
                     "M2vsM1.U")
all.sample.names.2 <-c( "M2vsM1.U")

for (runx in all.sample.names.2) {
  
    ### Move to Working Directory and prepare Genelists for gsea and e.gsea
    wd <- paste0('C:/Users/axi313/Documents/Bulk_RNAseq_UChicago_NLRP3/Pathway_Analysis/',gsub("\\..*", "", runx),'/',runx)
    setwd(wd)
    genelist <- get(paste0(runx,'.genelist'))
    egenelist <-get(paste0(runx,'.egenelist'))
    # Directory names
    directories <- c("GO", "KEGG", "KEGG_Pathview")
    
    # Create each directory if it doesn't already exist
    for (dir in directories) {
      if (!dir.exists(dir)) {
        dir.create(dir)
      } else {
        message(paste("Directory", dir, "already exists"))
      }
    }
    
    
    ############## GO Analysis ###########
    
    ### GO Enrichment (Overexpression Analysis)
    
    setwd(paste0(wd,'/GO'))
    
    ### GSEA Analysis
    go.gsea <- gseGO(geneList     = genelist,
                     OrgDb        = org.Mm.eg.db,
                     keyType = 'SYMBOL',
                     ont          = "All",
                     minGSSize    = 5,
                     maxGSSize    = 500,
                     pvalueCutoff = 0.05,
                     verbose      = TRUE)
    
    ### Results CSV
    
    df.gsea <- as.data.table(go.gsea)
    write.csv(df.gsea,paste0(runx,'_GO_GSEA_Enrichment.csv'))
    
    ### GO Dotplot
    
    dp.gsea <- dotplot(go.gsea, showCategory=30) + ggtitle("dotplot for GO GSEA")
    ggsave(paste0(runx,'_GO_Dotplot.png'),dpi=500, height =16, width = 13, dp.gsea)
    
    ### Gene-Concept Network
    
    go.gc.p1 <- cnetplot(go.gsea,color.params = list(foldChange = genelist))
    go.gc.p2 <- cnetplot(go.gsea, color.params = list(foldChange = genelist), circular = TRUE, colorEdge = TRUE) 
    go.gc.p3 <- cowplot::plot_grid(go.gc.p1, go.gc.p2, ncol=2, labels=LETTERS[1:2], rel_widths=c(.8, 1.2))
    ggsave(paste0(runx,'_GO_Cnetplot.png'),dpi=500, height =16, width = 21,bg='white', go.gc.p3)
    
    ### Heatplot
    hp.go.gsea <- heatplot(go.gsea, showCategory = 10, foldChange=genelist) 
    ggsave(paste0(runx,'_GO_Heatplot.png'),dpi=500, height =11, width = 28,bg='white', hp.go.gsea)
    
    ### Emap Plot
    go.gsea.tree <- pairwise_termsim(go.gsea)
    emap.go.gsea <- emapplot(go.gsea.tree)
    ggsave(paste0(runx,'_GO_Emmap.png'),dpi=500, height =11, width = 28,bg='white', emap.go.gsea)
    
    
    #### KEGG GSEA ####
    
    setwd(paste0(wd,'/KEGG'))
    
    kegg.gsea <- gseKEGG(
      geneList = egenelist,
      keyType = 'ncbi-geneid',
      organism = "mmu")
    
    ### Results CSV
    df.KEGG <- as.data.table(kegg.gsea)
    write.csv(df.KEGG,paste0(runx,'_KEGG_Enrichment.csv'))
    
    ### Dot Plot
    
    dp.kegg <- dotplot(kegg.gsea, showCategory=30) + ggtitle("dotplot for KEGG GSEA")
    ggsave(paste0(runx,'_KEGG_Dotplot.png'),dpi=500, height =16, width = 13, dp.kegg)
    
    ## Gene-Concept Network
    
    kegg.gsea.2 <- setReadable(kegg.gsea, 'org.Mm.eg.db', 'ENTREZID')
    kegg.p1 <- cnetplot(kegg.gsea.2,color.params = list(foldChange = genelist))
    kegg.p2 <- cnetplot(kegg.gsea.2, color.params = list(foldChange = genelist), circular = TRUE, colorEdge = TRUE) 
    kegg.p3 <- cowplot::plot_grid(kegg.p1, kegg.p2, ncol=2, labels=LETTERS[1:2], rel_widths=c(.8, 1.2))
    ggsave(paste0(runx,'_KEGG_Cnetplot.png'),dpi=500, height =16, width = 21,bg='white', kegg.p3)
    
    ### Heatplot
    
    hp.kegg.gsea <- heatplot(kegg.gsea.2, showCategory = 20, foldChange=genelist) 
    ggsave(paste0(runx,'_KEGG_Heatplot.png'),dpi=500, height =11, width = 28,bg='white', hp.kegg.gsea)
    
    ### Emap Plot
    go.kegg.tree <- pairwise_termsim(kegg.gsea.2)
    
    emap.kegg.gsea <- emapplot(go.kegg.tree)
    ggsave(paste0(runx,'_KEGG_Emmap.png'),dpi=500, height =11, width = 28,bg='white', emap.kegg.gsea)
    
    #### Pathway Overlay
    image_directory <- paste0(wd,'/KEGG_Pathview')
    setwd(image_directory)
    
    # Get Ids of pathways where adj Pval < 0.05
    
    KEGG.pathways <- df.KEGG[df.KEGG$p.adjust <0.05,]$ID
    
    for (KP in KEGG.pathways) {
      tryCatch({
        # Attempt to generate the plot with pathview
        pathview(
          gene.data = egenelist,
          pathway.id = KP,
          species = "mmu",
          limit = list(gene = max(abs(egenelist)), cpd = 1),
          new.signature = FALSE
        )
        dev.off()
      }, error = function(e) {
        # In case of an error, print the pathway ID and the error message
        message(sprintf("Error with pathway %s: %s", KP, e$message))
        # The loop will continue despite the error
      })
    }
    
    #Delete Extra Pathview files
    
    # List all files in the directory
    all_files <- list.files(image_directory, full.names = TRUE)
    
    # Identify the files that do not end with '.pathview.png'
    files_to_delete <- all_files[!grepl("\\.pathview\\.png$", all_files)]
    
    # Remove the identified files
    file.remove(files_to_delete)
    
}


