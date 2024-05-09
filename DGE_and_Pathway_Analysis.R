### Libraries
library(DESeq2)


### Read Files ####

setwd("C:/Users/axi313/Documents/Bulk_RNAseq_UChicago_NLRP3")
load.path <- "C:/Users/axi313/Documents/Bulk_RNAseq_UChicago_NLRP3/Raw_Data/"

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
