#RUNNING DEG ANALYSIS OF GOF_CF, GMF_CF, AND GMF_GOF using only female data from Macrogen sequenced dataset

#Removing all objects from the R environment
rm(list = ls())

getwd()
setwd("C:/R/2024/Differential_geneExpression_Sk/tsetse_DEGs")
library(rlang)
library(tidyverse)
library(dplyr)
#if (!require("BiocManager", quietly = TRUE))
#install.packages("BiocManager")
#BiocManager::install("edgeR")
library(edgeR)
library(ggplot2)
library(ggrepel)
library(Rcpp)
library(statmod)
library(limma)
#BiocManager::install("Glimma")
library(Glimma)
library(gplots)
library(RColorBrewer)
library(writexl)

#LOADING & PREPARING DATA

# Reading excel file into R
file_path <- "C:/R/2024/Differential_geneExpression_Sk/tsetse_DEGs/MacrogenFemOnlyCOUNTS.csv"
data <- read.csv(file_path)

# Print top list of the data to confirm correctness
head(data)
dim(data)

# Creating individual objects and saving them as files
for (i in 2:ncol(data)) {
  # Selecting the first and the ith columns
  subset_data <- data[, c(1, i)]
  
  # Getting the column name
  column_name <- colnames(data)[i]
  
  # Creating a dynamic object name
  object_name <- column_name
  
  # Assigning the subset data frame to a new object
  assign(object_name, subset_data)
  
  # Creating filename for each subset
  file_name <- paste0(column_name, ".xlsx")
  
  # Write the subset to a new Excel file
  write_xlsx(subset_data, file_name)
}

print(CF1)

alltranscripts <- data
head(alltranscripts)

#GENE COUNTS
#allmygenes <- read.csv("Expression_Profile.GpallidipesIAEA.gene.csv",sep = '\t', header = TRUE)
#head(allmygenes)
#print(allmygenes)

#TRANSCRIPTS COUNTS
#allmytranscripts <- read.csv("Expression_Profile.GpallidipesIAEA.transcript.csv",sep = '\t', header = TRUE)
#print(allmytranscripts)
#head(allmytranscripts)
#Reading the sample metadata into R
sample_metadata <- read.csv("metadata_FemaleOnly.csv", sep = ',', header = TRUE)
print(sample_metadata)

# Remove all rows with any NA values from alltranscript
#alltranscripts_clean <- na.omit(alltranscripts)
#(alltranscripts_clean)


#Replacing NAs with 0
alltranscripts[is.na(alltranscripts)] <- 0
print(alltranscripts)

#Save in alltranscripts as all_transcripts.csv in my directory to edit the duplictes
#write.csv(alltranscripts, "ALL_transcriptsORS.csv")
#re-import the edited file
#alltranscripts <- read.csv("all_transcripts.csv", sep = ',', header = TRUE)


#converting the transcript name column into row names
read_count <- alltranscripts %>% remove_rownames %>% column_to_rownames(var="Gene_ID")
dim(read_count)

# checking the output
print(read_count)


# Converting raw counts to DGEList object / Create 'DGEList'
dge <- DGEList(counts = read_count, group = sample_metadata$sample)
dim(dge)
view(dge)

# PRE-FILTERING & NORMALIZATION
# Filter lowly expressed genes
# Filter low count genes (adjust as per data characteristics)
keep <- filterByExpr(dge, group = sample_metadata$sample)
dge <- dge[keep, , keep.lib.sizes = FALSE]
dim(dge)
view(dge)

#save filtered dge# this are total genes with more than 2 counts across all the 9 groups
#write.csv(non_zero_counts, file = "C:/R/2024/Differential_geneExpression_Sk/tsetse_DEGs/Results_crickDEGS/9343_filtered_Counts2.csv")

# Converting counts to a logical matrix where non-zero values are TRUE
non_zero_matrix <- dge$counts != 0
# Sum the TRUE values across each row
non_zero_counts <- rowSums(non_zero_matrix)
print(non_zero_counts)

# Identifying rows where the sum of non-zero values is greater than 1
rows_to_keep <- non_zero_counts > 2 # keeps r=genes with atleast two counts across all the samples

print(rows_to_keep)

# Count the number of TRUE values
true_count <- sum(rows_to_keep)

print(true_count)

# Filter the DGEList object to keep only interesting rows
dge_filtered <- dge[rows_to_keep, , keep.lib.sizes = FALSE]
print(dge_filtered)
dim(dge_filtered)

#Normalize the data:
#Normalize the data to account for differences in library sizes using TMM method.
#Normalize for sequencing depth and RNA composition biases
dge_filtered <- calcNormFactors(dge_filtered, method = "TMM")
dim(dge_filtered)

# visualize normalized library sizes
print(dge_filtered$counts)  # Print normalized library size

# Estimatating dispersion (Fisher's exact test)
dge_filtered <- estimateDisp(dge_filtered)

# Normalize counts
normalized_counts <- cpm(dge_filtered, log=TRUE)

# Performing hierarchical clustering
dist_matrix <- dist(t(normalized_counts))
hc <- hclust(dist_matrix)

# Plot dendrogram
plot(hc)

#check for data dispersions
plotBCV(dge_filtered)

# Conducting the Fisher's exact test

# 1. GOF_CF
etGOF_CF <- exactTest(dge_filtered, pair = c("inf_gutOnly_F", "control_F"))
# View the top differentially expressed genes
topTags(etGOF_CF)
#Setting regulation threshold (significance and log-fold change) to determine upregulated and downregulated genes
p_value_threshold <- 0.05
logFC_threshold <- 1

# Extract results
results <- topTags(etGOF_CF, n = Inf, adjust.method="none")$table
# Adding columns for significance
results$Significance <- ifelse(results$PValue < p_value_threshold & abs(results$logFC) > logFC_threshold, "Significant", "Not significant")
# Identify upregulated and downregulated genes
upregulated_genesGOF_CF <- results[results$PValue < p_value_threshold & results$logFC > logFC_threshold, ]
downregulated_genesGOF_CF <- results[results$PValue < p_value_threshold & results$logFC < -logFC_threshold, ]

# Print results
cat("upregulated_genesGOF_CF:\n")
print(upregulated_genesGOF_CF)
cat("\ndownregulated_genes_CF_GOF:\n")
print(downregulated_genesGOF_CF)

#save these files
write.csv(results, file = "C:/2025/R_Analysis/Analysis_of_DEGs/VectorBase_FunctionalEnrichment/R_EdgeR_DEGs/GOF-CF_DEGs.csv")
write.csv(upregulated_genesGOF_CF, file = "C:/2025/R_Analysis/Analysis_of_DEGs/VectorBase_FunctionalEnrichment/R_EdgeR_DEGs/GOF-CF_UpregulatedGenes.csv")
write.csv(downregulated_genesGOF_CF, file = "C:/2025/R_Analysis/Analysis_of_DEGs/VectorBase_FunctionalEnrichment/R_EdgeR_DEGs/GOF-CF_DownregulatedGenes.csv")

# Categorize regulation
results$Regulation <- ifelse(
  results$PValue < 0.05 & results$logFC >= 1, "Upregulated",
  ifelse(results$PValue < 0.05 & results$logFC <= -1, "Downregulated", "Not significant")
)

# Plotting results
# Convert row names to a new column called Gene
results$Gene <- rownames(results)

# Plot
myplot1 <- ggplot(results, aes(x = logFC, y = -log10(PValue), color = Regulation)) +
  geom_point(alpha = 0.8, size = 1.5) +
  scale_color_manual(
    values = c(
      "Upregulated" = "red",
      "Downregulated" = "blue",
      "Not significant" = "gray"
    )
  ) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "blue", size = 0.5) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "blue", size = 0.5) +
  theme_minimal() +
  labs(
    title = "GOF-Control",
    x = "Log2 Fold Change",
    y = "-Log10 P-Value"
  ) +
  theme(legend.title = element_blank()) +
  # Label only significant genes
  geom_text_repel(
    data = subset(results, Regulation != "Not significant"),
    aes(label = Gene),   # <-- Make sure your data has a column named Gene
    size = 3,
    max.overlaps = 10
  )

# Print
print(myplot1)

# Output PDF file path
output_pdf1 <- "C:/R/2024/Differential_geneExpression_Sk/tsetse_DEGs/Results_crickDEGS/GOF_CF_Data/my_ggplot_GOF_CF27082025macrogenfemOnly.pdf"

# Creating the directory if it doesn't exist
dir.create(dirname(output_pdf1), recursive = TRUE, showWarnings = FALSE)

# Saving the ggplot copy as a PDF
ggsave(filename = output_pdf1,
       plot = myplot1,
       device = "pdf",
       width = 8, height = 6, units = "in")


#2. GMF_vs_CF

etGMF_CF <- exactTest(dge_filtered, pair = c("inf_gutMouth_F", "control_F"))

# View the top differentially expressed genes
topTags(etGMF_CF)

# Setting regulation thresholds
p_value_threshold <- 0.05
logFC_threshold <- 1

# Extract results
resultsGMF_CF <- topTags(etGMF_CF, n = Inf, adjust.method="none")$table
# Add columns for significance
resultsGMF_CF$Significance <- ifelse(resultsGMF_CF$PValue < p_value_threshold & abs(resultsGMF_CF$logFC) > logFC_threshold, "Significant", "Not significant")

# Identifying upregulated and downregulated genes
upregulated_genesGMF_CF <- resultsGMF_CF[resultsGMF_CF$PValue < p_value_threshold & resultsGMF_CF$logFC > logFC_threshold, ]
downregulated_genesGMF_CF <- resultsGMF_CF[resultsGMF_CF$PValue < p_value_threshold & resultsGMF_CF$logFC < -logFC_threshold, ]

# Print results
cat("upregulated_genesGMF_CF:\n")
print(upregulated_genesGMF_CF)
cat("\ndownregulated_genesGMF_CF:\n")
print(downregulated_genesGMF_CF)

#Saving the files
write.csv(resultsGMF_CF, file = "C:/2025/R_Analysis/Analysis_of_DEGs/VectorBase_FunctionalEnrichment/R_EdgeR_DEGs/GMF-CF_DEGs.csv")
write.csv(upregulated_genesGMF_CF, file = "C:/2025/R_Analysis/Analysis_of_DEGs/VectorBase_FunctionalEnrichment/R_EdgeR_DEGs/GMF-CF_UpregulatedGenes.csv")
write.csv(downregulated_genesGMF_CF, file = "C:/2025/R_Analysis/Analysis_of_DEGs/VectorBase_FunctionalEnrichment/R_EdgeR_DEGs/GMF-CF_DownregulatedGenes.csv")

# Categorizing regulation
resultsGMF_CF$Regulation <- ifelse(
  resultsGMF_CF$PValue < 0.05 & resultsGMF_CF$logFC >= 1, "Upregulated",
  ifelse(resultsGMF_CF$PValue < 0.05 & resultsGMF_CF$logFC <= -1, "Downregulated", "Not significant")
)

# Converting row names to a new column called Gene
resultsGMF_CF$Gene <- rownames(resultsGMF_CF)

# Plot
myplot2 <- ggplot(resultsGMF_CF, aes(x = logFC, y = -log10(PValue), color = Regulation)) +
  geom_point(alpha = 0.8, size = 1.5) +
  scale_color_manual(
    values = c(
      "Upregulated" = "Red",
      "Downregulated" = "blue",
      "Not significant" = "gray"
    )
  ) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "blue", size = 0.5) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "blue", size = 0.5) +
  theme_minimal() +
  labs(
    title = "GMF_Vs_Control",
    x = "Log2 Fold Change",
    y = "-Log10 P-Value"
  ) +
  theme(legend.title = element_blank()) +
  # Labelling only significant genes
  geom_text_repel(
    data = subset(resultsGMF_CF, Regulation != "Not significant"),
    aes(label = Gene),   # <-- Make sure your data has a column named Gene
    size = 3,
    max.overlaps = 10
  )

# Print
print(myplot2)

# Defining the path to save a copy as PDF file
output_pdf2 <- "C:/R/2024/Differential_geneExpression_Sk/tsetse_DEGs/Results_crickDEGS/GMF_CF_Data/my_ggplot_GMF_CF27082025macrogenfemOnly.pdf"

# Creating a directory if it doesn't exist
dir.create(dirname(output_pdf2), recursive = TRUE, showWarnings = FALSE)

# Save the ggplot as a PDF
ggsave(filename = output_pdf2,
       plot = myplot2,
       device = "pdf",
       width = 8, height = 6, units = "in")


#3. GMF_vs_GOF

etGMF_GOF <- exactTest(dge_filtered, pair = c("inf_gutMouth_F", "inf_gutOnly_F"))

# View the top differentially expressed genes
topTags(etGMF_GOF)

# Set thresholds
p_value_threshold <- 0.05
logFC_threshold <- 1

# Extract results
resultsGMF_GOF <- topTags(etGMF_GOF, n = Inf, adjust.method="none")$table
# Adding columns for significance
resultsGMF_GOF$Significance <- ifelse(resultsGMF_GOF$PValue < p_value_threshold & abs(resultsGMF_GOF$logFC) > logFC_threshold, "Significant", "Not significant")

# Identifying upregulated and downregulated genes
upregulated_GMF_GOF <- resultsGMF_GOF[resultsGMF_GOF$PValue < p_value_threshold & resultsGMF_GOF$logFC > logFC_threshold, ]
downregulated_GMF_GOF <- resultsGMF_GOF[resultsGMF_GOF$PValue < p_value_threshold & resultsGMF_GOF$logFC < -logFC_threshold, ]

# Print results
cat("upregulated_GMF_GOF:\n")
print(upregulated_GMF_GOF)
cat("\ndownregulated_GMF_GOF:\n")
print(downregulated_GMF_GOF)

#Saving results into a file
write.csv(resultsGMF_GOF, file = "C:/2025/R_Analysis/Analysis_of_DEGs/VectorBase_FunctionalEnrichment/R_EdgeR_DEGs/GMF-GOF_DEGs.csv")
write.csv(upregulated_GMF_GOF, file = "C:/2025/R_Analysis/Analysis_of_DEGs/VectorBase_FunctionalEnrichment/R_EdgeR_DEGs/GMF-GOF_UpregulatedGenes.csv")
write.csv(downregulated_GMF_GOF, file = "C:/2025/R_Analysis/Analysis_of_DEGs/VectorBase_FunctionalEnrichment/R_EdgeR_DEGs/GMF-GOF_DownregulatedGenes.csv")

# Categorizing regulation
resultsGMF_GOF$Regulation <- ifelse(
  resultsGMF_GOF$PValue < 0.05 & resultsGMF_GOF$logFC >= 1, "Upregulated",
  ifelse(resultsGMF_GOF$PValue < 0.05 & resultsGMF_GOF$logFC <= -1, "Downregulated", "Not significant")
)

# Convert row names to a new column called Gene
resultsGMF_GOF$Gene <- rownames(resultsGMF_GOF)

#Plotting
myplot3 <- ggplot(resultsGMF_GOF, aes(x = logFC, y = -log10(PValue), color = Regulation)) +
  geom_point(alpha = 0.8, size = 1.5) +
  scale_color_manual(
    values = c(
      "Upregulated" = "red",
      "Downregulated" = "blue",
      "Not significant" = "gray"
    )
  ) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "blue", size = 0.5) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "blue", size = 0.5) +
  theme_minimal() +
  labs(
    title = "GMF_Vs_GOF",
    x = "Log2 Fold Change",
    y = "-Log10 P-Value"
  ) +
  theme(legend.title = element_blank()) +
  # Labelling only significant genes
  geom_text_repel(
    data = subset(results, Regulation != "Not significant"),
    aes(label = Gene),   # <-- Make sure your data has a column named Gene
    size = 3,
    max.overlaps = 10
  )

# Print
print(myplot3)

# Defining a path where to save a PDF copy
output_pdf3 <- "C:/R/2024/Differential_geneExpression_Sk/tsetse_DEGs/Results_crickDEGS/GMF_GOF_Data/my_ggplot_GMF_GOF27082025macrogenfemOnly.pdf"

# Creating the directory if it doesn't exist
dir.create(dirname(output_pdf3), recursive = TRUE, showWarnings = FALSE)

# Save the ggplot file as a PDF
ggsave(filename = output_pdf3,
       plot = myplot3,
       device = "pdf",
       width = 8, height = 6, units = "in")

##All gene enrichment analysis were done using Vectorbase and results plotted in R

#Plotting Enrichment results in R
library(ggplot2)

#Enrichment analysis
setwd("C:/2025/R_Analysis/Analysis_of_DEGs/VectorBase_FunctionalEnrichment/R_EdgeR_DEGs/GMF-CF_VectorBaseAnnotation")

#Reading the enrichment results into R
Enriched <- read.csv("EnrichmentResult_GMF-CF_MF.csv", sep = ',', header = TRUE)
# Note: All the enrichment results for the other 2comparison groups were loaded in similar manner 

# Creating bubble plots
ggplot(Enriched, aes(x = Benjamini, y = Name, size = Gene_count, color = Regulation)) +
  geom_point(alpha = 0.7) +
  facet_grid(rows = vars(Category), scales = "free_y", space = "free") +  # Grouping by Category
  scale_size(range = c(3, 10)) +  # Adjusting bubble size
  scale_color_manual(values = c("Upregulated" = "red", "Downregulated" = "blue")) +  # Color for Up/Down
  theme_minimal() +
  labs(title = "GMF  vs Control",
       x = "adjusted P value", y = "Function", color = "Regulation", size = "Gene Count") +
  theme(
    axis.text.y = element_text(size = 10),  # Adjusting size of y-axis text
    panel.grid.major.y = element_line(color = "gray85"),
    panel.grid.major.x = element_line(color = "gray85"),  
    strip.background = element_rect(color = "black", fill = "white"),  # Dark line to separate categories
    panel.spacing = unit(0.8, "lines")  # Adjust panel spacing if needed
  )


##DONE##