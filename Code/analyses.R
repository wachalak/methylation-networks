## Load packages
library(dplyr)
library(huge)
library(conflicted)
library(tidyverse)
library(readr)
library(biomaRt)
library(clusterProfiler)
library(org.Ss.eg.db)
library(enrichplot)
library(GO.db)
library(ggplot2)
library(egg)
library(gridExtra)
library(gplots)
library(viridis) 
library(reshape2)
library(cluster)
library(rags2ridges)
library(huge)
library(igraph)

conflicts_prefer(base::setdiff)

################################## LOAD DATA ###################################

# Set working directory
setwd("C:/Users/PRACA/Desktop/RP/Methylation_networks/Data")

bed_files <- list.files(pattern = "\\.bed$")
bed_dataframes <- list()
for (file in bed_files) {
  # Define the dataframe name based on the file name (without extension)
  df_name <- sub("\\.bed$", "", file)
  # Read the BED file into a dataframe
  assign(df_name, read.table(file, header = TRUE))
  # Add the dataframe to the list
  bed_dataframes[[df_name]] <- get(df_name)
}

TPM <- read.table("reference_genes_TPM.tsv", sep="\t", 
                   header = TRUE, row.names = 1, 
                   stringsAsFactors = FALSE)


# Get a list of all dataframes starting with 'final_'
dataframe_list <- ls(pattern = '^final_.*')

# Loop through each dataframe and remove X chromosome
for (df_name in dataframe_list) {
  df <- get(df_name)
  assign(df_name, df[df$chromosome != 'X', ])
}

##################### DATA EXPLORATION - GET MEDIANS ###########################

calculate_medians_and_tidy <- function(df, tissue_name) {
  # Print column names to debug
  print(colnames(df))
  
  median_df <- df %>%
    group_by(gene) %>%
    summarise_at(c("meth_1", "meth_2", "meth_3"), median, na.rm = TRUE) %>%
    remove_rownames() %>%
    column_to_rownames(var = "gene")

  tidy_df <- assign(paste0(tissue_name, "_tidy"), median_df, envir = .GlobalEnv)
  return(tidy_df)
}

# Define a list of tissues
tissues <- c("liver_30dpf", "liver_70dpf", "liver_NB",
             "hindbrain_30dpf", "hindbrain_70dpf", "hindbrain_NB",
             "skin_30dpf", "skin_70dpf", "skin_NB",
             "lung_30dpf", "lung_70dpf", "lung_NB",
             "intestine_30dpf", "intestine_70dpf", "intestine_NB",
             "muscle_30dpf", "muscle_70dpf", "muscle_NB",
             "kidney_30dpf", "kidney_70dpf", "kidney_NB")

# Loop through tissues and perform data exploration
for (tissue in tissues) {
  # Get the dataframe for the specific tissue
  tissue_df <- get(paste0("final_", tissue))
  
  # Print tissue name for debugging
  print(paste("Processing tissue:", tissue))

  # Calculate medians and tidy the dataframe
  calculate_medians_and_tidy(tissue_df, tissue)
}

############### DATA PROCESSING - SCALE, TRANSPOSE & TRANSFORM DATA ####################

# Remove genes for which all methylation = 0 to avoid NAs
## Transpose to get samples in rows and variables in columns
### Scale so results from downstream analyses are not driven by differences in scale
#### Nonparanormal transformation to meet the assumptions of Gaussian Graphical Models - multivariate Gaussian

process_data <- function(df_tidy, tissue_name) {
  # Remove genes with all methylation = 0
  df_tidy <- df_tidy[rowSums(df_tidy != 0, na.rm = TRUE) > 0, ]
  
  # Scale data
  scaled_data <- scale(t(df_tidy))
  
  # Nonparanormal transformation
  npn_data <- huge.npn(scaled_data, npn.func = "shrinkage", npn.thresh = NULL, verbose = TRUE)
  
  # Assign the processed data to global environment
  assign(paste0(tissue_name, "_T"), npn_data, envir = .GlobalEnv)
}

# Loop through tissues and perform data processing
for (tissue in tissues) {
  # Get the tissue's scaled and transformed dataframe
  tissue_tidy <- get(paste0(tissue, "_tidy"))
  
  # Process data for the specific tissue
  process_data(tissue_tidy, tissue)
}

############ DATA EXPLORATION - INTERSECT TO OBTAIN UNION OF GENES #############

# List of tissues
tissues <- c("hindbrain", "liver", "skin", "intestine", "lung", "muscle", "kidney")

# Initialize an empty list to store unique gene sets for each tissue
gene_sets <- list()

# Loop through tissues
for (tissue in tissues) {
  # Get unique gene set for each tissue
  genes <- unique(colnames(get(paste0(tissue, "_30dpf_T"))))
  gene_sets[[paste0(tissue, "_30")]] <- genes
  
  genes <- unique(colnames(get(paste0(tissue, "_70dpf_T"))))
  gene_sets[[paste0(tissue, "_70")]] <- genes
  
  genes <- unique(colnames(get(paste0(tissue, "_NB_T"))))
  gene_sets[[paste0(tissue, "_NB")]] <- genes
}

# Get the union of all gene sets
all_genes <- unique(unlist(gene_sets))
write.csv(all_genes, "unique_genes.csv", row.names=TRUE)

####################### NETWORKS - SAMPLES INTERSECTION ########################
conflicts_prefer(base::intersect)    # base::intersect over any other package
## to avoid conflict

# ALL SAMPLES
col_extracted <- Reduce(intersect, list(colnames(liver_30dpf_T),
                                        colnames(liver_70dpf_T),
                                        colnames(liver_NB_T),
                                        colnames(lung_30dpf_T),
                                        colnames(lung_70dpf_T),
                                        colnames(lung_NB_T),
                                        colnames(intestine_30dpf_T),
                                        colnames(intestine_70dpf_T),
                                        colnames(intestine_NB_T),
                                        colnames(muscle_30dpf_T), 
                                        colnames(muscle_70dpf_T),
                                        colnames(muscle_NB_T),
                                        colnames(kidney_30dpf_T), 
                                        colnames(kidney_70dpf_T),
                                        colnames(kidney_NB_T),
                                        colnames(skin_30dpf_T),
                                        colnames(skin_70dpf_T),
                                        colnames(skin_NB_T),
                                        colnames(hindbrain_30dpf_T),
                                        colnames(hindbrain_70dpf_T),
                                        colnames(hindbrain_NB_T)))      # only 3 genes shared between all samples

# ENDODERM
col_extracted_endoderm <- Reduce(intersect, list(colnames(liver_30dpf_T), 
                                                 colnames(liver_70dpf_T),
                                                 colnames(liver_NB_T),
                                                 colnames(intestine_30dpf_T), 
                                                 colnames(intestine_70dpf_T),
                                                 colnames(intestine_NB_T),
                                                 colnames(lung_30dpf_T), 
                                                 colnames(lung_70dpf_T),
                                                 colnames(lung_NB_T)))

col_extracted_endoderm <- unique(col_extracted_endoderm)

write.csv(col_extracted_endoderm, "col_extracted_endoderm.csv", row.names=TRUE)

liver_30En <- liver_30dpf_T[, col_extracted_endoderm]
liver_70En <- liver_70dpf_T[, col_extracted_endoderm]
liver_NBEn <- liver_NB_T[, col_extracted_endoderm]

lung_30En <- lung_30dpf_T[, col_extracted_endoderm]
lung_70En <- lung_70dpf_T[, col_extracted_endoderm]
lung_NBEn <- lung_NB_T[, col_extracted_endoderm]

intestine_30En <- intestine_30dpf_T[, col_extracted_endoderm]
intestine_70En <- intestine_70dpf_T[, col_extracted_endoderm]
intestine_NBEn <- intestine_NB_T[, col_extracted_endoderm]

# Perform normality check
pdf("qq_plots.pdf")  # Open PDF device

for (i in 1:61) {
  # Open a new page for each plot
  if (i > 1)
    plot.new()
  
  # Set up the plot
  par(mfrow=c(1, 1))  # Set up a single plot per page
  qqnorm(test[,i], main = paste(i))
  qqline(test[,i])
}

dev.off()  # Close PDF device

# MESODERM
col_extracted_mesoderm <- Reduce(intersect, list(colnames(muscle_30dpf_T), 
                                                 colnames(muscle_70dpf_T),
                                                 colnames(muscle_NB_T),
                                                 colnames(kidney_30dpf_T), 
                                                 colnames(kidney_70dpf_T),
                                                 colnames(kidney_NB_T)))

write.csv(col_extracted_mesoderm, "col_extracted_mesoderm.csv", row.names=TRUE)

kidney_30Me <- kidney_30dpf_T[, col_extracted_mesoderm]
kidney_70Me <- kidney_70dpf_T[, col_extracted_mesoderm]
kidney_NBMe <- kidney_NB_T[, col_extracted_mesoderm]

muscle_30Me <- muscle_30dpf_T[, col_extracted_mesoderm]
muscle_70Me <- muscle_70dpf_T[, col_extracted_mesoderm]
muscle_NBMe <- muscle_NB_T[, col_extracted_mesoderm]

# ECTODERM
col_extracted_ectoderm <- Reduce(intersect, list(colnames(hindbrain_30dpf_T), 
                                                 colnames(hindbrain_70dpf_T),
                                                 colnames(hindbrain_NB_T),
                                                 colnames(skin_30dpf_T), 
                                                 colnames(skin_70dpf_T),
                                                 colnames(skin_NB_T)))

write.csv(col_extracted_ectoderm, "col_extracted_ectoderm.csv", row.names=TRUE)

hindbrain_30Ec <- hindbrain_30dpf_T[, col_extracted_ectoderm]
hindbrain_70Ec <- hindbrain_70dpf_T[, col_extracted_ectoderm]
hindbrain_NBEc <- hindbrain_NB_T[, col_extracted_ectoderm]

skin_30Ec <- skin_30dpf_T[,col_extracted_ectoderm]
skin_70Ec <- skin_70dpf_T[,col_extracted_ectoderm]
skin_NBEc <- skin_NB_T[,col_extracted_ectoderm]

############## DATA EXPLORATION - PRE-PROCESS DATA FOR CLUSTERING ##############

hind_30<-unique(colnames(hindbrain_30dpf_T))
hind_70<-unique(colnames(hindbrain_70dpf_T))
hind_NB<-unique(colnames(hindbrain_NB_T))
liv_30<-unique(colnames(liver_30dpf_T))
liv_70<-unique(colnames(liver_70dpf_T))
liv_NB<-unique(colnames(liver_NB_T))
skin_30<-unique(colnames(skin_30dpf_T))
skin_70<-unique(colnames(skin_70dpf_T))
skin_NB<-unique(colnames(skin_NB_T))
inte_30<-unique(colnames(intestine_30dpf_T))
inte_70<-unique(colnames(intestine_70dpf_T))
inte_NB<-unique(colnames(intestine_NB_T))
lung_30<-unique(colnames(lung_30dpf_T))
lung_70<-unique(colnames(lung_70dpf_T))
lung_NB<-unique(colnames(lung_NB_T))
musc_30<-unique(colnames(muscle_30dpf_T))
musc_70<-unique(colnames(muscle_70dpf_T))
musc_NB<-unique(colnames(muscle_NB_T))
kidn_30<-unique(colnames(kidney_30dpf_T))
kidn_70<-unique(colnames(kidney_70dpf_T))
kidn_NB<-unique(colnames(kidney_NB_T))


all <- c(hind_30, hind_70, hind_NB, liv_30, liv_70, liv_NB, skin_30, skin_70, 
         skin_NB, inte_30, inte_70, inte_NB, lung_30, lung_70, lung_NB, musc_30, 
         musc_70, musc_NB, kidn_30, kidn_70, kidn_NB)

all <- unique(all)

# Establish connection to Ensembl
ensembl <- useMart("ensembl")
ensembl_datasets <- listDatasets(ensembl)
ensembl <- useDataset("sscrofa_gene_ensembl", mart = ensembl)

# Get attributes and filters
filters <- listFilters(ensembl)
attributes <- listAttributes(ensembl)

# Retrieve gene information
t2g <- getBM(attributes = c('ensembl_gene_id', 'chromosome_name', 'start_position', 'end_position'), mart = ensembl)

# Assuming TPM is a data frame and ensembl_gene_id column exists
my_ids <- data.frame(row.names(TPM))
colnames(my_ids)[1] <- "ensembl_gene_id"

# Merge data frames
my_chrom <- merge(my_ids, t2g, by = 'ensembl_gene_id')

# Filter out rows containing "X" in column 2 of my_chrom
filtered_my_chrom <- subset(my_chrom, my_chrom$chromosome_name != "X")

# Match union with TPM file
conflicts_prefer(clusterProfiler::filter)

new_TPM <- TPM[row.names(TPM) %in% all, ]

# Filter new_TPM based on gene IDs present in filtered_my_chrom
new_TPM <- new_TPM[row.names(new_TPM) %in% filtered_my_chrom$ensembl_gene_id, ]
dim(new_TPM)

# Transform the data set for clustering
TPM_log <- log2(new_TPM + 1) # log transform and add 1
TPM_scaled <- data.frame(t(scale(t(TPM_log))))    # scale per gene
                                                  ## for plotting

# Convert df to a matrix
TPM_scaled <- as.matrix(TPM_scaled)   # for clustering 

################## DATA EXPLORATION - CLUSTERING ALL ###########################

# Cluster on samples
hc <- hclust(as.dist(1-cor(TPM_scaled, method="spearman")), method="complete")

# Cluster on genes
hr <- hclust(as.dist(1-cor(t(TPM_scaled), method="pearson")), method="complete") 

# Optimal number of clusters
set.seed(123)
sil <- rep(0, 20)
# repeat k-means for 1:20 and extract silhouette:
for(i in 2:20){
  k1to20 <- kmeans(t(TPM_scaled), centers = i, nstart = 25, iter.max = 20)
  ss <- silhouette(k1to20$cluster, dist(t(TPM_scaled)))
  sil[i] <- mean(ss[, 3])
}

# Plot the  average silhouette width
plot(1:20, sil, type = "b", pch = 19, xlab = "Number of clusters k", 
     ylab="Average silhouette width")
abline(v = which.max(sil), lty = 2)
cat("Average silhouette width optimal number of clusters:", 
    which.max(sil), "\n")

# Set the optimal number of cluster
hclustk16 = cutree(hc, k = 16)

# Visualize the bar indicating the clusters in combination with a heatmap
set.seed(123)
pdf(file = "CLUSTERING.pdf", width = 11, height = 11)

clustColBar <- viridis(16)
clustColBar <- clustColBar[as.vector(hclustk16)]

heatmap.2(TPM_scaled,
          Rowv = as.dendrogram(hr), 
          Colv = as.dendrogram(hc),
          col = inferno,
          scale = "none",
          margins = c(7, 7),
          cexCol = 0.6,
          labRow = F,
          main = NULL,
          trace = "none",
          ColSideColors=clustColBar,
          key = TRUE,
          key.title = 'Colour key')

dev.off()

########################## NETWORK EXTRACTION ENDODERM #########################

# Get correlation matrices
r_L30En <- cor(liver_30En)
r_L70En <- cor(liver_70En)
r_LNBEn <- cor(liver_NBEn)
r_Lu30En <- cor(lung_30En)
r_Lu70En <- cor(lung_70En)
r_LuNBEn <- cor(lung_NBEn)
r_I30En <- cor(intestine_30En)
r_I70En <- cor(intestine_70En)
r_INBEn <- cor(intestine_NBEn)

# Lists for constructing class-specific correlation matrices
Rlist <- list(r_L30En = r_L30En, r_L70En = r_L70En, r_LNBEn = r_LNBEn, 
              r_Lu30En = r_Lu30En, r_Lu70En = r_Lu70En, r_LuNBEn = r_LuNBEn, 
              r_I30En = r_I30En, r_I70En = r_I70En, r_INBEn = r_INBEn)

samps <- c(nrow(liver_30En), nrow(liver_70En), nrow(liver_NBEn),
           nrow(lung_30En), nrow(lung_70En), nrow(lung_NBEn),
           nrow(intestine_30En), nrow(intestine_70En), nrow(intestine_NBEn))

Tlist <- default.target.fused(Slist = Rlist, ns = samps, type = "DUPV")

Ylist <- list(L30Endata = liver_30En, L70Endata = liver_70En, 
              LNBEndata = liver_NBEn, Lu30Endata = lung_30En, 
              Lu70Endata = lung_70En, LuNBEndata = lung_NBEn,
              I30Endata = intestine_30En, I70Endata = intestine_70En,
              INBEndata = intestine_NBEn)

# Specify the matrix data
penalty_matrix_endoderm <- matrix(c(
    "ridge_liver_30", "fusion_L30-70", "0", "0", "0", "0", "0", "0", "0",
    "fusion_L30-70", "ridge_liver_70", "0", "0", "0", "0", "0", "0", "0",
    "0", "0", "ridge_liver_NB", "0", "0", "0", "0", "0", "0",
    "0", "0", "0", "ridge_lung_30", "0", "0", "fusion_Lu30-I30", "0", "0",
    "0", "0", "0", "0", "ridge_lung_70", "fusion_Lu70-NB", "0", "0", "0",
    "0", "0", "0", "0", "fusion_Lu70-NB", "ridge_lung_NB", "0", "0", "0",
    "0", "0", "0", "fusion_Lu30-I30", "0", "0", "ridge_intestine_30", "0", "0",
    "0", "0", "0", "0", "0", "0", "0", "ridge_intestine_70", "fusion_I70-NB",
    "0", "0", "0", "0", "0", "0", "0", "fusion_I70-NB", "ridge_intestine_NB"
), nrow = 9, byrow = TRUE)

# Specify the row and column names
rownames(penalty_matrix_endoderm) <- c("L30", "L70", "LNB", "Lu30", "Lu70", "LuNB", "I30", "I70", "INB")
colnames(penalty_matrix_endoderm) <- c("L30", "L70", "LNB", "Lu30", "Lu70", "LuNB", "I30", "I70", "INB")

# Print the matrix
print(penalty_matrix_endoderm)

# Find optimal regularized precision matrices
set.seed(123)

OPTf_endoderm <- optPenalty.fused(Ylist = Ylist, Tlist = Tlist,
                                  lambda = penalty_matrix_endoderm,
                                  cv.method = "LOOCV")

## Have a look at optimal penalties
OPTf_endoderm$lambda.unique
Plist <- ridgeP.fused(Slist=Rlist, ns=samps, lambda = OPTf_endoderm$lambda, maxit = 1000)
## CN plots
"do it separately for each data class - function needed"
CNplot(r_L30En, Iaids = TRUE, 
       lambdaMin = 1e-26, lambdaMax = 300, step = 100, vertical = TRUE, 
       value = OPTf_endoderm[["lambda.unique"]][["ridge_liver_30"]])


Plist <- unlist(Plist)
Plist <- as.matrix(Plist)
# Threshold the optimal penalty matrix
help(symm)
OPTf_endoderm$Plist$L30Endata <- symm(OPTf_endoderm$Plist$L30Endata)
OPTf_endoderm$Plist$L70Endata <- symm(OPTf_endoderm$Plist$L70Endata)
OPTf_endoderm$Plist$LNBEndata <- symm(OPTf_endoderm$Plist$LNBEndata)
OPTf_endoderm$Plist$Lu30Endata <- symm(OPTf_endoderm$Plist$Lu30Endata)
OPTf_endoderm$Plist$Lu70Endata <- symm(OPTf_endoderm$Plist$Lu70Endata)
OPTf_endoderm$Plist$LuNBEndata <- symm(OPTf_endoderm$Plist$LuNBEndata)
OPTf_endoderm$Plist$I30Endata <- symm(OPTf_endoderm$Plist$I30Endata)
OPTf_endoderm$Plist$I70Endata <- symm(OPTf_endoderm$Plist$I70Endata)
OPTf_endoderm$Plist$INBEndata <- symm(OPTf_endoderm$Plist$INBEndata)

P0s_endoderm <- sparsify.fused(OPTf_endoderm$Plist, threshold = 'localFDR', FDRcut = .99)

write.csv(P0s_endoderm, "P0s_endoderm.csv", row.names=TRUE)
write.csv(P0s_mesoderm, "P0s_mesoderm.csv", row.names=TRUE)
write.csv(P0s_ectoderm, "P0s_ectoderm.csv", row.names=TRUE)
conflicts_prefer(base::union)


# Take Union() of genes
TST_L.1.2_En <- Union(P0s_endoderm$L30Endata$sparseParCor, 
                      P0s_endoderm$L70Endata$sparseParCor)

TST_L.2.3_En <- Union(P0s_endoderm$L70Endata$sparseParCor, 
                      P0s_endoderm$LNBEndata$sparseParCor)

TST_Lu.1.2_En <- Union(P0s_endoderm$Lu30Endata$sparseParCor, 
                       P0s_endoderm$Lu70Endata$sparseParCor)

TST_Lu.2.3_En <- Union(P0s_endoderm$Lu70Endata$sparseParCor, 
                       P0s_endoderm$LuNBEndata$sparseParCor)

TST_I.1.2_En <- Union(P0s_endoderm$I30Endata$sparseParCor, 
                      P0s_endoderm$I70Endata$sparseParCor)

TST_I.2.3_En <- Union(P0s_endoderm$I70Endata$sparseParCor, 
                      P0s_endoderm$INBEndata$sparseParCor)

PC_L30_1.2_En <- TST_L.1.2_En$M1subset
PC_L70_1.2_En <- TST_L.1.2_En$M2subset
PC_L70_2.3_En <- TST_L.2.3_En$M1subset
PC_LNB_2.3_En <- TST_L.2.3_En$M2subset

PC_Lu30_1.2_En <- TST_Lu.1.2_En$M1subset
PC_Lu70_1.2_En <- TST_Lu.1.2_En$M2subset
PC_Lu70_2.3_En <- TST_Lu.2.3_En$M1subset
PC_LuNB_2.3_En <- TST_Lu.2.3_En$M2subset

PC_I30_1.2_En <- TST_I.1.2_En$M1subset
PC_I70_1.2_En <- TST_I.1.2_En$M2subset
PC_I70_2.3_En <- TST_I.2.3_En$M1subset
PC_INB_2.3_En <- TST_I.2.3_En$M2subset

########################## NETWORK VISUALIZATION ENDODERM ######################

# Obtain up-regulated and down-regulated genes
## LIVER 1.2

pattern_L_up_1.2_En <- subset(final_liver_30dpf, 
                              sign(final_liver_30dpf$logFC_1.2) ==  1, 
                              select = c('gene'))

pattern_L_1.2_UP_En <- unique(pattern_L_up_1.2_En[pattern_L_up_1.2_En$gene 
                                                  %in% rownames(PC_L30_1.2_En), 'gene'])

pattern_L_down_1.2_En <- subset(final_liver_30dpf, 
                                sign(final_liver_30dpf$logFC_1.2) ==  -1, 
                                select = c('gene'))

pattern_L_1.2_DOWN_En <- unique(pattern_L_down_1.2_En[pattern_L_down_1.2_En$gene 
                                                      %in% rownames(PC_L30_1.2_En), 'gene'])

Colors_L_1.2_En <- c()
color_L_1.2_up_En <- Colors_L_1.2_En[pattern_L_1.2_UP_En] <- "tomato3"
color_L_1.2_down_En <- Colors_L_1.2_En[pattern_L_1.2_DOWN_En]  <- 'seagreen3'

## LIVER 2.3
pattern_L_up_2.3_En <- subset(final_liver_30dpf, 
                              sign(final_liver_30dpf$logFC_2.3) ==  1, 
                              select = c('gene'))

pattern_L_2.3_UP_En <- unique(pattern_L_up_2.3_En[pattern_L_up_2.3_En$gene 
                                                  %in% rownames(PC_L70_2.3_En), 'gene'])

pattern_L_down_2.3_En <- subset(final_liver_30dpf, 
                                sign(final_liver_30dpf$logFC_2.3) ==  -1, 
                                select = c('gene'))

pattern_L_2.3_DOWN_En <- unique(pattern_L_down_2.3_En[pattern_L_down_2.3_En$gene 
                                                      %in% rownames(PC_L70_2.3_En), 'gene'])
Colors_L_2.3_En <- c()
color_L_2.3_up_En <- Colors_L_2.3_En[pattern_L_2.3_UP_En] <- "tomato3"
color_L_2.3_down_En <- Colors_L_2.3_En[pattern_L_2.3_DOWN_En]  <- 'seagreen3'

## LUNG 1.2
pattern_Lu_up_1.2_En <- subset(final_lung_30dpf, 
                               sign(final_lung_30dpf$logFC_1.2) ==  1, 
                               select = c('gene'))

pattern_Lu_1.2_UP_En <- unique(pattern_Lu_up_1.2_En[pattern_Lu_up_1.2_En$gene 
                                                    %in% rownames(PC_Lu30_1.2_En), 'gene'])

pattern_Lu_down_1.2_En <- subset(final_lung_30dpf, 
                                 sign(final_lung_30dpf$logFC_1.2) ==  -1, 
                                 select = c('gene'))

pattern_Lu_1.2_DOWN_En <- unique(pattern_Lu_down_1.2_En[pattern_Lu_down_1.2_En$gene 
                                                        %in% rownames(PC_Lu30_1.2_En), 'gene'])

Colors_Lu_1.2_En <- c()
color_Lu_1.2_up_En <- Colors_Lu_1.2_En[pattern_Lu_1.2_UP_En] <- "tomato3"
color_Lu_1.2_down_En <- Colors_Lu_1.2_En[pattern_Lu_1.2_DOWN_En]  <- 'seagreen3'

## LUNG 2.3

pattern_Lu_up_2.3_En <- subset(final_lung_30dpf, 
                               sign(final_lung_30dpf$logFC_2.3) ==  1, 
                               select = c('gene'))

pattern_Lu_2.3_UP_En <- unique(pattern_Lu_up_2.3_En[pattern_Lu_up_2.3_En$gene 
                                                    %in% rownames(PC_Lu70_2.3_En), 'gene'])

pattern_Lu_down_2.3_En <- subset(final_lung_30dpf, 
                                 sign(final_lung_30dpf$logFC_2.3) ==  -1, 
                                 select = c('gene'))

pattern_Lu_2.3_DOWN_En <- unique(pattern_Lu_down_2.3_En[pattern_Lu_down_2.3_En$gene 
                                                        %in% rownames(PC_Lu70_2.3_En), 'gene'])
Colors_Lu_2.3_En <- c()
color_Lu_2.3_up_En <- Colors_Lu_2.3_En[pattern_Lu_2.3_UP_En] <- "tomato3"
color_Lu_2.3_down_En <- Colors_Lu_2.3_En[pattern_Lu_2.3_DOWN_En]  <- 'seagreen3'

## INTESTINE 1.2
pattern_I_up_1.2_En <- subset(final_intestine_30dpf, 
                              sign(final_intestine_30dpf$logFC_1.2) ==  1, 
                              select = c('gene'))

pattern_I_1.2_UP_En <- unique(pattern_I_up_1.2_En[pattern_I_up_1.2_En$gene 
                                                  %in% rownames(PC_I30_1.2_En), 'gene'])

pattern_I_down_1.2_En <- subset(final_intestine_30dpf, 
                                sign(final_intestine_30dpf$logFC_1.2) ==  -1, 
                                select = c('gene'))

pattern_I_1.2_DOWN_En <- unique(pattern_I_down_1.2_En[pattern_I_down_1.2_En$gene 
                                                      %in% rownames(PC_I30_1.2_En), 'gene'])

Colors_I_1.2_En <- c()
color_I_1.2_up_En <- Colors_I_1.2_En[pattern_I_1.2_UP_En] <- "tomato3"
color_I_1.2_down_En <- Colors_I_1.2_En[pattern_I_1.2_DOWN_En]  <- 'seagreen3'

## INTESTINE 2.3

pattern_I_up_2.3_En <- subset(final_intestine_30dpf, 
                              sign(final_intestine_30dpf$logFC_2.3) ==  1, 
                              select = c('gene'))

pattern_I_2.3_UP_En <- unique(pattern_I_up_2.3_En[pattern_I_up_2.3_En$gene 
                                                  %in% rownames(PC_I70_2.3_En), 'gene'])

pattern_I_down_2.3_En <- subset(final_intestine_30dpf, 
                                sign(final_intestine_30dpf$logFC_2.3) ==  -1, 
                                select = c('gene'))

pattern_I_2.3_DOWN_En <- unique(pattern_I_down_2.3_En[pattern_I_down_2.3_En$gene 
                                                      %in% rownames(PC_I70_2.3_En), 'gene'])
Colors_I_2.3_En <- c()
color_I_2.3_up_En <- Colors_I_2.3_En[pattern_I_2.3_UP_En] <- "tomato3"
color_I_2.3_down_En <- Colors_I_2.3_En[pattern_I_2.3_DOWN_En]  <- 'seagreen3'

# Visualize networks with FR algorithm
## LIVER
set.seed(111213)
pdf(file = "NETWORKS_liver_endoderm.pdf", width = 11, height = 11)
opar <- par(mfrow = c(1, 1))

Coords_L_1.2_En <- Ugraph(PC_L30_1.2_En, type = "fancy", lay = "layout_with_fr",
                          prune = FALSE, Vsize = 15, Vcex = 0.5, 
                          Vcolor = 'lightblue', main = "Liver 30dpf")

Ugraph(PC_L70_1.2_En, type = "fancy", lay = NULL, coords = Coords_L_1.2_En,
       prune = FALSE, Vsize = 12, Vcex = 0.5, Vcolor = 'lightblue',
       main = "Liver 70dpf")

DiffGraph(PC_L30_1.2_En, PC_L70_1.2_En, lay = NULL, coords = Coords_L_1.2_En,
          Vsize = 12, Vcex = 0.5, Vcolor = Colors_L_1.2_En,
          P1color = 'magenta',
          P2color = 'sienna1',
          main = "Differential Network Liver 30-70dpf")

Coords_L_2.3_En <- Ugraph(PC_L70_2.3_En, type = "fancy", lay = "layout_with_fr",
                          prune = FALSE, Vsize = 12, Vcex = 0.5,
                          Vcolor = 'lightblue', main = "Liver 70dpf")

Ugraph(PC_LNB_2.3_En, type = "fancy", lay = NULL, 
       coords = Coords_L_2.3_En,
       prune = FALSE, Vsize = 12, Vcex = 0.5,
       Vcolor = 'lightblue',
       main = "Liver NB")

DiffGraph(PC_L70_2.3_En, PC_LNB_2.3_En, lay = NULL, coords = Coords_L_2.3_En,
          Vsize = 12, Vcex = 0.5, Vcolor = Colors_L_2.3_En,
          P1color = 'magenta',
          P2color = 'sienna1',
          main = "Differential Network Liver 70dpf-NB")
dev.off()

## LUNG
set.seed(111213)
pdf(file = "NETWORKS_lung_endoderm.pdf", width = 11, height = 11)
opar <- par(mfrow = c(1, 1))

Coords_Lu_1.2_En <- Ugraph(PC_Lu30_1.2_En, type = "fancy", lay = "layout_with_fr",
                           prune = FALSE, Vsize = 12, Vcex = 0.5, 
                           Vcolor = 'lightblue', main = "Lung 30dpf")

Ugraph(PC_Lu70_1.2_En, type = "fancy", lay = NULL, coords = Coords_Lu_1.2_En,
       prune = FALSE, Vsize = 12, Vcex = 0.5, Vcolor = 'lightblue',
       main = "Lung 70dpf")

DiffGraph(PC_Lu30_1.2_En, PC_Lu70_1.2_En, lay = NULL, coords = Coords_Lu_1.2_En,
          Vsize = 12, Vcex = 0.5, Vcolor = Colors_Lu_1.2_En,
          P1color = 'magenta',
          P2color = 'sienna1',
          main = "Differential Network Lung 30-70dpf")

Coords_Lu_2.3_En <- Ugraph(PC_Lu70_2.3_En, type = "fancy", lay = "layout_with_fr",
                           prune = FALSE, Vsize = 12, Vcex = 0.5,
                           Vcolor = 'lightblue', main = "Lung 70dpf")

Ugraph(PC_LuNB_2.3_En, type = "fancy", lay = NULL, 
       coords = Coords_Lu_2.3_En,
       prune = FALSE, Vsize = 12, Vcex = 0.5,
       Vcolor = 'lightblue',
       main = "Lung NB")

DiffGraph(PC_Lu70_2.3_En, PC_LuNB_2.3_En, lay = NULL, coords = Coords_Lu_2.3_En,
          Vsize = 12, Vcex = 0.5, Vcolor = Colors_Lu_2.3_En,
          P1color = 'magenta',
          P2color = 'sienna1',
          main = "Differential Network Lung 70dpf-NB")

par(opar)
dev.off()

## INTESTINE
set.seed(111213)
pdf(file = "NETWORKS_intestine_endoderm.pdf", width = 11, height = 11)
opar <- par(mfrow = c(1, 1))

Coords_I_1.2_En <- Ugraph(PC_I30_1.2_En, type = "fancy", lay = "layout_with_fr",
                          prune = FALSE, Vsize = 12, Vcex = 0.5, 
                          Vcolor = 'lightblue', main = "Intestine 30dpf")

Ugraph(PC_I70_1.2_En, type = "fancy", lay = NULL, coords = Coords_I_1.2_En,
       prune = FALSE, Vsize = 12, Vcex = 0.5, Vcolor = 'lightblue',
       main = "Intestine 70dpf")

DiffGraph(PC_I30_1.2_En, PC_I70_1.2_En, lay = NULL, coords = Coords_I_1.2_En,
          Vsize = 12, Vcex = 0.5, Vcolor = Colors_I_1.2_En,
          P1color = 'magenta',
          P2color = 'sienna1',
          main = "Differential Network Intestine 30-70dpf")

Coords_I_2.3_En <- Ugraph(PC_I70_2.3_En, type = "fancy", lay = "layout_with_fr",
                          prune = FALSE, Vsize = 12, Vcex = 0.5,
                          Vcolor = 'lightblue', main = "Intestine 70dpf")

Ugraph(PC_INB_2.3_En, type = "fancy", lay = NULL, 
       coords = Coords_I_2.3_En,
       prune = FALSE, Vsize = 12, Vcex = 0.5,
       Vcolor = 'lightblue',
       main = "Intestine NB")

DiffGraph(PC_I70_2.3_En, PC_INB_2.3_En, lay = NULL, coords = Coords_I_2.3_En,
          Vsize = 12, Vcex = 0.5, Vcolor = Colors_I_2.3_En,
          P1color = 'magenta',
          P2color = 'sienna1',
          main = "Differential Network Intestine 70dpf-NB")

par(opar)
dev.off()

########################## NETWORK EXTRACTION MESODERM #########################

# Get correlation matrices
r_M30Me <- cor(muscle_30Me)
r_M70Me <- cor(muscle_70Me)
r_MNBMe <- cor(muscle_NBMe)
r_K30Me <- cor(kidney_30Me)
r_K70Me <- cor(kidney_70Me)
r_KNBMe <- cor(kidney_NBMe)

# Lists for constructing class-specific correlation matrices
Rlist <- list(r_M30Me = r_M30Me, r_M70Me = r_M70Me, r_MNBMe = r_MNBMe, 
              r_K30Me = r_K30Me, r_K70Me = r_K70Me, r_KNBMe = r_KNBMe)

samps <- c(nrow(muscle_30Me), nrow(muscle_70Me), nrow(muscle_NBMe),
           nrow(kidney_30Me), nrow(kidney_70Me), nrow(kidney_NBMe))

Tlist <- default.target.fused(Slist = Rlist, ns = samps, type = "DUPV")

Ylist <- list(M30Medata = muscle_30Me, M70Medata = muscle_70Me, 
              MNBMedata = muscle_NBMe, K30Medata = kidney_30Me, 
              K70Medata = kidney_70Me, KNBMedata = kidney_NBMe)

# Specify the matrix data
penalty_matrix_mesoderm <- matrix(c(
    "ridge_kidney_30", "fusion_K30-70", "0", "0", "0", "0",
    "fusion_K30-70", "ridge_kidney_70", "0", "0", "0", "0",
    "0", "0", "ridge_kidney_NB", "0", "0", "0",
    "0", "0", "0", "ridge_muscle_30", "0", "0",
    "0", "0", "0", "0", "ridge_muscle_70", "fusion_M70-NB",
    "0", "0", "0", "0", "fusion_M70-NB", "ridge_muscle_NB"
), nrow = 6, byrow = TRUE)

# Specify the row and column names
rownames(penalty_matrix_mesoderm) <- c("K30", "K70", "KNB", "M30", "M70", "MNB")
colnames(penalty_matrix_mesoderm) <- c("K30", "K70", "KNB", "M30", "M70", "MNB")

# Print the matrix
print(penalty_matrix_mesoderm)

# Find optimal regularized precision matrices
set.seed(123)
OPTf_mesoderm <- optPenalty.fused(Ylist = Ylist, Tlist = Tlist,
                                  lambda = penalty_matrix_mesoderm,
                                  cv.method = "LOOCV")

## Have a look at optimal penalties
OPTf_mesoderm$lambda.unique

## CN plots
"do it separately for each data class - function needed"
CNplot(r_K30Me, Iaids = TRUE, 
       lambdaMin = 1e-5, lambdaMax = 20, step = 100, vertical = TRUE, 
       value = OPTf_mesoderm[["lambda.unique"]][["ridge_kidney_30"]])

Plist <- unlist(Plist)
Plist <- as.matrix(Plist)
# Threshold the optimal penalty matrix
OPTf_mesoderm$Plist$K30Medata <- symm(OPTf_mesoderm$Plist$K30Medata)
OPTf_mesoderm$Plist$K70Medata <- symm(OPTf_mesoderm$Plist$K70Medata)
OPTf_mesoderm$Plist$KNBMedata <- symm(OPTf_mesoderm$Plist$KNBMedata)
OPTf_mesoderm$Plist$M30Medata <- symm(OPTf_mesoderm$Plist$M30Medata)
OPTf_mesoderm$Plist$M70Medata <- symm(OPTf_mesoderm$Plist$M70Medata)
OPTf_mesoderm$Plist$MNBMedata <- symm(OPTf_mesoderm$Plist$MNBMedata)


Plist <- OPTf_mesoderm$Plist

# Threshold the optimal penalty matrix
P0s_mesoderm <- sparsify.fused(OPTf_mesoderm$Plist, threshold = 'localFDR', FDRcut = .99)

# Take Union() of genes
TST_M.1.2_Me <- Union(P0s_mesoderm$M30Medata$sparseParCor, 
                      P0s_mesoderm$M70Medata$sparseParCor)
TST_M.2.3_Me <- Union(P0s_mesoderm$M70Medata$sparseParCor, 
                      P0s_mesoderm$MNBMedata$sparseParCor)
TST_K.1.2_Me <- Union(P0s_mesoderm$K30Medata$sparseParCor, 
                      P0s_mesoderm$K70Medata$sparseParCor)
TST_K.2.3_Me <- Union(P0s_mesoderm$K70Medata$sparseParCor, 
                      P0s_mesoderm$KNBMedata$sparseParCor)

PC_M30_1.2_Me <- TST_M.1.2_Me$M1subset
PC_M70_1.2_Me <- TST_M.1.2_Me$M2subset
PC_M70_2.3_Me <- TST_M.2.3_Me$M1subset
PC_MNB_2.3_Me <- TST_M.2.3_Me$M2subset

PC_K30_1.2_Me <- TST_K.1.2_Me$M1subset
PC_K70_1.2_Me <- TST_K.1.2_Me$M2subset
PC_K70_2.3_Me <- TST_K.2.3_Me$M1subset
PC_KNB_2.3_Me <- TST_K.2.3_Me$M2subset


########################## NETWORK VISUALIZATION MESODERM ######################

# Obtain up-regulated and down-regulated genes
## MUSCLE 1.2
pattern_M_up_1.2_Me <- subset(final_muscle_30dpf, 
                              sign(final_muscle_30dpf$logFC_1.2) ==  1, 
                              select = c('gene'))

pattern_M_1.2_UP_Me <- unique(pattern_M_up_1.2_Me[pattern_M_up_1.2_Me$gene 
                                                  %in% rownames(PC_M30_1.2_Me), 'gene'])

pattern_M_down_1.2_Me <- subset(final_muscle_30dpf, 
                                sign(final_muscle_30dpf$logFC_1.2) ==  -1, 
                                select = c('gene'))

pattern_M_1.2_DOWN_Me <- unique(pattern_M_down_1.2_Me[pattern_M_down_1.2_Me$gene 
                                                      %in% rownames(PC_M30_1.2_Me), 'gene'])

Colors_M_1.2_Me <- c()
color_M_1.2_up_Me <- Colors_M_1.2_Me[pattern_M_1.2_UP_Me] <- "tomato3"
color_M_1.2_down_Me <- Colors_M_1.2_Me[pattern_M_1.2_DOWN_Me]  <- 'seagreen3'

## MUSCLE 2.3
pattern_M_up_2.3_Me <- subset(final_muscle_30dpf, 
                              sign(final_muscle_30dpf$logFC_2.3) ==  1, 
                              select = c('gene'))

pattern_M_2.3_UP_Me <- unique(pattern_M_up_2.3_Me[pattern_M_up_2.3_Me$gene 
                                                  %in% rownames(PC_M70_2.3_Me), 'gene'])

pattern_M_down_2.3_Me <- subset(final_muscle_30dpf, 
                                sign(final_muscle_30dpf$logFC_2.3) ==  -1, 
                                select = c('gene'))

pattern_M_2.3_DOWN_Me <- unique(pattern_M_down_2.3_Me[pattern_M_down_2.3_Me$gene 
                                                      %in% rownames(PC_M70_2.3_Me), 'gene'])
Colors_M_2.3_Me <- c()
color_M_2.3_up_Me <- Colors_M_2.3_Me[pattern_M_2.3_UP_Me] <- "tomato3"
color_M_2.3_down_Me <- Colors_M_2.3_Me[pattern_M_2.3_DOWN_Me]  <- 'seagreen3'

## KIDNEY 1.2
pattern_K_up_1.2_Me <- subset(final_kidney_30dpf, 
                              sign(final_kidney_30dpf$logFC_1.2) ==  1, 
                              select = c('gene'))

pattern_K_1.2_UP_Me <- unique(pattern_K_up_1.2_Me[pattern_K_up_1.2_Me$gene 
                                                  %in% rownames(PC_K30_1.2_Me), 'gene'])

pattern_K_down_1.2_Me <- subset(final_kidney_30dpf, 
                                sign(final_kidney_30dpf$logFC_1.2) ==  -1, 
                                select = c('gene'))

pattern_K_1.2_DOWN_Me <- unique(pattern_K_down_1.2_Me[pattern_K_down_1.2_Me$gene 
                                                      %in% rownames(PC_K30_1.2_Me), 'gene'])

Colors_K_1.2_Me <- c()
color_K_1.2_up_Me <- Colors_K_1.2_Me[pattern_K_1.2_UP_Me] <- "tomato3"
color_K_1.2_down_Me <- Colors_K_1.2_Me[pattern_K_1.2_DOWN_Me]  <- 'seagreen3'

## KIDNEY 2.3
pattern_K_up_2.3_Me <- subset(final_kidney_30dpf, 
                              sign(final_kidney_30dpf$logFC_2.3) ==  1, 
                              select = c('gene'))

pattern_K_2.3_UP_Me <- unique(pattern_K_up_2.3_Me[pattern_K_up_2.3_Me$gene 
                                                  %in% rownames(PC_K70_2.3_Me), 'gene'])

pattern_K_down_2.3_Me <- subset(final_kidney_30dpf, 
                                sign(final_kidney_30dpf$logFC_2.3) ==  -1, 
                                select = c('gene'))

pattern_K_2.3_DOWN_Me <- unique(pattern_K_down_2.3_Me[pattern_K_down_2.3_Me$gene 
                                                      %in% rownames(PC_K70_2.3_Me), 'gene'])
Colors_K_2.3_Me <- c()
color_K_2.3_up_Me <- Colors_K_2.3_Me[pattern_K_2.3_UP_Me] <- "tomato3"
color_K_2.3_down_Me <- Colors_K_2.3_Me[pattern_K_2.3_DOWN_Me]  <- 'seagreen3'

# Visualize networks with FR algorithm
## MUSCLE
set.seed(111213)
pdf(file = "NETWORKS_muscle_mesoderm.pdf", width = 11, height = 11)
opar <- par(mfrow = c(1, 1))

Coords_M_1.2_Me <- Ugraph(PC_M30_1.2_Me, type = "fancy", lay = "layout_with_fr",
                          prune = FALSE, Vsize = 12, Vcex = 0.5, 
                          Vcolor = 'lightblue', main = "Muscle 30dpf")

Ugraph(PC_M70_1.2_Me, type = "fancy", lay = NULL, coords = Coords_M_1.2_Me,
       prune = FALSE, Vsize = 12, Vcex = 0.5, Vcolor = 'lightblue',
       main = "Muscle 70dpf")

DiffGraph(PC_M30_1.2_Me, PC_M70_1.2_Me, lay = NULL, coords = Coords_M_1.2_Me,
          Vsize = 12, Vcex = 0.5, Vcolor = Colors_M_1.2_Me,
          P1color = 'magenta',
          P2color = 'sienna1',
          main = "Differential Network Muscle 30-70dpf")

Coords_M_2.3_Me <- Ugraph(PC_M70_2.3_Me, type = "fancy", lay = "layout_with_fr",
                          prune = FALSE, Vsize = 12, Vcex = 0.5,
                          Vcolor = 'lightblue', main = "Muscle 70dpf")

Ugraph(PC_MNB_2.3_Me, type = "fancy", lay = NULL, 
       coords = Coords_M_2.3_Me,
       prune = FALSE, Vsize = 12, Vcex = 0.5,
       Vcolor = 'lightblue',
       main = "Muscle NB")

DiffGraph(PC_M70_2.3_Me, PC_MNB_2.3_Me, lay = NULL, coords = Coords_M_2.3_Me,
          Vsize = 12, Vcex = 0.5, Vcolor = Colors_M_2.3_Me,
          P1color = 'magenta',
          P2color = 'sienna1',
          main = "Differential Network Muscle 70dpf-NB")

par(opar)
dev.off()

# Visualize networks with FR algorithm
## KIDNEY
set.seed(111213)
pdf(file = "NETWORKS_kidney_mesoderm.pdf", width = 11, height = 11)
opar <- par(mfrow = c(1, 1))

Coords_K_1.2_Me <- Ugraph(PC_K30_1.2_Me, type = "fancy", lay = "layout_with_fr",
                          prune = FALSE, Vsize = 12, Vcex = 0.5, 
                          Vcolor = 'lightblue', main = "Kidney 30dpf")

Ugraph(PC_K70_1.2_Me, type = "fancy", lay = NULL, coords = Coords_K_1.2_Me,
       prune = FALSE, Vsize = 12, Vcex = 0.5, Vcolor = 'lightblue',
       main = "Kidney 70dpf")

DiffGraph(PC_K30_1.2_Me, PC_K70_1.2_Me, lay = NULL, coords = Coords_K_1.2_Me,
          Vsize = 12, Vcex = 0.5, Vcolor = Colors_K_1.2_Me,
          P1color = 'magenta',
          P2color = 'sienna1',
          main = "Differential Network Kidney 30-70dpf")

Coords_K_2.3_Me <- Ugraph(PC_K70_2.3_Me, type = "fancy", lay = "layout_with_fr",
                          prune = FALSE, Vsize = 12, Vcex = 0.5,
                          Vcolor = 'lightblue', main = "Kidney 70dpf")

Ugraph(PC_KNB_2.3_Me, type = "fancy", lay = NULL, 
       coords = Coords_K_2.3_Me,
       prune = FALSE, Vsize = 12, Vcex = 0.5,
       Vcolor = 'lightblue',
       main = "Kidney NB")

DiffGraph(PC_K70_2.3_Me, PC_KNB_2.3_Me, lay = NULL, coords = Coords_K_2.3_Me,
          Vsize = 12, Vcex = 0.5, Vcolor = Colors_K_2.3_Me,
          P1color = 'magenta',
          P2color = 'sienna1',
          main = "Differential Network Kidney 70dpf-NB")

par(opar)
dev.off()

########################## NETWORK EXTRACTION ECTODERM #########################

# Get correlation matrices
r_H30Ec <- cor(hindbrain_30Ec)
r_H70Ec <- cor(hindbrain_70Ec)
r_HNBEc <- cor(hindbrain_NBEc)
r_S30Ec <- cor(skin_30Ec)
r_S70Ec <- cor(skin_70Ec)
r_SNBEc <- cor(skin_NBEc)

# Lists for constructing class-specific correlation matrices
Rlist <- list(r_H30Ec = r_H30Ec, r_H70Ec = r_H70Ec, r_HNBEc = r_HNBEc, 
              r_S30Ec = r_S30Ec, r_S70Ec = r_S70Ec, r_SNBEc = r_SNBEc)

samps <- c(nrow(hindbrain_30Ec), nrow(hindbrain_70Ec), nrow(hindbrain_NBEc),
           nrow(skin_30Ec), nrow(skin_70Ec), nrow(skin_NBEc))

Tlist <- default.target.fused(Slist = Rlist, ns = samps, type = "DUPV")

Ylist <- list(H30Ecdata = hindbrain_30Ec, H70Ecdata = hindbrain_70Ec, 
              HNBEcdata = hindbrain_NBEc, S30Ecdata = skin_30Ec, 
              S70Ecdata = skin_70Ec, SNBEcdata = skin_NBEc)

# Specify the matrix data
penalty_matrix_ectoderm <- matrix(c(
    "ridge_hindbrain_30", "fusion_H30-70", "0", "0", "0", "0",
    "fusion_H30-70", "ridge_hindbrain_70", "0", "0", "0", "0",
    "0", "0", "ridge_hindbrain_NB", "0", "0", "0",
    "0", "0", "0", "ridge_skin_30", "0", "0",
    "0", "0", "0", "0", "ridge_skin_70", "fusion_S70-NB",
    "0", "0", "0", "0", "fusion_S70-NB", "ridge_skin_NB"
), nrow = 6, byrow = TRUE)

# Specify the row and column names
rownames(penalty_matrix_ectoderm) <- c("H30", "H70", "HNB", "S30", "S70", "SNB")
colnames(penalty_matrix_ectoderm) <- c("H30", "H70", "HNB", "S30", "S70", "SNB")

# Print the matrix
print(penalty_matrix_ectoderm)

# Find optimal regularized precision matrices
set.seed(123)
OPTf_ectoderm <- optPenalty.fused(Ylist = Ylist, Tlist = Tlist,
                                  lambda = penalty_matrix_ectoderm,
                                  cv.method = "LOOCV")

## Have a look at optimal penalties
OPTf_ectoderm$lambda.unique

## CN plots
"do it separately for each data class - function needed"
CNplot(r_S30Ec, Iaids = TRUE, 
       lambdaMin = 1e-5, lambdaMax = 20, step = 100, vertical = TRUE, 
       value = OPTf_ectoderm[["lambda.unique"]][["ridge_skin_30"]])

Plist <- unlist(Plist)
Plist <- as.matrix(Plist)
# Threshold the optimal penalty matrix
OPTf_ectoderm$Plist$H30Ecdata <- symm(OPTf_ectoderm$Plist$H30Ecdata)
OPTf_ectoderm$Plist$H70Ecdata <- symm(OPTf_ectoderm$Plist$H70Ecdata)
OPTf_ectoderm$Plist$HNBEcdata <- symm(OPTf_ectoderm$Plist$HNBEcdata)
OPTf_ectoderm$Plist$S30Ecdata <- symm(OPTf_ectoderm$Plist$S30Ecdata)
OPTf_ectoderm$Plist$S70Ecdata <- symm(OPTf_ectoderm$Plist$S70Ecdata)
OPTf_ectoderm$Plist$SNBEcdata <- symm(OPTf_ectoderm$Plist$SNBEcdata)


Plist <- OPTf_mesoderm$Plist

# Threshold the optimal penalty matrix
P0s_ectoderm <- sparsify.fused(OPTf_ectoderm$Plist, threshold = 'localFDR', FDRcut = 0.99)

# Take Union() of genes

TST_H.1.2_Ec <- Union(P0s_ectoderm$H30Ecdata$sparseParCor, 
                       P0s_ectoderm$H70Ecdata$sparseParCor)
TST_H.2.3_Ec <- Union(P0s_ectoderm$H70Ecdata$sparseParCor, 
                       P0s_ectoderm$HNBEcdata$sparseParCor)
TST_S.1.2_Ec <- Union(P0s_ectoderm$S30Ecdata$sparseParCor, 
                       P0s_ectoderm$S70Ecdata$sparseParCor)
TST_S.2.3_Ec <- Union(P0s_ectoderm$S70Ecdata$sparseParCor, 
                       P0s_ectoderm$SNBEcdata$sparseParCor)

PC_H30_1.2_Ec <- TST_H.1.2_Ec$M1subset
PC_H70_1.2_Ec <- TST_H.1.2_Ec$M2subset
PC_H70_2.3_Ec <- TST_H.2.3_Ec$M1subset
PC_HNB_2.3_Ec <- TST_H.2.3_Ec$M2subset

PC_S30_1.2_Ec <- TST_S.1.2_Ec$M1subset
PC_S70_1.2_Ec <- TST_S.1.2_Ec$M2subset
PC_S70_2.3_Ec <- TST_S.2.3_Ec$M1subset
PC_SNB_2.3_Ec <- TST_S.2.3_Ec$M2subset


########################## NETWORK VISUALIZATION ECTODERM ######################

# Obtain up-regulated and down-regulated genes

## HINDBRAIN 1.2

pattern_H_up_1.2_Ec <- subset(final_hindbrain_30dpf, 
                           sign(final_hindbrain_30dpf$logFC_1.2) ==  1, 
                           select = c('gene'))

pattern_H_1.2_UP_Ec <- unique(pattern_H_up_1.2_Ec[pattern_H_up_1.2_Ec$gene 
                                     %in% rownames(PC_H30_1.2_Ec), 'gene'])

pattern_H_down_1.2_Ec <- subset(final_hindbrain_30dpf, 
                             sign(final_hindbrain_30dpf$logFC_1.2) ==  -1, 
                             select = c('gene'))

pattern_H_1.2_DOWN_Ec <- unique(pattern_H_down_1.2_Ec[pattern_H_down_1.2_Ec$gene 
                                         %in% rownames(PC_H30_1.2_Ec), 'gene'])

Colors_H_1.2_Ec <- c()
color_H_1.2_up_Ec <- Colors_H_1.2_Ec[pattern_H_1.2_UP_Ec] <- "tomato3"
color_H_1.2_down_Ec <- Colors_H_1.2_Ec[pattern_H_1.2_DOWN_Ec]  <- 'seagreen3'

## HINDBRAIN 2.3

pattern_H_up_2.3_Ec <- subset(final_hindbrain_30dpf, 
                           sign(final_hindbrain_30dpf$logFC_2.3) ==  1, 
                           select = c('gene'))

pattern_H_2.3_UP_Ec <- unique(pattern_H_up_2.3_Ec[pattern_H_up_2.3_Ec$gene 
                                     %in% rownames(PC_H70_2.3_Ec), 'gene'])

pattern_H_down_2.3_Ec <- subset(final_hindbrain_30dpf, 
                             sign(final_hindbrain_30dpf$logFC_2.3) ==  -1, 
                             select = c('gene'))

pattern_H_2.3_DOWN_Ec <- unique(pattern_H_down_2.3_Ec[pattern_H_down_2.3_Ec$gene 
                                         %in% rownames(PC_H70_2.3_Ec), 'gene'])
Colors_H_2.3_Ec <- c()
color_H_2.3_up_Ec <- Colors_H_2.3_Ec[pattern_H_2.3_UP_Ec] <- "tomato3"
color_H_2.3_down_Ec <- Colors_H_2.3_Ec[pattern_H_2.3_DOWN_Ec]  <- 'seagreen3'

## SKIN 1.2

pattern_S_up_1.2_Ec <- subset(final_skin_30dpf, 
                              sign(final_skin_30dpf$logFC_1.2) ==  1, 
                              select = c('gene'))

pattern_S_1.2_UP_Ec <- unique(pattern_S_up_1.2_Ec[pattern_S_up_1.2_Ec$gene 
                                                  %in% rownames(PC_S30_1.2_Ec), 'gene'])

pattern_S_down_1.2_Ec <- subset(final_skin_30dpf, 
                                sign(final_skin_30dpf$logFC_1.2) ==  -1, 
                                select = c('gene'))

pattern_S_1.2_DOWN_Ec <- unique(pattern_S_down_1.2_Ec[pattern_S_down_1.2_Ec$gene 
                                                      %in% rownames(PC_S30_1.2_Ec), 'gene'])

Colors_S_1.2_Ec <- c()
color_S_1.2_up_Ec <- Colors_S_1.2_Ec[pattern_S_1.2_UP_Ec] <- "tomato3"
color_S_1.2_down_Ec <- Colors_S_1.2_Ec[pattern_S_1.2_DOWN_Ec]  <- 'seagreen3'

## SKIN 2.3

pattern_S_up_2.3_Ec <- subset(final_skin_30dpf, 
                              sign(final_skin_30dpf$logFC_2.3) ==  1, 
                              select = c('gene'))

pattern_S_2.3_UP_Ec <- unique(pattern_S_up_2.3_Ec[pattern_S_up_2.3_Ec$gene 
                                                  %in% rownames(PC_S70_2.3_Ec), 'gene'])

pattern_S_down_2.3_Ec <- subset(final_skin_30dpf, 
                                sign(final_skin_30dpf$logFC_2.3) ==  -1, 
                                select = c('gene'))

pattern_S_2.3_DOWN_Ec <- unique(pattern_S_down_2.3_Ec[pattern_S_down_2.3_Ec$gene 
                                                      %in% rownames(PC_S70_2.3_Ec), 'gene'])
Colors_S_2.3_Ec <- c()
color_S_2.3_up_Ec <- Colors_S_2.3_Ec[pattern_S_2.3_UP_Ec] <- "tomato3"
color_S_2.3_down_Ec <- Colors_S_2.3_Ec[pattern_S_2.3_DOWN_Ec]  <- 'seagreen3'


# Visualize networks with FR algorithm
## HINDBRAIN
set.seed(111213)
pdf(file = "NETWORKS_hindbrain_ectoderm.pdf", width = 11, height = 11)
opar <- par(mfrow = c(1, 1))

Coords_H_1.2_Ec <- Ugraph(PC_H30_1.2_Ec, type = "fancy", lay = 'layout_with_fr',
                          prune = FALSE, Vsize = 12, Vcex = 0.5, 
                          Vcolor = 'lightblue', main = "Hindbrain 30dpf")

Ugraph(PC_H70_1.2_Ec, type = "fancy", lay = NULL, coords = Coords_H_1.2_Ec,
       prune = FALSE, Vsize = 12, Vcex = 0.5, Vcolor = 'lightblue',
       main = "Hindbrain 70dpf")

DiffGraph(PC_H30_1.2_Ec, PC_H70_1.2_Ec, lay = NULL, coords = Coords_H_1.2_Ec,
          Vsize = 12, Vcex = 0.5, Vcolor = Colors_H_1.2_Ec,
          P1color = 'magenta',
          P2color = 'sienna1',
          main = "Differential Network Hindbrain 30-70dpf")

Coords_H_2.3_Ec <- Ugraph(PC_H70_2.3_Ec, type = "fancy", lay = "layout_with_fr",
                          prune = FALSE, Vsize = 12, Vcex = 0.5,
                          Vcolor = 'lightblue', main = "Hindbrain 70dpf")

Ugraph(PC_HNB_2.3_Ec, type = "fancy", lay = NULL, 
       coords = Coords_H_2.3_Ec,
       prune = FALSE, Vsize = 12, Vcex = 0.5,
       Vcolor = 'lightblue',
       main = "Hindbrain NB")

DiffGraph(PC_H70_2.3_Ec, PC_HNB_2.3_Ec, lay = NULL, coords = Coords_H_2.3_Ec,
          Vsize = 12, Vcex = 0.5, Vcolor = Colors_H_2.3_Ec,
          P1color = 'magenta',
          P2color = 'sienna1',
          main = "Differential Network Hindbrain 70dpf-NB")

par(opar)
dev.off()

# Visualize networks with FR algorithm
## SKIN
set.seed(111213)
pdf(file = "NETWORKS_skin_ectoderm.pdf", width = 11, height = 11)
opar <- par(mfrow = c(1, 1))

Coords_S_1.2_Ec <- Ugraph(PC_S30_1.2_Ec, type = "fancy", lay = "layout_with_fr",
                          prune = FALSE, Vsize = 12, Vcex = 0.5, 
                          Vcolor = 'lightblue', main = "Skin 30dpf")

Ugraph(PC_S70_1.2_Ec, type = "fancy", lay = NULL, coords = Coords_S_1.2_Ec,
       prune = FALSE, Vsize = 12, Vcex = 0.5, Vcolor = 'lightblue',
       main = "Skin 70dpf")

DiffGraph(PC_S30_1.2_Ec, PC_S70_1.2_Ec, lay = NULL, coords = Coords_S_1.2_Ec,
          Vsize = 12, Vcex = 0.5, Vcolor = Colors_S_1.2_Ec,
          P1color = 'magenta',
          P2color = 'sienna1',
          main = "Differential Network Skin 30-70dpf")

Coords_S_2.3_Ec <- Ugraph(PC_S70_2.3_Ec, type = "fancy", lay = "layout_with_fr",
                          prune = FALSE, Vsize = 12, Vcex = 0.5,
                          Vcolor = 'lightblue', main = "Skin 70dpf")

Ugraph(PC_SNB_2.3_Ec, type = "fancy", lay = NULL, 
       coords = Coords_S_2.3_Ec,
       prune = FALSE, Vsize = 12, Vcex = 0.5,
       Vcolor = 'lightblue',
       main = "Skin NB")

DiffGraph(PC_S70_2.3_Ec, PC_SNB_2.3_Ec, lay = NULL, coords = Coords_S_2.3_Ec,
          Vsize = 12, Vcex = 0.5, Vcolor = Colors_S_2.3_Ec,
          P1color = 'magenta',
          P2color = 'sienna1',
          main = "Differential Network Skin 70dpf-NB")

par(opar)
dev.off()

######################### NETWORK ANALYSIS ENDODERM ############################

# Obtain simple network statistics
## LIVER
PCLlist_L.1.2_En <- list(PC_L30_1.2_En = PC_L30_1.2_En, 
                         PC_L70_1.2_En = PC_L70_1.2_En)

NwkSTATSList_L.1.2_En <- GGMnetworkStats.fused(PCLlist_L.1.2_En)
write.csv(NwkSTATSList_L.1.2_En, "NwkSTATSList_L.1.2_En.csv", row.names=TRUE)

## LUNG
PCLlist_Lu.2.3_En <- list(PC_Lu70_2.3_En = PC_Lu70_2.3_En,
                          PC_LuNB_2.3_En = PC_LuNB_2.3_En)

NwkSTATSList_Lu.2.3_En <- GGMnetworkStats.fused(PCLlist_Lu.2.3_En)
write.csv(NwkSTATSList_Lu.2.3_En, "NwkSTATSList_Lu.2.3_En.csv", row.names=TRUE)

######################### NETWORK ANALYSIS MESODERM ############################

## KIDNEY
PCLlist_K.1.2_Me <- list(PC_K30_1.2_Me = PC_K30_1.2_Me, 
                         PC_K70_1.2_Me = PC_K70_1.2_Me)

NwkSTATSList_K.1.2_Me <- GGMnetworkStats.fused(PCLlist_K.1.2_Me)
write.csv(NwkSTATSList_K.1.2_Me, "NwkSTATSList_K.1.2_Me.csv", row.names=TRUE)

######################### NETWORK ANALYSIS ECTODERM ############################

## SKIN
PCLlist_S.2.3_Ec <- list(PC_S70_2.3_Ec = PC_S70_2.3_Ec,
                         PC_SNB_2.3_Ec = PC_SNB_2.3_Ec)

NwkSTATSList_S.2.3_Ec <- GGMnetworkStats.fused(PCLlist_S.2.3_Ec)
write.csv(NwkSTATSList_S.2.3_Ec, "NwkSTATSList_S.2.3_Ec.csv", row.names=TRUE)

############################ COMMUNITY DETECTION ################################

# Function to transform negative weights into positive ones
transform_negative_weights <- function(adj_matrix) {
  adj_matrix[adj_matrix < 0] <- abs(adj_matrix[adj_matrix < 0])
  return(adj_matrix)
}

# Function to visualize networks
visualize_networks <- function(adj_matrix_stage1, adj_matrix_stage2, tissue_stage_name) {
  # Transform negative weights into positive ones
  adj_matrix_stage1 <- transform_negative_weights(adj_matrix_stage1)
  adj_matrix_stage2 <- transform_negative_weights(adj_matrix_stage2)
  
  # Convert adjacency matrices to graph objects
  graph_stage1 <- graph_from_adjacency_matrix(adj_matrix_stage1, mode = "undirected", weighted = TRUE)
  graph_stage2 <- graph_from_adjacency_matrix(adj_matrix_stage2, mode = "undirected", weighted = TRUE)
  
  # Remove self-loops from the graph using igraph::simplify
  graph_stage1 <- igraph::simplify(graph_stage1)
  graph_stage2 <- igraph::simplify(graph_stage2)
  
  # Perform community detection using cluster_fluid_communities
  communities_stage1 <- cluster_fluid_communities(graph_stage1, 2)
  communities_stage2 <- cluster_fluid_communities(graph_stage2, 2)
  
  coords_1 <- layout_with_fr(graph_stage1)
  coords_2 <- layout_with_fr(graph_stage2)
  
  # Get community membership for each node
  membership_stage1 <- membership(communities_stage1)
  membership_stage2 <- membership(communities_stage2)
  
  # Define unique node shapes for each community type
  node_shapes <- c("circle", "square")
  
  # Assign shapes based on community membership
  vertex_shapes_stage1 <- sapply(membership_stage1, function(x) node_shapes[x])
  vertex_shapes_stage2 <- sapply(membership_stage2, function(x) node_shapes[x])
  
  # Extract vertex names
  vertex_names <- unique(c(V(graph_stage1)$name, V(graph_stage2)$name))
  
  # Assign unique numbers to vertex names
  vertex_numbers <- seq_along(vertex_names)
  
  # Set vertex labels explicitly to vertex numbers
  V(graph_stage1)$label <- vertex_numbers[match(V(graph_stage1)$name, vertex_names)]
  V(graph_stage2)$label <- vertex_numbers[match(V(graph_stage2)$name, vertex_names)]
  
  # Create a legend to map numbers back to vertex names
  legend_text <- paste(vertex_numbers, vertex_names, sep = ": ")
  
  # Open a PDF device to save the plots
  pdf(paste("networks_with_legend_", tissue_stage_name, ".pdf", sep = ""))
  
  # Set up the plotting layout with smaller upper network and adjusted bottom margin
  set.seed(111213)
  layout(matrix(c(1, 2), nrow = 1, ncol = 2, byrow = FALSE), heights = c(1, 1), widths = c(3, 1))
  
  # Plot the graph for stage 1 with communities
  plot(communities_stage1, graph_stage1, 
       layout = coords_1, 
       edge.color = adjustcolor("#8d8181", alpha.f = 0.2),  # Adjust edge color with transparency
       vertex.label.cex = 0.7,
       vertex.size = 10,  # Set smaller node size
       vertex.shape = vertex_shapes_stage1)              # Set node shapes
  
  # Plot the legend on the right
  par(mar = c(0, 0, 0, 2))
  plot(0, 0, type = "n", axes = FALSE, xlab = "", ylab = "")
  legend("center", legend = legend_text, title = "Vertices", bty = "n", cex = 0.4, text.font = 0.7, horiz = FALSE)
  
  # Plot the graph for stage 2 with communities
  plot(communities_stage2, graph_stage2, 
       layout = coords_2, 
       edge.color = adjustcolor("#8d8181", alpha.f = 0.2),  # Adjust edge color with transparency
       vertex.label.cex = 0.7,
       vertex.size = 10,  # Set smaller node size
       vertex.shape = vertex_shapes_stage2)              # Set node shapes
  
  # Plot the legend on the right
  par(mar = c(0, 0, 0, 2))
  plot(0, 0, type = "n", axes = FALSE, xlab = "", ylab = "")
  legend("center", legend = legend_text, title = "Vertices", bty = "n", cex = 0.4, text.font = 0.7, horiz = FALSE)
  
  # Close the PDF device
  dev.off()
}

# Iterate over tissues and visualize networks
for (tissue_name in names(tissue_data)) {
  adj_matrices <- tissue_data[[tissue_name]]
  
  if (tissue_name == "endoderm") {
    # Liver networks
    visualize_networks(adj_matrices[["L30Endata"]], adj_matrices[["L70Endata"]], paste(tissue_name, "Liver", sep = "_"))
    
    # Lung networks
    visualize_networks(adj_matrices[["Lu70Endata"]], adj_matrices[["LuNBEndata"]], paste(tissue_name, "Lung", sep = "_"))
  } else {
    visualize_networks(adj_matrices[[1]], adj_matrices[[2]], tissue_name)
  }
}

