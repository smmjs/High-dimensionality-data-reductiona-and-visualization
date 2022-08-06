library(dplyr)
library(Seurat)
library(patchwork)
library(tcltk)
library(data.table)
library(scales)

#input feat thresh
max_feat <- 25

#input num nodes
num_nodes <- '250'
outcome_path <- "C:\\Users\\18452\\Downloads\\Vessel_GMM_All_Data\\Seurat_Analysis\\FPE_ordered.csv"

#define Min-Max normalization function
min_max_norm <- function(x) {
  (x - min(x)) / (max(x) - min(x))
}

#set working directory
#wd <- paste("C:\\Users\\18452\\Downloads\\Vessel_GMM_All_Data\\9_Results\\cross_sample#s\\center_lines\\",num_nodes,"\\PCA",sep='')
#setwd(wd)
setwd("C:\\Users\\18452\\Downloads\\Vessel_GMM_All_Data\\Seurat_Analysis")

#file = paste(wd,"\\",num_nodes,"_CL_All_PC_Scores.txt",sep='')
#file = "All_Sig_Feats_100_nodes_PC_FPE_Ordered.csv"
#file = "All_Sig_Feats_FPE_Ordered_wo_pvals.csv"
file = "All_Sig_Feats_FPE_Ordered_w_node_ids_wo_pvals.csv"

#read.table for csv files
#my_data <- read.delim(file,header=TRUE)
my_data <- read.csv(file,header=TRUE)
case_labels <- my_data$Case_ID
my_data <- subset(my_data,select=-c(Case_ID))
my_data <- my_data[,1:max_feat]


#my_data <- my_data[,10:45]
#my_data <- subset(my_data,select=-c(PC24,PC36))
#my_data <- subset(my_data,select=c(X75_12_12,X100_21_13,X50_8_12,X25_14_1,X125_11_13,
#                                   X125_6_13,X125_8_13,X125_61_12,X125_54_13,X50_30_5,
#                                   X100_54_13,X125_9_13,X25_15_1,X100_92_13,X125_7_13#))
my_data <- subset(my_data,select=c(X75_12_12,X100_21_13,X50_8_12,X25_14_1,X125_11_13))

my_feats <- colnames(my_data)
#my_data <-transpose(my_data)

#load outcome labels
outcome_labels <- read.csv(outcome_path,header=TRUE)

#reset outcome labels
outcome_labels[outcome_labels=='1'] <- 'FPE'
outcome_labels[outcome_labels=='0'] <- 'Failure'

#remove feature names from dataframe and set as rownames
my_data <- as.data.frame(t(my_data))
#row.names(my_data) <- my_feats
feat_names <- row.names(my_data)
row.names(my_data)<-feat_names

#normalize data
#my_data <- as.data.frame(t(my_data))
my_data <- as.data.frame(lapply(my_data,min_max_norm))
#my_data <- as.data.frame(t(my_data))

#create seurat object
seuratObject <- CreateSeuratObject(counts=data.matrix(my_data,rownames.force=TRUE),project='MT_Outcome',min.cells=25,min.features=1)

#identify variable features 
seuratObject <- FindVariableFeatures(seuratObject,selection.method='vst',nfeatures=)

#add labels to seurat object metadata
labels_df <- outcome_labels
rownames(labels_df) <- case_labels
seuratObject <- AddMetaData(seuratObject,labels_df)


all.genes <- rownames(seuratObject)
seuratObject <- ScaleData(seuratObject,features=all.genes)

seuratObject <- RunPCA(seuratObject,features=VariableFeatures(object=seuratObject))

VizDimLoadings(seuratObject, dims = 1:2, reduction = "pca")

seuratObject <- FindNeighbors(seuratObject, dims = 1:2)
seuratObject <- FindClusters(seuratObject, resolution = 0.5)

seuratObject <- RunUMAP(seuratObject, dims = 1:2)
my_plot <- DimPlot(seuratObject, reduction = "umap",pt.size=1)
my_plot + plot_annotation(title = '                         Unsupervised Output')

#extracting cluster identities
cids=Idents(seuratObject)

#reset labels to gender
Idents(seuratObject)=labels_df[,1]
my_plot <- DimPlot(seuratObject, reduction = "umap",pt.size=1)
my_plot + plot_annotation(title = '       Unsupervised Output - FPE Superimposed')

#extract specific features
FeaturePlot(seuratObject,features = '3')

VlnPlot(seuratObject, features = c("1","2","3","4","5"),cols=c("green","indianred"),idents=c("FPE","Failure"),pt.size=0.1)
