#============
# load required libraries #
library(dplyr)
library(Seurat)
library(patchwork)
library(data.table)
#install.packages("Seurat")

#============
# Load the cellchat.data dataset
cellchat.data <- read.table("cellchatdata.txt")
#transpose data frame
cellchat.data  <- transpose(cellchat.data )

# data preparation 
#========================
rownames(cellchat.data) <- cellchat.data$V1

cellchat.data<- cellchat.data[,-1]

colnames(cellchat.data) <- cellchat.data[1,]

cellchat.data <- cellchat.data[-1,]

orgdata.data <- cellchat.data

#drop 3 metadata columns
cellchat.data <- cellchat.data[,-1] # tumor
cellchat.data <- cellchat.data[,-1] #malignant
cellchat.data <- cellchat.data[,-1] #non-malignant

# Initialize the Seurat object with the raw (non-normalized data).
cellchat.suerat<- CreateSeuratObject(counts =  t(cellchat.data), project = "Proj4k", min.cells = 3, min.features = 10)
cellchat.suerat$nCount_RNA

# Transform the data 
cellchat.suerat <- log2(cellchat.suerat$nCount_RNA + 1)
View(cellchat.suerat)
