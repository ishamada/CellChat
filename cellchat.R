# run this command on ternimal
# conda install gfortran_linux-64
# install.packages('igraph','1.3.5')

# # Load the methods package
# library(methods)

library(data.table)

# log(as.matrix(CellChat::normalizeData(testobject@assays$RNA[1:5,1:2])) + 1)

#load(url("https://ndownloader.figshare.com/files/25950872")) # This is a combined data from two biological conditions: normal and diseases
#data.input = data_humanSkin$data # normalized data matrix


library(doParallel)
detectCores()
registerDoParallel(cores = 4)

future::plan("multiprocess", workers = 4) # do parallel

#Sys.sleep(2)

# Load the scdata dataset
df.data <- read.table("~/mel-data/scrdatafile.txt")

#Sys.sleep(2)
#transpose data frame
df.data <- transpose(df.data)

#Sys.sleep(2)

# data preparation 
#========================
rownames(df.data) <- df.data$V1

df.data <- df.data[,-1]

colnames(df.data) <- df.data[1,]

df.data <- df.data[-1,]

orgdata.data <- df.data
#========================
# unique(df.data$tumor)

# subset tumor number 79
mel79.data <- df.data[which(df.data[,1] == "79"),]

# which(mel79.data$`malignant(1=no,2=yes,0=unresolved)` == 2)

# set 7 to malignant cells
mel79.data[which(mel79.data$`malignant(1=no,2=yes,0=unresolved)` == 2),3] = 7

# which(mel79.data$`malignant(1=no,2=yes,0=unresolved)` == 0 & mel79.data$`non-malignant cell type (1=T,2=B,3=Macro.4=Endo.,5=CAF;6=NK)` == 0)

# set 9 to unresolved cells
mel79.data[which(mel79.data$`malignant(1=no,2=yes,0=unresolved)` == 0 & mel79.data$`non-malignant cell type (1=T,2=B,3=Macro.4=Endo.,5=CAF;6=NK)` == 0) , 3] = 9

# set 8 to non.malignant cells
mel79.data[which(mel79.data$`malignant(1=no,2=yes,0=unresolved)` == 1 & mel79.data$`non-malignant cell type (1=T,2=B,3=Macro.4=Endo.,5=CAF;6=NK)` == 0) , 3] = 8

# explore all cell groups
table(mel79.data$`non-malignant cell type (1=T,2=B,3=Macro.4=Endo.,5=CAF;6=NK)`)


# drop 2 columns

mel79.data <- mel79.data[,-1]
mel79.data <- mel79.data[,-1]


mel79.data[which(mel79.data$`non-malignant cell type (1=T,2=B,3=Macro.4=Endo.,5=CAF;6=NK)` == 1),1] = "Tcell"

mel79.data[which(mel79.data$`non-malignant cell type (1=T,2=B,3=Macro.4=Endo.,5=CAF;6=NK)` == 2),1] = "Bcell"

mel79.data[which(mel79.data$`non-malignant cell type (1=T,2=B,3=Macro.4=Endo.,5=CAF;6=NK)` == 3),1] = "Macro"

mel79.data[which(mel79.data$`non-malignant cell type (1=T,2=B,3=Macro.4=Endo.,5=CAF;6=NK)` == 4),1] = "Endo"

mel79.data[which(mel79.data$`non-malignant cell type (1=T,2=B,3=Macro.4=Endo.,5=CAF;6=NK)` == 5),1] = "CAF" 

mel79.data[which(mel79.data$`non-malignant cell type (1=T,2=B,3=Macro.4=Endo.,5=CAF;6=NK)` == 6),1] = "NK" 

mel79.data[which(mel79.data$`non-malignant cell type (1=T,2=B,3=Macro.4=Endo.,5=CAF;6=NK)` == 7),1] = "Mel79" 

mel79.data[which(mel79.data$`non-malignant cell type (1=T,2=B,3=Macro.4=Endo.,5=CAF;6=NK)` == 8),1] = "Nrml79" 

mel79.data[which(mel79.data$`non-malignant cell type (1=T,2=B,3=Macro.4=Endo.,5=CAF;6=NK)` == 9),1] = "Urslvd" 

# create meta data frame
df <- data.frame(rownames(mel79.data),  mel79.data[,1])

row.names(df) <-  df[,1]

colnames(df)[2]  <- "labels"

mel79.metadata <- subset(df, select = c("labels"))


# drop metadata column
mel79.data <- mel79.data[,-1] # normalized data matrix

# make a backup data
mel79.ready <- mel79.data

unique(mel79.metadata$labels) # check the cell labels

library(CellChat)
library(patchwork)
options(stringsAsFactors = FALSE)
library(dplyr)
library(Seurat)

# Initialize the Seurat object with the raw (non-normalized data).
seurat.mel79data <- CreateSeuratObject(counts =  t(mel79.data), project = "mel79", min.cells = 3, min.features = 10)

# seurat.mel79data <- SCTransform(seurat.mel79data,assay = "RNA",new.assay.name = "SCT")

# Create a CellChat object
cellchat <- createCellChat(object = seurat.mel79data, meta = mel79.metadata, group.by = "labels")


cellchat <- addMeta(cellchat, meta = mel79.metadata)
cellchat <- setIdent(cellchat, ident.use = "labels") # set "labels" as default cell identity
levels(cellchat@idents) # show factor levels of the cell labels
groupSize <- as.numeric(table(cellchat@idents)) # number of cells in each cell group

CellChatDB <- CellChatDB.human # use CellChatDB.mouse if running on mouse data
showDatabaseCategory(CellChatDB)

# Show the structure of the database
dplyr::glimpse(CellChatDB$interaction)

# use a subset of CellChatDB for cell-cell communication analysis
CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling") # use Secreted Signaling
# use all CellChatDB for cell-cell communication analysis
# CellChatDB.use <- CellChatDB # simply use the default CellChatDB

# set the used database in the object
cellchat@DB <- CellChatDB.use

# subset the expression data of signaling genes for saving computation cost
cellchat <- subsetData(cellchat) # This step is necessary even if using the whole database

cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)


# project gene expression data onto PPI (Optional: when running it, USER should set `raw.use = FALSE` in the function `computeCommunProb()` in order to use the projected data)
cellchat <- projectData(cellchat, PPI.human)

cellchat <- computeCommunProb(cellchat)

# Filter out the cell-cell communication if there are only few number of cells in certain cell groups
cellchat <- filterCommunication(cellchat, min.cells = 10)

df.net <- subsetCommunication(cellchat)

cellchat <- computeCommunProbPathway(cellchat)

cellchat <- aggregateNet(cellchat)

cellchat

groupSize <- as.numeric(table(cellchat@idents))
par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T,  label.edge= T, title.name = "Number of interactions")
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")

mat <- cellchat@net$weight
par(mfrow = c(3,4), xpd=TRUE)
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])
}

pathways.show <- c("CXCL") 
# Hierarchy plot
# Here we define `vertex.receive` so that the left portion of the hierarchy plot shows signaling to fibroblast and the right portion shows signaling to immune cells 
vertex.receiver = seq(1,4) # a numeric vector. 
netVisual_aggregate(cellchat, signaling = pathways.show,  vertex.receiver = vertex.receiver)
# Circle plot
par(mfrow=c(1,1))
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "circle")

# Chord diagram
par(mfrow=c(1,1))
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "chord")

# Heatmap
par(mfrow=c(1,1))
netVisual_heatmap(cellchat, signaling = pathways.show, color.heatmap = "Reds")
#> Do heatmap based on a single object

# Chord diagram
group.cellType <- c(rep("FIB", 4), rep("DC", 4), rep("TC", 4)) # grouping cell clusters into fibroblast, DC and TC cells
names(group.cellType) <- levels(cellchat@idents)
netVisual_chord_cell(cellchat, signaling = pathways.show, group = group.cellType, title.name = paste0(pathways.show, " signaling network"))
#> Plot the aggregated cell-cell communication network at the signaling pathway level

netAnalysis_contribution(cellchat, signaling = pathways.show)

pairLR.CXCL <- extractEnrichedLR(cellchat, signaling = pathways.show, geneLR.return = FALSE)
LR.show <- pairLR.CXCL[1,] # show one ligand-receptor pair
# Hierarchy plot
vertex.receiver = seq(1,4) # a numeric vector
netVisual_individual(cellchat, signaling = pathways.show,  pairLR.use = LR.show, vertex.receiver = vertex.receiver)
#> [[1]]
# Circle plot
netVisual_individual(cellchat, signaling = pathways.show, pairLR.use = LR.show, layout = "circle")

# Chord diagram
netVisual_individual(cellchat, signaling = pathways.show, pairLR.use = LR.show, layout = "chord")

# Access all the signaling pathways showing significant communications
pathways.show.all <- cellchat@netP$pathways
# check the order of cell identity to set suitable vertex.receiver
levels(cellchat@idents)
vertex.receiver = seq(1,4)
for (i in 1:length(pathways.show.all)) {
  # Visualize communication network associated with both signaling pathway and individual L-R pairs
  netVisual(cellchat, signaling = pathways.show.all[i], vertex.receiver = vertex.receiver, layout = "hierarchy")
  # Compute and visualize the contribution of each ligand-receptor pair to the overall signaling pathway
  gg <- netAnalysis_contribution(cellchat, signaling = pathways.show.all[i])
  ggsave(filename=paste0(pathways.show.all[i], "_L-R_contribution.pdf"), plot=gg, width = 3, height = 2, units = 'in', dpi = 300)
}

# show all the significant interactions (L-R pairs) from some cell groups (defined by 'sources.use') to other cell groups (defined by 'targets.use')
netVisual_bubble(cellchat, sources.use = 4, targets.use = c(5:11), remove.isolate = FALSE)
#> Comparing communications on a single object

# show all the significant interactions (L-R pairs) associated with certain signaling pathways
netVisual_bubble(cellchat, sources.use = 4, targets.use = c(5:11), signaling = c("CCL","CXCL"), remove.isolate = FALSE)
#> Comparing communications on a single object

# show all the significant interactions (L-R pairs) based on user's input (defined by `pairLR.use`)
pairLR.use <- extractEnrichedLR(cellchat, signaling = c("CCL","CXCL","FGF"))
netVisual_bubble(cellchat, sources.use = c(3,4), targets.use = c(5:8), pairLR.use = pairLR.use, remove.isolate = TRUE)
#> Comparing communications on a single object

# show all the significant interactions (L-R pairs) from some cell groups (defined by 'sources.use') to other cell groups (defined by 'targets.use')
# show all the interactions sending from Inflam.FIB
netVisual_chord_gene(cellchat, sources.use = 4, targets.use = c(5:11), lab.cex = 0.5,legend.pos.y = 30)

# show all the interactions received by Inflam.DC
netVisual_chord_gene(cellchat, sources.use = c(1,2,3,4), targets.use = 8, legend.pos.x = 15)

# show all the significant interactions (L-R pairs) associated with certain signaling pathways
netVisual_chord_gene(cellchat, sources.use = c(1,2,3,4), targets.use = c(5:11), signaling = c("CCL","CXCL"),legend.pos.x = 8)

# show all the significant signaling pathways from some cell groups (defined by 'sources.use') to other cell groups (defined by 'targets.use')
netVisual_chord_gene(cellchat, sources.use = c(1,2,3,4), targets.use = c(5:11), slot.name = "netP", legend.pos.x = 10)

plotGeneExpression(cellchat, signaling = "CXCL")

plotGeneExpression(cellchat, signaling = "CXCL", enriched.only = FALSE)

# Compute the network centrality scores
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP") # the slot 'netP' means the inferred intercellular communication network of signaling pathways
# Visualize the computed centrality scores using heatmap, allowing ready identification of major signaling roles of cell groups
netAnalysis_signalingRole_network(cellchat, signaling = pathways.show, width = 8, height = 2.5, font.size = 10)

# Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
gg1 <- netAnalysis_signalingRole_scatter(cellchat)
#> Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
# Signaling role analysis on the cell-cell communication networks of interest
gg2 <- netAnalysis_signalingRole_scatter(cellchat, signaling = c("CXCL", "CCL"))
#> Signaling role analysis on the cell-cell communication network from user's input
gg1 + gg2

# Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
ht1 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "outgoing")
ht2 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "incoming")
ht1 + ht2

# Signaling role analysis on the cell-cell communication networks of interest
ht <- netAnalysis_signalingRole_heatmap(cellchat, signaling = c("CXCL", "CCL"))

library(ggalluvial)

library(ggalluvial)

nPatterns = 3
cellchat <- identifyCommunicationPatterns(cellchat, pattern = "outgoing", k = nPatterns)

# river plot
netAnalysis_river(cellchat, pattern = "outgoing")
#> Please make sure you have load `library(ggalluvial)` when running this function

# dot plot
netAnalysis_dot(cellchat, pattern = "outgoing")

selectK(cellchat, pattern = "incoming")

nPatterns = 4
cellchat <- identifyCommunicationPatterns(cellchat, pattern = "incoming", k = nPatterns)

# river plot
netAnalysis_river(cellchat, pattern = "incoming")
#> Please make sure you have load `library(ggalluvial)` when running this function

# dot plot
netAnalysis_dot(cellchat, pattern = "incoming")

cellchat <- computeNetSimilarity(cellchat, type = "functional")
cellchat <- netEmbedding(cellchat, type = "functional")
#> Manifold learning of the signaling networks for a single dataset
cellchat <- netClustering(cellchat, type = "functional")
#> Classification learning of the signaling networks for a single dataset
# Visualization in 2D-space
netVisual_embedding(cellchat, type = "functional", label.size = 3.5)


cellchat <- computeNetSimilarity(cellchat, type = "structural")
cellchat <- netEmbedding(cellchat, type = "structural")
#> Manifold learning of the signaling networks for a single dataset
cellchat <- netClustering(cellchat, type = "structural")
#> Classification learning of the signaling networks for a single dataset
# Visualization in 2D-space
netVisual_embedding(cellchat, type = "structural", label.size = 3.5)

netVisual_embeddingZoomIn(cellchat, type = "structural", nCol = 2)

saveRDS(cellchat, file = "cellchat_humanSkin_LS.rds")