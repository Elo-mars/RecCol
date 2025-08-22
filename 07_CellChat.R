# example of cell-cell communications in the rectum HV
# R file

library("NMF")
library("circlize")
library("ComplexHeatmap")
library("CellChat")
library("patchwork")
library("Seurat")
library("future")
options(stringsAsFactors = FALSE)

PROJECT="/path/Figure_IBS"
setwd(PROJECT)
seuratObj=readRDS("Figure_IBS.rds")  

Idents(seuratObj)="Subtype"
table(seuratObj@active.ident)
sub=subset(seuratObj,idents="HV")
seuratObj=sub

dir.create("CellChat")
setwd(paste0(PROJECT,"/CellChat"))

Idents(seuratObj)= "Cluster"
# no space or special character allowed
seuratObj@meta.data$Cluster <- gsub(' ', '_', seuratObj@meta.data$Cluster)
seuratObj@meta.data$Cluster <- gsub('\\+', '_Pos', seuratObj@meta.data$Cluster)
seuratObj@meta.data$Cluster <- gsub('-', '_', seuratObj@meta.data$Cluster)
Idents(seuratObj)= "Cluster"
table(seuratObj@active.ident)

# need a column with patient IDs named "samples"
head(seuratObj@meta.data,2)
names(seuratObj@meta.data)[names(seuratObj@meta.data) == "orig.ident"] <- "samples"
head(seuratObj@meta.data,2)

# Part I: Data input & processing and initialization of CellChat object
# For the gene expression data matrix, genes should be in rows with rownames and cells in columns with colnames. Normalized data (e.g., library-size normalization and then log-transformed with a pseudocount of 1) is required as input for CellChat analysis. If user provides count data, we provide a normalizeData function to account for library size and then do log-transformed. For the cell group information, a dataframe with rownames is required as input for CellChat.In addition to taking a count data matrix as an input, we also provide instructions for how to prepare CellChat input files from other existing single-cell analysis toolkits, including Seurat, SingleCellExperiment and Scanpy. Please start to prepare the input data by following option A when the normalized count data and meta data are available, option B when the Seurat object is available, option C when the SingleCellExperiment object is available, and option D when the Anndata object is available. See details in the tutorial on Interface_with_other_single-cell_analysis_toolkits.

# Seurat 4 : @
data.input <- seuratObj[["RNA"]]@data # normalized data matrix
# For Seurat version >= “5.0.0”,: $
# get the normalized data via `seurat_object[["RNA"]]$data`
#data.input <- seuratObj[["RNA"]]$data # normalized data matrix

Idents(seuratObj)="Cluster"
labels <- Idents(seuratObj)
meta <- data.frame(labels = labels, row.names = names(labels)) # create a dataframe of the cell labels
head(meta)
cellChat <- createCellChat(object = seuratObj, group.by = "ident", assay = "RNA")

# Set the ligand-receptor interaction database
# Our database CellChatDB is a manually curated database of literature-supported ligand-receptor interactions in both human and mouse. CellChatDB v2 contains ~3,300 validated molecular interactions, including ~40% of secrete autocrine/paracrine signaling interactions, ~17% of extracellular matrix (ECM)-receptor interactions, ~13% of cell-cell contact interactions and ~30% non-protein signaling. Compared to CellChatDB v1, CellChatDB v2 adds more than 1000 protein and non-protein interactions such as metabolic and synaptic signaling. It should be noted that for molecules that are not directly related to genes measured in scRNA-seq, CellChat v2 estimates the expression of ligands and receptors using those molecules’ key mediators or enzymes for potential communication mediated by non-proteins.CellChatDB v2 also adds additional functional annotations of ligand-receptor pairs, such as UniProtKB keywords (including biological process, molecular function, functional class, disease, etc), subcellular location and relevance to neurotransmitter.
# When analyzing human samples, use the database CellChatDB.human; when analyzing mouse samples, use the database CellChatDB.mouse. CellChatDB categorizes ligand-receptor pairs into different types, including “Secreted Signaling”, “ECM-Receptor”, “Cell-Cell Contact” and “Non-protein Signaling”. By default, the “Non-protein Signaling” are not used.

CellChatDB <- CellChatDB.human # use CellChatDB.mouse if running on mouse data
pdf("CellChatDB.pdf", height=5, width=8)
showDatabaseCategory(CellChatDB)
dev.off()

# Show the structure of the database
dplyr::glimpse(CellChatDB$interaction)

# set the used database in the object
CellChatDB.use <- CellChatDB #  I want to use them all
cellChat@DB <- CellChatDB.use

# Preprocessing the expression data for cell-cell communication analysis
# To infer the cell state-specific communications, CellChat identifies over-expressed ligands or receptors in one cell group and then identifies over-expressed ligand-receptor interactions if either ligand or receptor are over-expressed.
# We also provide a function to project gene expression data onto protein-protein interaction (PPI) network. Specifically, a diffusion process is used to smooth genes’ expression values based on their neighbors’ defined in a high-confidence experimentally validated protein-protein network. This function is useful when analyzing single-cell data with shallow sequencing depth because the projection reduces the dropout effects of signaling genes, in particular for possible zero expression of subunits of ligands/receptors. One might be concerned about the possible artifact introduced by this diffusion process, however, it will only introduce very weak communications. By default CellChat uses the raw data (i.e., object@data.signaling) instead of the projected data. To use the projected data, users should run the function projectData before running computeCommunProb, and then set raw.use = FALSE when running computeCommunProb.

cellchat=cellChat
# subset the expression data of signaling genes for saving computation cost
cellchat <- subsetData(cellchat) # This step is necessary even if using the whole database
future::plan("multisession", workers = 60) # do parallel
options(future.globals.maxSize = 900000 * 1024^2) 

cellchat <- identifyOverExpressedGenes(cellchat,do.fast = FALSE)
cellchat <- identifyOverExpressedInteractions(cellchat)

# Part II: Inference of cell-cell communication network
# CellChat infers the biologically significant cell-cell communication by assigning each interaction with a probability value and peforming a permutation test. CellChat models the probability of cell-cell communication by integrating gene expression with prior known knowledge of the interactions between signaling ligands, receptors and their cofactors using the law of mass action.
# CAUTION: The number of inferred ligand-receptor pairs clearly depends on the method for calculating the average gene expression per cell group. By default, CellChat uses a statistically robust mean method called ‘trimean’, which produces fewer interactions than other methods. However, we find that CellChat performs well at predicting stronger interactions, which is very helpful for narrowing down on interactions for further experimental validations. In computeCommunProb, we provide an option for using other methods, such as 5% and 10% truncated mean, to calculating the average gene expression. Of note, ‘trimean’ approximates 25% truncated mean, implying that the average gene expression is zero if the percent of expressed cells in one group is less than 25%. To use 10% truncated mean, USER can set type = "truncatedMean" and trim = 0.1. To determine a proper value of trim, CellChat provides a function computeAveExpr, which can help to check the average expression of signaling genes of interest, e.g, computeAveExpr(cellchat, features = c("CXCL12","CXCR4"), type =  "truncatedMean", trim = 0.1). Therefore, if well-known signaling pathways in the studied biological process are not predicted, users can try truncatedMean with lower values of trim to change the method for calculating the average gene expression per cell group.
# When analyzing unsorted single-cell transcriptomes, under the assumption that abundant cell populations tend to send collectively stronger signals than the rare cell populations, CellChat can also consider the effect of cell proportion in each cell group in the probability calculation. USER can set population.size = TRUE.
# Compute the communication probability and infer cellular communication network
# type = "truncatedMean" and trim = 0.1
# population.size	
# whether consider the proportion of cells in each group across all sequenced cells. Set population.size = FALSE if analyzing sorting-enriched single cells, to remove the potential artifact of population size. Set population.size = TRUE if analyzing unsorted single-cell transcriptomes, with the reason that abundant cell populations tend to send collectively stronger signals than the rare cell populations.

cellchat <- computeCommunProb(cellchat, type = "truncatedMean", trim = 0.1,population.size = TRUE, nboot=20) # allows permutations
cellchat <- filterCommunication(cellchat, min.cells = 5) # standard=10
cellchat

df.net <- subsetCommunication(cellchat) 
head(df.net)
write.table(df.net, "dfnet.txt",col.names=NA, sep="\t")
library(readr)
write_csv(df.net, "dfnet.csv")

# returns a data frame consisting of all the inferred cell-cell communications at the level of ligands/receptors. Set slot.name = "netP" to access the the inferred communications at the level of signaling pathways
df.netPathways <- subsetCommunication(cellchat, slot.name = "netP") 
head(df.netPathways)
write.table(df.netPathways, "dfnetPathways.txt",col.names=NA, sep="\t")
write_csv(df.netPathways, "dfnetPathways.csv")

cellchat <- computeCommunProbPathway(cellchat)
# Calculate the aggregated cell-cell communication network
# CellChat calculates the aggregated cell-cell communication network by counting the number of links or summarizing the communication probability. Users can also calculate the aggregated network among a subset of cell groups by setting sources.use and targets.use.
cellchat <- aggregateNet(cellchat)
ptm = Sys.time()
execution.time = Sys.time() - ptm
print(as.numeric(execution.time, units = "secs"))

# CellChat can also visualize the aggregated cell-cell communication network. For example, showing the number of interactions or the total interaction strength (weights) between any two cell groups using circle plot.
groupSize <- as.numeric(table(cellchat@idents))
pdf("All_interactions.pdf", height=8, width=14)
par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")
dev.off()

# Due to the complicated cell-cell communication network, we can examine the signaling sent from each cell group. Here we also control the parameter edge.weight.max so that we can compare edge weights between differet networks.
pdf("Interactions_per_cluster.pdf", height=20, width=25)
mat <- cellchat@net$weight
par(mfrow = c(3,4), xpd=TRUE)
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])
}
dev.off()

#Part III: Visualization of cell-cell communication network
saveRDS(cellchat, "cellchat_2.rds")
# All the signaling pathways showing significant communications
print(cellchat@netP$pathways)

# Hierarchy plot
# Here we define `vertex.receive` so that the left portion of the hierarchy plot shows the index of the receiver pop
# how to see the index: 
# vertex.receiver = c(19,14) # a numeric vector. # give indexes of your clusters
# need more than 1, otherwise, matrix error
# Thicker edge line indicates a stronger signal

for (i in 1:length(cellchat@netP$pathways)){
#pathways.show <- c("TGFb") 
pathways.show <- cellchat@netP$pathways[i]
print(pathways.show) 

pdf(paste0(pathways.show,"_4pages.pdf"), height=15, width=15)
par(mfrow=c(1,1)) 
print(netVisual_aggregate(cellchat, signaling = pathways.show, layout = "circle"))
print(netVisual_aggregate(cellchat, signaling = pathways.show, layout = "chord"))
print(netVisual_heatmap(cellchat, signaling = pathways.show, color.heatmap = "Reds"))
print(netAnalysis_contribution(cellchat, signaling = pathways.show))
dev.off()
}

# examples
# Show all the significant interactions (L-R pairs) from some cell groups to other cell groups
pdf(paste0("Bubble_plots.pdf"), height=15, width=15)
netVisual_bubble(cellchat, sources.use = c(1:5), targets.use = c(5,8), remove.isolate = FALSE) 
dev.off()

# chord
pdf(paste0("Chord_diagram.pdf"), height=15, width=15)
netVisual_chord_gene(cellchat, sources.use = c(1:5), targets.use = c(5,8), lab.cex = 0.5,legend.pos.y = 30, small.gap = 0.1,big.gap = 0.2)
dev.off()

q()
n

