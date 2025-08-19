# R file
SAMPLE="GC110700_SI-GA-D1" # sample name of the directory where the data is
ORGANISM="Human" # or "Mouse"
name="GC110700_SI-GA-D1" # nice name you want to see in your plots

# load packages:
library("Seurat")                                                   
library("plyr")                                                                                                                         
library("dplyr")                                                                                                                        
library("future")                                    
library("ggplot2")                                                  
library("cowplot")                                                                                                                      
library("grid")                                                                                                                                                                                                                                                                  
library("gridExtra")                                                                                                                    
library("gtable")                                                                                                                       
library("readr")                                                                                                                        
library("beanplot")                                                                                                                     
library("gplots")                                                                                                                       
library("RColorBrewer")                                                                                                                 
library("devtools")                                                 
library("SoupX")
library("DoubletFinder")

mainDir<-"/staging/leuven/stg_00075/Project/211018_RecCol/Result_IM"
subDir<-paste0("Result_",SAMPLE)
dir.create(file.path(mainDir, subDir), showWarnings = FALSE)
opath_s<-paste0(mainDir,"/",subDir,"/")
opath<-opath_s

# ===== #
# SoupX #
# ===== #

testDir<-paste0("/staging/leuven/stg_00075/Project/211018_RecCol/Raw_Data/",SAMPLE,"/outs")
# testDir<-"/home/rstudio/host/UCMB_course/Hemsley_ECC/scRNA-seq/GC108011_H1/outs"
sc = load10X(testDir)
sc = autoEstCont(sc)
out = adjustCounts(sc)
sample= CreateSeuratObject(out)
saveRDS(sample, file =paste0(opath,name,"_After_SoupX.rds"))
# dim(sample)
### remove the genes expressed less than 3 cells
gene_sum<-Matrix::rowSums(sample@assays$RNA@data>=3)
gene_keep_index<-(which(gene_sum>0))
count.data <- GetAssayData(object = sample[["RNA"]])[rownames(sample)[gene_keep_index],]    
sample <- SetAssayData(
  object = subset(sample,features=rownames(sample)[gene_keep_index]),  
  new.data = count.data, 
  assay = "RNA"  
) 


saveRDS(sample, file =paste0(opath,name,"_After_SoupX_3cellgenes.rds"))
ggsave(file=paste0(opath,name,"Contamination_fraction.png"))



# ============================================== #
# 1) QC and selecting cells for further analysis #
# ============================================== #

# The number of unique genes detected in each cell.
# Low-quality cells or empty droplets will often have very few genes
# Cell doublets or multiplets may exhibit an aberrantly high gene count
# Similarly, the total number of molecules detected within a cell (correlates strongly with unique genes)
# The percentage of reads that map to the mitochondrial genome
# Low-quality / dying cells often exhibit extensive mitochondrial contamination
# We calculate mitochondrial QC metrics with the PercentageFeatureSet function, which calculates the percentage of counts originating from a set of features
# We use the set of all genes starting with MT- as a set of mitochondrial genes

# PCA
if(ORGANISM=="Mouse"){
  sample[["percent.mt"]] <- PercentageFeatureSet(sample, pattern = "^mt-")
} else {
  sample[["percent.mt"]] <- PercentageFeatureSet(sample, pattern = "^MT-")
}

pdf(paste0(opath,name,"_VlnPlot.pdf"), width = 15, height = 8)
VlnPlot(sample, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), pt.size = 0.5,ncol = 3)
dev.off()

# ====== #
# Filter #
# ====== #

# PCA
y = sample@meta.data %>% group_by(orig.ident) %>%
  filter(nFeature_RNA > -Inf, nFeature_RNA < +Inf, percent.mt > -Inf, percent.mt < +Inf, nCount_RNA > -Inf, nCount_RNA < +Inf) %>%
  summarise(n=n())
write.table(y, paste0(opath,name,"_beforeQC_numberofcells.txt"),col.names=NA, sep="\t")

x = sample@meta.data %>% group_by(orig.ident) %>%
  filter(nFeature_RNA > 400, percent.mt > -Inf, percent.mt < 40, nCount_RNA > 400, nCount_RNA < 50000) %>%
  summarise(n=n())
x
write.table(x, paste0(opath,name,"_QC_numberofcells.txt"),col.names=NA, sep="\t")

# features are genes and count_RNA are UMIs
sample <- subset(sample, subset = nFeature_RNA > 400 & nFeature_RNA < 6000  & percent.mt < 40 & nCount_RNA > 400 & nCount_RNA < 50000)

# sort(Matrix::rowSums(experiment.aggregate@data>=2))

pdf(paste0(opath,name,"_filtered_VlnPlot.pdf"), width = 15, height = 8)
VlnPlot(sample, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), pt.size = 0.5,ncol = 3)
dev.off()
sample


# ========================= #
# cell cycle effect(Elodie method)         #
# ========================= #
# remove gene cycle genes: 
sample<- NormalizeData(sample, verbose = FALSE)
gene.use=setdiff(rownames(sample),c(cc.genes.updated.2019$s.genes,cc.genes.updated.2019$g2m.genes))

FHF <- CellCycleScoring(sample, s.features = cc.genes.updated.2019$s.genes, g2m.features = cc.genes.updated.2019$g2m.genes, set.ident = TRUE)
#  install.packages("SeuratObject")                                                                                                                                                                         

count.data <- GetAssayData(object = FHF[["RNA"]])[gene.use,]                                                                                                                                                                      

# count.data <- as.matrix(x = count.data + 1)                                                                                                                                                                                       

FHF2 <- SetAssayData(
  object = subset(FHF,features=gene.use),  
  new.data = count.data, 
  assay = "RNA"  
)  

sample<-FHF2



# ==================== #
# Normalizing the data #
# ==================== #

# After removing unwanted cells from the dataset, the next step is to normalize the data.
# # PCA
Norm_sample <- NormalizeData(sample, normalization.method = "LogNormalize", scale.factor = 10000)

# ============================================================== #
# Identification of highly variable features (feature selection) #
# ============================================================== #

#PCA
FHF <- FindVariableFeatures(Norm_sample, selection.method = "vst", nfeatures = 2000)
FHF
# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(FHF), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(FHF)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
pdf(paste0(opath,name,"_Highly_Variable_genes.pdf"), width=10, height=5)
CombinePlots(plots = list(plot1, plot2))
dev.off()

all.genes=rownames(FHF)
length(all.genes)

# ===== #
# Scale #
# ===== #
saveRDS(FHF, file =paste0(opath,name,"_B4_regression.rds"))
# Scaled <- ScaleData(FHF, features = all.genes)
#saveRDS(Scaled, file =paste0(opath,name,"_Scaled.rds"))
Scaled <- ScaleData(FHF, features = gene.use, vars.to.regress = c("nCount_RNA", "percent.mt","S.Score", "G2M.Score")) # regress on cell cycle scores after remove them 

saveRDS(Scaled, file =paste0(opath,name,"_AFT_regression.rds"))




# ===== #
# Scale #
# ===== #

# Scaled <- ScaleData(FHF, features = all.genes)
#saveRDS(Scaled, file =paste0(opath,name,"_Scaled.rds"))
#Scaled <- ScaleData(FHF, features = all.genes, vars.to.regress = c("nCount_RNA", "percent.mt","S.Score", "G2M.Score")) # regress on cell cycle scores? or remove ? 

# ======= #
# Run PCA #
# ======= #

# library(future)
# plan("multiprocess", workers = 60) # 60 cores/CPU
# options(future.globals.maxSize = 10000 * 1024^2) # 10G RAM
#
Scaled <- RunPCA(Scaled, features = VariableFeatures(Scaled), npcs = 30, ndims.print = 1:5, nfeatures.print = 30)

pdf(paste0(opath,name,"_VizDim.pdf"))
VizDimLoadings(Scaled, dims = 1:2, reduction = "pca")
dev.off()

pdf(paste0(opath,name,"_Dimplot.pdf"))
DimPlot(Scaled, reduction = "pca")
dev.off()

pdf(paste0(opath,name,"_Dimheatmap.pdf"))
DimHeatmap(Scaled, dims = 1:30, cells = 500, balanced = TRUE)
dev.off()

# ======================================================== #
# Determine statistically significant principal components #
# ======================================================== #

# plan("multiprocess", workers = 60) # uses 45 CPU
# options(future.globals.maxSize = 10000 * 1024^2) # set it at 10Go of RAM)
Scaled <- JackStraw(Scaled, num.replicate = 100, dims = 30)
Scaled <- ScoreJackStraw(Scaled, dims = 1:30)

pdf(paste0(opath,name,"_JackStraw.pdf"))
JackStrawPlot(Scaled, dims = 1:30)
dev.off()
pdf(paste0(opath,name,"_Elbow.pdf"))
ElbowPlot(Scaled, ndims = 30, reduction = "pca")
dev.off()

NPC=30

# ==================== #
# 3) Cluster the cells #
# ==================== #
Scaled <- RunUMAP(Scaled, reduction = "pca", dims = 1:30)
FN <- FindNeighbors(Scaled, reduction = "pca", dims = 1:NPC, k.param = 30, compute.SNN = TRUE, nn.method = "rann", nn.eps = 0)
# plan("multiprocess", workers = 60) # uses 45 CPU
# options(future.globals.maxSize = 10000 * 1024^2) # set it at 10Go of RAM)

FC <- FindClusters(FN, resolution = reslist)

saveRDS(FC, paste0(opath,name,"_FC.rds"))
head(FC@meta.data)

# ================================================ #
# Run non-linear dimensional reduction (UMAP/tSNE) #
# ================================================ #

Idents(FC)

# Color by sample
Seurat_per_sample_Ident <- SetIdent(object = FC, value = 'orig.ident')
Seurat_per_sample <- RunUMAP(Seurat_per_sample_Ident, dims = 1:NPC)
Seurat_per_sample

# remove "plot.title" now it is a "+ labs"
pdf(paste0(opath,name,"_Plot_per_sample.pdf"))
UMAPPlot(Seurat_per_sample, reduction = "umap", pt.size=0.1,label=FALSE) +labs(title = "UMAP per sample")+ theme(plot.title = element_text(hjust = 0.5))
dev.off()

pdf(paste0(opath,name,"_Plot_per_sample_nolegend.pdf"))
UMAPPlot(Seurat_per_sample, reduction = "umap", pt.size=0.1,label=FALSE) +labs(title = "UMAP per sample")+ theme(plot.title = element_text(hjust = 0.5)) +NoLegend()
dev.off()

head(FC@meta.data)


pdf(paste0(opath,name,"_Clusters_diff_res.pdf"))
plotlist <- list()
for (i in (1:length(reslist))){
  res_i<-paste0('RNA_snn_res.',reslist[i])
  Seurat <- SetIdent(object = FC, value =res_i)
  Seurat <- RunUMAP(Seurat, dims = 1:NPC)
  # UMAPPlot(Seurat, reduction = "umap", pt.size=0.25,label=T) +labs(title = paste0("UMAP with res ",reslist[i]))+ theme(plot.title = element_text(hjust = 0.5))
  q<-DimPlot(Seurat,reduction = 'umap')+labs(title = paste0("UMAP with res ",reslist[i]))+ theme(plot.title = element_text(hjust = 0.5))
  plotlist[[i]]<-q
  print(plotlist[[i]])}
dev.off()


head(FC@meta.data)
u=FC@meta.data
write.table(u, paste0(opath,name,"_metadata.txt"),col.names=NA, sep="\t")


# tables - number of cells per cluster for each resolution
for (i in (1:length(reslist))){
  res_i<-paste0('RNA_snn_res.',reslist[i])
  
  Seurat_reso <- SetIdent(object = FC, value = res_i)
  cell_perC=table(Seurat_reso@active.ident)
  cell_perC
  write.table(cell_perC, paste0(opath,name,"Integrated_","_res0",reslist[i],"_cell_per_c.txt"),col.names=NA, sep="\t")
}


# FAM for each resolution
# with heatmaps with different #genes/cluster
FC<-readRDS(file=paste0(opath,name,'_FC.rds'))
library(future)
for (i in (1:length(reslist))){
  res_i<-paste0('RNA_snn_res.',reslist[i])
  
  Seurat_reso <- SetIdent(object = FC, value = res_i)
  Seurat.markers<-FindAllMarkers(Seurat_reso, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
  # Seurat.markers <- FindAllMarkers(Seurat_reso,test.use="MAST",  min.pct = 0, logfc.threshold = 0, max.cells.per.ident = Inf, latent.vars = NULL,min.cells.feature = 0, min.cells.group = 0)
  write.table(Seurat.markers, paste0(opath,name,"_res0",reslist[i],"_MAST_markers_unannotated_Clusters.txt"),col.names=NA, sep="\t")
  tsne.markers=Seurat.markers
  #DefaultAssay(Seurat_reso) <- 'integrated'
  topn = tsne.markers %>% group_by(cluster) %>% top_n(10, avg_log2FC)
  topn$gene
  pdf(paste0(opath,name,"_res0",reslist[i]*10,"_heatmap_top10genes.pdf"),width=25, height=20)
  q<-DoHeatmap(object = Seurat_reso, features=topn$gene, group.by = "ident", label = TRUE, angle=45) + ggtitle("top 10 genes expressed per cluster") + NoLegend()
  print(q)
  dev.off()
  topn = tsne.markers %>% group_by(cluster) %>% top_n(15, avg_log2FC)
  topn$gene
  pdf(paste0(opath,name,"_res0",reslist[i]*10,"_heatmap_top15genes.pdf"),width=25, height=25)
  q<-DoHeatmap(object = Seurat_reso, features=topn$gene, group.by = "ident", label = TRUE, angle=45) + ggtitle("top 15 genes expressed per cluster") + NoLegend()
  print(q)
  dev.off()
  topn = tsne.markers %>% group_by(cluster) %>% top_n(20, avg_log2FC)
  topn$gene
  pdf(paste0(opath,name,"_res0",reslist[i]*10,"_heatmap_top20genes.pdf"),width=25, height=30)
  q<-DoHeatmap(object = Seurat_reso, features=topn$gene, group.by = "ident", label = TRUE, angle=45) + ggtitle("top 20 genes expressed per cluster") + NoLegend()
  print(q)
  dev.off()
  topn = tsne.markers %>% group_by(cluster) %>% top_n(5, avg_log2FC)
  topn$gene
  pdf(paste0(opath,name,"_res0",reslist[i]*10,"_heatmap_top5genes.pdf"),width=25, height=15)
  q<-DoHeatmap(object = Seurat_reso, features=topn$gene, group.by = "ident", label = TRUE, angle=45) + ggtitle("top 5 genes expressed per cluster") + NoLegend()
  print(q)
  dev.off()
  topn = tsne.markers %>% group_by(cluster) %>% top_n(5, avg_log2FC)
  topn$gene
  topn
  write_csv(topn, paste0(opath,name,"_res0",reslist[i]*10,"_top5genespercluster_logFC.csv"))
  pdf(paste0(opath,name,"_res0",reslist[i]*10,"_VlnPlot.pdf"), width = 15, height = 8)
  VlnPlot(Seurat_reso, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), pt.size = 0.5,ncol = 3)
  dev.off()  
}
q()
n


