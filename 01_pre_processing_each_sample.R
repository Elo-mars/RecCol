# R file

SAMPLE="GC110700_SI-GA-D1" # sample name of the directory where the data is
ORGANISM="Human" # or "Mouse"
name="GC110700_SI-GA-D1" # nice name you want to see in your plots
reslist<-c(0.4,0.5,0.6,0.7,0.8,0.9,1,1.1,1.2,1.3,1.4)
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

testDir<-paste0("/path/",SAMPLE,"/outs")
sc = load10X(testDir)
sc = autoEstCont(sc)
out = adjustCounts(sc)
sample= CreateSeuratObject(out)
# remove the genes expressed less than 3 cells
gene_sum<-Matrix::rowSums(sample@assays$RNA@data>=3)
gene_keep_index<-(which(gene_sum>0))
count.data <- GetAssayData(object = sample[["RNA"]])[rownames(sample)[gene_keep_index],]    
sample <- SetAssayData(
  object = subset(sample,features=rownames(sample)[gene_keep_index]),  
  new.data = count.data, 
  assay = "RNA"  
) 

# ============================================== #
# 1) QC and selecting cells for further analysis #
# ============================================== #

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

y = sample@meta.data %>% group_by(orig.ident) %>%
  filter(nFeature_RNA > -Inf, nFeature_RNA < +Inf, percent.mt > -Inf, percent.mt < +Inf, nCount_RNA > -Inf, nCount_RNA < +Inf) %>%
  summarise(n=n())
write.table(y, paste0(opath,name,"_beforeQC_numberofcells.txt"),col.names=NA, sep="\t")

x = sample@meta.data %>% group_by(orig.ident) %>%
  filter(nFeature_RNA > 400, percent.mt > -Inf, percent.mt < 40, nCount_RNA > 400, nCount_RNA < 50000) %>%
  summarise(n=n())
x
write.table(x, paste0(opath,name,"_QC_numberofcells.txt"),col.names=NA, sep="\t")

sample <- subset(sample, subset = nFeature_RNA > 400 & nFeature_RNA < 6000  & percent.mt < 40 & nCount_RNA > 400 & nCount_RNA < 50000)

pdf(paste0(opath,name,"_filtered_VlnPlot.pdf"), width = 15, height = 8)
VlnPlot(sample, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), pt.size = 0.5,ncol = 3)
dev.off()
sample

# ================= #
# cell cycle effect #
# ================= #

sample<- NormalizeData(sample, verbose = FALSE)
gene.use=setdiff(rownames(sample),c(cc.genes.updated.2019$s.genes,cc.genes.updated.2019$g2m.genes))

FHF <- CellCycleScoring(sample, s.features = cc.genes.updated.2019$s.genes, g2m.features = cc.genes.updated.2019$g2m.genes, set.ident = TRUE)
count.data <- GetAssayData(object = FHF[["RNA"]])[gene.use,]                                                                                                                                                                      

FHF2 <- SetAssayData(
  object = subset(FHF,features=gene.use),  
  new.data = count.data, 
  assay = "RNA"  
)  

sample<-FHF2

# ==================== #
# Normalizing the data #
# ==================== #

Norm_sample <- NormalizeData(sample, normalization.method = "LogNormalize", scale.factor = 10000)

# ============================================================== #
# Identification of highly variable features (feature selection) #
# ============================================================== #

FHF <- FindVariableFeatures(Norm_sample, selection.method = "vst", nfeatures = 2000)
FHF
# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(FHF), 10)

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

Scaled <- ScaleData(FHF, features = gene.use, vars.to.regress = c("nCount_RNA", "percent.mt","S.Score", "G2M.Score"))

# ======= #
# Run PCA #
# ======= #

# library(future)
# plan("multiprocess", workers = 60)
# options(future.globals.maxSize = 10000 * 1024^2)

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

q()
n


