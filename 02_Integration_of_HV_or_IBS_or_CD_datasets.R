library("Seurat")
library("plyr")
library("dplyr")
library("future")
library("ggplot2")
library("cowplot")
library("tidyverse")
library("grid")
library("gridExtra")
library("gtable")
library("readr")
library("beanplot")
library("gplots") 
library("RColorBrewer")
library("devtools")
library("harmony")

sessionInfo()
project_name="Fig_IBS"
mainDir<-".../Result"
subDir<-paste0("Result_",PROJECT)
dir.create(file.path(mainDir, subDir), showWarnings = FALSE)
opath<-paste0(mainDir,"/",subDir,"/")
setwd(opath)
reslist<-c(0.4,0.5,0.6,0.7,0.8,0.9,1,1.1,1.2,1.3,1.4,1.5,1.6,3)

# ALL
IBS_list = c("GC1003333_C5","GC1003334_D5","GC1004260_D1","GC107954_H12","GC109183_SI-GA-D5","GC109184_SI-GA-E5","GC109290_SI-GA-G5","GC109291_SI-GA-H5","GC110306_SI-GA-H8","GC110307_SI-GA-A9","GC110365_SI-GA-C9","GC110366_SI-GA-D9","GC110704_SI-GA-H1","GC110705_SI-GA-A2","GC109185_SI-GA-F5","GC109691_SI-GA-D6","GC109692_SI-GA-E6","GC109902_SI-GA-F6","GC109903_SI-GA-G6","GC110698_SI-GA-B1","GC110699_SI-GA-C1","GC111007_SI-GA-G6","GC111008_SI-GA-H6","GC112065_SI-GA-G5","GC112066_SI-GA-H5", "GC108008_E1","GC108009_F1","GC108546_SI-GA-A3","GC108547_SI-GA-B3","GC110186_SI-GA-A8","GC110187_SI-GA-B8","GC110460_SI-GA-E9","GC110461_SI-GA-F9","GC110697_SI-GA-A1","GC110700_SI-GA-D1","GC110701_SI-GA-E1","GC111344_SI-GA-G7","GC111345_SI-GA-H7","GC112063_SI-GA-E5","GC112064_SI-GA-F5")   

sample_list<-list()
for (i in (1:length(IBS_list))) {
  sample=IBS_list[i]
  if (file.exists(paste0("/staging/leuven/stg_00075/Project/211018_RecCol/Result_IM/Result_",sample,"/",sample,"_FC.rds"))){
    so<-readRDS(paste0("/staging/leuven/stg_00075/Project/211018_RecCol/Result_IM/Result_",sample,"/",sample,"_FC.rds"))
    
    # Add metadata
    
    #Status: HV or IBS or CD
    #and subtype: HV, IBS-C or IBS-D or CD+IBS or CD-IBS
    if (sample %in% HV_list){
    so@meta.data$Status <- "HV"                                                                                                                                                                                                                                              
    so@meta.data$Subtype <- "HV"} 
    if (sample %in% IBS_C_list){
      so@meta.data$Status <- "IBS"                                                                                                                                                                                                                                                    so@meta.data$Subtype <- "IBS-C"} 
    if (sample %in% IBS_D_list){
      so@meta.data$Status <- "IBS"                                                                                                                                                                                                                                                   so@meta.data$Subtype <- "IBS-D"} 
    if (sample %in% IBS_list){
so@meta.data$orig.ident = toString(sample)}


so@meta.data$PIN <- ifelse(so@meta.data$orig.ident=="GC1003333_C5","PIN_1337",ifelse(so@meta.data$orig.ident=="GC1003334_D5","PIN_1337",ifelse(so@meta.data$orig.ident=="GC1004260_D1","PIN_1223",ifelse(so@meta.data$orig.ident=="GC107954_H12","PIN_1399",ifelse(so@meta.data$orig.ident=="GC108008_E1","PIN_1387",ifelse(so@meta.data$orig.ident=="GC108009_F1","PIN_1387",ifelse(so@meta.data$orig.ident=="GC108546_SI-GA-A3","PIN_1379",ifelse(so@meta.data$orig.ident=="GC108547_SI-GA-B3","PIN_1379",ifelse(so@meta.data$orig.ident=="GC109183_SI-GA-D5","PIN_1456",ifelse(so@meta.data$orig.ident=="GC109184_SI-GA-E5","PIN_1456",ifelse(so@meta.data$orig.ident=="GC109185_SI-GA-F5","PIN_1458",ifelse(so@meta.data$orig.ident=="GC109290_SI-GA-G5","PIN_1454",ifelse(so@meta.data$orig.ident=="GC109291_SI-GA-H5","PIN_1454",ifelse(so@meta.data$orig.ident=="GC109691_SI-GA-D6","PIN_1503",ifelse(so@meta.data$orig.ident=="GC109692_SI-GA-E6","PIN_1503",ifelse(so@meta.data$orig.ident=="GC109902_SI-GA-F6","PIN_1473",ifelse(so@meta.data$orig.ident=="GC109903_SI-GA-G6","PIN_1473",ifelse(so@meta.data$orig.ident=="GC110186_SI-GA-A8","PIN_1417",ifelse(so@meta.data$orig.ident=="GC110187_SI-GA-B8","PIN_1417",ifelse(so@meta.data$orig.ident=="GC110306_SI-GA-H8","PIN_1489",ifelse(so@meta.data$orig.ident=="GC110307_SI-GA-A9","PIN_1489",ifelse(so@meta.data$orig.ident=="GC110365_SI-GA-C9","PIN_1459",ifelse(so@meta.data$orig.ident=="GC110366_SI-GA-D9","PIN_1459",ifelse(so@meta.data$orig.ident=="GC110460_SI-GA-E9","PIN_1412",ifelse(so@meta.data$orig.ident=="GC110461_SI-GA-F9","PIN_1412",ifelse(so@meta.data$orig.ident=="GC110697_SI-GA-A1","PIN_1067",ifelse(so@meta.data$orig.ident=="GC110698_SI-GA-B1","PIN_1332",ifelse(so@meta.data$orig.ident=="GC110699_SI-GA-C1","PIN_1332",ifelse(so@meta.data$orig.ident=="GC110700_SI-GA-D1","PIN_1340",ifelse(so@meta.data$orig.ident=="GC110701_SI-GA-E1","PIN_1340",ifelse(so@meta.data$orig.ident=="GC110704_SI-GA-H1","PIN_1274",ifelse(so@meta.data$orig.ident=="GC110705_SI-GA-A2","PIN_1274",ifelse(so@meta.data$orig.ident=="GC111007_SI-GA-G6","PIN_1453",ifelse(so@meta.data$orig.ident=="GC111008_SI-GA-H6","PIN_1453",ifelse(so@meta.data$orig.ident=="GC111344_SI-GA-G7","PIN_1540",ifelse(so@meta.data$orig.ident=="GC111345_SI-GA-H7","PIN_1540",ifelse(so@meta.data$orig.ident=="GC112063_SI-GA-E5","PIN_1479",ifelse(so@meta.data$orig.ident=="GC112064_SI-GA-F5","PIN_1479",ifelse(so@meta.data$orig.ident=="GC112065_SI-GA-G5","PIN_1491","PIN_1491")))))))))))))))))))))))))))))))))))))))

      sample_list[[i]]<-so
}
}  
head(sample_list[[i]])


# 
# sample1="GC110527_SI-GA-E10/"
# sample2="GC110526_SI-GA-D10"
# sample3="GC1001745_E11"
# sample4="GC1001746_F11"
# # sample5=""
# #sample6=""
#sample7=""
#sample8=""

# =========== #
# Violin plot #
# =========== #

# All together for the pptx
GR <- merge(x = sample_list[[1]], 
            y = c(sample_list[[2]],sample_list[[3]],sample_list[[4]],sample_list[[5]],sample_list[[6]],sample_list[[7]],sample_list[[8]],sample_list[[9]],sample_list[[10]],sample_list[[11]],sample_list[[12]],sample_list[[13]],sample_list[[14]],sample_list[[15]],sample_list[[16]],sample_list[[17]],sample_list[[18]],sample_list[[19]],sample_list[[20]],sample_list[[21]],sample_list[[22]],sample_list[[23]],sample_list[[24]],sample_list[[25]],sample_list[[26]],sample_list[[27]],sample_list[[28]],sample_list[[29]],sample_list[[30]],sample_list[[31]],sample_list[[32]],sample_list[[33]],sample_list[[34]],sample_list[[35]],sample_list[[36]],sample_list[[37]],sample_list[[38]],sample_list[[39]],sample_list[[40]]),add.cell.ids = IBS_list, project = project_name)


GR
head(GR@meta.data)

GR <- SetIdent(object = GR, value = 'orig.ident')

# ===================== #
# filter more if needed # 
# ===================== #


# ====== #
# Violin #
# ====== #
pdf("ALL_VlnPlot.pdf",width = 15, height = 15)
VlnPlot(GR, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.ribo"), pt.size = 0.5,ncol = 2)
dev.off()

pdf("ALL_VlnPlot3.pdf",width = 15, height = 8)
VlnPlot(GR, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), pt.size = 0.5,ncol = 3)
dev.off()    

combined = GR
#combined <- merge(x =s1, y = c(s2,s3,s4), add.cell.ids = c(sample1,sample2,sample3,sample4,sample5), merge.data = TRUE, project = project_name)
tail(combined@meta.data)
# Yan
#samples.list <- c("GC1001745_E11" = SO_GC1001745,"GC1001746_F11" = SO_GC1001746,"GC1004261_D4" = GC1004261_D4_FC,"GC1004262_E4" =GC10004262_E4_FC,"GC107955_A1"=GC107955_A1,"GC107975_B1"=GC107975_B1,"GC108010_G1"=GC108010_G1,"GC108011_H1"=GC108011_H1)
#combined <- NormalizeData(combined) %>% FindVariableFeatures() %>% ScaleData(split.by = "sample", do.center = FALSE) %>% RunPCA(verbose = FALSE)
#combined <- RunHarmony(combined, group.by.vars = "sample")
#combined <- RunUMAP(combined, reduction = "harmony", dims = 1:30)
#combined <- FindNeighbors(combined, reduction = "harmony", dims = 1:30) %>% FindClusters()

# library(future)
# plan("multiprocess", workers = 60) 
# options(future.globals.maxSize = 50000 * 1024^2) 
combined <- NormalizeData(combined, normalization.method = "LogNormalize", scale.f4actor = 10000)
combined <- FindVariableFeatures(combined, selection.method = "vst", nfeatures = 2000)
combined <- ScaleData(combined,slip.by="orig.ident", do.center=FALSE,  features = VariableFeatures(object = combined)) #, vars.to.regress = c("nCount_RNA", "percent.mt"))

#combined <- NormalizeData(combined, normalization.method = "LogNormalize", scale.factor = 10000) %>% FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>% ScaleData(split.by = "orig.ident", do.center = FALSE, features =  VariableFeatures(object = combined) , vars.to.regress = c("nCount_RNA","percent.mt"))
combined =RunPCA(combined,features = VariableFeatures(object = combined), npcs = 35, ndims.print = 1:5, nfeatures.print = 30)
saveRDS(combined, "before_runharmony.rds")

# New with more RAM
#ombined = readRDS("combined.rds")
#lan("multiprocess", workers = 10) 
#ptions(future.globals.maxSize = 50000 *0 * 1024^2) 
#
# nofuture
#ombined = readRDS("combined.rds")
#etwd(paste0(PROJECT,"/210421_Harmony_5PIN_nofuture"))
#ombined <- RunHarmony(combined, c("nCount_RNA","percent.mt","orig.ident"))
#aveRDS(combined, "combinednofuture.rds")
#

# regress on orig only - no future - works!
# combined = readRDS("combined.rds")
# setwd(paste0(PROJECT,"/210421_Harmony_5PIN_regressonorigonly"))
combined <- RunHarmony(combined, c("orig.ident"))


# combined = readRDS("combineorigonly.rds")
# combined@meta.data[,10:23] = NULL
# head(combined@meta.data)

combined <- RunUMAP(combined, reduction = "harmony", dims = 1:35)
combined <- FindNeighbors(combined, reduction = "harmony", dims = 1:35, k.param = 30, compute.SNN = TRUE, nn.method = "rann", nn.eps = 0) 
combined <- FindClusters(combined, resolution = reslist)
saveRDS(combined, "after_runharmony.rds")
head(combined@meta.data)

# ============== #
# Doublet Finder #
# ============== #

# DOUBLET FINDER 3 
# define the function

doubletFinder_v3 <- function(seu, PCs, pN = 0.25, pK, nExp, reuse.pANN = FALSE, sct = FALSE) {
  require(Seurat); require(fields); require(KernSmooth)
  ## Generate new list of doublet classificatons from existing pANN vector to save time
  if (reuse.pANN != FALSE ) {
    pANN.old <- seu@meta.data[ , reuse.pANN]
    classifications <- rep("Singlet", length(pANN.old))
    classifications[order(pANN.old, decreasing=TRUE)[1:nExp]] <- "Doublet"
    seu@meta.data[, paste("DF.classifications",pN,pK,nExp,sep="_")] <- classifications
    return(seu)
  }
  
  if (reuse.pANN == FALSE) {
    ## Make merged real-artifical data
    real.cells <- rownames(seu@meta.data)
    data <- seu@assays$RNA@counts[, real.cells]
    n_real.cells <- length(real.cells)
    n_doublets <- round(n_real.cells/(1 - pN) - n_real.cells)
    print(paste("Creating",n_doublets,"artificial doublets...",sep=" "))
    real.cells1 <- sample(real.cells, n_doublets, replace = TRUE)
    real.cells2 <- sample(real.cells, n_doublets, replace = TRUE)
    doublets <- (data[, real.cells1] + data[, real.cells2])/2
    colnames(doublets) <- paste("X", 1:n_doublets, sep = "")
    data_wdoublets <- cbind(data, doublets)
    
    ## Store important pre-processing information
    orig.commands <- seu@commands
    
    ## Pre-process Seurat object
    if (sct == FALSE) {
      print("Creating Seurat object...")
      seu_wdoublets <- CreateSeuratObject(counts = data_wdoublets)
      
      print("Normalizing Seurat object...")
      seu_wdoublets <- NormalizeData(seu_wdoublets,
                                     normalization.method = orig.commands$NormalizeData.RNA@params$normalization.method,
                                     scale.factor = orig.commands$NormalizeData.RNA@params$scale.factor,
                                     margin = orig.commands$NormalizeData.RNA@params$margin)
      
      print("Finding variable genes...")
      seu_wdoublets <- FindVariableFeatures(seu_wdoublets,
                                            selection.method = orig.commands$FindVariableFeatures.RNA$selection.method,
                                            loess.span = orig.commands$FindVariableFeatures.RNA$loess.span,
                                            clip.max = orig.commands$FindVariableFeatures.RNA$clip.max,
                                            mean.function = orig.commands$FindVariableFeatures.RNA$mean.function,
                                            dispersion.function = orig.commands$FindVariableFeatures.RNA$dispersion.function,
                                            num.bin = orig.commands$FindVariableFeatures.RNA$num.bin,
                                            binning.method = orig.commands$FindVariableFeatures.RNA$binning.method,
                                            nfeatures = orig.commands$FindVariableFeatures.RNA$nfeatures,
                                            mean.cutoff = orig.commands$FindVariableFeatures.RNA$mean.cutoff,
                                            dispersion.cutoff = orig.commands$FindVariableFeatures.RNA$dispersion.cutoff)
      
      print("Scaling data...")
      seu_wdoublets <- ScaleData(seu_wdoublets,
                                 features = orig.commands$ScaleData.RNA$features,
                                 model.use = orig.commands$ScaleData.RNA$model.use,
                                 do.scale = orig.commands$ScaleData.RNA$do.scale,
                                 do.center = orig.commands$ScaleData.RNA$do.center,
                                 scale.max = orig.commands$ScaleData.RNA$scale.max,
                                 block.size = orig.commands$ScaleData.RNA$block.size,
                                 min.cells.to.block = orig.commands$ScaleData.RNA$min.cells.to.block)
      
      print("Running PCA...")
      seu_wdoublets <- RunPCA(seu_wdoublets,
                              features = orig.commands$ScaleData.RNA$features,
                              npcs = length(PCs),
                              rev.pca =  orig.commands$RunPCA.RNA$rev.pca,
                              weight.by.var = orig.commands$RunPCA.RNA$weight.by.var,
                              verbose=FALSE)
      pca.coord <- seu_wdoublets@reductions$pca@cell.embeddings[ , PCs]
      cell.names <- rownames(seu_wdoublets@meta.data)
      nCells <- length(cell.names)
      rm(seu_wdoublets); gc() # Free up memory
    }
    
    if (sct == TRUE) {
      require(sctransform)
      print("Creating Seurat object...")
      seu_wdoublets <- CreateSeuratObject(counts = data_wdoublets)
      
      print("Running SCTransform...")
      seu_wdoublets <- SCTransform(seu_wdoublets)
      
      print("Running PCA...")
      seu_wdoublets <- RunPCA(seu_wdoublets, npcs = length(PCs))
      pca.coord <- seu_wdoublets@reductions$pca@cell.embeddings[ , PCs]
      cell.names <- rownames(seu_wdoublets@meta.data)
      nCells <- length(cell.names)
      rm(seu_wdoublets); gc()
    }
    
    ## Compute PC distance matrix
    print("Calculating PC distance matrix...")
    dist.mat <- fields::rdist(pca.coord)
    
    ## Compute pANN
    print("Computing pANN...")
    pANN <- as.data.frame(matrix(0L, nrow = n_real.cells, ncol = 1))
    rownames(pANN) <- real.cells
    colnames(pANN) <- "pANN"
    k <- round(nCells * pK)
    for (i in 1:n_real.cells) {
      neighbors <- order(dist.mat[, i])
      neighbors <- neighbors[2:(k + 1)]
      neighbor.names <- rownames(dist.mat)[neighbors]
      pANN$pANN[i] <- length(which(neighbors > n_real.cells))/k
    }
    
    print("Classifying doublets..")
    classifications <- rep("Singlet",n_real.cells)
    classifications[order(pANN$pANN[1:n_real.cells], decreasing=TRUE)[1:nExp]] <- "Doublet"
    seu@meta.data[, paste("pANN",pN,pK,nExp,sep="_")] <- pANN[rownames(seu@meta.data), 1]
    seu@meta.data[, paste("DF.classifications",pN,pK,nExp,sep="_")] <- classifications
    return(seu)
  }
}

Object = FC
Object = SetIdent(object = Object, value = "orig.ident")
table(Object@active.ident)

#Idents(Object) <- "Sample"
#Object@meta.data$orig.ident<- "Sample"
head(Object@meta.data)
Object$pANN <- "NA"
Object$pANNPredictions <- "NA"
head(Object@meta.data)


for(sample in unique(Object$orig.ident)){
  sample.cluster <- subset(Object, idents = sample)
  print(paste0("sample:", sample))
  length(rownames(sample.cluster@meta.data))
  expected.doublets <- ceiling(0.039 * length(rownames(sample.cluster@meta.data)))
  # sample.cluster <- doubletFinder_v3(sample.cluster, PCs = 1:20, nExp = expected.doublets, pN = 0.25, pK = 0.01)
  sample.cluster  <- doubletFinder_v3(sample.cluster, PCs = 1:20, pN = 0.25, pK = 0.01, nExp = expected.doublets, reuse.pANN = FALSE, sct=TRUE)
  sample.cluster$pANN <- sample.cluster@meta.data[colnames(sample.cluster), paste("pANN_0.25_0.01", expected.doublets, sep = "_")]
  sample.cluster$pANNPredictions <- sample.cluster@meta.data[colnames(sample.cluster), paste("DF.classifications_0.25_0.01", expected.doublets, sep = "_")]
  Object$pANN[colnames(sample.cluster)] <- sample.cluster$pANN[colnames(sample.cluster)]
  Object$pANNPredictions[colnames(sample.cluster)] <- sample.cluster$pANNPredictions[colnames(sample.cluster)]
  sample.cluster <- NULL
}

head(Object@meta.data)

FC@meta.data = Object@meta.data
head(FC@meta.data)

x = FC@meta.data
write.table(x, "DF_metadata.txt",col.names=NA, sep="\t")

pdf("doublets.pdf")
DimPlot(Object,pt.size = 0.1,label=F, label.size = 0,reduction = "umap",group.by = "pANNPredictions" )+theme(aspect.ratio = 1)
dev.off()

Idents(FC)="RNA_snn_res.3"

table(FC@meta.data$RNA_snn_res.3)
# try this way:
Idents(FC)="pANNPredictions"
table(FC@active.ident)

#Singlet Doublet 
# 357027   14502 

#subset all the doublets
FC # 371,529 cells
FC2 = subset(FC, idents="Singlet")
FC2 #357,027 Cells

# umap of only these
pdf("singlets_per_PANN.pdf")
DimPlot(FC2,pt.size = 0.1,label=F, label.size = 0,reduction = "umap",group.by = "pANNPredictions" )+theme(aspect.ratio = 1)
dev.off()

pdf("singlets_per_cluster.pdf")
DimPlot(FC2,pt.size = 0.1,label=F, label.size = 0,reduction = "umap",group.by = "RNA_snn_res.1.4" )+theme(aspect.ratio = 1)
dev.off()

saveRDS(FC2, "FC_removed_doublet_cells.rds")




# p1 <- DimPlot(combined, reduction = "umap", group.by = "orig.ident", pt.size=0.1, label = F, label.size=3.5)
# # p2 <- DimPlot(combined, reduction = "umap", group.by = "Status",pt.size=0.1, label = F, label.size=3.5,  repel = TRUE) 
# # p4 <- DimPlot(combined, reduction = "umap", group.by = "Subtype",pt.size=0.1, label = F, label.size=3.5,  repel = TRUE) 
# 
# pdf("umap_combined.pdf", width=15, height=8)
# # plot_grid(p1, p2,  p4)
# p1
# dev.off()
# 
# 
# for (metadata in c("Status","Subtype","orig.ident")){
#   
#   combined <- SetIdent(object = combined, value = metadata)
#   pdf(paste0("Split_",metadata,".pdf"), width=30)
#   print(DimPlot(combined, reduction = "umap", split.by = metadata))
#   dev.off()
#   
#   pdf(paste0("Split_",metadata,"_NoLegend.pdf"), width=30)
#   print(DimPlot(combined, reduction = "umap", split.by = metadata) + NoLegend())
#   dev.off()
# }

p1 <- DimPlot(combined, reduction = "umap", group.by = "orig.ident", pt.size=0.1, label = F, label.size=3.5)                                                                                                                                                           
p2 <- DimPlot(combined, reduction = "umap", group.by = "Status",pt.size=0.1, label = F, label.size=3.5,  repel = TRUE)                                                                                                                                                 
p3 = DimPlot(combined, reduction = "umap", group.by = "PIN",pt.size=0.1, label = F, label.size=3.5,  repel = TRUE) 
p4 <- DimPlot(combined, reduction = "umap", group.by = "Subtype",pt.size=0.1, label = F, label.size=3.5,  repel = TRUE)                                                                                                                                                
                                                                                                                                                                                                                                                                       
pdf("umap_combined.pdf", width=15, height=15)                                                                                                                                                                                                                           
plot_grid(p1, p2,  p3,p4)                                                                                                                                                                                                                                                 
dev.off()                                                                                                                                                                                                                                                              
                                                                                                                                                                                                                                                                       
                                                                                                                                                                                                                                                                       
for (metadata in c("Status","Subtype","orig.ident","PIN")){                                                                                                                                                                                                                  
                                                                                                                                                                                                                                                                       
combined <- SetIdent(object = combined, value = metadata)                                                                                                                                                                                                              
pdf(paste0("Split_",metadata,".pdf"), width=40)                                                                                                                                                                                                                        
print(DimPlot(combined, reduction = "umap", split.by = metadata))                                                                                                                                                                                                      
dev.off()                                                                                                                                                                                                                                                              
                                                                                                                                                                                                                                                                       
pdf(paste0("Split_",metadata,"_NoLegend.pdf"), width=40)                                                                                                                                                                                                               
print(DimPlot(combined, reduction = "umap", split.by = metadata) + NoLegend())                                                                                                                                                                                         
dev.off()                                                                                                                                                                                                                                                              
}                                                                                                                                                                                                                                                                      
                                                                                                                                                                                                                                                                       
                                                                                                                                                                                                                                                                       
# ================================================ #                                                                                                                                                                                                                   
# Run non-linear dimensional reduction (UMAP/tSNE) #                                                                                                                                                                                                                   
# ================================================ #                                                                                                                                                                                                                   
                                                                                                                                                                                                                                                                       
# CCA only                                                                                                                                                                                                                                                             
NPC=35                                                                                                                                                                                                                                                                 
                                                                                                                                                                                                                                                                       
FC = combined                                                                                                                                                                                                                                                          
Idents(FC)                                                                                                                                                                                                                                                             
                                                                                                                                                                                                                                                                       
for (metadata in c("Status","Subtype","orig.ident","PIN")){                                                                                                                                                                                                                  
                                                                                                                                                                                                                                                                       
                                                                                                                                                                                                                                                                       
                                                                                                                                                                                                                                                                       
Seurat_per_sample_Ident <- SetIdent(object = FC, value = metadata)                                                                                                                                                                                                     
Seurat_per_sample <- RunUMAP(Seurat_per_sample_Ident, dims = 1:NPC, reduction="harmony")                                                                                                                                                                               
Seurat_per_sample                                                                                                                                                                                                                                                      
                                                                                                                                                                                                                                                                       
pdf(paste0("Plot_per_",metadata,".pdf"))                                                                                                                                                                                                                               
                                                                                                                                                                                                                                                                       
print(UMAPPlot(Seurat_per_sample, reduction = "umap", pt.size=0.1,label=FALSE) +labs(title = paste0("UMAP per ",metadata))+ theme(plot.title = element_text(hjust = 0.5)))                                                                                             
dev.off()                                                                                                                                                                                                                                                              
                                                                                                                                                                                                                                                                       
pdf(paste0("Plot_per_",metadata,"_NoLegend.pdf"))                                                                                                                                                                                                                      
print(UMAPPlot(Seurat_per_sample, reduction = "umap", pt.size=0.1,label=FALSE) +labs(title = paste0("UMAP per ",metadata))+ theme(plot.title = element_text(hjust = 0.5)) +NoLegend())                                                                                 
dev.off()                                                                                                                                                                                                                                                              
}                                                                                                                                                                                                                                                                      


head(FC@meta.data)
# rename new = old
#FC@meta.data = FC@meta.data %>% dplyr::rename(RNA_snn_res.0.4 = RNA_snn_res.0.4,RNA_snn_res.0.5 = RNA_snn_res.0.5,RNA_snn_res.0.6 = RNA_snn_res.0.6,RNA_snn_res.0.7 = RNA_snn_res.0.7, RNA_snn_res.0.8 = RNA_snn_res.0.8,RNA_snn_res.0.9 = RNA_snn_res.0.9, RNA_snn_re
s.1= RNA_snn_res.1 ,RNA_snn_res.1.1= RNA_snn_res.1.1,RNA_snn_res.1.2= RNA_snn_res.1.2, RNA_snn_res.1.4= RNA_snn_res.1.4, RNA_snn_res.1.6 = RNA_snn_res.1.6)

head(FC@meta.data)
saveRDS(FC, file = "combined_FC.rds")

pdf("Combined_Clusters_diff_res.pdf")                                                                                                                                                                                                                         [97/1964]
for (resolution in c("0.4","0.5","0.6","0.7","0.8","0.9","1","1.1","1.2","1.3","1.4","1.5","1.6","3")){
#pdf(paste0("res_",resolution,"Integrated_Clusters_diff_res.pdf")) 
Seurat <- SetIdent(object = FC, value = paste0('RNA_snn_res.',resolution))
Seurat <- RunUMAP(Seurat, dims = 1:NPC, reduction = "harmony")
print(UMAPPlot(Seurat, reduction = "umap", pt.size=0.1,label=T) +labs(title = paste0("UMAP with res ",resolution))+ theme(plot.title = element_text(hjust = 0.5)))
#dev.off()
}
dev.off()


head(FC@meta.data)
u=FC@meta.data
write.table(u,"Combined_metadata.txt",col.names=NA, sep="\t")


for (resolution in c("0.4","0.5","0.6","0.7","0.8","0.9","1","1.1","1.2","1.3","1.4","1.5","1.6","3")){
Seurat_reso <- SetIdent(object = FC, value = paste0('RNA_snn_res.',resolution))
cell_perC=table(Seurat_reso@active.ident)
cell_perC
write.table(cell_perC, paste0("Combined_res_",resolution,"_cell_per_c.txt"),col.names=NA, sep="\t")

cell_perC_perO=table(Seurat_reso@active.ident, Seurat_reso@meta.data$orig.ident)
cell_perC_perO
write.table(cell_perC_perO, paste0("Combined_res_",resolution,"_cell_per_c_per_Origin.txt"),col.names=NA, sep="\t")

cell_perC_perA=table(Seurat_reso@active.ident, Seurat_reso@meta.data$Subtype)
cell_perC_perA
write.table(cell_perC_perA, paste0("Combined_res_",resolution,"_cell_per_c_per_Subtype.txt"),col.names=NA, sep="\t")

cell_perC_perS=table(Seurat_reso@active.ident, Seurat_reso@meta.data$Status)
cell_perC_perS
write.table(cell_perC_perS, paste0("Combined_res_",resolution,"_cell_per_c_per_Status.txt"),col.names=NA, sep="\t")

cell_perC_perS=table(Seurat_reso@active.ident, Seurat_reso@meta.data$PIN)
cell_perC_perS
write.table(cell_perC_perS, paste0("Combined_res_",resolution,"_cell_per_c_per_PIN.txt"),col.names=NA, sep="\t")
# 
# cell_perC_perS=table(Seurat_reso@active.ident, Seurat_reso@meta.data$Status)
# cell_perC_perS
# write.table(cell_perC_perS, paste0("Combined_res_",resolution,"_cell_per_c_per_Status.txt"),col.names=NA, sep="\t")
}

for (resolution in c("0.4","0.5","0.6","0.7","0.8","0.9","1","1.1","1.2","1.3","1.4","1.5","1.6","3")){
Seurat_reso <- SetIdent(object = FC, value = paste0('RNA_snn_res.',resolution))
library(future)
plan("multiprocess", workers = 60)   
options(future.globals.maxSize = 10000 * 1024^2) 
Seurat.markers <- FindAllMarkers(Seurat_reso,test.use="MAST",  min.pct = 0.1, logfc.threshold = 0.25, max.cells.per.ident = Inf, latent.vars = NULL,min.cells.feature = 0, min.cells.group = 0)
write.table(Seurat.markers, paste0("Combined_res.",resolution,"_MAST_markers_unannotated_Clusters.txt"),col.names=NA, sep="\t")
tsne.markers=Seurat.markers         
#DefaultAssay(Seurat_reso) <- 'integrated'                     
topn = tsne.markers %>% group_by(cluster) %>% top_n(10, avg_log2FC)
topn$gene                                                     
pdf(paste0("Combined_res.",resolution,"_heatmap_top10genes.pdf"),width=25, height=20)                                              
   
print(DoHeatmap(object = Seurat_reso, features=topn$gene, group.by = "ident", label = TRUE, angle=45) + ggtitle("top 10 genes expressed per cluster") + NoLegend())
dev.off()
topn = tsne.markers %>% group_by(cluster) %>% top_n(15, avg_log2FC)
topn$gene
pdf(paste0("Combined_res.",resolution,"_heatmap_top15genes.pdf"),width=25, height=25)
print(DoHeatmap(object = Seurat_reso, features=topn$gene, group.by = "ident", label = TRUE, angle=45) + ggtitle("top 15 genes expressed per cluster") + NoLegend())
dev.off()
topn = tsne.markers %>% group_by(cluster) %>% top_n(20, avg_log2FC)
topn$gene                                                                                                                             
pdf(paste0("Combined_res.",resolution,"_heatmap_top20genes.pdf"),width=25, height=30)
print(DoHeatmap(object = Seurat_reso, features=topn$gene, group.by = "ident", label = TRUE, angle=45) + ggtitle("top 20 genes expressed per cluster") + NoLegend())
dev.off()
topn = tsne.markers %>% group_by(cluster) %>% top_n(5, avg_log2FC) 
topn$gene
pdf(paste0("Combined_res.",resolution,"_heatmap_top5genes.pdf"),width=25, height=15)
print(DoHeatmap(object = Seurat_reso, features=topn$gene, group.by = "ident", label = TRUE, angle=45) + ggtitle("top 5 genes expressed per cluster") + NoLegend())
dev.off()
topn = tsne.markers %>% group_by(cluster) %>% top_n(20, avg_log2FC) 
topn$gene
topn
write_csv(topn, paste0("Combined_res.",resolution,"_top20genespercluster_logFC.csv"))
}

q()
n
