# R file

# it always makes group 2 versus group 1 (logFC + in group 2 - blue)
# group 1 and group 2 are by alphabetical order
# eg: HV and IBS --> here HV will be group 1 and IBS will be group 2
# so HV will be group 1 - will have logFC negative - red
# IBS will be group 2 - will have positive logFC - in blue
# think about it if you want your HV to be blue --> then, they need to be group 2

library(miloR)
library(SingleCellExperiment)
library(dplyr)
library(patchwork)
library(Seurat)
library(ggplot2)
library(scater)

PROJECT="25xxxx_milor_HV_Rectum_vs_Colon"
mainDir<-"/path/Figure_1_HV"
subDir<-paste0("Result_",PROJECT)
dir.create(file.path(mainDir, subDir), showWarnings = FALSE)
opath<-paste0(mainDir,"/",subDir,"/")
rdspath<-paste0(mainDir,"/")
setwd(opath)

# if object is of seurat5 - BEFORE LAUNCHING THE SCRIPT
#        conda activate XXX # With Seurat 5
#        R
#        library("Seurat")
#        library("scCustomize")
#        PROJECT="/path/"
#        setwd(PROJECT)
#        FC<-readRDS("FINAL.rds")
#        obj <- Convert_Assay(seurat_object = FC, convert_to = "V3", assay="RNA")
#        saveRDS(obj, “FINAL_4.rds")


# Really needs not to be a Seurat 5 object!
FC<-readRDS(file=paste0(rdspath,'Figure_1_HV.rds')) 

condition="Colon_vs_rectum"
print(condition)

Idents(FC)<-"Tissue"
table(FC$Tissue)
FC$sample<-FC$orig.ident

# if you need to rename so the first group is red (logFC negative)
#FC@meta.data$Tissue[FC@meta.data$Tissue == "HV"] <- "2_HV"
#FC@meta.data$Tissue[FC@meta.data$Tissue == "IBS-D"] <- "1_IBS-D"
#FC@meta.data$Tissue[FC@meta.data$Tissue == "IBS-C"] <- "1_IBS-C"

Idents(FC)<-"Tissue"
table(FC@meta.data$orig.ident, FC@meta.data$Tissue)

# it's only one vs one, you will need to remove a group if you have more than 2
# sub_SO<-FC[,!(FC@active.ident %in% c("1_IBS-D"))] # the one you remove
# if you don't remove any
sub_SO<-FC

table(sub_SO$Tissue)
table(sub_SO$orig.ident)

sub_SO_se <- as.SingleCellExperiment(sub_SO)

# ## 2 visualize the data:
sub_SO_se$sample<-sub_SO_se$orig.ident

## 3 DA testing
### 3.1 create a milo object
sub_SO_milo <- Milo(sub_SO_se)
### 3.2 Construct KNN graph
sub_SO_milo <- buildGraph(sub_SO_milo, k = 35, d = 35, reduced.dim = "PCA")
### 3.3 Defining representative neighbourhoods on the KNN graph
sub_SO_milo <- makeNhoods(sub_SO_milo, prop = 0.1, k = 35, d=35, refined = TRUE, reduced_dims = "PCA")
### 3.4 Counting cells in neighbourhoods
sub_SO_milo <- countCells(sub_SO_milo, meta.data = as.data.frame(colData(sub_SO_milo)), sample="sample")
head(nhoodCounts(sub_SO_milo))
### 3.5 Defining experimental design
sub_SO_milo_design <- data.frame(colData(sub_SO_milo))[,c("sample", "Tissue")]
## Convert batch info from integer to factor
sub_SO_milo_design$Tissue <- as.factor(sub_SO_milo_design$Tissue) 
dim(sub_SO_milo_design)
sub_SO_milo_design <- distinct(sub_SO_milo_design)
dim(sub_SO_milo_design)
rownames(sub_SO_milo_design) <- sub_SO_milo_design$sample
sub_SO_milo_design # you need multiple samples (1st column) per group (2nd column)
## 3.6 Computing neighbourhood connectivity
sub_SO_milo <- calcNhoodDistance(sub_SO_milo, d=35, reduced.dim = "PCA")
sub_SO_milo 
        # class: Milo 
        # dim: 14858 21391 
        # metadata(0):
        # assays(2): counts logcounts
        # rownames(14858): Xkr4 Gm1992 ... Vamp7 Spry3
        # rowData names(0):
        # colnames(21391): GC141670_SI-GA-B12_AAACCCAAGCATGTTC-1
        #   GC141670_SI-GA-B12_AAACCCAAGGGATCGT-1 ...
        #   GC141671_SI-GA-C12_TTTGTTGTCCATTGTT-1
        #   GC141671_SI-GA-C12_TTTGTTGTCTAGTACG-1
        # colData names(26): orig.ident nCount_RNA ... sample ident
        # reducedDimNames(3): PCA INTEGRATED.HARMONY UMAP
        # mainExpName: RNA3
        # altExpNames(0):
        # nhoods dimensions(2): 21391 1357
        # nhoodCounts dimensions(2): 1357 4
        # nhoodDistances dimension(1): 1357
        # graph names(1): graph
        # nhoodIndex names(1): 1357
        # nhoodExpression dimension(2): 1 1
        # nhoodReducedDim names(0):
        # nhoodGraph names(0):
        # nhoodAdjacency dimension(2): 1 1
da_results <- testNhoods(sub_SO_milo, design = ~ Tissue, design.df = sub_SO_milo_design)
head(da_results)
da_results %>%
  arrange(SpatialFDR) %>%
  head() 

## 4 Inspecting DA testing results
pdf(paste0(condition,"_pvalue_histogram.pdf"))
ggplot(da_results, aes(PValue)) + geom_histogram(bins=50)
dev.off()

pdf(paste0(condition,"_pvalue_hline.pdf"))
ggplot(da_results, aes(logFC, -log10(SpatialFDR))) + 
  geom_point() +
  geom_hline(yintercept = 1)
dev.off()

sub_SO_milo <- buildNhoodGraph(sub_SO_milo)
umap_pl<-DimPlot(sub_SO, split.by='Tissue',group.by='Tissue',cols=c("#a31818", "#2171b5"),ncol=1)
nh_graph_pl <- plotNhoodGraphDA(sub_SO_milo, da_results, layout="UMAP",alpha=0.1, pt.size=1) 

pdf(paste0(condition,"_umap_nl_graph.pdf"),width=6, height=10)
umap_pl + nh_graph_pl + plot_layout(guides="collect")
dev.off()

# color groups
# cellular heterogeneity
da_results <- groupNhoods(sub_SO_milo, da_results, max.lfc.delta = 10, overlap = 1) #da.fdr=0.1 default - if no DA found, increase a bit !!! if change there, change below too with the DAbeeswarn
head(da_results)

pdf(paste0(condition,"_plotNhoodGroups.pdf"),height= 8, width=9)
plotNhoodGroups(sub_SO_milo, da_results, show_groups = NULL, alpha=0.05, delta=3) # delta changes the number of group --> delta 1 or delt 3 for example, delta 3 gives less groups
dev.off()

da_results <- annotateNhoods(sub_SO_milo, da_results, coldata_col = "Cluster")
head(da_results)

pdf(paste0(condition,"_Cluster_fraction.pdf"))
ggplot(da_results, aes(Cluster_fraction)) + geom_histogram(bins=50)
dev.off()

da_results$anno_celltype_cluster <- ifelse(da_results$Cluster_fraction < 0.7, "Mixed", da_results$Cluster)

pdf(paste0(condition,"_plotDAbeeswarm.pdf"),height= 15, width=17)
plotDAbeeswarm(da_results, group.by = "Cluster") + theme(text = element_text(size = 40)) # alpha=same value as da.fdr above
dev.off()


write.table(da_results, paste0(condition,"_da_results.txt"),col.names=NA, sep="\t")

library(dplyr)
MILO_results=da_results
head(MILO_results)
z <- MILO_results %>% group_by(anno_celltype_cluster) %>% summarise(across(everything(), mean))
z
write.table(z, paste0("Milo_",PROJECT,"_",condition,"_mean.txt"),col.names=NA, sep="\t")

