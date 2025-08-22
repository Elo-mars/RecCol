
# R file
library(Seurat)                                                                                                                         
library(devtools)
library(CytoTRACE2)
PROJECT="/path/Figure_HV"
setwd(PROJECT)
data=readRDS("Figure_HV.rds")
data

dir.create("CytoTrace")
setwd(paste0(PROJECT,"/CytoTraceB"))
getwd()

cytotrace2_result <- cytotrace2(data, species = "human",is_seurat = TRUE,slot_type="counts",full_model = TRUE,batch_size = 100000,smooth_batch_size = 10000,max_pcs = 200,seed = 14, ncores=60)
write.table(cytotrace2_result@meta.data, "cytotrace2_result.txt",col.names=FALSE, sep="\t")
head(cytotrace2_result)
saveRDS(cytotrace2_result, "Cyto.rds")
x=cytotrace2_result@meta.data
head(x)
write.table(x, "cytotrace2_result.txt",col.names=NA, sep="\t")

annotation <- data.frame(phenotype = data@meta.data$Cluster) %>% set_rownames(., colnames(data))
head(annotation)

plots <- plotData(cytotrace2_result,
                  annotation = annotation,
                  expression_data = NULL,
                  is_seurat = TRUE,
                  pc_dims = 35,
                  seed = 14)

pdf("Cytotrace_potency_score.pdf", width=8, height=8)
print(plots$CytoTRACE2_UMAP)
dev.off()

# Potency score distribution by phenotype
# A boxplot of predicted potency score separated by phenotype/group from the annotation file. Can be used to assess the distribution of predicted potency scores across different cell phenotypes.
pdf("Cytotrace_boxplot.pdf", width=8, height=8)
print(plots$CytoTRACE2_Boxplot_byPheno + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)))
dev.off()

#Potency category
#The UMAP embedding plot of predicted potency category reflects the discrete classification of cells into potency categories, taking possible values of Differentiated, Unipotent, Oligopotent, Multipotent, Pluripotent, and Totipotent.
pdf("Cytotrace_potency_category.pdf", width=8, height=8)
print(plots$CytoTRACE2_Potency_UMAP)
dev.off()

#Relative order
#UMAP embedding of predicted relative order, which is based on absolute predicted potency scores normalized to the range 0 (more differentiated) to 1 (less differentiated). Provides the relative ordering of cells by developmental potential
pdf("Cytotrace_relative_order.pdf", width=8, height=8)
print(plots$CytoTRACE2_Relative_UMAP)
dev.off()

#Phenotypes
#UMAP colored by phenotype annotation. Used to assess the distribution of cell phenotypes across the UMAP space.
pdf("Cytotrace_by_phenotype.pdf", width=8, height=8)
print(plots$Phenotype_UMAP)
dev.off()

q()
n
