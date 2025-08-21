# R file

library("Seurat")
library("tidyverse")
library("SCpubr")
library("viridis")
library("RColorBrewer")
library("devtools")
library("harmony")
library("ggplot2")
library("SCpubr")

packageVersion("Seurat")

setwd("/path/Figure3_CD")
Seu=readRDS("combined_FC.rds")
head(Seu@meta.data,2)

Idents(Seu)="RNA_snn_res.3" # go to the chosen resolution chosen after inspection of gene lists
levels(Seu)

library(future)
plan("multicore", workers = 15) 
options(future.globals.maxSize = 400000 * 1024^2) # 90Go


current.cluster.ids <- c("0","1","2",...)
new.cluster.ids <- c("TA","AbC-progenitor","SC","AQP8+ CEACAM7+ AbC",...)
names(new.cluster.ids) <- levels(Seu)
Seu <- RenameIdents(Seu, new.cluster.ids)
Seu <- StashIdent(Seu, save.name = "Cluster")
head(Seu@meta.data)
table(Seu@active.ident)

Seurat_Final=Seu
head(Seurat_Final@meta.data,2)

pdf("Figure_CD_SCpubR.pdf",width=15,height=13)
SCpubr::do_DimPlot(sample = Seurat_Final, label = TRUE, label.color = "white",label.fill=NULL, repel=TRUE,label.size=7, legend.icon.size=4,font.size=20)
dev.off()

# =============== #
# Final Gene list # 
# =============== #

# make sure you do not use Seurat 5 for that function, we do not want to use the presto package, we had weird results
DefaultAssay(Seurat_Final) <- "RNA"
Idents(Seurat_Final)="Cluster"
table(Seurat_Final@active.ident)
ORGANISM="Human"
plan("multicore", workers = 10) 
options(future.globals.maxSize = 10000 * 1024^2) 

if(ORGANISM=="Mouse"){                                       
    genes.use <- grep(pattern = "^Rp[sl][[:digit:]]|^Rplp[[:digit:]]|^Rpsa",rownames(Seurat_Final), value=TRUE, invert=TRUE) #get list of non-ribosomal genes 
} else { 
    genes.use <- grep(pattern = "^RP[SL][[:digit:]]|^RPLP[[:digit:]]|^RPSA",rownames(Seurat_Final), value=TRUE, invert=TRUE) #get list of non-ribosomal genes 
}
 
 
Seurat.markers <- FindAllMarkers(Seurat_Final,test.use="wilcox",  min.pct = 0.1, logfc.threshold = 0.25, features = genes.use, only.pos = FALSE)
write.table(Seurat.markers, "wilcox_markers_annotated_Clusters.txt",col.names=NA, sep="\t")
tsne.markers=Seurat.markers
topn = tsne.markers %>% group_by(cluster) %>% top_n(10, avg_log2FC)
topn$gene
pdf("heatmap_top10genes.pdf",width=20, height=25)
DoHeatmap(object = Seurat_Final, features=topn$gene, group.by = "ident", label = TRUE, angle=90) + ggtitle("top 10 genes expressed per cluster") + NoLegend() + theme(text = element_text(size = 20)) ##+ scale_fill_gradient(low = "deepskyblue2", high = "darkred")  
dev.off()

# ======== #
# Dotplots #
# ======== #

genes <- list(
"Junctional complex"
=c("OCLN","TJP1","TJP2","TJP3","CLDN1","CLDN2","CLDN3","CLDN4","CLDN5","CLDN6","CLDN7","CLDN8","CLDN9","CLDN10","CLDN11","CLDN12","CLDN13","CLDN14","CLDN15","CLDN16","CLDN17","CLDN18","CLDN19","CLDN20","CLDN21","CLDN22","CLDN23","CLDN24","CLDN25"),
"Anti-microbial peptides"
=c("DEFA5","DEFA6","DEFB1","CAMP","REG3G"),
"Mucins"
=c("MUC1","MUC2","MUC3A","MUC4","MUC5AC","MUC5B","MUC12","MUC13","MUC15","MUC17","MUC20"))

GROUP="Gut_Liver_axis"
genes <- list(
"Bile acid transporter"
=c("SLC51B","SLC51A","SLC10A2","SLC27A2","SLC5A9","ABCG5"),
"Bile acid receptors"
=c("NR1H4","GPBAR1","VDR","NR1I2","NR1I3","AHR","S1PR2","RORC","PPARA","PPARD","PPARG"))

GROUP="Nutrient_transporters"
genes <- list(
"Lipid"
=c("APOA1","APOC3","APOB","FABP6","FABP2","PLIN2","PLIN3","SAR1B","ABCG2","SLC27A4","APOM","APOC2","MGAT3","LPL","ACAT1","ACSL3","HMGCR","ACAT2","FFAR4"),
"Vitamin"
=c("RBP2","TCN2","CYP4F2","SLC23A1","SLC52A1","RDH5","CYP4F3","SLC23A3","SLC22A4","VNN1","BTD","BCO1","BCO2","CD320","DHR59","RBP4"),
"Amino acid"
=c("SLC7A9","SLC6A19","SLC3A1","SLC7A7","SLC3A2","SLC43A2","SLC1A1","SLC6A20","SLC6A6","SLC1A7","SLC25A39","SLC38A2","SLC38A1","SLC1A5","SLC25A13","SLC7A6","SLC38A5","SLC25A12","SLC1A4","SLC43A1"),
"Sugar"
=c("SLC5A1","SLC37A4","SLC2A5","SLC5A9","SLC2A2","SLC5A11","SLC50A1","SLC2A10","SLC2A1"),
"Nucleotide"
=c("SLC28A1","SLC25A36","SLC35A3","NTSE","SLC35A1","SLC35B3","SLC35B2","SLC25A33","SLC35D1","SLC35A2","SLC17A5","SLC28A2","SLC29A1","SLC29A2","SLC17A9"),
"Water"
=c("AQP1","AQP3","AQP7","AQP8","AQP11"),
"Metal ion"
=c("SLC39A4","SLC25A37","SLC39A5","SLC31A1","SLC31A2","SLC30A4","SLC39A7","SLC39A8","SLC30A9","SLC39A3","SLC39A9","SLC30A7","ATP2A3","SLC9A2","SCNN1A","SCNN1B","KCNS3","SLC30A5","SLC30A6","SLC25A28","KCTD14","KCTD17")
)

GROUP="inorganic_organic_solute_carriers"
genes <- list(
"Inorganic solute"
=c("SLC25A3","SLC20A1","NHERF1","SLC4A7","SLC13A1","SLC34A3","SLC34A2","SLC26A2","SLC4A4"),
"Organic solute"
=c("SLC13A2","SLC16A4","SLC16A5","SLC44A1","SLC44A3","SLC16A1","SLC5A3","SLC2A13","SLC22A23","SLC16A3","SLC26A6","SLC16A9","SLCO3A1")
)

GROUP="Ion_transporters_exchangers"
genes <- list(
"Ion exchanger"
=c("SLC26A1","SLC26A2","SLC26A3","SLC26A4","SLC26A5","SLC26A6","SLC26A7","SLC26A8","SLC26A9","SLC26A10","SLC26A11"),
"Na+/H+ exchangers"
=c("SLC9A1","SLC9A3","SLC9A4","SLC9A5","SLC9A6","SLC9A7","SLC9A8","SLC9A9","SLC9A10","SLC9A11"),
"CFTR pathway"
=c("CFTR","NHERF1","NHERF2","NHERF4","SH3","SHANK2","GOPC","ADRB2","LPAR2","PDZK1","MRP4"),
"Guanylate cyclase"
=c("GUCA2A","GUCA2B","GUCY2C"),
"Na+/K+ ATPase"
=c("ATP1A1","ATP1A2","ATP1A3","ATP1A4","ATP1B1","ATP1B2","ATP1B3","ATP1B4"),
"K+ channels"
=c("KCNE1","KCNE2","KCNE3","KCNE4","KCNE5","KCNN1","KCNN2","KCNN3","KCNN4","KCNQ1","KCNQ2","KCNQ3","KCNQ4","KCNQ5"),
"Na+K+Cl- T."
=c("SLC12A1","SLC12A2"),
"Epithelial Na+ channel"
=c("SCNNG")
)

GROUP="Genes_associated_with_diarrhea"
genes <- list(
"Congenital/infectious diarrhea"
=c("CFTR","SLC26A3","SLC9A3","GUCY2C"),
"Pri. bile acid diarr."
=c("SLC51B"),
"Gluc./galac. malabs."
=c("SLC5A1"))

Idents(Seurat_Final)="Subtype"
p1<-SCpubr::do_DotPlot(sample = Seurat_Final, features = genes, use_viridis = FALSE,na.value = "#dd1c77",dot_border = FALSE, font.size=10,legend.type="colorbar")
p3=p2+scale_colour_gradientn(colors = c("#2171b5","white","#a31818"))

#pdf(paste0("FIGURE3_",GROUP,"_Dotplot_per_subtype.pdf"), height=5, width=12)
#pdf(paste0("FIGURE3_",GROUP,"_Dotplot_per_subtype.pdf"), height=5, width=8)
#pdf(paste0("FIGURE3_",GROUP,"_Dotplot_per_subtype.pdf"), height=5, width=20)
#pdf(paste0("FIGURE3_",GROUP,"_Dotplot_per_subtype.pdf"), height=5, width=8)
#pdf(paste0("FIGURE3_",GROUP,"_Dotplot_per_subtype.pdf"), height=5, width=21)
pdf(paste0("FIGURE3_",GROUP,"_Dotplot_per_subtype.pdf"), height=5, width=10)
p3
dev.off()

Idents(Seurat_Final)="Cluster"
p1<-SCpubr::do_DotPlot(sample = Seurat_Final, features = genes, use_viridis = FALSE,na.value = "#dd1c77",dot_border = FALSE, font.size=10,legend.type="colorbar")    
p2=p1+theme(strip.background = element_rect(color="black", linetype="solid"))
p3=p2+scale_colour_gradientn(colors = c("#2171b5","white","#a31818"))

#pdf(paste0("FIGURE3_",GROUP,"_Dotplot_per_cluster.pdf"), height=8, width=15)
#pdf(paste0("FIGURE3_",GROUP,"_Dotplot_per_cluster.pdf"), height=8, width=8)
#pdf(paste0("FIGURE3_",GROUP,"_Dotplot_per_cluster.pdf"), height=8, width=20)
#pdf(paste0("FIGURE3_",GROUP,"_Dotplot_per_cluster.pdf"), height=8, width=8)
#pdf(paste0("FIGURE3_",GROUP,"_Dotplot_per_cluster.pdf"), height=8, width=21)
pdf(paste0("FIGURE3_",GROUP,"_Dotplot_per_cluster.pdf"), height=8, width=10)
p3
dev.off()

p1 <- DotPlot(Seurat_Final, features = genes, group.by="Cluster",split.by="Subtype", dot.scale=8, cols="RdBu") + theme(text = element_text(size = 20), axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
p2=p1+theme(strip.background = element_rect(color="black", linetype="solid"), axis.title= element_blank()) 
p3=p2+ scale_colour_gradientn(colors = c("#2171b5","white","#a31818")) 

#pdf(paste0("FIGURE3_",GROUP,"_Dotplot_per_cluster_per_subtype.pdf"), height=20, width=30)
#pdf(paste0("FIGURE3_",GROUP,"_Dotplot_per_cluster_per_subtype.pdf"), height=20, width=15)
#pdf(paste0("FIGURE3_",GROUP,"_Dotplot_per_cluster_per_subtype.pdf"), height=25, width=30)
#pdf(paste0("FIGURE3_",GROUP,"_Dotplot_per_cluster_per_subtype.pdf"), height=20, width=15)
#pdf(paste0("FIGURE3_",GROUP,"_Dotplot_per_cluster_per_subtype.pdf"), height=20, width=47)
pdf(paste0("FIGURE3_",GROUP,"_Dotplot_per_cluster_per_subtype.pdf"), height=20, width=18)
p3
dev.off()

# === #
# DGE #
# === #

Seurat_Final=SetIdent(Seurat_Final, value="orig.ident")
head(Seurat_Final@meta.data)
table(Seurat_Final@active.ident)
GRP1=grep(pattern = "GC1001745_E11|GC1001746_F11|GC107955_A1|GC107975_B1|GC108010_G1|GC108011_H1|GC108673_SI-GA-B5|GC117610_SI-GA-C1|GC117611_SI-GA-D1|GC117612_SI-GA-G2|GC117613_SI-GA-H2", colnames(x = Seurat_Final), value=TRUE) # HV 
GRP2=grep(pattern = "GC111886_SI-GA-C5|GC111887_SI-GA-D5|GC112519_SI-GA-H7|GC112520_SI-GA-A8|GC117608_SI-GA-A1|GC117609_SI-GA-B1|GC128024_SI-GA-D8|GC128025_SI-GA-E8|GC130099_SI-GA-C12|GC110188_SI-GA-C8|GC110189_SI-GA-D8|GC110308_SI-GA-B9|GC110526_SI-GA-D10|GC110527_SI-GA-E10|GC111001_SI-GA-E6|GC111002_SI-GA-F6|GC112070_SI-GA-C6|GC112071_SI-GA-D6|GC114506_SI-GA-C2|GC114507_SI-GA-D2", colnames(x = Seurat_Final), value=TRUE) # CD
GRP3=grep(pattern = "GC110188_SI-GA-C8|GC110189_SI-GA-D8|GC110308_SI-GA-B9|GC110526_SI-GA-D10|GC110527_SI-GA-E10|GC111001_SI-GA-E6|GC111002_SI-GA-F6|GC112070_SI-GA-C6|GC112071_SI-GA-D6|GC114506_SI-GA-C2|GC114507_SI-GA-D2", colnames(x = Seurat_Final), value=TRUE) # CD+IBS 
GRP4=grep(pattern = "GC111886_SI-GA-C5|GC111887_SI-GA-D5|GC112519_SI-GA-H7|GC112520_SI-GA-A8|GC117608_SI-GA-A1|GC117609_SI-GA-B1|GC128024_SI-GA-D8|GC128025_SI-GA-E8|GC130099_SI-GA-C12", colnames(x = Seurat_Final), value=TRUE) # CD-IBS  
    
length(GRP1) 
length(GRP2)
length(GRP3)
length(GRP4)

Seurat_Final=SetIdent(Seurat_Final, value="Cluster")

for (CLUSTER in c("AbC-immature","GC-immature","CA1+ CEACAM7+ AbC-immature",...)){

    plan("multicore", workers = 15) 
    options(future.globals.maxSize = 10000 * 1024^2) 
    GRP1=grep(pattern = "GC1001745_E11|GC1001746_F11|GC107955_A1|GC107975_B1|GC108010_G1|GC108011_H1|GC108673_SI-GA-B5|GC117610_SI-GA-C1|GC117611_SI-GA-D1|GC117612_SI-GA-G2|GC117613_SI-GA-H2", colnames(x = Seurat_Final), value=TRUE) # HV 
    GRP2=grep(pattern = "GC111886_SI-GA-C5|GC111887_SI-GA-D5|GC112519_SI-GA-H7|GC112520_SI-GA-A8|GC117608_SI-GA-A1|GC117609_SI-GA-B1|GC128024_SI-GA-D8|GC128025_SI-GA-E8|GC130099_SI-GA-C12|GC110188_SI-GA-C8|GC110189_SI-GA-D8|GC110308_SI-GA-B9|GC110526_SI-GA-D10|GC110527_SI-GA-E10|GC111001_SI-GA-E6|GC111002_SI-GA-F6|GC112070_SI-GA-C6|GC112071_SI-GA-D6|GC114506_SI-GA-C2|GC114507_SI-GA-D2", colnames(x = Seurat_Final), value=TRUE) # CD
    GRP3=grep(pattern = "GC110188_SI-GA-C8|GC110189_SI-GA-D8|GC110308_SI-GA-B9|GC110526_SI-GA-D10|GC110527_SI-GA-E10|GC111001_SI-GA-E6|GC111002_SI-GA-F6|GC112070_SI-GA-C6|GC112071_SI-GA-D6|GC114506_SI-GA-C2|GC114507_SI-GA-D2", colnames(x = Seurat_Final), value=TRUE) # CD+IBS 
    GRP4=grep(pattern = "GC111886_SI-GA-C5|GC111887_SI-GA-D5|GC112519_SI-GA-H7|GC112520_SI-GA-A8|GC117608_SI-GA-A1|GC117609_SI-GA-B1|GC128024_SI-GA-D8|GC128025_SI-GA-E8|GC130099_SI-GA-C12", colnames(x = Seurat_Final), value=TRUE) # CD-IBS  
    length(GRP1) 
    length(GRP2) 
    length(GRP3) 
    length(GRP4) 
    
    CLUSTER
    GP1 = WhichCells(Seurat_Final, ident=CLUSTER, cells = GRP1)
    length(GP1)
    GP2 = WhichCells(Seurat_Final, ident=CLUSTER, cells = GRP2)
    length(GP2)
    GP3 = WhichCells(Seurat_Final, ident=CLUSTER, cells = GRP3)
    length(GP3)
    GP4 = WhichCells(Seurat_Final, ident=CLUSTER, cells = GRP4)
    length(GP4)

  if(ORGANISM=="Mouse"){
    genes.use <- grep(pattern = "^Rp[sl][[:digit:]]|^Rplp[[:digit:]]|^Rpsa|^Igkc|^Jchain|^Iglc2|^Igha1|^Igha2|^mt-",rownames(sub), value=TRUE, invert=TRUE) #get list of non-ribosomal genes
} else {
    genes.use <- grep(pattern = "^RP[SL][[:digit:]]|^RPLP[[:digit:]]|^RPSA|^IGKC|^JCHAIN|^IGLC2|^IGHA1|^IGHA2|^MT-",rownames(sub), value=TRUE, invert=TRUE) #get list of non-ribosomal genes 
}
    
    DE=FindMarkers(Seurat_Final, test.use = "wilcox", ident.1=GP1, ident.2=GP2, only.pos = F, logfc.threshold = 0, min.pct = 0, features=genes.use)
    write.table(DE, paste0("DE_HV_vs_allCD_in_",CLUSTER,".txt"),col.names=NA, sep="\t")
    # HV vs IBS-C
    DE=FindMarkers(Seurat_Final, test.use = "wilcox", ident.1=GP1, ident.2=GP3, only.pos = F, logfc.threshold = 0, min.pct = 0, features=genes.use)
    write.table(DE, paste0("DE_HV_vs_CD+IBS_in_",CLUSTER,".txt"),col.names=NA, sep="\t") 
    # HV vs IBS-D
    DE=FindMarkers(Seurat_Final, test.use = "wilcox", ident.1=GP1, ident.2=GP4, only.pos = F, logfc.threshold = 0, min.pct = 0, features=genes.use)
    write.table(DE, paste0("DE_HV_vs_CD-IBS_in_",CLUSTER,".txt"),col.names=NA, sep="\t")  
    #IBS-C vs IBS-D
    DE=FindMarkers(Seurat_Final, test.use = "wilcox", ident.1=GP3, ident.2=GP4, only.pos = F, logfc.threshold = 0, min.pct = 0, features=genes.use)
    write.table(DE, paste0("DE_CD+IBS_vs_CD-IBS_in_",CLUSTER,".txt"),col.names=NA, sep="\t")  
} 

# ============= #
# Volcano plots #
# ============= #

for (CLUSTER in c("AbC-immature","GC-immature","CA1+ CEACAM7+ AbC-immature",...)){
    for (THRESH in c(0.5,1,1.5,2,3)){
    
G1="HV"
G2="all_CD"
G3="CD+IBS"
G4="CD-IBS"

de=read.table(paste0("DE_HV_vs_allCD_in_",CLUSTER,".txt"), sep="\t", header=T)
#de=read.table(paste0("DE_HV_vs_CD+IBS_in_",CLUSTER,".txt"), sep="\t", header=T)
#de=read.table(paste0("DE_HV_vs_CD-IBS_in_",CLUSTER,".txt"), sep="\t", header=T)
#de=read.table(paste0("DE_CD+IBS_vs_CD-IBS_in_",CLUSTER,".txt"), sep="\t", header=T)

print(CLUSTER)
print(THRESH)
group1=G1
group2=G2
## add a column of NAs
de$diffexpressed <- "NO"
de$diffexpressed[de$avg_log2FC > THRESH & de$p_val_adj < 0.05] <- "UP"
de$diffexpressed[de$avg_log2FC < -THRESH & de$p_val_adj < 0.05] <- "DOWN"
de$delabel <- NA
colnames(de)[1]="name"
de$delabel[de$diffexpressed != "NO"] <- as.character(de$name[de$diffexpressed != "NO"])
library(ggrepel)
g =ggplot(data=de, aes(x=avg_log2FC, y=-log10(p_val_adj), col=diffexpressed, label=delabel)) + 
    geom_point(size=2) + 
    theme_minimal() +
    geom_text_repel(show_guide  = F , size=4) +
	ggtitle(paste0(group1," vs ",group2," in cluster ",CLUSTER)) +
	scale_colour_manual("Genes", labels=c(paste0("upregulated in ",group2),paste0("not sign. and below threshold of logFC ±", THRESH),paste0("upregulated in ",group1)),
    values=c("#2471A3", "gray80", "#F39C12"))+
    guides(colour = guide_legend(override.aes = list(size=4))) +
    geom_vline(xintercept=c(-THRESH,THRESH), linetype="dotted") +
	theme(plot.title = element_text(color="black", size=12, face="bold.italic"),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
panel.background = element_blank(), axis.line = element_line(colour = "black"))
ggsave(paste0("Volcano_",group1,"_vs_",group2,"_in_",CLUSTER,"_threshold_",THRESH,".pdf"),g, width=8, height=5)	
}}

q()
n
