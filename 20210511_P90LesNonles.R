##################### QC and  preliminary analysis ####################
####################### lesional vs non lesional #####################

set.seed(2)

setwd("~/Derma/CTCL/LES_vs_NONLES_vs_HC")
#setwd("../Derma/LES_vs_NONLES_vs_HC") #on windows

###load packages ###
#devtools::install_github(repo = "satijalab/seurat", ref = "develop")
library(Seurat)
library(ggplot2)
library(cowplot)
library(magrittr)
library(dplyr)

Pat<-90
L.data<-Read10X(paste("./data/P",Pat,"_resequenced/lesional", sep="")) 
NL.data<-Read10X(paste("./data/P",Pat,"_resequenced/nonlesional", sep=""))

#Create Seurat Object
#MF312agg<-CreateSeuratObject(MF312agg.data) #inspect: no metadata about origin (thin/thick) present, dont go on with this
Les<-CreateSeuratObject(L.data) 
NonLes<-CreateSeuratObject(NL.data) 

#clean workspace
rm(list = ls(pattern = ".data"))

######################## VDJ specific ###########################################################
#load VDJ data from csv file supplied
L_clone<-read.csv(paste("./data/P",Pat,"_resequenced/lesional/TCR/filtered_contig_annotations.csv", sep=""))
NL_clone<-read.csv(paste("./data/P",Pat,"_resequenced/nonlesional/TCR/filtered_contig_annotations.csv", sep=""))

#check gpr -1 at the end of the barcode in sample
head(colnames(Les)) #no -1
head(colnames(NonLes)) #no -1


# Remove the -1 at the end of each barcode.
L_clone$barcode <- gsub("-1", "", L_clone$barcode)
NL_clone$barcode <- gsub("-1", "", NL_clone$barcode)

# Subset so only the first line of each barcode is kept,
# as each entry for given barcode will have same clonotype.
L_clone <- L_clone[!duplicated(L_clone$barcode), ]
NL_clone <- NL_clone[!duplicated(NL_clone$barcode), ]


# Only keep the barcode and clonotype columns.
# We'll get additional clonotype info from the clonotype table.
L_clone <- L_clone[,c("barcode", "raw_clonotype_id")]
names(L_clone)[names(L_clone) == "raw_clonotype_id"] <- "clonotype_id"

NL_clone <- NL_clone[,c("barcode", "raw_clonotype_id")]
names(NL_clone)[names(NL_clone) == "raw_clonotype_id"] <- "clonotype_id"

# Clonotype-centric info.
clono_L <- read.csv(paste("./data/P",Pat,"_resequenced/lesional/TCR/clonotypes.csv", sep=""))
clono_NL <- read.csv(paste("./data/P",Pat,"_resequenced/nonlesional/TCR/clonotypes.csv", sep=""))

# Slap the AA sequences onto our original table by clonotype_id.
L_clone<- merge(L_clone, clono_L[, c("clonotype_id", "cdr3s_aa")])
NL_clone<- merge(NL_clone, clono_NL[, c("clonotype_id", "cdr3s_aa")])


# Reorder so barcodes are first column and set them as rownames.
L_clone <- L_clone[, c(2,1,3)]
rownames(L_clone) <- L_clone[,1]
L_clone[,1] <- NULL

NL_clone <- NL_clone[, c(2,1,3)]
rownames(NL_clone) <- NL_clone[,1]
NL_clone[,1] <- NULL

# Add to the Seurat object's metadata.
Les <- AddMetaData(object=Les, metadata=L_clone)
NonLes <- AddMetaData(object=NonLes, metadata=NL_clone)

tail(Les$cdr3s_aa, 20)
tail(NonLes$cdr3s_aa, 200)
head(clono_L)
head(clono_NL)

sum(grepl(pattern ="TRA:CILRDLNNDMRF;TRB:CSARGWGADTGELFF",Les$cdr3s_aa)) #42
sum(grepl(pattern ="TRA:CILRDLNNDMRF;TRB:CSARGWGADTGELFF",NonLes$cdr3s_aa)) #29

rm(clono_L)
rm(clono_NL)
rm(L_clone)
rm(NL_clone)

##############################End of VDJ specific #########################################


# Add additional metadata
Les <- AddMetaData(Les, "Lesional",  col.name= "tissue")
NonLes <- AddMetaData(NonLes, "Nonlesional", col.name= "tissue")

Les <- AddMetaData(Les, paste0("Pat",Pat),  col.name= "patient")
NonLes <- AddMetaData(NonLes, paste0("Pat",Pat), col.name= "patient")

Les <- AddMetaData(Les, paste0("Pat",Pat, "_lesional"),  col.name= "sample")
NonLes <- AddMetaData(NonLes, paste0("Pat",Pat, "_nonlesional"), col.name= "sample")

#### quality control and filter ####

Les[["percent.mt"]] <- PercentageFeatureSet(Les, pattern = "^MT-")
head(Les@meta.data, 15)
NonLes[["percent.mt"]] <- PercentageFeatureSet(NonLes, pattern = "^MT-")
head(NonLes@meta.data, 5)



# QC plots
pdf(paste0("./output/QC/VlnPlot_Pat", Pat, "_Les.pdf"),
    width=14, height=7)
VlnPlot(Les, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.off()

pdf(paste0("./output/QC/VlnPlot_Pat", Pat, "_NonLes.pdf"),
    width=14, height=7)
VlnPlot(NonLes, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.off()

Les[["percent.ribo"]] <- PercentageFeatureSet(Les, pattern = "^RP[SL][[:digit:]]")
NonLes[["percent.ribo"]] <- PercentageFeatureSet(NonLes, pattern = "^RP[SL][[:digit:]]")


VlnPlot(Les, features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.ribo"), ncol = 4, pt.size = 0)
VlnPlot(NonLes, features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.ribo"), ncol = 4, pt.size = 0)

####  ------------------------------------------ code to get to sce1 object
Les.sce = as.SingleCellExperiment(Les)
NonLes.scr = as.SingleCellExperiment(NonLes)


## ----------------------------------------------code for doublet estimation (from Matthias Wielscher)
set.seed(2)
library(scran)
doublet_score=scran::doubletCells(Les.sce, k = 200,force.match=TRUE)  # k can be reduced to 50 -- works faster 
doub_score=log10(doublet_score+1)
thresh=median(doub_score)+3*mad(doub_score)



## plotting results
hist(doub_score,breaks=25,main="Distribution of Doublet Scores",xlab="log(doublet scores)",ylab=c("number of cells"))
abline(v=thresh,col="blue",lty=2,lwd=2)
legend("topright", legend=c("Doublet exclusion threshold"),
       col=c("blue"), lty=c(2),lwd=1.5, cex=0.8)
#Add doublet score to object
Les.sce$DoubletScore<-doub_score

exclude_doublet=table(doub_score<thresh)[1]
sum(exclude_doublet)# cells (doublets) to exclude =14

Les.sce2=Les.sce[,doub_score<thresh]  #### exluding the doublets from the data set

#redo for second sample
doublet_score1=scran::doubletCells(NonLes.scr, k = 200,force.match=TRUE)  # k can be reduced to 50 -- works faster 
doub_score1=log10(doublet_score1+1)
thresh1=median(doub_score1)+3*mad(doub_score1) #very low 1,8977

## plotting results
hist(doub_score1,breaks=25,main="Distribution of Doublet Scores",xlab="log(doublet scores)",ylab=c("number of cells"))
abline(v=thresh1,col="blue",lty=2,lwd=2)
legend("topright", legend=c("Doublet exclusion threshold"),
       col=c("blue"), lty=c(2),lwd=1.5, cex=0.8)
#Add doublet score to object
NonLes.scr$DoubletScore<-doub_score1
#exclude doublets
exclude_doublet1=table(doub_score1<thresh1)[1]
sum(exclude_doublet1)# cells (doublets) to exclude =30

NonLes.sce2=NonLes.scr[,doub_score1<thresh1] 


## ----------------------------------------------------- finished
rm(list=ls(pattern="score"))
rm(list=ls(pattern="exclude_"))
rm(list=ls(pattern="thresh"))

####  convert back to seurat object 
#important: use further for making a list for integration:
Les2 = as.Seurat(Les.sce2) #9147 (was 9161) cells
NonLes2 = as.Seurat(NonLes.sce2) #1494 (was 1524) cells



#compute 2 SD up and down from mean of total mRNAs (nCount_RNA)
upper_bound <- 10^(mean(log10((Les)$nCount_RNA)) +  #8662
                     2*sd(log10((Les)$nCount_RNA)))
lower_bound <- 10^(mean(log10((Les)$nCount_RNA)) -
                     2*sd(log10((Les)$nCount_RNA))) #189

#plot data
qplot(nCount_RNA, data = Les@meta.data, color = "tissue", geom =
        "density") +
  geom_vline(xintercept = lower_bound) +
  geom_vline(xintercept = upper_bound)

#nFeature
#compute 2 SD up and down from mean of total mRNAs (nFeature_RNA)
upper_bound <- 10^(mean(log10((Les2)$nFeature_RNA)) +  #2458
                     2*sd(log10((Les2)$nFeature_RNA)))
lower_bound <- 10^(mean(log10((Les2)$nFeature_RNA)) -
                     2*sd(log10((Les2)$nFeature_RNA))) #175

upper_bound <- 10^(mean(log10((NonLes2)$nFeature_RNA)) +  #5827
                     2*sd(log10((NonLes2)$nFeature_RNA)))
lower_bound <- 10^(mean(log10((NonLes2)$nFeature_RNA)) -
                     2*sd(log10((NonLes2)$nFeature_RNA))) #183

#plot data
qplot(nFeature_RNA, data = NonLes2@meta.data, color = "tissue", geom =
        "density") +
  geom_vline(xintercept = lower_bound) +
  geom_vline(xintercept = upper_bound)

VlnPlot(subset(Les2, nFeature_RNA<1200), features = "nFeature_RNA", pt.size = 0) #250
VlnPlot(subset(NonLes2, nFeature_RNA<1200), features = "nFeature_RNA") #250
VlnPlot(subset(Les2, percent.mt<40), features = "percent.mt",pt.size=0) #12 cave two pop
VlnPlot(subset(NonLes2, percent.mt<40), features = "percent.mt", pt.size = 0) #12

plot1 <- FeatureScatter(Les2, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(Les2, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") #4200 nfeature
plot3 <- FeatureScatter(Les2, feature1 = "nCount_RNA", feature2 = "percent.ribo")
plot1 + plot2 + plot3
FeatureScatter(Les2, feature1 = "nFeature_RNA", feature2 = "percent.mt")

plot1 <- FeatureScatter(NonLes2, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(NonLes2, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")#4000 nf
plot3 <- FeatureScatter(NonLes2, feature1 = "nCount_RNA", feature2 = "percent.ribo")
plot1 + plot2 + plot3
FeatureScatter(NonLes2, feature1 = "nFeature_RNA", feature2 = "percent.mt")




#filter each sample seperately
Les3 <- subset(Les2, subset = nFeature_RNA > 250 & nFeature_RNA < 4000 & percent.mt<12)
               
NonLes3 <- subset(NonLes2, subset = nFeature_RNA > 250 & nFeature_RNA < 4000 &  percent.mt<12)

##how many genes per cell after filtering
median(Les$nFeature_RNA)  #513
median(Les3$nFeature_RNA)  #514

median(NonLes$nFeature_RNA)  #1164
median(NonLes3$nFeature_RNA)  #1336

#Save
saveRDS(Les3, paste0("./data/Pat", Pat, "_Lesional_filtered.RDS"))
saveRDS(NonLes3, paste0("./data/Pat", Pat, "_Nonlesional_filtered.RDS"))

###### 

Les3<-readRDS( paste0("./data/Pat", Pat, "_Lesional_filtered.RDS"))
NonLes3<-readRDS(paste0("./data/Pat", Pat, "_Nonlesional_filtered.RDS"))




###SCT Transformation
#Make list for SCT normalization
MF3.list<-c(Les3, NonLes3)

#for 2.5 Gigabyte:
4300*1024^2# = 4508876800 
#future.globals.maxSize
options(future.globals.maxSize= 4508876800)

for (i in 1:length(MF3.list)) {
  MF3.list[[i]] <- SCTransform(MF3.list[[i]], verbose = FALSE,  vars.to.regress = "percent.mt", 
                               conserve.memory = F, return.only.var.genes = F)
}

#save SCT transformed objects
saveRDS(MF3.list, paste0("./data/Pat",Pat, "_list_sctransform.RDS"))

####until here 11.5.2021

#### Instead Merge:
comb <- merge(Les3, y = NonLes3, add.cell.ids = c("lesional", "nonlesional"), project = paste0("Patient",Pat))
comb <- NormalizeData(comb, normalization.method = "LogNormalize", scale.factor = 10000)
comb <- FindVariableFeatures(comb, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(comb)
comb <- ScaleData(comb, vars.to.regress = "percent.mt")#scales only variable features (if you want to scale all use features = all.genes)
comb<-RunPCA(comb, features = VariableFeatures(object = comb))
PCHeatmap(comb)
ElbowPlot(comb, ndims = 50 )
comb<-FindNeighbors(comb, dims=1:25)
comb<-FindClusters(comb, resolution = 0.6)
comb<-RunUMAP(comb, dims=1:25)
DimPlot(comb, reduction = "umap", label = TRUE)
DimPlot(comb, split.by = "tissue", label = TRUE, pt.size = 0.5)
DimPlot(comb, group.by = "tissue", label = TRUE, pt.size = 0.5)


#
#DefaultAssay(MF.combined)<-"integrated"
#MF.combined<-ScaleData(MF.combined)
#MF.combined<-RunPCA(MF.combined)
#ElbowPlot(MF.combined,ndims = 50 ) #22
#MF.combined<-FindNeighbors(MF.combined, dims=1:22)
#MF.combined<-FindClusters(MF.combined, resolution = 0.6)
#MF.combined<-RunUMAP(MF.combined, dims=1:22)

#### UMAPPlots, unnamed clusters ####
DimPlot(comb, label=T, reduction = "umap")
DimPlot(MF.combined, label=T, reduction = "umap", group.by = "tissue")


Idents(comb)<-"seurat_clusters"
pdf("./output/Pat65_UMAPPlot.pdf",
    width=7, height=5)
DimPlot(comb, label=T, reduction = "umap")
dev.off()

pdf("./output/Pat65_UMAPPlot_splitbytissue.pdf",
    width=10, height=5)
DimPlot(comb, label=T, reduction = "umap")
dev.off()



##### Save combined object
saveRDS(comb, "./data/Pat65.combined.RDS" )
comb<-readRDS("./data/Pat65.combined.RDS" )
#### Clustermarker- plots ####
#DefaultAssay(MF.combined)<-"RNA"

#Cluster markers acc to NatMed Publ: Nat Med. 2019 Aug; 25(8): 1251–1259. 31359002
Clustermarker<-c("CD3D", "CD8A", "GZMA", "CD4", "FOXP3", "KLRC1", "KLRC3", "CD19", "CD79A",
                 "SLAMF7", "IGKC", "FCGR2A", "CSF1R", "FLT3", "COL1A2", "MCAM", "MYLK", "PECAM1", "VWF",
                 "FAP", "MKI67", "KRT5", "TPSAB1", "MRC1", "PMEL")
#CD3G, CD3D, CD3E, CD2 (T cells), 
#CD8A, GZMA (CD8+ T cells), 
#CD4, FOXP3 (CD4+ T cells/Tregs),
#KLRC1, KLRC3 (NK cells),
#CD19, CD79A (B cells), #cluster12
#SLAMF7, IGKC (Plasma cells), #cluster 22 
#FCGR2A, CSF1R (Macrophages), 
#FLT3 (Dendritic cells),
#CLEC4C (Plasmacytoid Dendritic cells), 
#COL1A2 (Fibroblasts), 
#MCAM, MYLK (Myofibroblasts),
#FAP, PDPN (CAFs), #
#PECAM1, VWF (Endothelial cells), 
#TPSAB1 (Mast cells)
#MRC1 -M2 macrophage
pdf(paste0("./output/Pat" ,Pat, "_FeaturePlot_clustermarker.pdf"),
    width=21, height=21)
FeaturePlot(comb, features =Clustermarker , cols = c("lightgrey", "red"), ncol = 5, order=T)
dev.off()

FeaturePlot(comb, features= "DoubletScore")
VlnPlot(comb, features= "DoubletScore", split.by = "tissue")
########### subset for clonotypes ##########
clonotype_seq<-table(comb$cdr3s_aa,comb$tissue )
write.csv(table(comb$cdr3s_aa, comb$tissue), "./output/Pat65_cdr3s_aa_frequency.csv")

#add metadata about TRA/TRB chains
comb$tcr <- ifelse(grepl(pattern ="TRA:CILRDLNNDMRF;TRB:CSARGWGADTGELFF",comb$cdr3s_aa), "mal_clono", 
                          ifelse(!is.na(comb$cdr3s_aa), "other_clono", "n/a"))
head(comb$tcr ) 
sum(grepl(pattern ="CILRDLNNDMRF",comb$cdr3s_aa)) #71 cells
sum(grepl(pattern ="CSARGWGADTGELFF",comb$cdr3s_aa)) #109 cells

Idents(comb) <- "tcr"
pdf(paste0("./output/Pat", Pat, "_UMAPPlot_clono.pdf"),
    width=7, height=7)
DimPlot(comb, label=T, reduction = "umap")
dev.off()

pdf(paste0("./output/Pat", Pat, "_UMAPPlot_clono_splitbytissue.pdf"),
    width=10, height=5)
DimPlot(comb, label=T, reduction = "umap", split.by = "tissue")
dev.off()




##### Save combined object
saveRDS(comb, paste0("./data/Pat", Pat,"_combined.RDS"))
