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

L.data<-Read10X("./data/MF309_resequenced/P76_followup") 

#Create Seurat Object
#MF312agg<-CreateSeuratObject(MF312agg.data) #inspect: no metadata about origin (thin/thick) present, dont go on with this
Les<-CreateSeuratObject(L.data) 

#clean workspace
rm(list = ls(pattern = ".data"))

######################## VDJ specific ###########################################################
#load VDJ data from csv file supplied
L_clone<-read.csv("./data/MF309_resequenced/P76_followup/TCR/filtered_contig_annotations.csv")

#check gpr -1 at the end of the barcode in sample
head(colnames(Les)) #no -1

# Remove the -1 at the end of each barcode.
L_clone$barcode <- gsub("-1", "", L_clone$barcode)

# Subset so only the first line of each barcode is kept,
# as each entry for given barcode will have same clonotype.
L_clone <- L_clone[!duplicated(L_clone$barcode), ]


# Only keep the barcode and clonotype columns.
# We'll get additional clonotype info from the clonotype table.
L_clone <- L_clone[,c("barcode", "raw_clonotype_id")]
names(L_clone)[names(L_clone) == "raw_clonotype_id"] <- "clonotype_id"



# Clonotype-centric info.
clono_L <- read.csv("./data/MF309_resequenced/P76_followup/TCR/clonotypes.csv")

# Slap the AA sequences onto our original table by clonotype_id.
L_clone<- merge(L_clone, clono_L[, c("clonotype_id", "cdr3s_aa")])


# Reorder so barcodes are first column and set them as rownames.
L_clone <- L_clone[, c(2,1,3)]
rownames(L_clone) <- L_clone[,1]
L_clone[,1] <- NULL

# Add to the Seurat object's metadata.
Les <- AddMetaData(object=Les, metadata=L_clone)

rm(clono_thick)
rm(clono_thin)
rm(Les_clone)
rm(NonLes_clone)

tail(Les$cdr3s_aa, 20)

sum(grepl(pattern ="TRA:CAVGTSGSRLTF;TRB:CASSLALSGGAYNEQFF",Les$cdr3s_aa)) #5360

##############################End of VDJ specific #########################################


# Add additional metadata
Les <- AddMetaData(Les, "followup",  col.name= "tissue")

Les <- AddMetaData(Les, "MF309",  col.name= "patient")
NonLes <- AddMetaData(NonLes, "MF309", col.name= "patient")

Les <- AddMetaData(Les, "MF309_followup",  col.name= "sample")

#### quality control and filter ####

Les[["percent.mt"]] <- PercentageFeatureSet(Les, pattern = "^MT-")
head(Les@meta.data, 15)
NonLes[["percent.mt"]] <- PercentageFeatureSet(NonLes, pattern = "^MT-")
head(NonLes@meta.data, 5)



# QC plots
pdf("./output/QC/VlnPlot_MF309_followup.pdf",
    width=14, height=7)
VlnPlot(Les, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.off()



Les[["percent.ribo"]] <- PercentageFeatureSet(Les, pattern = "^RP[SL][[:digit:]]")


VlnPlot(Les, features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.ribo"), ncol = 4, pt.size = 0)
VlnPlot(NonLes, features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.ribo"), ncol = 4, pt.size = 0)

####  ------------------------------------------ code to get to sce1 object
Les.sce = as.SingleCellExperiment(Les)
NonLes.scr = as.SingleCellExperiment(NonLes)


## ----------------------------------------------code for doublet estimation (from Matthias Wielscher)
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
Les.sce2=Les.sce[,doub_score<thresh]  #### exluding the doublets from the data set




## ----------------------------------------------------- finished
rm(list=ls(pattern="score"))
rm(list=ls(pattern="exclude_"))
rm(list=ls(pattern="thresh"))

####  convert back to seurat object 
#important: use further for making a list for integration:
Les2 = as.Seurat(Les.sce2) #17476 (was 9734) cells
NonLes2 = as.Seurat(NonLes.sce2) #7057 (was 7106) cells

####until here 14.28

#compute 2 SD up and down from mean of total mRNAs (nCount_RNA)
upper_bound <- 10^(mean(log10((Les)$nCount_RNA)) +  #19705
                     2*sd(log10((Les)$nCount_RNA)))
lower_bound <- 10^(mean(log10((Les)$nCount_RNA)) -
                     2*sd(log10((Les)$nCount_RNA))) #396

#plot data
qplot(nCount_RNA, data = Les@meta.data, color = "tissue", geom =
        "density") +
  geom_vline(xintercept = lower_bound) +
  geom_vline(xintercept = upper_bound)

#nFeature
#compute 2 SD up and down from mean of total mRNAs (nFeature_RNA)
upper_bound <- 10^(mean(log10((Les2)$nFeature_RNA)) +  #4713
                     2*sd(log10((Les2)$nFeature_RNA)))
lower_bound <- 10^(mean(log10((Les2)$nFeature_RNA)) -
                     2*sd(log10((Les2)$nFeature_RNA))) #223
#plot data
qplot(nCount_RNA, data = Les@meta.data, color = "tissue", geom =
        "density") +
  geom_vline(xintercept = lower_bound) +
  geom_vline(xintercept = upper_bound)


VlnPlot(subset(Les2, nFeature_RNA<500), features = "nFeature_RNA", pt.size = 0) #250
VlnPlot(subset(Les2, nFeature_RNA>2200), features = "nFeature_RNA", pt.size = 0) #4700

VlnPlot(subset(Les2, percent.mt<40), features = "percent.mt", pt.size = 0) #12
VlnPlot(subset(NonLes2, percent.mt<40), features = "percent.mt", pt.size = 0) #12

plot1 <- FeatureScatter(Les2, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(Les2, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot3 <- FeatureScatter(Les2, feature1 = "nCount_RNA", feature2 = "percent.ribo")
plot1 + plot2 + plot3
FeatureScatter(Les2, feature1 = "nFeature_RNA", feature2 = "percent.mt")



#filter each sample seperately
Les3 <- subset(Les2, subset = nFeature_RNA > 250 & nFeature_RNA < 4700 & percent.mt<12)

##how many genes per cell after filtering
median(Les$nFeature_RNA)  #1147
median(Les3$nFeature_RNA)  #1403


#Save
saveRDS(Les3, "./data/MF309_followup_filtered.RDS")


FU<-readRDS( "./data/MF309_followup_filtered.RDS")


NonLes3<-readRDS("./data/MF309_NonLesional_filtered.RDS")




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
saveRDS(MF3.list, "./data/Pat65_list_sctransform.RDS")

####until here 1.3.2021

#### Instead Merge:
comb <- merge(Les3, y = NonLes3, add.cell.ids = c("lesional", "nonlesional"), project = "Patient65")
comb <- NormalizeData(comb, normalization.method = "LogNormalize", scale.factor = 10000)
comb <- FindVariableFeatures(comb, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(comb)
comb <- ScaleData(comb, vars.to.regress = "percent.mt")#scales only variable features (if you want to scale all use features = all.genes)
comb<-RunPCA(comb, features = VariableFeatures(object = comb))
PCHeatmap(comb)
ElbowPlot(comb, ndims = 50 )
comb<-FindNeighbors(comb, dims=1:23)
comb<-FindClusters(comb, resolution = 0.6)
comb<-RunUMAP(comb, dims=1:23)
DimPlot(comb, reduction = "umap", label = TRUE)
DimPlot(comb, split.by = "tissue", label = TRUE, pt.size = 0.5)
DimPlot(comb, group.by = "tissue", label = TRUE, pt.size = 0.5)

