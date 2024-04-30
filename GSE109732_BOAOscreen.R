library(patchwork)
library(Seurat)
library(dplyr)
library(future)
library(ggplot2)
library(ggpubr)
library(RColorBrewer)
options(warn=-1)
set.seed(1)
library(Hmisc)

GCB1.data<-Read10X("D:/GCB1_2.2.filtered_feature_bc_matrix/")
GCB2.data<-Read10X("D:/GCB2_2.2.filtered_feature_bc_matrix/")

GCB1.data<-as.data.frame(GCB1.data)
GCB2.data<-as.data.frame(GCB2.data)

# GCB1.metadata<-read.table("d:/GCB1metadata.txt",sep = "\t",header = TRUE,row.names = 1)
# GCB2.metadata<-read.table("d:/GCB2metadata.txt",sep = "\t",header = TRUE,row.names = 1)
# 
# GCB1.data<-GCB1.data[,rownames(GCB1.metadata)]
# GCB2.data<-GCB2.data[,rownames(GCB2.metadata)]

for (i in 1:length(colnames(GCB1.data))) {
  colnames(GCB1.data)[i] <- paste(colnames(GCB1.data)[i],"GCB1",i,sep = "-")  
}

for (i in 1:length(colnames(GCB2.data))) {
  colnames(GCB2.data)[i] <- paste(colnames(GCB2.data)[i],"GCB2",i,sep = "-")  
}

# rownames(GCB1.metadata)<-colnames(GCB1.data)
# rownames(GCB2.metadata)<-colnames(GCB2.data)
# 
# GCB1.metadata$group<-rep("GCB1",length(colnames(GCB1.data)))
# GCB2.metadata$group<-rep("GCB2",length(colnames(GCB2.data)))

# pbmc.metadata<-rbind(GCB1.metadata,GCB2.metadata)

pbmc.data<-cbind(GCB1.data,GCB2.data)
pbmc.data <- as.data.frame(pbmc.data)
pbmc.data_t<-as.data.frame(t(pbmc.data))
pbmc.data_t<-filter(pbmc.data_t,Ighd <=0)
pbmc.data<-as.data.frame(t(pbmc.data_t))
pbmc <- CreateSeuratObject(counts = pbmc.data, project = "GCB",min.cells = 3, min.features = 200)
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^mt-")
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
library(tidyverse)
library(SummarizedExperiment)
library(scuttle)
ref.data<-read.table("D:/CJJ immunity/GSE109732/GSE109732_GeneExpression_DZ_Fraction1-3_GC-PB.txt",sep = "\t",header = TRUE)
ref.data<-ref.data[!duplicated(ref.data$Gene),]
rownames(ref.data)<-ref.data$Gene
ref.data<-ref.data[,-1]
ref.metadata<-data.frame(colnames(ref.data))
ref.metadata$Clustername<-ref.metadata[,1]
colnames(ref.metadata)[1]<-"ID"
rownames(ref.metadata)<-ref.metadata[,1]
ds_seurat<-SummarizedExperiment(assays=list(counts=ref.data),colData = ref.metadata[,c("ID","Clustername")])
ds_seurat<-logNormCounts(ds_seurat)
library(SingleR)
remove(GCB1.data,GCB2.data,pbmc.data,pbmc.data_t,ref.data,ref.metadata,i)
gc()
pbmc_raw<-pbmc
remove(pbmc)
gc()
library(ggridges)
library(ggplot2)
setwd("e:/CJJscreen_BOAO/")
for (s in c(2500,3000,3500,4000,4500,5000,5500,6000,6500,7000)) {
  for (t in c(5,10,15,20,25,30,35,40,45)) {
    pbmc <- subset(pbmc_raw, subset = nFeature_RNA > 200 & nFeature_RNA < s & percent.mt < t)
    pbmc <- NormalizeData(pbmc)
    pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000, verbose=F)
    all.genes <- rownames(pbmc)
    pbmc <- ScaleData(pbmc, features = all.genes, verbose=F)
    pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc), verbose=F,npcs = 100)
    pred.hesc <- SingleR(test = pbmc@assays$RNA@counts, ref = ds_seurat, assay.type.test=1,
                         labels = ds_seurat$Clustername)
    
    Idents(pbmc)<-pred.hesc@listData$labels
    pbmc@meta.data$Cluster<-Idents(pbmc)
    
    
    Idents(pbmc)<-pbmc$Cluster
    a<-levels(Idents(pbmc))
    new.cluster.ids <-a
    for (q in 1:length(a)) {
      if (a[q] == "GC.PB_1"){
        new.cluster.ids[q] <- "prePB"
      }
      if (a[q] == "GC.PB_2"){
        new.cluster.ids[q] <- "prePB"
      }
      if (a[q] == "GC.PB_3"){
        new.cluster.ids[q] <- "prePB"
      }
      if (a[q] == "Fraction1_1"){
        new.cluster.ids[q] <- "prePB"
      }
      if (a[q] == "Fraction1_2"){
        new.cluster.ids[q] <- "prePB"
      }
      if (a[q] == "Fraction1_3"){
        new.cluster.ids[q] <- "prePB"
      }
    }
    names(new.cluster.ids) <- levels(pbmc)
    pbmc <- RenameIdents(pbmc, new.cluster.ids)
    pbmc@meta.data$newcluster<-Idents(pbmc)
    pbmc<-RunLDA(pbmc,labels = pbmc$newcluster,reduction.name = "lda3")
    pbmc<-RunUMAP(pbmc,reduction = "lda3",reduction.name = "lda3umap",dims = 1:9)
    pbmc<-RunTSNE(pbmc,reduction = "lda3",reduction.name = "lda3tsne",dims = 1:9)
    Idents(pbmc)<-pbmc$newcluster
    Glutamine_1<-read.table("d:/GSE60927.txt",sep = "\t")
    Glutamine_1_new<-Glutamine_1[1:length(rownames(Glutamine_1)),]
    Glutamine_1_new<-intersect(Glutamine_1_new,rownames(pbmc))
    Glutamine1<-FetchData(pbmc,vars = Glutamine_1_new)
    for (i in 1:length(rownames(Glutamine1))) {
      Glutamine1$Average[i]<-mean(as.numeric(Glutamine1[i,]))
    }
    
    Total<-FetchData(pbmc,vars = rownames(pbmc))
    Data<-data.frame()
    remove(i)
    for (j in 1:length(colnames(Total))) {
      Data[1,j]<-mean(Total[,j])
    }
    average_Glutamine1<-mean(Glutamine1$Average)
    colnames(Data)<-colnames(Total)
    Data_T<-as.data.frame(t(Data))
    remove(j)
    rownamesDataT<-as.data.frame(rownames(Data_T))
    Data_T<-cbind(Data_T,rownamesDataT)
    colnames(Data_T)<-c("S","Gene")
    Sortdata<-Data_T[order(Data_T$S,decreasing = TRUE),]
    
    Sortdatahi<-filter(Sortdata,S > average_Glutamine1)
    Sortdatalo<-filter(Sortdata,S < average_Glutamine1)
    higene<-rownames(Sortdatahi[(length(rownames(Sortdatahi))-(10*length((Glutamine_1_new)))/2+1):length(rownames(Sortdatahi)),])
    logene<-rownames(Sortdatalo[1:(10*length((Glutamine_1_new))/2),])
    k<-c(higene,logene)
    control<-FetchData(pbmc,vars = k )
    for (p in 1:length(rownames(control))) {
      control$Average[p]<-mean(as.numeric(control[p,]))
    }
    
    pbmc@meta.data$plasmacell<-Glutamine1$Average-control$Average
    remove(p,control,Data,Data_T,Glutamine_1,Glutamine1,rownamesDataT,Sortdata,Sortdatahi,Sortdatalo,Total,average_Glutamine1,Glutamine_1_new,higene,logene,k)
    pdf(paste("Kdm6b",s,t,".pdf",sep = "_"))
    n1<-DotPlot(pbmc,features = c("Kdm6b"))
    print(n1)
    dev.off()
    pdf(paste("Ly75",s,t,".pdf",sep = "_"))
    n1<-DotPlot(pbmc,features = c("Ly75"))
    print(n1)
    dev.off()
    pdf(paste("Irf4",s,t,".pdf",sep = "_"))
    n1<-DotPlot(pbmc,features = c("Irf4"))
    print(n1)
    dev.off()
    pdf(paste("plasmacell",s,t,".pdf",sep = "_"))
    n1<-DotPlot(pbmc,features = c("plasmacell"))
    print(n1)
    dev.off()
  }
}
