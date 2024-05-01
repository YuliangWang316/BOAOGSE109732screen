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


pbmc.data <- Read10X(data.dir = "D:/B18 cell Transfer&prdm1 reporter single cell sequencing Rawdata/CJJS1/outs/filtered_feature_bc_matrix/")

pbmc <- CreateSeuratObject(counts = pbmc.data, project = "pbmc3k", min.cells = 3, min.features = 200)
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^mt-")
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,pt.size = 0)
pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 7500 & percent.mt < 25)
pbmc <- NormalizeData(pbmc)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(pbmc)
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
s.genes<-capitalize(tolower(s.genes))
g2m.genes<-capitalize(tolower(g2m.genes))
pbmc <- CellCycleScoring(pbmc, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
pbmc <- ScaleData(pbmc, features = all.genes, verbose=F)
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
pbmc <- FindNeighbors(pbmc, dims = 1:10)
pbmc <- FindClusters(pbmc, resolution = 0.5)
pbmc <- RunUMAP(pbmc, dims = 1:10)
pbmc <- RunTSNE(pbmc, dims = 1:10)
DimPlot(pbmc,reduction = "umap",label = TRUE)
DimPlot(pbmc,reduction = "tsne",label = TRUE)
remove(all.genes,pbmc.data)
remove(g2m.genes,s.genes)
gc()
pbmc2<-readRDS("e:/CJJGSE109732.rds")
FeaturePlot(pbmc,features = c("Prdm1","Irf4"),order = TRUE,label = TRUE)
Idents(pbmc)<-pbmc$seurat_clusters
pbmc_pc<-subset(pbmc,idents = c("0","1","6","7"))
pbmc_new<-merge(x=pbmc2,y = pbmc_pc)
pbmc_new <- NormalizeData(pbmc_new)
pbmc_new <- FindVariableFeatures(pbmc_new, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(pbmc_new)
pbmc_new <- ScaleData(pbmc_new, features = all.genes)
pbmc_new <- RunPCA(pbmc_new, features = VariableFeatures(object = pbmc_new))
pbmc_new$Cluster_merge<-Idents(pbmc_new)
pbmc_new <- FindNeighbors(pbmc_new, dims = 1:30)
pbmc_new <- FindClusters(pbmc_new, resolution = 1.2)
pbmc_new <- RunUMAP(pbmc_new, dims = 1:30)
pbmc_new <- RunTSNE(pbmc_new, dims = 1:30)
DimPlot(pbmc_new,reduction = "umap",label = TRUE)
DimPlot(pbmc_new,reduction = "tsne",label = TRUE)
Idents(pbmc_new)<-pbmc_new$Cluster_merge
remove(pbmc,pbmc_pc,pbmc2,all.genes)
gc()

FeaturePlot(pbmc_new,features = "Prdm1",order = TRUE)
new.cluster.ids <- c("PC", "PC", "PC", "PC", "DZ", "DZ",
                     "DZ", "LZ", "LZ","LZ","LZ","LZ","LZ","prePB")
names(new.cluster.ids) <- levels(pbmc_new)
pbmc_new <- RenameIdents(pbmc_new, new.cluster.ids)
pbmc_new$Cluster_merge_NEW<-Idents(pbmc_new)

pbmc_new<-RunLDA(pbmc_new,labels = pbmc_new$Cluster_merge_NEW)
pbmc_new<-RunUMAP(pbmc_new,reduction = "lda",reduction.name = "ldaumap",dims = 1:3)
pbmc_new<-RunTSNE(pbmc_new,reduction = "lda",reduction.name = "ldatsne",dims = 1:3)
DimPlot(pbmc_new,reduction = "ldaumap",label = TRUE)
DimPlot(pbmc_new,reduction = "ldatsne",label = TRUE)

gc()
saveRDS(pbmc_new,"E:/boao_b18.rds")
remove(new.cluster.ids)
gc()
library(monocle)
trace('project2MST',edit = T,where = asNamespace("monocle"))
data<-as.sparse(pbmc_new@assays$RNA@counts)
pd <-pbmc_new@meta.data
pd <- new('AnnotatedDataFrame', data = pd)  
fData<-data.frame(gene_short_name = row.names(data),row.names = row.names(data))
fd <- new('AnnotatedDataFrame', data = fData)
monocle_cds <- newCellDataSet(data, phenoData = pd,featureData = fd,lowerDetectionLimit = 0.1,
                              expressionFamily = VGAM::negbinomial.size())
remove(data,fd,fData,pbmc_new,pd)
gc()
monocle_cds <- estimateSizeFactors(monocle_cds)
monocle_cds <- estimateDispersions(monocle_cds)
monocle_cds <- detectGenes(monocle_cds, min_expr = 0.1)
# 
# 
# 
pbmc.marker<-FindAllMarkers(pbmc_new_new,only.pos = TRUE,min.pct = 0)
diff_test_res <- differentialGeneTest(monocle_cds,fullModelFormulaStr = "~Cluster_merge_NEW",cores = 20)

pbmcmarkers_new<-pbmc.marker[which(pbmc.marker$p_val_adj < 0.05 ),]
pbmcmarkers_new$filter<-pbmcmarkers_new$pct.1 - pbmcmarkers_new$pct.2
# ordering_genes <- diff_test_res$gene_short_name
ordering_genes <- row.names (subset(diff_test_res, qval < 1e-300))
monocle_cds <-setOrderingFilter(monocle_cds,ordering_genes = ordering_genes)
monocle_cds <-reduceDimension(monocle_cds,reduction_method = "DDRTree",max_components = 2)
gc()
saveRDS(monocle_cds,"e:/boao_b18_monocle_raw.rds")
diff_test_res<-cbind(rownames(diff_test_res),diff_test_res)
colnames(diff_test_res)[1]<-"Gene_Symbol"
write.table(diff_test_res,"e:/diff_test_Res_boao_b18.txt",sep = "\t",row.names = FALSE,quote = FALSE)
monocle_cds <-orderCells(monocle_cds)
plot_cell_trajectory(monocle_cds, color_by = "State",cell_size = 0.75) 
plot_cell_trajectory(monocle_cds, color_by = "Pseudotime",cell_size = 0.75)
plot_cell_trajectory(monocle_cds, color_by = "Cluster_merge_NEW",cell_size = 0.75,)
plot_cell_trajectory(monocle_cds, color_by = "Cluster_merge_NEW",cell_size = 0.75,)+ facet_wrap(~Cluster_merge_NEW, nrow = 1)

my_genes <- row.names(subset(fData(monocle_cds),
                             gene_short_name %in% c("Irf4","Ly75","Kdm6b","Prdm1")))
cds_subset <- monocle_cds[my_genes,]
# plot_genes_in_pseudotime(cds_subset, color_by = "group")
plotdf=pData(monocle_cds)
library(ggridges)
mycolor<-c("#619CFF","#00BA38","#F8766D","Black")
ggplot(plotdf, aes(x=Pseudotime,y=Cluster_merge_NEW,fill=Cluster_merge_NEW))+
  geom_density_ridges(scale=1) +
  geom_vline(xintercept = c(5,10),linetype=2)+
  scale_y_discrete("")+
  theme_minimal()+
  theme(
    panel.grid = element_blank()
  )+scale_fill_manual(values = mycolor)
saveRDS(monocle_cds,file = "e:/boao_b18_monocle.rds")
