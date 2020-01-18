library(dplyr)
library(Seurat)
library(Matrix)
library(reticulate)
py_config()
py_install("umap-learn")

sample5.adt <- Read10X(data.dir = "/Users/yunzhehuang/Downloads/Antibody-seq/0005_wl/umi_count/",gene.column = 1)
sample5.rna <- Read10X(data.dir = "/Users/yunzhehuang/Downloads/RNA-seq-Cellranger_out/19CT019936/filtered_feature_bc_matrix/")
sample5 <- CreateSeuratObject(counts = sample5.rna, project = "sample5")
sample5 <- NormalizeData(sample5)
sample5 <- FindVariableFeatures(sample5)
sample5 <- ScaleData(sample5)
sample5[["ADT"]] <- CreateAssayObject(counts = sample5.adt)

sample6.adt <- as.sparse(read.csv(file = "/Users/yunzhehuang/Downloads/Antibody-seq/0006_wl/umi_count/final.csv",sep = "\t", header = TRUE, row.names = 1))
sample6.rna <- Read10X(data.dir = "/Users/yunzhehuang/Downloads/RNA-seq-Cellranger_out/19CT019938/filtered_feature_bc_matrix")
sample6 <- CreateSeuratObject(counts = sample6.rna, project = "sample6")
sample6 <- NormalizeData(sample6)
sample6 <- FindVariableFeatures(sample6)
sample6 <- ScaleData(sample6)
sample6[["ADT"]] <- CreateAssayObject(counts = sample6.adt)
sample7.adt <- as.sparse(read.csv(file = "/Users/yunzhehuang/Downloads/Antibody-seq/0007_wl/umi_count/final.csv",sep = "\t", header = TRUE, row.names = 1))
sample7.rna <- Read10X(data.dir = "/Users/yunzhehuang/Downloads/RNA-seq-Cellranger_out/19CT019939/filtered_feature_bc_matrix")
sample7 <- CreateSeuratObject(counts = sample7.rna, project = "sample7")
sample7 <- NormalizeData(sample7)
sample7 <- FindVariableFeatures(sample7)
sample7 <- ScaleData(sample7)
sample7[["ADT"]] <- CreateAssayObject(counts = sample7.adt)


sample.merged <- merge(sample5, y = c(sample6, sample7), add.cell.ids = c("5", "6", "7"), project = "Merged", merge.data = TRUE)
head(colnames(sample.merged))


sample.merged[["percent.mt"]] <- PercentageFeatureSet(sample.merged, pattern = "^MT-")
sample.merged <- subset(sample.merged, percent.mt < 5)
sample.merged <- FindVariableFeatures(sample.merged)
sample.merged <- ScaleData(sample.merged)

sample.merged <- RunPCA(sample.merged, features = VariableFeatures(object = sample.merged))
PCAPlot(object = sample.merged, dim.1 = 1, dim.2 = 2)
PCHeatmap(sample.merged, use.pcs = 2)

#T-reg
FeaturePlot(sample.merged, features = c("adt_CD19", "adt_CD10", "adt_CD24"))

FeaturePlot(sample.merged, features = c("S100A10"))
FeaturePlot(sample.merged, features = c("PRDM1", "PAX5", "XBP1","CD38", "ITGB7"))
grep('FOXP3', rownames(sample.merged@raw.data))
which(FetchData(sample.merged, vars = 'FOXP3', slot = 'data') > 1)
my.data <- FetchData(sample.merged,c("FOXP3"))
tail(my.data,5)

Treg1[["percent.FOXP3"]] <- PercentageFeatureSet(sample.merged, pattern = "FOXP3")
Treg2 <- subset(sample.merged,subset = percent.FOXP3 > 50)


sample.merged <- RunUMAP(sample.merged, npcs = 20, reduction = "pca",verbose = FALSE, dims = 1:10, resolution = 0.8)
sampleplot <- DimPlot(sample.merged, reduction = "umap") 
sampleplot + DarkTheme()

sample.merged <- FindNeighbors(sample.merged, dims = 1:10)
sample.merged <- FindClusters(sample.merged, redection.type = "pca", dims.use = 1:20,resolution = 0.8)
sample.merged <- RunTSNE(sample.merged, dims = 1:20, method = "FIt-SNE")
TSNEPlot(object = sample.merged, do.label = TRUE, group.by = "orig.ident")

sample.merged.rna.markers <- FindAllMarkers(sample.merged, max.cells.per.ident = 100, min.diff.pct = 0.3, only.pos = TRUE)
new.cluster.ids <- c("Memory CD4 T", "CD14+ Mono", "Naive CD4 T", "NK", "CD14+ Mono","B", 
                     "CD8 T", "CD16+ Mono", "T/Mono doublets", "NK", "CD34+", "Multiplets", "Eryth", "Mk", 
                     "DC", "pDCs","Mouse", "HS_Stress","FCGR3A+ Monocytes","Dendritic Cells","Megakaryocytes","CD14+ Monocytes")
names(new.cluster.ids) <- levels(sample.merged)
sample.merged <- RenameIdents(sample.merged, new.cluster.ids)  
cluster.markers <- FindAllMarkers(object = sample.merged, ident.1 = 1, min.pct = 0.25)
TSNEPlot(object = sample.merged, do.label = TRUE)
AddModuleScore(sample.merged, name = "Cluster")

cluster.markers <- FindAllMarkers(object = sample.merged, ident.1 = 1, min.pct = 0.25)
sample.merged <- RenameIdents(sample.merged, new.cluster.ids)
DimPlot(sample.merged, reductiion = "umap")
DoHeatmap(sample.merged, features = unlist(TopFeatures(sample.merged[["pca"]], balanced = TRUE)), size = 3)
WhichCells(sample.merged, idents = "B")
B.raw.data <- as.matrix(GetAssayData(sample.merged, slot = "counts")[, WhichCells(sample.merged, ident = "B")])
sample.B <- subset(sample.merged, idents = c("B"))
sample.B <- FindVariableFeatures(sample.B)
sample.B <- ScaleData(sample.B)
sample.B <- RunPCA(sample.B, features = VariableFeatures(object = sample.B))

sample.B <- RunUMAP(sample.B, npcs = 20, reduction = "pca",verbose = FALSE, dims = 1:10, resolution = 0.8)
DimPlot(sample.B, reduction = "umap")
sample.B <- FindNeighbors(sample.B, dims = 1:10)
sample.B <- FindClusters(sample.B, redection.type = "pca", dims.use = 1:20,resolution = 0.8)
sample.B <- RunTSNE(sample.B, dims = 1:20, method = "FIt-SNE")
TSNEPlot(object = sample.B, do.label = TRUE)

sample.B.rna.markers <- FindAllMarkers(sample.B, max.cells.per.ident = 100, min.diff.pct = 0.3, only.pos = TRUE)
new.cluster.ids <- c("Memory B", "Naive B", "Plasma", "Unknown")
names(new.cluster.ids) <- levels(sample.merged)
sample.B <- RenameIdents(sample.B, new.cluster.ids)  
cluster.markers <- FindAllMarkers(object = sample.B, ident.1 = 1, min.pct = 0.25)
TSNEPlot(object = sample.B, do.label = TRUE)
VlnPlot(sample.B, features = c("WARS"))
FeaturePlot(sample.B, features = c("WARS", "NME3"))

clustersample.markers <- FindMarkers(sample.merged, ident.1 = 16)
head(clustersample.markers, n = 10)

