library(dplyr)
library(Seurat)
library(Matrix)
library(reticulate)
py_config()
py_install("umap-learn")



sample5.adt <- Read10X(data.dir = "/Users/yunzhehuang/Downloads/Antibody-seq/0005_wl/umi_count/",gene.column = 1)
sample5.adt[1:21,1:3]
sample5.rna <- Read10X(data.dir = "/Users/yunzhehuang/Downloads/RNA-seq-Cellranger_out/19CT019936/filtered_feature_bc_matrix/")
sample5.rna[1:10,1:5]
sample5 <- CreateSeuratObject(counts = sample5.rna)
sample5[["ADT"]] <- CreateAssayObject(counts = sample5.adt)
sample5[["percent.mt"]] <- PercentageFeatureSet(sample5, pattern = "^MT-")
sample5 <- subset(sample5, subset = percent.mt < 15)
sample5 <- NormalizeData(sample5)
sample5 <- FindVariableFeatures(sample5, selection.method = "vst")
sample5 <- ScaleData(sample5)

RidgePlot(sample6,features=c("adt_CD38","adt_CD3", "adt_CD34", "adt_CD14"), ncol = 2)
sample5 <- NormalizeData(sample5, assay = "ADT", normalization.method = "CLR", margin = 2)
sample5 <- ScaleData(sample5, assay = "ADT")
sample5 <- RunPCA(sample5, features = VariableFeatures(object = sample5))
print(sample5[["pca"]], dims = 1:5, nfeatures = 5)
PCAPlot(object = sample5, dim.1 = 1, dim.2 = 2)
sample5 <- FindNeighbors(sample5, dims = 1:25)
sample5 <- FindClusters(sample5, redection.type = "pca", dims.use = 1:10,resolution = 0.8)
sample5 <- RunUMAP(sample5, assay = "RNA", npcs = 20, reduction = "pca",verbose = FALSE, dims = 1:10, resolution = 0.8)
DimPlot(sample5, assay = "RNA", reduction = "umap", label = TRUE)
RidgePlot(sample5,assay = "ADT",features=c("adt_CD34","adt_CD24", "adt_CD138", "adt_CD3"), ncol = 2)
RidgePlot(sample5, assay = "ADT", features=c("adt_AnnexinV","adt_CD16","adt_IgD","adt_IgM"), ncol = 2)
RidgePlot(sample5, assay = "ADT", features=c("adt_CD56", "adt_CD127","adt_CD38", "adt_CD20"), ncol = 2)
RidgePlot(sample5,assay = "ADT", features=c("adt_CD27", "adt_CD14", "adt_CD19","adt_CD5"), ncol = 2)
RidgePlot(sample5, assay = "ADT", features=c("adt_CD43","adt_CD4", "adt_CD8", "adt_CD25"), ncol = 2)

FeaturePlot(sample5, features = c("adt_AnnexinV","adt_CD19","adt_IgD","adt_IgM"))


sample5 <- RunTSNE(sample5, dims = 1:25, method = "FIt-SNE")
sample5.rna.markers <- FindAllMarkers(sample5, max.cells.per.ident = 100, min.diff.pct = 0.3, only.pos = TRUE)
new.cluster.ids <- c("Memory CD4 T", "CD14+ Mono", "Naive CD4 T", "NK", "CD14+ Mono", "Mouse", "B", 
                     "CD8 T", "CD16+ Mono", "T/Mono doublets", "NK", "CD34+", "Multiplets", "Mouse", "Eryth", "Mk", 
                     "Mouse", "DC", "pDCs")
names(new.cluster.ids) <- levels(sample5)
sample5 <- RenameIdents(sample5, new.cluster.ids)
FeatureScatter(sample5, feature1 ="adt_AnnexinV", feature2 = "MT-ND4")

FeatureScatter(sample5, feature1 = "adt_CD4", feature2 = "adt_CD8")
FeatureScatter(sample5, feature1 = "adt_CD14", feature2 = "adt_CD16")

