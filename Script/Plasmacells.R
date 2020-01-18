library(dplyr)
library(Seurat)
library(Matrix)
library(reticulate)
library(devtools)
library(ggplot2)
library(ggridges)
py_config()
py_install("umap-learn")

options(future.globals.maxSize = 4000 * 1024^2)


#Loading Blood sample
BEI.A<- Read10X(data.dir = "/Users/yunzhehuang/Downloads/BEI/",gene.column = 2)
BEI.RNA <- BEI.A[["Gene Expression"]]
BEI <- CreateSeuratObject(counts = BEI.RNA, project = "Blood")

BEI[["ADT"]] <- CreateAssayObject(counts = BEI.A[["Antibody Capture"]])






BEI <- PercentageFeatureSet(BEI, pattern = "^MT-", col.name = "percent.mt")

# run sctransform
BEI <- SCTransform(BEI, vars.to.regress = "percent.mt", verbose = FALSE, variable.features.n = 4000)
BEI <- RunPCA(BEI, verbose = FALSE)
BEI <- FindNeighbors(BEI, dims = 1:15)
BEI <- FindClusters(BEI, resolution = 0.8)
BEI <- RunTSNE(BEI, dims = 1:15, method = "FIt-SNE")

TSNEPlot(Bt, methods = "FIt-SNE", label = TRUE)

FeaturePlot(Bt, features = "CD79A", ncol = 2)


BEI[["percent.CD79A"]] <- PercentageFeatureSet(BEI, pattern = "CD79A")
Bt <- subset(BEI, subset = percent.CD79A > 0)

Bcells <- subset(Bt, idents = c('0','2','3','8','9','13','11','6','7'))

Bcells <- SCTransform(Bcells, verbose = FALSE, variable.features.n = 4000)
BEI <- RunPCA(BEI, verbose = FALSE)
BEI <- FindNeighbors(BEI, dims = 1:15)
BEI <- FindClusters(BEI, resolution = 0.8)
BEI <- RunTSNE(BEI, dims = 1:15, method = "FIt-SNE")







FeaturePlot(BEI,features = "", ncol = 2)

custom_colours <- c("#ffff33", "#fed976", "#ffd700", "#feb24c", "#fe9929", "#fd8d3c", "#ff7f00", "#f16913", "#ec7014", "#d94801", "#fc4e2a", "#e31a1c", "#e41a1c", "#bd0026", "#800026", "#7f0000", "#67000d", "#000000")

FeaturePlot(Bt, features = c("CD79A"), ncol = 2, min.cutoff = 0, max.cutoff = 5, cols=custom_colours)

"SDC1", "XBP1", "IRF4", "CD9", "CD38", "ITGB7", "VCAM1"


which(BEI@raw.d["CD79A",]>1
B <- FilterCells(BEI, BEI@assays[["SCT"]]["CD79A"], low.thresholds = 0.1 )

matrix_dir = "/Users/yunzhehuang/Downloads/samplebm/umi_count/"
barcode.path <- paste0(matrix_dir, "barcodes.tsv.gz")
features.path <- paste0(matrix_dir, "features.tsv.gz")
matrix.path <- paste0(matrix_dir, "matrix.mtx.gz")
mat <- readMM(file = matrix.path)
feature.names = read.delim(features.path, header = FALSE, stringsAsFactors = FALSE)
barcode.names = read.delim(barcode.path, header = FALSE, stringsAsFactors = FALSE)
colnames(mat) = barcode.names$V1
rownames(mat) = feature.names$V1
mat.rna <- Read10X(data.dir = "/Users/yunzhehuang/Downloads/BMRNA/expression/")

BMM <- CreateSeuratObject(counts = mat.rna, project = "BoneMarrow")
BMM[["ADT"]] <- CreateAssayObject(counts = mat)

BMM <- PercentageFeatureSet(BMM, pattern = "^MT-", col.name = "percent.mt")

# run sctransform
BMM <- SCTransform(BMM, vars.to.regress = "percent.mt", verbose = FALSE, variable.features.n = 4000)
BMM <- RunPCA(BMM, verbose = FALSE)
BMM <- FindNeighbors(BMM, dims = 1:15)
BMM <- FindClusters(BMM, resolution = 0.8)
BMM <- RunTSNE(BMM, dims = 1:15, method = "FIt-SNE")



BMMBcells <- subset(BMM, idents = c('7','0','5','4','9','6','1','2'))
BMMBcells <- SCTransform(BMMBcells, verbose = FALSE, variable.features.n = 4000)
 <- RunPCA(BEI, verbose = FALSE)
BEI <- FindNeighbors(BEI, dims = 1:15)
BEI <- FindClusters(BEI, resolution = 0.8)
BEI <- RunTSNE(BEI, dims = 1:15, method = "FIt-SNE")

