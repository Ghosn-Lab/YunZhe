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

#without merging
FL.B<- Read10X(data.dir = "/Users/yunzhehuang/Downloads/Liversample2/",gene.column = 2)
LiverSeurat <- CreateSeuratObject(counts = FL.B, project = "Liver2")

LiverSeurat <- PercentageFeatureSet(LiverSeurat, pattern = "^MT-", col.name = "percent.mt")

# run sctransform
LiverSeurat <- SCTransform(LiverSeurat, vars.to.regress = "percent.mt", verbose = FALSE)


LiverSeurat <- RunPCA(LiverSeurat, verbose = FALSE)

LiverSeurat <- RunTSNE(LiverSeurat, dims = 1:15, method = "FIt-SNE")
LiverSeurat <- FindNeighbors(LiverSeurat, dims = 1:15)
LiverSeurat <- FindClusters(LiverSeurat, resolution = 0.8)
TSNEPlot(LiverSeurat, methods = "FIt-SNE")

FeaturePlot(LiverSeurat, features = c("NR4A1", "CD9"), ncol = 2)
FeaturePlot(LiverSeurat, features = "NR4A1", ncol = 2)

Blood.B <- Read10X(data.dir = "/Users/yunzhehuang/Downloads/Bloodsample2/",gene.column = 2)
BloodSeurat <- CreateSeuratObject(counts = Blood.B, project = "Blood2")
BloodSeurat <- PercentageFeatureSet(BloodSeurat, pattern = "^MT-", col.name = "percent.mt")

# run sctransform
BloodSeurat <- SCTransform(BloodSeurat, vars.to.regress = "percent.mt", verbose = FALSE)


BloodSeurat <- RunPCA(BloodSeurat, verbose = FALSE)

BloodSeurat <- RunTSNE(BloodSeurat, dims = 1:15, method = "FIt-SNE")
BloodSeurat <- FindNeighbors(BloodSeurat, dims = 1:15)
BloodSeurat <- FindClusters(BloodSeurat, resolution = 0.8)
TSNEPlot(BloodSeurat, methods = "FIt-SNE")

FeaturePlot(BloodSeurat, features = c("NR4A1", "CD9"), ncol = 2)

FeaturePlot(BloodSeurat, features = "CD19", ncol = 2)



BoneMarrow.B  <- Read10X(data.dir = "/Users/yunzhehuang/Downloads/Bonemarrowsample2/",gene.column = 2)
BoneMarrowSeurat <- CreateSeuratObject(counts = BoneMarrow.B, project = "BoneMarrow2")

BoneMarrowSeurat <- PercentageFeatureSet(BoneMarrowSeurat, pattern = "^MT-", col.name = "percent.mt")

# run sctransform
BoneMarrowSeurat <- SCTransform(BoneMarrowSeurat, vars.to.regress = "percent.mt", verbose = FALSE)


BoneMarrowSeurat <- RunPCA(BoneMarrowSeurat, verbose = FALSE)

BoneMarrowSeurat <- RunTSNE(BoneMarrowSeurat, dims = 1:15, method = "FIt-SNE")
BoneMarrowSeurat <- FindNeighbors(BoneMarrowSeurat, dims = 1:15)
BoneMarrowSeurat <- FindClusters(BoneMarrowSeurat, resolution = 0.8)
TSNEPlot(BoneMarrowSeurat, methods = "FIt-SNE")

FeaturePlot(BoneMarrowSeurat, features = "NR4A1", ncol = 2)

BEI.A<- Read10X(data.dir = "/Users/yunzhehuang/Downloads/BEI/",gene.column = 2)
BEI.RNA <- BEI.A[["Gene Expression"]]
BEI <- CreateSeuratObject(counts = BEI.RNA, project = "Blood")
BEI <- PercentageFeatureSet(BEI, pattern = "^MT-", col.name = "percent.mt")

# run sctransform
BEI <- SCTransform(BEI, vars.to.regress = "percent.mt", verbose = FALSE)


BEI <- RunPCA(BEI, verbose = FALSE)

BEI <- RunTSNE(BEI, dims = 1:15, method = "FIt-SNE")
BEI <- FindNeighbors(BEI, dims = 1:15)
BEI <- FindClusters(BEI, resolution = 0.8)
TSNEPlot(BEI, methods = "FIt-SNE")


custom_colours <- c("#ffff33", "#fed976", "#ffd700", "#feb24c", "#fe9929", "#fd8d3c", "#ff7f00", "#f16913", "#ec7014", "#d94801", "#fc4e2a", "#e31a1c", "#e41a1c", "#bd0026", "#800026", "#7f0000", "#67000d", "#000000")

FeaturePlot(BEI, features = c("CD19"), ncol = 2, min.cutoff = 0, max.cutoff = 5, cols=custom_colours)




FeaturePlot(BEI, features = "CD19", ncol = 2)

FL.A<- Read10X(data.dir = "/Users/yunzhehuang/Downloads/FLADT/",gene.column = 2)
FL.RNA <- FL.A[["Gene Expression"]]
FL <- CreateSeuratObject(counts = FL.RNA, project = "Liver")

FL <- PercentageFeatureSet(FL, pattern = "^MT-", col.name = "percent.mt")

# run sctransform
FL <- SCTransform(FL, vars.to.regress = "percent.mt", verbose = FALSE)


FL <- RunPCA(FL, verbose = FALSE)

FL <- RunTSNE(FL, dims = 1:15, method = "FIt-SNE")
FL <- FindNeighbors(FL, dims = 1:15)
FL <- FindClusters(FL, resolution = 0.8)
TSNEPlot(FL, methods = "FIt-SNE")



custom_colours <- c("#ffff33", "#fed976", "#ffd700", "#feb24c", "#fe9929", "#fd8d3c", "#ff7f00", "#f16913", "#ec7014", "#d94801", "#fc4e2a", "#e31a1c", "#e41a1c", "#bd0026", "#800026", "#7f0000", "#67000d", "#000000")

FeaturePlot(FL, features = c("NR4A1"), ncol = 2, min.cutoff = 0, max.cutoff = 5, cols=custom_colours)





FeaturePlot(FL, features = "CD19", ncol = 2)


mat.rna <- Read10X(data.dir = "/Users/yunzhehuang/Downloads/BMRNA/expression/")
BMM <- CreateSeuratObject(counts = mat.rna, project = "BoneMarrow")

BMM <- PercentageFeatureSet(BMM, pattern = "^MT-", col.name = "percent.mt")

# run sctransform
BMM <- SCTransform(BMM, vars.to.regress = "percent.mt", verbose = FALSE)


BMM <- RunPCA(BMM, verbose = FALSE)

BMM <- RunTSNE(BMM, dims = 1:15, method = "FIt-SNE")
BMM <- FindNeighbors(BMM, dims = 1:15)
BMM <- FindClusters(BMM, resolution = 0.8)
TSNEPlot(BMM, methods = "FIt-SNE")



custom_colours <- c("#ffff33", "#fed976", "#ffd700", "#feb24c", "#fe9929", "#fd8d3c", "#ff7f00", "#f16913", "#ec7014", "#d94801", "#fc4e2a", "#e31a1c", "#e41a1c", "#bd0026", "#800026", "#7f0000", "#67000d", "#000000")

FeaturePlot(BMM, features = c("CD9"), ncol = 2, min.cutoff = 0, max.cutoff = 5, cols=custom_colours)








FeaturePlot(BMM, features = "CD19", ncol = 2)

Bloodmerged <- merge(BEI, y = c(FL, BMM), add.cell.ids = c("Blood", "Liver", "BoneMarrow"), project = "BloodSample", merge.data = TRUE)
Bloodmerged <- PercentageFeatureSet(Bloodmerged, pattern = "^MT-", col.name = "percent.mt")

Bloodmerged2 <- merge(LiverSeurat, y = c(BloodSeurat, BoneMarrowSeurat), add.cell.ids = c("Blood", "Liver", "BoneMarrow"), project = "BloodSample", merge.data = TRUE)



Bloodmerged <- SCTransform(Bloodmerged, vars.to.regress = "percent.mt", verbose = FALSE)
Bloodmerged2 <- SCTransform(Bloodmerged2, verbose = FALSE)

BloodList <- SplitObject(Bloodmerged2)

for (i in 1:length(x = BloodList)) {
  BloodList[[i]] <- SCTransform(object = BloodList[[i]], verbose = FALSE)
  
}

Blood.features <- SelectIntegrationFeatures(object.list = BloodList, nfeatures = 4000)

BloodList <- PrepSCTIntegration(object.list = BloodList, anchor.features = Blood.features, 
                                verbose = FALSE)





BloodAnchors <- FindIntegrationAnchors(object.list = BloodList, dims = 1:20, normalization.method = "SCT", anchor.features = Blood.features)

BloodIntegrated2 <- IntegrateData(anchorset = BloodAnchors, normalization.method = "SCT",dims = 1:20, verbose = FALSE)

DefaultAssay(object = BloodIntegrated) <- "integrated"


BloodIntegrated2 <- RunPCA(BloodIntegrated2, verbose = FALSE)
BloodIntegrated2 <- RunTSNE(BloodIntegrated2, dims = 1:15, method = "FIt-SNE")
BloodIntegrated2 <- FindNeighbors(BloodIntegrated2, dims = 1:15)
BloodIntegrated2 <- ScaleData(BloodIntegrated2)
BloodIntegrated2 <- FindClusters(BloodIntegrated2, resolution = 0.8)
TSNEPlot(BloodIntegrated2, methods = "FIt-SNE", group.by = "orig.ident")
TSNEPlot(Bloodmerged, methods = "FIt-SNE")







Bloodmerged <- RunPCA(Bloodmerged, verbose = FALSE)
Bloodmerged <- RunTSNE(Bloodmerged, dims = 1:15, method = "FIt-SNE")
Bloodmerged <- FindNeighbors(Bloodmerged, dims = 1:15)
Bloodmerged <- FindClusters(Bloodmerged, resolution = 0.8)
TSNEPlot(Bloodmerged, methods = "FIt-SNE", group.by = "orig.ident")
TSNEPlot(Bloodmerged, methods = "FIt-SNE")



custom_colours <- c("#ffff33", "#fed976", "#ffd700", "#feb24c", "#fe9929", "#fd8d3c", "#ff7f00", "#f16913", "#ec7014", "#d94801", "#fc4e2a", "#e31a1c", "#e41a1c", "#bd0026", "#800026", "#7f0000", "#67000d", "#000000")

FeaturePlot(B2, features = c("CD9"), ncol = 2, min.cutoff = 0, max.cutoff = 5, cols=custom_colours)





FeaturePlot(BloodIntegrated2, features = "NR4A1", ncol = 2)







BloodIntegrated <- RunPCA(BloodIntegrated, verbose = FALSE)
BloodIntegrated <- RunTSNE(BloodIntegrated, dims = 1:15, method = "FIt-SNE")
BloodIntegrated <- FindNeighbors(BloodIntegrated, dims = 1:15)
BloodIntegrated <- FindClusters(BloodIntegrated, resolution = 0.8)
TSNEPlot(BloodIntegrated, methods = "FIt-SNE", group.by = "orig.ident")
TSNEPlot(BloodIntegrated, methods = "FIt-SNE", label = TRUE)

FeaturePlot(BloodIntegrated, features = "CD9", ncol = 2)

DimPlot(projectB.merged, assay = "RNA", reduction = "FIt-SNE", label = FALSE, group.by = "orig.ident", cols = c('BoneMarrow2' = 'green', 'Liver2' = 'grey', 'Blood2' = 'grey'), pt.size = 2 )

B2 <- SubsetData(BloodIntegrated, ident.use = c("0","1","6","20","14","11"))

B2 <- SCTransform(B2, verbose = FALSE)
B2 <- RunPCA(B2, verbose = FALSE)
B2 <- FindNeighbors(B2, dims = 1:15)
B2 <- FindClusters(B2, resolution = 0.8)
B2 <- RunTSNE(B2, dims = 1:15, method = "FIt-SNE")



TSNEPlot(B2, methods = "FIt-SNE", group.by = "orig.ident")
TSNEPlot(B2, methods = "FIt-SNE", label = TRUE)
FeaturePlot(B2, features = "NR4A1", ncol = 2)
