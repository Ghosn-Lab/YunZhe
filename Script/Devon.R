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

#Reading GFP ADT data
GFP.adt <- as.sparse(read.csv(file = "/Users/yunzhehuang/Desktop/Devon/GFP_ADT/umi_count/final.csv",sep = "\t", header = TRUE, row.names = 1))
#Reading RNA data
GFP.rna <- Read10X(data.dir = "/Users/yunzhehuang/Desktop/Devon/DE_GFP_GEX",gene.column = 2 )
GFP <- CreateSeuratObject(counts = GFP.rna, project = "GFP")
GFP[["ADT"]] <- CreateAssayObject(counts = GFP.adt)




#Reading Total data
Total.raw <- Read10X(data.dir = "/Users/yunzhehuang/Desktop/Devon/Total/",gene.column = 2)
Total.RNA <- Total.raw[["Gene Expression"]]
Total <- CreateSeuratObject(counts = Total.RNA, project = "Total")

Total[["ADT"]] <- CreateAssayObject(counts = Total.raw[["Antibody Capture"]])

Merged <- merge(GFP, y = c(Total), add.cell.ids = c("GFP", "Total"), project = "Merged", merge.data = TRUE)

Merged <- SCTransform(Merged, verbose = TRUE)

Merged <- NormalizeData(Merged, assay = "ADT", normalization.method = "CLR", margin = 2)
Merged <- ScaleData(Merged, assay = "ADT")


Merged <- RunPCA(Merged, verbose = TRUE)
ElbowPlot(Merged)

DimHeatmap(Merged, dims = 1:10)

Merged <- JackStraw(Merged, num.replicate = 100)
Merged <- ScoreJackStraw(Merged, dims = 1:20)

JackStrawPlot(Merged, dims = 1:15)

Merged <- FindNeighbors(Merged, dims = 1:7)
Merged <- FindClusters(Merged, resolution = 0.5)
Merged <- RunTSNE(Merged, dims = 1:7, method = "FIt-SNE")
TSNEPlot(Merged, methods = "FIt-SNE", group.by = "orig.ident")
TSNEPlot(Merged, methods = "FIt-SNE", label = TRUE)

custom_colours <- c("#ffff33", "#fed976", "#ffd700", "#feb24c", "#fe9929", "#fd8d3c", "#ff7f00", "#f16913", "#ec7014", "#d94801", "#fc4e2a", "#e31a1c", "#e41a1c", "#bd0026", "#800026", "#7f0000", "#67000d", "#000000")

FeaturePlot(Merged, features = c("CD45"), ncol = 2, min.cutoff = 0, max.cutoff = 5, cols=custom_colours)

FeaturePlot(Merged, features = c("Ccr2", "Olr1"), min.cutoff = 0, max.cutoff = 5)
FeaturePlot(Merged, features = c("Il1b", "Cxcl2"), min.cutoff = 0, max.cutoff = 5)
FeaturePlot(Merged, features = c("Lrp2", "Pirb"), min.cutoff = 0, max.cutoff = 5)
FeaturePlot(Merged, features = c("Gipc1", "Cd163"), min.cutoff = 0, max.cutoff = 5)
FeaturePlot(Merged, features = c("Mrc1", "Itgax"), min.cutoff = 0, max.cutoff = 5)
FeaturePlot(Merged, features = c("Sirpa", "Xcr1"), min.cutoff = 0, max.cutoff = 5)
FeaturePlot(Merged, features = c("Cd74", "Ly6c1"), min.cutoff = 0, max.cutoff = 5)
FeaturePlot(Merged, features = c("Ptprc", "Cx3cr1"), min.cutoff = 0, max.cutoff = 5)
FeaturePlot(Merged, features = c("Adgre1", "Cd80"), min.cutoff = 0, max.cutoff = 5)
FeaturePlot(Merged, features = c("Cd86", "Cd40"), min.cutoff = 0, max.cutoff = 5)
FeaturePlot(Merged, features = c("Dpp4", "Cd164"), min.cutoff = 0, max.cutoff = 5)
FeaturePlot(Merged, features = c("Cybb", "Car4"), min.cutoff = 0, max.cutoff = 5)
FeaturePlot(Merged, features = c("Mafb", "Tlr4"), min.cutoff = 0, max.cutoff = 5)
FeaturePlot(Merged, features = c("Lamp1", "Fcgr1"), min.cutoff = 0, max.cutoff = 5)
FeaturePlot(Merged, features = c("Siglecf", "Siglech"), min.cutoff = 0, max.cutoff = 5)
FeaturePlot(Merged, features = c("Ccr7", "Bcl2"), min.cutoff = 0, max.cutoff = 5)
FeaturePlot(Merged, features = c("Mcl1", "Fadd"), min.cutoff = 0, max.cutoff = 5)
FeaturePlot(Merged, features = c("Cflar", "Casp1"), min.cutoff = 0, max.cutoff = 5)
FeaturePlot(Merged, features = c("Casp3", "Casp8"), min.cutoff = 0, max.cutoff = 5)
FeaturePlot(Merged, features = c("Casp9", "Pidd1"), min.cutoff = 0, max.cutoff = 5)
FeaturePlot(Merged, features = c("Tnfsf10", "Bbc3"), min.cutoff = 0, max.cutoff = 5)
FeaturePlot(Merged, features = c("Mx1", "Pawr"), min.cutoff = 0, max.cutoff = 5)
FeaturePlot(Merged, features = c("Rela", "Relb"), min.cutoff = 0, max.cutoff = 5)
FeaturePlot(Merged, features = c("Tnf"), min.cutoff = 0, max.cutoff = 5)
FeaturePlot(Merged, features = c("Pira2"), min.cutoff = 0, max.cutoff = 5)


#Classifying T cells(CD3 ADT, CD8 genes)
FeaturePlot(Merged, features = c("CD3"), min.cutoff = 0, max.cutoff = 5)
#Classifying CD8 and Cd4 T cells
FeaturePlot(Merged, features = c("Cd8a", "Cd4"), min.cutoff = 0, max.cutoff = 5)

#NKT
FeaturePlot(Merged, features = c("Klrb1c"), min.cutoff = 0, max.cutoff = 5)

#NK cells
FeaturePlot(Merged, features = c("Klrb1b"), min.cutoff = 0, max.cutoff = 5)

#Plasma cells
FeaturePlot(Merged, features = c("Itgb7", "Cd38"), min.cutoff = 0, max.cutoff = 5)

FeaturePlot(Merged, features = c("CD19", "Cd79a"), min.cutoff = 0, max.cutoff = 5)

new.cluster.ids <- c("CD4+ Tcells", "CD4+ Tcells", "Myeloid Cells", "CD8a+ Tcells", "B Cells/Plasma Cells", "NKT + CD8a+ NK Cells", "Myeloid Cells", "Myeloid Cells", "T cells")
names(new.cluster.ids) <- levels(Merged)
Merged <- RenameIdents(Merged, new.cluster.ids)

TSNEPlot(Merged, methods = "FIt-SNE", label = TRUE)






MyeloidCellTest <- SubsetData(Merged, ident.use = c("2","6","7"))
MyeloidCell <- SCTransform(MyeloidCell, verbose = FALSE, variable.features.n = 4000)
MyeloidCellTest <- RunPCA(MyeloidCellTest, verbose = FALSE)

MyeloidCellTest <- JackStraw(MyeloidCellTest, num.replicate = 100)
MyeloidCellTest <- ScoreJackStraw(MyeloidCellTest, dims = 1:20)

JackStrawPlot(MyeloidCellTest, dims = 1:15)

MyeloidCellTest <- FindNeighbors(MyeloidCellTest, dims = 1:11)
MyeloidCellTest <- FindClusters(MyeloidCellTest, resolution = 0.8)
MyeloidCellTest <- RunTSNE(MyeloidCellTest, dims = 1:11, method = "FIt-SNE")
TSNEPlot(MyeloidCellTest, methods = "FIt-SNE", group.by = "orig.ident")
TSNEPlot(MyeloidCellTest, methods = "FIt-SNE", label = TRUE)

FeaturePlot(MyeloidCellTest, features = c("CD3", "CD86"), min.cutoff = 0, max.cutoff = 5)
FeaturePlot(MyeloidCellTest, features = c("Gr-1", "CD11b"), min.cutoff = 0, max.cutoff = 5)
FeaturePlot(MyeloidCellTest, features = c("CD11c", "I-A/I-E"), min.cutoff = 0, max.cutoff = 5)
FeaturePlot(MyeloidCellTest, features = c("CD193", "CD19"), min.cutoff = 0, max.cutoff = 5)
FeaturePlot(MyeloidCellTest, features = c("F4/80", "CD115"), min.cutoff = 0, max.cutoff = 5)
FeaturePlot(MyeloidCellTest, features = c("CD69", "CD335"), min.cutoff = 0, max.cutoff = 5)
FeaturePlot(MyeloidCellTest, features = c("CD45", "Ly6C"), min.cutoff = 0, max.cutoff = 5)
FeaturePlot(MyeloidCellTest, features = c("CD5"), min.cutoff = 0, max.cutoff = 5)


FeaturePlot(MyeloidCellTest, features = c("Ccr2", "Olr1"), min.cutoff = 0, max.cutoff = 5)
FeaturePlot(MyeloidCellTest, features = c("Il1b", "Cxcl2"), min.cutoff = 0, max.cutoff = 5)
FeaturePlot(MyeloidCellTest, features = c("Lrp2", "Pirb"), min.cutoff = 0, max.cutoff = 5)
FeaturePlot(MyeloidCellTest, features = c("Gipc1", "Cd163"), min.cutoff = 0, max.cutoff = 5)
FeaturePlot(MyeloidCellTest, features = c("Mrc1", "Itgax"), min.cutoff = 0, max.cutoff = 5)
FeaturePlot(MyeloidCellTest, features = c("Sirpa", "Xcr1"), min.cutoff = 0, max.cutoff = 5)
FeaturePlot(MyeloidCellTest, features = c("Cd74", "Ly6c1"), min.cutoff = 0, max.cutoff = 5)
FeaturePlot(MyeloidCellTest, features = c("Ptprc", "Cx3cr1"), min.cutoff = 0, max.cutoff = 5)
FeaturePlot(MyeloidCellTest, features = c("Adgre1", "Cd80"), min.cutoff = 0, max.cutoff = 5)
FeaturePlot(MyeloidCellTest, features = c("Cd86", "Cd40"), min.cutoff = 0, max.cutoff = 5)
FeaturePlot(MyeloidCellTest, features = c("Dpp4", "Cd164"), min.cutoff = 0, max.cutoff = 5)
FeaturePlot(MyeloidCellTest, features = c("Cybb", "Car4"), min.cutoff = 0, max.cutoff = 5)
FeaturePlot(MyeloidCellTest, features = c("Mafb", "Tlr4"), min.cutoff = 0, max.cutoff = 5)
FeaturePlot(MyeloidCellTest, features = c("Lamp1", "Fcgr1"), min.cutoff = 0, max.cutoff = 5)
FeaturePlot(MyeloidCellTest, features = c("Siglecf", "Siglech"), min.cutoff = 0, max.cutoff = 5)
FeaturePlot(MyeloidCellTest, features = c("Ccr7", "Bcl2"), min.cutoff = 0, max.cutoff = 5)
FeaturePlot(MyeloidCellTest, features = c("Mcl1", "Fadd"), min.cutoff = 0, max.cutoff = 5)
FeaturePlot(MyeloidCellTest, features = c("Cflar", "Casp1"), min.cutoff = 0, max.cutoff = 5)
FeaturePlot(MyeloidCellTest, features = c("Casp3", "Casp8"), min.cutoff = 0, max.cutoff = 5)
FeaturePlot(MyeloidCellTest, features = c("Casp9", "Pidd1"), min.cutoff = 0, max.cutoff = 5)
FeaturePlot(MyeloidCellTest, features = c("Tnfsf10","Bbc3"), min.cutoff = 0, max.cutoff = 5)
FeaturePlot(MyeloidCellTest, features = c("Mx1", "Pawr"), min.cutoff = 0, max.cutoff = 5)
FeaturePlot(MyeloidCellTest, features = c("Rela", "Relb"), min.cutoff = 0, max.cutoff = 5)
FeaturePlot(MyeloidCellTest, features = c("Tnf"), min.cutoff = 0, max.cutoff = 5)
FeaturePlot(MyeloidCellTest, features = c("Pira2"), min.cutoff = 0, max.cutoff = 5)









