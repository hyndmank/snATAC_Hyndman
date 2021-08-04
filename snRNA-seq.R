#this is the workflow for our analyses of the snRNA-seq in 2021 using LogNormalization and DESeq2 for determination of differentially expressed genes. 
library(Seurat)#v4.0.0
library(dplyr)
library(harmony)
library(ggplot2)
library(patchwork)
set.seed(1234)                
library(tidyr)
library(future)
options(future.globals.maxSize = 50000 * 1024^2)
library(DESeq2)

#import the data
Sample67c.data <- Read10X(data.dir = "/data/user/hyndmank/RNAseq/sample67c/outs/filtered_feature_bc_matrix/")
Sample66c.data <- Read10X(data.dir = "/data/user/hyndmank/RNAseq/sample66c/outs/filtered_feature_bc_matrix/")
Sample76c.data <- Read10X(data.dir = "/data/user/hyndmank/RNAseq/sample76c/outs/filtered_feature_bc_matrix/")
Sample77c.data <- Read10X(data.dir = "/data/user/hyndmank/RNAseq/sample77c/outs/filtered_feature_bc_matrix/")

#create seurat objects for each sample
Sample67c <- CreateSeuratObject(counts = Sample67c.data, project = "Sample67c", min.cells = 3, min.features = 500)
Sample67c <- PercentageFeatureSet(Sample67c, pattern = "^mt-", col.name = "percent.mt")
Sample67c <- subset(Sample67c, subset = nFeature_RNA > 500  & percent.mt < 2)
Sample67c

Sample66c <- CreateSeuratObject(counts = Sample66c.data, project = "Sample66c", min.cells = 3, min.features = 500)
Sample66c <- PercentageFeatureSet(Sample66c, pattern = "^mt-", col.name = "percent.mt")
Sample66c <- subset(Sample66c, subset = nFeature_RNA > 500  & percent.mt < 2)
Sample66c

Sample76c <- CreateSeuratObject(counts = Sample76c.data, project = "Sample76c", min.cells = 3, min.features = 500)
Sample76c <- PercentageFeatureSet(Sample76c, pattern = "^mt-", col.name = "percent.mt")
Sample76c <- subset(Sample76c, subset = nFeature_RNA > 500  & percent.mt < 2)
Sample76c

Sample77c <- CreateSeuratObject(counts = Sample77c.data, project = "Sample77c", min.cells = 3, min.features = 500)
Sample77c <- PercentageFeatureSet(Sample77c, pattern = "^mt-", col.name = "percent.mt")
Sample77c <- subset(Sample77c, subset = nFeature_RNA > 500  & percent.mt < 2)
Sample77c

#add genotype and sex to each file
Sample66c$groupid <- "ko"
Sample67c$groupid <- "con"
Sample76c$groupid <- "ko"
Sample77c$groupid <- "con"
Sample66c$sex <- "female"
Sample67c$sex <- "female"
Sample76c$sex <- "male"
Sample77c$sex <- "male"

#merge to create one object
everyone <- merge(
  x = Sample66c,
  y = list(Sample67c, Sample76c, Sample77c))
everyone[["RNA"]]

#log normalize, PCA etc
everyone <-NormalizeData(everyone, normalization.method = "LogNormalize", scale.factor = 10000)
everyone <- FindVariableFeatures(everyone, selection.method = "vst", nfeatures = 2000)
top10 <- head(VariableFeatures(everyone), 10)
all.genes <- rownames(everyone)
everyone <- ScaleData(everyone, features = all.genes)
everyone <- RunPCA(everyone, features = VariableFeatures(object = everyone))
DimHeatmap(everyone, dims = 1:15, cells = 500, balanced = TRUE)
everyone <- JackStraw(everyone, num.replicate = 100)
everyone <- ScoreJackStraw(everyone, dims = 1:20)
JackStrawPlot(everyone, dims = 1:20)
everyone <- FindNeighbors(everyone, dims = 1:20)
everyone <- FindClusters(everyone, resolution = 0.4)
everyone <- RunUMAP(everyone, dims = 1:20)

#Harmony batch correction
combined.harmony <- RunHarmony(
  object = everyone,
  group.by.vars = 'groupid',
  reduction = 'pca',
  assay.use = 'RNA',
  project.dim = FALSE)

combined.harmony <- RunUMAP(combined.harmony, dims = 2:30, reduction = 'harmony')
combined.harmony <- FindNeighbors(combined.harmony, dims = 1:20)
combined.harmony <- FindClusters(combined.harmony, resolution = 0.4)
p1 <- DimPlot(everyone, group.by = 'groupid', pt.size = 0.1) + ggplot2::ggtitle("Unintegrated")
p2 <- DimPlot(combined.harmony, group.by = 'groupid', pt.size = 0.1) + ggplot2::ggtitle("Harmony integration")
p1 + p2
p3 <-DimPlot(combined.harmony, label = TRUE, pt.size = 0.1)
p3

#find positive markers for each cluster
everyone.marker <- FindAllMarkers(combined.harmony, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
everyone.marker %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)
write.csv(everyone.marker, file = "/data/scratch/RNA2021markers.csv")

#to rename idents based on the markers
combined.harmony <- RenameIdents(combined.harmony, '0' = "TAL", '1' = "EC1", '2' = "mPT", '3' = "fPT", '4' = "DCT2", '5' = "Stromal", '6' = "fPT", '7' = "DCT1", '8' = "mPT", '9' = "fPT", '10' = "PT5", '11' = "fPT", '12' = "mPT",  '13' = "Monocytes", '14' = "PT", '15' = "Stromal", '16' = "PT", '17' = "EC2", '18' = "Novel", '19' = "ICa", '20' = "mPT", '21' = "ICb", '22' = "Podocyte", '23' = "PT", '24' = "fPT", '25' = "CDPC", '26' = "Stromal")

#Toredorder clusters levels(combined.harmony)
levels(combined.harmony) <- c("Podocyte", "EC1", "EC2", "Stromal", "mPT","fPT","PT", "PT5", "TAL", "DCT1", "DCT2", "CDPC", "ICa","ICb", "Monocytes", "Novel")
combined.harmony$celltype <- Idents(combined.harmony)
DimPlot(combined.harmony,  label = TRUE) #772X500
DimPlot(combined.harmony, reduction = "umap", split.by = "groupid", label = TRUE) #772X500
DimPlot(combined.harmony, reduction = "umap", group.by = "sex", label = FALSE) #772X500

#create idents based on genotype or sex
Idents(combined.harmony) <- "celltype"
combined.harmony$celltype.groupid <- paste(Idents(combined.harmony), combined.harmony$groupid, sep = "_")
Idents(combined.harmony) <- "celltype.groupid"
combined.harmony$celltype.groupid.sex <- paste(Idents(combined.harmony), combined.harmony$sex, sep = "_")
Idents(combined.harmony) <- "celltype.groupid.sex"

#cell counts etc.
table(Idents(combined.harmony))
table(combined.harmony$groupid)
table(combined.harmony$sex)
table(Idents(combined.harmony), combined.harmony$groupid)
###Export normalized counts
counts.matrix <-as.matrix(combined.harmony [["RNA"]]@data)
write.csv(counts.matrix, file = "/data/scratch/RNA2021counts.csv")

#Calculate DEGS within a cluster between Con and KO using DESeq2.
#use write csv because write xlsx doesn't give gene name, only pct and pvalues. output are indivudal files for each cluster.

output <- "/data/scratch/DEG_"
Idents(combined.harmony) <- "celltype"
for (i in (Idents(combined.harmony)))({
  try({
    ident1 <- paste0(i,"_Con")
    ident2 <- paste0(i,"_KO")
    Idents(combined.harmony) <- "celltype.groupid"
    deg <- FindMarkers(combined.harmony, ident.1 = ident1, ident.2=ident2, test.use = "DESeq2", verbose = TRUE) 
    write.csv(deg, file=paste0(output,i,".csv"))
  })
})

#export all barcodes
RNA2021_harm <-saveRDS(combined.harmony, file = "/data/scratch/RNA2021_harm.rds")
Idents(RNA2021_harm) <- "celltype"
RNA_barcodes <-Idents(RNA2021_harm)
write.csv(RNA_barcodes, file = "/data/scratch/Final/RNAbarcodes.csv")
save.image("/scratch/RNA2021.RData")
