
library(Signac)#v1.2.0
library(Seurat)#v4.0.0
library(GenomeInfoDb)
library(EnsDb.Mmusculus.v79)
library(ggplot2)
library(patchwork)
set.seed(1234)                
library(magrittr)
library(readr)
library(Matrix)
library(tidyr)
library(dplyr)
library(GenomicRanges)
library(future)
options(future.globals.maxSize = 50000 * 1024^2)
library(cowplot)
library(harmony)
library(ChIPpeakAnno)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(tibble)
library(openxlsx)

#key to the samples
CONF - 67
CONM - 77
KOF - 66
KOM - 76

#Import files for each sample
peaks.mko <- read.table(
  file = "/data/user/hyndmank/Chromatin/CellRanger/khATAC2/76/outs/filtered_peak_bc_matrix/peaks.bed",
  col.names = c("chr", "start", "end"))
peaks.mcon <- read.table(
  file = "/data/user/hyndmank/Chromatin/CellRanger/khATAC2/77/outs/filtered_peak_bc_matrix/peaks.bed",
  col.names = c("chr", "start", "end"))

peaks.fko <- read.table(
  file = "/data/user/hyndmank/Chromatin/CellRanger/khATAC2/66/outs/filtered_peak_bc_matrix/peaks.bed",
  col.names = c("chr", "start", "end"))
peaks.fcon <- read.table(
  file = "/data/user/hyndmank/Chromatin/CellRanger/khATAC2/67/outs/filtered_peak_bc_matrix/peaks.bed",
  col.names = c("chr", "start", "end"))


# convert to genomic ranges
gr.mko <- makeGRangesFromDataFrame(peaks.mko)
gr.mcon <- makeGRangesFromDataFrame(peaks.mcon)
gr.fko <- makeGRangesFromDataFrame(peaks.fko)
gr.fcon <- makeGRangesFromDataFrame(peaks.fcon)

# Create a unified set of peaks to quantify in each dataset
combined.peaks <- reduce(x = c(gr.mko, gr.mcon,gr.fko, gr.fcon ))

# Filter out bad peaks based on length
peakwidths <- width(combined.peaks)
combined.peaks <- combined.peaks[peakwidths  < 10000 & peakwidths > 20]
combined.peaks

# load metadata
metadata.mko <- read.table(
  file = "/data/user/hyndmank/Chromatin/CellRanger/khATAC2/76/outs/singlecell.csv",
  stringsAsFactors = FALSE,
  sep = ",",
  header = TRUE,
  row.names = 1
)[-1, ]

metadata.mcon <- read.table(
  file = "/data/user/hyndmank/Chromatin/CellRanger/khATAC2/77/outs/singlecell.csv",
  stringsAsFactors = FALSE,
  sep = ",",
  header = TRUE,
  row.names = 1
)[-1, ]

metadata.fko <- read.table(
  file = "/data/user/hyndmank/Chromatin/CellRanger/khATAC2/66/outs/singlecell.csv",
  stringsAsFactors = FALSE,
  sep = ",",
  header = TRUE,
  row.names = 1
)[-1, ]

metadata.fcon <- read.table(
  file = "/data/user/hyndmank/Chromatin/CellRanger/khATAC2/67/outs/singlecell.csv",
  stringsAsFactors = FALSE,
  sep = ",",
  header = TRUE,
  row.names = 1
)[-1, ]

# perform an initial filtering of low count cells
metadata.mko <- metadata.mko[metadata.mko$passed_filters > 500, ]
metadata.mcon <- metadata.mcon[metadata.mcon$passed_filters > 500, ]
metadata.fko <- metadata.fko[metadata.fko$passed_filters > 500, ]
metadata.fcon <- metadata.fcon[metadata.fcon$passed_filters > 500, ]

# create fragment objects
frags.mko <- CreateFragmentObject(
  path = "/data/user/hyndmank/Chromatin/CellRanger/khATAC2/76/outs/fragments.tsv.gz",
  cells = rownames(metadata.mko))

frags.mcon <- CreateFragmentObject(
  path = "/data/user/hyndmank/Chromatin/CellRanger/khATAC2/77/outs/fragments.tsv.gz",
  cells = rownames(metadata.mcon))

frags.fko <- CreateFragmentObject(
  path = "/data/user/hyndmank/Chromatin/CellRanger/khATAC2/66/outs/fragments.tsv.gz",
  cells = rownames(metadata.fko))

frags.fcon <- CreateFragmentObject(
  path = "/data/user/hyndmank/Chromatin/CellRanger/khATAC2/67/outs/fragments.tsv.gz",
  cells = rownames(metadata.fcon))

########## Create count matrix
mko.counts <- FeatureMatrix(
  fragments = frags.mko,
  features = combined.peaks,
  cells = rownames(metadata.mko))

mcon.counts <- FeatureMatrix(
  fragments = frags.mcon,
  features = combined.peaks,
  cells = rownames(metadata.mcon))

fko.counts <- FeatureMatrix(
  fragments = frags.fko,
  features = combined.peaks,
  cells = rownames(metadata.fko))

fcon.counts <- FeatureMatrix(
  fragments = frags.fcon,
  features = combined.peaks,
  cells = rownames(metadata.fcon))

######## chromatin assays
chrom_assay_mko <- CreateChromatinAssay(
  counts = mko.counts,
  sep = c(":", "-"),
  fragments = frags.mko,
  min.cells = 10,
  min.features = 200)

mko <- CreateSeuratObject(
  counts = chrom_assay_mko,
  assay = "peaks",
  meta.data = metadata.mko)

chrom_assay_mcon <- CreateChromatinAssay(
  counts = mcon.counts,
  sep = c(":", "-"),
  fragments = frags.mcon,
  min.cells = 10,
  min.features = 200)

mcon <- CreateSeuratObject(
  counts = chrom_assay_mcon,
  assay = "peaks",
  meta.data = metadata.mcon)

chrom_assay_fcon <- CreateChromatinAssay(
  counts = fcon.counts,
  sep = c(":", "-"),
  fragments = frags.fcon,
  min.cells = 10,
  min.features = 200)

fcon <- CreateSeuratObject(
  counts = chrom_assay_fcon,
  assay = "peaks",
  meta.data = metadata.fcon)

chrom_assay_fko <- CreateChromatinAssay(
  counts = fko.counts,
  sep = c(":", "-"),
  fragments = frags.fko,
  min.cells = 10,
  min.features = 200)

fko <- CreateSeuratObject(
  counts = chrom_assay_fko,
  assay = "peaks",
  meta.data = metadata.fko)

# add information to identify dataset of origin
mko$groupid <- 'ko'
mcon$groupid <- 'con'
fko$groupid <- 'ko'
fcon$groupid <- 'con'
mko$sex <- 'male'
mcon$sex <- 'male'
fko$sex <- 'female'
fcon$sex <- 'female'

# extract gene annotations from EnsDb
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79)
seqlevelsStyle(annotations) <- 'UCSC'
genome(annotations) <- "mm10"
Annotation(mko) <- annotations
Annotation(mcon) <- annotations
Annotation(fko) <- annotations
Annotation(fcon) <- annotations

#QC checks
mko <- NucleosomeSignal(object = mko)
mko$nucleosome_group <- ifelse(mko$nucleosome_signal > 4, 'NS > 4', 'NS < 4')
FragmentHistogram(object = mko, group.by = 'nucleosome_group', region = 'chr1-1-10000000')
mko <- TSSEnrichment(mko, fast = FALSE)
mko$pct_reads_in_peaks <- mko$peak_region_fragments / mko$passed_filters * 100
mko$blacklist_ratio <- mko$blacklist_region_fragments / mko$peak_region_fragments
mko$high.tss <- ifelse(mko$TSS.enrichment > 2, 'High', 'Low')
TSSPlot(mko, group.by = 'high.tss') + NoLegend()

VlnPlot(
  object = mko,
  features = c('pct_reads_in_peaks', 'peak_region_fragments',
               'TSS.enrichment', 'blacklist_ratio', 'nucleosome_signal'),
  pt.size = 0.1,
  ncol = 5
)

mko <- subset(
  x = mko,
  subset = peak_region_fragments > 2500 &
    peak_region_fragments < 25000 &
    pct_reads_in_peaks > 15 &
    blacklist_ratio < 0.02 &
    nucleosome_signal < 4 &
    TSS.enrichment > 2)
mko

#An object of class Seurat 
#199266 features across 3537 samples within 1 assay 
#Active assay: peaks (199266 features, 0 variable features)

mcon <- NucleosomeSignal(object = mcon)
mcon$nucleosome_group <- ifelse(mcon$nucleosome_signal > 4, 'NS > 4', 'NS < 4')
FragmentHistogram(object = mcon, group.by = 'nucleosome_group', region = 'chr1-1-10000000')
mcon <- TSSEnrichment(mcon, fast = FALSE)
mcon$pct_reads_in_peaks <- mcon$peak_region_fragments / mcon$passed_filters * 100
mcon$blacklist_ratio <- mcon$blacklist_region_fragments / mcon$peak_region_fragments
mcon$high.tss <- ifelse(mcon$TSS.enrichment > 2, 'High', 'Low')
TSSPlot(mcon, group.by = 'high.tss') + NoLegend()

VlnPlot(
  object = mcon,
  features = c('pct_reads_in_peaks', 'peak_region_fragments',
               'TSS.enrichment', 'blacklist_ratio', 'nucleosome_signal'),
  pt.size = 0.1,
  ncol = 5
)

mcon <- subset(
  x = mcon,
  subset = peak_region_fragments > 2500 &
    peak_region_fragments < 25000 &
    pct_reads_in_peaks > 15 &
    blacklist_ratio < 0.02 &
    nucleosome_signal < 4 &
    TSS.enrichment > 2)
mcon

#An object of class Seurat 
#182502 features across 4432 samples within 1 assay 
#Active assay: peaks (182502 features, 0 variable features)

fko <- NucleosomeSignal(object = fko)
fko$nucleosome_group <- ifelse(fko$nucleosome_signal > 4, 'NS > 4', 'NS < 4')
FragmentHistogram(object = fko, group.by = 'nucleosome_group', region = 'chr1-1-10000000')
fko <- TSSEnrichment(fko, fast = FALSE)
fko$pct_reads_in_peaks <- fko$peak_region_fragments / fko$passed_filters * 100
fko$blacklist_ratio <- fko$blacklist_region_fragments / fko$peak_region_fragments
fko$high.tss <- ifelse(fko$TSS.enrichment > 2, 'High', 'Low')
TSSPlot(fko, group.by = 'high.tss') + NoLegend()

fko <- subset(
  x = fko,
  subset = peak_region_fragments > 2500 &
    peak_region_fragments < 25000 &
    pct_reads_in_peaks > 15 &
    blacklist_ratio < 0.02 &
    nucleosome_signal < 4 &
    TSS.enrichment > 2)
fko

#An object of class Seurat 
#196015 features across 4314 samples within 1 assay 
#Active assay: peaks (196015 features, 0 variable features)

fcon <- NucleosomeSignal(object = fcon)
fcon$nucleosome_group <- ifelse(fcon$nucleosome_signal > 4, 'NS > 4', 'NS < 4')
FragmentHistogram(object = fcon, group.by = 'nucleosome_group', region = 'chr1-1-10000000')
fcon <- TSSEnrichment(fcon, fast = FALSE)
fcon$pct_reads_in_peaks <- fcon$peak_region_fragments / fcon$passed_filters * 100
fcon$blacklist_ratio <- fcon$blacklist_region_fragments / fcon$peak_region_fragments
fcon$high.tss <- ifelse(fcon$TSS.enrichment > 2, 'High', 'Low')
TSSPlot(fcon, group.by = 'high.tss') + NoLegend()

fcon <- subset(
  x = fcon,
  subset = peak_region_fragments > 2500 &
    peak_region_fragments < 25000 &
    pct_reads_in_peaks > 15 &
    blacklist_ratio < 0.02 &
    nucleosome_signal < 4 &
    TSS.enrichment > 2)
fcon

#An object of class Seurat 
#178245 features across 3584 samples within 1 assay 
#Active assay: peaks (178245 features, 0 variable features)

# merge all datasets, adding a cell ID to make sure cell names are unique
combined <- merge(
  x = mko,
  y = list(mcon, fko, fcon))
combined[["peaks"]]

#ChromatinAssay data with 203918 features for 15868 cells
#Variable features: 0 
#Genome: mm10 
#Annotation present: TRUE 
#Motifs present: FALSE 
#Fragment files: 4 

#TF-IDF normalization, non-linear reduction and clustering
combined <- RunTFIDF(combined)
combined <- FindTopFeatures(combined, min.cutoff = 20)
combined <- RunSVD(combined)
combined <- RunUMAP(combined, dims = 2:50, reduction = 'lsi')


#Harmony batch correction
combined.harmony <- RunHarmony(
  object = combined,
  group.by.vars = 'groupid',
  reduction = 'lsi',
  assay.use = 'peaks',
  project.dim = FALSE)

# re-compute the UMAP using corrected LSI embeddings
combined.harmony <- RunUMAP(combined.harmony, dims = 2:30, reduction = 'harmony')

combined.harmony <- FindNeighbors(
  object = combined.harmony,
  reduction = 'harmony',
  dims = 2:30)

combined.harmony <- FindClusters(
  object = combined.harmony,
  algorithm = 3,
  resolution = 1.2,
  verbose = FALSE
)

p4 <- DimPlot(combined, group.by = 'groupid', pt.size = 0.1) + ggplot2::ggtitle("Unintegrated")
p5 <- DimPlot(combined.harmony, group.by = 'groupid', pt.size = 0.1) + ggplot2::ggtitle("Harmony integration")
p4 + p5 #1220x588

DimPlot(object = combined.harmony, label = TRUE) + NoLegend()

#add gene activities
gene.activities <- GeneActivity(combined.harmony)
combined.harmony[['RNA']] <- CreateAssayObject(counts = gene.activities)
combined.harmony <- NormalizeData(
  object = combined.harmony,
  assay = 'RNA',
  normalization.method = 'LogNormalize',
  scale.factor = median(combined.harmony$nCount_RNA))

#save combined.harmony as ATAC2021.rds
saveRDS(combined.harmony, file = "/data/scratch/hyndmank/Humphreys/test/ATAC2021.rds")
save.image("/scratch/hyndmank/Humphreys/test/ATAC.RData")

#Integrating with scRNA-seq data, Load the pre-processed scRNA-seq data
ATAC2021 <- readRDS("/data/scratch/hyndmank/Humphreys/test/ATAC2021.rds")
RNA2021_harm <- readRDS("/data/scratch/hyndmank/Humphreys/test/RNA2021_harm.rds")
DefaultAssay(ATAC2021) <- 'RNA'
DefaultAssay(RNA2021_harm) <- 'RNA'

transfer.anchors <- FindTransferAnchors(
  reference = RNA2021_harm,
  query = ATAC2021,
  reduction = 'cca')

predicted.labels <- TransferData(
  anchorset = transfer.anchors,
  refdata = RNA2021_harm$celltype,
  weight.reduction = ATAC2021[['lsi']],
  dims = 2:30)

ATAC2021 <- AddMetaData(object = ATAC2021, metadata = predicted.labels)

p6 <- DimPlot(
  object = RNA2021_harm,
  group.by = 'celltype',
  label = TRUE,
  repel = TRUE) + NoLegend() + ggtitle('scRNA-seq')

p7 <- DimPlot(
  object = ATAC2021,
  group.by = 'predicted.id',
  label = TRUE,
  repel = TRUE) + NoLegend() + ggtitle('scATAC-seq')

figS3a <-p6 + p7 #1220x588
figS3a

#FindAllMarkers for positive markers only in each cluster, must be in assay 'RNA'
combined.markers <- FindAllMarkers(ATAC2021, only.pos = TRUE, min.pct = 0.5, logfc.threshold = 0.58)
write.csv(combined.markers, "/data/scratch/hyndmank/Humphreys/test/khATAC2mergeb.markers.csv", quote = F)

#rename the clusters
ATAC2021 <- RenameIdents(ATAC2021, '0' = "mTAL", '1' = "mPTS1", '2' = "fPTS1", '3' = "cTAL", '4' = "EC1", '5' = "DCT1", '6' = "mPTS2", '7' = "fPTS2", '8' = "mPTS2", '9' = "fPTS2", '10' = "DCT2", '11' = "fPTS1", '12' = "Stromal",  '13' = "PT", '14' = "CDPC", '15' = "fPTS3", '16' = "mPTS3", '17' = "EC2", '18' = "mPTS3", '19' = "IC", '20' = "fPTS3", '21' = "PT5", '22' = "Monocytes", '23' = "Novel",'24' = "Podocyte", '25' = "Stromal", '26' = "Monocytes")
ATAC2021$celltype <- Idents(ATAC2021)
DimPlot(object = ATAC2021, label = TRUE, repel = TRUE)

#Toredorder clusters levels(ATAC2021)
levels(ATAC2021) <- c("Podocyte", "EC1", "EC2", "Stromal", "mPTS1","fPTS1","mPTS2","fPTS2", "mPTS3", "fPTS3", "PT", "PT5", "mTAL","cTAL", "DCT1", "DCT2", "CDPC", "IC","Monocytes", "Novel")

#export barcode ids
Idents(ATAC2021) <- "celltype"
ATAC_barcodes <-Idents(ATAC2021)
write.csv(ATAC_barcodes, file = "/data/scratch/hyndmank/Humphreys/test/Final/ATACbarcodes.csv")

#Dimplots for snRNA and snATAC
p8 <- DimPlot(object = RNA2021_harm, label = FALSE, pt.size = 0.1) + NoLegend() + ggplot2::ggtitle("snRNA-seq")
p9 <- DimPlot(object = ATAC2021, label = FALSE, pt.size = 0.1)+ NoLegend() + ggplot2::ggtitle("snATAC-seq")
fig2ab <-p8 + p9
fig2ab#1220x588

p10 <- DimPlot(object = RNA2021_harm, label = FALSE, pt.size = 0.1, group.by = "groupid") + NoLegend() + ggplot2::ggtitle("snRNA-seq")
p11 <- DimPlot(object = ATAC2021, label = FALSE, pt.size = 0.1, group.by = "groupid")+ NoLegend() + ggplot2::ggtitle("snATAC-seq")
fig4ab <-p10 + p11
fig4ab#1220x588

p12 <- DimPlot(object = RNA2021_harm, label = FALSE, pt.size = 0.1, group.by = "sex") + NoLegend() + ggplot2::ggtitle("snRNA-seq")
p13 <- DimPlot(object = ATAC2021, label = FALSE, pt.size = 0.1, group.by = "sex")+ NoLegend() + ggplot2::ggtitle("snATAC-seq")
figS2ab <-p12 + p13
figS2ab#1220x588


#to change idents to separate con and KO. Original Clusters are in meta$seurat_clusters.  Then we can separate based on sex.
ATAC2021$celltype.groupid <- paste(Idents(ATAC2021), ATAC2021$groupid, sep = "_")
Idents(ATAC2021) <- "celltype.groupid"
ATAC2021$celltype.groupid.sex <- paste(Idents(ATAC2021), ATAC2021$sex, sep = "_")
Idents(ATAC2021) <- "celltype"
Idents(ATAC2021) <- "celltype.groupid.sex"
Idents(ATAC2021) <-"seurat_clusters"

#cell counts etc.
table(Idents(ATAC2021))
table(ATAC2021$groupid)
table(ATAC2021$sex)
table(Idents(ATAC2021), ATAC2021$groupid)
table(Idents(ATAC2021), ATAC2021$sex)
prop.table(table(Idents(ATAC2021)))
prop.table(table(Idents(ATAC2021), ATAC2021$groupid), margin = 2)

#dotplot of cluster markers
DefaultAssay(ATAC2021) <- 'RNA'
DefaultAssay(RNA2021_harm) <- 'RNA'
markers.to.plot <- c("Bcat1", "Cd74", "Atp6v1g3", "Aqp2","Scnn1b","Calb1","Slc12a3", "Slc12a1","Plcl2", "Slco1a6","Slc5a1","Slc5a2", "Lrp2","Pdgfrb", "Pecam1","Nphs2")
fig2c <-DotPlot(RNA2021_harm, features = rev(markers.to.plot), cols = c("lightblue", "purple"), scale.min =0, dot.scale = 6) + RotatedAxis() 
fig2c #685x434
fig2d <-DotPlot(ATAC2021, features = rev(markers.to.plot), cols = c("lightblue", "purple"), scale.min =0, dot.scale = 6) + RotatedAxis() 
fig2d #685x434

#Differentially accessible chromatin within each cluster between genotypes.
DefaultAssay(ATAC2021) <- 'peaks'
Idents(ATAC2021) <- "celltype.groupid"

DiffGenesPod <- FindMarkers(ATAC2021, ident.1 = "Podocyte_con", ident.2 = "Podocyte_ko", min.pct = 0.2, test.use = 'LR', logfc.threshold = 0.25, latent.vars = 'peak_region_fragments')
open <- rownames(DiffGenesPod)  
cf <- ClosestFeature(ATAC2021, regions = open)
Pod <-(cbind(DiffGenesPod, gene=cf$gene_name, distance=cf$distance))
write.xlsx(Pod, file = "/data/scratch/hyndmank/Humphreys/test/Genotype_DAC/Pod.xlsx")

DiffGenesCDPC <- FindMarkers(ATAC2021, ident.1 = "CDPC_con", ident.2 = "CDPC_ko", min.pct = 0.2, test.use = 'LR', logfc.threshold = 0.25, latent.vars = 'peak_region_fragments')
open <- rownames(DiffGenesCDPC)  
cf <- ClosestFeature(ATAC2021, regions = open)
CDPC <-(cbind(DiffGenesCDPC, gene=cf$gene_name, distance=cf$distance))
write.xlsx(CDPC, file = "/data/scratch/hyndmank/Humphreys/test/Genotype_DAC/CDPC.xlsx")

DiffGenesmPTS1 <- FindMarkers(ATAC2021, ident.1 = "mPTS1_con", ident.2 = "mPTS1_ko", min.pct = 0.2, test.use = 'LR', logfc.threshold = 0.25, latent.vars = 'peak_region_fragments')
open <- rownames(DiffGenesmPTS1)  
cf <- ClosestFeature(ATAC2021, regions = open)
mPTS1 <-(cbind(DiffGenesmPTS1, gene=cf$gene_name, distance=cf$distance))
write.xlsx(mPTS1, file = "/data/scratch/hyndmank/Humphreys/test/Genotype_DAC/mPTS1.xlsx")

DiffGenesmPTS2 <- FindMarkers(ATAC2021, ident.1 = "mPTS2_con", ident.2 = "mPTS2_ko", min.pct = 0.2, test.use = 'LR', logfc.threshold = 0.25, latent.vars = 'peak_region_fragments')
open <- rownames(DiffGenesmPTS2)  
cf <- ClosestFeature(ATAC2021, regions = open)
mPTS2 <-(cbind(DiffGenesmPTS2, gene=cf$gene_name, distance=cf$distance))
write.xlsx(mPTS2, file = "/data/scratch/hyndmank/Humphreys/test/Genotype_DAC/mPTS2.xlsx")

DiffGenesmPTS3 <- FindMarkers(ATAC2021, ident.1 = "mPTS3_con", ident.2 = "mPTS3_ko", min.pct = 0.2, test.use = 'LR', logfc.threshold = 0.25, latent.vars = 'peak_region_fragments')
open <- rownames(DiffGenesmPTS3)  
cf <- ClosestFeature(ATAC2021, regions = open)
mPTS3 <-(cbind(DiffGenesmPTS3, gene=cf$gene_name, distance=cf$distance))
write.xlsx(mPTS3, file = "/data/scratch/hyndmank/Humphreys/test/Genotype_DAC/mPTS3.xlsx")

DiffGenesfPTS1 <- FindMarkers(ATAC2021, ident.1 = "fPTS1_con", ident.2 = "fPTS1_ko", min.pct = 0.2, test.use = 'LR', logfc.threshold = 0.25, latent.vars = 'peak_region_fragments')
open <- rownames(DiffGenesfPTS1)  
cf <- ClosestFeature(ATAC2021, regions = open)
fPTS1 <-(cbind(DiffGenesfPTS1, gene=cf$gene_name, distance=cf$distance))
write.xlsx(fPTS1, file = "/data/scratch/hyndmank/Humphreys/test/Genotype_DAC/fPTS1.xlsx")

DiffGenesfPTS2 <- FindMarkers(ATAC2021, ident.1 = "fPTS2_con", ident.2 = "fPTS2_ko", min.pct = 0.2, test.use = 'LR', logfc.threshold = 0.25, latent.vars = 'peak_region_fragments')
open <- rownames(DiffGenesfPTS2)  
cf <- ClosestFeature(ATAC2021, regions = open)
fPTS2 <-(cbind(DiffGenesfPTS2, gene=cf$gene_name, distance=cf$distance))
write.xlsx(fPTS2, file = "/data/scratch/hyndmank/Humphreys/test/Genotype_DAC/fPTS2.xlsx")

DiffGenesfPTS3 <- FindMarkers(ATAC2021, ident.1 = "fPTS3_con", ident.2 = "fPTS3_ko", min.pct = 0.2, test.use = 'LR', logfc.threshold = 0.25, latent.vars = 'peak_region_fragments')
open <- rownames(DiffGenesfPTS3)  
cf <- ClosestFeature(ATAC2021, regions = open)
fPTS3 <-(cbind(DiffGenesfPTS3, gene=cf$gene_name, distance=cf$distance))
write.xlsx(fPTS3, file = "/data/scratch/hyndmank/Humphreys/test/Genotype_DAC/fPTS3.xlsx")

DiffGenesPT5 <- FindMarkers(ATAC2021, ident.1 = "PT5_con", ident.2 = "PT5_ko", min.pct = 0.2, test.use = 'LR', logfc.threshold = 0.25, latent.vars = 'peak_region_fragments')
open <- rownames(DiffGenesPT5)  
cf <- ClosestFeature(ATAC2021, regions = open)
PT5 <-(cbind(DiffGenesPT5, gene=cf$gene_name, distance=cf$distance))
write.xlsx(PT5, file = "/data/scratch/hyndmank/Humphreys/test/Genotype_DAC/PT5.xlsx")

DiffGenesEC1 <- FindMarkers(ATAC2021, ident.1 = "EC1_con", ident.2 = "EC1_ko", min.pct = 0.2, test.use = 'LR', logfc.threshold = 0.25, latent.vars = 'peak_region_fragments')
open <- rownames(DiffGenesEC1)  
cf <- ClosestFeature(ATAC2021, regions = open)
EC1 <-(cbind(DiffGenesEC1, gene=cf$gene_name, distance=cf$distance))
write.xlsx(EC1, file = "/data/scratch/hyndmank/Humphreys/test/Genotype_DAC/EC1.xlsx")

DiffGenesEC2 <- FindMarkers(ATAC2021, ident.1 = "EC2_con", ident.2 = "EC2_ko", min.pct = 0.2, test.use = 'LR', logfc.threshold = 0.25, latent.vars = 'peak_region_fragments')
open <- rownames(DiffGenesEC2)  
cf <- ClosestFeature(ATAC2021, regions = open)
EC2 <-(cbind(DiffGenesEC2, gene=cf$gene_name, distance=cf$distance))
write.xlsx(EC2, file = "/data/scratch/hyndmank/Humphreys/test/Genotype_DAC/EC2.xlsx")

DiffGenesDCT1 <- FindMarkers(ATAC2021, ident.1 = "DCT1_con", ident.2 = "DCT1_ko", min.pct = 0.2, test.use = 'LR', logfc.threshold = 0.25, latent.vars = 'peak_region_fragments')
open <- rownames(DiffGenesDCT1)  
cf <- ClosestFeature(ATAC2021, regions = open)
DCT1 <-(cbind(DiffGenesDCT1, gene=cf$gene_name, distance=cf$distance))
write.xlsx(DCT1, file = "/data/scratch/hyndmank/Humphreys/test/Genotype_DAC/DCT1.xlsx")

DiffGenesDCT2 <- FindMarkers(ATAC2021, ident.1 = "DCT2_con", ident.2 = "DCT2_ko", min.pct = 0.2, test.use = 'LR', logfc.threshold = 0.25, latent.vars = 'peak_region_fragments')
open <- rownames(DiffGenesDCT2)  
cf <- ClosestFeature(ATAC2021, regions = open)
DCT2 <-(cbind(DiffGenesDCT2, gene=cf$gene_name, distance=cf$distance))
write.xlsx(DCT2, file = "/data/scratch/hyndmank/Humphreys/test/Genotype_DAC/DCT2.xlsx")

DiffGenescTAL <- FindMarkers(ATAC2021, ident.1 = "cTAL_con", ident.2 = "cTAL_ko", min.pct = 0.2, test.use = 'LR', logfc.threshold = 0.25, latent.vars = 'peak_region_fragments')
open <- rownames(DiffGenescTAL)  
cf <- ClosestFeature(ATAC2021, regions = open)
cTAL <-(cbind(DiffGenescTAL, gene=cf$gene_name, distance=cf$distance))
write.xlsx(cTAL, file = "/data/scratch/hyndmank/Humphreys/test/Genotype_DAC/cTAL.xlsx")

DiffGenesmTAL <- FindMarkers(ATAC2021, ident.1 = "mTAL_con", ident.2 = "mTAL_ko", min.pct = 0.2, test.use = 'LR', logfc.threshold = 0.25, latent.vars = 'peak_region_fragments')
open <- rownames(DiffGenesmTAL)  
cf <- ClosestFeature(ATAC2021, regions = open)
mTAL <-(cbind(DiffGenesmTAL, gene=cf$gene_name, distance=cf$distance))
write.xlsx(mTAL, file = "/data/scratch/hyndmank/Humphreys/test/Genotype_DAC/mTAL.xlsx")

DiffGenesIC <- FindMarkers(ATAC2021, ident.1 = "IC_con", ident.2 = "IC_ko", min.pct = 0.2, test.use = 'LR', logfc.threshold = 0.25, latent.vars = 'peak_region_fragments')
open <- rownames(DiffGenesIC)  
cf <- ClosestFeature(ATAC2021, regions = open)
IC <-(cbind(DiffGenesIC, gene=cf$gene_name, distance=cf$distance))
write.xlsx(IC, file = "/data/scratch/hyndmank/Humphreys/test/Genotype_DAC/IC.xlsx")

DiffGenesstromal <- FindMarkers(ATAC2021, ident.1 = "Stromal_con", ident.2 = "Stromal_ko", min.pct = 0.2, test.use = 'LR', logfc.threshold = 0.25, latent.vars = 'peak_region_fragments')
open <- rownames(DiffGenesstromal)  
cf <- ClosestFeature(ATAC2021, regions = open)
Stromal <-(cbind(DiffGenesstromal, gene=cf$gene_name, distance=cf$distance))
write.xlsx(Stromal, file = "/data/scratch/hyndmank/Humphreys/test/Genotype_DAC/Stromal.xlsx")

DiffGenesNovel <- FindMarkers(ATAC2021, ident.1 = "Novel_con", ident.2 = "Novel_ko", min.pct = 0.2, test.use = 'LR', logfc.threshold = 0.25, latent.vars = 'peak_region_fragments')
open <- rownames(DiffGenesNovel)  
cf <- ClosestFeature(ATAC2021, regions = open)
Novel <-(cbind(DiffGenesNovel, gene=cf$gene_name, distance=cf$distance))
write.xlsx(Novel, file = "/data/scratch/hyndmank/Humphreys/test/Genotype_DAC/Novel.xlsx")

DiffGenesMonocytes <- FindMarkers(ATAC2021, ident.1 = "Monocytes_con", ident.2 = "Monocytes_ko", min.pct = 0.2, test.use = 'LR', logfc.threshold = 0.25, latent.vars = 'peak_region_fragments')
open <- rownames(DiffGenesMonocytes)  
cf <- ClosestFeature(ATAC2021, regions = open)
Monocytes <-(cbind(DiffGenesMonocytes, gene=cf$gene_name, distance=cf$distance))
write.xlsx(Monocytes, file = "/data/scratch/hyndmank/Humphreys/test/Genotype_DAC/Monocytes.xlsx")

saveRDS(ATAC2021, file = "/data/scratch/hyndmank/Humphreys/test/Final/ATAC2021.rds")
saveRDS(RNA2021_harm, file = "/data/scratch/hyndmank/Humphreys/test/Final/RNA2021_harm.rds")
#combine all xlsx into one file
#now moved on to DACvs DEG.R code
