library(Signac) 
library(Seurat) 
library(GenomeInfoDb)
library(harmony)
library(EnsDb.Mmusculus.v79)
library(ggplot2)
library(here)
library(tibble)
library(dplyr)
library(openxlsx)
set.seed(1234)
library(ChIPseeker)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene
library(clusterProfiler)
library(org.Mm.eg.db)

#t5n insertion like Muto et al. Figure 2.  did for control and ko separately
#ATAC2021 separate controls from KOs
ATAC2021 <- readRDS("/data/scratch/Final/ATAC2021.rds")
KO <-subset(ATAC2021, subset = groupid == "ko")
CON <-subset(ATAC2021, subset = groupid == "con")
DefaultAssay(CON) <- 'peaks'
DefaultAssay(KO) <- 'peaks'
Idents(CON) <- "celltype"

#open the DAC file generated with the DAC_vs_DEG code.
dacfile <- ("/data/scratch/Final/dac.con.xlsx")
idents <- getSheetNames(dacfile)
list.dac <- lapply(idents, function(x) {
  df <- read.xlsx(dacfile, sheet = x, rowNames = T) %>%
    rownames_to_column(var = "coord") %>%
    dplyr::mutate(celltype = x) 
})

# identify all unique cell-type-specific peaks and filter for log2fc > 0
all_dac <- bind_rows(list.dac) %>%
  dplyr::filter(avg_log2FC > 0) %>%
  dplyr::select("coord") %>%
  dplyr::distinct()

dac_aver <- AverageExpression(CON, features = all_dac$coord, assays = "peaks")
figxa <- pheatmap::pheatmap(dac_aver[["peaks"]], scale = 'row', cluster_rows = F, cluster_cols = F,
                            show_rownames = FALSE) #520x630

all_dac.gr <- StringToGRanges(all_dac$coord, sep = c(":","-"))
list.dac.gr <- lapply(seq(list.dac), function(x) {
  df <- list.dac[[x]]
  gr <- StringToGRanges(df$coord, sep = c(":","-"))
  return(gr)
})
names(list.dac.gr) <- idents

list.peakAnno <- lapply(list.dac.gr, annotatePeak, TxDb = txdb,
                        tssRegion = c(-3000, 3000), verbose = FALSE)
all.peakAnno <- annotatePeak(all_dac.gr, TxDb = txdb,
                             tssRegion = c(-3000, 3000), verbose = FALSE)

fig3a <- plotAnnoPie(all.peakAnno) #total T5n insertion in the dataset
fig3b <- plotAnnoBar(list.peakAnno) #celltype-specific analysis
figS5a <- plotDistToTSS(list.peakAnno) #Distribution of transcription factor-binding loci relative to TSS

#to do statistics on the distribution export as xlsx
write.xlsx(list.peakAnno, file = ("/data/scratch/Final/Con_peakanno.xlsx"))

#same thing but with the KO files
Idents(KO) <- "celltype"
levels(KO) <- c("Podocyte", "EC1", "EC2", "Stromal", "mPTS1","fPTS1","mPTS2","fPTS2", "mPTS3", "fPTS3", "PT", "PT5", "mTAL","cTAL", "DCT1", "DCT2", "CDPC", "IC","Monocytes", "Novel")
dacfile <- ("/data/scratch/Final/dac.ko.xlsx")
idents <- getSheetNames(dacfile)
list.dac <- lapply(idents, function(x) {
  df <- read.xlsx(dacfile, sheet = x, rowNames = T) %>%
    rownames_to_column(var = "coord") %>%
    dplyr::mutate(celltype = x) 
})

all_dac <- bind_rows(list.dac) %>%
  dplyr::filter(avg_log2FC > 0) %>%
  dplyr::select("coord") %>%
  dplyr::distinct()

dac_aver <- AverageExpression(KO, features = all_dac$coord, assays = "peaks")
figxb <- pheatmap::pheatmap(dac_aver[["peaks"]], scale = 'row', cluster_rows = F, cluster_cols = F,
                            show_rownames = FALSE) #520x630

all_dac.gr <- StringToGRanges(all_dac$coord, sep = c(":","-"))
list.dac.gr <- lapply(seq(list.dac), function(x) {
  df <- list.dac[[x]]
  gr <- StringToGRanges(df$coord, sep = c(":","-"))
  return(gr)
})
names(list.dac.gr) <- idents
# annotate the list of GRanges dac for each cell type
list.peakAnno <- lapply(list.dac.gr, annotatePeak, TxDb = txdb,
                        tssRegion = c(-3000, 3000), verbose = FALSE)
all.peakAnno <- annotatePeak(all_dac.gr, TxDb = txdb,
                             tssRegion = c(-3000, 3000), verbose = FALSE)

fig3c <- plotAnnoPie(all.peakAnno) 
fig3d <- plotAnnoBar(list.peakAnno) 
figS5b <- plotDistToTSS(list.peakAnno)
fig3c
fig3d
figS5b
