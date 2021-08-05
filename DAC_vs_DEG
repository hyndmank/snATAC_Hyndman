library(Signac)
library(Seurat)
library(EnsDb.Mmusculus.v79)
library(openxlsx)
library(readxl)
library(dplyr)

#We will compare cluster markers differentially accessible chromatin (DAC) with differentially expressed genes (DEG)
ATAC2021 <- readRDS("/data/scratch/Final/ATAC2021.rds")
RNA2021_harm <- readRDS("/data/scratch/Final/RNA2021_harm.rds")
Idents(ATAC2021) <- "celltype"
Idents(RNA2021_harm) <- "celltype"
DefaultAssay(ATAC2021) <- 'peaks'
DefaultAssay(RNA2021_harm) <- 'RNA'

#subset the control and KO animals.
KO <-subset(ATAC2021, subset = groupid == "ko")
CON <-subset(ATAC2021, subset = groupid == "con")
DefaultAssay(CON) <- 'peaks'
DefaultAssay(KO) <- 'peaks'

#look for dacs within CON
GetMarkers <- function(cluster, seurat_aggregate) {
  print(paste0("Finding dac for: ",cluster))
  CON <- seurat_aggregate
  dac <- FindMarkers(CON, 
                     ident.1 = cluster,    
                     test.use = 'LR', 
                     latent.vars = "peak_region_fragments",
                     min.pct = 0.2)
  open <- rownames(dac)  
  cf <- ClosestFeature(CON, regions = open)
  return(cbind(dac, gene=cf$gene_name, distance=cf$distance))}

idents <- levels(CON@meta.data$celltype)
list.cluster.dac <- lapply(idents, function(x) {GetMarkers(x, seurat_aggregate = CON)})
write.xlsx(list.cluster.dac, file = "/data/scratch/Final/dac.con.xlsx", sheetName = idents, rowNames = T)

#look for dacs within KO
GetMarkers <- function(cluster, seurat_aggregate) {
  print(paste0("Finding dac for: ",cluster))
  KO <- seurat_aggregate
  dac <- FindMarkers(KO, 
                     ident.1 = cluster,    
                     test.use = 'LR', 
                     latent.vars = "peak_region_fragments",
                     min.pct = 0.2)
  open <- rownames(dac)  
  cf <- ClosestFeature(KO, regions = open)
  return(cbind(dac, gene=cf$gene_name, distance=cf$distance))}

idents <- levels(KO@meta.data$celltype)
list.cluster.dac <- lapply(idents, function(x) {GetMarkers(x, seurat_aggregate = KO)})
write.xlsx(list.cluster.dac, file = "/data/scratch/Final/dac.KO.xlsx", sheetName = idents, rowNames = T)

#DEG
KO <-subset(RNA2021_harm, subset = groupid == "KO")
CON <-subset(RNA2021_harm, subset = groupid == "Con")

# FindMarkers and write to an xlsx file with default parameters
GetMarkers <- function(cluster, seurat_aggregate) {
  print(paste0("Finding DEG for: ",cluster))
  CON <- seurat_aggregate
  deg <- FindMarkers(CON, 
                     ident.1 = cluster,    
                     min.pct = 0.2) 
  return(deg)}
idents <- levels(CON@meta.data$celltype)
list.cluster.deg <- lapply(idents, function(x) {GetMarkers(x, seurat_aggregate = CON)})
write.xlsx(list.cluster.deg, file = "/data/scratch/Final/deg.con.xlsx", sheetName = idents, rowNames = T)

GetMarkers <- function(cluster, seurat_aggregate) {
  print(paste0("Finding DEG for: ",cluster))
  KO <- seurat_aggregate
  deg <- FindMarkers(KO, 
                     ident.1 = cluster,    
                     min.pct = 0.2) 
  return(deg)}
idents <- levels(KO@meta.data$celltype)
list.cluster.deg <- lapply(idents, function(x) {GetMarkers(x, seurat_aggregate = KO)})
write.xlsx(list.cluster.deg, file = "/data/scratch/Final/deg.ko.xlsx", sheetName = idents, rowNames = T)

#kept only adjusted p<0.05, import excel where each tab is a different cluster.  make sure tab cluster names are the same in the dac and deg.
sheet = excel_sheets("/data/scratch/Final/dac.sig.con.xlsx")
df = lapply(setNames(sheet, sheet), function(x) read_excel("/data/scratch/Final/dac.sig.con.xlsx", sheet=x))
df = bind_rows(df, .id="Sheet")
write.xlsx(df, file = ("/data/scratch/Final/Con_DAC.xlsx"))

sheet = excel_sheets("/data/scratch/data/scratch/Final/deg.sig.con.xlsx")
df = lapply(setNames(sheet, sheet), function(x) read_excel("/data/scratch/Final/deg.sig.con.xlsx", sheet=x))
df = bind_rows(df, .id="Sheet")
write.xlsx(df, file = ("/data/scratch/Final/Con_DEG.xlsx"))

sheet = excel_sheets("/data/scratch/Final/dac.sig.KO.xlsx")
df = lapply(setNames(sheet, sheet), function(x) read_excel("/data/scratch/Final/dac.sig.KO.xlsx", sheet=x))
df = bind_rows(df, .id="Sheet")
write.xlsx(df, file = ("/data/scratch/Final/KO_DAC.xlsx"))

sheet = excel_sheets("/data/scratch/Final/deg.sig.ko.xlsx")
df = lapply(setNames(sheet, sheet), function(x) read_excel("/data/scratch/Final/deg.sig.ko.xlsx", sheet=x))
df = bind_rows(df, .id="Sheet")
write.xlsx(df, file = ("/data/scratch/Final/KO_DEG.xlsx"))

DACCON <- read_excel("/data/scratch/Final/Con_DAC.xlsx")
DACCONdf <- DACCON[order(DACCON$gene, -abs(DACCON$avg_log2FC) ), ]### sort first
DACCONdf <- DACCONdf %>% distinct(Sheet, gene, .keep_all = TRUE) #filter distinct genes in each celltype

DEGCON <- read_excel("/data/scratch/Final/Con_DEG.xlsx")
DEGCONdf <- DEGCON[order(DEGCON$gene, -abs(DEGCON$avg_log2FC) ), ]### sort first

#join the two dataframes based upon the Sheet (cluster) and gene.
CON <-left_join(DEGCON, DACCONdf, by=c("Sheet","gene"))
write.xlsx(CON, file = ("/data/scratch/Final/CON_DACDEG.xlsx")) 

DACKO <- read_excel("/data/scratch/Final/KO_DAC.xlsx")
DACKOdf <- DACKO[order(DACKO$gene, -abs(DACKO$avg_log2FC) ), ]### sort first
DACKOdf <- DACKOdf %>% distinct(Sheet, gene, .keep_all = TRUE) #filter distinct genes in each celltype

DEGKO <- read_excel("/data/scratch/Final/KO_DEG.xlsx")
DEGKOdf <- DEGKO[order(DEGKO$gene, -abs(DEGKO$avg_log2FC) ), ]### sort first

#join the two dataframes based upon the Sheet (cluster) and gene.
KO <-left_join(DEGKO, DACKOdf, by=c("Sheet","gene"))
write.xlsx(KO, file = ("/data/scratch/Final/KO_DACDEG.xlsx")) 

#downloaded files and removed rows that were not shared between the datasets. Used Prism to do correlations reported in Figs 4CD.

#ggplot2 the data if you would like
p <- ggplot(CON, aes(x=avg_log2FC.y, y=avg_log2FC.x, label=gene,xmin=-2.5, xmax = 6.5, ymin =-5, ymax=8))
p + geom_point()

quad_count_CON <- CON %>%
  # Count how many with each combination of X and Y being positive
  count(right = avg_log2FC.y > 0, top = avg_log2FC.x > 0) %>%
  # TRUE = 1, FALSE = 0, so these map the TRUE to +1 and FALSE to -1
  mutate(avg_log2FC.y = 2 * (right - 0.5), avg_log2FC.x = 2 * (top - 0.5))
quad_count_CON

quad_count_KO <- KO %>%
  # Count how many with each combination of X and Y being positive
  count(right = avg_log2FC.y > 0, top = avg_log2FC.x > 0) %>%
  # TRUE = 1, FALSE = 0, so these map the TRUE to +1 and FALSE to -1
  mutate(avg_log2FC.y = 2 * (right - 0.5), avg_log2FC.x = 2 * (top - 0.5))
quad_count_KO
