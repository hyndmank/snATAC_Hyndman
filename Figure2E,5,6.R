library(Signac)
library(Seurat)
library(ggplot2)
library(patchwork)
library(magrittr)
library(readr)
library(Matrix)
library(tidyr)
library(dplyr)
library(cowplot)
library(ChIPpeakAnno)
library(tibble)
library(ChIPseeker)
library(hrbrthemes)
library(clusterProfiler)
library(tidyverse)
library(ggpubr)
library(ggrepel)
library(ggfittext)

RNA2021_harm <- readRDS("/data/scratch/Final/RNA2021_harm.rds")
ATAC2021 <- readRDS("/data/scratch/Final/ATAC2021.rds")
DefaultAssay(ATAC2021) <- 'peaks'
Idents(ATAC2021) <- "celltype"
DefaultAssay(RNA2021_harm) <-"RNA"
Idents(RNA2021_harm) <- "celltype"

#coveragebrowser 
DefaultAssay(ATAC2021) <- 'peaks'
p.list <-CoverageBrowser(ATAC2021, region = "Aqp4")

#coverage plots in figure 2E. Exported and arranged in powerpoint.
CoveragePlot(
  object = ATAC2021,
  region = "chr1-156310727-156328035",
  features = "Nphs2",
  annotation = FALSE,
  peaks = FALSE,
  tile = FALSE,
  links = FALSE
)
CoveragePlot(
  object = ATAC2021,
  region = "chr11-106654217-106750628",
  features = "Pecam1",
  annotation = FALSE,
  peaks = TRUE,
  tile = FALSE,
  links = FALSE
)
CoveragePlot(
  object = ATAC2021,
  region = "chr2-173153048-173159273",
  features = "Pck1",
  annotation = FALSE,
  peaks = FALSE,
  tile = FALSE,
  links = FALSE
)
CoveragePlot(
  object = ATAC2021,
  region = "chr7-128265657-128272430",
  features = "Slc5a2",
  annotation = FALSE,
  peaks = FALSE,
  tile = FALSE,
  links = FALSE
)
CoveragePlot(
  object = ATAC2021,
  region = "chr2-125152505-125230002",
  features = "Slc12a1",
  annotation = FALSE,
  peaks = FALSE,
  tile = FALSE,
  links = FALSE
)
CoveragePlot(
  object = ATAC2021,
  region = "chr8-94329192-94366213",
  features = "Slc12a3",
  annotation = FALSE,
  peaks = FALSE,
  tile = FALSE,
  links = FALSE
)
CoveragePlot(
  object = ATAC2021,
  region = "chr15-99577550-99585545",
  features = "Aqp2",
  annotation = FALSE,
  peaks = FALSE,
  tile = FALSE,
  links = FALSE
)

#For figures 5 and 6 to compared DAC and DEG we need cluster names to be the same, so we combined PTs and IC to match the DAC and DEG cluster ids.
ATAC2021 <- RenameIdents(ATAC2021, `mTAL` = "TAL", `mPTS1` = "mPT", 
                         `fPTS1` = "fPT", `cTAL` = "TAL", 
                         'mPTS2' = "mPT", `fPTS2` = "fPT", `mPTS2` = "mPT", 
                         `fPTS2` = "fPT", 
                         'fPTS3' = "fPT", 'mPTS3' = "mPT", 
                         'mPTS3' = "mPT", 'fPTS3' = "fPT")
ATAC2021$PTgroup <- paste(Idents(ATAC2021))
levels(ATAC2021) <- c("Podocyte", "EC1", "EC2", "Stromal", "mPT","fPT","PT", "PT5", "TAL", 
                      "DCT1", "DCT2", "CDPC", "IC","Monocytes", "Novel")
Idents(ATAC2021) <- "PTgroup"
RNA2021_harm <-RenameIdents(RNA2021_harm, 'ICa' = "IC", "ICb" = "IC")
levels(RNA2021_harm) <- c("Podocyte", "EC1", "EC2", "Stromal", "mPT","fPT","PT", "PT5", "TAL", 
                      "DCT1", "DCT2", "CDPC", "IC","Monocytes", "Novel")

#Figure 5 Coverage plots, snATAC and snRNA expression plots.  Figure 5D, H, L are qPCR results plotted with prism.
#gene plot and coverage plot
#Aqp1 Figure 5A-C

coverage_plot <-CoveragePlot(
  object = ATAC2021,
  region = "chr6-55336432-55348555",
  group.by = "groupid",
  annotation = TRUE,
  peaks = TRUE,
  tile = FALSE,
  links = FALSE
)
coverage_plot #625x484

#RNA expression plots 200x300
DefaultAssay(ATAC2021) <- 'RNA'
VlnPlot <-VlnPlot(RNA2021_harm, features = "Aqp1", split.by = "groupid", idents = c("mPT", "fPT", "PT"))
VlnPlot
VlnPlot + theme(axis.title.y = element_blank(), axis.title = element_blank(), 
                plot.title = element_blank(), legend.position = "None", 
                axis.text.x = element_text(angle = 90), axis.text.y = element_text(angle = 90))

VlnPlot <-VlnPlot(ATAC2021, features = "Aqp1", split.by = "groupid", idents = c("mPT", "fPT", "PT"))
VlnPlot
VlnPlot + theme(axis.title.y = element_blank(), axis.title = element_blank(), 
                plot.title = element_blank(), legend.position = "None", 
                axis.text.x = element_text(angle = 90), axis.text.y = element_text(angle = 90))

#Aqp2 Figure 5E-G
DefaultAssay(ATAC2021) <- 'peaks'
coverage_plot <-CoveragePlot(
  object = ATAC2021,
  region = "chr15-99577550-99585545",
  group.by = "groupid",
  annotation = TRUE,
  peaks = TRUE,
  tile = FALSE,
  links = FALSE
)
coverage_plot 
VlnPlot <-VlnPlot(RNA2021_harm, features = "Aqp2", split.by = "groupid", idents = "CDPC")
VlnPlot
VlnPlot + theme(axis.title.y = element_blank(), axis.title = element_blank(), 
                plot.title = element_blank(), legend.position = "None", 
                axis.text.x = element_text(angle = 90), axis.text.y = element_text(angle = 90))
DefaultAssay(ATAC2021) <- 'RNA'
VlnPlot <-VlnPlot(ATAC2021, features = "Aqp2", split.by = "groupid", idents = "CDPC")
VlnPlot
VlnPlot + theme(axis.title.y = element_blank(), axis.title = element_blank(), 
                plot.title = element_blank(), legend.position = "None", 
                axis.text.x = element_text(angle = 90), axis.text.y = element_text(angle = 90))

#Aqp4 Figure 5I-K
DefaultAssay(ATAC2021) <- 'peaks'
coverage_plot <-CoveragePlot(
  object = ATAC2021,
  region = "chr18-15389394-15404684",
  group.by = "groupid",
  annotation = TRUE,
  peaks = TRUE,
  tile = FALSE,
  links = FALSE
)
coverage_plot 

VlnPlot <-VlnPlot(RNA2021_harm, features = "Aqp4", split.by = "groupid", idents = "CDPC")
VlnPlot
VlnPlot + theme(axis.title.y = element_blank(), axis.title = element_blank(), 
                plot.title = element_blank(), legend.position = "None", 
                axis.text.x = element_text(angle = 90), axis.text.y = element_text(angle = 90))
DefaultAssay(ATAC2021) <- 'RNA'
VlnPlot <-VlnPlot(ATAC2021, features = "Aqp4", split.by = "groupid", idents = "CDPC")
VlnPlot
VlnPlot + theme(axis.title.y = element_blank(), axis.title = element_blank(), 
                plot.title = element_blank(), legend.position = "None", 
                axis.text.x = element_text(angle = 90), axis.text.y = element_text(angle = 90))

#Abi3bp Figure 6A-C, 1452x 988
DefaultAssay(ATAC2021) <- 'peaks'
CoveragePlot(
  object = ATAC2021,
  region = "chr16-56475846-56690128",
  group.by = "groupid",
  annotation = TRUE,
  peaks = TRUE,
  tile = FALSE,
  links = FALSE
)

gene_plot <- AnnotationPlot(
  object = ATAC2021,
  region = "chr16-56475846-56690128")
gene_plot 

VlnPlot <-VlnPlot(RNA2021_harm, features = "Abi3bp", split.by = "groupid",idents = c("Stromal", "mPT", "fPT", "PT", "PT5"))
VlnPlot
VlnPlot + theme(axis.title.y = element_blank(), axis.title = element_blank(), 
                plot.title = element_blank(), legend.position = "None", 
                axis.text.x = element_text(angle = 90), axis.text.y = element_text(angle = 90))

DefaultAssay(ATAC2021) <-'RNA'
VlnPlot <-VlnPlot(ATAC2021, features = "Abi3bp", split.by = "groupid", idents = c("Stromal", "mPT", "fPT", "PT", "PT5"))
VlnPlot
VlnPlot + theme(axis.title.y = element_blank(), axis.title = element_blank(), 
                plot.title = element_blank(), legend.position = "None", 
                axis.text.x = element_text(angle = 90), axis.text.y = element_text(angle = 90))

#Acsm2 1452x 988
DefaultAssay(ATAC2021) <- 'peaks'
CoveragePlot(
  object = ATAC2021,
  region = "chr7-119554469-119600694",
  group.by = "groupid",
  annotation = TRUE,
  peaks = TRUE,
  tile = FALSE,
  links = FALSE
)

gene_plot <- AnnotationPlot(
  object = ATAC2021,
  region = "chr7-119554469-119600694")
gene_plot 

VlnPlot <-VlnPlot(RNA2021_harm, features = "Acsm2", split.by = "groupid",idents = c("mPT", "fPT", "PT", "PT5"))
VlnPlot
VlnPlot + theme(axis.title.y = element_blank(), axis.title = element_blank(), 
                plot.title = element_blank(), legend.position = "None", 
                axis.text.x = element_text(angle = 90), axis.text.y = element_text(angle = 90))
DefaultAssay(ATAC2021) <-'RNA'
VlnPlot <-VlnPlot(ATAC2021, features = "Acsm2", split.by = "groupid", idents = c("mPT", "fPT", "PT", "PT5"))
VlnPlot
VlnPlot + theme(axis.title.y = element_blank(), axis.title = element_blank(), 
                plot.title = element_blank(), legend.position = "None", 
                axis.text.x = element_text(angle = 90), axis.text.y = element_text(angle = 90))
