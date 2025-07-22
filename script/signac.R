# script/signac_excn.R

#  Libraries 

library(Seurat)
library(Signac)
library(GenomicRanges)
library(data.table)
library(BSgenome.Mmusculus.UCSC.mm10)

set.seed(1234)

#  File paths 
seurat_path <- "data/seurat_chr1to3_25k_multiome_new.rds"
granges_path <- "data/gr_ranges.rds"
output_file <- "output/signac/linked_peaks_ExcN.tsv"

#  Load Seurat object 
multiome <- readRDS(seurat_path)

#  Subset Excitatory Neurons (ExcN) 
ExcN <- multiome[, multiome$cell_type_manual == "ExcN"]
cat("   →", ncol(ExcN), "cells\n")
saveRDS(ExcN, "data/seurat_ExcN.rds")


#  Load gene coordinates 
geneCoords <- readRDS(granges_path)
geneCoordDt <- as.data.table(geneCoords)
geneCoordDt <- subset(geneCoordDt, seqnames %in% c(1:3))  # limit to chr1–3
geneCoordDt[, seqnames := paste0("chr", seqnames)]

geneCoords <- makeGRangesFromDataFrame(as.data.frame(geneCoordDt), keep.extra.columns = TRUE)
geneCoords <- geneCoords[geneCoords$gene_id %in% rownames(ExcN[["RNA"]])]
names(geneCoords) <- geneCoords$gene_id

#  Region stats (optional, but safe) 
#head(ExcN[["ATAC"]]@misc$region_stats)
head(seqlevels(ExcN[["ATAC"]]))


#  Annotate ATAC assay 
Annotation(ExcN[["ATAC"]]) <- geneCoords

#  Run LinkPeaks 
linked <- LinkPeaks(
  object = ExcN[,1:100], # this is the subset part to try the small amount
  peak.assay = "ATAC",
  expression.assay = "RNA",
  peak.slot = "counts",
  expression.slot = "counts",
  gene.coords = geneCoords,
  method = "pearson",
  distance = 5e5,
  min.distance = NULL,
  min.cells = 10,
  n_sample = 200,
  pvalue_cutoff = 0.05,
  score_cutoff = 0.05,
  gene.id = TRUE,
  verbose = TRUE
)

#  Save output 
dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)

df <- as.data.frame(linked[["ATAC"]]@links)
out_df <- df[, c("peak", "gene", "score")]

write.table(out_df, file = output_file, sep = "\t", quote = FALSE, row.names = FALSE)

