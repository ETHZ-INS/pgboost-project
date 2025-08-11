#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Signac)
  library(Seurat)
  library(data.table)
  library(GenomicRanges)
  library(BiocParallel)
})

# ---- CLI Arguments ----
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4) {
  stop("Usage:\n  Rscript generate_candidate_links_batch.R <seurat_file> <gene_coords_file> <n_links_per_file> <celltypes.tsv>")
}

seurat_file <- args[1]
gene_coords_file <- args[2]
number_candidate_links_per_file <- as.numeric(args[3])
celltypes_tsv <- args[4]

# ---- Load Data ----
cat("Reading Seurat object...\n")
seurat <- readRDS(seurat_file)

cat("Reading gene coordinates...\n")
geneCoords <- readRDS(gene_coords_file)
geneCoordDt <- as.data.table(geneCoords)
geneCoordDt[, seqnames := paste0('chr', seqnames)]
geneCoords <- makeGRangesFromDataFrame(as.data.frame(geneCoordDt), keep.extra.columns = TRUE)
names(geneCoords) <- geneCoords$gene_id

# ---- Prepare ATAC peaks ----
atac_mtx <- GetAssayData(seurat, assay = "ATAC", slot = "counts")
atac_mtx <- 1 * (atac_mtx > 0)
atac_mtx <- atac_mtx[rowSums(atac_mtx) / ncol(atac_mtx) > 0.05, ]

peaks <- data.frame(peak = rownames(atac_mtx))
peaks$chr <- sapply(strsplit(peaks$peak, "[:-]"), `[`, 1)
peaks$start <- as.numeric(sapply(strsplit(peaks$peak, "[:-]"), `[`, 2))
peaks$end <- as.numeric(sapply(strsplit(peaks$peak, "[:-]"), `[`, 3))
peaks$center <- (peaks$start + peaks$end) / 2

# ---- Function to generate links per celltype ----
generate_links_for_celltype <- function(celltype) {
  cat("Processing:", celltype, "\n")
  
  # RNA filtering
  rna_mtx <- GetAssayData(seurat, assay = "RNA", slot = "counts")
  rna_mtx <- rna_mtx[rowSums(rna_mtx) / ncol(rna_mtx) > 0.05, ]
  
  # Gene info
  gene_info <- data.frame(
    gene = geneCoords$gene_id,
    chr = as.character(seqnames(geneCoords)),
    strand = as.character(strand(geneCoords)),
    tss = ifelse(strand(geneCoords) == "+", start(geneCoords), end(geneCoords))
  )
  gene_info <- gene_info[gene_info$gene %in% rownames(rna_mtx), ]
  gene_info$tss_minus_500kb <- gene_info$tss - 500000
  gene_info$tss_plus_500kb <- gene_info$tss + 500000
  
  # Generate candidate links
  candidate_links <- data.frame()
  for (i in seq_len(nrow(gene_info))) {
    gene <- gene_info$gene[i]
    chr <- gene_info$chr[i]
    tss_min <- gene_info$tss_minus_500kb[i]
    tss_max <- gene_info$tss_plus_500kb[i]
    
    hits <- peaks[peaks$chr == chr & peaks$center >= tss_min & peaks$center <= tss_max, ]
    if (nrow(hits) > 0) {
      df <- data.frame(gene = gene, peak = hits$peak)
      candidate_links <- rbind(candidate_links, df)
    }
  }
  
  cat("Found", nrow(candidate_links), "links for", celltype, "\n")
  
  # Chunking and writing
  candidate_links <- candidate_links[sample(nrow(candidate_links)), ]
  n_chunks <- ceiling(nrow(candidate_links) / number_candidate_links_per_file)
  
  output_dir <- file.path("output/feature_scores/SCENT/scent_candidate_links", celltype)
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  
  for (i in seq_len(n_chunks)) {
    start <- (i - 1) * number_candidate_links_per_file + 1
    end <- min(i * number_candidate_links_per_file, nrow(candidate_links))
    chunk <- candidate_links[start:end, ]
    outfile <- file.path(output_dir, sprintf("chunk_%s_%d.txt", celltype, i))
    write.table(chunk, outfile, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
    cat("Saved:", outfile, "\n")
  }
}

# ---- Load Celltypes ----
df <- read.delim(celltypes_tsv, stringsAsFactors = FALSE)
celltypes <- df$celltype_key

# ---- Run in parallel ----
param <- MulticoreParam(workers = min(4, length(celltypes)))  
bplapply(celltypes, generate_links_for_celltype, BPPARAM = param)

cat("All jobs complete.\n")
