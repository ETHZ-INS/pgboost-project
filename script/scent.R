library(Signac)
library(Seurat)
library(data.table)
library(GenomicRanges)


seurat <- readRDS("data/seurat_ExcN.rds")
seurat

genes <- rownames(seurat[["RNA"]])
genes

geneCoords <- readRDS("data/gr_ranges.rds")

geneCoordDt <- as.data.table(geneCoords)
geneCoordDt <- subset(geneCoordDt, seqnames %in% c(1:3))  # or full genome
geneCoordDt[, seqnames := paste0('chr', seqnames)]

geneCoords <- makeGRangesFromDataFrame(as.data.frame(geneCoordDt), keep.extra.columns = TRUE)
names(geneCoords) <- geneCoords$gene_id


print(geneCoords[1:3])

output_dir <- "output/scent_candidate_links"
number_candidate_links_per_file <- 70000

if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

rna_mtx <- GetAssayData(seurat, assay = "RNA", slot = "counts")
atac_mtx <- GetAssayData(seurat, assay = "ATAC", slot = "counts")
atac_mtx <- 1 * (atac_mtx > 0)  # Binarize


rna_mtx <- rna_mtx[rowSums(rna_mtx) / ncol(rna_mtx) > 0.05, ]
atac_mtx <- atac_mtx[rowSums(atac_mtx) / ncol(atac_mtx) > 0.05, ]


gene_info <- data.frame(
  gene = geneCoords$gene_id,
  chr = as.character(seqnames(geneCoords)),
  strand = as.character(strand(geneCoords)),
  tss = ifelse(strand(geneCoords) == "+", start(geneCoords), end(geneCoords))
)

gene_info <- gene_info[gene_info$gene %in% rownames(rna_mtx), ]
gene_info$tss_minus_500kb <- gene_info$tss - 500000
gene_info$tss_plus_500kb <- gene_info$tss + 500000

peaks <- data.frame(peak = rownames(atac_mtx))
peaks$chr <- sapply(strsplit(peaks$peak, "[:-]"), `[`, 1)
peaks$start <- as.numeric(sapply(strsplit(peaks$peak, "[:-]"), `[`, 2))
peaks$end <- as.numeric(sapply(strsplit(peaks$peak, "[:-]"), `[`, 3))
peaks$center <- (peaks$start + peaks$end) / 2



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

candidate_links <- candidate_links[sample(nrow(candidate_links)), ]
n_chunks <- ceiling(nrow(candidate_links) / number_candidate_links_per_file)


for (i in 1:n_chunks) {
  start <- (i - 1) * number_candidate_links_per_file + 1
  end <- min(i * number_candidate_links_per_file, nrow(candidate_links))
  chunk <- candidate_links[start:end, ]
  outfile <- file.path(output_dir, sprintf("chunk%d.txt", i))
  write.table(chunk, outfile, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
}




