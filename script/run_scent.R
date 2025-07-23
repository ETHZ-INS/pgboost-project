library(SCENT)
library(Signac)
library(Seurat)
library(Matrix)

set.seed(1234)

seurat_path <- "data/seurat_ExcN.rds"
link_file <- "output/scent_candidate_links/chunk1.txt"
output_file <- "output/scent/chunk1_output.tsv"
skip_bootstrap <- TRUE
max_bootstrap_iter <- 50000

seurat <- readRDS(seurat_path)

rna_mtx <- GetAssayData(seurat, assay = "RNA", slot = "counts")
atac_mtx <- GetAssayData(seurat, assay = "ATAC", slot = "counts")
atac_mtx <- 1 * (atac_mtx > 0)  # Binarize ATAC


seurat@meta.data$nUMI <- seurat@meta.data$nCount_RNA
seurat@meta.data$celltype <- "focal"

mito.genes <- grep("^MT-", rownames(rna_mtx), value = TRUE)

percent.mito <- if (length(mito.genes) > 0) {
  colSums(rna_mtx[mito.genes, , drop = FALSE]) / colSums(rna_mtx)
} else {
  rep(0, ncol(rna_mtx))
}
seurat <- AddMetaData(seurat, percent.mito, col.name = "percent.mito")


meta <- seurat@meta.data
meta$log_nUMI <- log1p(meta$nUMI)
meta$cell <- rownames(meta)
meta <- meta[, c("cell", "nUMI", "log_nUMI", "percent.mito", "celltype")]

links <- read.table(link_file, sep = "\t", header = FALSE)
colnames(links) <- c("gene", "peak")


scent_obj <- CreateSCENTObj(
  rna = rna_mtx,
  atac = atac_mtx,
  meta.data = meta,
  peak.info = links,
  covariates = c("log_nUMI", "percent.mito"),
  celltypes = "celltype"
)

scent_obj <- SCENT_algorithm(
  object = scent_obj,
  celltype = "focal",
  ncores = 1,
  boot = !skip_bootstrap,
  maxboot = if (!skip_bootstrap) max_bootstrap_iter else NULL
)

head(scent_obj@SCENT.result)

# just write part is this
write.table(scent_obj@SCENT.result,
            "output/scent/chunk1_output.tsv",
            sep = "\t", quote = FALSE, row.names = FALSE)



