suppressPackageStartupMessages({
  library(Signac)
  library(Seurat)
  library(monocle3)
  library(monocle)
  library(cicero)
  library(scuttle)
  library(Matrix)
  library(BSgenome.Mmusculus.UCSC.mm10)
  library(dplyr)
  library(BiocParallel)
})

# ---- Function: run Cicero for one cell type ----
run_cicero_on_celltype <- function(
    seurat_file,
    celltype_key = "glia",
    cellsubtype_key = "astrocytes",
    leiden_col = "leiden",
    genome_size_file = "data/mm10.chrom.sizes"
) {
  output_dir <- file.path("output/feature_scores/cicero")
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  
  cat("Reading Seurat object...\n")
  data <- readRDS(seurat_file)
  data$cell_type_manual <- as.character(data$cell_type_manual)
  data$celltype <- as.character(data$celltype)
  
  # Subset to relevant cells
  if (is.null(cellsubtype_key)) {
    subset <- data[, data$cell_type_manual == celltype_key]
    cat(paste0("Subsetted to celltype: ", celltype_key, "\n"))
  } else {
    subset <- data[, data$cell_type_manual == celltype_key & data$celltype == cellsubtype_key]
    cat(paste0("Subsetted to celltype: ", celltype_key, " and subtype: ", cellsubtype_key, "\n"))
  }
  
  obj <- subset
  atac_counts <- GetAssayData(obj, assay = "ATAC", slot = "data")
  
  cat("Number of cells: ", ncol(obj), "\n")
  cat("Number of peaks (raw): ", nrow(atac_counts), "\n")
  
  cluster_ids <- obj[[leiden_col]][, 1]
  meta_counts <- scuttle::sumCountsAcrossCells(atac_counts, cluster_ids)
  expr_mat <- 1 * (assay(meta_counts, "sum") > 0)
  
  # Filter to valid peak names
  peak_names <- rownames(expr_mat)
  valid <- grepl("^chr.+-\\d+-\\d+$", peak_names)
  expr_mat <- expr_mat[valid, ]
  peak_names <- rownames(expr_mat)
  cat("Number of peaks (valid): ", nrow(expr_mat), "\n")
  
  parsed_coords <- strcapture(
    "^chr(\\w+)-(\\d+)-(\\d+)$",
    peak_names,
    data.frame(chr = character(), bp1 = integer(), bp2 = integer(), stringsAsFactors = FALSE)
  )
  
  # Build CDS
  cat("Build CDS...\n")
  cellinfo <- data.frame(cells = colnames(expr_mat), row.names = colnames(expr_mat))
  peakinfo <- data.frame(peaks = peak_names, gene_short_name = peak_names, row.names = peak_names)
  cds <- new_cell_data_set(expr_mat, cell_metadata = cellinfo, gene_metadata = peakinfo)
  cds <- cds[Matrix::rowSums(exprs(cds)) != 0, ]
  parsed_coords <- parsed_coords[peak_names %in% rownames(cds), , drop = FALSE]
  
  cat("Detect Genes, estimate size factor, preprocess, reduce dimension \n")
  cds <- detect_genes(cds)
  cds <- estimate_size_factors(cds)
  cds <- preprocess_cds(cds, method = "LSI")
  cds <- reduce_dimension(cds, reduction_method = "UMAP", preprocess_method = "LSI")
  
  # Add coordinates
  rowdata <- as.data.frame(rowData(cds))
  rowdata <- cbind(rowdata, parsed_coords)
  rowData(cds) <- rowdata
  
  fd <- new("AnnotatedDataFrame", data = rowdata)
  pd <- new("AnnotatedDataFrame", data = as.data.frame(colData(cds)))
  cat("Creating newCellDataSet...\n")
  cicero_cds <- newCellDataSet(expr_mat[rownames(cds), ],
                               phenoData = pd,
                               featureData = fd,
                               expressionFamily = binomialff(),
                               lowerDetectionLimit = 0)
  umap_coords <- reducedDims(cds)$UMAP
  reducedDimA(cicero_cds) <- as.matrix(umap_coords[colnames(cicero_cds), , drop = FALSE])
  
  # Run Cicero
  cat("1")
  chromosome_length <- read.table(genome_size_file)
  colnames(chromosome_length) <- c("chr", "length")
  cat("2")
  suffix <- if (is.null(cellsubtype_key)) celltype_key else paste(celltype_key, cellsubtype_key, sep = "_")
  outfile <- file.path(output_dir, paste0("cicero_connections_", suffix, ".tsv"))
  
  cat("Running Cicero...\n")
  cicero_cds <- detectGenes(cicero_cds, min_expr = 0.1)
  cicero_cds <- cicero:::estimateSizeFactorsSimp(cicero_cds)
  Biobase::exprs(cicero_cds) <- t(t(Biobase::exprs(cicero_cds)) / Biobase::pData(cicero_cds)$Size_Factor)
  
  conns <- run_cicero(cicero_cds, genomic_coords = chromosome_length)
  conns <- conns[, c("Peak1", "Peak2", "coaccess")]
  
  if (nrow(conns) > 0) {
    write.table(conns, outfile, sep = "\t", quote = FALSE, row.names = FALSE)
    cat(sprintf("Saved Cicero output to: %s\n", outfile))
  } else {
    warning("Cicero returned no connections.")
  }
}

# ---- CLI interface ----
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: Rscript run_cicero_batch.R <seurat_file> <celltype_tsv>")
}

seurat_file <- args[1]
celltype_tsv <- args[2]

df <- read.delim(celltype_tsv, stringsAsFactors = FALSE)
df$cellsubtype_key[df$cellsubtype_key == "NA"] <- NA

n_workers <- nrow(df)
param <- MulticoreParam(workers = n_workers)
register(param)

results <- bplapply(seq_len(n_workers), function(i) {
  x <- df[i, ]
  message("\n===== Running Cicero for: ", x$celltype_key,
          if (!is.na(x$cellsubtype_key)) paste0(" / ", x$cellsubtype_key) else "", " =====")
  
  tryCatch({
    run_cicero_on_celltype(
      seurat_file = seurat_file,
      celltype_key = x$celltype_key,
      cellsubtype_key = if (!is.na(x$cellsubtype_key)) x$cellsubtype_key else NULL
    )
    TRUE
  }, error = function(e) {
    message("Error: ", e$message)
    FALSE
  })
}, BPPARAM = param)

message("All jobs complete.")
