run_cicero_on_celltype <- function(
    seurat_file,
    celltype_key = "glia",
    cellsubtype_key = "astrocytes",
    leiden_col = "leiden",
    genome_size_file = "/mnt/germain/work/tomohl/pgBoost/mm10.chrom.sizes",
    output_dir = "Cicero/Cicero_OutputDir"
) {
  # Load required libraries 
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
  
  # Output Folder
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  
  # Load Seurat Object
  cat("Reading Seurat object...\n")
  data <- readRDS(seurat_file)
  data$cell_type_manual <- as.character(data$cell_type_manual)
  data$celltype <- as.character(data$celltype)
  cat(unique(data$cell_type_manual), unique(data$celltype))
  # Subset to Celltype
  if (is.null(cellsubtype_key)) {
    subset <- data[, data$cell_type_manual == celltype_key]
    cat(paste0("Subsetted to celltype:", celltype_key, "\n"))
  } else {
    subset <- data[, data$cell_type_manual == celltype_key & data$celltype == cellsubtype_key]
    cat(paste0("Subsetted to celltype:", celltype_key, " and cellsubtype:", cellsubtype_key, "\n"))
  }
  
  obj <- subset
  
  # Aggregate into Meta Cells
  cat("Aggregate into Meta Cells \n")
  atac_counts <- GetAssayData(obj, assay = "ATAC", slot = "data")
  
  cat("Number of cells in subset: ", ncol(obj), "\n")
  cat("Number of peaks (before filtering): ", nrow(atac_counts), "\n")
  
  cluster_ids <- obj$leiden
  meta_counts <- scuttle::sumCountsAcrossCells(atac_counts, cluster_ids)
  expr_mat <- 1 * (assay(meta_counts, "sum") > 0)
  
  # Filter to valid Peak names
  peak_names <- rownames(expr_mat)
  valid <- grepl("^chr.+-\\d+-\\d+$", peak_names)
  expr_mat <- expr_mat[valid, ]
  peak_names <- rownames(expr_mat)
  cat("Number of peaks (valid format): ", nrow(expr_mat), "\n")
  
  # parse coordinates
  parsed_coords <- strcapture(
    "^chr(\\w+)-(\\d+)-(\\d+)$",
    peak_names,
    data.frame(chr = character(), bp1 = integer(), bp2 = integer(), stringsAsFactors = FALSE)
  )
  
  # build CDS
  cat("Build CDS \n")
  cellinfo <- data.frame(cells = colnames(expr_mat), row.names = colnames(expr_mat))
  peakinfo <- data.frame(peaks = peak_names, gene_short_name = peak_names, row.names = peak_names)
  cds <- new_cell_data_set(expr_mat, cell_metadata = cellinfo, gene_metadata = peakinfo)
  cds <- cds[Matrix::rowSums(exprs(cds)) != 0, ]
  parsed_coords <- parsed_coords[peak_names %in% rownames(cds), , drop = FALSE]
  
  cds <- detect_genes(cds)
  cds <- estimate_size_factors(cds)
  cds <- preprocess_cds(cds, method = "LSI")
  cds <- reduce_dimension(cds, reduction_method = "UMAP", preprocess_method = "LSI")
  
  # add coordinates to CDS
  rowdata <- as.data.frame(rowData(cds))
  rowdata <- cbind(rowdata, parsed_coords)
  rowData(cds) <- rowdata
  
  # Convert to Cicero CDS 
  fd <- new("AnnotatedDataFrame", data = rowdata)
  pd <- new("AnnotatedDataFrame", data = as.data.frame(colData(cds)))
  cicero_cds <- newCellDataSet(expr_mat[rownames(cds), ],
                               phenoData = pd,
                               featureData = fd,
                               expressionFamily = binomialff(),
                               lowerDetectionLimit = 0)
  umap_coords <- reducedDims(cds)$UMAP
  umap_coords <- umap_coords[colnames(cicero_cds), , drop = FALSE]
  reducedDimA(cicero_cds) <- as.matrix(umap_coords)
  
  # Load genome size 
  chromosome_length <- read.table(genome_size_file)
  colnames(chromosome_length) <- c("chr", "length")
  
  # Save output path early
  suffix <- if (is.null(cellsubtype_key)) celltype_key else paste(celltype_key, cellsubtype_key, sep = "_")
  outfile <- file.path(output_dir, paste0("cicero_connections_", suffix, ".tsv"))
  
  #  Run Cicero 
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
    warning("Cicero returned no connections. No file written.")
  }
  
  write.table(conns, outfile, sep = "\t", quote = FALSE, row.names = FALSE)
  cat(sprintf("Saved Cicero output to: %s\n", outfile))
}

# --- CLI interface ---
args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 2) {
  stop("Usage: Rscript run_cicero_batch.R <seurat_file> <output_dir> [celltype_tsv]")
}

seurat_file <- args[1]
output_dir <- args[2]
celltype_tsv <- ifelse(length(args) >= 3, args[3], NULL)

# Load BiocParallel
suppressPackageStartupMessages(library(BiocParallel))

if (!is.null(celltype_tsv)) {
  # Load table of celltypes to process
  df <- read.delim(celltype_tsv, stringsAsFactors = FALSE)
  df$cellsubtype_key[df$cellsubtype_key == "NA"] <- NA  # convert string "NA" to actual NA
} else {
  # fallback: single run
  stop("You must provide a TSV with celltype_key and cellsubtype_key for parallel processing.")
}

# Define number of workers = number of rows
n_workers <- nrow(df)
param <- MulticoreParam(workers = n_workers)
register(param)

# Run in parallel
results <- bplapply(seq_len(n_workers), function(i) {
  x <- df[i, ]
  message("\n===== Running Cicero for: ", x$celltype_key,
          if (!is.na(x$cellsubtype_key)) paste0(" / ", x$cellsubtype_key) else "", " =====")
  
  tryCatch({
    run_cicero_on_celltype(
      seurat_file = seurat_file,
      celltype_key = x$celltype_key,
      cellsubtype_key = if (!is.na(x$cellsubtype_key)) x$cellsubtype_key else NULL,
      output_dir = output_dir
    )
    TRUE
  }, error = function(e) {
    message("Error: ", e$message)
    FALSE
  })
}, BPPARAM = param)

# Done
message("All jobs complete.")