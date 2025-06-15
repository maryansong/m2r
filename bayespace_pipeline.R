
Sys.setenv(OMP_NUM_THREADS = "1")  # Prevent OpenMP threading errors on macOS
knitr::opts_chunk$set(echo = TRUE)

# Install BiocManager if not already present
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

# Install required packages using BiocManager
required_pkgs <- c("zellkonverter", "SpatialExperiment", "BayesSpace", "SingleCellExperiment", "fastmap")

for (pkg in required_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    BiocManager::install(pkg, ask = FALSE)
  }
}

# Update fastmap if below minimum version
if (packageVersion("fastmap") < "1.2.0") install.packages("fastmap")

# Load all required libraries
library(zellkonverter)
library(SpatialExperiment)
library(BayesSpace)
library(SingleCellExperiment)
library(ggplot2)

# Function to process a single sample through BayesSpace
run_bayespace_pipeline <- function(h5ad_file, pca_csv, sample_id, q = 7, d = 15, fast = FALSE, out_pdf = TRUE) {
  cat(paste0("\nProcessing sample: ", sample_id, "\n"))
  sce <- readH5AD(h5ad_file)
  spe <- as(sce, "SpatialExperiment")
  pca <- read.csv(pca_csv, row.names = 1)

  # Clean names to R-compatible form and ensure valid row/col names
  colnames(spe) <- make.names(colnames(spe), unique = TRUE)
  rownames(pca) <- make.names(rownames(pca), unique = TRUE)

  colnames(spe)[is.na(colnames(spe)) | colnames(spe) == ""] <- paste0("Spot", seq_len(sum(is.na(colnames(spe)) | colnames(spe) == "")))
  rownames(pca)[is.na(rownames(pca)) | rownames(pca) == ""] <- paste0("Spot", seq_len(sum(is.na(rownames(pca)) | rownames(pca) == "")))

  # Check for matching barcodes and align in same order
  common <- intersect(colnames(spe), rownames(pca))
  if (length(common) == 0) {
    cat("No matching spot names between PCA and SPE\n")
    cat("colnames(spe):\n")
    print(head(colnames(spe)))
    cat("rownames(pca):\n")
    print(head(rownames(pca)))
    stop("Aborting due to mismatched names.")
  }

  spe <- spe[, common]
  pca <- pca[common, , drop = FALSE]

  if (ncol(spe) != nrow(pca)) {
    stop(paste("Dimension mismatch:", ncol(spe), "columns in spe vs", nrow(pca), "rows in PCA"))
  }

  reducedDims(spe)$PCA <- as.matrix(pca)

  if (!("row" %in% colnames(colData(spe))) || !("col" %in% colnames(colData(spe)))) {
    row_candidates <- grep("row", colnames(colData(spe)), value = TRUE)
    col_candidates <- grep("col", colnames(colData(spe)), value = TRUE)
    if (length(row_candidates) > 0 && length(col_candidates) > 0) {
      colData(spe)$row <- colData(spe)[[row_candidates[1]]]
      colData(spe)$col <- colData(spe)[[col_candidates[1]]]
    } else {
      stop("Could not find valid 'row' and 'col' columns in colData.")
    }
  }

  # BayesSpace clustering
  spe <- spatialCluster(
    spe,
    q = q,
    d = d,
    platform = "Visium",
    use.dimred = "PCA",
    nrep = ifelse(fast, 1000, 50000),
    burn.in = ifelse(fast, 100, 1000)
  )

  df <- as.data.frame(colData(spe))
  p <- ggplot(df, aes(x = col, y = row, fill = factor(spatial.cluster))) +
    geom_tile() +
    coord_equal() +
    scale_y_reverse() +
    labs(title = paste("BayesSpace Clustering:", sample_id), fill = "Cluster") +
    theme_minimal()

  ggsave(paste0(sample_id, "_BayesSpace_cluster.png"), plot = p, width = 6, height = 6)
  if (out_pdf) {
    ggsave(paste0(sample_id, "_BayesSpace_cluster.pdf"), plot = p, width = 6, height = 6)
  }

  saveRDS(spe, file = paste0(sample_id, "_BayesSpace_result.rds"))
  cat("Done.\n")
}

run_bayespace_pipeline("151673.h5ad", "151673_pca_embeddings.csv", "151673")
