suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(SingleCellExperiment)
  library(Matrix)
  library(zellkonverter)
  library(ggplot2)
  library(patchwork)
})

options(future.globals.maxSize = 3e9)
set.seed(0)

cosmx_h5ad <- "colon_adata_for_seurat.h5ad"
sc1_h5ad <- "combined_adata_first_half.h5ad"
sc2_h5ad <- "combined_adata_second_half.h5ad"

nfeatures <- 3000
npcs <- 25
resolution <- 0.8

out_integrated_rds <- "cosmx_scRNA_rpca_integrated.rds"
out_cosmx_rds <- "cosmx_rpca_reclustered.rds"
out_umap_png <- "cosmx_umap_rpca.png"

cosmx_counts_mtx <- "colon_adata_for_seurat_counts.mtx"
cosmx_cells_tsv <- "colon_adata_for_seurat_cells.tsv"
cosmx_genes_tsv <- "colon_adata_for_seurat_genes.tsv"
cosmx_obs_csv <- "colon_adata_for_seurat_obs.csv"

load_h5ad <- function(path) {
  if (!file.exists(path)) {
    stop("Missing file: ", path)
  }
  sce <- zellkonverter::readH5AD(path)
  assay_names <- as.character(SummarizedExperiment::assayNames(sce))
  if (length(assay_names) == 0) {
    stop("No assays found in: ", path)
  }

  mat_name <- if ("X" %in% assay_names) "X" else assay_names[1]
  mat <- SummarizedExperiment::assay(sce, mat_name)
  if (!inherits(mat, "dgCMatrix")) {
    mat <- as(mat, "dgCMatrix")
  }

  # detect if mat is transposed (cells x genes) vs expected genes x cells
  sce_rn <- rownames(sce)
  sce_cn <- colnames(sce)
  if (!is.null(sce_rn) && !is.null(sce_cn)) {
    if (nrow(mat) == length(sce_cn) && ncol(mat) == length(sce_rn)) {
      mat <- Matrix::t(mat)
    }
  }

  rn <- rownames(mat)
  if (is.null(rn) || anyNA(rn) || any(rn == "")) {
    rn <- rownames(sce)
  }
  if (is.null(rn) || anyNA(rn) || any(rn == "")) {
    rn <- paste0("gene_", seq_len(nrow(mat)))
  }

  cn <- colnames(mat)
  if (is.null(cn) || anyNA(cn) || any(cn == "")) {
    cn <- colnames(sce)
  }
  if (is.null(cn) || anyNA(cn) || any(cn == "")) {
    cdat <- as.data.frame(SummarizedExperiment::colData(sce))
    if ("barcodes" %in% colnames(cdat)) {
      cn <- as.character(cdat$barcodes)
    } else {
      cn <- paste0("cell_", seq_len(ncol(mat)))
    }
  }

  rownames(mat) <- make.unique(as.character(rn))
  colnames(mat) <- make.unique(as.character(cn))

  meta <- as.data.frame(SummarizedExperiment::colData(sce))
  if (nrow(meta) != ncol(mat)) {
    meta <- data.frame(row.names = colnames(mat))
  } else {
    rownames(meta) <- colnames(mat)
    # Fix type mismatches: convert numeric columns that are whole numbers to integer
    for (col in colnames(meta)) {
      if (is.numeric(meta[[col]]) && !is.integer(meta[[col]])) {
        if (all(meta[[col]] == floor(meta[[col]]), na.rm = TRUE)) {
          meta[[col]] <- as.integer(meta[[col]])
        }
      }
    }
  }

  obj <- CreateSeuratObject(
    counts = mat,
    meta.data = meta,
    assay = "RNA",
    project = tools::file_path_sans_ext(basename(path))
  )

  DefaultAssay(obj) <- "RNA"
  obj
}

ensure_sct <- function(obj) {
  if ("SCT" %in% Assays(obj)) {
    return(obj)
  }
  assay_name <- DefaultAssay(obj)
  SCTransform(obj, assay = assay_name, verbose = FALSE)
}

load_sc_half <- function(path, prefix) {
  obj <- load_h5ad(path)
  obj <- RenameCells(obj, add.cell.id = prefix)
  obj
}

merge_sc_halves <- function(obj1, obj2) {
  # Avoid problematic Seurat::merge() by manually concatenating
  # Get counts matrices
  mat1 <- GetAssayData(obj1, assay = "RNA", layer = "counts")
  mat2 <- GetAssayData(obj2, assay = "RNA", layer = "counts")
  
  # Ensure same genes and order
  common_genes <- intersect(rownames(mat1), rownames(mat2))
  if (length(common_genes) == 0) {
    stop("No common genes between the two datasets")
  }
  
  mat1 <- mat1[common_genes, ]
  mat2 <- mat2[common_genes, ]
  
  # Concatenate matrices (cbind = add columns/cells)
  mat_combined <- Matrix::cbind2(mat1, mat2)
  
  # Concatenate metadata
  meta1 <- obj1@meta.data
  meta2 <- obj2@meta.data
  
  # Ensure metadata columns are compatible
  all_cols <- union(colnames(meta1), colnames(meta2))
  for (col in all_cols) {
    if (!(col %in% colnames(meta1))) meta1[[col]] <- NA
    if (!(col %in% colnames(meta2))) meta2[[col]] <- NA
  }
  meta_combined <- rbind(meta1[all_cols], meta2[all_cols])
  
  # Create merged Seurat object
  obj_merged <- CreateSeuratObject(
    counts = mat_combined,
    meta.data = meta_combined,
    assay = "RNA",
    project = "scRNA_reference"
  )
  
  DefaultAssay(obj_merged) <- "RNA"
  obj_merged
}

load_cosmx_export <- function(counts_mtx, cells_tsv, genes_tsv, obs_csv = NULL) {
  required_files <- c(counts_mtx, cells_tsv, genes_tsv)
  missing_files <- required_files[!file.exists(required_files)]
  if (length(missing_files) > 0) {
    stop("Missing CosMx export files: ", paste(missing_files, collapse = ", "))
  }

  mat <- Matrix::readMM(counts_mtx)
  if (!inherits(mat, "dgCMatrix")) {
    mat <- as(mat, "dgCMatrix")
  }
  cells <- read.delim(cells_tsv, header = FALSE, stringsAsFactors = FALSE)
  genes <- read.delim(genes_tsv, header = FALSE, stringsAsFactors = FALSE)

  cells <- as.character(cells[[1]])
  genes <- as.character(genes[[1]])

  if (nrow(mat) != length(genes) || ncol(mat) != length(cells)) {
    stop(
      "CosMx export dimensions do not match names: matrix is ",
      nrow(mat), " x ", ncol(mat),
      ", genes = ", length(genes),
      ", cells = ", length(cells)
    )
  }

  rownames(mat) <- make.unique(genes)
  colnames(mat) <- make.unique(cells)

  meta <- data.frame(row.names = colnames(mat))
  if (!is.null(obs_csv) && file.exists(obs_csv)) {
    obs <- read.csv(obs_csv, stringsAsFactors = FALSE)
    if ("cell_id" %in% colnames(obs)) {
      obs$cell_id <- as.character(obs$cell_id)
      rownames(obs) <- obs$cell_id
      obs$cell_id <- NULL
    } else if ("cell" %in% colnames(obs)) {
      obs$cell <- as.character(obs$cell)
      rownames(obs) <- obs$cell
      obs$cell <- NULL
    }
    if (nrow(obs) == ncol(mat)) {
      # Ensure consistent numeric types: convert whole-number doubles to integers
      for (col in colnames(obs)) {
        if (is.numeric(obs[[col]]) && !is.integer(obs[[col]])) {
          if (all(obs[[col]] == floor(obs[[col]]), na.rm = TRUE)) {
            obs[[col]] <- as.integer(obs[[col]])
          }
        }
      }
      meta <- obs[colnames(mat), , drop = FALSE]
    }
  }

  obj <- CreateSeuratObject(
    counts = mat,
    meta.data = meta,
    assay = "RNA",
    project = "colon_adata_for_seurat"
  )

  DefaultAssay(obj) <- "RNA"
  obj
}

sc_ref <- merge_sc_halves(
  load_sc_half(sc1_h5ad, "scRNA1"),
  load_sc_half(sc2_h5ad, "scRNA2")
)
gc()

query <- load_cosmx_export(
  counts_mtx = cosmx_counts_mtx,
  cells_tsv = cosmx_cells_tsv,
  genes_tsv = cosmx_genes_tsv,
  obs_csv = cosmx_obs_csv
)
query <- RenameCells(query, add.cell.id = "CosMx")

sc_ref <- ensure_sct(sc_ref)
query <- ensure_sct(query)

obj_list <- list(sc_ref, query)

################################################################################

features <- SelectIntegrationFeatures(object.list = obj_list, nfeatures = nfeatures)

obj_list <- lapply(obj_list, PrepSCTIntegration, anchor.features = features, verbose = FALSE)
obj_list <- lapply(obj_list, RunPCA, features = features, npcs = npcs, verbose = FALSE)

anchors <- FindIntegrationAnchors(
  object.list = obj_list,
  normalization.method = "SCT",
  anchor.features = features,
  reduction = "rpca",
  reference = 1,
  dims = 1:npcs,
  verbose = FALSE
)

integrated <- IntegrateData(
  anchorset = anchors,
  normalization.method = "SCT",
  verbose = FALSE
)

DefaultAssay(integrated) <- "integrated"
integrated <- RunPCA(integrated, npcs = npcs, verbose = FALSE)

integrated$dataset <- ifelse(grepl("^CosMx_", colnames(integrated)), "CosMx", "scRNA")
cosmx <- subset(integrated, cells = colnames(integrated)[integrated$dataset == "CosMx"])

cosmx <- RunUMAP(cosmx, dims = 1:npcs, reduction = "pca", verbose = FALSE)
cosmx <- FindNeighbors(cosmx, dims = 1:npcs, verbose = FALSE)
cosmx <- FindClusters(cosmx, resolution = resolution, verbose = FALSE)

p1 <- DimPlot(cosmx, group.by = "seurat_clusters", label = TRUE) + NoLegend()

png(out_umap_png, width = 1200, height = 900, res = 200)
print(p1)
dev.off()

saveRDS(integrated, file = out_integrated_rds)
saveRDS(cosmx, file = out_cosmx_rds)
message("Saved integrated object: ", out_integrated_rds)
message("Saved CosMx object: ", out_cosmx_rds)
message("Saved CosMx UMAP: ", out_umap_png)
