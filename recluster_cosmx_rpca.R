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

sc_ref <- merge(
  load_sc_half(sc1_h5ad, "scRNA1"),
  y = load_sc_half(sc2_h5ad, "scRNA2"),
  project = "scRNA_reference"
)
gc()

query <- load_h5ad(cosmx_h5ad)
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
