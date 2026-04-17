# ============================================================
# run_fastreseg_all_fov.R
# ============================================================
# Full FastReseg run across all FOVs using the same workflow
# that worked for single FOV runs.
#
# Outputs per FOV:
# - FastReseg_full_fov_<FOV>.rds
# - FastReseg_summary_FOV_<FOV>.txt
# - FastReseg_changed_transcripts_FOV_<FOV>.csv
# - FastReseg_cell_rename_summary_FOV_<FOV>.csv
# - tutorial_style_before_<cell>.png
# - tutorial_style_after_<cell>.png
# - top_changed_before_<cell>.png
# - top_changed_after_<cell>.png
#
# Combined outputs:
# - all_fov_run_status.csv
# - all_fov_summary_metrics.csv
# - all_fov_changed_transcripts.csv
# - all_fov_cell_rename_summary.csv
# ============================================================

# -------------------------------
# 1) Load packages
# -------------------------------
library(data.table)
library(Matrix)
library(FastReseg)
library(ggplot2)
library(scales)

# -------------------------------
# 2) Define input files
# -------------------------------
expr_file   <- "Raw Data/Colon_TMA_exprMat_file.csv.gz"
meta_file   <- "Raw Data/Colon_TMA_metadata_file.csv.gz"
tx_file     <- "Raw Data/Colon_TMA_tx_file.csv.gz"
poly_file   <- "Raw Data/Colon_TMA-polygons.csv.gz"
fovpos_file <- "Raw Data/Colon_TMA_fov_positions_file.csv.gz"

cluster_col <- "RNA_RNA.QC.Module_Cell.Typing.InSituType.1_1_clusters"

# -------------------------------
# 3) Global settings
# -------------------------------
pixel_size_um <- 0.12028
zstep_size_um <- 0.8

percent_cores <- 0.25

# Set to NULL to run all FOVs
fov_subset <- NULL
# example:
# fov_subset <- c(7, 232)

# -------------------------------
# 4) Output folders
# -------------------------------
root_out <- "FastReseg_allFOV_same_workflow"
dir.create(root_out, recursive = TRUE, showWarnings = FALSE)

status_dir  <- file.path(root_out, "00_status")
fov_root    <- file.path(root_out, "01_per_fov")
summary_dir <- file.path(root_out, "02_combined_outputs")

dir.create(status_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fov_root, recursive = TRUE, showWarnings = FALSE)
dir.create(summary_dir, recursive = TRUE, showWarnings = FALSE)

cat("Root output folder:\n", normalizePath(root_out), "\n\n")

# -------------------------------
# 5) Read all metadata once
# -------------------------------
cat("Reading metadata file once...\n")
metadata_all <- fread(meta_file)

if (!"cell" %in% colnames(metadata_all)) {
  metadata_all[, cell := paste0("c_1_", fov, "_", cell_ID)]
}

cat("Metadata rows:", nrow(metadata_all), "\n")
cat("Metadata columns:", ncol(metadata_all), "\n\n")

# -------------------------------
# 6) Discover FOV list
# -------------------------------
fov_positions <- fread(fovpos_file)
all_fovs <- sort(unique(as.integer(fov_positions$FOV)))

if (!is.null(fov_subset)) {
  all_fovs <- intersect(all_fovs, as.integer(fov_subset))
}

cat("Total FOVs to run:", length(all_fovs), "\n")
cat("First 20 FOVs:\n")
print(head(all_fovs, 20))
cat("\n")

# -------------------------------
# 7) Containers for combined outputs
# -------------------------------
run_status_list <- list()
summary_metrics_list <- list()
all_changed_transcripts_list <- list()
all_cell_rename_summary_list <- list()

# -------------------------------
# 8) Loop over FOVs
# -------------------------------
for (idx in seq_along(all_fovs)) {
  
  test_fov <- all_fovs[idx]
  
  cat("\n====================================================\n")
  cat(sprintf("[%d / %d] Running FOV %d\n", idx, length(all_fovs), test_fov))
  cat("====================================================\n")
  
  out_dir <- file.path(fov_root, paste0("FOV_", test_fov))
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  
  start_file <- file.path(status_dir, paste0("FOV_", test_fov, "_STARTED.txt"))
  writeLines(
    c(
      paste("FOV", test_fov, "started"),
      paste("Time:", as.character(Sys.time()))
    ),
    start_file
  )
  
  fov_result <- tryCatch({
    
    # -------- 3) Read raw metadata for one FOV only --------
    meta_fov <- copy(metadata_all[fov == test_fov])
    
    cat("Metadata rows in FOV", test_fov, ":", nrow(meta_fov), "\n")
    cat("Metadata columns:", ncol(meta_fov), "\n\n")
    
    if (nrow(meta_fov) == 0) stop("meta_fov has 0 rows")
    
    if (!"cell" %in% colnames(meta_fov)) {
      meta_fov[, cell := paste0("c_1_", fov, "_", cell_ID)]
    }
    
    # -------- 4) Read raw expression matrix for one FOV only --------
    expr_fov <- fread(
      cmd = sprintf(
        "gzip -dc %s | awk -F',' 'NR==1 || $1==%d'",
        shQuote(expr_file), test_fov
      )
    )
    
    cat("Expression rows in FOV", test_fov, ":", nrow(expr_fov), "\n")
    cat("Expression columns:", ncol(expr_fov), "\n\n")
    
    if (nrow(expr_fov) == 0) stop("expr_fov has 0 rows")
    
    expr_fov[, cell := paste0("c_1_", fov, "_", cell_ID)]
    
    # -------- 5) Match metadata cells and expression cells --------
    common_cells <- intersect(expr_fov$cell, meta_fov$cell)
    
    expr_fov <- expr_fov[cell %in% common_cells]
    meta_fov <- meta_fov[cell %in% common_cells]
    
    setkey(expr_fov, cell)
    setkey(meta_fov, cell)
    expr_fov <- expr_fov[meta_fov$cell]
    
    cat("Matched cells between metadata and expression:", length(common_cells), "\n\n")
    
    if (length(common_cells) == 0) stop("No matched cells between metadata and expression")
    
    # -------- 6) Build cell x gene counts matrix --------
    id_cols <- c("fov", "cell_ID", "cell")
    gene_cols <- setdiff(colnames(expr_fov), id_cols)
    
    counts_fov <- as.matrix(expr_fov[, ..gene_cols])
    rownames(counts_fov) <- expr_fov$cell
    storage.mode(counts_fov) <- "numeric"
    counts_fov <- Matrix(counts_fov, sparse = TRUE)
    
    cat("Counts matrix dim (cells x genes):", dim(counts_fov)[1], "x", dim(counts_fov)[2], "\n\n")
    
    # -------- 7) Extract cluster labels for each cell --------
    clust_fov <- as.character(meta_fov[[cluster_col]])
    names(clust_fov) <- meta_fov$cell
    
    cat("Number of unique clusters in pilot:", length(unique(clust_fov)), "\n")
    print(table(clust_fov))
    cat("\n")
    
    # -------- 8) Create a one-FOV transcript file --------
    tx_fov_file <- file.path(out_dir, paste0("Colon_TMA_tx_fov_", test_fov, ".csv"))
    
    cmd_tx <- sprintf(
      "gzip -dc %s | awk -F',' 'NR==1 || $1==%d' > %s",
      shQuote(tx_file), test_fov, shQuote(tx_fov_file)
    )
    
    system(cmd_tx)
    
    cat("Transcript pilot file written to:\n", normalizePath(tx_fov_file), "\n\n")
    
    # -------- 9) Read one-FOV transcript file --------
    tx_fov <- fread(tx_fov_file)
    
    cat("Transcript rows in pilot FOV:", nrow(tx_fov), "\n")
    cat("Transcript columns:\n")
    print(colnames(tx_fov))
    cat("\n")
    
    if (nrow(tx_fov) == 0) stop("tx_fov has 0 rows")
    
    # -------- 10) Define extracellular / unassigned transcript bucket --------
    extracellular_cellID <- paste0("c_1_", test_fov, "_0")
    
    cat("Extracellular cell label:", extracellular_cellID, "\n")
    print(table(tx_fov$cell == extracellular_cellID, useNA = "ifany"))
    cat("\n")
    
    # -------- 11) Keep only genes shared between counts and transcript file --------
    genes_tx <- unique(tx_fov$target)
    common_genes <- intersect(colnames(counts_fov), genes_tx)
    
    counts_fov <- counts_fov[, common_genes, drop = FALSE]
    
    cat("Genes in transcript file:", length(genes_tx), "\n")
    cat("Genes shared with counts matrix:", length(common_genes), "\n\n")
    
    if (length(common_genes) == 0) stop("No common genes between counts and transcript file")
    
    # -------- 12) Keep only cells that appear in transcript file --------
    tx_cells_assigned <- unique(tx_fov$cell)
    tx_cells_assigned <- setdiff(tx_cells_assigned, extracellular_cellID)
    
    pilot_cells <- intersect(rownames(counts_fov), tx_cells_assigned)
    
    counts_fov <- counts_fov[pilot_cells, , drop = FALSE]
    clust_fov <- clust_fov[pilot_cells]
    
    cat("Cells used in FastReseg pilot:", length(pilot_cells), "\n")
    cat("Final counts matrix dim:", dim(counts_fov)[1], "x", dim(counts_fov)[2], "\n\n")
    
    if (length(pilot_cells) == 0) stop("No pilot cells remained after transcript filtering")
    
    # -------- 13) Run FULL FastReseg pipeline --------
    set.seed(1)
    
    res_full <- fastReseg_full_pipeline(
      counts = counts_fov,
      clust = clust_fov,
      refProfiles = NULL,
      
      pixel_size = pixel_size_um,
      zstep_size = zstep_size_um,
      
      transcript_df = tx_fov,
      transID_coln = NULL,
      transGene_coln = "target",
      cellID_coln = "cell",
      spatLocs_colns = c("x_local_px", "y_local_px", "z"),
      invert_y = TRUE,
      
      extracellular_cellID = extracellular_cellID,
      
      flagModel_TransNum_cutoff = 50,
      flagCell_lrtest_cutoff = 2,
      svmClass_score_cutoff = -2,
      
      molecular_distance_cutoff = NULL,
      cellular_distance_cutoff = NULL,
      score_baseline = NULL,
      lowerCutoff_transNum = NULL,
      higherCutoff_transNum = NULL,
      
      groupTranscripts_method = "dbscan",
      spatialMergeCheck_method = "geometryDiff",
      cutoff_spatialMerge = 0.5,
      
      imputeFlag_missingCTs = TRUE,
      
      path_to_output = out_dir,
      transDF_export_option = 2,
      save_intermediates = TRUE,
      return_perCellData = TRUE,
      combine_extra = FALSE,
      
      seed_process = 1,
      percentCores = percent_cores
    )
    
    # -------- 14) Inspect returned objects --------
    cat("=== RESULT OBJECT NAMES ===\n")
    print(names(res_full))
    cat("\n")
    
    updated_transDF <- NULL
    if ("updated_transDF_list" %in% names(res_full) && length(res_full$updated_transDF_list) >= 1) {
      updated_transDF <- as.data.table(res_full$updated_transDF_list[[1]])
    }
    
    if (is.null(updated_transDF) || nrow(updated_transDF) == 0) {
      stop("updated_transDF not found in returned object or is empty")
    }
    
    cat("=== updated_transDF preview ===\n")
    print(dim(updated_transDF))
    print(head(updated_transDF))
    cat("\nColumns:\n")
    print(colnames(updated_transDF))
    cat("\n")
    
    # -------- 15) Tutorial-style before/after transcript assignment plots --------
    poly_fov <- fread(
      cmd = sprintf(
        "gzip -dc %s | awk -F',' 'NR==1 || $1==%d'",
        shQuote(poly_file), test_fov
      )
    )
    
    changed_transcripts <- updated_transDF[UMI_cellID != updated_cellID]
    changed_cells <- unique(changed_transcripts$UMI_cellID)
    
    if (length(changed_cells) == 0) {
      changed_cells <- unique(updated_transDF$UMI_cellID)[1:3]
    }
    
    target_cells <- unique(changed_cells)[1:min(3, length(unique(changed_cells)))]
    
    cat("Target cells selected for tutorial-style plots:\n")
    print(target_cells)
    cat("\n")
    
    poly_fov[, x_um := x_local_px * pixel_size_um]
    poly_fov[, y_um := y_local_px * pixel_size_um]
    
    tx_before <- copy(updated_transDF)
    tx_before[, plot_cellID := UMI_cellID]
    tx_before[, x_um := x]
    tx_before[, y_um := y]
    
    tx_after <- copy(updated_transDF)
    tx_after[, plot_cellID := updated_cellID]
    tx_after[, x_um := x]
    tx_after[, y_um := y]
    
    get_plot_limits <- function(df, pad = 5) {
      list(
        xlim = range(df$x_um, na.rm = TRUE) + c(-pad, pad),
        ylim = range(df$y_um, na.rm = TRUE) + c(-pad, pad)
      )
    }
    
    get_local_neighborhood <- function(target_cell, radius_px = 120) {
      target_row <- meta_fov[cell == target_cell]
      if (nrow(target_row) == 0) return(NULL)
      
      cx <- target_row$CenterX_local_px[1]
      cy <- target_row$CenterY_local_px[1]
      
      meta_local <- copy(meta_fov)
      meta_local[, dist_px := sqrt((CenterX_local_px - cx)^2 + (CenterY_local_px - cy)^2)]
      near_cells <- meta_local[dist_px <= radius_px, cell]
      
      poly_local <- copy(poly_fov[cell %in% near_cells])
      tx_before_local <- copy(tx_before[UMI_cellID %in% near_cells | updated_cellID %in% near_cells])
      tx_after_local  <- copy(tx_after[UMI_cellID %in% near_cells | updated_cellID %in% near_cells])
      
      list(
        target_cell = target_cell,
        poly_local = poly_local,
        tx_before_local = tx_before_local,
        tx_after_local = tx_after_local
      )
    }
    
    for (tc in target_cells) {
      nb <- get_local_neighborhood(tc, radius_px = 120)
      if (is.null(nb)) next
      
      lims_before <- get_plot_limits(nb$tx_before_local)
      
      p_before <- ggplot() +
        geom_path(
          data = nb$poly_local,
          aes(x = x_um, y = y_um, group = cell),
          color = "grey60",
          linewidth = 0.3,
          alpha = 0.7
        ) +
        geom_point(
          data = nb$tx_before_local,
          aes(x = x_um, y = y_um, color = plot_cellID),
          size = 1.2,
          alpha = 0.85
        ) +
        coord_fixed(xlim = lims_before$xlim, ylim = lims_before$ylim) +
        labs(
          title = paste("Original transcript assignment near", tc),
          x = "x (micron)",
          y = "y (micron)",
          color = "original cell ID"
        ) +
        theme_minimal(base_size = 14) +
        theme(
          legend.position = "right",
          plot.title = element_text(size = 16)
        )
      
      ggsave(
        filename = file.path(out_dir, paste0("tutorial_style_before_", tc, ".png")),
        plot = p_before,
        width = 10,
        height = 8,
        dpi = 300
      )
      
      lims_after <- get_plot_limits(nb$tx_after_local)
      
      p_after <- ggplot() +
        geom_path(
          data = nb$poly_local,
          aes(x = x_um, y = y_um, group = cell),
          color = "grey60",
          linewidth = 0.3,
          alpha = 0.7
        ) +
        geom_point(
          data = nb$tx_after_local,
          aes(x = x_um, y = y_um, color = plot_cellID),
          size = 1.2,
          alpha = 0.85
        ) +
        coord_fixed(xlim = lims_after$xlim, ylim = lims_after$ylim) +
        labs(
          title = paste("Updated transcript assignment near", tc),
          x = "x (micron)",
          y = "y (micron)",
          color = "updated cell ID"
        ) +
        theme_minimal(base_size = 14) +
        theme(
          legend.position = "right",
          plot.title = element_text(size = 16)
        )
      
      ggsave(
        filename = file.path(out_dir, paste0("tutorial_style_after_", tc, ".png")),
        plot = p_after,
        width = 10,
        height = 8,
        dpi = 300
      )
    }
    
    saveRDS(res_full, file.path(out_dir, paste0("FastReseg_full_fov_", test_fov, ".rds")))
    cat("Saved full FastReseg result to:\n", file.path(out_dir, paste0("FastReseg_full_fov_", test_fov, ".rds")), "\n")
    
    # -------- 16) Plot TOP 5 most changed cells --------
    updated_transDF[, changed := UMI_cellID != updated_cellID]
    
    cell_change_counts <- updated_transDF[changed == TRUE, .N, by = UMI_cellID]
    setorder(cell_change_counts, -N)
    
    top_cells <- head(cell_change_counts$UMI_cellID, 5)
    
    cat("Top 5 cells with most changed transcripts:\n")
    print(top_cells)
    cat("\n")
    
    if (!"x_um" %in% colnames(poly_fov)) poly_fov[, x_um := x_local_px * pixel_size_um]
    if (!"y_um" %in% colnames(poly_fov)) poly_fov[, y_um := y_local_px * pixel_size_um]
    
    tx_before <- copy(updated_transDF)
    tx_before[, plot_cellID := UMI_cellID]
    tx_before[, x_um := x]
    tx_before[, y_um := y]
    
    tx_after <- copy(updated_transDF)
    tx_after[, plot_cellID := updated_cellID]
    tx_after[, x_um := x]
    tx_after[, y_um := y]
    
    get_plot_limits <- function(df, pad = 5) {
      list(
        xlim = range(df$x_um, na.rm = TRUE) + c(-pad, pad),
        ylim = range(df$y_um, na.rm = TRUE) + c(-pad, pad)
      )
    }
    
    make_local_color_map <- function(before_ids, after_ids) {
      before_ids <- sort(unique(as.character(before_ids)))
      after_ids  <- sort(unique(as.character(after_ids)))
      
      shared_ids   <- intersect(before_ids, after_ids)
      before_only  <- setdiff(before_ids, after_ids)
      after_only   <- setdiff(after_ids, before_ids)
      
      ordered_ids <- c(shared_ids, before_only, after_only)
      
      n_ids <- length(ordered_ids)
      cols <- hue_pal()(max(n_ids, 1))
      names(cols) <- ordered_ids
      
      cols
    }
    
    for (tc in top_cells) {
      target_row <- meta_fov[cell == tc]
      if (nrow(target_row) == 0) next
      
      cx <- target_row$CenterX_local_px[1]
      cy <- target_row$CenterY_local_px[1]
      
      meta_fov[, dist := sqrt((CenterX_local_px - cx)^2 + (CenterY_local_px - cy)^2)]
      near_cells <- meta_fov[dist <= 120, cell]
      
      poly_local <- poly_fov[cell %in% near_cells]
      tx_before_local <- tx_before[UMI_cellID %in% near_cells | updated_cellID %in% near_cells]
      tx_after_local  <- tx_after[UMI_cellID %in% near_cells | updated_cellID %in% near_cells]
      
      local_colors <- make_local_color_map(
        before_ids = tx_before_local$plot_cellID,
        after_ids  = tx_after_local$plot_cellID
      )
      
      lims_before <- get_plot_limits(tx_before_local)
      
      p_before <- ggplot() +
        geom_path(
          data = poly_local,
          aes(x = x_um, y = y_um, group = cell),
          color = "grey55",
          linewidth = 0.35,
          alpha = 0.8
        ) +
        geom_point(
          data = tx_before_local,
          aes(x = x_um, y = y_um, color = plot_cellID),
          size = 1.4,
          alpha = 0.9
        ) +
        scale_color_manual(values = local_colors, drop = FALSE) +
        coord_fixed(xlim = lims_before$xlim, ylim = lims_before$ylim) +
        labs(
          title = paste("BEFORE -", tc),
          x = "x (micron)",
          y = "y (micron)",
          color = "original cell ID"
        ) +
        theme_minimal(base_size = 14) +
        theme(
          plot.title = element_text(size = 16),
          legend.position = "right"
        )
      
      ggsave(
        filename = file.path(out_dir, paste0("top_changed_before_", tc, ".png")),
        plot = p_before,
        width = 10,
        height = 8,
        dpi = 300
      )
      
      lims_after <- get_plot_limits(tx_after_local)
      
      p_after <- ggplot() +
        geom_path(
          data = poly_local,
          aes(x = x_um, y = y_um, group = cell),
          color = "grey55",
          linewidth = 0.35,
          alpha = 0.8
        ) +
        geom_point(
          data = tx_after_local,
          aes(x = x_um, y = y_um, color = plot_cellID),
          size = 1.4,
          alpha = 0.9
        ) +
        scale_color_manual(values = local_colors, drop = FALSE) +
        coord_fixed(xlim = lims_after$xlim, ylim = lims_after$ylim) +
        labs(
          title = paste("AFTER -", tc),
          x = "x (micron)",
          y = "y (micron)",
          color = "updated cell ID"
        ) +
        theme_minimal(base_size = 14) +
        theme(
          plot.title = element_text(size = 16),
          legend.position = "right"
        )
      
      ggsave(
        filename = file.path(out_dir, paste0("top_changed_after_", tc, ".png")),
        plot = p_after,
        width = 10,
        height = 8,
        dpi = 300
      )
    }
    
    # -------- 17) Clean summary + write to file --------
    updated_transDF[, changed := UMI_cellID != updated_cellID]
    
    total_transcripts <- nrow(updated_transDF)
    changed_transcripts_n <- sum(updated_transDF$changed)
    pct_changed <- changed_transcripts_n / total_transcripts * 100
    
    cell_changes <- updated_transDF[, .(
      total = .N,
      changed = sum(changed)
    ), by = UMI_cellID]
    
    cell_changes[, pct_changed := changed / total * 100]
    
    mean_pct <- mean(cell_changes$pct_changed)
    median_pct <- median(cell_changes$pct_changed)
    max_pct <- max(cell_changes$pct_changed)
    
    cells_gt5 <- sum(cell_changes$pct_changed > 5)
    cells_gt10 <- sum(cell_changes$pct_changed > 10)
    
    summary_text <- paste0(
      "FastReseg Results Summary (FOV ", test_fov, ")\n",
      "========================================\n\n",
      "Dataset size:\n",
      "- Cells: ", nrow(meta_fov), "\n",
      "- Transcripts: ", total_transcripts, "\n\n",
      "Resegmentation:\n",
      "- Transcripts changed: ", changed_transcripts_n, "\n",
      "- Percent changed: ", round(pct_changed, 4), "%\n\n",
      "Per-cell changes:\n",
      "- Mean % change: ", round(mean_pct, 4), "\n",
      "- Median % change: ", round(median_pct, 4), "\n",
      "- Max % change: ", round(max_pct, 4), "\n\n",
      "Cells affected:\n",
      "- >5% changed: ", cells_gt5, "\n",
      "- >10% changed: ", cells_gt10, "\n\n",
      "Interpretation:\n",
      if (pct_changed < 0.5 && max_pct < 10) {
        "Segmentation is mostly correct.\nOnly minor local corrections were needed.\n"
      } else if (pct_changed < 2) {
        "Moderate corrections detected.\nSome local improvements present.\n"
      } else {
        "Significant segmentation issues detected.\nResegmentation recommended.\n"
      }
    )
    
    summary_file <- file.path(out_dir, paste0("FastReseg_summary_FOV_", test_fov, ".txt"))
    writeLines(summary_text, summary_file)
    
    cat(summary_text)
    cat("\nSummary saved to:\n", summary_file, "\n")
    
    # -------- 18) Save transcript-level old -> new assignments --------
    changed_transcripts_dt <- copy(updated_transDF)[UMI_cellID != updated_cellID]
    
    keep_cols <- intersect(
      c(
        "target",
        "UMI_cellID",
        "updated_cellID",
        "updated_celltype",
        "score_updated_celltype",
        "x", "y", "z",
        "transComp"
      ),
      colnames(changed_transcripts_dt)
    )
    
    changed_transcripts_dt <- changed_transcripts_dt[, ..keep_cols]
    
    setnames(
      changed_transcripts_dt,
      old = c("UMI_cellID", "updated_cellID"),
      new = c("old_cell", "new_cell"),
      skip_absent = TRUE
    )
    
    changed_transcripts_dt[, fov := test_fov]
    setcolorder(changed_transcripts_dt, c("fov", setdiff(colnames(changed_transcripts_dt), "fov")))
    
    changed_transcripts_file <- file.path(
      out_dir,
      paste0("FastReseg_changed_transcripts_FOV_", test_fov, ".csv")
    )
    
    fwrite(changed_transcripts_dt, changed_transcripts_file)
    
    cat("Saved transcript-level changed assignments to:\n", changed_transcripts_file, "\n")
    cat("Number of changed transcripts:", nrow(changed_transcripts_dt), "\n\n")
    
    # -------- 19) Save cell-level old -> new cell name summary --------
    cell_rename_summary <- copy(updated_transDF)[
      UMI_cellID != updated_cellID,
      .(
        n_changed_transcripts = .N,
        n_genes_changed = uniqueN(target)
      ),
      by = .(old_cell = UMI_cellID, new_cell = updated_cellID)
    ]
    
    setorder(cell_rename_summary, -n_changed_transcripts, old_cell, new_cell)
    
    old_cell_totals <- updated_transDF[, .(old_cell_total_transcripts = .N), by = .(old_cell = UMI_cellID)]
    
    cell_rename_summary <- merge(
      cell_rename_summary,
      old_cell_totals,
      by = "old_cell",
      all.x = TRUE
    )
    
    cell_rename_summary[, pct_of_old_cell_transcripts := 100 * n_changed_transcripts / old_cell_total_transcripts]
    cell_rename_summary[, fov := test_fov]
    setcolorder(cell_rename_summary, c("fov", setdiff(colnames(cell_rename_summary), "fov")))
    setorder(cell_rename_summary, -n_changed_transcripts, -pct_of_old_cell_transcripts)
    
    cell_rename_file <- file.path(
      out_dir,
      paste0("FastReseg_cell_rename_summary_FOV_", test_fov, ".csv")
    )
    
    fwrite(cell_rename_summary, cell_rename_file)
    
    cat("Saved cell-level rename summary to:\n", cell_rename_file, "\n")
    cat("Number of old->new cell mappings:", nrow(cell_rename_summary), "\n\n")
    
    # Return objects needed for combined outputs
    list(
      status = "success",
      summary_metrics = data.table(
        fov = test_fov,
        n_cells = nrow(meta_fov),
        n_transcripts = total_transcripts,
        n_changed_transcripts = changed_transcripts_n,
        pct_changed = pct_changed,
        mean_pct_cell_change = mean_pct,
        median_pct_cell_change = median_pct,
        max_pct_cell_change = max_pct,
        cells_gt5pct = cells_gt5,
        cells_gt10pct = cells_gt10
      ),
      changed_transcripts_dt = changed_transcripts_dt,
      cell_rename_summary = cell_rename_summary
    )
    
  }, error = function(e) {
    err_msg <- conditionMessage(e)
    cat("FOV", test_fov, "FAILED:\n", err_msg, "\n\n")
    
    fail_file <- file.path(status_dir, paste0("FOV_", test_fov, "_FAILED.txt"))
    writeLines(
      c(
        paste("FOV", test_fov, "failed"),
        paste("Time:", as.character(Sys.time())),
        paste("Error:", err_msg)
      ),
      fail_file
    )
    
    list(
      status = "failed",
      error_message = err_msg
    )
  })
  
  # Save status row
  if (identical(fov_result$status, "success")) {
    run_status_list[[length(run_status_list) + 1]] <- data.table(
      fov = test_fov,
      status = "success",
      error_message = NA_character_
    )
    
    summary_metrics_list[[length(summary_metrics_list) + 1]] <- fov_result$summary_metrics
    all_changed_transcripts_list[[length(all_changed_transcripts_list) + 1]] <- fov_result$changed_transcripts_dt
    all_cell_rename_summary_list[[length(all_cell_rename_summary_list) + 1]] <- fov_result$cell_rename_summary
    
    success_file <- file.path(status_dir, paste0("FOV_", test_fov, "_SUCCESS.txt"))
    writeLines(
      c(
        paste("FOV", test_fov, "completed successfully"),
        paste("Time:", as.character(Sys.time()))
      ),
      success_file
    )
  } else {
    run_status_list[[length(run_status_list) + 1]] <- data.table(
      fov = test_fov,
      status = "failed",
      error_message = fov_result$error_message
    )
  }
  
  gc()
}

# -------------------------------
# 9) Combine and save all outputs
# -------------------------------
cat("\nCombining outputs across all FOVs...\n")

run_status_dt <- rbindlist(run_status_list, fill = TRUE)
summary_metrics_dt <- rbindlist(summary_metrics_list, fill = TRUE)
all_changed_transcripts_dt <- rbindlist(all_changed_transcripts_list, fill = TRUE)
all_cell_rename_summary_dt <- rbindlist(all_cell_rename_summary_list, fill = TRUE)

fwrite(run_status_dt, file.path(summary_dir, "all_fov_run_status.csv"))
fwrite(summary_metrics_dt, file.path(summary_dir, "all_fov_summary_metrics.csv"))
fwrite(all_changed_transcripts_dt, file.path(summary_dir, "all_fov_changed_transcripts.csv"))
fwrite(all_cell_rename_summary_dt, file.path(summary_dir, "all_fov_cell_rename_summary.csv"))

saveRDS(run_status_dt, file.path(summary_dir, "all_fov_run_status.rds"))
saveRDS(summary_metrics_dt, file.path(summary_dir, "all_fov_summary_metrics.rds"))
saveRDS(all_changed_transcripts_dt, file.path(summary_dir, "all_fov_changed_transcripts.rds"))
saveRDS(all_cell_rename_summary_dt, file.path(summary_dir, "all_fov_cell_rename_summary.rds"))

# -------------------------------
# 10) Final run report
# -------------------------------
n_success <- sum(run_status_dt$status == "success", na.rm = TRUE)
n_failed  <- sum(run_status_dt$status == "failed", na.rm = TRUE)

final_report <- c(
  "FastReseg all-FOV run complete",
  "================================",
  paste("Time:", as.character(Sys.time())),
  paste("Total FOVs attempted:", length(all_fovs)),
  paste("Successful FOVs:", n_success),
  paste("Failed FOVs:", n_failed),
  "",
  "Combined output files:",
  paste("-", file.path(summary_dir, "all_fov_run_status.csv")),
  paste("-", file.path(summary_dir, "all_fov_summary_metrics.csv")),
  paste("-", file.path(summary_dir, "all_fov_changed_transcripts.csv")),
  paste("-", file.path(summary_dir, "all_fov_cell_rename_summary.csv"))
)

writeLines(final_report, file.path(summary_dir, "FINAL_REPORT.txt"))

cat("\n========================================\n")
cat("ALL DONE\n")
cat("========================================\n")
cat("Successful FOVs:", n_success, "\n")
cat("Failed FOVs:", n_failed, "\n")
cat("Combined outputs saved in:\n", normalizePath(summary_dir), "\n")