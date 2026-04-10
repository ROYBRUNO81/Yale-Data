# ============================================================
# FASTRESEG ALL-FOV WORKFLOW USING RAW FILES ONLY
# ============================================================
# WHAT THIS SCRIPT DOES
# 1) Reads raw Yale files directly (no adata object)
# 2) Runs FastReseg FLAGGING on every FOV
# 3) Saves a dataset-wide flag summary (csv + txt + rds)
# 4) Suggests a threshold for choosing FOVs to refine
# 5) Runs FULL FastReseg only on selected FOVs
# 6) Saves per-FOV full objects, summaries, and plots
# 7) Saves a master object so you can revisit key FOVs later
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
expr_file <- "Raw Data/Colon_TMA_exprMat_file.csv.gz"
meta_file <- "Raw Data/Colon_TMA_metadata_file.csv.gz"
tx_file   <- "Raw Data/Colon_TMA_tx_file.csv.gz"
poly_file <- "Raw Data/Colon_TMA-polygons.csv.gz"
fovpos_file <- "Raw Data/Colon_TMA_fov_positions_file.csv.gz"

cluster_col <- "RNA_RNA.QC.Module_Cell.Typing.InSituType.1_1_clusters"

# -------------------------------
# 3) Define output folders
# -------------------------------
root_out <- "FastReseg_allFOV"
dir.create(root_out, recursive = TRUE, showWarnings = FALSE)

flag_out_dir    <- file.path(root_out, "01_flagging")
full_out_dir    <- file.path(root_out, "02_full_refinement")
plot_out_dir    <- file.path(root_out, "03_plots")
summary_out_dir <- file.path(root_out, "04_summaries")
object_out_dir  <- file.path(root_out, "05_objects")

dir.create(flag_out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(full_out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(plot_out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(summary_out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(object_out_dir, recursive = TRUE, showWarnings = FALSE)

cat("Root output folder:\n", normalizePath(root_out), "\n\n")

# -------------------------------
# 4) User-tunable settings
# -------------------------------
pixel_size_um <- 0.12028
zstep_size_um <- 0.8

flagModel_TransNum_cutoff <- 50
flagCell_lrtest_cutoff <- 2
svmClass_score_cutoff <- -2

percentCores_flag <- 0.25
percentCores_full <- 0.25

# If not NULL, overrides auto-selected FOVs
manual_selected_fovs <- NULL
max_refine_fovs <- 20

# Plot settings for refined FOVs
n_top_changed_cells_to_plot <- 3
n_average_changed_cells_to_plot <- 2
plot_radius_px <- 120
plot_pad_um <- 5

# -------------------------------
# 5) Discover all FOVs
# -------------------------------
fov_positions <- fread(fovpos_file)
all_fovs <- sort(unique(as.integer(fov_positions$FOV)))

cat("Number of FOVs detected:", length(all_fovs), "\n")
cat("First 10 FOVs:", paste(head(all_fovs, 10), collapse = ", "), "\n\n")

# -------------------------------
# 6) Helper functions
# -------------------------------

read_meta_fov <- function(fov_id) {
  meta_fov <- fread(
    cmd = sprintf(
      "gzip -dc %s | awk -F',' 'NR==1 || $3==%d'",
      shQuote(meta_file), fov_id
    )
  )
  if (!"cell" %in% colnames(meta_fov)) {
    meta_fov[, cell := paste0("c_1_", fov, "_", cell_ID)]
  }
  meta_fov
}

read_expr_fov <- function(fov_id) {
  expr_fov <- fread(
    cmd = sprintf(
      "gzip -dc %s | awk -F',' 'NR==1 || $1==%d'",
      shQuote(expr_file), fov_id
    )
  )
  expr_fov[, cell := paste0("c_1_", fov, "_", cell_ID)]
  expr_fov
}

read_tx_fov <- function(fov_id) {
  fread(
    cmd = sprintf(
      "gzip -dc %s | awk -F',' 'NR==1 || $1==%d'",
      shQuote(tx_file), fov_id
    )
  )
}

read_poly_fov <- function(fov_id) {
  fread(
    cmd = sprintf(
      "gzip -dc %s | awk -F',' 'NR==1 || $1==%d'",
      shQuote(poly_file), fov_id
    )
  )
}

build_counts_and_clusters <- function(meta_fov, expr_fov, tx_fov, cluster_col, extracellular_cellID) {
  if (nrow(meta_fov) == 0) stop("meta_fov has 0 rows")
  if (nrow(expr_fov) == 0) stop("expr_fov has 0 rows")
  if (nrow(tx_fov) == 0) stop("tx_fov has 0 rows")
  if (!(cluster_col %in% colnames(meta_fov))) stop("cluster_col not found in meta_fov")
  
  common_cells <- intersect(expr_fov$cell, meta_fov$cell)
  if (length(common_cells) == 0) stop("No common cells between expr_fov and meta_fov")
  
  expr_fov <- expr_fov[cell %in% common_cells]
  meta_fov <- meta_fov[cell %in% common_cells]
  
  setkey(expr_fov, cell)
  setkey(meta_fov, cell)
  expr_fov <- expr_fov[meta_fov$cell]
  
  id_cols <- c("fov", "cell_ID", "cell")
  gene_cols <- setdiff(colnames(expr_fov), id_cols)
  if (length(gene_cols) == 0) stop("No gene columns found in expr_fov")
  
  counts_fov <- as.matrix(expr_fov[, ..gene_cols])
  rownames(counts_fov) <- expr_fov$cell
  storage.mode(counts_fov) <- "numeric"
  
  clust_fov <- as.character(meta_fov[[cluster_col]])
  names(clust_fov) <- meta_fov$cell
  
  genes_tx <- unique(tx_fov$target)
  common_genes <- intersect(colnames(counts_fov), genes_tx)
  if (length(common_genes) == 0) stop("No common genes between counts_fov and tx_fov")
  
  counts_fov <- counts_fov[, common_genes, drop = FALSE]
  
  tx_cells_assigned <- unique(tx_fov$cell)
  tx_cells_assigned <- setdiff(tx_cells_assigned, extracellular_cellID)
  
  pilot_cells <- intersect(rownames(counts_fov), tx_cells_assigned)
  if (length(pilot_cells) == 0) stop("No usable assigned cells remain after transcript filtering")
  
  counts_fov <- counts_fov[pilot_cells, , drop = FALSE]
  clust_fov <- clust_fov[pilot_cells]
  
  if (is.null(dim(counts_fov))) stop("counts_fov lost matrix dimensions")
  if (nrow(counts_fov) < 2) stop("counts_fov has fewer than 2 cells")
  if (ncol(counts_fov) < 2) stop("counts_fov has fewer than 2 genes")
  if (length(clust_fov) != nrow(counts_fov)) stop("length(clust_fov) does not match nrow(counts_fov)")
  if (length(clust_fov) == 0) stop("clust_fov is empty")
  
  counts_fov <- Matrix::Matrix(counts_fov, sparse = TRUE)
  
  list(
    meta_fov = meta_fov,
    expr_fov = expr_fov,
    counts_fov = counts_fov,
    clust_fov = clust_fov,
    common_genes = common_genes,
    pilot_cells = pilot_cells
  )
}

extract_flagged_cells <- function(res_obj) {
  flagged_cells <- character(0)
  
  if ("combined_flaggedCells" %in% names(res_obj)) {
    x <- res_obj$combined_flaggedCells
    
    if (is.null(x)) {
      flagged_cells <- character(0)
    } else if (is.vector(x) || is.factor(x)) {
      flagged_cells <- as.character(x)
    } else if (is.list(x) && length(x) >= 1) {
      flagged_cells <- unique(as.character(unlist(x)))
    } else if (is.data.frame(x)) {
      possible_cols <- c("cell", "cellID", "cell_id", "cell_ID", "UMI_cellID")
      hit <- intersect(possible_cols, colnames(x))
      if (length(hit) > 0) {
        flagged_cells <- unique(as.character(x[[hit[1]]]))
      }
    }
  }
  
  unique(flagged_cells)
}

safe_nrow <- function(x) {
  if (is.null(x)) return(NA_integer_)
  if (is.data.frame(x) || is.matrix(x)) return(nrow(x))
  NA_integer_
}

make_local_color_map <- function(before_ids, after_ids) {
  before_ids <- sort(unique(as.character(before_ids)))
  after_ids  <- sort(unique(as.character(after_ids)))
  
  shared_ids  <- intersect(before_ids, after_ids)
  before_only <- setdiff(before_ids, after_ids)
  after_only  <- setdiff(after_ids, before_ids)
  
  ordered_ids <- c(shared_ids, before_only, after_only)
  n_ids <- length(ordered_ids)
  
  cols <- hue_pal()(max(n_ids, 1))
  names(cols) <- ordered_ids
  cols
}

get_plot_limits <- function(df, pad = 5) {
  list(
    xlim = range(df$x_um, na.rm = TRUE) + c(-pad, pad),
    ylim = range(df$y_um, na.rm = TRUE) + c(-pad, pad)
  )
}

save_local_before_after_plots <- function(fov_id, meta_fov, poly_fov, updated_transDF,
                                          out_dir,
                                          n_top = 3,
                                          n_avg = 2,
                                          radius_px = 120,
                                          pixel_size_um = 0.12028,
                                          pad_um = 5) {
  
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  
  updated_transDF[, changed := UMI_cellID != updated_cellID]
  
  cell_changes <- updated_transDF[, .(
    total_transcripts = .N,
    changed_transcripts = sum(changed)
  ), by = UMI_cellID]
  cell_changes[, pct_changed := 100 * changed_transcripts / total_transcripts]
  
  setorder(cell_changes, -changed_transcripts, -pct_changed)
  
  top_cells <- head(cell_changes$UMI_cellID, n_top)
  
  nonzero_cells <- cell_changes[pct_changed > 0]
  avg_cells <- character(0)
  if (nrow(nonzero_cells) > 0) {
    target_pct <- median(nonzero_cells$pct_changed)
    nonzero_cells[, dist_to_median := abs(pct_changed - target_pct)]
    setorder(nonzero_cells, dist_to_median, -changed_transcripts)
    avg_cells <- head(setdiff(nonzero_cells$UMI_cellID, top_cells), n_avg)
  }
  
  selected_cells <- unique(c(top_cells, avg_cells))
  
  poly_fov <- copy(poly_fov)
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
  
  selected_info <- cell_changes[UMI_cellID %in% selected_cells]
  fwrite(selected_info, file.path(out_dir, sprintf("FOV_%03d_selected_cells_for_plots.csv", fov_id)))
  
  for (tc in selected_cells) {
    target_row <- meta_fov[cell == tc]
    if (nrow(target_row) == 0) next
    
    cx <- target_row$CenterX_local_px[1]
    cy <- target_row$CenterY_local_px[1]
    
    meta_local <- copy(meta_fov)
    meta_local[, dist_px := sqrt((CenterX_local_px - cx)^2 + (CenterY_local_px - cy)^2)]
    near_cells <- meta_local[dist_px <= radius_px, cell]
    
    poly_local <- poly_fov[cell %in% near_cells]
    tx_before_local <- tx_before[UMI_cellID %in% near_cells | updated_cellID %in% near_cells]
    tx_after_local  <- tx_after[UMI_cellID %in% near_cells | updated_cellID %in% near_cells]
    
    local_colors <- make_local_color_map(
      before_ids = tx_before_local$plot_cellID,
      after_ids  = tx_after_local$plot_cellID
    )
    
    lims_before <- get_plot_limits(tx_before_local, pad = pad_um)
    lims_after  <- get_plot_limits(tx_after_local,  pad = pad_um)
    
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
        title = paste("BEFORE -", tc, sprintf("(FOV %d)", fov_id)),
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
      filename = file.path(out_dir, sprintf("FOV_%03d_BEFORE_%s.png", fov_id, tc)),
      plot = p_before,
      width = 10,
      height = 8,
      dpi = 300
    )
    
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
        title = paste("AFTER -", tc, sprintf("(FOV %d)", fov_id)),
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
      filename = file.path(out_dir, sprintf("FOV_%03d_AFTER_%s.png", fov_id, tc)),
      plot = p_after,
      width = 10,
      height = 8,
      dpi = 300
    )
  }
  
  invisible(selected_cells)
}

# -------------------------------
# 7) ALL-FOV FLAGGING PASS
# -------------------------------
flag_summary_list <- vector("list", length(all_fovs))
names(flag_summary_list) <- as.character(all_fovs)

cat("Starting all-FOV flagging pass...\n\n")

for (i in seq_along(all_fovs)) {
  fov_id <- all_fovs[i]
  cat(sprintf("[%d/%d] Flagging FOV %d\n", i, length(all_fovs), fov_id))
  
  result_row <- tryCatch({
    extracellular_cellID <- paste0("c_1_", fov_id, "_0")
    
    meta_fov <- read_meta_fov(fov_id)
    expr_fov <- read_expr_fov(fov_id)
    tx_fov   <- read_tx_fov(fov_id)
    
    prep <- build_counts_and_clusters(
      meta_fov = meta_fov,
      expr_fov = expr_fov,
      tx_fov = tx_fov,
      cluster_col = cluster_col,
      extracellular_cellID = extracellular_cellID
    )
    
    counts_fov <- prep$counts_fov
    clust_fov  <- prep$clust_fov
    
    fov_flag_dir <- file.path(flag_out_dir, sprintf("FOV_%03d", fov_id))
    dir.create(fov_flag_dir, recursive = TRUE, showWarnings = FALSE)
    
    res_flag <- fastReseg_flag_all_errors(
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
      
      flagModel_TransNum_cutoff = flagModel_TransNum_cutoff,
      flagCell_lrtest_cutoff = flagCell_lrtest_cutoff,
      svmClass_score_cutoff = svmClass_score_cutoff,
      
      path_to_output = fov_flag_dir,
      transDF_export_option = 0,
      return_trimmed_perCell = FALSE,
      combine_extra = FALSE,
      
      percentCores = percentCores_flag,
      seed_transError = 1
    )
    
    flagged_cells <- extract_flagged_cells(res_flag)
    n_flagged <- length(flagged_cells)
    n_eval <- safe_nrow(res_flag$combined_modStats_ToFlagCells)
    if (is.na(n_eval)) n_eval <- length(prep$pilot_cells)
    
    data.table(
      fov = fov_id,
      status = "ok",
      error_message = NA_character_,
      cells_total = nrow(meta_fov),
      cells_evaluated = n_eval,
      flagged_cells = n_flagged,
      pct_flagged = ifelse(n_eval > 0, 100 * n_flagged / n_eval, NA_real_),
      transcripts_total = nrow(tx_fov),
      extracellular_transcripts = sum(tx_fov$cell == extracellular_cellID, na.rm = TRUE),
      pct_extracellular = 100 * sum(tx_fov$cell == extracellular_cellID, na.rm = TRUE) / max(nrow(tx_fov), 1),
      common_genes = length(prep$common_genes),
      cells_used = length(prep$pilot_cells)
    )
    
  }, error = function(e) {
    data.table(
      fov = fov_id,
      status = "failed",
      error_message = as.character(e$message),
      cells_total = NA_integer_,
      cells_evaluated = NA_integer_,
      flagged_cells = NA_integer_,
      pct_flagged = NA_real_,
      transcripts_total = NA_integer_,
      extracellular_transcripts = NA_integer_,
      pct_extracellular = NA_real_,
      common_genes = NA_integer_,
      cells_used = NA_integer_
    )
  })
  
  flag_summary_list[[as.character(fov_id)]] <- result_row
  gc()
}

flag_summary <- rbindlist(flag_summary_list, fill = TRUE)
setorder(flag_summary, -pct_flagged, -flagged_cells)

fwrite(flag_summary, file.path(summary_out_dir, "allFOV_flag_summary.csv"))
saveRDS(flag_summary, file.path(object_out_dir, "allFOV_flag_summary.rds"))

# -------------------------------
# 8) AUTO-SELECT FOVS TO REFINE
# -------------------------------
ok_flag_summary <- flag_summary[status == "ok" & !is.na(pct_flagged)]

pct_threshold_auto <- as.numeric(quantile(ok_flag_summary$pct_flagged, 0.90, na.rm = TRUE))
nflag_threshold_auto <- max(
  20L,
  as.integer(quantile(ok_flag_summary$flagged_cells, 0.75, na.rm = TRUE))
)

selected_fovs_auto <- ok_flag_summary[
  pct_flagged >= pct_threshold_auto &
    flagged_cells >= nflag_threshold_auto
][order(-pct_flagged, -flagged_cells)]$fov

selected_fovs_auto <- head(selected_fovs_auto, max_refine_fovs)

if (!is.null(manual_selected_fovs)) {
  selected_fovs <- sort(unique(as.integer(manual_selected_fovs)))
  selection_mode <- "manual"
} else {
  selected_fovs <- selected_fovs_auto
  selection_mode <- "automatic"
}

cat("\nFlagging pass complete.\n")
cat("Suggested refinement threshold based on flag burden:\n")
cat("- pct_flagged >=", round(pct_threshold_auto, 2), "%\n")
cat("- flagged_cells >=", nflag_threshold_auto, "\n")
cat("Selection mode:", selection_mode, "\n")
cat("Selected FOVs for full refinement:\n")
print(selected_fovs)
cat("\n")

flag_summary_txt <- c(
  "FastReseg all-FOV flagging summary",
  "================================",
  "",
  sprintf("Total FOVs processed: %d", length(all_fovs)),
  sprintf("Successful FOVs: %d", nrow(flag_summary[status == 'ok'])),
  sprintf("Failed FOVs: %d", nrow(flag_summary[status == 'failed'])),
  sprintf("Automatic review threshold for pct_flagged: %.2f%%", pct_threshold_auto),
  sprintf("Automatic review threshold for flagged_cells: %d", nflag_threshold_auto),
  sprintf("Selection mode: %s", selection_mode),
  sprintf(
    "Selected FOVs for full refinement (%d): %s",
    length(selected_fovs),
    ifelse(length(selected_fovs) > 0, paste(selected_fovs, collapse = ", "), "none")
  ),
  "",
  "Top 15 successful FOVs by pct_flagged:",
  paste(
    apply(
      as.matrix(head(ok_flag_summary[order(-pct_flagged, -flagged_cells), .(fov, pct_flagged, flagged_cells, cells_evaluated)], 15)),
      1,
      function(x) sprintf("FOV %s: pct_flagged=%s, flagged_cells=%s, cells_evaluated=%s", x[1], x[2], x[3], x[4])
    ),
    collapse = "\n"
  ),
  "",
  "Failed FOVs:",
  if (nrow(flag_summary[status == "failed"]) > 0) {
    paste(
      apply(
        as.matrix(flag_summary[status == "failed", .(fov, error_message)]),
        1,
        function(x) sprintf("FOV %s failed: %s", x[1], x[2])
      ),
      collapse = "\n"
    )
  } else {
    "None"
  }
)
writeLines(flag_summary_txt, file.path(summary_out_dir, "allFOV_flagging_summary.txt"))

# -------------------------------
# 9) RUN FULL REFINEMENT ONLY ON SELECTED FOVS
# -------------------------------
refine_summary_list <- list()

if (length(selected_fovs) == 0) {
  cat("No FOVs selected for full refinement under current heuristic.\n")
} else {
  cat("Starting full refinement on selected FOVs...\n\n")
  
  for (j in seq_along(selected_fovs)) {
    fov_id <- selected_fovs[j]
    cat(sprintf("[%d/%d] Full refinement for FOV %d\n", j, length(selected_fovs), fov_id))
    
    tryCatch({
      extracellular_cellID <- paste0("c_1_", fov_id, "_0")
      
      meta_fov <- read_meta_fov(fov_id)
      expr_fov <- read_expr_fov(fov_id)
      tx_fov   <- read_tx_fov(fov_id)
      poly_fov <- read_poly_fov(fov_id)
      
      prep <- build_counts_and_clusters(
        meta_fov = meta_fov,
        expr_fov = expr_fov,
        tx_fov = tx_fov,
        cluster_col = cluster_col,
        extracellular_cellID = extracellular_cellID
      )
      
      counts_fov <- prep$counts_fov
      clust_fov  <- prep$clust_fov
      
      fov_full_dir <- file.path(full_out_dir, sprintf("FOV_%03d", fov_id))
      dir.create(fov_full_dir, recursive = TRUE, showWarnings = FALSE)
      
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
        
        flagModel_TransNum_cutoff = flagModel_TransNum_cutoff,
        flagCell_lrtest_cutoff = flagCell_lrtest_cutoff,
        svmClass_score_cutoff = svmClass_score_cutoff,
        
        molecular_distance_cutoff = NULL,
        cellular_distance_cutoff = NULL,
        score_baseline = NULL,
        lowerCutoff_transNum = NULL,
        higherCutoff_transNum = NULL,
        
        imputeFlag_missingCTs = TRUE,
        groupTranscripts_method = "dbscan",
        spatialMergeCheck_method = "leidenCut",
        cutoff_spatialMerge = 0.5,
        
        path_to_output = fov_full_dir,
        transDF_export_option = 2,
        save_intermediates = TRUE,
        return_perCellData = TRUE,
        combine_extra = FALSE,
        
        seed_process = 1,
        percentCores = percentCores_full
      )
      
      updated_transDF <- NULL
      if ("updated_transDF_list" %in% names(res_full) && length(res_full$updated_transDF_list) >= 1) {
        updated_transDF <- as.data.table(res_full$updated_transDF_list[[1]])
      }
      
      if (!is.null(updated_transDF)) {
        updated_transDF[, changed := UMI_cellID != updated_cellID]
        
        total_transcripts <- nrow(updated_transDF)
        changed_transcripts <- sum(updated_transDF$changed)
        pct_changed <- 100 * changed_transcripts / max(total_transcripts, 1)
        
        cell_changes <- updated_transDF[, .(
          total_transcripts = .N,
          changed_transcripts = sum(changed)
        ), by = UMI_cellID]
        cell_changes[, pct_changed := 100 * changed_transcripts / total_transcripts]
        
        mean_pct <- mean(cell_changes$pct_changed)
        median_pct <- median(cell_changes$pct_changed)
        max_pct <- max(cell_changes$pct_changed)
        cells_gt5 <- sum(cell_changes$pct_changed > 5)
        cells_gt10 <- sum(cell_changes$pct_changed > 10)
        
        fwrite(cell_changes, file.path(fov_full_dir, sprintf("FOV_%03d_perCell_change_summary.csv", fov_id)))
        
        fov_plot_dir <- file.path(plot_out_dir, sprintf("FOV_%03d", fov_id))
        selected_cells_for_plots <- save_local_before_after_plots(
          fov_id = fov_id,
          meta_fov = meta_fov,
          poly_fov = poly_fov,
          updated_transDF = updated_transDF,
          out_dir = fov_plot_dir,
          n_top = n_top_changed_cells_to_plot,
          n_avg = n_average_changed_cells_to_plot,
          radius_px = plot_radius_px,
          pixel_size_um = pixel_size_um,
          pad_um = plot_pad_um
        )
        
        interpretation_line <- if (pct_changed < 0.5 && max_pct < 10) {
          "Interpretation: segmentation is mostly correct; only minor local corrections were needed."
        } else if (pct_changed < 2) {
          "Interpretation: moderate corrections detected; some local improvements are present."
        } else {
          "Interpretation: stronger segmentation issues detected; resegmentation appears worthwhile."
        }
        
        fov_summary_text <- c(
          sprintf("FastReseg Results Summary (FOV %d)", fov_id),
          "========================================",
          "",
          "Dataset size:",
          sprintf("- Cells: %d", nrow(meta_fov)),
          sprintf("- Transcripts: %d", total_transcripts),
          sprintf(
            "- Extracellular transcripts removed upstream: %d (%.3f%% of transcript rows in this FOV)",
            sum(tx_fov$cell == extracellular_cellID, na.rm = TRUE),
            100 * sum(tx_fov$cell == extracellular_cellID, na.rm = TRUE) / max(nrow(tx_fov), 1)
          ),
          "",
          "Resegmentation:",
          sprintf("- Transcripts changed: %d", changed_transcripts),
          sprintf("- Percent changed: %.4f%%", pct_changed),
          "",
          "Per-cell changes:",
          sprintf("- Mean %% change: %.4f", mean_pct),
          sprintf("- Median %% change: %.4f", median_pct),
          sprintf("- Max %% change: %.4f", max_pct),
          "",
          "Cells affected:",
          sprintf("- >5%% changed: %d", cells_gt5),
          sprintf("- >10%% changed: %d", cells_gt10),
          "",
          interpretation_line,
          "",
          sprintf(
            "Cells plotted for review: %s",
            ifelse(length(selected_cells_for_plots) > 0, paste(selected_cells_for_plots, collapse = ", "), "none")
          )
        )
        
        writeLines(
          fov_summary_text,
          file.path(fov_full_dir, sprintf("FOV_%03d_FastReseg_summary.txt", fov_id))
        )
        
        refine_summary_list[[as.character(fov_id)]] <- data.table(
          fov = fov_id,
          status = "ok",
          error_message = NA_character_,
          cells_total = nrow(meta_fov),
          transcripts_total = total_transcripts,
          transcripts_changed = changed_transcripts,
          pct_changed = pct_changed,
          mean_pct_change_per_cell = mean_pct,
          median_pct_change_per_cell = median_pct,
          max_pct_change_per_cell = max_pct,
          cells_gt5pct = cells_gt5,
          cells_gt10pct = cells_gt10,
          selected_for_refinement = TRUE
        )
        
      } else {
        refine_summary_list[[as.character(fov_id)]] <- data.table(
          fov = fov_id,
          status = "failed",
          error_message = "updated_transDF not found in returned object",
          cells_total = nrow(meta_fov),
          transcripts_total = NA_integer_,
          transcripts_changed = NA_integer_,
          pct_changed = NA_real_,
          mean_pct_change_per_cell = NA_real_,
          median_pct_change_per_cell = NA_real_,
          max_pct_change_per_cell = NA_real_,
          cells_gt5pct = NA_integer_,
          cells_gt10pct = NA_integer_,
          selected_for_refinement = TRUE
        )
      }
      
      saveRDS(
        res_full,
        file.path(object_out_dir, sprintf("FOV_%03d_FastReseg_full_result.rds", fov_id))
      )
      
      rm(meta_fov, expr_fov, tx_fov, poly_fov, prep, counts_fov, clust_fov, res_full, updated_transDF)
      gc()
      
    }, error = function(e) {
      refine_summary_list[[as.character(fov_id)]] <<- data.table(
        fov = fov_id,
        status = "failed",
        error_message = as.character(e$message),
        cells_total = NA_integer_,
        transcripts_total = NA_integer_,
        transcripts_changed = NA_integer_,
        pct_changed = NA_real_,
        mean_pct_change_per_cell = NA_real_,
        median_pct_change_per_cell = NA_real_,
        max_pct_change_per_cell = NA_real_,
        cells_gt5pct = NA_integer_,
        cells_gt10pct = NA_integer_,
        selected_for_refinement = TRUE
      )
    })
  }
}

# -------------------------------
# 10) SAVE MASTER OUTPUTS
# -------------------------------
if (length(refine_summary_list) > 0) {
  refine_summary <- rbindlist(refine_summary_list, fill = TRUE)
  setorder(refine_summary, -pct_changed)
  fwrite(refine_summary, file.path(summary_out_dir, "selectedFOV_full_refinement_summary.csv"))
  saveRDS(refine_summary, file.path(object_out_dir, "selectedFOV_full_refinement_summary.rds"))
} else {
  refine_summary <- data.table()
}

master_obj <- list(
  all_fovs = all_fovs,
  selected_fovs = selected_fovs,
  flag_summary = flag_summary,
  refine_summary = refine_summary,
  settings = list(
    pixel_size_um = pixel_size_um,
    zstep_size_um = zstep_size_um,
    flagModel_TransNum_cutoff = flagModel_TransNum_cutoff,
    flagCell_lrtest_cutoff = flagCell_lrtest_cutoff,
    svmClass_score_cutoff = svmClass_score_cutoff,
    percentCores_flag = percentCores_flag,
    percentCores_full = percentCores_full,
    pct_threshold_auto = pct_threshold_auto,
    nflag_threshold_auto = nflag_threshold_auto,
    selection_mode = selection_mode,
    manual_selected_fovs = manual_selected_fovs,
    max_refine_fovs = max_refine_fovs
  )
)

saveRDS(master_obj, file.path(object_out_dir, "FastReseg_master_object.rds"))

final_summary_txt <- c(
  "FastReseg Yale Colon all-FOV workflow complete",
  "==============================================",
  "",
  sprintf("Total FOVs processed in flagging: %d", length(all_fovs)),
  sprintf("Successful FOVs in flagging: %d", nrow(flag_summary[status == 'ok'])),
  sprintf("Failed FOVs in flagging: %d", nrow(flag_summary[status == 'failed'])),
  sprintf("FOVs selected for full refinement: %d", length(selected_fovs)),
  sprintf("Selected FOV list: %s", ifelse(length(selected_fovs) > 0, paste(selected_fovs, collapse = ", "), "none")),
  "",
  sprintf("- Flag summary csv: %s", file.path(summary_out_dir, "allFOV_flag_summary.csv")),
  sprintf("- Flag summary txt: %s", file.path(summary_out_dir, "allFOV_flagging_summary.txt")),
  sprintf("- Full refinement summary csv: %s", file.path(summary_out_dir, "selectedFOV_full_refinement_summary.csv")),
  sprintf("- Master object rds: %s", file.path(object_out_dir, "FastReseg_master_object.rds"))
)
writeLines(final_summary_txt, file.path(summary_out_dir, "RUN_COMPLETE.txt"))

cat("\nWorkflow finished.\n")
cat("Main summary files are in:\n", normalizePath(summary_out_dir), "\n")
cat("Main objects are in:\n", normalizePath(object_out_dir), "\n")