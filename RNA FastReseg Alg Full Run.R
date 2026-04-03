# ============================================================
# FASTRESEG DATASET-WIDE SCREEN + OPTIONAL TARGETED REFINEMENT
# ============================================================

# ----------------------------
# 1) Load packages
# ----------------------------
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(scales)
  library(zellkonverter)
  library(SingleCellExperiment)
  library(DelayedArray)
  library(FastReseg)
})

# ----------------------------
# 2) User configuration
# ----------------------------
# Core inputs
adata_file         <- "colon_adata_clustered.h5ad"              # H5AD with cell x gene clustering object
tx_all_file        <- "Raw Data/Colon_TMA_tx_file.csv.gz"       # raw transcript file
meta_file          <- "Raw Data/Colon_TMA_metadata_file.csv.gz" # raw metadata
fov_pos_file       <- "Raw Data/Colon_TMA_fov_positions_file.csv.gz"

# Output root
run_root           <- "FastReseg_allFOV_screen"

# Cluster column to use from H5AD
cluster_col        <- "RNA_RNA.QC.Module_Cell.Typing.InSituType.1_1_clusters"

# Spatial scaling
pixel_size_um      <- 0.12028
zstep_size_um      <- 0.8

# Screening parameters (conservative first pass)
# FastReseg docs list flagCell_lrtest_cutoff default as 5 and svmClass_score_cutoff default as -2.
flagModel_TransNum_cutoff <- 50
flagCell_lrtest_cutoff    <- 5
svmClass_score_cutoff     <- -2

# Parallelism
percentCores_screen <- 0.50
percentCores_full   <- 0.25

# Thresholding strategy after screening
# If NULL, code will use 90th percentile of % flagged and count flagged.
refine_threshold_pctflag <- NULL
refine_threshold_nflag   <- NULL

# Full refinement controls
run_full_refine    <- FALSE     # FIRST RUN: leave FALSE. After reviewing summary txt, set TRUE if desired.
max_refine_fovs    <- 10        # cap number of FOVs to refine in one run
top_changed_n      <- 3         # top-changed cells to save plots for each refined FOV
avg_changed_n      <- 2         # average-like changed cells to save plots for each refined FOV
neighbor_radius_px <- 120

# Save objects
master_rds_file    <- file.path(run_root, "FastReseg_master_screen_object.rds")

# ----------------------------
# 3) Create output folders
# ----------------------------
dir.create(run_root, recursive = TRUE, showWarnings = FALSE)

tx_split_dir       <- file.path(run_root, "transcripts_by_fov")
flag_dir           <- file.path(run_root, "flagging_outputs")
summary_dir        <- file.path(run_root, "summary_outputs")
refine_dir         <- file.path(run_root, "refinement_outputs")
plot_dir           <- file.path(run_root, "saved_plots")

dir.create(tx_split_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(flag_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(summary_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(refine_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

cat("Run root:", normalizePath(run_root), "\n")

# ----------------------------
# 4) Read H5AD counts + clusters using disk-backed HDF5 mode
#    This is safer than materializing the full raw expression CSV in RAM.
# ----------------------------
cat("\n[Phase 0] Reading H5AD counts/clusters in HDF5-backed mode ...\n")

sce <- readH5AD(
  adata_file,
  use_hdf5 = TRUE,
  obs = c("fov", cluster_col),
  uns = FALSE,
  obsm = FALSE,
  varm = FALSE,
  obsp = FALSE,
  varp = FALSE,
  raw = FALSE
)

# H5AD assays are gene x cell; FastReseg expects cells x genes.
# DelayedArray transpose keeps this disk-backed.
x <- assay(sce, "X")
counts <- DelayedArray::t(x)

if (is.null(rownames(counts))) rownames(counts) <- colnames(sce)
if (is.null(colnames(counts))) colnames(counts) <- rownames(sce)

clust <- as.character(colData(sce)[[cluster_col]])
names(clust) <- colnames(sce)

fov_vec <- as.integer(as.character(colData(sce)$fov))
cell_counts_per_fov <- as.data.table(table(fov_vec))
setnames(cell_counts_per_fov, c("fov", "total_cells"))
cell_counts_per_fov[, fov := as.integer(as.character(fov))]

cat("H5AD loaded.\n")
cat("Counts dimensions (cells x genes):", nrow(counts), "x", ncol(counts), "\n")
cat("Unique clusters:", length(unique(clust)), "\n")
cat("Unique FOVs in H5AD:", length(unique(fov_vec)), "\n")

# ----------------------------
# 5) Split the giant transcript file into per-FOV CSVs (one-time step)
#    FastReseg multi-FOV wrappers are designed around transDF_fileInfo.
# ----------------------------
cat("\n[Phase 0] Preparing per-FOV transcript files ...\n")

existing_split_files <- list.files(tx_split_dir, pattern = "^FOV_[0-9]{3}_tx\\.csv$", full.names = TRUE)

if (length(existing_split_files) < 300) {
  cat("Per-FOV transcript files not found or incomplete. Splitting the raw transcript file now ...\n")
  
  # One-pass split of the raw gz file into per-FOV CSVs.
  # Requires a Unix-like shell (macOS/Linux).
  split_cmd <- sprintf(
    "mkdir -p %s && gzip -dc %s | awk -F',' 'NR==1{header=$0; next} {outfile=sprintf(\"%s/FOV_%%03d_tx.csv\", $1); if (!(outfile in seen)) {print header > outfile; seen[outfile]=1} print >> outfile}'",
    shQuote(tx_split_dir),
    shQuote(tx_all_file),
    tx_split_dir
  )
  
  status <- system(split_cmd)
  if (status != 0) stop("Transcript split command failed. Stop and inspect shell availability / file paths.")
} else {
  cat("Per-FOV transcript files already present. Reusing them.\n")
}

split_files <- list.files(tx_split_dir, pattern = "^FOV_[0-9]{3}_tx\\.csv$", full.names = TRUE)
cat("Per-FOV transcript files available:", length(split_files), "\n")

# ----------------------------
# 6) Read FOV positions and create transDF_fileInfo
#    FastReseg expects stage offsets in MICRONS for multi-FOV stitching.
# ----------------------------
cat("\n[Phase 0] Building transDF_fileInfo ...\n")

fov_pos <- fread(fov_pos_file)
setnames(fov_pos, old = c("FOV"), new = c("fov"))

fov_pos[, stage_X := x_global_mm * 1000]
fov_pos[, stage_Y := y_global_mm * 1000]
fov_pos[, slide := 1L]

file_info <- data.table(
  file_path = split_files,
  fov = as.integer(sub("^FOV_([0-9]{3})_tx\\.csv$", "\\1", basename(split_files)))
)

file_info <- merge(file_info, fov_pos[, .(fov, stage_X, stage_Y, slide)], by = "fov", all.x = TRUE)
setorder(file_info, fov)

# Keep only FOVs that exist in the H5AD object too
file_info <- file_info[fov %in% unique(fov_vec)]

cat("transDF_fileInfo rows:", nrow(file_info), "\n")
cat("FOV range:", min(file_info$fov), "to", max(file_info$fov), "\n")

# ----------------------------
# 7) Run dataset-wide flagging across all FOVs
# ----------------------------
cat("\n[Phase 1] Running FastReseg flagging across all FOVs ...\n")
cat("This is the multi-FOV screening pass.\n")
cat("Conservative screening settings:\n")
cat("- flagCell_lrtest_cutoff =", flagCell_lrtest_cutoff, "\n")
cat("- svmClass_score_cutoff  =", svmClass_score_cutoff, "\n")
cat("- flagModel_TransNum_cutoff =", flagModel_TransNum_cutoff, "\n\n")

set.seed(1)

res_flag_all <- fastReseg_flag_all_errors(
  counts = counts,
  clust = clust,
  refProfiles = NULL,
  
  transDF_fileInfo = as.data.frame(file_info),
  filepath_coln = "file_path",
  prefix_colns = c("slide", "fov"),
  fovOffset_colns = c("stage_X", "stage_Y"),
  
  pixel_size = pixel_size_um,
  zstep_size = zstep_size_um,
  
  transcript_df = NULL,
  transID_coln = NULL,
  transGene_coln = "target",
  cellID_coln = "cell_ID",
  spatLocs_colns = c("x_local_px", "y_local_px", "z"),
  invert_y = TRUE,
  
  extracellular_cellID = c(0),
  
  flagModel_TransNum_cutoff = flagModel_TransNum_cutoff,
  flagCell_lrtest_cutoff = flagCell_lrtest_cutoff,
  svmClass_score_cutoff = svmClass_score_cutoff,
  
  path_to_output = flag_dir,
  transDF_export_option = 2,
  return_trimmed_perCell = FALSE,
  combine_extra = FALSE,
  
  seed_transError = 1,
  percentCores = percentCores_screen
)

# ----------------------------
# 8) Summarize flagging across all FOVs
# ----------------------------
cat("\n[Phase 1] Summarizing FOV-level flagging ...\n")

flagged_list <- res_flag_all$combined_flaggedCells
flagged_counts <- data.table(
  fov = file_info$fov,
  flagged_cells = vapply(flagged_list, length, integer(1))
)

flag_summary <- merge(cell_counts_per_fov, flagged_counts, by = "fov", all.x = TRUE)
flag_summary[is.na(flagged_cells), flagged_cells := 0L]
flag_summary[, pct_flagged := 100 * flagged_cells / pmax(total_cells, 1)]

setorder(flag_summary, -pct_flagged, -flagged_cells)

# Suggested threshold based on upper tail
auto_pct_thresh <- as.numeric(quantile(flag_summary$pct_flagged, probs = 0.90, na.rm = TRUE))
auto_n_thresh   <- as.integer(round(as.numeric(quantile(flag_summary$flagged_cells, probs = 0.90, na.rm = TRUE))))

if (is.null(refine_threshold_pctflag)) refine_threshold_pctflag <- auto_pct_thresh
if (is.null(refine_threshold_nflag))   refine_threshold_nflag   <- max(10L, auto_n_thresh)

selected_fovs_dt <- flag_summary[
  pct_flagged >= refine_threshold_pctflag & flagged_cells >= refine_threshold_nflag
]

setorder(selected_fovs_dt, -pct_flagged, -flagged_cells)

# Save CSV
fwrite(flag_summary, file.path(summary_dir, "flag_summary_all_fovs.csv"))
fwrite(selected_fovs_dt, file.path(summary_dir, "selected_fovs_for_refinement.csv"))

# Save summary plots (do not display)
p_top_flagged <- ggplot(head(flag_summary, 30), aes(x = reorder(as.factor(fov), pct_flagged), y = pct_flagged)) +
  geom_col() +
  coord_flip() +
  theme_minimal(base_size = 12) +
  labs(
    title = "Top 30 FOVs by percent flagged cells",
    x = "FOV",
    y = "% flagged cells"
  )

ggsave(
  filename = file.path(plot_dir, "top30_fovs_by_pct_flagged.png"),
  plot = p_top_flagged,
  width = 8,
  height = 8,
  dpi = 300
)

p_hist_flagged <- ggplot(flag_summary, aes(x = pct_flagged)) +
  geom_histogram(bins = 30) +
  theme_minimal(base_size = 12) +
  labs(
    title = "Distribution of flagged-cell percentage across FOVs",
    x = "% flagged cells",
    y = "Number of FOVs"
  )

ggsave(
  filename = file.path(plot_dir, "histogram_pct_flagged_across_fovs.png"),
  plot = p_hist_flagged,
  width = 8,
  height = 6,
  dpi = 300
)

# Write screening summary TXT
screening_txt <- c(
  paste0("FastReseg multi-FOV flagging summary"),
  paste0("===================================="),
  "",
  paste0("Number of FOVs screened: ", nrow(flag_summary)),
  paste0("Total cells across screened FOVs: ", sum(flag_summary$total_cells)),
  paste0("Total flagged cells across screened FOVs: ", sum(flag_summary$flagged_cells)),
  paste0("Median % flagged across FOVs: ", round(median(flag_summary$pct_flagged), 3)),
  paste0("90th percentile % flagged across FOVs: ", round(auto_pct_thresh, 3)),
  paste0("Suggested refinement threshold (% flagged): ", round(refine_threshold_pctflag, 3)),
  paste0("Suggested refinement threshold (flagged cells): ", refine_threshold_nflag),
  "",
  paste0("Selected FOVs for full refinement under current threshold: ", nrow(selected_fovs_dt)),
  paste0("Selected FOV list: ", paste(selected_fovs_dt$fov, collapse = ", ")),
  "",
  "Interpretation suggestion:",
  "- Use this screening output to focus full refinement on the most affected FOVs.",
  "- If very few FOVs are selected, the dataset may already be largely stable under these conservative settings.",
  "- If you want a more sensitive second pass later, lower flagCell_lrtest_cutoff from 5 to 2."
)

screening_txt_file <- file.path(summary_dir, "FastReseg_screening_summary.txt")
writeLines(screening_txt, screening_txt_file)

cat("\nDataset-wide flagging completed.\n")
cat("Screening summary saved to:\n", screening_txt_file, "\n")
cat("Suggested threshold now in use:\n")
cat("- pct_flagged >=", round(refine_threshold_pctflag, 3), "\n")
cat("- flagged_cells >=", refine_threshold_nflag, "\n")

if (!run_full_refine) {
  cat("\nrun_full_refine is FALSE.\n")
  cat("Review the CSV/TXT summary first, then set run_full_refine <- TRUE to refine selected FOVs.\n")
}

# ----------------------------
# 9) Helper functions for per-FOV full refinement and saved plots
# ----------------------------
make_local_color_map <- function(before_ids, after_ids) {
  before_ids <- sort(unique(as.character(before_ids)))
  after_ids  <- sort(unique(as.character(after_ids)))
  
  shared_ids  <- intersect(before_ids, after_ids)
  before_only <- setdiff(before_ids, after_ids)
  after_only  <- setdiff(after_ids, before_ids)
  
  ordered_ids <- c(shared_ids, before_only, after_only)
  cols <- hue_pal()(max(length(ordered_ids), 1))
  names(cols) <- ordered_ids
  cols
}

get_plot_limits <- function(df, pad = 5) {
  list(
    xlim = range(df$x_um, na.rm = TRUE) + c(-pad, pad),
    ylim = range(df$y_um, na.rm = TRUE) + c(-pad, pad)
  )
}

save_before_after_plots_for_fov <- function(updated_transDF, meta_fov, poly_fov, out_dir_one_fov,
                                            top_n = 3, avg_n = 2, radius_px = 120,
                                            pixel_size_um = 0.12028) {
  dir.create(out_dir_one_fov, recursive = TRUE, showWarnings = FALSE)
  
  upd <- copy(as.data.table(updated_transDF))
  upd[, changed := UMI_cellID != updated_cellID]
  
  # Top-changed cells
  cell_change_counts <- upd[changed == TRUE, .N, by = UMI_cellID]
  setorder(cell_change_counts, -N)
  
  top_cells <- head(cell_change_counts$UMI_cellID, top_n)
  
  # Average-like changed cells: choose cells whose changed count is closest to median among nonzero changed cells
  avg_cells <- character(0)
  if (nrow(cell_change_counts) > 0) {
    med_n <- median(cell_change_counts$N)
    avg_cells <- cell_change_counts[order(abs(N - med_n))][1:min(avg_n, .N), UMI_cellID]
  }
  
  chosen_cells <- unique(c(top_cells, avg_cells))
  if (length(chosen_cells) == 0) return(invisible(NULL))
  
  # Build plotting tables
  poly_local_all <- copy(poly_fov)
  poly_local_all[, x_um := x_local_px * pixel_size_um]
  poly_local_all[, y_um := y_local_px * pixel_size_um]
  
  tx_before <- copy(upd)
  tx_before[, plot_cellID := UMI_cellID]
  tx_before[, x_um := x]
  tx_before[, y_um := y]
  
  tx_after <- copy(upd)
  tx_after[, plot_cellID := updated_cellID]
  tx_after[, x_um := x]
  tx_after[, y_um := y]
  
  for (tc in chosen_cells) {
    target_row <- meta_fov[cell == tc]
    if (nrow(target_row) == 0) next
    
    cx <- target_row$CenterX_local_px[1]
    cy <- target_row$CenterY_local_px[1]
    
    meta_tmp <- copy(meta_fov)
    meta_tmp[, dist_px := sqrt((CenterX_local_px - cx)^2 + (CenterY_local_px - cy)^2)]
    near_cells <- meta_tmp[dist_px <= radius_px, cell]
    
    poly_local <- poly_local_all[cell %in% near_cells]
    tx_before_local <- tx_before[UMI_cellID %in% near_cells | updated_cellID %in% near_cells]
    tx_after_local  <- tx_after[UMI_cellID %in% near_cells | updated_cellID %in% near_cells]
    
    local_colors <- make_local_color_map(
      before_ids = tx_before_local$plot_cellID,
      after_ids  = tx_after_local$plot_cellID
    )
    
    # BEFORE
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
      theme(plot.title = element_text(size = 16), legend.position = "right")
    
    ggsave(
      filename = file.path(out_dir_one_fov, paste0("before_", tc, ".png")),
      plot = p_before,
      width = 10,
      height = 8,
      dpi = 300
    )
    
    # AFTER
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
      theme(plot.title = element_text(size = 16), legend.position = "right")
    
    ggsave(
      filename = file.path(out_dir_one_fov, paste0("after_", tc, ".png")),
      plot = p_after,
      width = 10,
      height = 8,
      dpi = 300
    )
  }
}

# ----------------------------
# 10) Optional Phase 2: full refinement only for selected FOVs
# ----------------------------
refine_manifest <- list()

if (run_full_refine) {
  cat("\n[Phase 2] Running full refinement on selected FOVs ...\n")
  
  if (nrow(selected_fovs_dt) == 0) {
    cat("No FOVs passed the current refinement threshold. Skipping full refinement.\n")
  } else {
    fovs_to_refine <- head(selected_fovs_dt$fov, max_refine_fovs)
    cat("FOVs to refine:", paste(fovs_to_refine, collapse = ", "), "\n")
    
    for (this_fov in fovs_to_refine) {
      cat("\n--- Refining FOV", this_fov, "---\n")
      
      one_fov_dir <- file.path(refine_dir, paste0("FOV_", this_fov))
      dir.create(one_fov_dir, recursive = TRUE, showWarnings = FALSE)
      
      # one-row file info so wrapper still prefixes IDs consistently
      file_info_one <- file_info[fov == this_fov]
      
      # run full pipeline for this one FOV
      set.seed(1)
      res_full_one <- fastReseg_full_pipeline(
        counts = counts,
        clust = clust,
        refProfiles = NULL,
        
        transDF_fileInfo = as.data.frame(file_info_one),
        filepath_coln = "file_path",
        prefix_colns = c("slide", "fov"),
        fovOffset_colns = c("stage_X", "stage_Y"),
        
        pixel_size = pixel_size_um,
        zstep_size = zstep_size_um,
        
        transcript_df = NULL,
        transID_coln = NULL,
        transGene_coln = "target",
        cellID_coln = "cell_ID",
        spatLocs_colns = c("x_local_px", "y_local_px", "z"),
        invert_y = TRUE,
        
        extracellular_cellID = c(0),
        
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
        
        path_to_output = one_fov_dir,
        transDF_export_option = 2,
        save_intermediates = TRUE,
        return_perCellData = TRUE,
        combine_extra = FALSE,
        
        seed_process = 1,
        percentCores = percentCores_full
      )
      
      # Save full result object for this FOV
      fov_rds <- file.path(one_fov_dir, paste0("FastReseg_full_FOV_", this_fov, ".rds"))
      saveRDS(res_full_one, fov_rds)
      
      # Pull updated transcript table
      updated_transDF <- NULL
      if ("updated_transDF_list" %in% names(res_full_one) && length(res_full_one$updated_transDF_list) >= 1) {
        updated_transDF <- as.data.table(res_full_one$updated_transDF_list[[1]])
      }
      
      # Read raw metadata / polygons for this FOV for plotting and summaries
      meta_fov <- fread(
        cmd = sprintf(
          "gzip -dc %s | awk -F',' 'NR==1 || $3==%d'",
          shQuote(meta_file), this_fov
        )
      )
      if (!"cell" %in% colnames(meta_fov)) {
        meta_fov[, cell := paste0("c_1_", fov, "_", cell_ID)]
      }
      
      poly_file <- "Raw Data/Colon_TMA-polygons.csv.gz"
      poly_fov <- fread(
        cmd = sprintf(
          "gzip -dc %s | awk -F',' 'NR==1 || $1==%d'",
          shQuote(poly_file), this_fov
        )
      )
      
      # Per-FOV summary text
      if (!is.null(updated_transDF)) {
        upd <- copy(updated_transDF)
        upd[, changed := UMI_cellID != updated_cellID]
        
        total_transcripts <- nrow(upd)
        changed_transcripts <- sum(upd$changed)
        pct_changed <- 100 * changed_transcripts / pmax(total_transcripts, 1)
        
        cell_changes <- upd[, .(total = .N, changed = sum(changed)), by = UMI_cellID]
        cell_changes[, pct_changed := 100 * changed / pmax(total, 1)]
        
        mean_pct <- mean(cell_changes$pct_changed)
        median_pct <- median(cell_changes$pct_changed)
        max_pct <- max(cell_changes$pct_changed)
        cells_gt5 <- sum(cell_changes$pct_changed > 5)
        cells_gt10 <- sum(cell_changes$pct_changed > 10)
        
        fov_txt <- c(
          paste0("FastReseg Results Summary (FOV ", this_fov, ")"),
          "========================================",
          "",
          "Dataset size:",
          paste0("- Cells: ", nrow(meta_fov)),
          paste0("- Transcripts: ", total_transcripts),
          "",
          "Resegmentation:",
          paste0("- Transcripts changed: ", changed_transcripts),
          paste0("- Percent changed: ", round(pct_changed, 4), "%"),
          "",
          "Per-cell changes:",
          paste0("- Mean % change: ", round(mean_pct, 4)),
          paste0("- Median % change: ", round(median_pct, 4)),
          paste0("- Max % change: ", round(max_pct, 4)),
          "",
          "Cells affected:",
          paste0("- >5% changed: ", cells_gt5),
          paste0("- >10% changed: ", cells_gt10),
          "",
          "Interpretation:",
          if (pct_changed < 0.5 && max_pct < 10) {
            "Segmentation is mostly correct. Only minor local corrections were needed."
          } else if (pct_changed < 2) {
            "Moderate corrections detected. Some local improvements present."
          } else {
            "Significant segmentation issues detected. Resegmentation likely useful."
          }
        )
        
        fov_txt_file <- file.path(one_fov_dir, paste0("FastReseg_summary_FOV_", this_fov, ".txt"))
        writeLines(fov_txt, fov_txt_file)
        
        # Save before/after plots for top-changed and average-like changed cells
        plot_out_dir <- file.path(one_fov_dir, "key_cell_plots")
        save_before_after_plots_for_fov(
          updated_transDF = upd,
          meta_fov = meta_fov,
          poly_fov = poly_fov,
          out_dir_one_fov = plot_out_dir,
          top_n = top_changed_n,
          avg_n = avg_changed_n,
          radius_px = neighbor_radius_px,
          pixel_size_um = pixel_size_um
        )
        
        refine_manifest[[as.character(this_fov)]] <- list(
          fov = this_fov,
          rds_file = fov_rds,
          summary_txt = fov_txt_file,
          plot_dir = plot_out_dir,
          transcripts_changed = changed_transcripts,
          pct_changed = pct_changed
        )
      } else {
        refine_manifest[[as.character(this_fov)]] <- list(
          fov = this_fov,
          rds_file = fov_rds,
          summary_txt = NA_character_,
          plot_dir = NA_character_,
          transcripts_changed = NA_integer_,
          pct_changed = NA_real_
        )
      }
      
      # free some memory each loop
      rm(res_full_one, updated_transDF, meta_fov, poly_fov)
      gc()
    }
  }
}

# ----------------------------
# 11) Save master screen object for later revisit
# ----------------------------
master_obj <- list(
  config = list(
    adata_file = adata_file,
    tx_all_file = tx_all_file,
    meta_file = meta_file,
    fov_pos_file = fov_pos_file,
    cluster_col = cluster_col,
    pixel_size_um = pixel_size_um,
    zstep_size_um = zstep_size_um,
    flagModel_TransNum_cutoff = flagModel_TransNum_cutoff,
    flagCell_lrtest_cutoff = flagCell_lrtest_cutoff,
    svmClass_score_cutoff = svmClass_score_cutoff,
    refine_threshold_pctflag = refine_threshold_pctflag,
    refine_threshold_nflag = refine_threshold_nflag
  ),
  file_info = file_info,
  flag_summary = flag_summary,
  selected_fovs = selected_fovs_dt,
  screening_summary_txt = screening_txt_file,
  flagging_rds = file.path(run_root, "FastReseg_flagging_allFOVs.rds"),
  refinement_manifest = refine_manifest
)

# save main flagging result separately
saveRDS(res_flag_all, file.path(run_root, "FastReseg_flagging_allFOVs.rds"))
saveRDS(master_obj, master_rds_file)

cat("\nAll requested screen outputs are saved.\n")
cat("Master object:\n", master_rds_file, "\n")
cat("Flagging result object:\n", file.path(run_root, "FastReseg_flagging_allFOVs.rds"), "\n")
cat("Screening TXT summary:\n", screening_txt_file, "\n")