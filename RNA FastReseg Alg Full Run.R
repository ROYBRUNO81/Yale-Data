# ============================================================
# FASTRESEG ALL-FOV SCREENING + OPTIONAL TARGETED REFINEMENT
# PURE RAW MODE (NO H5AD) - ROBUST COLUMN-CLEANING VERSION
# ============================================================

# ----------------------------
# 1) Load packages
# ----------------------------
suppressPackageStartupMessages({
  library(data.table)
  library(Matrix)
  library(ggplot2)
  library(scales)
  library(FastReseg)
})

# ----------------------------
# 2) User settings
# ----------------------------
expr_file   <- "Raw Data/Colon_TMA_exprMat_file.csv.gz"
meta_file   <- "Raw Data/Colon_TMA_metadata_file.csv.gz"
tx_file     <- "Raw Data/Colon_TMA_tx_file.csv.gz"
poly_file   <- "Raw Data/Colon_TMA-polygons.csv.gz"
fovpos_file <- "Raw Data/Colon_TMA_fov_positions_file.csv.gz"

cluster_col <- "RNA_RNA.QC.Module_Cell.Typing.InSituType.1_1_clusters"

run_root <- "FastReseg_allFOV_raw_mode"
dir.create(run_root, recursive = TRUE, showWarnings = FALSE)

tx_split_dir  <- file.path(run_root, "tx_by_fov")
screen_dir    <- file.path(run_root, "screen_outputs")
refine_dir    <- file.path(run_root, "refine_outputs")
summary_dir   <- file.path(run_root, "summary_outputs")
plot_dir      <- file.path(run_root, "plot_outputs")

dir.create(tx_split_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(screen_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(refine_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(summary_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

pixel_size_um <- 0.12028
zstep_size_um <- 0.8

flagModel_TransNum_cutoff <- 50
flagCell_lrtest_cutoff    <- 5
svmClass_score_cutoff     <- -2

percentCores_screen <- 0.50
percentCores_full   <- 0.25

run_full_refine <- FALSE
max_refine_fovs <- 10

top_changed_n <- 3
avg_changed_n <- 2
neighbor_radius_px <- 120

master_rds <- file.path(run_root, "FastReseg_master_raw_mode.rds")

cat("Run root:\n", normalizePath(run_root), "\n\n")

# ----------------------------
# 3) Helper functions
# ----------------------------
clean_names_dt <- function(dt) {
  nm <- names(dt)
  nm <- trimws(nm)
  names(dt) <- nm
  dt
}

force_colname <- function(dt, wanted) {
  nm <- names(dt)
  if (wanted %in% nm) return(dt)
  
  hit <- nm[tolower(nm) == tolower(wanted)]
  if (length(hit) == 1) {
    setnames(dt, hit, wanted)
    return(dt)
  }
  
  stop(sprintf("Could not find required column '%s'. Available columns include: %s",
               wanted, paste(head(nm, 20), collapse = ", ")))
}

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

# ----------------------------
# 4) Read raw metadata and expression safely
# ----------------------------
cat("[1/10] Reading raw metadata and expression...\n")

meta_all <- as.data.table(fread(meta_file))
expr_all <- as.data.table(fread(expr_file))

meta_all <- clean_names_dt(meta_all)
expr_all <- clean_names_dt(expr_all)

# Force critical column names
meta_all <- force_colname(meta_all, "fov")
meta_all <- force_colname(meta_all, "cell_ID")

expr_all <- force_colname(expr_all, "fov")
expr_all <- force_colname(expr_all, "cell_ID")

cat("Metadata columns (first 15):\n")
print(names(meta_all)[1:min(15, length(names(meta_all)))])
cat("\nExpression columns (first 15):\n")
print(names(expr_all)[1:min(15, length(names(expr_all)))])
cat("\n")

# Build cell labels
if (!"cell" %in% names(meta_all)) {
  meta_all[, cell := paste0("c_1_", fov, "_", cell_ID)]
}
if (!"cell" %in% names(expr_all)) {
  expr_all[, cell := paste0("c_1_", fov, "_", cell_ID)]
}

cat("Metadata rows:", nrow(meta_all), "\n")
cat("Expression rows:", nrow(expr_all), "\n\n")

# ----------------------------
# 5) Match metadata and expression
# ----------------------------
cat("[2/10] Matching metadata and expression rows...\n")

common_cells <- intersect(meta_all[["cell"]], expr_all[["cell"]])

meta_all <- meta_all[meta_all[["cell"]] %in% common_cells, ]
expr_all <- expr_all[expr_all[["cell"]] %in% common_cells, ]

setkey(meta_all, cell)
setkey(expr_all, cell)
expr_all <- expr_all[meta_all$cell]

cat("Matched cells:", length(common_cells), "\n\n")

# ----------------------------
# 6) Build counts matrix and cluster vector
# ----------------------------
cat("[3/10] Building counts matrix and cluster vector...\n")

if (!(cluster_col %in% names(meta_all))) {
  stop(sprintf("Cluster column '%s' not found in metadata.", cluster_col))
}

id_cols <- intersect(c("fov", "cell_ID", "cell"), names(expr_all))
gene_cols <- setdiff(names(expr_all), id_cols)

counts_all <- as.matrix(expr_all[, ..gene_cols])
rownames(counts_all) <- expr_all$cell
storage.mode(counts_all) <- "numeric"
counts_all <- Matrix(counts_all, sparse = TRUE)

clust_all <- as.character(meta_all[[cluster_col]])
names(clust_all) <- meta_all$cell

cell_counts_per_fov <- as.data.table(table(meta_all[["fov"]]))
setnames(cell_counts_per_fov, c("fov", "total_cells"))
cell_counts_per_fov[, fov := as.integer(as.character(fov))]
cell_counts_per_fov[, total_cells := as.integer(total_cells)]
setorder(cell_counts_per_fov, fov)

cat("Counts matrix dim:", dim(counts_all)[1], "x", dim(counts_all)[2], "\n")
cat("Unique clusters:", length(unique(clust_all)), "\n")
cat("Unique FOVs:", length(unique(meta_all[["fov"]])), "\n\n")

# ----------------------------
# 7) Split transcript file by FOV
# ----------------------------
cat("[4/10] Preparing per-FOV transcript files...\n")

existing_split <- list.files(tx_split_dir, pattern = "^FOV_[0-9]{3}_tx\\.csv$", full.names = TRUE)

if (length(existing_split) < 300) {
  cat("Splitting raw transcript file by FOV...\n")
  
  split_cmd <- sprintf(
    "mkdir -p %s && gzip -dc %s | awk -F',' 'NR==1{header=$0; next} {outfile=sprintf(\"%s/FOV_%%03d_tx.csv\", $1); if (!(outfile in seen)) {print header > outfile; seen[outfile]=1} print >> outfile}'",
    shQuote(tx_split_dir),
    shQuote(tx_file),
    tx_split_dir
  )
  
  status <- system(split_cmd)
  if (status != 0) stop("Transcript split failed.")
} else {
  cat("Transcript split files already present. Reusing.\n")
}

tx_files <- list.files(tx_split_dir, pattern = "^FOV_[0-9]{3}_tx\\.csv$", full.names = TRUE)
cat("Transcript FOV files:", length(tx_files), "\n\n")

# ----------------------------
# 8) Build transDF_fileInfo
# ----------------------------
cat("[5/10] Building transDF_fileInfo...\n")

fov_pos <- as.data.table(fread(fovpos_file))
fov_pos <- clean_names_dt(fov_pos)
fov_pos <- force_colname(fov_pos, "fov")

# Handle FOV file if it came as uppercase
if ("FOV" %in% names(fov_pos) && !("fov" %in% names(fov_pos))) {
  setnames(fov_pos, "FOV", "fov")
}

required_fovpos <- c("fov", "x_global_mm", "y_global_mm")
missing_fovpos <- setdiff(required_fovpos, names(fov_pos))
if (length(missing_fovpos) > 0) {
  stop(sprintf("Missing required columns in FOV positions file: %s",
               paste(missing_fovpos, collapse = ", ")))
}

fov_pos[, stage_X := x_global_mm * 1000]
fov_pos[, stage_Y := y_global_mm * 1000]
fov_pos[, slide := 1L]

file_info <- data.table(
  file_path = tx_files,
  fov = as.integer(sub("^FOV_([0-9]{3})_tx\\.csv$", "\\1", basename(tx_files)))
)

file_info <- merge(file_info, fov_pos[, .(fov, stage_X, stage_Y, slide)], by = "fov", all.x = TRUE)
file_info <- file_info[file_info[["fov"]] %in% unique(meta_all[["fov"]]), ]
setorder(file_info, fov)

cat("transDF_fileInfo rows:", nrow(file_info), "\n\n")

# ----------------------------
# 9) Run dataset-wide flagging
# ----------------------------
cat("[6/10] Running dataset-wide FastReseg flagging...\n")
cat("This is the screening pass across all FOVs.\n")
cat("Suggested first-pass settings:\n")
cat("- flagCell_lrtest_cutoff =", flagCell_lrtest_cutoff, "\n")
cat("- svmClass_score_cutoff  =", svmClass_score_cutoff, "\n")
cat("- flagModel_TransNum_cutoff =", flagModel_TransNum_cutoff, "\n\n")

set.seed(1)

res_flag_all <- fastReseg_flag_all_errors(
  counts = counts_all,
  clust = clust_all,
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
  
  path_to_output = screen_dir,
  transDF_export_option = 2,
  return_trimmed_perCell = FALSE,
  combine_extra = FALSE,
  
  seed_transError = 1,
  percentCores = percentCores_screen
)

saveRDS(res_flag_all, file.path(run_root, "FastReseg_flagging_allFOVs_raw_mode.rds"))

# ----------------------------
# 10) Summarize flagging across FOVs
# ----------------------------
cat("[7/10] Building all-FOV flag summary...\n")

flagged_list <- res_flag_all$combined_flaggedCells
flagged_counts <- data.table(
  fov = file_info$fov,
  flagged_cells = vapply(flagged_list, length, integer(1))
)

flag_summary <- merge(cell_counts_per_fov, flagged_counts, by = "fov", all.x = TRUE)
flag_summary[is.na(flagged_cells), flagged_cells := 0L]
flag_summary[, pct_flagged := 100 * flagged_cells / pmax(total_cells, 1)]
setorder(flag_summary, -pct_flagged, -flagged_cells)

fwrite(flag_summary, file.path(summary_dir, "flag_summary_all_fovs.csv"))

pct_thresh <- as.numeric(quantile(flag_summary$pct_flagged, probs = 0.90, na.rm = TRUE))
n_thresh   <- as.integer(round(as.numeric(quantile(flag_summary$flagged_cells, probs = 0.90, na.rm = TRUE))))
n_thresh   <- max(10L, n_thresh)

selected_fovs <- flag_summary[pct_flagged >= pct_thresh & flagged_cells >= n_thresh]
fwrite(selected_fovs, file.path(summary_dir, "selected_fovs_for_refinement.csv"))

summary_txt <- c(
  "FastReseg all-FOV screening summary (pure raw mode)",
  "===================================================",
  "",
  paste0("Number of FOVs screened: ", nrow(flag_summary)),
  paste0("Total cells screened: ", sum(flag_summary$total_cells)),
  paste0("Total flagged cells: ", sum(flag_summary$flagged_cells)),
  "",
  paste0("Median % flagged across FOVs: ", round(median(flag_summary$pct_flagged), 3)),
  paste0("90th percentile % flagged: ", round(pct_thresh, 3)),
  paste0("Suggested threshold for full refinement: pct_flagged >= ", round(pct_thresh, 3),
         " and flagged_cells >= ", n_thresh),
  "",
  paste0("Selected FOV count: ", nrow(selected_fovs)),
  paste0("Selected FOVs: ", paste(selected_fovs$fov, collapse = ", ")),
  "",
  "Interpretation guidance:",
  "- Start full refinement only on the selected FOVs.",
  "- If too many FOVs are selected, raise the threshold slightly.",
  "- If too few are selected, lower the threshold or use lrtest cutoff = 2 in a second pass."
)

summary_txt_file <- file.path(summary_dir, "FastReseg_screening_summary.txt")
writeLines(summary_txt, summary_txt_file)

p_top <- ggplot(head(flag_summary, 30), aes(x = reorder(as.factor(fov), pct_flagged), y = pct_flagged)) +
  geom_col() +
  coord_flip() +
  theme_minimal(base_size = 12) +
  labs(title = "Top 30 FOVs by % flagged cells", x = "FOV", y = "% flagged")

ggsave(
  file.path(plot_dir, "top30_fovs_by_pct_flagged.png"),
  p_top, width = 8, height = 8, dpi = 300
)

p_hist <- ggplot(flag_summary, aes(x = pct_flagged)) +
  geom_histogram(bins = 30) +
  theme_minimal(base_size = 12) +
  labs(title = "Distribution of % flagged across FOVs", x = "% flagged", y = "FOV count")

ggsave(
  file.path(plot_dir, "histogram_pct_flagged_across_fovs.png"),
  p_hist, width = 8, height = 6, dpi = 300
)

cat("Screening summary saved to:\n", summary_txt_file, "\n")
cat("Suggested full-refinement threshold:\n")
cat("- pct_flagged >=", round(pct_thresh, 3), "\n")
cat("- flagged_cells >=", n_thresh, "\n\n")

# ----------------------------
# 11) Optional full refinement for selected FOVs
# ----------------------------
refine_manifest <- list()

if (run_full_refine) {
  cat("[8/10] Running targeted full refinement...\n")
  
  if (nrow(selected_fovs) == 0) {
    cat("No FOVs passed threshold. Skipping full refinement.\n")
  } else {
    fovs_to_refine <- head(selected_fovs$fov, max_refine_fovs)
    cat("Refining FOVs:\n")
    print(fovs_to_refine)
    cat("\n")
    
    for (this_fov in fovs_to_refine) {
      cat("Refining FOV", this_fov, "...\n")
      
      one_fov_dir <- file.path(refine_dir, paste0("FOV_", this_fov))
      dir.create(one_fov_dir, recursive = TRUE, showWarnings = FALSE)
      
      meta_fov <- meta_all[meta_all[["fov"]] == this_fov, ]
      expr_fov <- expr_all[expr_all[["fov"]] == this_fov, ]
      
      common_cells_fov <- intersect(meta_fov[["cell"]], expr_fov[["cell"]])
      meta_fov <- meta_fov[meta_fov[["cell"]] %in% common_cells_fov, ]
      expr_fov <- expr_fov[expr_fov[["cell"]] %in% common_cells_fov, ]
      
      setkey(meta_fov, cell)
      setkey(expr_fov, cell)
      expr_fov <- expr_fov[meta_fov$cell]
      
      id_cols <- intersect(c("fov", "cell_ID", "cell"), names(expr_fov))
      gene_cols <- setdiff(names(expr_fov), id_cols)
      
      counts_fov <- as.matrix(expr_fov[, ..gene_cols])
      rownames(counts_fov) <- expr_fov$cell
      storage.mode(counts_fov) <- "numeric"
      counts_fov <- Matrix(counts_fov, sparse = TRUE)
      
      clust_fov <- as.character(meta_fov[[cluster_col]])
      names(clust_fov) <- meta_fov$cell
      
      file_info_one <- file_info[file_info[["fov"]] == this_fov, ]
      
      set.seed(1)
      res_full <- fastReseg_full_pipeline(
        counts = counts_fov,
        clust = clust_fov,
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
      
      fov_rds <- file.path(one_fov_dir, paste0("FastReseg_full_FOV_", this_fov, ".rds"))
      saveRDS(res_full, fov_rds)
      
      updated_transDF <- NULL
      if ("updated_transDF_list" %in% names(res_full) && length(res_full$updated_transDF_list) >= 1) {
        updated_transDF <- as.data.table(res_full$updated_transDF_list[[1]])
      }
      
      if (!is.null(updated_transDF)) {
        updated_transDF[, changed := UMI_cellID != updated_cellID]
        
        total_transcripts <- nrow(updated_transDF)
        changed_transcripts <- sum(updated_transDF$changed)
        pct_changed <- 100 * changed_transcripts / pmax(total_transcripts, 1)
        
        cell_changes <- updated_transDF[, .(total = .N, changed = sum(changed)), by = UMI_cellID]
        cell_changes[, pct_changed := 100 * changed / pmax(total, 1)]
        
        mean_pct <- mean(cell_changes$pct_changed)
        median_pct <- median(cell_changes$pct_changed)
        max_pct <- max(cell_changes$pct_changed)
        cells_gt5 <- sum(cell_changes$pct_changed > 5)
        cells_gt10 <- sum(cell_changes$pct_changed > 10)
        
        summary_lines <- c(
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
        
        summary_file <- file.path(one_fov_dir, paste0("FastReseg_summary_FOV_", this_fov, ".txt"))
        writeLines(summary_lines, summary_file)
        
        poly_fov <- as.data.table(fread(
          cmd = sprintf(
            "gzip -dc %s | awk -F',' 'NR==1 || $1==%d'",
            shQuote(poly_file), this_fov
          )
        ))
        
        poly_fov <- clean_names_dt(poly_fov)
        if (!"cell" %in% names(poly_fov)) stop("Column 'cell' not found in polygon file.")
        
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
        
        cell_change_counts <- updated_transDF[changed == TRUE, .N, by = UMI_cellID]
        setorder(cell_change_counts, -N)
        top_cells <- head(cell_change_counts$UMI_cellID, top_changed_n)
        
        avg_cells <- character(0)
        if (nrow(cell_change_counts) > 0) {
          med_n <- median(cell_change_counts$N)
          avg_cells <- cell_change_counts[order(abs(N - med_n))][1:min(avg_changed_n, .N), UMI_cellID]
        }
        
        chosen_cells <- unique(c(top_cells, avg_cells))
        plot_save_dir <- file.path(one_fov_dir, "key_cell_plots")
        dir.create(plot_save_dir, recursive = TRUE, showWarnings = FALSE)
        
        for (tc in chosen_cells) {
          target_row <- meta_fov[meta_fov[["cell"]] == tc, ]
          if (nrow(target_row) == 0) next
          
          cx <- target_row$CenterX_local_px[1]
          cy <- target_row$CenterY_local_px[1]
          
          meta_fov[, dist := sqrt((CenterX_local_px - cx)^2 + (CenterY_local_px - cy)^2)]
          near_cells <- meta_fov[dist <= neighbor_radius_px, cell]
          
          poly_local <- poly_fov[poly_fov[["cell"]] %in% near_cells, ]
          tx_before_local <- tx_before[tx_before[["UMI_cellID"]] %in% near_cells | tx_before[["updated_cellID"]] %in% near_cells, ]
          tx_after_local  <- tx_after[tx_after[["UMI_cellID"]] %in% near_cells | tx_after[["updated_cellID"]] %in% near_cells, ]
          
          local_colors <- make_local_color_map(tx_before_local$plot_cellID, tx_after_local$plot_cellID)
          
          lims_before <- get_plot_limits(tx_before_local)
          p_before <- ggplot() +
            geom_path(data = poly_local, aes(x = x_um, y = y_um, group = cell),
                      color = "grey55", linewidth = 0.35, alpha = 0.8) +
            geom_point(data = tx_before_local,
                       aes(x = x_um, y = y_um, color = plot_cellID),
                       size = 1.4, alpha = 0.9) +
            scale_color_manual(values = local_colors, drop = FALSE) +
            coord_fixed(xlim = lims_before$xlim, ylim = lims_before$ylim) +
            labs(title = paste("BEFORE -", tc), x = "x (micron)", y = "y (micron)", color = "original cell ID") +
            theme_minimal(base_size = 14) +
            theme(plot.title = element_text(size = 16), legend.position = "right")
          
          ggsave(
            filename = file.path(plot_save_dir, paste0("before_", tc, ".png")),
            plot = p_before, width = 10, height = 8, dpi = 300
          )
          
          lims_after <- get_plot_limits(tx_after_local)
          p_after <- ggplot() +
            geom_path(data = poly_local, aes(x = x_um, y = y_um, group = cell),
                      color = "grey55", linewidth = 0.35, alpha = 0.8) +
            geom_point(data = tx_after_local,
                       aes(x = x_um, y = y_um, color = plot_cellID),
                       size = 1.4, alpha = 0.9) +
            scale_color_manual(values = local_colors, drop = FALSE) +
            coord_fixed(xlim = lims_after$xlim, ylim = lims_after$ylim) +
            labs(title = paste("AFTER -", tc), x = "x (micron)", y = "y (micron)", color = "updated cell ID") +
            theme_minimal(base_size = 14) +
            theme(plot.title = element_text(size = 16), legend.position = "right")
          
          ggsave(
            filename = file.path(plot_save_dir, paste0("after_", tc, ".png")),
            plot = p_after, width = 10, height = 8, dpi = 300
          )
        }
        
        refine_manifest[[as.character(this_fov)]] <- list(
          fov = this_fov,
          rds_file = fov_rds,
          summary_txt = summary_file,
          plot_dir = plot_save_dir,
          transcripts_changed = changed_transcripts,
          pct_changed = pct_changed
        )
      }
      
      rm(res_full, updated_transDF, counts_fov, clust_fov, poly_fov)
      gc()
    }
  }
} else {
  cat("run_full_refine is FALSE. Screening only was completed.\n")
}

# ----------------------------
# 12) Save master object
# ----------------------------
master_obj <- list(
  config = list(
    expr_file = expr_file,
    meta_file = meta_file,
    tx_file = tx_file,
    fovpos_file = fovpos_file,
    cluster_col = cluster_col,
    pixel_size_um = pixel_size_um,
    zstep_size_um = zstep_size_um,
    flagModel_TransNum_cutoff = flagModel_TransNum_cutoff,
    flagCell_lrtest_cutoff = flagCell_lrtest_cutoff,
    svmClass_score_cutoff = svmClass_score_cutoff,
    run_full_refine = run_full_refine
  ),
  flag_summary = flag_summary,
  selected_fovs = selected_fovs,
  screening_summary_txt = summary_txt_file,
  refinement_manifest = refine_manifest
)

saveRDS(master_obj, master_rds)

cat("\n[10/10] Done.\n")
cat("Master object saved to:\n", master_rds, "\n")
cat("Screening summary saved to:\n", summary_txt_file, "\n")