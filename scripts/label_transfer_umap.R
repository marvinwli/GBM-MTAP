# label_transfer_and_umap_parallel.R

library(Seurat)
library(dplyr)
library(ggplot2)
library(future.apply)

# PARAMETERS
reference_path <- "/home/mli110/GBM-MTAP/rds objects/reference_merged_normalized_collapsed.RDS"
query_dir_root <- "/home/mli110/GBM-MTAP/data/GSE182109_RAW"
output_dir <- "/home/mli110/GBM-MTAP/figures"
data_output_dir <- "/home/mli110/GBM-MTAP/data/"
umap_pdf_file <- file.path(output_dir, "query_umap_plots.pdf")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(data_output_dir, showWarnings = FALSE, recursive = TRUE)

# --- Set parallel plan ---
# You can set this from a SLURM job, e.g. --cpus-per-task=8
plan("multicore", workers = 44) # Adjust number of workers to match your HPC allocation

# LOAD REFERENCE & PICK PCs
reference <- readRDS(reference_path)
reference <- RunPCA(reference, npcs = 40)
stdevs <- reference[["pca"]]@stdev
pc_diffs <- diff(stdevs)
drop_threshold <- -0.1
min_pcs <- 15
max_pcs <- 30
elbow_pc <- which(pc_diffs > drop_threshold)[1]
if (is.na(elbow_pc)) {
  chosen_pcs <- max_pcs
} else if (elbow_pc < min_pcs) {
  chosen_pcs <- min_pcs
} else if (elbow_pc > max_pcs) {
  chosen_pcs <- max_pcs
} else {
  chosen_pcs <- elbow_pc
}
dims_to_use <- 1:chosen_pcs
cat("Using dims 1:", chosen_pcs, "for all queries.\n")

# LIST QUERY SAMPLES
sample_dirs <- list.dirs(path = query_dir_root, recursive = FALSE)

# PROCESSING FUNCTION
process_sample <- function(sample_dir) {
  sample_name <- basename(sample_dir)
  cat("Processing:", sample_name, "\n")
  
  query_data <- Read10X(sample_dir)
  query <- CreateSeuratObject(query_data, project = sample_name)
  query <- SCTransform(query, verbose = FALSE)
  VariableFeatures(query) <- intersect(VariableFeatures(reference), rownames(query))
  
  # Label transfer
  anchors <- FindTransferAnchors(reference = reference, query = query, normalization.method = "SCT", dims = dims_to_use)
  predictions <- TransferData(anchorset = anchors, refdata = reference$celltype_collapsed, dims = dims_to_use)
  query <- AddMetaData(query, predictions)
  
  # UMAP
  query <- RunPCA(query, npcs = max(dims_to_use))
  query <- RunUMAP(query, dims = dims_to_use)
  umap_plot <- DimPlot(query, group.by = "predicted.id", label = TRUE, repel = TRUE) +
    ggtitle(paste0(sub("^[^_]*_", "", sample_name), " - Transferred Labels (", length(dims_to_use), " PCs)"))
  
  # Cell type counts (optional)
  cell_type_counts <- as.data.frame(table(query$predicted.id))
  colnames(cell_type_counts) <- c("Cell_Type", "Count")
  cell_type_counts$Sample <- sample_name
  
  list(
    Sample = sample_name,
    Cell_Type_Counts = cell_type_counts,
    umap_plot = umap_plot,
    query = query
  )
}

# PROCESS ALL SAMPLES IN PARALLEL!
results_list <- future_lapply(sample_dirs, process_sample)
names(results_list) <- basename(sample_dirs)

# SAVE ALL UMAPs TO A SINGLE PDF
pdf(umap_pdf_file, width = 8, height = 6)
for (res in results_list) {
  print(res$umap_plot)
}
dev.off()

# SAVE THE MASTER RESULTS LIST (for future any analysis)
saveRDS(results_list, file = file.path(data_output_dir, "results_list.rds"))

# SAVE ALL CELL TYPE COUNTS (optional, helpful for future plotting)
all_counts <- bind_rows(lapply(results_list, function(x) x$Cell_Type_Counts))
write.csv(all_counts, file = file.path(data_output_dir, "results_celltype_counts_long.csv"), row.names = FALSE)

cat("All results saved.\n")

