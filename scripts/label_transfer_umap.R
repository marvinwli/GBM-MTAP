# label_transfer_and_umap_parallel.R

library(Seurat)
library(dplyr)
library(ggplot2)
library(future.apply)

# PARAMETERS
reference_path <- "/home/mli110/GBM-MTAP/rds objects/reference_merged_normalized_collapsed.RDS"
query_dir_root <- "/home/mli110/GBM-MTAP/data/GSE182109_RAW/"
output_dir <- "/home/mli110/GBM-MTAP/figures"
data_output_dir <- "/home/mli110/GBM-MTAP/data/"
umap_pdf_file <- file.path(output_dir, "query_umap_plots.pdf")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(data_output_dir, showWarnings = FALSE, recursive = TRUE)
options(future.globals.maxSize = 300 * 1024^3)
# --- Set parallel plan ---
# You can set this from a SLURM job, e.g. --cpus-per-task=8
plan("multisession", workers = 4) # Adjust number of workers to match your HPC allocation

# LOAD REFERENCE & PICK PCs
reference <- readRDS(reference_path)

dims_to_use <- 1:21
cat("Using dims 1:", "21", "for all queries.\n")

# LIST QUERY SAMPLES
sample_dirs <- list.dirs(path = query_dir_root, recursive = FALSE)

# PROCESSING FUNCTION
process_sample <- function(sample_dir) {
  tryCatch({
    sample_name <- basename(sample_dir)
    cat("Processing:", sample_name, "\n")
    
    # 1) LOAD raw 10X counts and make Seurat object
    query_data <- Read10X(sample_dir)
    query      <- CreateSeuratObject(query_data, project = sample_name)
    
    # 2) LOG‐NORMALIZE + VARIABLE FEATURES + SCALE (⟵ REPLACED the SCTransform block)
    DefaultAssay(query) <- "RNA"
    query <- NormalizeData(
      query,
      normalization.method = "LogNormalize",
      scale.factor        = 10000,
      verbose             = FALSE
    )
    query <- FindVariableFeatures(
      query,
      selection.method = "vst",
      nfeatures        = length(VariableFeatures(reference))  # or pick 2000 
    )
    query <- ScaleData(
      query,
      features        = VariableFeatures(query),
      vars.to.regress = c("nCount_RNA", "percent.mt"),  # if you’d like
      verbose         = FALSE
    )
    
    # 3) LABEL TRANSFER (use LogNormalize on both reference + query)
    anchors <- FindTransferAnchors(
      reference         = reference,
      query             = query,
      normalization.method = "LogNormalize",        # ⟵ CHANGED
      reference.assay   = "RNA",                    # ⟵ CHANGED
      query.assay       = "RNA",                    # ⟵ CHANGED
      dims              = dims_to_use
    )
    predictions <- TransferData(
      anchorset = anchors,
      refdata   = reference$celltype_collapsed,
      dims      = dims_to_use
    )
    query <- AddMetaData(query, predictions)
    
    # 4) PCA + UMAP on the query (using the “RNA” assay since that’s what we scaled)
    query <- RunPCA(query, npcs = max(dims_to_use), verbose = FALSE)
    query <- RunUMAP(query, dims = dims_to_use, verbose = FALSE)
    umap_plot <- DimPlot(
      query,
      reduction = "umap",
      group.by  = "predicted.id",
      label     = TRUE,
      repel     = TRUE
    ) + ggtitle(
      paste0(
        sub("^[^_]*_", "", sample_name),
        " – Transferred Labels (", length(dims_to_use), " PCs)"
      )
    )
    
    # 5) Cell type counts
    cell_type_counts <- as.data.frame(table(query$predicted.id), stringsAsFactors = FALSE)
    colnames(cell_type_counts) <- c("Cell_Type", "Count")
    cell_type_counts$Sample <- sample_name
    
    list(
      Sample           = sample_name,
      Cell_Type_Counts = cell_type_counts,
      umap_plot        = umap_plot,
      query            = query
    )
  }, error = function(e) {
    message("Error processing ", sample_dir, ": ", e$message)
    return(NULL)
  })
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

