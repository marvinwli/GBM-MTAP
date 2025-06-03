library(Seurat)
library(dplyr)
library(ggplot2)
options(future.globals.maxSize = 700 * 1024^3)

# Load and prepare reference
reference_list <- list(
  GBM030_tumor   = readRDS("/home/mli110/GBM-MTAP/data/reference/GBM030tumor.RDS"),
  GBM049_tumor   = readRDS("/home/mli110/GBM-MTAP/data/reference/GBM049tumor.RDS"),
  GBM030_lymph   = readRDS("/home/mli110/GBM-MTAP/data/reference/GBM030lymph.RDS"),
  GBM030_myeloid = readRDS("/home/mli110/GBM-MTAP/data/reference/GBM030myeloid.RDS"),
  GBM049_lymph   = readRDS("/home/mli110/GBM-MTAP/data/reference/GBM049lymph.RDS"),
  GBM049_myeloid = readRDS("/home/mli110/GBM-MTAP/data/reference/GBM049myeloid.RDS")
)

reference_list$GBM030_tumor$orig.ident   <- "GBM030tumor"
reference_list$GBM049_tumor$orig.ident   <- "GBM049tumor"
reference_list$GBM030_lymph$orig.ident   <- "GBM030lym"
reference_list$GBM030_myeloid$orig.ident <- "GBM030myl"
reference_list$GBM049_lymph$orig.ident   <- "GBM049lym"
reference_list$GBM049_myeloid$orig.ident <- "GBM049myl"

reference_list$GBM030_tumor$celltype   <- "tumor"
reference_list$GBM049_tumor$celltype   <- "tumor"
reference_list$GBM030_lymph$celltype   <- as.character(Idents(reference_list$GBM030_lymph))
reference_list$GBM030_myeloid$celltype <- as.character(Idents(reference_list$GBM030_myeloid))
reference_list$GBM049_lymph$celltype   <- as.character(Idents(reference_list$GBM049_lymph))
reference_list$GBM049_myeloid$celltype <- as.character(Idents(reference_list$GBM049_myeloid))

reference <- merge(
  reference_list[[1]],
  y = reference_list[-1],
  add.cell.ids = names(reference_list)
)

# Get all variable features from all references
all_var_features <- unique(unlist(lapply(reference_list, VariableFeatures)))
length(all_var_features)  # Just to check
VariableFeatures(reference) <- all_var_features

print(unique(reference$celltype))

label_map <- c(
  # Group CD8 T cell subtypes
  " CD8 EFF" = "CD8 T", "CD8 IFNG" = "CD8 T", "CD8 TEFF" = "CD8 T", 
  "CD8 RM-GZMK" = "CD8 T", "CD8 RM-XCL1" = "CD8 T", "CD8 Naive" = "CD8 T", 
  "MAIT T" = "CD8 T", "CD8 EX" = "CD8 T",
  # Group CD4 T cell subtypes
  "CD4 CM" = "CD4 T", "CD4 EM" = "CD4 T", "Th17" = "CD4 T", "CD4 Naive" = "CD4 T",
  # T-reg and Gamma-Delta separate
  "TREG" = "Treg", "GD T" = "GD T",        
  # Group microglia
  "MCG1" = "MCG", "MCG2" = "MCG", "MCG3" = "MCG", "MCG4" = "MCG", "MCG5" = "MCG",
  # Group MDSC
  "M-MDSC" = "MDSC","PMN-MDSC" = "MDSC", "E-MDSC" = "MDSC",
  # NK cells
  "NK1" = "NK", "NK2" = "NK", 
  # B cells
  "B cells" = "B", "Plasma" = "B",
  # Cycling Cells
  "CYC1" = "Cycling", "CYC2" = "Cycling",
  # Macrophages
  "MAC1" = "MAC1", "MAC2" = "MAC2",
  # Myeloid Dendritic Cells
  "MDC1" = "MDC", "MDC2" = "MDC",
  # Neutrophil
  "Neut"="Neut",
  # Tumor
  "tumor" = "Tumor"
)
reference$celltype_collapsed <- dplyr::recode(reference$celltype, !!!label_map)
table(reference$celltype_collapsed, useNA = "always")

DefaultAssay(reference) <- "RNA"

# 5a) Percentage of mitochondrial genes (for potential regression)
#     Adjust the pattern if using non-human gene names
reference[["percent.mt"]] <- PercentageFeatureSet(reference, pattern = "^MT-")

# 5b) Normalize counts per cell
reference <- NormalizeData(
  reference,
  normalization.method = "LogNormalize",
  scale.factor        = 10000,
  verbose             = TRUE
)

# 5c) Identify variable features (on log-normalized data)
reference <- FindVariableFeatures(
  reference,
  selection.method = "vst",
  nfeatures        = length(all_var_features)  # we already assigned them, but recalc if desired
)
# If you prefer to keep the union from step 3, skip the previous two lines and rely on VariableFeatures(reference)=all_var_features.

# 5d) Scale data and regress out technical covariates
reference <- ScaleData(
  reference,
  features       = VariableFeatures(reference),
  vars.to.regress = c("nCount_RNA", "percent.mt"),
  verbose        = TRUE
)

reference <- RunPCA(reference, npcs = 30)
elbow_plot <- ElbowPlot(reference, ndims = 30, reduction = "pca")
ggsave("/home/mli110/GBM-MTAP/figures/reference_elbow_plot.png", elbow_plot, width = 6, height = 5, dpi = 300)
reference <- RunUMAP(reference, dims = 1:11)
umap_collapsed <- DimPlot(reference, group.by = "celltype_collapsed", label = TRUE, repel = TRUE) +
  ggtitle("Reference UMAP by Collapsed Cell Type")
ggsave("/home/mli110/GBM-MTAP/figures/reference_umap_collapsed_celltype.png", umap_collapsed, width = 8, height = 6, dpi = 300)
umap_origident <- DimPlot(reference, group.by = "orig.ident", label = TRUE, repel = TRUE) + 
  ggtitle("Reference UMAP by Source Batch")
ggsave("/home/mli110/GBM-MTAP/figures/reference_umap_origident.png", umap_origident, width = 8, height = 6, dpi = 300)

saveRDS(reference, file = "/home/mli110/GBM-MTAP/rds objects/reference_merged_normalized_collapsed.RDS")