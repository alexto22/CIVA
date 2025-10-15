# -------- Code for cell type distribution --------------
options(repr.plot.width = 8, repr.plot.height = 7)
all_subset <- subset(all, subset = orig.ident %in% c("Flex1_assem", "Flex2_VIO"))
cell_dist <- table(all_subset@meta.data$orig.ident, all_subset@meta.data$cell_type)
cell_dist_df <- as.data.frame.table(cell_dist)
colnames(cell_dist_df) <- c("Sample", "CellType", "Count")

# Calculate proportions
cell_dist_df <- cell_dist_df %>% 
  group_by(Sample) %>% 
  mutate(Proportion = Count / sum(Count))

# Plot with predefined palette
ggplot(cell_dist_df, aes(x = Sample, y = Proportion, fill = CellType)) +
  geom_bar(stat = "identity", width = 0.6) +
  scale_fill_manual(values = pal) +  # Use your predefined palette
  theme_minimal() +
  labs(title = "Cell Type Distribution: Assembloids vs VIO", 
       x = "Sample", y = "Proportion of Cells", fill = "Cell Type") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "right",
        legend.key.size = unit(1, "cm"))

# -------- Code for cell types distribution across niches  --------------
# Load necessary libraries
library(Seurat)
library(ggplot2)

# Define the custom color palette (artistic style with immune highlighted, tumors dimmed brown-based)
pal <- c(
  "Fibro" = "#A6CEE3",        # Light blue (Fibro group: blues)
  "Fibro_POSTN" = "#1F78B4",  # Dark blue (Fibro group)
  "Fibro_B7-H4" = "#6495ED",  # Medium blue (Fibro group)
  "Tumor" = "#8B4513",        # Muted dark brown (Tumor group: browns, dimmed)
  "Tumor_PDL1" = "#D2B48C",   # Light tan brown (Tumor group, faded)
  "Macrophages" = "#FF4500",  # Vivid orange-red (highlighted immune)
  "Endothelial" = "#228B22",  # Warm peach (vascular)
  "ECM" = "#808080",          # Medium gray (neutral matrix)
  "NA" = "#D3D3D3"            # Light gray (subtle unknown)
)
options(repr.plot.width = 5, repr.plot.height = 5)

# Visualize subclusters
# Filter indices for orig.ident == "Assembloid"
idx <- combined$orig.ident == "Assembloid"

# Create a contingency table of niches vs. third_annotations_filtered
count_table <- table(combined$niches[idx], combined$third_annotations_filtered[idx])

# Convert table to data.frame for ggplot
count_df <- as.data.frame(count_table)
colnames(count_df) <- c("Niche", "CellType", "Count")

# Determine the unique cell types and their order (matching Seurat's behavior for color assignment)
if (is.factor(combined$third_annotations_filtered)) {
  unique_cell_types <- levels(combined$third_annotations_filtered)
} else {
  unique_cell_types <- sort(unique(combined$third_annotations_filtered))
}

# Assign colors from pal (ensures match to unique_cell_types)
colors <- pal[unique_cell_types]  # Subset pal to match your data's unique types if needed

# Create a stacked bar plot showing cell type proportions per niche
ggplot(count_df, aes(x = Niche, y = Count, fill = CellType)) +
  geom_bar(stat = "identity", position = "fill") + # Use "fill" for proportions; change to "stack" for counts
  scale_fill_manual(values = colors) + # Use custom colors matching previous plots
  labs(title = "Cell Type Distributions per Niche",
       x = "Niche",
       y = "Proportion",
       fill = "Cell Type") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) # Rotate x-axis labels if needed for readability


# ----------- code for calculating cell type correlations ---------
library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(RColorBrewer)
library(tibble)
library(clusterProfiler)
library(pheatmap)

VIO <- readRDS("/storage/riosugimuralab/alexto/A2_AssemFlex/VIO_reanalysed.rds")
Zemin_lung<-load("/storage/riosugimuralab/alexto/06_Public/250527_ZeminZhang2021/IM06_ALL_cells_Immune_AND_nonImmune_annotations.RData")

desired_annotations <- c("MF-Monocytes", "fibroblast", "endothelial", "Mast-cells", "neuro.e cell", "tumor")
Zemin.subset.celltypes <- subset(Zemin_lung, subset = Final_annotation %in% desired_annotations)


VIO.order <- rev(c("MonoMac", "Granulocyte","Endothelium","Fibroblast","Myofibroblast",
          "Neurons","Yolk","Trophoblast"))
Zemin.order <- rev(c("MF-Monocytes", "Mast-cells", "endothelial", "fibroblast", 
                    "neuro.e cell", "tumor"))


get_corr_mat<-function(obj, ref, 
                   obj.group.by, ref.group.by,
                      obj.order, ref.order){
    obj <- NormalizeData(obj) %>% FindVariableFeatures() %>% ScaleData()
    ref <- NormalizeData(ref) %>% FindVariableFeatures() %>% ScaleData()
    
    ref@meta.data[[ref.group.by]] <- factor(ref@meta.data[[ref.group.by]], levels = sort(unique(ref@meta.data[[ref.group.by]])))
    
    avg_ref <- as.matrix(AverageExpression(ref, group.by = ref.group.by, assays = "RNA", slot = "scale.data")$RNA)
    avg_obj <- as.matrix(AverageExpression(obj, group.by = obj.group.by, assays = "RNA", slot = "scale.data")$RNA)
        
    # Step 3: Find common genes (intersect rownames)
    common_genes <- intersect(rownames(avg_ref), rownames(avg_obj))
        # Check if there are common genes
    if (length(common_genes) < 2) {
      stop("Insufficient common genes for correlation analysis. Check your data.")
    }

        
    # Step 4: Compute correlation matrix (Spearman is robust for expression data; can change to "pearson")
    cor_mat <- cor(avg_obj[common_genes, ], avg_ref[common_genes, ], method = "spearman")

    'if (obj.order == FALSE){
    # Define the desired all order (using hyphens as per data)
    obj_order <- rev(unique(obj@meta.data[[obj.group.by]]))
    }'
    
    # Reorder the rows of cor_mat to reverse the desired order (to account for pheatmap's y-axis plotting from bottom to top)
    cor_mat <- cor_mat[rev(obj.order), ]
    cor_mat <- cor_mat[, rev(ref.order)]
    
    return(cor_mat)
        
}

options(repr.plot.height=7, repr.plot.width=4)

plot_corr_mat <- function(corr_mat, font=12,
                          colorset = rev(brewer.pal(n = 11, name = "BrBG"))) {
  p <- pheatmap(corr_mat,
                cluster_rows = FALSE,
                cluster_cols = FALSE,
                main = "Correlation Heatmap",
                fontsize = font,
                color = colorset)
  print(p)
}

cor_mat <- get_corr_mat(obj=VIO, ref=Zemin.subset.celltypes, 
                        obj.group.by = "cell_type", ref.group.by="Final_annotation",
                        obj.order = VIO.order, ref.order = Zemin.order
                        )

options(repr.plot.height=5, repr.plot.width=3.5)
library(scCustomize)
plot_corr_mat(corr_mat = cor_mat)
