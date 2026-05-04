############################################################
## Ribosomal protein abundance heatmap in CPTAC COAD
## Data source: LinkedOmics CPTAC Colon
## https://linkedomics.org/cptac-colon/
##
## Reference:
## Vasaikar et al., Cell, 2019
############################################################

## 1. Load packages
library(ComplexHeatmap)
library(circlize)
library(grid)

## 2. Set working directory
setwd("CJX_workspace/RPs/")

## 3. Load CPTAC COAD tumor-normal proteome log2FC data
exp <- read.table(
  "Human__CPTAC_COAD__PNNL__Proteome__TMT__03_01_2017__BCM__Gene__Tumor_Normal_log2FC.cct",
  header = TRUE,
  row.names = 1,
  sep = "\t",
  check.names = FALSE
)

exp[1:3, 1:3]
mat <- as.matrix(exp)

## 4. Load ribosomal protein gene list
rp_genes <- read.csv(
  "../NPM1-BIOID2/NPM1_ribosomal_proteins_subset.csv",
  stringsAsFactors = FALSE
)

rp_genes <- rp_genes$Gene
rp_genes <- unique(rp_genes)

## 5. Extract ribosomal proteins from CPTAC matrix
matched_genes <- intersect(rp_genes, rownames(mat))
missing_genes <- setdiff(rp_genes, rownames(mat))

message("Matched genes: ", length(matched_genes))
message("Missing genes: ", length(missing_genes))

if (length(missing_genes) > 0) {
  message("Missing genes:")
  print(missing_genes)
}

data <- mat[matched_genes, , drop = FALSE]

## 6. Remove invalid values
data[is.na(data) | is.nan(data) | is.infinite(data)] <- 0

## 7. Limit extreme values for heatmap visualization
data_plot <- data
data_plot[data_plot > 2] <- 2
data_plot[data_plot < -2] <- -2

## 8. Calculate mean log2 ratio for each protein
protein_mean <- rowMeans(data, na.rm = TRUE)

## 9. Sort proteins by mean log2 ratio from high to low
row_order <- order(protein_mean, decreasing = TRUE)

data_ordered <- data_plot[row_order, , drop = FALSE]
mean_ordered <- protein_mean[row_order]

## 10. Define heatmap color scale
col_fun <- colorRamp2(
  c(-2, -0.5, 0, 0.5, 2),
  c("#FFF5F0", "#FCAE91", "#FB6A4A", "#9E1F63", "#3B0F3F")
)

## 11. Add right-side barplot showing mean log2 ratio
right_anno <- rowAnnotation(
  `Mean\nLog2 ratio` = anno_barplot(
    mean_ordered,
    gp = gpar(fill = "#67000D", col = NA),
    border = FALSE,
    width = unit(3, "cm")
  ),
  annotation_name_gp = gpar(fontsize = 7)
)

## 12. Draw heatmap
ht <- Heatmap(
  data_ordered,
  name = "Log2 ratio",
  col = col_fun,
  
  cluster_rows = FALSE,
  cluster_columns = TRUE,
  
  row_names_side = "right",
  show_column_names = FALSE,
  
  row_names_gp = gpar(fontsize = 4),
  
  width = unit(5, "cm"),
  height = unit(10, "cm")
)

ht_final <- ht + right_anno

## 13. Save as PDF
pdf(
  "heatmap_sorted_all_mean.pdf",
  width = 8,
  height = 10,
  useDingbats = FALSE
)

draw(ht_final)

dev.off()
