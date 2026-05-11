############################################################
## Functional module heatmap of NPM1 interactors
## in CPTAC COAD proteome data
##
## Purpose:
## 1. Define functional modules from GO enrichment of
##    NPM1 high-confidence interactors.
## 2. Match these interactors to CPTAC COAD tumor/normal
##    proteome log2FC data.
## 3. Visualize module-level and gene-level protein abundance
##    changes in colorectal cancer.
############################################################

## =========================================================
## 0. Load packages
## =========================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(ComplexHeatmap)
  library(circlize)
  library(grid)
  library(ggplot2)
})

## =========================================================
## 1. Set working directory and create output folders
## =========================================================

setwd("/home/xxm_xxm/CJX_workspace/RPs/")

dir.create("results", showWarnings = FALSE)
dir.create("figures", showWarnings = FALSE)

## =========================================================
## 2. Load CPTAC COAD tumor/normal proteome log2FC matrix
## =========================================================

exp <- read.table(
  "Human__CPTAC_COAD__PNNL__Proteome__TMT__03_01_2017__BCM__Gene__Tumor_Normal_log2FC.cct",
  header = TRUE,
  row.names = 1,
  sep = "\t",
  check.names = FALSE
)

mat <- as.matrix(exp)

## =========================================================
## 3. Load NPM1 high-confidence interactors
## =========================================================

npm1_interactors <- read.csv(
  "../NPM1-BIOID2/NPM1_high_confidence_interactors.csv",
  stringsAsFactors = FALSE
)

gene_list <- unique(npm1_interactors$Gene)

## =========================================================
## 4. GO enrichment analysis
## =========================================================

gene_entrez <- clusterProfiler::bitr(
  gene_list,
  fromType = "SYMBOL",
  toType = "ENTREZID",
  OrgDb = org.Hs.eg.db
)

ego_bp <- clusterProfiler::enrichGO(
  gene          = gene_entrez$ENTREZID,
  OrgDb         = org.Hs.eg.db,
  keyType       = "ENTREZID",
  ont           = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.2,
  readable      = TRUE
)

ego_bp_simplified <- clusterProfiler::simplify(
  ego_bp,
  cutoff = 0.7,
  by = "p.adjust",
  select_fun = min
)

## =========================================================
## 5. Select representative GO terms as functional modules
## =========================================================

go_show <- ego_bp_simplified@result %>%
  dplyr::arrange(p.adjust) %>%
  dplyr::slice_head(n = 10)

write.csv(
  go_show,
  "results/NPM1_interactor_GO_modules_used_for_heatmap.csv",
  row.names = FALSE
)

## =========================================================
## 6. Build GO term-gene relationship
## =========================================================

go_gene_df <- go_show %>%
  dplyr::select(ID, Description, geneID, p.adjust, Count) %>%
  tidyr::separate_rows(geneID, sep = "/") %>%
  dplyr::rename(
    GO_ID = ID,
    Functional_module = Description,
    Gene = geneID
  )

## If one gene belongs to multiple GO terms,
## assign it to the most significant GO term.
gene_module <- go_gene_df %>%
  dplyr::group_by(Gene) %>%
  dplyr::slice_min(
    order_by = p.adjust,
    n = 1,
    with_ties = FALSE
  ) %>%
  dplyr::ungroup() %>%
  dplyr::select(Gene, Functional_module)

## =========================================================
## 7. Add GO-module annotation to NPM1 interactors
## =========================================================

module_df <- npm1_interactors %>%
  dplyr::left_join(gene_module, by = "Gene") %>%
  dplyr::mutate(
    Functional_module = ifelse(
      is.na(Functional_module),
      "Other NPM1 interactors",
      Functional_module
    )
  )

## =========================================================
## 8. Match NPM1 interactors to CPTAC proteome matrix
## =========================================================

matched_genes <- intersect(module_df$Gene, rownames(mat))
missing_genes <- setdiff(module_df$Gene, rownames(mat))

message("Matched NPM1 interactors in CPTAC: ", length(matched_genes))
message("Missing NPM1 interactors in CPTAC: ", length(missing_genes))

data <- mat[matched_genes, , drop = FALSE]
data[is.na(data) | is.nan(data) | is.infinite(data)] <- 0

## =========================================================
## 9. Build annotation dataframe
## =========================================================

anno_df <- module_df %>%
  dplyr::filter(Gene %in% matched_genes) %>%
  dplyr::distinct(
    Gene,
    Functional_module,
    log2FC_NPM1_VEC,
    X.Spec.NPM1
  )

anno_df <- anno_df[match(rownames(data), anno_df$Gene), ]

anno_df$CPTAC_mean_log2FC <- rowMeans(data, na.rm = TRUE)

## =========================================================
## 10. Remove unclassified proteins
## =========================================================

anno_df_focus <- anno_df %>%
  dplyr::filter(Functional_module != "Other NPM1 interactors")

## =========================================================
## 11. Summarize CPTAC changes by functional module
## =========================================================

module_summary <- anno_df_focus %>%
  dplyr::group_by(Functional_module) %>%
  dplyr::summarise(
    n = dplyr::n(),
    mean_log2FC = mean(CPTAC_mean_log2FC, na.rm = TRUE),
    median_log2FC = median(CPTAC_mean_log2FC, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::arrange(desc(mean_log2FC))

write.csv(
  module_summary,
  "results/NPM1_interactor_GO_module_CPTAC_summary.csv",
  row.names = FALSE
)

print(module_summary)

## =========================================================
## 12. Keep all classified genes
## =========================================================

anno_df_plot <- anno_df_focus

## =========================================================
## 13. Order functional modules and genes
## =========================================================

## Order modules by their average CPTAC tumor/normal log2FC
module_order <- module_summary$Functional_module

anno_df_plot$Functional_module <- factor(
  anno_df_plot$Functional_module,
  levels = module_order
)

## Within each module, order genes by CPTAC mean log2FC
row_order <- anno_df_plot %>%
  dplyr::arrange(
    Functional_module,
    desc(CPTAC_mean_log2FC)
  ) %>%
  dplyr::pull(Gene)

data_ordered <- data[row_order, , drop = FALSE]

anno_ordered <- anno_df_plot[
  match(row_order, anno_df_plot$Gene),
]

## =========================================================
## 14. Limit extreme values for visualization
## =========================================================

data_plot <- data_ordered
data_plot[data_plot > 2] <- 2
data_plot[data_plot < -2] <- -2

## =========================================================
## 15. Define heatmap color scale
## =========================================================

col_fun <- circlize::colorRamp2(
  c(-2, -1, 0, 1, 2),
  c("#3B4CC0", "#AFCBFF", "white", "#FDBE85", "#B40426")
)

## =========================================================
## 16. Define functional module colors
## =========================================================

module_levels <- levels(anno_ordered$Functional_module)

module_palette <- c(
  "#E64B35FF",
  "#4DBBD5FF",
  "#00A087FF",
  "#3C5488FF",
  "#F39B7FFF",
  "#8491B4FF",
  "#91D1C2FF",
  "#7E6148FF",
  "#B09C85FF",
  "#6A3D9AFF"
)

module_cols <- module_palette[seq_along(module_levels)]
names(module_cols) <- module_levels

## =========================================================
## 17. Row annotations
## =========================================================

left_anno <- rowAnnotation(
  Module = anno_ordered$Functional_module,
  col = list(Module = module_cols),
  show_annotation_name = FALSE
)

right_anno <- rowAnnotation(
  `Mean\nlog2FC` = anno_barplot(
    anno_ordered$CPTAC_mean_log2FC,
    gp = gpar(fill = "#B40426", col = NA),
    border = FALSE,
    width = unit(2.2, "cm")
  ),
  annotation_name_gp = gpar(fontsize = 8)
)

## =========================================================
## 18. Draw heatmap
## =========================================================

ht <- Heatmap(
  data_plot,
  name = "Tumor/Normal\nlog2FC",
  col = col_fun,
  
  cluster_rows = FALSE,
  cluster_columns = TRUE,
  
  row_split = anno_ordered$Functional_module,
  row_title_gp = gpar(fontsize = 8, fontface = "bold"),
  
  show_column_names = FALSE,
  row_names_side = "right",
  row_names_gp = gpar(fontsize = 6),
  
  width = unit(5, "cm"),
  height = unit(12, "cm")
)

ht_final <- left_anno + ht + right_anno

pdf(
  "figures/Fig_NPM1_interactor_GO_module_CPTAC_heatmap.pdf",
  width = 9,
  height = 12,
  useDingbats = FALSE
)

draw(ht_final)

dev.off()

png(
  "figures/Fig_NPM1_interactor_GO_module_CPTAC_heatmap.png",
  width = 2700,
  height = 3600,
  res = 300
)

draw(ht_final)

dev.off()

## =========================================================
## 19. Export ordered gene list
## =========================================================

module_gene_list <- anno_ordered %>%
  dplyr::select(
    Functional_module,
    Gene,
    CPTAC_mean_log2FC,
    log2FC_NPM1_VEC,
    X.Spec.NPM1
  )

write.csv(
  module_gene_list,
  "results/NPM1_interactor_GO_module_CPTAC_ordered_gene_list.csv",
  row.names = FALSE
)

## =========================================================
## 20. Save processed objects
## =========================================================

save(
  npm1_interactors,
  ego_bp,
  ego_bp_simplified,
  go_show,
  module_df,
  anno_df,
  anno_df_focus,
  module_summary,
  module_gene_list,
  file = "results/NPM1_interactor_GO_module_CPTAC_analysis.RData"
)

message("Analysis completed.")
message("Displayed classified genes: ", nrow(anno_ordered))
message("Functional modules: ", length(unique(anno_ordered$Functional_module)))
