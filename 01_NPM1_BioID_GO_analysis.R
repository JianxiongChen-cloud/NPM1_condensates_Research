## =========================================================
## NPM1 BioID-MS analysis
## Identify NPM1-associated proteins, perform GO enrichment,
## visualize functional modules, and extract ribosomal proteins
## =========================================================

## 0. Set working directory
setwd("/home/xxm_xxm/CJX_workspace/NPM1-BIOID2/")

## 1. Load packages
suppressPackageStartupMessages({
  library(tidyverse)
  library(stringr)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(enrichplot)
  library(igraph)
  library(ggraph)
  library(ggplot2)
  library(ggrepel)
})

## 2. Create output folders
dir.create("results", showWarnings = FALSE)
dir.create("figures", showWarnings = FALSE)

## 3. Read protein table
df <- read.csv("proteins.csv", check.names = FALSE)

## 4. Extract gene symbols and calculate enrichment
df <- df %>%
  dplyr::mutate(
    Gene = stringr::str_extract(Description, "(?<=GN=)[^ ]+"),
    log2FC_NPM1_VEC = log2((`Area NPM1` + 1) / (`Area VEC` + 1)),
    log2FC_RPL_VEC  = log2((`Area RPL` + 1) / (`Area VEC` + 1))
  )

## 5. Remove common contaminants
contaminants <- c(
  "KRT1", "KRT2", "KRT9", "KRT10", "KRT14",
  "ACTB", "ACTG1",
  "TUBA1A", "TUBA1B", "TUBB",
  "HSPA1A", "HSPA1B", "HSP90AA1"
)

df_clean <- df %>%
  dplyr::filter(!Gene %in% contaminants) %>%
  dplyr::filter(!is.na(Gene))

## 6. Define high-confidence NPM1 interactors
npm1_interactors <- df_clean %>%
  dplyr::filter(
    log2FC_NPM1_VEC > 1,
    `#Unique` >= 2,
    `#Spec NPM1` >= 5,
    `#Spec NPM1` > `#Spec VEC`
  ) %>%
  dplyr::arrange(desc(log2FC_NPM1_VEC))

write.csv(
  npm1_interactors,
  "results/NPM1_high_confidence_interactors.csv",
  row.names = FALSE
)

## 7. Export top candidates
top_candidates <- npm1_interactors %>%
  dplyr::select(
    Gene, Accession, log2FC_NPM1_VEC,
    `#Peptides`, `#Unique`,
    `#Spec NPM1`, `#Spec VEC`,
    Description
  ) %>%
  dplyr::slice_head(n = 30)

write.csv(
  top_candidates,
  "results/Top30_NPM1_interactors.csv",
  row.names = FALSE
)

## 8. Scatter plot for NPM1 enrichment
plot_df <- df_clean %>%
  dplyr::mutate(
    group = dplyr::case_when(
      Gene %in% npm1_interactors$Gene ~ "High-confidence interactors",
      TRUE ~ "Other proteins"
    )
  )

p_scatter <- ggplot(
  plot_df,
  aes(x = log2FC_NPM1_VEC, y = `#Spec NPM1`)
) +
  geom_point(aes(size = `#Unique`, color = group), alpha = 0.75) +
  geom_vline(xintercept = 1, linetype = "dashed") +
  scale_color_manual(
    values = c(
      "High-confidence interactors" = "#E64B35FF",
      "Other proteins" = "grey75"
    )
  ) +
  theme_classic(base_size = 14) +
  labs(
    x = "log2FC NPM1 / VEC",
    y = "Spectral count in NPM1 BioID",
    color = "",
    size = "Unique peptides"
  )

ggsave(
  "figures/Fig_IPMS_NPM1_enrichment_scatter.pdf",
  p_scatter,
  width = 6,
  height = 5,
  device = cairo_pdf
)

ggsave(
  "figures/Fig_IPMS_NPM1_enrichment_scatter.png",
  p_scatter,
  width = 6,
  height = 5,
  dpi = 600
)

## 9. GO enrichment analysis
gene_list <- unique(npm1_interactors$Gene)

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

write.csv(
  as.data.frame(ego_bp),
  "results/GO_BP_NPM1_interactors.csv",
  row.names = FALSE
)

## 10. Remove redundant GO terms
ego_bp_simplified <- clusterProfiler::simplify(
  ego_bp,
  cutoff = 0.7,
  by = "p.adjust",
  select_fun = min
)

write.csv(
  as.data.frame(ego_bp_simplified),
  "results/GO_BP_NPM1_interactors_simplified.csv",
  row.names = FALSE
)

## 11. GO enrichment map
ego_bp_focus <- ego_bp_simplified
ego_bp_focus@result <- ego_bp_simplified@result %>%
  dplyr::arrange(p.adjust) %>%
  dplyr::slice_head(n = 15)

ego_bp_focus <- enrichplot::pairwise_termsim(ego_bp_focus)

p_emap <- enrichplot::emapplot(
  ego_bp_focus,
  showCategory = 15,
  color = "p.adjust",
  layout = "kk",
  label_format = 28
) +
  scale_color_gradientn(
    colors = c("#D95F5F", "#F2B6A0", "#D9E6F2", "#4C78A8"),
    name = "Adjusted\nP value"
  ) +
  scale_size_continuous(
    range = c(3.5, 8),
    name = "Gene count"
  ) +
  theme_void(base_size = 12) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 13),
    legend.position = "right",
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 9),
    plot.margin = margin(10, 10, 10, 10)
  ) +
  labs(
    title = "Functional modules of NPM1-associated proteins"
  )

ggsave(
  "figures/Fig_GO_BP_emapplot_NPM1_interactors.pdf",
  p_emap,
  width = 7,
  height = 5.8,
  device = cairo_pdf
)

ggsave(
  "figures/Fig_GO_BP_emapplot_NPM1_interactors.png",
  p_emap,
  width = 7,
  height = 5.8,
  dpi = 600
)

## 12. GO term-gene network
go_show <- ego_bp_simplified@result %>%
  dplyr::arrange(p.adjust) %>%
  dplyr::slice_head(n = 10)

edge_df <- go_show %>%
  dplyr::select(ID, Description, geneID, p.adjust, Count) %>%
  tidyr::separate_rows(geneID, sep = "/") %>%
  dplyr::rename(
    GO_ID = ID,
    GO_term = Description,
    Gene = geneID
  )

gene_module <- edge_df %>%
  dplyr::group_by(Gene) %>%
  dplyr::slice_min(order_by = p.adjust, n = 1, with_ties = FALSE) %>%
  dplyr::ungroup() %>%
  dplyr::select(Gene, module = GO_term)

key_genes <- npm1_interactors %>%
  dplyr::arrange(desc(log2FC_NPM1_VEC)) %>%
  dplyr::slice_head(n = 80) %>%
  dplyr::pull(Gene)

go_nodes <- go_show %>%
  dplyr::transmute(
    name = Description,
    type = "GO term",
    module = Description,
    size = Count,
    pvalue_score = -log10(p.adjust),
    label = Description
  )

gene_nodes <- edge_df %>%
  dplyr::distinct(Gene) %>%
  dplyr::left_join(gene_module, by = "Gene") %>%
  dplyr::left_join(
    npm1_interactors %>%
      dplyr::select(Gene, log2FC_NPM1_VEC, `#Spec NPM1`),
    by = "Gene"
  ) %>%
  dplyr::transmute(
    name = Gene,
    type = "Gene",
    module = module,
    size = ifelse(is.na(`#Spec NPM1`), 3, `#Spec NPM1`),
    pvalue_score = NA_real_,
    label = ifelse(Gene %in% key_genes, Gene, "")
  )

nodes <- dplyr::bind_rows(go_nodes, gene_nodes)

edges <- edge_df %>%
  dplyr::transmute(
    from = GO_term,
    to = Gene
  )

g <- igraph::graph_from_data_frame(
  d = edges,
  vertices = nodes,
  directed = FALSE
)

module_cols <- c(
  "#4DBBD5FF",
  "#00A087FF",
  "#3C5488FF",
  "#8491B4FF",
  "#91D1C2FF",
  "#7E6148FF",
  "#B09C85FF",
  "#6A3D9AFF",
  "#1F78B4",
  "#33A02C"
)

names(module_cols) <- unique(go_show$Description)[seq_along(module_cols)]

p_go_gene_net <- ggraph(g, layout = "fr") +
  geom_edge_link(
    color = "grey75",
    linewidth = 0.3,
    alpha = 0.3
  ) +
  geom_node_point(
    data = function(x) dplyr::filter(x, type == "Gene"),
    aes(size = size, fill = module),
    shape = 21,
    color = "white",
    stroke = 0.15,
    alpha = 0.9
  ) +
  geom_node_point(
    data = function(x) dplyr::filter(x, type == "GO term"),
    aes(size = size, color = pvalue_score),
    shape = 16,
    alpha = 0.95
  ) +
  scale_fill_manual(
    values = module_cols,
    name = "Functional module"
  ) +
  scale_color_gradient(
    low = "#FCAE91",
    high = "#CB181D",
    name = "-log10(adj. P)"
  ) +
  scale_size_continuous(
    range = c(1.5, 8),
    name = "Size"
  ) +
  geom_node_text(
    aes(label = label),
    repel = TRUE,
    size = 2.6,
    color = "black",
    max.overlaps = 150,
    segment.size = 0.2,
    segment.alpha = 0.4
  ) +
  theme_void(base_size = 12) +
  theme(
    legend.position = "right",
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 8),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 13),
    plot.margin = margin(10, 10, 10, 10)
  ) +
  labs(
    title = "GO term-gene network of NPM1-associated proteins"
  )

ggsave(
  "figures/Fig_GO_gene_network_NPM1_interactors.pdf",
  p_go_gene_net,
  width = 9,
  height = 7,
  units = "in",
  device = cairo_pdf
)

ggsave(
  "figures/Fig_GO_gene_network_NPM1_interactors.png",
  p_go_gene_net,
  width = 9,
  height = 7,
  dpi = 600
)

## 13. Extract ribosomal proteins
ribosomal_gene_set <- c(
  "RPSA", "RPS9", "RPS8", "RPS7", "RPS6KA3", "RPS6", "RPS5", "RPS4X",
  "RPS3A", "RPS3", "RPS29", "RPS28", "RPS27", "RPS26", "RPS25",
  "RPS24", "RPS23", "RPS21", "RPS20", "RPS2", "RPS19BP1", "RPS19",
  "RPS18", "RPS17", "RPS16", "RPS15A", "RPS15", "RPS14", "RPS13",
  "RPS12", "RPS11", "RPS10",
  "RPLP2", "RPLP1", "RPLP0", "RPL9P9", "RPL8", "RPL7L1", "RPL7A",
  "RPL7", "RPL6", "RPL5", "RPL4", "RPL39P5", "RPL39", "RPL38",
  "RPL37A", "RPL37", "RPL36A", "RPL36", "RPL35A", "RPL35", "RPL34",
  "RPL32", "RPL31", "RPL30", "RPL3", "RPL29", "RPL28", "RPL27A",
  "RPL27", "RPL26L1", "RPL26", "RPL24", "RPL23A", "RPL23", "RPL22",
  "RPL21", "RPL19", "RPL18A", "RPL18", "RPL17", "RPL15", "RPL14",
  "RPL13A", "RPL13", "RPL12", "RPL11", "RPL10A", "RPL10",
  "RP9", "UBA52", "FAU"
)

ribosomal_subset <- npm1_interactors %>%
  dplyr::filter(Gene %in% ribosomal_gene_set) %>%
  dplyr::arrange(desc(log2FC_NPM1_VEC))

write.csv(
  ribosomal_subset,
  "results/NPM1_interacting_ribosomal_proteins.csv",
  row.names = FALSE
)

npm1_interactors <- npm1_interactors %>%
  dplyr::mutate(
    Ribosome = ifelse(Gene %in% ribosomal_gene_set, "Ribosomal protein", "Other")
  )

top_rp <- npm1_interactors %>%
  dplyr::filter(Ribosome == "Ribosomal protein") %>%
  dplyr::slice_max(order_by = log2FC_NPM1_VEC, n = 10)

p_rp <- ggplot(
  npm1_interactors,
  aes(x = log2FC_NPM1_VEC, y = `#Spec NPM1`)
) +
  geom_point(color = "grey80", alpha = 0.6, size = 2.2) +
  geom_point(
    data = dplyr::filter(npm1_interactors, Ribosome == "Ribosomal protein"),
    color = "#E64B35",
    size = 2.6
  ) +
  geom_text_repel(
    data = top_rp,
    aes(label = Gene),
    size = 3,
    max.overlaps = 50
  ) +
  theme_classic(base_size = 13) +
  labs(
    x = "log2FC NPM1 / VEC",
    y = "Spectral count in NPM1 BioID",
    title = "Ribosomal proteins are enriched in the NPM1 interactome"
  )

ggsave(
  "figures/Fig_NPM1_ribosomal_proteins_scatter.pdf",
  p_rp,
  width = 6,
  height = 5,
  device = cairo_pdf
)

ggsave(
  "figures/Fig_NPM1_ribosomal_proteins_scatter.png",
  p_rp,
  width = 6,
  height = 5,
  dpi = 600
)

## 14. Save processed R objects
save(
  df,
  df_clean,
  npm1_interactors,
  ego_bp,
  ego_bp_simplified,
  ribosomal_subset,
  file = "results/NPM1_BioID_GO_analysis_results.RData"
)

## 15. Print summary
cat("Analysis completed.\n")
cat("High-confidence NPM1 interactors:", nrow(npm1_interactors), "\n")
cat("Ribosomal proteins among interactors:", length(unique(ribosomal_subset$Gene)), "\n")
