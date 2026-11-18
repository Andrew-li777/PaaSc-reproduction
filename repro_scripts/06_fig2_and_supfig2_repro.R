# ==============================================================================
# PaaSc Paper Main Figure 2 & Supplementary Figure 2A Reproduction
# 论文主图 Fig 2（多方法基准、降采样鲁棒性、CITE-seq蛋白验证）与附图 SupFig 2A（回归t值去偏原理）
# 
# 数据来源：
# 1. SupFig 2A: 基于 data/pbmc3k_filtered_gene_bc_matrices/hg19 现场重算 doRegression
# 2. Fig 2B: zenodo_data/PaaSc_data/Fig2/2B.txt (8种算法 x 7种PBMC亚型 AUC)
# 3. Fig 2C: zenodo_data/PaaSc_data/Fig2/2C.txt (20%-80% 降采样 Dropout 鲁棒性)
# 4. Fig 2D: zenodo_data/PaaSc_data/Fig2/2D.txt (CITE-seq 蛋白真值 F1/Precision/Recall)
#
# 输出文件：
# - repro_results/fig3_to_fig9/Fig2_and_SupFig2_benchmark.png
# - repro_results/fig3_to_fig9/fig2_verification_summary.txt
# ==============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(PaaSc)
  library(tidyverse)
  library(patchwork)
})

base_dir <- "d:/project/gemini/reproduction/PaaSc"
data_dir <- file.path(base_dir, "data")
zen_dir  <- file.path(base_dir, "zenodo_data", "PaaSc_data")
out_dir  <- file.path(base_dir, "repro_results", "fig3_to_fig9")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

cat("=== 1. Calculating SupFig 2A: MCA Regression t-values on pbmc3k ===\n")
pbmc_mat <- Read10X(file.path(data_dir, "pbmc3k_filtered_gene_bc_matrices", "hg19"))
pbmc_obj <- CreateSeuratObject(counts = pbmc_mat, min.cells = 3, min.features = 200)
pbmc_obj <- computeMCA(pbmc_obj, nmcs = 20)

bg_gmt <- file.path(data_dir, "kegg_legacy.gmt")
bcell_pt <- file.path(data_dir, "B_cell_activation.txt")
gene_rate <- getGeneRate(background.geneset = bg_gmt, pathway.geneset = list(B_cell_activation = readLines(bcell_pt)))
rdata <- doRegression(pbmc_obj, gene.rate = gene_rate)

n_dims <- ncol(rdata)
df_s2a <- tibble(
  Dimension = rep(seq_len(n_dims), 2),
  DimName   = rep(colnames(rdata), 2),
  tvalue    = c(as.numeric(rdata["B_cell_activation.t", ]), as.numeric(rdata["B_cell_activation.background.t", ])),
  pvalue    = c(as.numeric(rdata["B_cell_activation.pvalue", ]), as.numeric(rdata["B_cell_activation.background.pvalue", ])),
  Group     = c(rep("Pathway (B-cell)", n_dims), rep("Background (KEGG)", n_dims))
)

p_2a <- ggplot(df_s2a, aes(x = factor(Dimension), y = tvalue, fill = Group)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.75, alpha = 0.9) +
  geom_hline(yintercept = c(-1.96, 1.96), linetype = "dashed", color = "grey50", linewidth = 0.6) +
  scale_fill_manual(values = c("Pathway (B-cell)" = "#E64B35", "Background (KEGG)" = "#7F8C8D")) +
  labs(
    title = "SupFig 2A: PaaSc Linear Regression t-values Across MCA Dimensions",
    subtitle = "De-biasing pathway loadings against background gene frequencies (pbmc3k)",
    x = "MCA Dimension",
    y = "Regression t-statistic",
    fill = NULL
  ) +
  theme_minimal(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(color = "grey35", hjust = 0.5, size = 9),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
    panel.grid.minor = element_blank(),
    legend.position = "top"
  )

cat("=== 2. Loading & Plotting Fig 2B: Multi-Method AUC on 7 PBMC Cell Types ===\n")
f2b <- read.table(file.path(zen_dir, "Fig2", "2B.txt"), header = TRUE, check.names = FALSE)
colnames(f2b)[1] <- "Method"

df_2b <- f2b %>%
  pivot_longer(-Method, names_to = "CellType", values_to = "AUC") %>%
  mutate(Method = factor(Method, levels = c("PaaSc", "CelliD", "GSDensity", "AUCell", "VAM", "GSEA", "GSVA", "scGSEA")))

method_cols <- c(
  "PaaSc" = "#E89DA0", "CelliD" = "#E6CECF", "GSDensity" = "#F6C8A8", "AUCell" = "#9FBA95",
  "VAM" = "#B2D3A4", "GSEA" = "#80C1C4", "GSVA" = "#B696B6", "scGSEA" = "#88CEE6"
)

p_2b <- ggplot(df_2b, aes(x = Method, y = AUC, fill = Method)) +
  geom_boxplot(alpha = 0.8, outlier.shape = NA, width = 0.6) +
  geom_jitter(width = 0.15, size = 2.2, alpha = 0.8) +
  scale_fill_manual(values = method_cols) +
  coord_cartesian(ylim = c(0.5, 1.02)) +
  labs(
    title = "Fig 2B: PBMC Cell-Type Identification AUC (7 Subtypes)",
    x = NULL, y = "ROC AUC"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    axis.text.x = element_text(angle = 40, hjust = 1, face = "bold", color = "black"),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
    legend.position = "none"
  )

cat("=== 3. Loading & Plotting Fig 2C: Dropout Subsampling Robustness (20% to 100%) ===\n")
f2c <- read.table(file.path(zen_dir, "Fig2", "2C.txt"), header = TRUE)
df_2c <- f2c %>%
  mutate(
    Percent = case_when(
      Group.2 == "full.marker.set" ~ 100,
      TRUE ~ as.numeric(sub("_percent.ds", "", Group.2))
    ),
    Method = case_when(
      tolower(Group.1) == "paasc" ~ "PaaSc",
      tolower(Group.1) == "cellid" ~ "CelliD",
      tolower(Group.1) == "gsdensity" ~ "GSDensity",
      tolower(Group.1) == "aucell" ~ "AUCell",
      tolower(Group.1) == "vam" ~ "VAM",
      tolower(Group.1) == "gsea" ~ "GSEA",
      tolower(Group.1) == "gsva" ~ "GSVA",
      tolower(Group.1) == "scgsea" ~ "scGSEA",
      TRUE ~ as.character(Group.1)
    )
  )

p_2c <- ggplot(df_2c, aes(x = Percent, y = x, color = Method, group = Method)) +
  geom_line(linewidth = 1.0) +
  geom_point(size = 2.2) +
  scale_color_manual(values = method_cols) +
  scale_x_continuous(breaks = seq(20, 100, by = 10), labels = paste0(seq(20, 100, by = 10), "%")) +
  labs(
    title = "Fig 2C: Downsampling Robustness Across Dropout Levels",
    x = "Subsampled Sequencing Depth",
    y = "Mean AUC"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
    panel.grid.minor = element_blank()
  )

cat("=== 4. Loading & Plotting Fig 2D: CITE-seq Ground Truth Benchmark ===\n")
f2d <- read.table(file.path(zen_dir, "Fig2", "2D.txt"), header = TRUE, sep = "\t")
df_2d <- f2d %>%
  filter(metrics == "F1") %>%
  mutate(
    Method = case_when(
      methods == "GDSdensity" ~ "GSDensity",
      methods == "CellID"     ~ "CelliD",
      methods == "ssGSEA"     ~ "GSEA",
      TRUE ~ as.character(methods)
    ),
    Method = factor(Method, levels = c("PaaSc", "CelliD", "GSDensity", "AUCell", "VAM", "GSEA", "GSVA", "scGSEA"))
  )

p_2d <- ggplot(df_2d, aes(x = Method, y = value, fill = Method)) +
  geom_boxplot(alpha = 0.8, outlier.shape = NA, width = 0.55) +
  geom_jitter(width = 0.15, size = 2.2, alpha = 0.8) +
  scale_fill_manual(values = method_cols) +
  labs(
    title = "Fig 2D: CITE-seq Surface Protein Ground Truth (F1 Score)",
    x = NULL, y = "F1 Score"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    axis.text.x = element_text(angle = 35, hjust = 1, face = "bold", color = "black"),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
    legend.position = "none"
  )

p_all <- (p_2a | p_2b) / (p_2c | p_2d) +
  plot_annotation(
    title = "PaaSc Paper Figure 2 & Supplementary Figure 2A Comprehensive Benchmark",
    subtitle = "Methodological Comparison, Subsampling Robustness, CITE-seq Validation, and Regression De-biasing",
    theme = theme(
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 10, hjust = 0.5, color = "grey35")
    )
  )

fig_path <- file.path(out_dir, "Fig2_and_SupFig2_benchmark.png")
ggsave(fig_path, p_all, width = 14, height = 10, dpi = 300, bg = "white")
cat("Saved composite plot to:", fig_path, "\n")

# Summary text output
sink(file.path(out_dir, "fig2_verification_summary.txt"))
cat("PaaSc Figure 2 & SupFig 2A Verification Summary\n")
cat("==============================================\n")
cat("SupFig 2A (pbmc3k doRegression):\n")
cat("  Evaluated MCA Dimensions:", n_dims, "\n")
cat("  Max Pathway t-value:", sprintf("%.4f", max(df_s2a$tvalue[df_s2a$Group == "Pathway (B-cell)"])), "\n")
cat("  Significant Dimensions (|t| > 1.96):", sum(abs(df_s2a$tvalue[df_s2a$Group == "Pathway (B-cell)"]) > 1.96), "/", n_dims, "\n\n")

cat("Fig 2B (PBMC 7 Subtypes AUC):\n")
print(f2b)
cat("\nMean AUC by Method across 7 subtypes:\n")
mean_auc_2b <- df_2b %>% group_by(Method) %>% summarise(Mean_AUC = mean(AUC), .groups = "drop") %>% arrange(desc(Mean_AUC))
print(as.data.frame(mean_auc_2b), row.names = FALSE)

cat("\nFig 2C (Downsampling Robustness at 20% vs 80%):\n")
f2c_summary <- df_2c %>% filter(Percent %in% c(20, 80)) %>% arrange(Percent, desc(x))
print(as.data.frame(f2c_summary %>% select(Method, Percent, Mean_AUC = x)), row.names = FALSE)

cat("\nFig 2D (CITE-seq Ground Truth F1 Score):\n")
f2d_summary <- df_2d %>% group_by(methods) %>% summarise(Mean_F1 = mean(value), Median_F1 = median(value), .groups = "drop") %>% arrange(desc(Mean_F1))
print(as.data.frame(f2d_summary), row.names = FALSE)
sink()

cat("=== Fig2 & SupFig2 Reproduction Complete! ===\n")
