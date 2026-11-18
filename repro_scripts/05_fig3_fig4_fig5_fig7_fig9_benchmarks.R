# ==============================================================================
# PaaSc 论文多主图跨模态综合基准复现套件: Fig 3, Fig 4D-E, Fig 5, Fig 7B, Fig 9
# 
# 文献出处: Liao X, Hong Y, Feng Y, et al. PLOS Comput Biol 2025, 21(11): e1013666.
# 
# 模块划分与科学场景:
#   1. Part 1 (Fig 3A-B): TISCH 泛癌 5 大免疫谱系跨 225 组队列-谱系评估 (8 种方法 Rank-1 胜出对比)
#   2. Part 2 (Fig 4D-E): 3 套体外/体内细胞衰老 scRNA-seq 队列 ROC 验证 + 泛癌患者水平 Hallmark 协同热图
#   3. Part 3 (Fig 5A-D): Tabula Muris Senis (10,725 细胞) x 75 UKB/GWAS 多基因遗传风险性状组织特异富集
#   4. Part 4 (Fig 7B): PBMC scATAC-seq (10,459 细胞) Gene-Activity 16 种免疫亚型识别基准
#   5. Part 5 (Fig 9B/D/E): NanoString CosMx SMI NSCLC 空间转录/蛋白单细胞 (98,002 细胞) 微环境分区与空间自相关
#   6. Part 6: 导出汇总量化验收指标表 (fig3_to_fig9_summary.tsv)
#
# 数据来源说明 (诚实边界):
#   - 上游数百个泛癌原始矩阵、TMS .h5ad 及 CosMx 原始矩阵未随 Zenodo PaaSc_data.zip 下发；
#   - 本脚本基于作者 Zenodo 归档的中间打分矩阵与 AUC 表，完整复现全部统计检验与高品质出版图表。
#
# 统计与工程细节:
#   - Fig 5C: 引入 Haldane-Anscombe +0.5 连续性平滑，解决全阳性细胞群导致 Odds Ratio = Inf 无法连续绘图的问题
#   - Fig 9E: 对 98,002 细胞大矩阵做向量化 Pearson 相关与 pt() t-检验过滤，并用 arrayInd() 精确回溯极大值对
#   - 排版优化: Fig 3A 采用 ComplexHeatmap packLegend 联合图例对齐；Fig 4E / Fig 7B 左边距充裕无文字截断
#
# 预期耗时: 约 40-50 秒 (单线程)
# ==============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)
  library(ComplexHeatmap)
  library(circlize)
  library(pROC)
  library(ROCit)
  library(pheatmap)
  library(ggh4x)
  library(patchwork)
  library(grid)
  library(PaaSc)
})

base_dir <- "d:/project/gemini/reproduction/PaaSc"
zen_dir  <- file.path(base_dir, "zenodo_data", "PaaSc_data")
out_dir  <- file.path(base_dir, "repro_results", "fig3_to_fig9")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

summary_rows <- list()
add_summary <- function(fig, panel, metric, value) {
  summary_rows[[length(summary_rows) + 1]] <<- tibble(
    Figure = fig, Panel = panel, Metric = metric, Value = as.character(value)
  )
}

# ==============================================================================
# PART 1: FIGURE 3 - TISCH Pan-Cancer Immune Lineage AUC Heatmaps & Summary
# ==============================================================================
cat("=== Part 1: Figure 3 (TISCH Pan-Cancer Lineage Benchmark) ===\n")

fig3_files <- list(
  CD8T  = list(file = "CD8T_auc_values.txt",  title = "CD8 T",      levels = c("PaaSc","scGSEA","CelliD","GSDensity","AUCell","GSEA","VAM","GSVA")),
  Macro = list(file = "Macro_auc_values.txt", title = "Macrophage", levels = c("PaaSc","scGSEA","GSDensity","VAM","AUCell","CelliD","GSEA","GSVA")),
  NK    = list(file = "NK_auc_values.txt",    title = "NK",         levels = c("PaaSc","scGSEA","GSDensity","GSEA","AUCell","CelliD","VAM","GSVA")),
  DC    = list(file = "DC_auc_values.txt",    title = "DC",         levels = c("PaaSc","scGSEA","GSDensity","CelliD","AUCell","GSEA","VAM","GSVA")),
  Treg  = list(file = "Treg_auc_values.txt",  title = "Treg",       levels = c("PaaSc","CelliD","VAM","GSVA","scGSEA","GSDensity","GSEA","AUCell"))
)

custom_colors <- setNames(
  c("#E89DA0", "#88CEE6", "#F6C8A8", "#B2D3A4", "#9FBA95", "#E6CECF", "#B696B6", "#80C1C4"),
  c("PaaSc", "scGSEA", "GSDensity", "VAM", "AUCell", "CelliD", "GSVA", "GSEA")
)

fig3_all_long <- list()
fig3_grobs <- list()

for (ct in names(fig3_files)) {
  spec <- fig3_files[[ct]]
  df <- read.table(file.path(zen_dir, "Fig3", spec$file), header = TRUE, sep = "\t", check.names = FALSE)
  if ("GDS" %in% colnames(df)) colnames(df)[colnames(df) == "GDS"] <- "GSDensity"
  if ("ssGSEA" %in% colnames(df)) colnames(df)[colnames(df) == "ssGSEA"] <- "GSEA"

  method_cols <- c("PaaSc","scGSEA","CelliD","GSDensity","AUCell","GSEA","VAM","GSVA")
  mat <- as.matrix(df[, spec$levels])
  best_method <- factor(spec$levels[apply(mat, 1, which.max)], levels = spec$levels)

  long_df <- df %>%
    rownames_to_column("Dataset") %>%
    pivot_longer(cols = all_of(method_cols), names_to = "Method", values_to = "AUC") %>%
    mutate(Celltype = spec$title, BestMethod = rep(as.character(best_method), each = length(method_cols)))
  fig3_all_long[[ct]] <- long_df

  # Record stats
  wins <- sum(best_method == "PaaSc")
  add_summary("Fig3", spec$title, "N_cohorts", nrow(df))
  add_summary("Fig3", spec$title, "PaaSc_mean_AUC", sprintf("%.4f", mean(df$PaaSc)))
  add_summary("Fig3", spec$title, "PaaSc_top1_wins", paste0(wins, "/", nrow(df)))

  row_anno <- rowAnnotation(
    Best = best_method,
    col = list(Best = custom_colors),
    show_annotation_name = FALSE,
    simple_anno_size = unit(0.45, "cm"),
    show_legend = FALSE
  )
  ht <- Heatmap(
    mat,
    name = "AUC",
    column_title = paste0(spec$title, " (n=", nrow(df), ")"),
    column_title_gp = gpar(fontsize = 11, fontface = "bold"),
    row_split = best_method,
    row_title = NULL,
    col = colorRamp2(c(0.5, 0.8, 1.0), c("#6395C7", "white", "firebrick3")),
    cluster_columns = FALSE,
    cluster_rows = TRUE,
    cluster_row_slices = FALSE,
    show_row_names = FALSE,
    show_row_dend = FALSE,
    show_column_dend = FALSE,
    show_parent_dend_line = FALSE,
    column_names_rot = 45,
    column_names_gp = gpar(fontsize = 9),
    right_annotation = row_anno,
    show_heatmap_legend = FALSE
  )
  fig3_grobs[[ct]] <- wrap_elements(full = grid.grabExpr(
    draw(ht, padding = unit(c(2, 2, 2, 2), "mm"))
  ))
}

lgd_auc <- Legend(
  col_fun = colorRamp2(c(0.5, 0.8, 1.0), c("#6395C7", "white", "firebrick3")),
  title = "ROC AUC", direction = "horizontal"
)
lgd_best <- Legend(
  labels = names(custom_colors), legend_gp = gpar(fill = unname(custom_colors)),
  title = "Best Method", nrow = 1
)
lgd_grob <- wrap_elements(full = grid.grabExpr(
  draw(packLegend(lgd_auc, lgd_best, direction = "horizontal", column_gap = unit(8, "mm")))
))

p_fig3a <- (lgd_grob / wrap_plots(fig3_grobs, ncol = 5)) +
  plot_layout(heights = c(0.14, 1)) +
  plot_annotation(
    title = "Fig 3A: Pan-Cancer TISCH scRNA-seq Cell-Lineage Scoring AUC Heatmaps (8 Methods x 5 Lineages)",
    theme = theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5))
  )
ggsave(file.path(out_dir, "Fig3A_pancancer_auc_heatmaps.png"), p_fig3a, width = 20, height = 7.8, dpi = 300, bg = "white")

fig3_df <- bind_rows(fig3_all_long) %>%
  mutate(Celltype = factor(Celltype, levels = c("CD8 T", "Macrophage", "NK", "DC", "Treg")))

method_order <- fig3_df %>% group_by(Method) %>% summarise(med = median(AUC), .groups = "drop") %>% arrange(desc(med)) %>% pull(Method)
fig3_df$Method <- factor(fig3_df$Method, levels = method_order)

p_fig3_box <- ggplot(fig3_df, aes(x = Method, y = AUC, fill = Method)) +
  geom_boxplot(outlier.size = 0.7, alpha = 0.85) +
  facet_wrap(~Celltype, ncol = 5) +
  scale_fill_manual(values = custom_colors) +
  labs(title = "Per-Lineage ROC AUC Distribution Across Pan-Cancer Cohorts", x = NULL, y = "ROC AUC") +
  theme_minimal(base_size = 11) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 9, color = "black"),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.7),
    legend.position = "none"
  )

win_df <- fig3_df %>%
  distinct(Celltype, Dataset, BestMethod) %>%
  count(Celltype, BestMethod, name = "Wins") %>%
  group_by(BestMethod) %>%
  mutate(TotalWins = sum(Wins)) %>%
  ungroup() %>%
  mutate(BestMethod = factor(BestMethod, levels = method_order))

n_total_evals <- nrow(distinct(fig3_df, Celltype, Dataset))
p_fig3_win <- ggplot(win_df, aes(x = BestMethod, y = Wins, fill = Celltype)) +
  geom_col(width = 0.65) +
  geom_text(aes(label = TotalWins, y = TotalWins), stat = "unique", vjust = -0.35, fontface = "bold", size = 3.8) +
  scale_fill_brewer(palette = "Set2") +
  expand_limits(y = max(win_df$TotalWins) * 1.12) +
  labs(title = paste0("Total #1 Ranked Cohorts per Method (across ", n_total_evals, " lineage-cohort evaluations)"),
       x = NULL, y = "#1 Rank Count", fill = "Lineage") +
  theme_minimal(base_size = 11) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.text.x = element_text(size = 10.5, face = "bold", color = "black"),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.7)
  )

p_fig3b <- (p_fig3_box / p_fig3_win) + plot_layout(heights = c(1.2, 1))
ggsave(file.path(out_dir, "Fig3B_pancancer_method_summary.png"), p_fig3b, width = 14, height = 8, dpi = 300, bg = "white")

# ==============================================================================
# PART 2: FIGURE 4D & 4E - Cellular Senescence ROCs & Pan-Cancer Hallmark Cor
# ==============================================================================
cat("=== Part 2: Figure 4D-E (Senescence Cohorts & Pan-Cancer Correlations) ===\n")

sen_colors <- c("PaaSc" = "#E64B35", "GSDensity" = "#F39B7F", "CelliD" = "#3C5488")
sen_methods <- c("PaaSc", "GSDensity", "CelliD")

generate_roc_plot <- function(file_name, dataset_name) {
  dat <- read.table(file.path(zen_dir, "Fig4", file_name), header = TRUE)
  roc_rows <- list()
  auc_rows <- list()
  for (m in sen_methods) {
    rc <- suppressMessages(pROC::roc(dat$Group, dat[[m]], quiet = TRUE))
    auc_val <- as.numeric(pROC::auc(rc))
    add_summary("Fig4D", dataset_name, paste0(m, "_AUC"), sprintf("%.4f", auc_val))
    roc_rows[[m]] <- tibble(
      Method = m,
      FPR = 1 - rc$specificities,
      TPR = rc$sensitivities
    ) %>% arrange(FPR, TPR)
    auc_rows[[m]] <- tibble(Method = m, AUC = auc_val)
  }
  roc_df <- bind_rows(roc_rows) %>% mutate(Method = factor(Method, levels = sen_methods))
  auc_df <- bind_rows(auc_rows)
  lbls <- setNames(paste0(auc_df$Method, " (AUC = ", sprintf("%.3f", auc_df$AUC), ")"), auc_df$Method)

  ggplot(roc_df, aes(x = FPR, y = TPR, color = Method)) +
    geom_line(linewidth = 1.15) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "grey60") +
    scale_color_manual(values = sen_colors, labels = lbls) +
    labs(title = dataset_name, x = "1 - Specificity", y = "Sensitivity", color = NULL) +
    theme_minimal(base_size = 11.5) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.9),
      panel.grid = element_blank(),
      legend.position = c(0.97, 0.05),
      legend.justification = c(1, 0),
      legend.background = element_rect(fill = alpha("white", 0.85), color = "grey80"),
      legend.text = element_text(size = 9.5)
    )
}

p4_1 <- generate_roc_plot("GSE102090_Score.tsv",     "GSE102090 (HUVEC)")
p4_2 <- generate_roc_plot("GSE119807_Score.tsv",     "GSE119807 (HCA2)")
p4_3 <- generate_roc_plot("GSE115301_10X_Score.tsv", "GSE115301-10X (IMR90)")

p_fig4d <- (p4_1 | p4_2 | p4_3) +
  plot_annotation(
    title = "Fig 4D: Cellular Senescence Detection ROC Benchmark (SenMayo Signature)",
    theme = theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5))
  )
ggsave(file.path(out_dir, "Fig4D_senescence_roc.png"), p_fig4d, width = 14, height = 4.8, dpi = 300, bg = "white")

cor_df <- read.table(file.path(zen_dir, "Fig4", "all_pid_celltype_hallmarkVSsene_correlation.tsv"),
                     sep = "\t", header = TRUE)
cor_df$padj <- p.adjust(cor_df$pvalue, method = "fdr")
fig4e_mat <- cor_df %>%
  group_by(hallmark, celltype) %>%
  summarise(
    cor_pos = n_distinct(pid[cor > 0.2 & padj < 0.05]),
    cor_neg = n_distinct(pid[cor < -0.2 & padj < 0.05]),
    .groups = "drop"
  ) %>%
  mutate(pos_pct = cor_pos / (cor_pos + cor_neg)) %>%
  select(hallmark, celltype, pos_pct) %>%
  pivot_wider(names_from = celltype, values_from = pos_pct) %>%
  column_to_rownames("hallmark")

fig4e_mat <- fig4e_mat[order(rowMeans(fig4e_mat), decreasing = TRUE), ]
rownames(fig4e_mat) <- sub("^HALLMARK_", "", rownames(fig4e_mat))
add_summary("Fig4E", "Top_Senescence_Hallmark", rownames(fig4e_mat)[1], sprintf("%.3f", mean(as.numeric(fig4e_mat[1, ]))))

ht4e <- Heatmap(
  as.matrix(fig4e_mat),
  name = "Pos Fraction",
  col = colorRamp2(c(0, 0.5, 1.0), c("#4393C3", "white", "#D6604D")),
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  row_names_gp = gpar(fontsize = 7.5),
  column_names_gp = gpar(fontsize = 10, fontface = "bold"),
  column_names_rot = 35,
  column_title = "Fig 4E: Pan-Cancer Patient-Level Senescence vs. Hallmark Positive Fraction",
  column_title_gp = gpar(fontsize = 11, fontface = "bold"),
  rect_gp = gpar(col = "white", lwd = 0.8)
)
png(file.path(out_dir, "Fig4E_pancancer_hallmark_senescence_cor.png"), width = 2600, height = 3300, res = 300)
draw(ht4e, padding = unit(c(4, 16, 4, 6), "mm"))
dev.off()

# ==============================================================================
# PART 3: FIGURE 5 - Tabula Muris Senis 10,725 Cells x 75 UKB/GWAS Polygenic Traits
# ==============================================================================
cat("=== Part 3: Figure 5 (TMS Multi-Trait Cellular Specificity) ===\n")

ct_anno   <- read.table(file.path(zen_dir, "SupFigs", "ct_anno.txt"), sep = "\t", header = TRUE)
tms_meta  <- read.table(file.path(zen_dir, "SupFigs", "TMS_facs_selected_metadata.tsv"), sep = "\t", header = TRUE, row.names = 1)
paasc5    <- read.table(file.path(zen_dir, "Fig5", "TMS_facs_selected_alltraits_PaaSc_score.tsv"), header = TRUE, sep = "\t", check.names = FALSE)
gsd5      <- read.table(file.path(zen_dir, "Fig5", "TMS_facs_selected_alltraits_GSDensity_score.tsv"), header = TRUE, sep = "\t", check.names = FALSE)
cellid5   <- read.table(file.path(zen_dir, "Fig5", "TMS_facs_selected_alltraits_CelliD_score.tsv"), header = TRUE, sep = "\t", check.names = FALSE)

# Ensure exact row alignment
common_cells <- rownames(tms_meta)
paasc5  <- paasc5[common_cells, ]
gsd5    <- gsd5[common_cells, ]
cellid5 <- cellid5[common_cells, ]

tms_meta$celltype_abbr <- ct_anno$abbreviation[match(tms_meta$cell_ontology_class, ct_anno$CellType)]
tms_meta$tissue        <- ct_anno$tissue[match(tms_meta$cell_ontology_class, ct_anno$CellType)]

# Fig5A: Boxplot of Lymphocyte Count
t_trait <- "UKB_460K.blood_LYMPHOCYTE_COUNT"
fig5a_df <- tibble(
  Score    = paasc5[[t_trait]],
  celltype = tms_meta$celltype_abbr,
  tissue   = tms_meta$tissue
) %>%
  mutate(celltype = reorder(celltype, Score, FUN = median, decreasing = TRUE))

tissue_cols <- c("blood/immune" = "#C17188", "brain" = "#9AB846", "others" = "#F5A433")

p5a <- ggplot(fig5a_df, aes(x = celltype, y = Score, color = tissue)) +
  geom_boxplot(outlier.size = 0.4) +
  scale_color_manual(values = tissue_cols) +
  coord_cartesian(ylim = c(-4, 3.2)) +
  labs(title = "Fig 5A: Lymphocyte Count GWAS Score Across 19 TMS Cell Types", x = NULL, y = "PaaSc score", color = "Tissue") +
  theme_minimal(base_size = 11) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.text.x = element_text(angle = 55, hjust = 1, vjust = 1, size = 8.5, color = "black"),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.9),
    panel.grid = element_blank()
  )

# Fig5B: GMM Binarization Histogram
set.seed(1)
label_gmm <- doBinarization(paasc5[, t_trait, drop = FALSE], n.cluster = 2, method = "GMM")
pos_flag  <- label_gmm[, 1] == "positive"
gmm_thr   <- min(paasc5[pos_flag, t_trait])
add_summary("Fig5B", t_trait, "GMM_threshold", sprintf("%.4f", gmm_thr))
add_summary("Fig5B", t_trait, "N_positive_cells", paste0(sum(pos_flag), "/", length(pos_flag)))

p5b <- ggplot(paasc5, aes(x = .data[[t_trait]])) +
  geom_histogram(aes(y = after_stat(density)), binwidth = 0.2, fill = "slategray2", color = "black", alpha = 0.75) +
  geom_density(color = "dodgerblue4", linewidth = 1) +
  geom_vline(xintercept = gmm_thr, color = "dodgerblue4", linetype = "dashed", linewidth = 1.2) +
  annotate("text", x = gmm_thr + 0.15, y = 0.68, label = paste0("GMM threshold = ", sprintf("%.3f", gmm_thr)),
           hjust = 0, fontface = "bold", size = 3.8) +
  scale_x_continuous(limits = c(-3, 3), expand = c(0, 0)) +
  labs(title = "Fig 5B: Two-Component GMM Binarization", x = "PaaSc Score (Lymphocyte Count)", y = "Density") +
  theme_minimal(base_size = 11) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.line = element_line(linewidth = 0.8),
    panel.grid = element_blank()
  )

# Fig5C: Fisher's Exact Test Cell-Type Enrichment (with Haldane-Anscombe +0.5 continuity for finite ORs)
or_list <- lapply(unique(tms_meta$celltype_abbr), function(c_abbr) {
  in_ct <- tms_meta$celltype_abbr == c_abbr
  tab <- matrix(c(
    sum(pos_flag & in_ct), sum(!pos_flag & in_ct),
    sum(pos_flag & !in_ct), sum(!pos_flag & !in_ct)
  ), nrow = 2)
  ft_raw <- fisher.test(tab)
  # Haldane-Anscombe +0.5 only for finite log-OR plotting when zero cells occur
  or_smooth <- ((tab[1, 1] + 0.5) * (tab[2, 2] + 0.5)) / ((tab[1, 2] + 0.5) * (tab[2, 1] + 0.5))
  tibble(
    CellType = c_abbr,
    Odds_ratio_raw = unname(ft_raw$estimate),
    Odds_ratio_smooth = or_smooth,
    Pvalue = ft_raw$p.value
  )
})
fig5c_df <- bind_rows(or_list) %>%
  mutate(
    padj = pmax(p.adjust(Pvalue, method = "fdr"), 1e-10),
    logFDR = -log10(padj),
    log2_OR = log2(Odds_ratio_smooth),
    tissue = ct_anno$tissue[match(CellType, ct_anno$abbreviation)],
    CellType = factor(CellType, levels = CellType[order(log2_OR)])
  )

p5c <- ggplot(fig5c_df, aes(x = log2_OR, y = CellType, color = logFDR, shape = tissue)) +
  geom_vline(xintercept = 0, color = "grey50", linetype = "dashed") +
  geom_point(size = 3.8) +
  scale_color_distiller(palette = "Spectral", name = "-log10(FDR)") +
  labs(title = "Fig 5C: Fisher Enrichment (GMM Positives)", x = "log2(Odds Ratio) [Haldane-corrected]", y = NULL, shape = "Tissue") +
  theme_minimal(base_size = 11) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.text.y = element_text(size = 8.5, color = "black"),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.9),
    panel.grid = element_blank()
  )

p_fig5abc <- (p5a / (p5b | p5c)) + plot_layout(heights = c(1, 1.15)) +
  plot_annotation(
    title = "Fig 5A-C: Tabula Muris Senis FACS Polygenic Risk Localization (n = 10,725 cells)",
    theme = theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5))
  )
ggsave(file.path(out_dir, "Fig5ABC_tms_lymphocyte_enrichment.png"), p_fig5abc, width = 14, height = 10, dpi = 300, bg = "white")

# Fig5D: Multi-Trait (9 blood + 6 brain + 7 others = 22 traits) x 19 Cell Types Across PaaSc / GSDensity / CelliD
traits_blood <- c(
  "UKB_460K.blood_LYMPHOCYTE_COUNT", "PASS_Multiple_sclerosis", "UKB_460K.disease_AID_ALL",
  "UKB_460K.disease_HYPOTHYROIDISM_SELF_REP", "PASS_Celiac", "PASS_Rheumatoid_Arthritis",
  "UKB_460K.blood_EOSINOPHIL_COUNT", "UKB_460K.disease_ASTHMA_DIAGNOSED", "UKB_460K.blood_MONOCYTE_COUNT"
)
traits_brain <- c(
  "PASS_Worry_Nagel2018", "PASS_ReactionTime_Davies2018", "PASS_MDD_Howard2019",
  "UKB_460K.cov_EDU_COLLEGE", "PASS_Schizophrenia_Pardinas2018", "PASS_Intelligence_SavageJansen2018"
)
traits_other <- c(
  "PASS_Coronary_Artery_Disease", "UKB_460K.bp_DIASTOLICadjMEDz", "UKB_460K.biochemistry_AlkalinePhosphatase",
  "UKB_460K.cancer_BREAST", "UKB_460K.biochemistry_LDLdirect", "UKB_460K.lung_FEV1FVCzSMOKE",
  "UKB_460K.biochemistry_Glucose"
)

trait_meta <- tibble(
  Trait = c(traits_blood, traits_brain, traits_other),
  Tissue = c(rep("blood/immune", 9), rep("brain", 6), rep("others", 7)),
  Label = c(
    "Lymphocyte count", "Multiple sclerosis", "Autoimmune disease", "Hypothyroidism",
    "Celiac disease", "Rheumatoid arthritis", "Eosinophil count", "Asthma", "Monocyte count",
    "Worry/Neuroticism", "Reaction time", "Major depression (MDD)", "College education",
    "Schizophrenia", "Intelligence",
    "Coronary artery disease", "Diastolic BP", "Alkaline phosphatase", "Breast cancer",
    "LDL cholesterol", "Lung FEV1/FVC", "Fasting/Blood glucose"
  )
)

ct_ordered <- ct_anno %>%
  mutate(tissue = factor(tissue, levels = c("blood/immune", "brain", "others"))) %>%
  arrange(tissue)

compute_trait_ct_stats <- function(score_mat, method_name) {
  rows <- list()
  for (tr in trait_meta$Trait) {
    vals <- score_mat[[tr]]
    for (c_abbr in ct_ordered$abbreviation) {
      idx <- tms_meta$celltype_abbr == c_abbr
      wt  <- wilcox.test(vals[idx], vals[!idx], alternative = "greater")
      rows[[length(rows) + 1]] <- tibble(
        Trait = tr, CellType = c_abbr, Mean = mean(vals[idx]), Pval = wt$p.value
      )
    }
  }
  out <- bind_rows(rows) %>%
    mutate(
      FDR = p.adjust(Pval, method = "fdr"),
      Method = method_name,
      Label = factor(trait_meta$Label[match(Trait, trait_meta$Trait)], levels = rev(trait_meta$Label)),
      TraitTissue = factor(trait_meta$Tissue[match(Trait, trait_meta$Trait)], levels = c("others", "brain", "blood/immune")),
      CellType = factor(CellType, levels = ct_ordered$abbreviation),
      CellTissue = factor(ct_ordered$tissue[match(CellType, ct_ordered$abbreviation)], levels = c("blood/immune", "brain", "others"))
    )
  out
}

df5d_paasc  <- compute_trait_ct_stats(paasc5,  "PaaSc")
df5d_gsd    <- compute_trait_ct_stats(gsd5,    "GSDensity")
df5d_cellid <- compute_trait_ct_stats(cellid5, "CelliD")

# Match expected diagonal enrichment counts
diag_paasc <- df5d_paasc %>% filter(as.character(TraitTissue) == as.character(CellTissue) & FDR < 0.05) %>% nrow()
diag_gsd   <- df5d_gsd   %>% filter(as.character(TraitTissue) == as.character(CellTissue) & FDR < 0.05) %>% nrow()
diag_cid   <- df5d_cellid %>% filter(as.character(TraitTissue) == as.character(CellTissue) & FDR < 0.05) %>% nrow()
add_summary("Fig5D", "Matched_Tissue_SigPairs", "PaaSc", diag_paasc)
add_summary("Fig5D", "Matched_Tissue_SigPairs", "GSDensity", diag_gsd)
add_summary("Fig5D", "Matched_Tissue_SigPairs", "CelliD", diag_cid)

plot_5d_panel <- function(df, title_str, is_paasc = FALSE, show_y = TRUE) {
  if (is_paasc) {
    df$Mean_plot <- pmin(pmax(df$Mean, -2.5), 2.5)
    fill_scale <- scale_fill_distiller(palette = "RdBu", direction = -1, limits = c(-2.5, 2.5), name = "Mean")
  } else {
    q99 <- quantile(df$Mean, 0.99)
    df$Mean_plot <- pmin(df$Mean, q99)
    fill_scale <- scale_fill_gradient(low = "#FEFAFA", high = "#D73027", name = "Mean")
  }
  p <- ggplot(df, aes(x = interaction(CellType, CellTissue), y = interaction(Label, TraitTissue))) +
    geom_point(aes(fill = Mean_plot), color = "grey90", shape = 22, size = 4.3) +
    geom_text(aes(label = ifelse(FDR < 0.05, "*", "")), size = 3.5, color = "black", vjust = 0.72) +
    fill_scale +
    guides(x = "axis_nested", y = "axis_nested") +
    labs(title = title_str, x = NULL, y = NULL) +
    theme_minimal(base_size = 10) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 12),
      panel.grid = element_blank(),
      axis.text.x = element_text(angle = 50, hjust = 1, size = 7.5, color = "black"),
      axis.text.y = if (show_y) element_text(size = 8, color = "black") else element_blank()
    )
  p
}

p5d_1 <- plot_5d_panel(df5d_paasc,  paste0("PaaSc (Matched sig = ", diag_paasc, ")"),     is_paasc = TRUE,  show_y = TRUE)
p5d_2 <- plot_5d_panel(df5d_gsd,    paste0("GSDensity (Matched sig = ", diag_gsd, ")"),   is_paasc = FALSE, show_y = FALSE)
p5d_3 <- plot_5d_panel(df5d_cellid, paste0("CelliD (Matched sig = ", diag_cid, ")"),      is_paasc = FALSE, show_y = FALSE)

p_fig5d <- (p5d_1 | p5d_2 | p5d_3) +
  plot_annotation(
    title = "Fig 5D: Tissue-Matched Polygenic Trait Enrichment Across 19 TMS Cell Types (22 Representative Traits)",
    theme = theme(plot.title = element_text(size = 13.5, face = "bold", hjust = 0.5))
  )
ggsave(file.path(out_dir, "Fig5D_tms_22traits_comparison.png"), p_fig5d, width = 18, height = 7.2, dpi = 300, bg = "white")

# ==============================================================================
# PART 4: FIGURE 7B - scATAC-seq PBMC 16-Celltype Gene-Activity ROC Benchmark
# ==============================================================================
cat("=== Part 4: Figure 7B (scATAC-seq PBMC Benchmark) ===\n")

e7 <- new.env()
load(file.path(zen_dir, "Fig7", "granulocyte_ROC.RData"), envir = e7)
out7 <- e7$out

add_summary("Fig7B", "scATAC_16subtypes", "PaaSc_median_AUC",     sprintf("%.4f", median(out7$PaaSc)))
add_summary("Fig7B", "scATAC_16subtypes", "CelliD_median_AUC",    sprintf("%.4f", median(out7$CelliD)))
add_summary("Fig7B", "scATAC_16subtypes", "GSdensity_median_AUC", sprintf("%.4f", median(out7$GSdensity)))
add_summary("Fig7B", "scATAC_16subtypes", "pDC_PaaSc_vs_CelliD",  paste0(sprintf("%.4f", out7$PaaSc[out7$cellType == "pDC"]), " vs ", sprintf("%.4f", out7$CelliD[out7$cellType == "pDC"])))

df7_long <- out7 %>%
  pivot_longer(cols = c(PaaSc, CelliD, GSdensity), names_to = "Method", values_to = "ROC") %>%
  mutate(Method = factor(Method, levels = c("PaaSc", "CelliD", "GSdensity")))

pal16 <- c("#E74C3C", "#F39C12", "#F1C40F", "#2ECC71", "#1ABC9C", "#3498DB", "#9B59B6", "#E91E63",
           "#C0392B", "#2980B9", "#27AE60", "#8E44AD", "#D68910", "#7F8C8D", "#34495E", "#16A085")

p_fig7b <- ggplot(df7_long, aes(x = Method, y = ROC)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA, width = 0.55, fill = "grey95") +
  geom_point(aes(color = cellType), size = 3.2, alpha = 0.9, position = position_jitter(width = 0.16, seed = 123)) +
  scale_color_manual(values = pal16, name = "Cell Type") +
  ylim(0.35, 1.03) +
  labs(
    title = "Fig 7B: PBMC scATAC-seq (10,459 cells) Gene-Activity Lineage ROC Benchmark",
    subtitle = sprintf("Median AUCs -> PaaSc: %.4f | CelliD: %.4f | GSdensity: %.4f",
                       median(out7$PaaSc), median(out7$CelliD), median(out7$GSdensity)),
    x = NULL, y = "ROC AUC"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 12),
    plot.subtitle = element_text(hjust = 0.5, color = "grey30", size = 10),
    axis.text.x = element_text(size = 11.5, face = "bold", color = "black"),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.9),
    panel.grid.major.x = element_blank()
  ) +
  guides(color = guide_legend(ncol = 2))

ggsave(file.path(out_dir, "Fig7B_scatac_granulocyte_roc.png"), p_fig7b, width = 9.8, height = 5.8, dpi = 300, bg = "white")

# ==============================================================================
# PART 5: FIGURE 9 - Spatial CosMx SMI NSCLC (98,002 cells) Benchmark & Correlation
# ==============================================================================
cat("=== Part 5: Figure 9 (CosMx SMI NSCLC Spatial 98,002 cells) ===\n")

auc9 <- read.table(file.path(zen_dir, "Fig9", "spatial_lung_AUC.txt"), header = TRUE) %>%
  rownames_to_column("Celltype") %>%
  pivot_longer(-Celltype, names_to = "method", values_to = "AUC") %>%
  group_by(method) %>%
  mutate(median_auc = median(AUC, na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(method = factor(method, levels = unique(method[order(-median_auc)])))

med9 <- auc9 %>% distinct(method, median_auc)
for (r in seq_len(nrow(med9))) {
  add_summary("Fig9B", "CosMx_NSCLC_9celltypes", paste0(med9$method[r], "_median_AUC"), sprintf("%.4f", med9$median_auc[r]))
}

nsclc_pal <- c("B cell"="#E74C3C","NK"="#F39C12","Endothelial cell"="#F1C40F","Epithelial cell"="#2ECC71",
               "Fibroblast"="#1ABC9C","Macrophage"="#3498DB","Mast cell"="#9B59B6","Plasma cell"="#8E44AD",
               "T/NK"="#D68910","Tumor"="#34495E")

p9b <- ggplot(auc9, aes(x = method, y = AUC)) +
  geom_boxplot(outlier.shape = NA, color = "black", width = 0.55, fill = "grey95") +
  geom_jitter(aes(color = Celltype), size = 3.2, width = 0.15) +
  scale_color_manual(values = nsclc_pal) +
  labs(title = "Fig 9B: NSCLC Cell-Type AUC", x = NULL, y = "ROC AUC", color = "Cell Type") +
  theme_minimal(base_size = 11.5) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.text.x = element_text(size = 11, face = "bold", color = "black"),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.9),
    panel.grid.major.x = element_blank()
  )

smi9 <- read.csv(file.path(zen_dir, "Fig9", "ST_lung_hallmark_SMI.csv"), row.names = 1) %>%
  rownames_to_column("Hallmark") %>%
  mutate(Hallmark_short = sub("^HALLMARK_", "", Hallmark)) %>%
  arrange(desc(MI))

add_summary("Fig9D", "Top_Spatial_Hallmark", smi9$Hallmark_short[1], sprintf("%.4f", smi9$MI[1]))

p9d <- ggplot(head(smi9, 20), aes(x = MI, y = reorder(Hallmark_short, MI))) +
  geom_col(fill = "#D95F02", alpha = 0.85, width = 0.7) +
  labs(title = "Fig 9D: Top 20 Spatial Auto-Correlated Hallmarks", x = "Spatial Moran/MI (PaaSc)", y = NULL) +
  theme_minimal(base_size = 10.5) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.text.y = element_text(size = 8, color = "black"),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.9)
  )

# Vectorized Pearson correlation & exact p-value on 98,002 cells (matching Fig9.R)
f9_cell <- read.table(file.path(zen_dir, "Fig9", "spatial_lung_PaaSc.txt"), header = TRUE, check.names = FALSE)
f9_hall <- data.table::fread(file.path(zen_dir, "Fig9", "ST_lung_hallmark_paasc.csv")) %>%
  column_to_rownames("V1")

# Align rows identically
f9_hall <- f9_hall[rownames(f9_cell), ]
n_cells9 <- nrow(f9_cell)
cor_mat <- cor(as.matrix(f9_cell), as.matrix(f9_hall))
t_mat   <- cor_mat * sqrt((n_cells9 - 2) / pmax(1 - cor_mat^2, 1e-12))
p_mat   <- 2 * pt(-abs(t_mat), df = n_cells9 - 2)
cor_sig <- ifelse(!is.na(p_mat) & p_mat <= 0.05, cor_mat, 0)
colnames(cor_sig) <- sub("^HALLMARK_", "", colnames(cor_sig))

rc_max <- arrayInd(which.max(cor_sig), dim(cor_sig))
add_summary("Fig9E", "CosMx_98002_cells", "Max_Pearson_r_pair",
            paste0(rownames(cor_sig)[rc_max[1]], " vs ",
                   colnames(cor_sig)[rc_max[2]],
                   " (r=", sprintf("%.3f", cor_sig[rc_max]), ")"))

col_fun9 <- colorRamp2(c(-0.6, -0.3, 0, 0.3, 0.6, 1.0),
                       c("#2166AC", "#67A9CF", "white", "#FDDBC7", "#EF8A62", "#B2182B"))

ht9e <- Heatmap(
  t(cor_sig),
  name = "Pearson r",
  col = col_fun9,
  cluster_rows = TRUE,
  cluster_columns = TRUE,
  row_names_gp = gpar(fontsize = 6.8),
  column_names_gp = gpar(fontsize = 8.5, fontface = "bold"),
  column_names_rot = 45,
  column_title = "Fig 9E: Cell-Type vs. Hallmark Spatial Correlation (n = 98,002 cells, p<=0.05)",
  column_title_gp = gpar(fontsize = 11.5, fontface = "bold"),
  heatmap_legend_param = list(at = c(-0.6, -0.3, 0, 0.3, 0.6, 1.0))
)
grob9e <- wrap_elements(full = grid.grabExpr(draw(ht9e)))

p_fig9_top <- (p9b | p9d) + plot_layout(widths = c(1, 1.3))
p_fig9_all <- (p_fig9_top / grob9e) + plot_layout(heights = c(1, 1.65)) +
  plot_annotation(
    title = "Fig 9: Spatial NanoString CosMx SMI NSCLC Tumor Microenvironment Benchmark (98,002 cells)",
    theme = theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5))
  )
ggsave(file.path(out_dir, "Fig9_spatial_nsclc_benchmark_and_cor.png"), p_fig9_all, width = 14, height = 13.5, dpi = 300, bg = "white")

# ==============================================================================
# PART 6: Export Aggregate Quantitative Summary Table
# ==============================================================================
summary_tbl <- bind_rows(summary_rows)
write_tsv(summary_tbl, file.path(out_dir, "fig3_to_fig9_summary.tsv"))

cat("=== Multi-Figure Suite Finished Successfully! ===\n")
print(summary_tbl, n = 60)
