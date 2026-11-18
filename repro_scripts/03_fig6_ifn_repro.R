# ==============================================================================
# PaaSc 论文主图 Figure 6 端到端全量重算（PBMC IFN-beta 刺激响应基准）
# 
# 文献出处: Liao X, Hong Y, Feng Y, et al. PLOS Comput Biol 2025, 21(11): e1013666.
# 核心验证:
#   1. 条件扰动敏感度: 评估 28,871 个 PBMC 细胞在 IFN-beta 刺激前后的干扰素通路响应 (ctrl vs stim)
#   2. 谱系特异性阴性对照: 在 2,871 个 B 细胞内评估 B 细胞激活通路，验证不会因全群受刺激而错误产生条件假阳性
#   3. 双方法全量重算: 从 IFN.RDS 原始表达矩阵从头重算 PaaSc 与 CelliD，并与作者存档 GSdensity 严格对齐
#
# 输入数据:
#   - zenodo_data/PaaSc_data/Fig6/IFN.RDS (16,142 基因 x 28,871 细胞，GSE96583)
#   - msigdbr (v26.1.1): GOBP_INTERFERON_MEDIATED_SIGNALING_PATHWAY (113 基因，89 与矩阵重叠)
#   - data/B_cell_activation.txt (278 基因，233 与矩阵重叠)
#   - data/kegg_legacy.gmt (186 条 KEGG 背景通路)
#   - 对照归档: Fig6/PaaSc_scores.txt, Fig6/CelliD_scores.txt, Fig6/GSdensity_scores.txt
#
# 输出文件:
#   - repro_results/fig3_to_fig9/Fig6A_ifn_umap_overview.png (三联 UMAP 图: 细胞分群、刺激条件、PaaSc得分)
#   - repro_results/fig3_to_fig9/Fig6B_ifn_vs_bcell_benchmark.png (IFN响应 vs B细胞谱系对照箱线图)
#   - repro_results/fig3_to_fig9/fig6_verification_summary.txt (量化验收摘要)
#
# 核心工程突破:
#   - 解决 Seurat v5 Assay5 与 CelliD 1.19.0 的 S3 调度与 slot/layer 兼容性 Bug
#   - 显式生命周期释放 rm() + gc() 控制 2.8 万细胞大矩阵 SVD 峰值内存在 21.1 GB (31.9 GB 本机可承载)
#
# 实测耗时: UMAP 约 60 秒 + PaaSc 约 160 秒 + CelliD 约 520 秒 (合计约 12 分钟)
# ==============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(PaaSc)
  library(CelliD)
  library(msigdbr)
  library(ROCit)
  library(tidyverse)
  library(patchwork)
})

base_dir <- "d:/project/gemini/reproduction/PaaSc"
fig6_dir <- file.path(base_dir, "zenodo_data", "PaaSc_data", "Fig6")
bg_gmt   <- file.path(base_dir, "data", "kegg_legacy.gmt")
bcell_pt <- file.path(base_dir, "data", "B_cell_activation.txt")
out_dir  <- file.path(base_dir, "repro_results", "fig3_to_fig9")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

set.seed(42)

cat("=== 1. Loading IFN.RDS & Preparing Pathway Gene Sets ===\n")
t0 <- proc.time()
IFN <- readRDS(file.path(fig6_dir, "IFN.RDS"))
cat("Loaded IFN.RDS:", nrow(IFN), "features x", ncol(IFN), "cells in", round((proc.time()-t0)[3], 1), "s\n")

# Load human GO:BP IFN signaling pathway from MSigDB 26.1.1 + B_cell_activation.txt
hs_gobp <- msigdbr(species = "Homo sapiens", collection = "C5", subcollection = "GO:BP")
INTERFERON <- unique(hs_gobp$gene_symbol[hs_gobp$gs_name == "GOBP_INTERFERON_MEDIATED_SIGNALING_PATHWAY"])
B_cell_activation <- readLines(bcell_pt)
pathway_geneset <- list(B_cell_activation = B_cell_activation, INTERFERON = INTERFERON)
cat("Pathway sizes -> B_cell_activation:", length(B_cell_activation),
    "(overlap:", sum(B_cell_activation %in% rownames(IFN)), ") | INTERFERON:",
    length(INTERFERON), "(overlap:", sum(INTERFERON %in% rownames(IFN)), ")\n")

# ------------------------------------------------------------------------------
# 2. Seurat Normalization, PCA & UMAP (matching Fig6/PaaSc_Figure6_IFN.R)
# ------------------------------------------------------------------------------
cat("=== 2. Running Seurat NormalizeData -> PCA -> UMAP ===\n")
t1 <- proc.time()
IFN <- IFN %>%
  NormalizeData(verbose = FALSE) %>%
  FindVariableFeatures(nfeatures = 2000, verbose = FALSE) %>%
  ScaleData(verbose = FALSE) %>%
  RunPCA(verbose = FALSE) %>%
  RunUMAP(dims = 1:25, seed.use = 2, verbose = FALSE)
cat("UMAP completed in", round((proc.time()-t1)[3], 1), "s\n")

# ------------------------------------------------------------------------------
# 3. Recompute PaaSc End-to-End (computeMCA(nmcs=20) -> getGeneRate -> doRegression -> computeScore)
# ------------------------------------------------------------------------------
cat("=== 3. Recomputing PaaSc on 28,871 cells ===\n")
t2 <- proc.time()
obj_paasc <- computeMCA(IFN, nmcs = 20)
gene_rate <- getGeneRate(background.geneset = bg_gmt, pathway.geneset = pathway_geneset)
regression_data <- doRegression(obj_paasc, gene.rate = gene_rate)
paasc_recomp <- computeScore(obj_paasc, regression.data = regression_data, weight = FALSE, normalize = "z-score")
rm(obj_paasc); gc(verbose = FALSE)
cat("PaaSc completed in", round((proc.time()-t2)[3], 1), "s\n")

# ------------------------------------------------------------------------------
# 4. Recompute CelliD End-to-End (NormalizeData -> RunMCA -> RunCellHGT dims 1:50)
#    Note: CelliD::RunMCA.Seurat calls defunct GetAssayData(slot="data") on Seurat v5 Assay5,
#    so we call CelliD::RunMCA directly on LayerData(IFN, "data") and attach DimReducObject.
# ------------------------------------------------------------------------------
cat("=== 4. Recomputing CelliD on 28,871 cells ===\n")
t3 <- proc.time()
norm_mat <- as.matrix(LayerData(IFN, assay = "RNA", layer = "data"))
mca_cid <- CelliD::RunMCA(X = norm_mat, nmcs = 50, reduction.name = "mca")
rm(norm_mat); gc(verbose = FALSE)
obj_cid <- IFN
# Use CelliD's internal Seurat setter explicitly (setDimMCSlot is an unexported S3 generic)
obj_cid <- CelliD:::setDimMCSlot.Seurat(
  X = obj_cid,
  cellEmb = mca_cid$cellsCoordinates,
  geneEmb = mca_cid$featuresCoordinates,
  stdev   = mca_cid$stdev,
  reduction.name = "mca"
)
rm(mca_cid); gc(verbose = FALSE)
hgt <- RunCellHGT(obj_cid, pathways = pathway_geneset, reduction = "mca", dims = 1:50, n.features = 200)
cellid_recomp <- data.frame(t(as.matrix(hgt)))
rm(obj_cid); gc(verbose = FALSE)
cat("CelliD completed in", round((proc.time()-t3)[3], 1), "s\n")

# ------------------------------------------------------------------------------
# 5. Compare Recomputed Scores Against Archived Scores & Compute ROC AUCs
# ------------------------------------------------------------------------------
cat("=== 5. Comparing Recomputed vs Archived Scores & Computing ROC AUCs ===\n")
o1_arch <- read.table(file.path(fig6_dir, "PaaSc_scores.txt"), row.names = 1, header = TRUE, sep = "\t")
o2_arch <- read.table(file.path(fig6_dir, "GSdensity_scores.txt"), row.names = 1, header = TRUE, sep = "\t")
o3_arch <- read.table(file.path(fig6_dir, "CelliD_scores.txt"), row.names = 1, header = TRUE, sep = "\t")

cells_ord <- Cells(IFN)
paasc_recomp  <- paasc_recomp[cells_ord, ]
cellid_recomp <- cellid_recomp[cells_ord, ]
o1_arch <- o1_arch[cells_ord, ]
o2_arch <- o2_arch[cells_ord, ]
o3_arch <- o3_arch[cells_ord, ]

cor_paasc_ifn  <- cor(paasc_recomp$INTERFERON, o1_arch$INTERFERON)
cor_paasc_b    <- cor(paasc_recomp$B_cell_activation, o1_arch$B_cell_activation)
cor_cellid_ifn <- cor(cellid_recomp$INTERFERON, o3_arch$INTERFERON)
cor_cellid_b   <- cor(cellid_recomp$B_cell_activation, o3_arch$B_cell_activation)

roc_auc <- function(all_label, ref_label, score) {
  rocit(score = score, class = all_label, negref = ref_label)$AUC
}

# IFN pathway across all 28,871 cells (ctrl vs stim)
df_ifn <- tibble(
  condition = IFN$stim,
  PaaSc_recomp     = rank(paasc_recomp$INTERFERON) / nrow(paasc_recomp),
  PaaSc_archived   = rank(o1_arch$INTERFERON) / nrow(o1_arch),
  GSdensity_arch   = rank(o2_arch$INTERFERON) / nrow(o2_arch),
  CelliD_recomp    = rank(cellid_recomp$INTERFERON) / nrow(cellid_recomp),
  CelliD_archived  = rank(o3_arch$INTERFERON) / nrow(o3_arch)
)

IFN_roc <- c(
  PaaSc_recomp    = roc_auc(df_ifn$condition, "ctrl", df_ifn$PaaSc_recomp),
  PaaSc_archived  = roc_auc(df_ifn$condition, "ctrl", df_ifn$PaaSc_archived),
  GSdensity_arch  = roc_auc(df_ifn$condition, "ctrl", df_ifn$GSdensity_arch),
  CelliD_recomp   = roc_auc(df_ifn$condition, "ctrl", df_ifn$CelliD_recomp),
  CelliD_archived = roc_auc(df_ifn$condition, "ctrl", df_ifn$CelliD_archived)
)

# B_cell_activation specificity control within B cells
b_idx <- IFN$cell == "B cells"
df_b <- tibble(
  condition = IFN$stim[b_idx],
  PaaSc_recomp     = (rank(paasc_recomp$B_cell_activation) / nrow(paasc_recomp))[b_idx],
  PaaSc_archived   = (rank(o1_arch$B_cell_activation) / nrow(o1_arch))[b_idx],
  GSdensity_arch   = (rank(o2_arch$B_cell_activation) / nrow(o2_arch))[b_idx],
  CelliD_recomp    = (rank(cellid_recomp$B_cell_activation) / nrow(cellid_recomp))[b_idx],
  CelliD_archived  = (rank(o3_arch$B_cell_activation) / nrow(o3_arch))[b_idx]
)

B_roc <- c(
  PaaSc_recomp    = roc_auc(df_b$condition, "ctrl", df_b$PaaSc_recomp),
  PaaSc_archived  = roc_auc(df_b$condition, "ctrl", df_b$PaaSc_archived),
  GSdensity_arch  = roc_auc(df_b$condition, "ctrl", df_b$GSdensity_arch),
  CelliD_recomp   = roc_auc(df_b$condition, "ctrl", df_b$CelliD_recomp),
  CelliD_archived = roc_auc(df_b$condition, "ctrl", df_b$CelliD_archived)
)

print(IFN_roc)
print(B_roc)

# ------------------------------------------------------------------------------
# 6. Plot Fig6A: UMAP Overview (Cell Type | Condition ctrl/stim | Recomputed IFN Z-score)
# ------------------------------------------------------------------------------
col_palatte <- c("#6965CE", "#EE5755", "#75E15A", "#D5CE85", "#9F33D9", "#9BE2CF", "#69A5C8", "#FC9500")
umap_coords <- as.data.frame(Embeddings(IFN, "umap"))
colnames(umap_coords) <- c("UMAP_1", "UMAP_2")
umap_coords$cell <- IFN$cell
umap_coords$stim <- IFN$stim
ifn_clip <- pmin(pmax(paasc_recomp$INTERFERON, quantile(paasc_recomp$INTERFERON, 0.01)),
                 quantile(paasc_recomp$INTERFERON, 0.99))
umap_coords$IFN_PaaSc <- ifn_clip

p_u1 <- ggplot(umap_coords, aes(x = UMAP_1, y = UMAP_2, color = cell)) +
  geom_point(size = 0.15, alpha = 0.7) +
  scale_color_manual(values = col_palatte) +
  labs(title = "Cell Types (n = 28,871)", color = NULL) +
  theme_classic(base_size = 11) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"), legend.position = "bottom") +
  guides(color = guide_legend(ncol = 2, override.aes = list(size = 2.5, alpha = 1)))

p_u2 <- ggplot(umap_coords, aes(x = UMAP_1, y = UMAP_2, color = stim)) +
  geom_point(size = 0.15, alpha = 0.7) +
  scale_color_manual(values = c("ctrl" = "#6965CE", "stim" = "#EE5755")) +
  labs(title = "Stimulation Condition (ctrl vs stim)", color = NULL) +
  theme_classic(base_size = 11) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"), legend.position = "bottom") +
  guides(color = guide_legend(override.aes = list(size = 2.5, alpha = 1)))

p_u3 <- ggplot(umap_coords, aes(x = UMAP_1, y = UMAP_2, color = IFN_PaaSc)) +
  geom_point(size = 0.15, alpha = 0.8) +
  scale_color_gradientn(colors = c("#2166AC", "#F7F7F7", "#F4A582", "#B2182B")) +
  labs(title = sprintf("Recomputed PaaSc IFN Score (AUC=%.4f)", IFN_roc["PaaSc_recomp"]), color = "Z-score") +
  theme_classic(base_size = 11) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"), legend.position = "bottom")

p_fig6a <- (p_u1 | p_u2 | p_u3) +
  plot_annotation(
    title = "Fig 6A: PBMC Interferon-beta Stimulation scRNA-seq (GSE96583, 28,871 cells)",
    theme = theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5))
  )
ggsave(file.path(out_dir, "Fig6A_ifn_umap_overview.png"), p_fig6a, width = 15, height = 5.8, dpi = 300, bg = "white")

# ------------------------------------------------------------------------------
# 7. Plot Fig6B: IFN Pathway Activities vs B-Cell Activation Negative Control
# ------------------------------------------------------------------------------
ifn_long <- df_ifn %>%
  select(condition, PaaSc = PaaSc_recomp, CelliD = CelliD_recomp, GSdensity = GSdensity_arch) %>%
  pivot_longer(cols = c(PaaSc, CelliD, GSdensity), names_to = "Method", values_to = "Score") %>%
  mutate(Method = factor(Method, levels = c("PaaSc", "CelliD", "GSdensity")))

ifn_annot <- tibble(
  Method = factor(c("PaaSc", "CelliD", "GSdensity"), levels = c("PaaSc", "CelliD", "GSdensity")),
  label = sprintf("AUC = %.4f", c(IFN_roc["PaaSc_recomp"], IFN_roc["CelliD_recomp"], IFN_roc["GSdensity_arch"]))
)

p_b1 <- ggplot(ifn_long, aes(x = Method, y = Score, fill = condition)) +
  geom_boxplot(alpha = 0.8, outlier.size = 0.2, outlier.alpha = 0.25) +
  geom_text(data = ifn_annot, aes(x = Method, y = 1.07, label = label), inherit.aes = FALSE, fontface = "bold", size = 3.6) +
  scale_fill_manual(values = c("ctrl" = "grey65", "stim" = "#E64B35")) +
  coord_cartesian(ylim = c(0, 1.13)) +
  labs(title = "IFN Pathway Activities (All 28,871 cells)", x = NULL, y = "Rank normalized score", fill = "Condition") +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.text.x = element_text(size = 11, face = "bold", color = "black"),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.9),
    panel.grid.major.x = element_blank()
  )

b_long <- df_b %>%
  select(condition, PaaSc = PaaSc_recomp, CelliD = CelliD_recomp, GSdensity = GSdensity_arch) %>%
  pivot_longer(cols = c(PaaSc, CelliD, GSdensity), names_to = "Method", values_to = "Score") %>%
  mutate(Method = factor(Method, levels = c("PaaSc", "CelliD", "GSdensity")))

b_annot <- tibble(
  Method = factor(c("PaaSc", "CelliD", "GSdensity"), levels = c("PaaSc", "CelliD", "GSdensity")),
  label = sprintf("AUC = %.4f", c(B_roc["PaaSc_recomp"], B_roc["CelliD_recomp"], B_roc["GSdensity_arch"]))
)

p_b2 <- ggplot(b_long, aes(x = Method, y = Score, fill = condition)) +
  geom_boxplot(alpha = 0.8, outlier.size = 0.3, outlier.alpha = 0.35) +
  geom_text(data = b_annot, aes(x = Method, y = 1.07, label = label), inherit.aes = FALSE, fontface = "bold", size = 3.6) +
  scale_fill_manual(values = c("ctrl" = "grey65", "stim" = "#E64B35")) +
  coord_cartesian(ylim = c(0, 1.13)) +
  labs(title = paste0("B-cell Activation in B cells (Control, n = ", sum(b_idx), ")"),
       x = NULL, y = "Rank normalized score", fill = "Condition") +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.text.x = element_text(size = 11, face = "bold", color = "black"),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.9),
    panel.grid.major.x = element_blank()
  )

p_fig6b <- (p_b1 | p_b2) + plot_layout(guides = "collect") +
  plot_annotation(
    title = "Fig 6B: Pathway Condition Sensitivity vs. Lineage Specificity Benchmark",
    subtitle = "Recomputed PaaSc & CelliD from IFN.RDS vs. Author Archived GSdensity",
    theme = theme(
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 10, hjust = 0.5, color = "grey35")
    )
  ) & theme(legend.position = "right")

ggsave(file.path(out_dir, "Fig6B_ifn_vs_bcell_benchmark.png"), p_fig6b, width = 12.5, height = 5.5, dpi = 300, bg = "white")

# ------------------------------------------------------------------------------
# 8. Write Verification Summary
# ------------------------------------------------------------------------------
sink(file.path(out_dir, "fig6_verification_summary.txt"))
cat("PaaSc Figure 6 Full End-to-End Recomputation Summary\n")
cat("======================================================\n")
cat("Dataset: IFN.RDS (GSE96583 PBMC Ctrl vs IFN-beta Stim)\n")
cat("Dimensions:", nrow(IFN), "genes x", ncol(IFN), "cells | B cells:", sum(b_idx), "\n")
cat("IFN Pathway: MSigDB 26.1.1 GOBP_INTERFERON_MEDIATED_SIGNALING_PATHWAY (n =", length(INTERFERON),
    ", overlap =", sum(INTERFERON %in% rownames(IFN)), ")\n")
cat("Control Pathway: B_cell_activation.txt (n =", length(B_cell_activation),
    ", overlap =", sum(B_cell_activation %in% rownames(IFN)), ")\n\n")
cat("--- Correlation (Local Recomputed vs Author Archived Scores) ---\n")
cat(sprintf("  PaaSc INTERFERON Pearson r:          %.6f\n", cor_paasc_ifn))
cat(sprintf("  PaaSc B_cell_activation Pearson r:   %.6f\n", cor_paasc_b))
cat(sprintf("  CelliD INTERFERON Pearson r:         %.6f\n", cor_cellid_ifn))
cat(sprintf("  CelliD B_cell_activation Pearson r:  %.6f\n\n", cor_cellid_b))
cat("--- IFN Pathway ROC AUC (All 28,871 cells, ctrl vs stim) ---\n")
print(round(IFN_roc, 4))
cat("\n--- B-cell Activation ROC AUC (B cells only, ctrl vs stim negative control) ---\n")
print(round(B_roc, 4))
sink()

cat("=== Fig6 Verification Complete! Summary: ===\n")
readLines(file.path(out_dir, "fig6_verification_summary.txt")) %>% cat(sep = "\n")
