# ==============================================================================
# PaaSc 官方教程 Part 1: 单细胞 scRNA-seq 通路活性识别与验证 (pbmc3k)
# 
# 文献出处: Liao X, Hong Y, Feng Y, et al. PLOS Comput Biol 2025, 21(11): e1013666.
# 算法核心: MCA 细胞-基因联合降维 -> 背景频率线性回归去偏 -> GMM 双峰二分类 -> ROC-AUC 验证
#
# 输入数据:
#   - data/pbmc3k_filtered_gene_bc_matrices/hg19 (10x Genomics 原始稀疏计数)
#   - data/kegg_legacy.gmt (KEGG Legacy 背景通路集，186条通路)
#   - data/B_cell_activation.txt (目标验证通路: B细胞激活，278个基因)
#   - data/pbmc3k_seurat_annotation.txt (独立参考注释: 344个B细胞 vs 2,356个非B细胞)
#
# 输出文件:
#   - repro_results/01_umap_B_cell_activation_score.png (UMAP 通路活性评分图)
#   - repro_results/02_hist_B_cell_activation_gmm.png (GMM 二分类双峰直方图及阈值线)
#   - repro_results/03_umap_gmm_vs_celltypes.png (PaaSc GMM 状态 vs 参考细胞类型对比图)
#   - repro_results/04_boxplot_score_by_celltype.png (各细胞类型通路活性分布箱线图)
#   - repro_results/05_roc_B_cell_activation.png (ROC 曲线，AUC = 0.9981)
#   - repro_results/verification_summary.txt (量化验收摘要)
#
# 实测环境: R 4.6.0 + Seurat 5.6.0 + SeuratObject 5.4.0 + PaaSc 1.0.1
# 预期耗时: 约 25-35 秒 (单线程)
# ==============================================================================

set.seed(42)

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(patchwork)
  library(PaaSc)
  library(ROCit)
})

# 创建结果输出目录
out_dir <- "repro_results"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# ------------------------------------------------------------------------------
# 步骤 1: 读取 10x Genomics 原始稀疏矩阵并构建 Seurat 基础对象
# ------------------------------------------------------------------------------
cat("[1/6] Loading pbmc3k 10x matrices...\n")
data_mat <- Read10X(data.dir = "data/pbmc3k_filtered_gene_bc_matrices/hg19")
object <- CreateSeuratObject(counts = data_mat, min.cells = 3, min.features = 200)
print(object) # 13,714 features x 2,700 cells

# ------------------------------------------------------------------------------
# 步骤 2: 多重对应分析 (MCA) 联合降维
# 原理: 将细胞与基因同时投影到同一低维空间 (nmcs = 20)，以距离衡量表达关联
# ------------------------------------------------------------------------------
cat("[2/6] Running computeMCA (nmcs = 20)...\n")
t0 <- Sys.time()
object <- computeMCA(object, nmcs = 20)
cat("computeMCA elapsed:", round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 2), "sec\n")

# ------------------------------------------------------------------------------
# 步骤 3: 背景通路基因频率统计与多元线性回归去偏
# 原理: 统计目标通路和全部背景通路集合中的基因出现频率 (getGeneRate)；
# 以各维度 MCA 基因 loading 为因变量，以背景频率与通路成员标志为自变量做回归 (doRegression)，
# 提取剔除背景高频看家基因偏倚后的显著特征维度 t 统计量。
# ------------------------------------------------------------------------------
cat("[3/6] Calculating gene rates & MCA regression...\n")
background_geneset <- "data/kegg_legacy.gmt"
B_cell_activation <- readLines("data/B_cell_activation.txt")
pathway_geneset <- list(B_cell_activation = B_cell_activation)
gene_rate <- getGeneRate(background.geneset = background_geneset, pathway.geneset = pathway_geneset)
regression_data <- doRegression(object, gene.rate = gene_rate)
print(round(regression_data[, 1:5], 4))

# ------------------------------------------------------------------------------
# 步骤 4: 计算单细胞通路活性评分与高斯混合模型 (GMM) 自动二分类
# 原理: computeScore 仅使用回归显著维度 (p < 0.05) 进行加权投影与 z-score 标准化；
# doBinarization 使用 mixtools::normalmixEM 拟合双峰分布，阈值取两个高斯均值的中点。
# ------------------------------------------------------------------------------
cat("[4/6] Computing pathway activity scores & GMM binarization...\n")
score_data <- computeScore(object, regression.data = regression_data, weight = FALSE, normalize = "z-score")
label_data <- doBinarization(score_data, n.cluster = 2, method = "GMM")
print(table(label_data$B_cell_activation.label))

# ------------------------------------------------------------------------------
# 步骤 5: Seurat 标准 UMAP 聚类与结果图绘制
# ------------------------------------------------------------------------------
cat("[5/6] Running UMAP & generating figures...\n")
object <- FindVariableFeatures(object, verbose = FALSE)
object <- ScaleData(object, verbose = FALSE)
object <- RunPCA(object, verbose = FALSE)
object <- FindNeighbors(object, dims = 1:20, verbose = FALSE)
object <- RunUMAP(object, dims = 1:20, verbose = FALSE)

# 极值平滑截断 (1% 至 99% 分位数) 便于热图连续展示
cscore <- score_data[colnames(object), 1]
cscore[cscore > quantile(cscore, 0.99)] <- quantile(cscore, 0.99)
cscore[cscore < quantile(cscore, 0.01)] <- quantile(cscore, 0.01)
object@meta.data$B_cell_activation <- cscore
object@meta.data$B_cell_activation_label <- label_data[colnames(object), ]

# 读取独立参考注释文件 (pbmc3k_seurat_annotation.txt)
annotation_file <- "data/pbmc3k_seurat_annotation.txt"
seurat_annotations <- read.table(annotation_file, header = TRUE, row.names = 1, sep = "\t")
object@meta.data$seurat_annotations <- seurat_annotations[colnames(object), ]

# 图 1: UMAP 通路活性评分图
p1 <- FeaturePlot(object, features = "B_cell_activation", cols = c("white", "darkred")) +
  ggtitle("PaaSc B-cell Activation Score (pbmc3k)")
ggsave(file.path(out_dir, "01_umap_B_cell_activation_score.png"), p1, width = 6.5, height = 5.5, dpi = 200)

# 图 2: GMM 二分类双峰直方图及判定阈值线
threshold <- min(score_data[rownames(label_data)[label_data[["B_cell_activation.label"]] == "positive"], "B_cell_activation"])
p2 <- ggplot(score_data, aes(x = B_cell_activation)) +
  geom_histogram(bins = 50, fill = "grey70", alpha = 0.8) +
  geom_vline(xintercept = threshold, color = "red", linewidth = 0.9) +
  annotate("text", x = threshold, y = 50, label = paste0("threshold = ", signif(threshold, 3)),
           hjust = -0.05, color = "black", size = 4.5) +
  xlab("Activity score (z-score)") + ylab("Cell frequency") +
  ggtitle("GMM Binarization Threshold") + theme_minimal()
ggsave(file.path(out_dir, "02_hist_B_cell_activation_gmm.png"), p2, width = 6.5, height = 4.5, dpi = 200)

# 图 3: PaaSc GMM 状态 vs 参考细胞类型对比图
p3a <- DimPlot(object, group.by = "B_cell_activation_label") + ggtitle("PaaSc GMM State")
p3b <- DimPlot(object, group.by = "seurat_annotations") + ggtitle("Reference Cell Types")
ggsave(file.path(out_dir, "03_umap_gmm_vs_celltypes.png"), p3a + p3b, width = 12, height = 5, dpi = 200)

# 图 4: 各细胞类型通路活性分布箱线图
score_data_withanno <- cbind(score_data, seurat_annotations = seurat_annotations[rownames(score_data), , drop = FALSE])
p4 <- ggplot(score_data_withanno, aes(x = seurat_annotations, y = B_cell_activation, color = seurat_annotations)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.2, size = 0.2, alpha = 0.5) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        legend.position = "none",
        axis.text.x = element_text(vjust = 1, hjust = 1, angle = 45, size = 11)) +
  ylab("B_cell_activation score") + xlab("Cell Type") +
  ggtitle("Pathway Score Distribution Across Cell Types")
ggsave(file.path(out_dir, "04_boxplot_score_by_celltype.png"), p4, width = 7, height = 5, dpi = 200)

# ------------------------------------------------------------------------------
# 步骤 6: 对照独立细胞注释计算经验 ROC-AUC
# ------------------------------------------------------------------------------
cat("[6/6] Validating ROC-AUC against B cell annotations...\n")
label <- rep("others", nrow(score_data_withanno))
label[score_data_withanno$seurat_annotations == "B"] <- "B"
roc_empirical <- rocit(score = score_data_withanno$B_cell_activation, class = label, negref = "others")
print(summary(roc_empirical))

# 图 5: ROC 曲线
png(file.path(out_dir, "05_roc_B_cell_activation.png"), width = 1200, height = 1000, res = 200)
plot(roc_empirical, YIndex = FALSE, values = FALSE, col = c(2, 1), main = "ROC: B-cell Activation Score vs B Cells")
text(0.25, 0.55, paste0("AUC = ", round(roc_empirical$AUC, 4)), cex = 1.2)
dev.off()

# ------------------------------------------------------------------------------
# 步骤 7: 输出标准化量化验收摘要
# ------------------------------------------------------------------------------
summary_lines <- c(
  paste0("Features: ", nrow(object)),
  paste0("Cells: ", ncol(object)),
  paste0("B_cells_positive: ", sum(label == "B")),
  paste0("Others_negative: ", sum(label == "others")),
  paste0("GMM_threshold: ", signif(threshold, 4)),
  paste0("ROC_AUC: ", round(roc_empirical$AUC, 4)),
  paste0("PaaSc_version: ", as.character(packageVersion("PaaSc"))),
  paste0("Seurat_version: ", as.character(packageVersion("Seurat"))),
  paste0("SeuratObject_version: ", as.character(packageVersion("SeuratObject")))
)
writeLines(summary_lines, file.path(out_dir, "verification_summary.txt"))
writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo.txt"))
cat("=== Reproduction Completed Successfully ===\n")
cat(paste(summary_lines, collapse = "\n"), "\n")
