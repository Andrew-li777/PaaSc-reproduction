# ==============================================================================
# PaaSc 官方教程 Part 2: 空间转录组通路活性推断与空间互信息 (Spatial MI) 检验
# 
# 文献出处: Liao X, Hong Y, Feng Y, et al. PLOS Comput Biol 2025, 21(11): e1013666.
# 核心方法: SCTransform 空间预处理 -> computeMCA 降维 -> doRegression 背景去偏 -> calculateSpatialMI 3D 互信息置换检验
#
# 输入数据:
#   - SeuratData: stxBrain (type = "posterior1", 10x Visium 小鼠脑后部切片，31,053 基因 x 3,353 spots)
#   - msigdbr (v26.1.1): Mus musculus C5 GO:BP 全量生物过程背景基因集 (7,535 条通路)
#
# 输出文件:
#   - repro_results/06_spatial_clusters_umap_and_slice.png (切片组织聚类与 UMAP 分群对比图)
#   - repro_results/07_spatial_cocaine_pathway.png (可卡因响应通路皮层/纹状体空间特异性图)
#   - repro_results/08_spatial_cerebellar_pathways.png (小脑皮层与浦肯野细胞发育空间定位图)
#   - repro_results/spatial_MI_results.tsv (各通路空间互信息 MI、p-value 与 q-value 汇总表)
#   - repro_results/spatial_verification_summary.txt (空间量化验收摘要)
#
# 实测环境: R 4.6.0 + Seurat 5.6.0 + SeuratData 0.2.2 + msigdbr 26.1.1 + PaaSc 1.0.1
# 预期耗时: 约 60-80 秒 (SCT 约 38 秒，MCA 约 18 秒，MI 置换约 5 秒)
# ==============================================================================

set.seed(42)

suppressPackageStartupMessages({
  library(ggplot2)
  library(ggrepel)
  library(patchwork)
  library(msigdbr)
  library(Seurat)
  library(SeuratData)
  library(PaaSc)
})

# 创建结果输出目录
out_dir <- "repro_results"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# ------------------------------------------------------------------------------
# 步骤 1: 加载 10x Visium 小鼠脑后部组织切片 (posterior1)
# ------------------------------------------------------------------------------
cat("[1/6] Loading stxBrain (posterior1) dataset...\n")
brain <- LoadData("stxBrain", type = "posterior1")
print(brain) # 31,053 features across 3,353 samples within 1 assay

# ------------------------------------------------------------------------------
# 步骤 2: SCTransform 空间方差稳定化变换、PCA 降维、邻近图构建与分群
# ------------------------------------------------------------------------------
cat("[2/6] Running SCTransform, PCA, clustering & UMAP...\n")
t0 <- Sys.time()
brain <- SCTransform(brain, assay = "Spatial", verbose = FALSE)
brain <- RunPCA(brain, assay = "SCT", verbose = FALSE)
brain <- FindNeighbors(brain, reduction = "pca", dims = 1:30, verbose = FALSE)
brain <- FindClusters(brain, verbose = FALSE)
brain <- RunUMAP(brain, reduction = "pca", dims = 1:30, verbose = FALSE)
cat("Spatial preprocessing elapsed:", round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 2), "sec\n")

# 图 6: 切片组织聚类与 UMAP 分群并排图
p_umap <- DimPlot(brain, reduction = "umap", label = TRUE) + ggtitle("stxBrain UMAP Clusters")
p_slice <- SpatialDimPlot(brain, label = TRUE, label.size = 4) + ggtitle("stxBrain Tissue Clusters")
ggsave(file.path(out_dir, "06_spatial_clusters_umap_and_slice.png"), p_umap + p_slice, width = 12, height = 5.5, dpi = 200)

# ------------------------------------------------------------------------------
# 步骤 3: 提取 MSigDB 小鼠全量 GO:BP 背景通路库
# 注: msigdbr 26.x 参数为 collection = 'C5', subcollection = 'GO:BP'
# ------------------------------------------------------------------------------
cat("[3/6] Fetching MSigDB Mus musculus GO:BP gene sets...\n")
mm_c5_gobp <- msigdbr(species = "Mus musculus", collection = "C5", subcollection = "GO:BP")
background_geneset <- split(mm_c5_gobp$gene_symbol, mm_c5_gobp$gs_name)
background_geneset <- lapply(background_geneset, unique)
cat("Total GO:BP pathways loaded:", length(background_geneset), "\n") # 7,535 条通路

# ------------------------------------------------------------------------------
# 步骤 4: 基于 SCT assay 表达矩阵运行 PaaSc MCA 联合降维
# ------------------------------------------------------------------------------
cat("[4/6] Running computeMCA (nmcs = 20, assay = 'SCT')...\n")
t1 <- Sys.time()
brain <- computeMCA(brain, nmcs = 20, assay = "SCT")
cat("computeMCA elapsed:", round(as.numeric(difftime(Sys.time(), t1, units = "secs")), 2), "sec\n")

# ------------------------------------------------------------------------------
# 步骤 5: 选取 8 条代表性测试通路进行去偏评分
# 涵盖: 纹状体/皮层神经通路 (可卡因响应)、小脑结构发育通路、通用代谢通路及阴性对照 (生物生长)
# ------------------------------------------------------------------------------
cat("[5/6] Scoring 8 test pathways with PaaSc...\n")
example_pathways <- c(
  "GOBP_RESPONSE_TO_COCAINE",
  "GOBP_CEREBELLAR_CORTEX_FORMATION",
  "GOBP_CEREBELLAR_CORTEX_DEVELOPMENT",
  "GOBP_CEREBELLAR_CORTEX_MORPHOGENESIS",
  "GOBP_CEREBELLAR_PURKINJE_CELL_LAYER_DEVELOPMENT",
  "GOBP_CGMP_MEDIATED_SIGNALING",
  "GOBP_PROTEIN_KINASE_B_SIGNALING",
  "GOBP_MULTICELLULAR_ORGANISM_GROWTH"
)
# 使用 intersect 防御部分历史版本 MSigDB 通路命名更迭
pathway_geneset <- background_geneset[intersect(example_pathways, names(background_geneset))]
gene_rate <- getGeneRate(background.geneset = background_geneset, pathway.geneset = pathway_geneset)
regression_data <- doRegression(brain, gene.rate = gene_rate)
score_data <- computeScore(brain, regression.data = regression_data, pvalue = 0.05, weight = FALSE, normalize = "z-score")
score_data <- score_data[, colSums(is.na(score_data)) == 0, drop = FALSE]
cat("Valid scored pathways:", ncol(score_data), "\n")

# ------------------------------------------------------------------------------
# 步骤 6: 空间互信息 (Spatial Mutual Information) 置换检验
# 原理: 将空间二维坐标 (x, y) 与通路活性得分 z 进行多维分箱，
# 计算三维互信息 I(X, Y; Z) 并执行 1,000 次空间坐标重抽样置换，评估空间非均匀聚集显著性。
# 关键防护: 强制使用 coords[rownames(score_data), ] 进行行名显式对齐，杜绝错行风险。
# ------------------------------------------------------------------------------
cat("[6/6] Computing spatial mutual information (calculateSpatialMI)...\n")
coords <- GetTissueCoordinates(brain[["posterior1"]])
coords <- coords[rownames(score_data), , drop = FALSE]
spatialMI <- calculateSpatialMI(x = coords$x, y = coords$y, z = score_data)
spatialMI <- data.frame(pathway = rownames(spatialMI), spatialMI, qvalue = p.adjust(spatialMI$pvalue))
print(spatialMI)
write.table(spatialMI, file.path(out_dir, "spatial_MI_results.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

# ------------------------------------------------------------------------------
# 步骤 7: 空间原位切片通路活性图渲染
# ------------------------------------------------------------------------------
for (pw in colnames(score_data)) {
  brain@meta.data[[pw]] <- score_data[colnames(brain), pw]
}

# 图 7: 可卡因响应通路空间分布图 (显著富集于皮层和纹状体区域，MI = 1.6655, q = 0)
if ("GOBP_RESPONSE_TO_COCAINE" %in% colnames(score_data)) {
  p_coc <- SpatialFeaturePlot(brain, features = "GOBP_RESPONSE_TO_COCAINE", pt.size.factor = 2.5) +
    ggtitle("GOBP_RESPONSE_TO_COCAINE")
  ggsave(file.path(out_dir, "07_spatial_cocaine_pathway.png"), p_coc, width = 6.5, height = 5.5, dpi = 200)
}

# 图 8: 小脑皮层形成与浦肯野细胞发育通路 (高特异性局限在后部小脑叶，MI = 1.5349, q = 0)
cerebellar_pws <- intersect(
  c("GOBP_CEREBELLAR_CORTEX_FORMATION", "GOBP_CEREBELLAR_PURKINJE_CELL_LAYER_DEVELOPMENT"),
  colnames(score_data)
)
if (length(cerebellar_pws) > 0) {
  p_cer <- SpatialFeaturePlot(brain, features = cerebellar_pws, pt.size.factor = 2.5)
  ggsave(file.path(out_dir, "08_spatial_cerebellar_pathways.png"), p_cer, width = 12, height = 5.5, dpi = 200)
}

# ------------------------------------------------------------------------------
# 步骤 8: 空间量化验收摘要输出
# ------------------------------------------------------------------------------
summary_lines <- c(
  paste0("Spatial_features: ", nrow(brain)),
  paste0("Spatial_spots: ", ncol(brain)),
  paste0("GOBP_pathways_total: ", length(background_geneset)),
  paste0("Valid_scored_pathways: ", ncol(score_data)),
  paste0("Significant_spatial_pathways_q0.05: ", sum(spatialMI$qvalue < 0.05)),
  paste0("Top_MI_pathway: ", spatialMI$pathway[which.max(spatialMI$MI)], " (MI=", round(max(spatialMI$MI), 4), ")")
)
writeLines(summary_lines, file.path(out_dir, "spatial_verification_summary.txt"))
cat("=== Spatial Reproduction Completed Successfully ===\n")
cat(paste(summary_lines, collapse = "\n"), "\n")
