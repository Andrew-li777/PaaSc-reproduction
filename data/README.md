# PaaSc 内置与教程数据集说明文档 (`data/`)

本目录存放 GitHub 官方仓库 `yoyoong/PaaSc` 自带的示例数据与预定义基因集文件，主要用于支持官方教程 `01_pbmc3k_repro.R` 与通用通路去偏评分。

---

## 1. 目录文件清单与元信息

| 文件/子目录 | 格式 | 规模 / 维度 | 生物学来源与背景 | 在复现工程中的用途 |
|---|---|---|---|---|
| [`pbmc3k_filtered_gene_bc_matrices/hg19/`](./pbmc3k_filtered_gene_bc_matrices/hg19) | 10x 稀疏矩阵 (`mtx` + `tsv`) | 原始过滤后 `13,714` 基因 × `2,700` 单细胞 | 10x Genomics 官方标准 PBMC 3k 测序数据（健康人外周血单核细胞） | 用于 `01_pbmc3k_repro.R` 和 `06_fig2_and_supfig2_repro.R` 的单细胞表达矩阵输入 |
| ├── `matrix.mtx` | Matrix Market | 30.47 MB | 稀疏表达矩阵（UMI counts） | `Read10X()` 解析计数主体 |
| ├── `genes.tsv` | TSV (无表头) | 830 KB (`32,738` 基因 Ensembl + Symbol) | 人类基因参考注释 (hg19) | 基因行名对齐 |
| └── `barcodes.tsv` | TSV (无表头) | 47 KB (`2,700` 细胞条形码) | 测序液滴细胞 Barcode | 细胞列名对齐 |
| [`pbmc3k_seurat_annotation.txt`](./pbmc3k_seurat_annotation.txt) | TSV 表格 (带表头) | 72 KB (`2,701` 行，含表头) | Seurat 官方标准聚类注释（9 种亚型：B 细胞 `344`、CD4 T `1,144`、CD8 T `316`、CD14+ Mono `480`、FCGR3A+ Mono `150`、NK `155`、DC `32`、Platelet `14`、Dendritic `66`） | 作为独立真值标签，用于在 `01_pbmc3k_repro.R` 中验证 B 细胞激活通路识别能力（`AUC = 0.9981`） |
| [`B_cell_activation.txt`](./B_cell_activation.txt) | 纯文本换行 | 1.94 KB (`278` 个基因 Symbol) | 人类 B 细胞激活（B cell activation）通路标记基因集合 | 目标功能通路基因集，与 `pbmc3k` 重叠 233 个基因，与 `IFN.RDS` 重叠 233 个基因 |
| [`kegg_legacy.gmt`](./kegg_legacy.gmt) | GMT 标准基因集格式 | 98 KB (`186` 条通路) | MSigDB KEGG Legacy 经典代谢与信号通路库 | 作为 `getGeneRate()` 的全量背景通路集合，用于计算每个基因的泛在出现频率，支撑多元线性回归去偏 |

---

## 2. 外部动态调用的辅助数据集

除本目录下的静态文件外，`PaaSc` 复现工程还在运行时动态调用了以下两项官方扩展数据：

1. **`stxBrain` 小鼠脑空间转录组数据集**：
   - **加载方式**：通过 `SeuratData::LoadData("stxBrain", type = "posterior1")` 调用（依赖包 `stxBrain.SeuratData v0.1.2`）。
   - **数据特征**：10x Visium 小鼠脑矢状后部切片，SCTransform 预处理后包含 `17,668` 个高变特征与 `3,353` 个组织 spots。
   - **使用脚本**：[02_spatial_stxBrain_repro.R](../repro_scripts/02_spatial_stxBrain_repro.R)。
2. **`MSigDB` 全量 GO:BP 通路库**：
   - **加载方式**：通过 `msigdbr::msigdbr(species = "Mus musculus", collection = "C5", subcollection = "GO:BP")` 提取。
   - **数据特征**：涵盖小鼠全部 `7,535` 条 Biological Process 通路。
   - **使用脚本**：[02_spatial_stxBrain_repro.R](../repro_scripts/02_spatial_stxBrain_repro.R)（空间互信息检验）与 [03_fig6_ifn_repro.R](../repro_scripts/03_fig6_ifn_repro.R)（人类 `GOBP_INTERFERON_MEDIATED_SIGNALING_PATHWAY`，113 个基因）。
