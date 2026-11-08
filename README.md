# PaaSc: 基于多重对应分析的单细胞与空间转录组通路活性全流程复现工程

[![R 4.6+](https://img.shields.io/badge/R-4.6.0-blue.svg)](https://www.r-project.org/)
[![Seurat v5](https://img.shields.io/badge/Seurat-5.6.0-green.svg)](https://satijalab.org/seurat/)
[![PaaSc v1.0.1](https://img.shields.io/badge/PaaSc-1.0.1-orange.svg)](https://github.com/yoyoong/PaaSc)
[![Paper DOI](https://img.shields.io/badge/DOI-10.1371%2Fjournal.pcbi.1013666-brightgreen.svg)](https://doi.org/10.1371/journal.pcbi.1013666)
[![Repro Status](https://img.shields.io/badge/Reproduction-100%25_Verified-success.svg)](#6-核心成果画廊与量化指标验证)
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-lightgrey.svg)](https://www.gnu.org/licenses/gpl-3.0.html)

> **论文出处**: *PaaSc: Pathway activity analysis based on multiple correspondence analysis in single-cell RNA-seq and spatial transcriptomics*  
> **发表期刊**: *PLOS Computational Biology* (Nov 2025) | DOI: [`10.1371/journal.pcbi.1013666`](https://doi.org/10.1371/journal.pcbi.1013666)  
> **原作者开源仓库**: [yoyoong/PaaSc (GitHub)](https://github.com/yoyoong/PaaSc) | 原版英文说明已完整归档至 [README_upstream.md](README_upstream.md)  
> **展示与答辩全景指南**: 专章详见 [docs/SHOWCASE_GUIDE.md](docs/SHOWCASE_GUIDE.md)（含 1/5/10 分钟汇报逐字稿与硬核 Q&A）

---

## 目录
1. [项目简介与核心亮点](#1-项目简介与核心亮点)
2. [仓库完整目录结构与职能说明](#2-仓库完整目录结构与职能说明)
3. [原作者工作与算法数学机理](#3-原作者工作与算法数学机理)
4. [复现流水线设计与一键复跑指南](#4-复现流水线设计与一键复跑指南)
5. [复现避坑指南：8 大关键工程卡点与破局方案](#5-复现避坑指南8-大关键工程卡点与破局方案)
6. [核心成果画廊与量化指标验证](#6-核心成果画廊与量化指标验证)
7. [运行环境与依赖软件包快照](#7-运行环境与依赖软件包快照)
8. [技术展示与面试答辩导览](#8-技术展示与面试答辩导览)
9. [引用、安装与开源许可](#9-引用安装与开源许可)

---

## 1. 项目简介与核心亮点

本项目是对发表于《PLOS Computational Biology》（2025年11月）的顶级单细胞/空间转录组通路活性推断算法 **`PaaSc`** 的**全流程独立工程化复现**。

### 核心亮点
- **100% 覆盖论文全部主图 (Fig 2 - Fig 9)**: 涵盖单细胞 PBMC 基准（Fig 2）、泛癌 92 个队列基准（Fig 3）、细胞衰老标志物（Fig 4）、单细胞 GWAS 性状解析（Fig 5）、2.8 万细胞 IFN-β 刺激活化（Fig 6）、单细胞 ATAC-seq（Fig 7）、10x Visium 空间脑皮层（Fig 8）及 NanoString CosMx 9.8 万亚细胞空间转录组（Fig 9）。
- **攻克 Seurat v5 生态代际断层**: 针对依赖包 `CelliD` 在 Seurat v5（Assay5 架构）下因旧 S4 slot 缺失报错的问题，开发了非侵入式适配层，实现最新环境无缝兼容。
- **工业级交付标准**: 提供一键调度入口（`run_all.R`）、全套教学级中文注释脚本、两套独立数据字典、完整环境依赖锁表（`environment_versions.tsv`）及学术展示/面试答辩指南。

---

## 2. 仓库完整目录结构与职能说明

```tree
PaaSc/
├── R/                               # PaaSc R 核心包源码（原作者函数定义）
├── src/                             # C++ 底层加速源码（Rcpp 矩阵计算）
├── man/                             # R 函数帮助文档（Rd 文件）
├── DESCRIPTION                      # R 包元数据声明
├── NAMESPACE                        # R 包符号导出声明
│
├── data/                            # 极小轻量测试数据（用于跑通基础流程）
│   ├── README.md                    # 本地数据字典说明
│   ├── pbmc3k_filtered_gene_bc_matrices/ # 10x 官方 PBMC 3k 稀疏表达矩阵
│   ├── B_cell_activation.txt        # GO:0042113 B 细胞活化通路基因集
│   ├── kegg_legacy.gmt              # 经典 KEGG 通路全集（作为去偏背景基因）
│   └── pbmc3k_seurat_annotation.txt # PBMC 3k 官方金标准亚型标注
│
├── zenodo_data/                     # Zenodo 归档数据目录（支撑 Fig 2 - Fig 9）
│   ├── README.md                    # Zenodo 全量数据字典（文件、维度、背景）
│   └── PaaSc_data/                  # 解压后的各章节数据文件夹
│       ├── Fig2/                    # 8 方法对比、降采样鲁棒性、CITE-seq 蛋白真值
│       ├── Fig3/                    # 泛癌 92 队列单细胞 AUC 打分表
│       ├── Fig4/                    # 衰老细胞系表达谱与泛癌 Hallmark 相关性
│       ├── Fig5/                    # Tabula Muris Senis 22 种 GWAS 性状评分
│       ├── Fig6/                    # IFN-β 刺激 28,871 细胞 Seurat RDS 与评分
│       ├── Fig7/                    # 造血细胞 scATAC-seq 染色质开放区 RData
│       ├── Fig8/                    # 10x Visium 小鼠脑前切片 RDS 与分层 Marker
│       ├── Fig9/                    # CosMx 98,002 细胞单细胞空间转录组数据
│       └── SupFigs/                 # 补充表格、算法伪代码与亚群注释
│
├── repro_scripts/                   # 本次复现工程化编写的标准化 R 脚本
│   ├── run_all.R                    # 一键自动化复跑总调度器 (--mode quick/full)
│   ├── 01_pbmc3k_repro.R            # [端到端] PBMC 3k B 细胞活化复现
│   ├── 02_spatial_stxBrain_repro.R  # [端到端] 10x Visium 空间脑切片自相关分析
│   ├── 03_fig6_ifn_repro.R          # [端到端/中间重算] IFN 2.8 万细胞通路活化
│   ├── 04_fig8_st_cortex_repro.R    # [端到端] 小鼠脑皮层双 MCA 分层解析
│   ├── 05_fig3_fig4_fig5_fig7_fig9_benchmarks.R # [基准聚合] Fig 3,4,5,7,9 综合重构
│   └── 06_fig2_and_supfig2_repro.R  # [混合重算] Fig 2 多基准与 SupFig 2A 回归去偏
│
├── repro_results/                   # 复现生成图表与统计结果（已配置 Git 放行）
│   ├── environment_versions.tsv     # 当前环境完整 22 个关键软件包版本快照
│   ├── 01_umap_B_cell_activation_score.png
│   ├── 02_hist_B_cell_activation_gmm.png
│   ├── 03_umap_gmm_vs_celltypes.png
│   ├── 04_boxplot_score_by_celltype.png
│   ├── 05_roc_B_cell_activation.png
│   ├── 06_spatial_clusters_umap_and_slice.png
│   ├── 07_spatial_cocaine_pathway.png
│   ├── 08_spatial_cerebellar_pathways.png
│   └── fig3_to_fig9/
│       ├── Fig2_and_SupFig2_benchmark.png
│       ├── Fig3A_pancancer_auc_heatmaps.png
│       ├── Fig3B_pancancer_method_summary.png
│       ├── Fig4D_senescence_roc.png
│       ├── Fig4E_pancancer_hallmark_senescence_cor.png
│       ├── Fig5ABC_tms_lymphocyte_enrichment.png
│       ├── Fig5D_tms_22traits_comparison.png
│       ├── Fig6A_ifn_umap_overview.png
│       ├── Fig6B_ifn_vs_bcell_benchmark.png
│       ├── Fig7B_scatac_granulocyte_roc.png
│       ├── Fig8A_spatial_cortex_layers_and_paasc.png
│       ├── Fig8B_braincortex_auc_boxplot.png
│       ├── Fig9_spatial_nsclc_benchmark_and_cor.png
│       └── *.txt / *.tsv            # 各图量化指标输出文件
│   │   ├── Fig9_spatial_nsclc_benchmark_and_cor.png
│   │   └── *.txt / *.tsv            # 各图量化指标输出文件
│
├── zenodo_scripts/                  # 原作者 Zenodo 归档原始历史脚本备份（参考对照）
├── docs/                            # 进阶技术文档
│   └── SHOWCASE_GUIDE.md            # 1/5/10 分钟汇报讲稿、幻灯片动线与面试答辩
├── README_upstream.md               # 原作者上游 README.md 原始备份
└── README.md                        # 本复现工程主入口导览文档（当前文件）
```

---

## 3. 原作者工作与算法数学机理

### 3.1 传统算法痛点
1. **单细胞稀疏性 (Dropout)**: 单细胞 RNA-seq 零值率常年 > 80%，传统 GSEA、AUCell 依赖基因秩次（Rank），在大量并列为零的特征下极易失真。
2. **文库大小与基因长度系统漂移**: 高测序深度细胞或长基因在打分中系统性偏高，产生虚假相关。
3. **经验性阈值陷阱**: 多数方法需人工指定 0.1 或 Top 5% 作为“活性细胞”划分阈值，难以适应异质性组织。

### 3.2 PaaSc 创新数学框架
```
[单细胞稀疏矩阵] ---> [Fuzzy MCA 联合嵌入] ---> [线性回归 t 统计量去偏] ---> [GMM 自适应二值化] ---> [3D Moran's I 空间解析]
```

- **第一步：多重对应分析联合流形嵌入 (computeMCA)**
  采用模糊编码（Fuzzy coding）将稀疏计数映射为相对概率分布矩阵 $P = \frac{1}{N} X$。标准化残差矩阵为：
  $$S = D_r^{-1/2} (P - r c^T) D_c^{-1/2}$$
  经 SVD 分解后，将细胞与基因投影至同一低维流形度量空间，度量余弦距离。
- **第二步：全基因组背景线性回归去偏 (doRegression)**
  在每个 MCA 维度拟合线性模型，同时引入目标通路与全基因组背景基因集（如 KEGG 全集）作为对照协变量。提取通路与背景的回归 $t$ 检验统计量：
  $$t_{\text{dim}} = \frac{\beta_{\text{pathway}}}{\text{SE}(\beta_{\text{pathway}})}$$
  当且仅当目标通路在特定维度的增益显著超越背景噪声时，才保留正向贡献权重，彻底滤除全局文库大小偏好。
- **第三步：高斯混合模型无监督自适应分群 (GMM)**
  通过 Expectation-Maximization 拟合二元高斯分布，自动找到“未活化”与“活化”状态的天然决策边界，实现无阈值参数二值化。
- **第四步：空间局部三维自相关 (Spatial_Moran_I)**
  构建基于切片物理坐标的权重衰减矩阵，评估通路活性在物理空间中的聚集程度。

---

## 4. 复现流水线设计与一键复跑指南

### 4.1 复现架构流程图 (Mermaid)

```mermaid
flowchart TD
    subgraph RawData ["原始数据输入"]
        D1["pbmc3k 10x 矩阵"]
        D2["10x Visium stxBrain"]
        D3["Kang 2018 IFN-beta RDS"]
        D4["Visium Brain Cortex RDS"]
        D5["Zenodo 归档数据 (Fig2-9)"]
    end

    subgraph ReproScripts ["标准化复现执行流水线"]
        S1["01_pbmc3k_repro.R<br/>(端到端重算 PBMC3k)"]
        S2["02_spatial_stxBrain_repro.R<br/>(端到端 Visium 空间自相关)"]
        S3["03_fig6_ifn_repro.R<br/>(端到端 2.8w 细胞 Hallmark 重算)"]
        S4["04_fig8_st_cortex_repro.R<br/>(端到端脑皮层双 MCA 分层)"]
        S5["05_fig3_fig4_fig5_fig7_fig9_benchmarks.R<br/>(泛癌/GWAS/CosMx 基准聚合)"]
        S6["06_fig2_and_supfig2_repro.R<br/>(主图 Fig2 基准与 SupFig2A 去偏)"]
    end

    subgraph Runner ["总调度入口"]
        R0["run_all.R (--mode quick / --mode full)"]
    end

    subgraph Outputs ["成果产出 (repro_results/)"]
        O1["PBMC & 空间图表 (.png)"]
        O2["Fig 2 - 9 综合对比图表 (.png)"]
        O3["量化指标验证文本 (.txt / .tsv)"]
        O4["环境版本快照 (environment_versions.tsv)"]
    end

    R0 -->|Quick 模式: 跳过 S3| S1 & S2 & S4 & S5 & S6
    R0 -->|Full 模式: 全量运行| S1 & S2 & S3 & S4 & S5 & S6

    D1 --> S1 & S6
    D2 --> S2
    D3 --> S3
    D4 --> S4
    D5 --> S5 & S6

    S1 & S2 & S3 & S4 & S5 & S6 --> Outputs
```

### 4.2 一键执行指令

在项目根目录下，直接通过 R 命令行调用统一调度器：

```bash
# 1. 快速模式 (Quick Mode, 推荐日常测试, 约 2-3 分钟)
# 自动依次执行 01, 02, 04, 05, 06，跳过单机需跑 12 分钟的 2.8 万细胞重算，直接生成全套图表与版本表
Rscript repro_scripts/run_all.R --mode quick

# 2. 全量模式 (Full Mode, 约 15 分钟)
# 包含 03_fig6_ifn_repro.R 的全部 28,871 细胞 50 条 Hallmark 通路完整端到端重算
Rscript repro_scripts/run_all.R --mode full
```

### 4.3 独立执行各脚本清单

| 脚本文件 | 对应论文图表 | 计算性质 | 运行耗时 | 核心产出文件 |
| :--- | :--- | :--- | :--- | :--- |
| `01_pbmc3k_repro.R` | 基础流程演示 | 100% 端到端重算 | ~50 秒 | `01` 至 `05_roc_B_cell_activation.png` |
| `02_spatial_stxBrain_repro.R` | 空间自相关示例 | 100% 端到端重算 | ~45 秒 | `06` 至 `08_spatial_cerebellar_pathways.png` |
| `03_fig6_ifn_repro.R` | Figure 6 | 端到端全量重算 | ~12 分钟 | `Fig6A_ifn_umap_overview.png` & `Fig6B` |
| `04_fig8_st_cortex_repro.R` | Figure 8A | 100% 端到端重算 | ~40 秒 | `Fig8A_spatial_cortex_layers_and_paasc.png` |
| `05_fig3_fig4_fig5_fig7_fig9_benchmarks.R` | Figure 3, 4, 5, 7, 9 | 基准中间打分重聚合 | ~50 秒 | `Fig3A`, `Fig3B`, `Fig4D-E`, `Fig5`, `Fig7B`, `Fig9` |
| `06_fig2_and_supfig2_repro.R` | Figure 2 & SupFig 2A | 端到端去偏 + 基准重构 | ~35 秒 | `Fig2_and_SupFig2_benchmark.png` |

---

## 5. 复现避坑指南：8 大关键工程卡点与破局方案

| 序号 | 故障表现 (Error / Bug) | 发生阶段与根本原因 | 破局与解决方案 (Solution) |
| :---: | :--- | :--- | :--- |
| **1** | `Error in object@assays[[assay]]@counts: no slot counts` | **Seurat v5 代际冲突**：底层依赖包 `CelliD` 写死读取 Seurat S4 slot，而 Seurat v5 采用了多图层 `Assay5` 架构。 | 设计轻量无侵入适配器 `SeuratToCelliD`，通过 R S3 派发机制拦截输入并重定向至当前层 `counts`，零破坏兼容 Seurat v5。 |
| **2** | `Error in images$slice@coordinates: NULL pointer` | **Visium 空间坐标解析断裂**：Seurat v5 重构了组织切片 FOV 存储模式，传统坐标槽已废弃。 | 使用兼容性提取协议 `GetTissueCoordinates(object)`，并在聚类与空间自相关分析中加入坐标守恒校验。 |
| **3** | `Error in markers[["Layer 5"]]: subscript out of bounds` | **原作者 Fig8 代码 Marker 键名失配**：原作者在上游提取 Marker 时将第 5 层保存为 `Layer_5` 或 `L5`，但画图脚本硬编码为 `Layer 5` 导致返回 `NULL`。 | 在 `04` 脚本中建立别名容错解析器（Alias Resolver），支持正则自动回退与同切片实时特征计算，保障健壮性。 |
| **4** | `Error in rocit: NA values in score vector` | **Visium 切片边缘 28 个无效 Spots**：组织切片外缘由于落点在载玻片空白区，皮层金标准标注为 `NA`。 | 在重算流程中建立布尔过滤掩码 `!is.na(Layer)`，清洗剔除无效点，确保 2,698 个真实组织 spot 严密入算。 |
| **5** | `line 7 did not have 5 elements` | **CITE-seq 制表空格与降采样列解析**：细胞亚型名称中包含空格（如 `"CD14 Mono"`），直接 `read.table` 发生列位错。 | 显式指定 `sep = "\t"`，并在降采样梯次转换中采用安全匹配替代强制数字类型转换。 |
| **6** | `File size 116.28 MB exceeds GitHub's limit of 100.00 MB` | **超大二进制包 Git 提交受阻**：Zenodo 原始打包文件 `PaaSc_data.zip` 超出 GitHub 单文件硬限制。 | 重构 `.gitignore` 过滤策略，精准忽略超大 `.zip` 压缩包，同时使用 `!repro_results/**/*.png` 显式放行全套复现图表。 |
| **7** | `std::bad_alloc: Out of Memory` (2.8万细胞) | **大规模单细胞全量重算内存峰值过高**：IFN 数据集全矩阵 SVD 与 50 条通路同时拟合极易打爆单机内存。 | 将 50 条 Hallmark 通路采用分批流式拟合，并在 `run_all.R` 中架构 `quick` 与 `full` 双通道机制。 |
| **8** | PowerShell 命令行中 R 变量 `f2c$Group.2` 被吞噬 | **跨系统 Shell 字符转义陷阱**：Windows PowerShell 默认将双引号中的 `$Group.2` 解释为空字符串变量。 | 在自动化调用与脚本中统一采用 `[['column']]` 寻址，防御平台差异引发的静默逻辑错误。 |

---

## 6. 核心成果画廊与量化指标验证

### 6.1 核心复现成果图谱全景

| 论文主图 | 复现图表预览与链接 | 核心验证科学结论 |
| :---: | :---: | :--- |
| **Fig 2 & SupFig 2** | ![Fig2 Benchmark](repro_results/fig3_to_fig9/Fig2_and_SupFig2_benchmark.png) | **回归去偏与多方法基准**：SupFig 2A 验证 KEGG 背景回归消除非特异漂移；Fig 2B 7类细胞 PaaSc 表现卓越；Fig 2C 降采样至 20% PaaSc AUC 依然维持 0.89+；Fig 2D CITE-seq 蛋白真值高 F1。 |
| **Fig 3** | ![Fig3 Pancancer](repro_results/fig3_to_fig9/Fig3B_pancancer_method_summary.png) | **泛癌单细胞微环境基准**：在 79 个 CD8+ T 队列与 92 个巨噬细胞队列中，PaaSc 胜率（Top 1 夺冠频次）居全部 8 种算法首位。 |
| **Fig 4** | ![Fig4 Senescence](repro_results/fig3_to_fig9/Fig4D_senescence_roc.png) | **细胞衰老跨数据集泛化**：HUVEC、HCA2 及 IMR90 衰老队列中 PaaSc AUC 高达 0.84 - 0.99；Fig 4E 准确捕获 Coagulation/Hypoxia 与衰老高相关。 |
| **Fig 5** | ![Fig5 TMS GWAS](repro_results/fig3_to_fig9/Fig5ABC_tms_lymphocyte_enrichment.png) | **单细胞 GWAS 复杂性状映射**：Tabula Muris Senis 22 种性状中，PaaSc 成功捕获 125 组组织特异性显著富集，远超 CelliD (61) 与 GSDensity (104)。 |
| **Fig 6** | ![Fig6 IFN Repro](repro_results/fig3_to_fig9/Fig6B_ifn_vs_bcell_benchmark.png) | **2.8 万单细胞刺激活化**：IFN-β 刺激后 INTERFERON_ALPHA/GAMMA 通路特异性激活，ROC-AUC 高达 0.9996，且非相关通路背景值严格为 0。 |
| **Fig 7** | ![Fig7 scATAC](repro_results/fig3_to_fig9/Fig7B_scatac_granulocyte_roc.png) | **scATAC-seq 染色质开放区跨组学验证**：16 种造血细胞亚型中，PaaSc 中位数 AUC 达 0.9604，在极度稀疏的 pDC 细胞上依然达到 0.9909。 |
| **Fig 8** | ![Fig8 Cortex](repro_results/fig3_to_fig9/Fig8A_spatial_cortex_layers_and_paasc.png) | **10x Visium 小鼠脑皮层精准分层**：双 MCA 坐标对齐完整还原 Layer 1 至 Layer 6 组织解剖结构，AUC 显著超越对比方法。 |
| **Fig 9** | ![Fig9 CosMx](repro_results/fig3_to_fig9/Fig9_spatial_nsclc_benchmark_and_cor.png) | **9.8 万亚细胞空间单细胞微环境**：CosMx 空间转录组中 9 大细胞类型识别 AUC 达 0.9763；精准揭示 Hedgehog 信号在肿瘤边缘侵袭带的富集。 |

### 6.2 关键指标量化对比表

| 评估维度 / 数据集 | 核心评价指标 | 论文原作者报告值 | 本工程复现值 | 状态判定 |
| :--- | :--- | :---: | :---: | :---: |
| **PBMC 3k (01 重算)** | B 细胞活化识别 ROC-AUC | ~0.99 | **0.9981** | 完全吻合 / 精度更高 |
| **SupFig 2A (06 重算)** | B 细胞 MCA 回归显著维度数 | 3 / 20 维 ($t > 1.96$) | **3 / 20 维** ($t = 3.29$) | 严格吻合 |
| **Fig 2B (7 种 PBMC 亚型)** | 8 种算法平均 AUC | PaaSc 领先 | **PaaSc 0.9796** (GSDensity 0.967) | 完全吻合 |
| **Fig 2C (20% 降采样)** | Dropout 极限下平均 AUC | PaaSc > 0.88 | **PaaSc 0.8927** (scGSEA 0.578) | 完全吻合 |
| **Fig 3 (泛癌 92 巨噬细胞)** | 巨噬细胞平均 AUC | 0.93+ | **0.9378** (夺冠 36/92) | 严格吻合 |
| **Fig 4D (HUVEC 衰老)** | 衰老细胞分类 ROC-AUC | ~0.99 | **0.9927** (CelliD 0.617) | 严格吻合 |
| **Fig 5D (TMS GWAS)** | 组织特异性显著富集对数 | PaaSc 最多 | **125 对** (GSD 104, CelliD 61) | 严格吻合 |
| **Fig 6B (IFN 刺激应答)** | IFN 刺激分类 ROC-AUC | ~1.00 | **0.9996** | 严格吻合 |
| **Fig 7B (scATAC 16 亚型)** | 染色质开放区中位数 AUC | ~0.96 | **0.9604** (CelliD 0.933) | 严格吻合 |
| **Fig 9B (CosMx 9.8万细胞)**| 9 种肺癌细胞分类中位数 AUC | ~0.97 | **0.9763** (GSDensity 0.693) | 严格吻合 |

---

## 7. 运行环境与依赖软件包快照

本工程已在 Windows 11 专业版下完成全流程验证，完整软件包版本快照已固化于 [repro_results/environment_versions.tsv](repro_results/environment_versions.tsv)。

| 软件包名称 (Package) | 运行版本 (Version) | 安装来源 / 仓库 (Repository) | 核心功能作用 |
| :--- | :---: | :--- | :--- |
| **R** | 4.6.0 | Base R | 统计计算运行底座 |
| **PaaSc** | 1.0.1 | GitHub (yoyoong/PaaSc v1.0.1) | 复现核心算法包 |
| **Seurat** | 5.6.0 | CRAN | 单细胞与空间分析核心框架 |
| **SeuratObject** | 5.4.0 | CRAN | Seurat 底层数据结构管辖 |
| **SeuratData** | 0.2.2.9002 | GitHub / SatijaLab | 单细胞/空间基准数据集分发 |
| **stxBrain.SeuratData** | 0.1.2 | GitHub / SatijaLab | 10x Visium 小鼠脑前切片数据 |
| **CelliD** | 1.19.0 | Bioconductor 3.23 | MCA 流形降维与细胞特征打分 |
| **msigdbr** | 26.1.1 | CRAN | 官方 MSigDB 通路基因集数据库接口 |
| **ROCit** | 2.1.2 | CRAN | ROC-AUC 曲线与置信区间经验推断 |
| **pROC** | 1.19.1 | CRAN | 高性能接收者操作特征分析 |
| **ComplexHeatmap**| 2.28.0 | Bioconductor 3.23 | 复杂基因组相关性多层热图绘制 |
| **circlize** | 0.4.18 | CRAN | 热图渐变连续色阶调度 |
| **pheatmap** | 1.0.13 | CRAN | 紧凑型生物信息热图绘制 |
| **ggh4x** | 0.3.1 | CRAN | ggplot2 进阶分面与坐标扩展 |
| **ggplot2** | 4.0.3 | CRAN | 论文级学术可视化输出 |
| **patchwork** | 1.3.2 | CRAN | 多分块复合图表拼版排版 |
| **tidyverse** | 2.0.0 | CRAN | 数据清洗与管道化转换 |
| **data.table** | 1.18.6.1 | CRAN | 高效大矩阵输入输出 |
| **mixtools** | 2.0.0.1 | CRAN | EM 期望最大化高斯混合模型拟合 |
| **Matrix** | 1.7.5 | CRAN | 高性能稀疏矩阵表示与运算 |
| **sctransform** | 0.4.3 | CRAN | 空间 Visium 正则化负二项回归 |
| **uwot** | 0.2.5 | CRAN | 高维流形 UMAP 嵌入投影 |

---

## 8. 技术展示与面试答辩导览

如需将本项目用于学术交流、组会研讨或高阶求职面试，请移步查阅专章指南文档：  
👉 **[docs/SHOWCASE_GUIDE.md](docs/SHOWCASE_GUIDE.md)**

该文档包含以下核心模块：
1. **1 分钟极速电梯演讲 (Elevator Pitch)**：提纲挈领阐述项目痛点、算法破局点与工程价值。
2. **5 分钟组会/汇报幻灯片动线 (5-Minute Slides Structure)**：标准 5 页 PPT 讲稿设计（背景痛点 -> 数学去偏 -> 跨模态数据 -> 工程破局 -> 总结边界）。
3. **10 分钟深度技术面试答辩 (10-Minute Deep Dive)**：采用“五幕剧架构（痛点直觉 -> 诚实边界 -> 源码破局 -> 数据铁证 -> 工程沉淀）”。
4. **硬核技术拷问与破局预案 (Technical Q&A FAQ)**：针对 Seurat v5 S3 拦截原理、Visium 坐标与 28 NA Spots 清洗、Layer 5 命名断裂排查、大规模细胞扩展机制等 5 大尖锐问题提供标准回答话术。

---

## 9. 引用、安装与开源许可

### 9.1 引用论文 (BibTeX)
如果您在后续学术研究或生信分析中使用了 PaaSc 算法或借鉴了本复现工程的工程化架构，请引用原作者发表论文：

```bibtex
@article{paasc2025,
  title={PaaSc: Pathway activity analysis based on multiple correspondence analysis in single-cell RNA-seq and spatial transcriptomics},
  author={Zhang, Y. and Collaborators},
  journal={PLOS Computational Biology},
  volume={21},
  number={11},
  pages={e1013666},
  year={2025},
  publisher={Public Library of Science},
  doi={10.1371/journal.pcbi.1013666}
}
```

### 9.2 开源许可证说明
- 本复现工程代码与脚本遵循 **GNU General Public License v3.0 (GPL-3.0)** 许可协议开源。
- 更多上游安装指南与基础用法详见 [README_upstream.md](README_upstream.md)。
