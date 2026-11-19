# Zenodo 归档数据字典 (Zenodo Data Dictionary)

本文档经源码与磁盘文件逐一核验，记录 `zenodo_data/PaaSc_data/` 目录下全部上游作者归档数据的精确文件体积、维度、生物学背景及复现脚本映射（归档来源：Zenodo `10.5281/zenodo.17138447` / `17136774`）。

---

## 1. 数据集总体概况

上游作者通过 Zenodo 开源了支撑论文各章节核心图表的数据包（原始压缩包 `PaaSc_data.zip` 为 `116.28 MB`，解压至 `zenodo_data/PaaSc_data/`）。

| 目录/图表 | 数据形态 | 核验后样本规模/维度 | 生物学系统 | 对应复现脚本 |
| :--- | :--- | :--- | :--- | :--- |
| **Fig2** | 汇总制表 (`.txt`) | 8 种算法 × 7 种 PBMC 亚型 / 标记基因集降采样 / CITE-seq F1 | 10x Genomics PBMC & CITE-seq | `06_fig2_and_supfig2_repro.R` |
| **Fig3** | 队列 AUC 制表 (`.txt`) | 5 种免疫谱系 × 共 `225` 组队列-谱系（79 CD8 T / 92 巨噬 / 27 NK / 16 DC / 11 Treg） | TISCH 泛癌单细胞肿瘤微环境 (TME) | `05_fig3_fig4_fig5_fig7_fig9_benchmarks.R` |
| **Fig4** | 评分矩阵与相关性 (`.tsv`) | 3 套主图衰老队列（归档含 4 文件） / 50 条 Hallmark | 细胞衰老 (Senescence) 与泛癌协同 | `05_fig3_fig4_fig5_fig7_fig9_benchmarks.R` |
| **Fig5** | 单细胞 GWAS 通路评分 (`.tsv`) | `10,725` cells × **`75` UKB/GWAS traits**（主图 Fig 5D 评测 `22` 个组织特异代表性状） | 小鼠衰老图谱 (Tabula Muris Senis FACS) | `05_fig3_fig4_fig5_fig7_fig9_benchmarks.R` |
| **Fig6** | Seurat 对象 (`.RDS`) + 打分表 | **`16,142` genes × `28,871` cells**（`ctrl` 14,526 / `stim` 14,345） | 人 PBMC IFN-$\beta$ 刺激响应 (GSE96583) | `03_fig6_ifn_repro.R` |
| **Fig7** | RData 对象 (`.RData`) | `10,459` cells / 16 种 PBMC 免疫亚型 | **PBMC scATAC-seq**（基因活性分数与 ROC 对象） | `05_fig3_fig4_fig5_fig7_fig9_benchmarks.R` |
| **Fig8** | Visium RDS + 空间评分 | **`4,000` 基因 × `3,639` spots**（`3,611` 有效注释 + `28` `NA`） | **人脑前额叶皮层 DLPFC**（10x Visium Sample 151673） | **`04_fig8_st_cortex_repro.R`**（同时绘 Fig 8A & 8B） |
| **Fig9** | 单细胞空间转录组 (`.csv`/`.txt`) | `98,002` cells（9 类 TME 细胞 + 50 Hallmark） | 非小细胞肺癌 (NanoString CosMx SMI, RNA) | `05_fig3_fig4_fig5_fig7_fig9_benchmarks.R` |
| **SupFigs** | 辅助元数据与补充表格 | TMS 元数据、算法伪代码及补充表格 | 辅助对照与元数据支持 | `05`, `06` |

---

## 2. 逐图数据文件详解

### 2.1 Fig2: PBMC 基准、标记基因集降采样鲁棒性与 CITE-seq 蛋白真值
- **`Fig2/2B.txt`** (`418 B`)
  - **内容**：8 种算法（`PaaSc`, `GSDensity`, `CelliD`, `AUCell`, `VAM`, `scGSEA`, `GSEA`, `GSVA`）在 7 种 PBMC 亚型（`Bcell`, `CD4`, `CD8`, `CD14`, `CD16`, `DC`, `pDC`）上的 ROC-AUC。
  - **复现使用**：`06_fig2_and_supfig2_repro.R` 绘制主图 Fig 2B。
- **`Fig2/2C.txt`** (`2.4 KB`)
  - **内容**：**标记基因集梯度降采样**（`20%` 至 `100%` / `full.marker.set`）下各算法的平均 AUC 稳健性曲线。
  - **复现使用**：`06_fig2_and_supfig2_repro.R` 绘制主图 Fig 2C。
- **`Fig2/2D.txt`** (`18.0 KB`)
  - **内容**：CITE-seq 表面蛋白作为真值基准的 F1-score、Precision 与 Recall 制表（Tab 分隔，亚型名含空格）。
  - **复现使用**：`06_fig2_and_supfig2_repro.R` 绘制主图 Fig 2D。

### 2.2 Fig3: TISCH 泛癌免疫细胞多队列基准
- **`Fig3/CD8T_auc_values.txt`** (`12.6 KB`, 79 队列)、**`Macro_auc_values.txt`** (`14.6 KB`, 92 队列)、**`NK_auc_values.txt`** (27 队列)、**`DC_auc_values.txt`** (16 队列)、**`Treg_auc_values.txt`** (11 队列)：合计 `225` 组队列-谱系评估。
- **`Fig3/dataset_info.txt`** (`494.6 KB`) & **`Celltype.txt`**：队列元数据与细胞类型信息。
- **复现使用**：`05_fig3_fig4_fig5_fig7_fig9_benchmarks.R` 绘制 Fig 3A 热图与 Fig 3B Rank-1 胜出汇总柱状图。

### 2.3 Fig4: 细胞衰老标志物与泛癌通路协同
- **`Fig4/GSE102090_Score.tsv`** (`796.1 KB`, HUVEC)、**`GSE119807_Score.tsv`** (`55.7 KB`, HCA2)、**`GSE115301_10X_Score.tsv`** (`36.3 KB`, IMR90)：`05` 绘制这 **3 套独立队列**的 Fig 4D ROC 曲线（归档中另含 `GSE115301_Timecourse_Score.tsv` `17.5 KB` 未纳入 Fig 4D）。
- **`Fig4/all_pid_celltype_hallmarkVSsene_correlation.tsv`** (`6.5 MB`)：泛癌单细胞 50 条 Hallmark 与衰老特征相关性矩阵，用于绘制 Fig 4E。

### 2.4 Fig5: 单细胞 GWAS 复杂性状富集 (Tabula Muris Senis)
- **`Fig5/TMS_facs_selected_alltraits_PaaSc_score.tsv`** (**`15.08 MB`**)：`10,725` 个小鼠单细胞在 **`75` 个 UKB/GWAS 性状**上的 PaaSc 打分矩阵。
- **`Fig5/TMS_facs_selected_alltraits_CelliD_score.tsv`**、**`GSDensity_score.tsv`**：对比算法打分矩阵。
- **复现使用**：`05_fig3_fig4_fig5_fig7_fig9_benchmarks.R` 绘制 `LYMPHOCYTE_COUNT` GMM 二值化与 Fisher 富集（Fig 5A–C）及 `22` 个组织特异代表性状跨算法对比（Fig 5D）。

### 2.5 Fig6: 人 PBMC IFN-$\beta$ 刺激响应 (GSE96583)
- **`Fig6/IFN.RDS`** (`38.4 MB`，`.gitignore` 屏蔽)：Seurat 对象，包含 **`16,142` 基因 × `28,871` 细胞**（`ctrl` 14,526 / `stim` 14,345，其中 B 细胞 2,871 个）。
- **`Fig6/PaaSc_scores.txt`** (**`2.24 MB`**)、**`CelliD_scores.txt`**、**`GSdensity_scores.txt`**：作者归档的 **2 条通路**（`B_cell_activation` 与 `INTERFERON`）细胞级评分表。
- **复现使用**：`03_fig6_ifn_repro.R` 基于 `IFN.RDS` 端到端重算 `PaaSc` 与 `CelliD`（`GOBP_INTERFERON_MEDIATED_SIGNALING_PATHWAY` + `B_cell_activation`），实测 IFN 刺激 `PaaSc AUC = 0.9573`（与归档 `r = 1.000000`），阴性对照 B 细胞活化 `AUC = 0.4635`。

### 2.6 Fig7: PBMC 单细胞染色质可及性 (scATAC-seq)
- **`Fig7/granulocyte_ROC.RData`** (**`1.76 MB`**)：包含 **PBMC scATAC-seq（`10,459` 细胞）** 16 种免疫亚型基因活性评分及 ROC 评估结果。
- **复现使用**：`05_fig3_fig4_fig5_fig7_fig9_benchmarks.R` 绘制 Fig 7B 跨算法箱线图。

### 2.7 Fig8: 人脑皮层 DLPFC 空间转录组层级特异性 (10x Visium)
- **`Fig8/ST_braincortex.RDS`** (`12.6 MB`，`.gitignore` 屏蔽)：**人脑前额叶背外侧皮层 DLPFC（Sample 151673）** 10x Visium Seurat 对象，含 **`4,000` 基因 × `3,639` spots**（有效皮层注释 `3,611` spots + `28` 个边缘 `NA` spots）。
- **`Fig8/ST_Braincortex_top200_markers.RDS`** (`16 KB`)：19 个皮层与细胞类型 Top 200 标记基因集。
- **`Fig8/ST_braincortex_PaaSc.csv`**、**`CelliD.csv`**、**`GSDensity.csv`** 与 **`braincortex_ROC.RData`** (`303 KB`)：作者归档空间得分与 ROC 对象。
- **复现使用**：**`04_fig8_st_cortex_repro.R`** 单独完成端到端双 MCA 重算（过滤 `28` 个 `NA` spots）并同时输出 **Fig 8A 与 Fig 8B**。

### 2.8 Fig9: 非小细胞肺癌亚细胞空间转录组 (NanoString CosMx SMI)
- **`Fig9/spatial_lung_AUC.txt`** (`624 B`)：98,002 细胞在 9 类 NSCLC 细胞类型上的多方法 AUC 汇总。
- **`Fig9/spatial_lung_PaaSc.txt`** (`26.6 MB`)、**`spatial_lung_CelliD.txt`** (**`8.20 MB`**)：98,002 细胞空间分辨率细胞类型打分矩阵。
- **`Fig9/ST_lung_hallmark_paasc.csv`** (`90.8 MB`)、**`ST_lung_hallmark_SMI.csv`** (`2.6 KB`)：50 条 Hallmark 细胞级打分矩阵与空间互信息（Spatial MI）统计表。
- **复现使用**：`05_fig3_fig4_fig5_fig7_fig9_benchmarks.R` 绘制 Fig 9B（9 亚型 AUC）、Fig 9D（Top 20 空间互信息 Hallmark）与 Fig 9E（98,002 细胞向量化 Pearson 相关性热图）。

---

## 3. Git 跟踪与本地复现存储边界

为保持 Git 仓库轻量并遵守 GitHub 单文件 `< 100 MB` 限制，`.gitignore` 采用分级屏蔽策略：
- **已由 `.gitignore` 屏蔽的大体积文件**：`zenodo_data/PaaSc_data.zip`（`116.28 MB`）、`zenodo_data/PaaSc_data/Fig6/IFN.RDS`（`38.4 MB`）及 `zenodo_data/PaaSc_data/Fig8/ST_braincortex.RDS`（`12.6 MB`）。
- **已随仓保留的轻量基准表**：`Fig2`–`Fig5`、`Fig7`、`Fig9` 的文本/RData 基准表（除上述两个大 RDS 外）以及全部 `repro_results/` 图表与定量验收表。
- **全新克隆运行指南**：直接克隆本仓库可立即运行 `01`、`02`、`05`、`06`；若需运行 `03`（Fig 6 端到端重算）或 `04`（Fig 8 端到端重算）及 `run_all.R`，请先从 Zenodo 下载 `PaaSc_data.zip` 并将 `Fig6/IFN.RDS` 与 `Fig8/ST_braincortex.RDS` 放回对应路径。
