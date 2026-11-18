# Zenodo 归档数据字典 (Zenodo Data Dictionary)

本文档详细记录 `zenodo_data/PaaSc_data/` 目录下全部上游作者归档数据的来源、文件结构、生物学背景及对应的复现脚本映射。

---

## 1. 数据集总体概况

上游作者通过 Zenodo 开源了支撑论文各章节核心图表的数据包（压缩包 `PaaSc_data.zip` 约 116.3 MB，已在 `.gitignore` 中忽略，解压后目录为 `zenodo_data/PaaSc_data/`）。

| 目录/图表 | 数据形态 | 样本规模/维度 | 生物学系统 | 对应复现脚本 |
| :--- | :--- | :--- | :--- | :--- |
| **Fig2** | 汇总制表 (`.txt`) | 8种算法 × 7种亚型 / 降采样梯次 | 10x Genomics PBMC 3k & CITE-seq | `06_fig2_and_supfig2_repro.R` |
| **Fig3** | 队列AUC制表 (`.txt`) | 5种主要免疫细胞 × 79-92个癌症队列 | 泛癌单细胞肿瘤微环境 (TME) | `05_fig3_fig4_fig5_fig7_fig9_benchmarks.R` |
| **Fig4** | 评分矩阵与相关性 (`.tsv`) | 4个细胞系数据集 / 50条Hallmark | 细胞衰老 (Senescence) 与泛癌相关性 | `05_fig3_fig4_fig5_fig7_fig9_benchmarks.R` |
| **Fig5** | 单细胞GWAS通路评分 (`.tsv`) | 10,725 cells × 22 GWAS traits | 小鼠衰老图谱 (Tabula Muris Senis) | `05_fig3_fig4_fig5_fig7_fig9_benchmarks.R` |
| **Fig6** | Seurat对象 (`.RDS`) + 打分表 | 28,871 cells × 35,635 genes | 人外周血单核细胞 IFN-β 刺激响应 | `03_fig6_ifn_repro.R` |
| **Fig7** | RData 对象 (`.RData`) | 16种造血分化细胞亚型 | 骨髓 scATAC-seq 染色质开放区 | `05_fig3_fig4_fig5_fig7_fig9_benchmarks.R` |
| **Fig8** | Visium RDS + 空间评分 | 2,698 spots × 6皮层区域 | 小鼠脑前脑冠状切片 (10x Visium) | `04_fig8_st_cortex_repro.R` & `05` |
| **Fig9** | 单细胞空间转录组 (`.csv`/`.txt`)| 98,002 cells × 960 基因 | 非小细胞肺癌 (NanoString CosMx SMI) | `05_fig3_fig4_fig5_fig7_fig9_benchmarks.R` |
| **SupFigs** | 辅助元数据与补充表格 | 元数据、伪代码、B细胞亚群AUC | TMS 元数据、算法伪代码及补充材料 | `05`, `06` |

---

## 2. 逐图数据文件详解

### 2.1 Fig2: PBMC基准、降采样鲁棒性与CITE-seq蛋白真值
- **`Fig2/2B.txt`** (418 B)
  - **内容**: 8种通路活性评估算法（`PaaSc`, `GSDensity`, `CellID`, `AUCell`, `VAM`, `scGSEA`, `ssGSEA`, `GSVA`）在 PBMC 7 种主要免疫亚型（B, CD4 T, CD8 T, CD14 Mono, FCGR3A Mono, NK, DC）上的分类 ROC-AUC。
  - **复现使用**: `06_fig2_and_supfig2_repro.R` 绘制主图 Fig 2B 多方法箱线图。
- **`Fig2/2C.txt`** (2.4 KB)
  - **内容**: 单细胞测序深度梯度降采样（20%, 30%, 40%, 50%, 60%, 70%, 80%, 90%, 100%/full.marker.set）下各算法的 AUC 衰减曲线。
  - **复现使用**: `06_fig2_and_supfig2_repro.R` 绘制主图 Fig 2C Dropout 鲁棒性折线图。
- **`Fig2/2D.txt`** (18.0 KB)
  - **内容**: 10x Genomics CITE-seq 表面蛋白作为细胞真值标志，评估转录组通路打分的 F1-score、Precision 和 Recall。
  - **复现使用**: `06_fig2_and_supfig2_repro.R` 绘制主图 Fig 2D 蛋白真值基准对比。

### 2.2 Fig3: 泛癌免疫细胞表征多中心基准
- **`Fig3/CD8T_auc_values.txt`** (12.6 KB): 79 个不同癌种/临床队列中 CD8+ T 细胞的 8 种算法 AUC 分值。
- **`Fig3/Macro_auc_values.txt`** (14.6 KB): 92 个癌症队列中巨噬细胞 (Macrophage) 的算法 AUC 分值。
- **`Fig3/NK_auc_values.txt`**, **`DC_auc_values.txt`**, **`Treg_auc_values.txt`**: 对应自然杀伤细胞 (27队列)、树突状细胞 (16队列) 和调节性 T 细胞 (11队列) 的 AUC 评测。
- **`Fig3/dataset_info.txt`** (494.6 KB) & **`Celltype.txt`**: 样本元数据（GEO 号、癌种、细胞数）。
- **复现使用**: `05_fig3_fig4_fig5_fig7_fig9_benchmarks.R` 汇总重构 5 大细胞类型泛癌热图与胜率汇总柱状图 (Fig 3A & 3B)。

### 2.3 Fig4: 细胞衰老标志物与泛癌通路关联
- **`Fig4/GSE102090_Score.tsv`** (796.1 KB): HUVEC 内皮细胞衰老模型打分 (PaaSc, CelliD, GSDensity)。
- **`Fig4/GSE119807_Score.tsv`** (55.7 KB): HCA2 人皮肤成纤维细胞衰老打分。
- **`Fig4/GSE115301_10X_Score.tsv`** (36.3 KB) & **`GSE115301_Timecourse_Score.tsv`** (17.5 KB): IMR90 肺成纤维细胞时间序列衰老打分。
- **`Fig4/all_pid_celltype_hallmarkVSsene_correlation.tsv`** (6.5 MB): 泛癌单细胞中 50 条 MSigDB Hallmark 与衰老特征的相关性矩阵。
- **复现使用**: `05_fig3_fig4_fig5_fig7_fig9_benchmarks.R` 绘制 4 个独立队列的衰老 ROC 曲线 (Fig 4D) 及泛癌 Hallmark-衰老相关性气泡/热图 (Fig 4E)。

### 2.4 Fig5: 单细胞 GWAS 复杂性状富集 (Tabula Muris Senis)
- **`Fig5/TMS_facs_selected_alltraits_PaaSc_score.tsv`** (1.5 MB): 10,725 个小鼠单细胞在 22 个人类复杂表型/疾病 GWAS 基因集上的 PaaSc 活性评分。
- **`Fig5/TMS_facs_selected_alltraits_CelliD_score.tsv`**, **`GSDensity_score.tsv`**: 对应对比算法评分。
- **复现使用**: `05_fig3_fig4_fig5_fig7_fig9_benchmarks.R` 绘制淋巴细胞计数性状的 GMM 双峰二值化 (Fig 5B) 与组织特异性显著富集热图 (Fig 5D)。

### 2.5 Fig6: 单细胞 IFN-β 刺激免疫应答 (Kang et al. 2018)
- **`Fig6/IFN.RDS`** (38.4 MB): 完整 Seurat 对象，包含 28,871 个外周血免疫单细胞（分为 CONTROL 与 STIM 组），含 8 类免疫亚型。
- **`Fig6/PaaSc_scores.txt`** (22.4 MB), **`CelliD_scores.txt`**, **`GSdensity_scores.txt`**: 预存的 50 条 MSigDB Hallmark 评分表。
- **复现使用**: `03_fig6_ifn_repro.R` 支持直接载入 `IFN.RDS` 进行**端到端全量重算** (耗时约 12 分钟)，验证 INTERFERON_ALPHA/GAMMA_RESPONSE 活性的特异性激活与 ROC-AUC（PaaSc AUC 达 0.9996）。

### 2.6 Fig7: 单细胞染色质可及性 (scATAC-seq 造血分化)
- **`Fig7/granulocyte_ROC.RData`** (17.6 MB): 包含 16 种不同造血前体及分化阶段免疫细胞的染色质开放区 Peak 矩阵与通路 AUC 曲线数据。
- **复现使用**: `05_fig3_fig4_fig5_fig7_fig9_benchmarks.R` 验证 PaaSc 在表观遗传稀疏单细胞数据上的跨组学适用性 (Fig 7B)。

### 2.7 Fig8: 空间转录组脑皮层层级特异性 (10x Visium)
- **`Fig8/ST_braincortex.RDS`** (12.6 MB): 10x Genomics 官方小鼠大脑冠状切片 Visium 数据 (2,698 spots)。
- **`Fig8/ST_Braincortex_top200_markers.RDS`** (16 KB): 大脑皮层各层 (Layer 1 至 Layer 6 及白质) 前 200 个特征 Marker 基因集。
- **`Fig8/ST_braincortex_PaaSc.csv`**, **`CelliD.csv`**, **`GSDensity.csv`**: 空间点通路活性矩阵。
- **`Fig8/braincortex_ROC.RData`** (303 KB): 空间分层识别 ROC 曲线数据。
- **复现使用**:
  - `04_fig8_st_cortex_repro.R` 从原始 counts 矩阵进行端到端重算（含双 MCA 空间分析、28个空斑点清洗与皮层特征映射）。
  - `05_fig3_fig4_fig5_fig7_fig9_benchmarks.R` 绘制全方法脑皮层 AUC 对比箱线图 (Fig 8B)。

### 2.8 Fig9: 亚细胞空间转录组微环境互作 (NanoString CosMx SMI)
- **`Fig9/spatial_lung_AUC.txt`** (624 B): 98,002 个细胞在 9 种非小细胞肺癌 (NSCLC) 细胞类型上的 AUC 评估指标。
- **`Fig9/spatial_lung_PaaSc.txt`** (26.6 MB), **`CelliD.txt`** (82.0 MB): 9.8 万单细胞空间分辨率通路活性矩阵。
- **`Fig9/ST_lung_hallmark_paasc.csv`** (90.8 MB), **`ST_lung_hallmark_SMI.csv`** (2.6 KB): 空间临近区与肿瘤侵袭前沿的 Hallmark 通路共表达矩阵。
- **复现使用**: `05_fig3_fig4_fig5_fig7_fig9_benchmarks.R` 绘制 9.8 万细胞基准条形图与肿瘤微环境 Hallmark 空间互作热图 (Fig 9B, 9D, 9E)。

---

## 3. 存储与分发说明

由于 GitHub 仓库存在 **100 MB 单文件传输限制**：
- 原始压缩包 `zenodo_data/PaaSc_data.zip` (116.28 MB) 已被 `.gitignore` 规则屏蔽。
- 解压后的所有中间打分矩阵、RDS 对象及制表文件均严格控制在 100 MB 以内，可随 Git 仓库顺畅协同。
- 如需重新下载最原始压缩包，请参阅上游 Zenodo 归档页面或运行数据下载指引。
