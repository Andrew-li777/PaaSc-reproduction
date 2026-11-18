
library('PaaSc')
library('Seurat')
library('Matrix')
library('dplyr')
library('gsdensity')
library('CelliD')
library("ROCit")
library('tidyr')
library('ggplot2')

# loading mtx files
    # matrix_data <- readMM("data/IFN/matrix.mtx")
    # obs = read.table("data/IFN/GSE96583_batch2.total.tsne.df.tsv", row.names = 1, header = TRUE, sep="\t")
    # genef = read.table("data/IFN/GSE96583_batch2.genes.tsv", header=FALSE, row.names = NULL, sep="\t")
    # dimnames(matrix_data) <- list(make.unique(genef$V2), rownames(obs))

    # IFN <- CreateSeuratObject(counts = matrix_data, project = "IFN", min.cells = 3, min.features = 200)
    # tmp <- obs[Cells(IFN), ]
    # IFN[[]] <- data.frame(IFN[[]], tmp)
    # saveRDS(IFN, file = 'data/IFN.RDS')

IFN <- readRDS(file = 'data/IFN.RDS')

# Seurat pipeline
IFN <- IFN %>% NormalizeData %>% FindVariableFeatures(nfeatures = 2000) %>% 
    ScaleData %>% RunPCA %>% RunUMAP(dims = 1:25, seed.use=2)

pdf('out/IFN_PCA.pdf', width = 6, height = 5)
    ElbowPlot(IFN, ndims = 40, reduction = "pca")
dev.off()

col_palatte <- c('#6965CE', '#EE5755', '#75E15A', '#D5CE85', '#9F33D9', '#9BE2CF', '#69A5C8', 
    '#FC9500', '#D846AF', '#CE95C4', '#0018F9')

pdf('out/IFN_DimPlot_cellType.pdf', width = 6, height = 5)
    DimPlot(IFN, group.by = 'cell', label = F, pt.size = 0.01, cols = col_palatte)
dev.off()

pdf('out/IFN_DimPlot_batch.pdf', width = 6, height = 5)
    DimPlot(IFN, group.by = 'stim', label = F, pt.size = 0.01, cols = col_palatte)
dev.off()

# IFN and B-cell activation pathways
object <- computeMCA(IFN, nmcs = 20)
background_geneset <- "data/kegg_legacy.gmt"
INTERFERON <- strsplit(readLines("data/GOBP_INTERFERON_MEDIATED_SIGNALING_PATHWAY.v2025.1.Hs.gmt"), '\t')[[1]]
INTERFERON <- INTERFERON[3:length(INTERFERON)]

B_cell_activation <- readLines("data/B_cell_activation.txt")
pathway_geneset <- list(B_cell_activation = B_cell_activation, INTERFERON = INTERFERON)
gene_rate <- getGeneRate(background.geneset = background_geneset, pathway.geneset = pathway_geneset)
regression_data <- doRegression(object, gene.rate = gene_rate)
score_data <- computeScore(object, regression.data = regression_data, weight = FALSE, normalize = "z-score")
out <- data.frame(score_data, cellType = IFN$cell, condition = IFN$stim)
write.table(out, file = 'out/PaaSc_scores.txt', row.names = TRUE, col.names = TRUE, sep = '\t')

# GSdensity
IFN <- readRDS(file = 'data/IFN.RDS')
cells <- colnames(IFN)
ce <- compute.mca(object = IFN)
el <- compute.nn.edges(coembed = ce, nn.use = 300)
cv.df <- run.rwr.list(el = el, gene_set_list = pathway_geneset, cells = cells)
out <- data.frame(cv.df, cellType = IFN$cell, condition = IFN$stim)
write.table(out, file = 'out/GSdensity_scores.txt', row.names = TRUE, col.names = TRUE, sep = '\t')

# MCA
IFN <- readRDS(file = 'data/IFN.RDS')
IFN <- NormalizeData(IFN)
IFN <- RunMCA(IFN)
HGT <- RunCellHGT(IFN, pathways = pathway_geneset, dims = 1:50, n.features = 200)
score_data <- data.frame(t(HGT))
out <- data.frame(score_data, cellType = IFN$cell, condition = IFN$stim)
write.table(out, file = 'out/CelliD_scores.txt', row.names = TRUE, col.names = TRUE, sep = '\t')

# summary
o1 <- read.table(file = 'out/PaaSc_scores.txt', row.names = 1, header = TRUE, sep = '\t')
o2 <- read.table(file = 'out/GSdensity_scores.txt', row.names = 1, header = TRUE, sep = '\t')
o3 <- read.table(file = 'out/CelliD_scores.txt', row.names = 1, header = TRUE, sep = '\t')

roc_summary <- function(all_label, ref_label, score){
    roc_empirical <- rocit(score = score, class = all_label, negref = ref_label)
    tmp <- strsplit(as.character(summary(roc_empirical)[4,]), ':')[[1]][2]
    tmp <- as.numeric(gsub(' ', '', tmp))
    return(tmp)
}

# 比较IFN信号
data <- data.frame(condition = IFN$stim, 
                   PaaSc = rank(o1$INTERFERON)/nrow(o1),
                   GSdensity = rank(o2$INTERFERON)/nrow(o2),
                   CelliD = rank(o3$INTERFERON)/nrow(o3))
IFN_roc <- c(PaaSc = roc_summary(all_label = data$condition, ref_label = 'ctrl', score = data$PaaSc),
             GSdensity = roc_summary(all_label = data$condition, ref_label = 'ctrl', score = data$GSdensity),
             CelliD = roc_summary(all_label = data$condition, ref_label = 'ctrl', score = data$CelliD))

data_long <- data %>%
  pivot_longer(cols = c(PaaSc, GSdensity, CelliD), 
               names_to = "variable", 
               values_to = "value")

p <- ggplot(data_long, aes(x = variable, y = value, fill = condition)) +
  geom_boxplot(alpha = 0.7, outlier.alpha = 0.6) +
  scale_fill_manual(values = c("ctrl" = "grey60", "stim" = "red")) +
  labs(x = "Variables", 
       y = "Rank normalized scores", 
       fill = "Condition",
       title = "IFN pathway activities") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(hjust = 0.5))

pdf('out/Score_IFN.pdf', width = 6, height = 4)
    print(p)
dev.off()

# B-cell activation
data <- data.frame(condition = IFN$stim, 
                   PaaSc = rank(o1$B_cell_activation)/nrow(o1),
                   GSdensity = rank(o2$B_cell_activation)/nrow(o2),
                   CelliD = rank(o3$B_cell_activation)/nrow(o3))
data <- data[IFN$cell == 'B cells', ]

B_roc <- c(PaaSc = roc_summary(all_label = data$condition, ref_label = 'ctrl', score = data$PaaSc),
             GSdensity = roc_summary(all_label = data$condition, ref_label = 'ctrl', score = data$GSdensity),
             CelliD = roc_summary(all_label = data$condition, ref_label = 'ctrl', score = data$CelliD))

data_long <- data %>%
  pivot_longer(cols = c(PaaSc, GSdensity, CelliD), 
               names_to = "variable", 
               values_to = "value")

p <- ggplot(data_long, aes(x = variable, y = value, fill = condition)) +
  geom_boxplot(alpha = 0.7, outlier.alpha = 0.6) +
  scale_fill_manual(values = c("ctrl" = "grey60", "stim" = "red")) +
  labs(x = "Variables", 
       y = "Rank normalized scores", 
       fill = "Condition",
       title = "B-cell activation in B cells") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(hjust = 0.5))

pdf('out/Score_Bcells.pdf', width = 6, height = 4)
    print(p)
dev.off()

# ROC
> IFN_roc
    PaaSc GSdensity    CelliD 
   0.9573    0.5575    0.6037 

> B_roc
    PaaSc GSdensity    CelliD 
   0.4635    0.5562    0.5122

