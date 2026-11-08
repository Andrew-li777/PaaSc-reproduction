# source /sibcb/program/install/r-4.4/profile
# source /sibcb/program/install/hdf5-1.10/profile

library('PaaSc')
library('gsdensity')
library('Signac')
library('GenomicRanges')
library('Seurat')
library('dplyr')
library('CelliD')
library('Matrix')
library('rGREAT')
library('ROCit')
library('ggplot2')
library('AnnotationHub')
library('EnsDb.Hsapiens.v86')
library('SeuratDisk')
library('tidyr')
library('VennDiagram')
library('BSgenome.Hsapiens.UCSC.hg38')

load('scATACseq_granulocyte_GA.Rata')

roc_summary <- function(all_label, ref_label, score){
    roc_empirical <- rocit(score = score, class = all_label, negref = ref_label)
    tmp <- strsplit(as.character(summary(roc_empirical)[4,]), ':')[[1]][2]
    tmp <- as.numeric(gsub(' ', '', tmp))
    return(tmp)
}

# find markers
reference <- LoadH5Seurat("pbmc_multimodal.h5seurat", assays = list("SCT" = "counts"), reductions = 'spca')
reference <- UpdateSeuratObject(reference)

newtype <- sapply(strsplit(as.character(Idents(reference)), '_'), function(z) z[1])
newtype <- sub("^(\\S+(?:\\s+\\S+)?).*", "\\1", newtype)
newtype[newtype == 'Treg Memory'] = 'Treg'
newtype[newtype == 'Treg Naive'] = 'Treg'
Idents(reference) = newtype

ktype <- levels(unique(Idents(pbmc)))
marker_list <- list()

for(i in 1:length(ktype)){
	tag <- ktype[i]
	cat('Processing', tag, '\n')
	tmp <- FindMarkers(reference, ident.1 = tag, only.pos = TRUE, test.use = "wilcox") %>% dplyr::filter(p_val_adj < 0.05)
	marker_list[[tag]] <- rownames(tmp)[1:min(200, nrow(tmp))]
}
names(marker_list) <- gsub(' ', '_', names(marker_list))

# PaaSc
GA <- CreateSeuratObject(counts = gene.activities, project = "GA", min.cells = 3, min.features = 200)
GA <- computeMCA(GA, nmcs = 20)
background_geneset <- "data/kegg_legacy.gmt"
gene_rate <- getGeneRate(background.geneset = background_geneset, pathway.geneset = marker_list)	
regression_data <- doRegression(GA, gene.rate = gene_rate)
score_data <- computeScore(GA, regression.data = regression_data, weight = FALSE, normalize = "z-score")

roc_PaaSc <- rep(0, ncol(score_data))
names(roc_PaaSc) <- gsub('_', ' ', colnames(score_data))

for(i in 1:ncol(score_data)){
	all_label <- as.character(Idents(pbmc)[rownames(score_data)])
	all_label[all_label != names(roc_PaaSc)[i]] = 'others'
	roc_PaaSc[i] <- roc_summary(all_label = all_label, ref_label = 'others', score = score_data[,i])
}
score_PaaSc <- score_data

# CelliD
GA <- CreateSeuratObject(counts = gene.activities, project = "GA", min.cells = 3, min.features = 200)
GA <- NormalizeData(GA)
GA <- RunMCA(GA)
HGT <- RunCellHGT(GA, pathways = marker_list, dims = 1:50, n.features = 200)
score_data <- data.frame(t(HGT))

roc_CelliD <- rep(0, ncol(score_data))
names(roc_CelliD) <- gsub('_', ' ', colnames(score_data))

for(i in 1:ncol(score_data)){
	all_label <- as.character(Idents(pbmc)[rownames(score_data)])
	all_label[all_label != names(roc_CelliD)[i]] = 'others'
	roc_CelliD[i] <- roc_summary(all_label = all_label, ref_label = 'others', score = score_data[,i])
}
score_CelliD <- score_data

# GSdensity
GA <- CreateSeuratObject(counts = gene.activities, project = "GA", min.cells = 3, min.features = 200)
cells <- colnames(GA)
ce <- compute.mca(object = GA)
el <- compute.nn.edges(coembed = ce, nn.use = 300)
score_data <- run.rwr.list(el = el, gene_set_list = marker_list, cells = cells)

roc_GSdensity <- rep(0, ncol(score_data))
names(roc_GSdensity) <- gsub('_', ' ', colnames(score_data))

for(i in 1:ncol(score_data)){
	all_label <- as.character(Idents(pbmc)[rownames(score_data)])
	all_label[all_label != names(roc_GSdensity)[i]] = 'others'
	roc_GSdensity[i] <- roc_summary(all_label = all_label, ref_label = 'others', score = score_data[,i])
}
score_GSdensity <- score_data

# plot
out <- data.frame(cellType = names(marker_list), PaaSc = roc_PaaSc, CelliD = roc_CelliD, GSdensity = roc_GSdensity)
rownames(out) <- NULL

#         cellType  PaaSc CelliD GSdensity
#1       CD4_Naive 0.8982 0.9285    0.6705
#2         CD4_TCM 0.9617 0.9115    0.6402
#3       CD8_Naive 0.9735 0.9519    0.4872
#4       CD16_Mono 0.9544 0.9354    0.5223
#5              NK 0.9908 0.9815    0.9497
#6            Treg 0.7031 0.8430    0.5525
#7       CD14_Mono 0.9232 0.7518    0.7048
#8  B_intermediate 0.9740 0.9758    0.9508
#9            cDC2 0.8599 0.6457    0.5202
#10       B_memory 0.9669 0.9604    0.9482
#11        CD4_TEM 0.9215 0.9001    0.5604
#12           MAIT 0.9779 0.9726    0.8891
#13        CD8_TEM 0.9590 0.9773    0.9375
#14        B_naive 0.9746 0.9758    0.6921
#15            gdT 0.8753 0.8679    0.7696
#16            pDC 0.9909 0.4499    0.4595

# plot
df_long <- out %>%
  pivot_longer(cols = c(PaaSc, CelliD, GSdensity), 
               names_to = "Method", 
               values_to = "Accuracy") %>%
  mutate(Method = factor(Method, levels = c("PaaSc", "CelliD", "GSdensity")))


color_pat <- c("#E74C3C", "#F39C12", "#F1C40F", "#2ECC71", "#1ABC9C", "#3498DB", "#9B59B6", "#E91E63", "#C0392B", "#2980B9", "#27AE60", "#8E44AD", "#D68910", "#BDC3C7", "#34495E", "#2C3E50")

p1 <- ggplot(df_long, aes(x = Method, y = Accuracy)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +
  geom_point(aes(color = cellType), 
             size = 3, 
             alpha = 0.8,
             position = position_jitter(width = 0.2, seed = 123)) +
  scale_color_manual(values = color_pat, name = "Cell Type") +
  labs(title = "",
       x = "",
       y = "ROC") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5, size = 12),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 11),
        axis.title = element_text(size = 12),
        legend.title = element_text(size = 11),
        legend.text = element_text(size = 9)) +
  ylim(0, 1)

pdf('granulocyte_ROC2.pdf', width = 6, height = 5)
	print(p1)
dev.off()

save(score_PaaSc, score_CelliD, score_GSdensity, out, file = 'granulocyte_ROC.RData')

## Venn diagram
tmpList <- list(score_PaaSc, score_CelliD, score_GSdensity)
tmpPred <- matrix(0, ncol(pbmc), 3)
for(i in 1:3){
  MM <- tmpList[[i]]
  tmp <- colnames(MM)[apply(MM, 1, which.max)]
  tmp <- gsub('_', ' ', tmp)
  tmpPred[,i] <- (as.character(Idents(pbmc)) == tmp) + 0
}
dimnames(tmpPred) <- list(Cells(pbmc), c('PaaSc', 'CelliD', 'GSdensity'))

PaaSc_correct = which(tmpPred[, "PaaSc"] == 1)
CelliD_correct = which(tmpPred[, "CelliD"] == 1)
GSdensity_correct = which(tmpPred[, "GSdensity"] == 1)

venn.plot <- venn.diagram(
  x = list(
    PaaSc = PaaSc_correct,
    CelliD = CelliD_correct,
    GSdensity = GSdensity_correct
  ),
  category.names = c("PaaSc", "CelliD", "GSdensity"),
  filename = NULL,  
  output = TRUE,
  
  
  fill = c("lightblue", "lightgreen", "lightcoral"),
  alpha = 0.6,
  
 
  cex = 1.2,
  fontfamily = "serif",
  cat.cex = 1.2,
  cat.fontfamily = "serif",
  cat.default.pos = "outer",
  cat.pos = c(-27, 27, 135),
  cat.dist = c(0.055, 0.055, 0.085),
  
  
  lwd = 2,
  lty = 'blank'
)

pdf('granulocyte_Venn.pdf', width = 6, height = 6)
  grid.draw(venn.plot)
dev.off()


#### Fig7D

suppressPackageStartupMessages(library(msigdbr))
hallmark_genes <- msigdbr(species = "Homo sapiens", category = "H")
hallmark_list <- split(hallmark_genes$gene_symbol, hallmark_genes$gs_name)

gene_rate <- getGeneRate(background.geneset = background_geneset, pathway.geneset = hallmark_list)
GA=computeMCA(GA,nmcs = 20)
regression_data <- doRegression(GA, gene.rate = gene_rate)
score_data <- computeScore(GA, regression.data = regression_data, weight = FALSE, normalize = "z-score")



cc=score_data
                  
cc=cc%>%filter(Type %in% c('CD8 TEM','CD8 Naive'))
cc$Type[cc$Type == "CD8 Naive"] = 'Naive'
cc$Type[cc$Type == "CD8 TEM"] = 'Tem'

wilcox_results <- lapply(1:50, function(i) {
  test <- wilcox.test(cc[cc$Type == "Naive", i],
                      cc[cc$Type == "Tem", i])
  data.frame(
    Feature = colnames(cc)[i],
    p.value = test$p.value,
    W = test$statistic
  )
})

wilcox_df <- do.call(rbind, wilcox_results)
wilcox_df$FDR <- p.adjust(wilcox_df$p.value, method = "fdr")


wilcox_df <- wilcox_df %>%
  mutate(
    direction = ifelse(
      apply(cc[cc$Type == "Tem", Feature], 2, median) >
        apply(cc[cc$Type == "Naive", Feature], 2, median),
      "Tem_higher", "Naive_higher"
    )
  )

wilcox_df <- wilcox_df %>%
  mutate(
    logFDR = -log10(FDR),
    logFDR_signed = ifelse(direction == "Naive_higher",
                           -logFDR, logFDR)
  )


plot_df = wilcox_df[order(wilcox_df$FDR),]



ggplot(plot_df[c(1:10,41:50),], aes(x = reorder(Feature, logFDR_signed),
                    y = logFDR_signed,
                    fill = direction)) +
  geom_col() +
  coord_flip() +
  labs(x = "Feature", y = "log10FDR",
       title = "TEM vs Naive (Wilcox)") +
  scale_fill_manual(values = c("Tem_higher" = "firebrick",
                               "Naive_higher" = "steelblue")) +
  theme_minimal(base_size = 14)+
  theme(axis.text.x = element_text(color = "black", angle = 0, hjust = 0.5, size = 15),
          axis.text.y = element_text(color = "black", size = 15),
          legend.title = element_text(size = 15),
          legend.text = element_text(size = 15),
          panel.background = element_blank(),
          panel.border = element_rect(colour = "black", fill = NA, linewidth = 1.2))





