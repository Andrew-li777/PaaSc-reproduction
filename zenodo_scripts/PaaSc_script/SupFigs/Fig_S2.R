library(Seurat)
library(PaaSc)
library(tidyverse)
library(ROCit)

# Load the PBMC dataset
pbmc.data <- Read10X(data.dir = "/sibcb1/bioinformatics/shijiantao/rasc/pbmc3k/filtered_gene_bc_matrices/hg19/")
object <- CreateSeuratObject(counts = pbmc.data, project = "pbmc3k", min.cells = 3, min.features = 200)


#---SupFig2A---#
object <- computeMCA(object, nmcs = 50)
background_geneset <- "/sibcb1/bioinformatics/shijiantao/rasc/Tres.kegg"
B_cell_activation <- readLines("/sibcb1/bioinformatics/shijiantao/rasc/B_cell_activation.txt")
pathway_geneset <- list(B_cell_activation = B_cell_activation)
gene_rate <- getGeneRate(background.geneset = background_geneset, pathway.geneset = pathway_geneset)
rdata <- doRegression(object, gene.rate = gene_rate)
out <- data.frame(tvalue = c(rdata[1,], rdata[2,]),
				  pvalue = c(rdata[3,], rdata[4,]),
				  mca = colnames(rdata),
				  Factor = as.integer(gsub('mca_', '', colnames(rdata))),
				  group = c(rep('pathway', ncol(rdata)), rep('background', ncol(rdata)))
				  )

SupFig2A <- ggplot(out, aes(x = Factor, y = tvalue, fill = group)) + 
		 geom_bar(stat="identity", position = position_dodge(width = 0.9), width = 0.9) + 
		 theme_bw() + scale_fill_manual(values = c("grey", "red")) + 
	     theme(panel.grid = element_blank(),panel.background=element_blank(),
	        	axis.ticks.length = unit(5,'pt'),axis.line=element_line(color="black"),
	            axis.title.x = element_blank(),
	            axis.text.y = element_text(color="black",size=10),
	            axis.ticks = element_line(color="black"))
SupFig2A
#ggsave(file = 'SupFig1A.pdf',SupFig1A, width = 6, height = 3)

###compute PaaSc score
object <- computeMCA(object, nmcs = 20)
regression_data <- doRegression(object, gene.rate = gene_rate)
score_data <- computeScore(object, regression.data = regression_data, weight = FALSE, normalize = "z-score")
head(score_data)
# Ceiling score for better visualization
cscore <- score_data[colnames(object), 1]
cscore[cscore > quantile(cscore, 0.99)] = quantile(cscore, 0.99)
cscore[cscore < quantile(cscore, 0.01)] = quantile(cscore, 0.01)
object@meta.data$B_cell_activation <- cscore


# pre-process
object <- FindVariableFeatures(object, verbose = FALSE)
object <- ScaleData(object, verbose = FALSE)
object <- RunPCA(object, verbose = FALSE)
object <- FindNeighbors(object, verbose = FALSE)
object <- RunUMAP(object, dims = 1:20, verbose = FALSE)
annotation_file <- "/sibcb1/bioinformatics/shijiantao/rasc/pbmc3k/pbmc3k_seurat_annotation.txt.gz"
seurat_annotations <- read.table(annotation_file, header = TRUE, row.names = 1, sep = "\t")
object@meta.data$seurat_annotations <- seurat_annotations[sub("-.*", "", colnames(object)), ]


#---SupFig2B---#
colors=c('#E64825','#F18C25','#CAC0E1','#F6C0CC','#ABDAEC','#FAECA8','#FDD5C0','#6A9ACE','#1E803D','grey')
SupFig2B <- DimPlot(object, group.by = "seurat_annotations",cols = colors,size=5)+
                  labs(title = NULL)  + 
                  theme(aspect.ratio = 1,legend.position = 'top')
SupFig2B
#ggsave(file = 'SupFig2B.pdf',SupFig2B, width = 3, height = 3)


#---SupFig1C---#
SupFig2C <- FeaturePlot(object, features = "B_cell_activation", cols = c('white', '#E64825'))+ 
  labs(title = NULL)  + 
theme(aspect.ratio = 1,legend.position = 'top')
SupFig2C
#ggsave(file = 'SupFig2C.pdf',SupFig2C, width = 3, height = 3)


set.seed(1)
label_data <- doBinarization(score_data, n.cluster = 2, method = "GMM")
threshold <- min(score_data[rownames(label_data[label_data['B_cell_activation.label'] == 'positive', ]), 'B_cell_activation'])
#---SupFig1D---#
SupFig2D <- ggplot(score_data, aes(x = B_cell_activation)) + geom_density(fill = "grey", alpha = 0.7) + theme_bw() +
	theme(panel.grid = element_blank(),aspect.ratio = 1,
		  axis.line = element_line(color = "black"),
	      axis.text = element_text(color = "black")) +
	geom_vline(xintercept = threshold, color = "red", linetype = "dashed") + xlim(-5,5) +
	labs(x = "Value", y = "Density")
SupFig2D
#ggsave(file = 'SupFig2D.pdf',SupFig2D, width = 3, height = 3)


object@meta.data$B_cell_activation_label <- label_data[colnames(object), ]
#---SupFig2E---#
SupFig2E <- DimPlot(object, group.by = "B_cell_activation_label",cols=c('#BE0E23','#073E7F'))+
  labs(title = NULL)  + 
theme(aspect.ratio = 1,legend.position = 'top')
SupFig2E


score_data_withanno <- cbind(score_data, seurat_annotations)
#---SupFig2F---#
SupFig2F<- ggplot(score_data_withanno, aes(x = seurat_annotations, y = B_cell_activation, color = seurat_annotations)) + 
    scale_color_manual(values = colors)+ ylim(-6,6) +
    labs(x=NULL,y='PaaSc score')+
    geom_boxplot(outlier.shape = NA) + geom_jitter(width = 0.2, size = 0.8,shape = 1) + theme_bw() + 
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
    theme(plot.title = element_text(hjust = 0.5)) + theme(legend.title=element_blank()) + 
    theme(legend.position="none") +  
    theme(aspect.ratio = 1,axis.text.x = element_text(vjust = 1, hjust = 1, angle = 90)
         )
SupFig2F
#ggsave(file = 'SupFig2F.pdf',SupFig2F, width = 3, height = 3)


label <- rep('others', nrow(score_data_withanno))
label[score_data_withanno$seurat_annotations == "B"] <- "B"
roc_empirical <- rocit(score = score_data_withanno$B_cell_activation, class = label, negref = "others")

auc_value <- as.numeric(strsplit(as.character(summary(roc_empirical)[4,]), ":")[[1]][2])
auc_label <- paste0("AUC: ", round(auc_value, 3))  

roc_df <- data.frame(
  TPR = roc_empirical$TPR,  
  FPR = roc_empirical$FPR   
)
#---SupFig2G---#
SupFig2G <- ggplot(roc_df, aes(x = FPR, y = TPR)) +
  geom_line(color = "red", linewidth = 1) +  # 画 ROC 曲线
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black") +  
  annotate("text", x = 0.2, y = 0.6, label = auc_label, size = 5, color = "blue") +  
  labs(title = "ROC Curve", x = "1 - Specificity", y = "Sensitivity") +
  theme_minimal()+ theme(aspect.ratio = 1,panel.border = element_rect(color = "black", fill = NA, linewidth = 1))
SupFig2G
#ggsave(file = 'SupFig2G.pdf',SupFig2G, width = 3, height = 3)