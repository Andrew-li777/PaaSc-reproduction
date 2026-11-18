#source /sibcb/program/install/r-4.4/profile
#source /sibcb/program/install/hdf5-1.10/profile


suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(PaaSc))
suppressPackageStartupMessages(library(CelliD))
suppressPackageStartupMessages(library(gsdensity))
suppressPackageStartupMessages(library(ROCit))


library(SeuratDisk)

Convert(
  source = "/sibcb1/bioinformatics/liaoxiqi/project/PaaSc/MouseBrain/Reproducibility/ProcessedData/SMI_Lung/Lung5-3.h5ad",
  dest = "~/Lung5-3.h5seurat", 
  assay = "RNA",
  overwrite = TRUE
)


library(Seurat)
adata <- LoadH5Seurat("~/Lung5-3.h5seurat")

library(readxl)

#marker
obj=readRDS('/sibcb1/bioinformatics/liaoxiqi/PaaSc_Script/test/ST_Lung.RDS')
#marker=read_xlsx('/sibcb1/bioinformatics/liaoxiqi/PaaSc_Script/test/data.xlsx')
#marker = read.csv('/sibcb1/bioinformatics/fengyan/Project/fig9/st4.csv')
#gene_sets <- map(marker, ~ unique(na.omit(.x)))

#library(dplyr)
#obj$cell_type_group <- obj$cell_type
#obj$cell_type_group <- case_when(
#  obj$cell_type %in% c("T CD4 naive", "T CD4 memory", "Treg") ~ "CD4.T.cells",
#  obj$cell_type %in% c("T CD8 naive", "T CD8 memory") ~ "CD8.T.cells",
#  grepl("^tumor", obj$cell_type) ~ "Tumor",
#  obj$cell_type == "NK" ~ "NK.cells",
#  obj$cell_type == "Neutrophils" | obj$cell_type == "neutrophil" ~ "Neutrophils",
#  obj$cell_type == "CD14.Monocytes" | obj$cell_type == "monocyte" ~ "CD14.Monocytes",
#  obj$cell_type == "mDC" ~ "cDC",
#  obj$cell_type == "pDC" ~ "pDC",
#  obj$cell_type == "macrophage" ~ "Macrophage",
# obj$cell_type == "mast" ~ "Mast.cell",
#  obj$cell_type == "B-cell" ~ "B.cell",
#  obj$cell_type == "epithelial" ~ "Epithelial.cell",
#  obj$cell_type == "endothelial" ~ "Endothelial.cell",
#  obj$cell_type == "fibroblast" ~ "Fibroblast",
#  obj$cell_type == "plasmablast" ~ "Plasma.cell",
#  TRUE ~ as.character(obj$cell_type)  # 其它原样保留
#)
#table(obj$cell_type_group)


load("marker_list.rdata")
gene_sets = marker_list

###PaaSc
#bgList <- loadGMT('/sibcb1/bioinformatics/shijiantao/scATACseq/data/Tres.kegg')
gmt <- gmtPathways("/sibcb1/bioinformatics/shijiantao/scATACseq/data/Tres.kegg")
bgList = gmt
gene_rate = getGeneRate(background.geneset = bgList, pathway.geneset = gene_sets)
regression_data = doRegression(obj, gene.rate = gene_rate)
PaaSc.res = computeScore(obj, regression.data = regression_data, pvalue = 0.05, weight = FALSE, normalize = "z-score")
colnames(PaaSc.res)=names( gene_sets)
write.table(PaaSc.res,'spatial_lung_PaaSc.txt')



###CelliD
set.seed(123)
    HGT <-
        RunCellHGT(
          X = obj,
            reduction = "mca",
            pathways = gene_sets,
            log.trans = T, minSize = 5,
            dims=1:20
        )
CelliD.res = HGT%>%t()%>%as.data.frame()
write.table(CelliD.res,'spatial_lung_CelliD.txt')



###GSDensity
ce <- compute.mca(object = obj)
cells <- colnames(obj)
el <- compute.nn.edges(coembed = ce, nn.use = 300)

# Seurat annotations
Idents(obj)=obj@meta.data$cell_type_group
refLabel <- as.character(Idents(obj))
names(refLabel) <- Cells(obj)



roc <- data.frame()
for(i in 1:length(gene_sets)){
	tag = names(gene_sets)[i]
    print(tag)
	gs = gene_sets[[i]]
    GSD <- run.rwr(el = el, gene_set = gs, cells = cells)
    GSDensity.res <- as.data.frame(GSD)%>%setnames('GSDensity')
    check = sum(refLabel == tag) > 20
    print(check)
    #print(check)
    if(check){
		label = rep('others', length(PaaSc.res))
		label[rownames(PaaSc.res) %in% names(refLabel)[refLabel == tag]] = tag
        label[is.na(label)] <- "others"
        #print(label)
        print(length(PaaSc.res[,tag]))
        print(length(label))
		PaaSc.roc_empirical <- rocit(score = PaaSc.res[,tag], class = label, negref = "others")
        roc[tag,'PaaSc'] <- PaaSc.roc_empirical$AUC

        CelliD.roc_empirical <- rocit(score = CelliD.res[,tag], class = label, negref = "others")
        roc[tag,'CelliD'] <- CelliD.roc_empirical$AUC
     
        GSDensity.roc_empirical <- rocit(score = GSDensity.res$GSDensity, class = label, negref = "others")
        roc[tag,'GSDensity'] <- GSDensity.roc_empirical$AUC
	}

}

###PLOT
data=roc%>%rownames_to_column(var='Celltype')%>%pivot_longer(.,col=-1,names_to = 'method',values_to = 'AUC')

data <- data %>%
  group_by(method) %>%
  mutate(median_auc = median(AUC, na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(method = factor(method, levels = unique(method[order(-median_auc)])))


# cell type
celltypes <- c(
  "B.cell", "NK.cells", "Endothelial.cell", "Epithelial.cell",
  "Fibroblast", "cDC", "Macrophage", "Mast.cell",
  "CD14.Monocytes", "Neutrophils", "Plasma.cell",
  "CD4.T.cells", "CD8.T.cells", "Tumor"
)

colors <- c(
  "#E74C3C", "#F39C12", "#F1C40F", "#2ECC71",
  "#1ABC9C", "#3498DB", "#9B59B6", "#E91E63",
  "#C0392B", "#2980B9", "#8E44AD",
  "#D68910", "#BDC3C7", "#34495E"
)
celltype_colors <- setNames(colors, celltypes)


ggplot(data, aes(x = method, y = AUC)) +
    geom_boxplot(outlier.shape = NA,color='black',width=0.6, position = position_dodge(width = 0.5)) +
    geom_jitter(aes(color = Celltype),size=2.3,width=0.2) +
    guides(title = NULL) +
    scale_color_manual(values=celltype_colors)+
    labs(x = NULL ,y = 'AUC', title = 'NSCLC') +
    theme_minimal() +
    theme(
        axis.text.x = element_text(hjust = 1, vjust = 1, size = 10),
        axis.title = element_text(size = 12, face = "bold"),
        plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
        legend.position = "right",
        panel.background = element_blank(),
        plot.background = element_blank(),
        panel.border = element_rect(fill = NA, linewidth = 1),
        panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank()
    ) 




####hallmark
hallmark_genes <- msigdbr(species = "Homo Sapiens", category = "H")
hallmark_list <- split(hallmark_genes$gene_symbol, hallmark_genes$gs_name)

gene_rate = getGeneRate(background.geneset = bgList, pathway.geneset = hallmark_list)
regression_data = doRegression(obj, gene.rate = gene_rate)
hallmark.res = computeScore(obj, regression.data = regression_data, pvalue = 0.05, weight = FALSE, normalize = "z-score")

write.csv(hallmark.res,'ST_lung_hallmark_paasc.csv')


#SMI
spatial=as.data.frame(adata$obsm)
x <- spatial$spatial.1
y <- spatial$spatial.2
b <- calculateSpatialMI(x = x, y = y, z = hallmark.res)
bb <- calculateSpatialMI(x = x, y = y, z = PaaSc.res)


##cor
# b is the HALLMARK PaaSc score, bb is the celltype marker PaaSc score
res <- matrix(0, nrow = ncol(bb), ncol = ncol(b))  

for (i in 1:ncol(bb)) {
  for (j in 1:ncol(b)) {
    tmp <- suppressWarnings(cor.test(bb[, i], b[, j], method = "pearson")) 
    if (!is.na(tmp$p.value) && tmp$p.value <= 0.05) {
      res[i, j] <- tmp$estimate  
    } else {
      res[i, j] <- 0
    }
  }
}


res_df <- as.data.frame(res)
rownames(res_df) <- colnames(bb)
colnames(res_df) <- colnames(b)


Heatmap(res,
        name = "score",
        col = col_fun,
        cluster_rows = FALSE,
        cluster_columns = FALSE,
        show_row_names = TRUE,
        show_column_names = TRUE,
        height = unit(9, "cm"),
        width = unit(5, "cm"),
        row_names_gp = gpar(fontsize = 4),
        column_names_gp = gpar(fontsize = 8),
        heatmap_legend_param = list(
            at = c(-0.6, -0.3, 0, 0.3, 0.6, 1),
            labels = c("-0.6","-0.3","0","0.3","0.6","1")
        )
)




