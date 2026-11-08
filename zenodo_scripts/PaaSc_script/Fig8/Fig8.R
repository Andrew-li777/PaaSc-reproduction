suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(PaaSc))
suppressPackageStartupMessages(library(CelliD))
suppressPackageStartupMessages(library(gsdensity))
suppressPackageStartupMessages(library(ROCit))
suppressPackageStartupMessages(library(BiocParallel))


obj=readRDS('ST_braincortex.RDS')
obj


#marker=readRDS('/sibcb1/bioinformatics/liaoxiqi/PaaSc_Script/test/ST_mouseBrain_markers.RDS')
Cellname=read.csv('GSE144136_CellNames.csv.gz',col.names = 1)
Genenames=read.csv('GSE144136_GeneNames.csv.gz',col.names = 1)
library(Matrix)
counts=readMM('GSE144136_GeneBarcodeMatrix_Annotated.mtx.gz')%>%as.data.frame()
rownames(counts)=Genenames$X1
colnames(counts)=Cellname$X1
obj=CreateSeuratObject(counts=counts)

extract_part <- function(s) {
  dot_pos <- regexpr("\\.", s)
  if (dot_pos == -1) return(s) 
  
  prefix <- substr(s, 1, dot_pos - 1) 
  underscores <- gregexpr("_", prefix)[[1]]  
  
  if (length(underscores) >= 2) {
    start_pos <- underscores[2] + 1
    substr(prefix, start_pos, nchar(prefix))
  } else {
    prefix
  }
}

results <- sapply(rownames(obj@meta.data), function(s) {
  res <- extract_part(s)
  #if (is.na(res)) "" else res  # 将 NA 替换为空字符串
})
names(results) <- NULL

results[results %in% c('OPCs_1', 'OPCs_2','OPCs_3')] = 'OPCs'
results[results %in% c('Mix_1', 'Mix_2','Mix_3','Mix_4','Mix_5')] = 'others'
results[results %in% c('Oligos_1', 'Oligos_2','Oligos_3')] = 'Olig'
results[results == 'VIP'] = 'inhibitory Vip'
results[results == 'SST'] = 'inhibitory Sst'
results[results == 'L5'] = 'Excitatory L5'
results[results == 'L_6'] = 'Excitatory L6'
results[results == 'L2_4'] = 'Excitatory L2/3'
results[results == 'L5_6'] = 'Excitatory L5'
results[results %in% c('L4_6','L4_5')] = 'Excitatory L4'
results[results %in% c('Astros_3','Astros_2','Astros_1')] = 'Astro'
results[results == 'PVALB'] = 'Pvalb'
results[results %in% c('Inhib_5','Inhib_1')] = 'inhibitory Other'

Idents(obj) <- results
UC <- as.character(unique(results))
DefaultAssay(obj) <- "RNA"
obj <- NormalizeData(obj, verbose=F)
marker <- list()
for(i in 1:length(UC)){
	tag <- UC[i]
	NS <- sum(results == tag)
    print(tag)
    print(NS)
	if(NS > 20){
		m1 = FindMarkers(obj, ident.1 = tag, only.pos = TRUE, test.use = "wilcox")		
		marker[[tag]] <- m1
	}
}


# geneset
genesets_top200 <- lapply(marker, function(genevec) {
  genes=rownames(genevec[1:200,]) 
return(genes)
})

saveRDS(genesets_top200, file = 'ST_Braincortex_top200_markers.RDS')



###CelliD
set.seed(123)
    HGT <-
        RunCellHGT(
          X = obj,
            reduction = "mca",
            pathways = genesets_top200,
            log.trans = T, minSize = 5,
            dims=1:20
        )
CelliD.res = HGT%>%t()%>%as.data.frame()
write.csv(CelliD.res,'/sibcb1/bioinformatics/liaoxiqi/PaaSc_Script/Fig8/ST_braincortex_CelliD.csv')


###PaaSc
bgList <- 'Tres.kegg'
gene_rate = getGeneRate(background.geneset = bgList, pathway.geneset = genesets_top200)
regression_data = doRegression(obj, gene.rate = gene_rate)
PaaSc.res = computeScore(obj, regression.data = regression_data, pvalue = 0.05, weight = FALSE, normalize = "z-score")
colnames(PaaSc.res)=names(genesets_top200)
write.csv(PaaSc.res,'/sibcb1/bioinformatics/liaoxiqi/PaaSc_Script/Fig8/ST_braincortex_PaaSc.csv')

###GSDensity
ce <- compute.mca(object = obj)
cells <- colnames(obj)
el <- compute.nn.edges(coembed = ce, nn.use = 300)
GSDensity.res=list()
for(i in 1:length(genesets_top200)){
	tag = names(genesets_top200)[i]
    print(tag)
	gs = genesets_top200[[i]]
    GSD <- run.rwr(el = el, gene_set = gs, cells = cells)
    GSDensity.res[[tag]]=GSD
	}


GSDensity.res=bind_cols(GSDensity.res)%>%as.data.frame()
rownames(GSDensity.res)=colnames(obj)
head(GSDensity.res)
write.csv(GSDensity.res,'/sibcb1/bioinformatics/liaoxiqi/PaaSc_Script/Fig8/ST_braincortex_GSDensity.csv')

refLabel <- obj@meta.data$layer%>%as.character()
names(refLabel) <- Cells(obj)
PaaSc.res=PaaSc.res[names(refLabel),c('Excitatory L2/3','Excitatory L4','Excitatory L5','Excitatory L6')]
names(PaaSc.res)=c('L2/3','L4','L5','L6')
CelliD.res=CelliD.res[names(refLabel),c('Excitatory L2/3','Excitatory L4','Excitatory L5','Excitatory L6')]
names(CelliD.res)=c('L2/3','L4','L5','L6')
GSDensity.res=GSDensity.res[names(refLabel),c('Excitatory L2/3','Excitatory L4','Excitatory L5','Excitatory L6')]
names(GSDensity.res)=c('L2/3','L4','L5','L6')


refLabel[refLabel %in% c('Layer2','Layer3')] = 'L2/3'
refLabel[refLabel == 'Layer5'] = 'L5'
refLabel[refLabel == 'Layer6'] = 'L6'
refLabel[refLabel == 'Layer3'] = 'L3'
refLabel[refLabel == 'Layer4'] = 'L4'

a=genesets_top200[c('Excitatory L2/3','Excitatory L4','Excitatory L5','Excitatory L6')]
names(a)=c('L2/3','L4','L5','L6')
refLabel=na.omit(refLabel)


roc <- data.frame()
for(i in 1:length(a)){
	tag = names(a)[i]
    check = sum(refLabel == tag) > 20
    if(check){
		label = rep('others', length(refLabel))
		label[rownames(PaaSc.res) %in% names(refLabel)[refLabel == tag]] = tag
        label[is.na(label)] <- "others"
      
		PaaSc.roc_empirical <- rocit(score = PaaSc.res[,tag], class = label, negref = "others")
        roc[tag,'PaaSc'] <- PaaSc.roc_empirical$AUC

        CelliD.roc_empirical <- rocit(score = CelliD.res[,tag], class = label, negref = "others")
        roc[tag,'CelliD'] <- CelliD.roc_empirical$AUC
     
        GSDensity.roc_empirical <- rocit(score = GSDensity.res[,tag], class = label, negref = "others")
        roc[tag,'GSDensity'] <- GSDensity.roc_empirical$AUC
	}
}
roc


data <- roc %>%
  rownames_to_column("Layer") %>%
  pivot_longer(-Layer, names_to = "method", values_to = "AUC") %>%
  group_by(method) %>%
  mutate(median_auc = median(AUC, na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(method = factor(method, levels = unique(method[order(-median_auc)])))

save(PaaSc.res,CelliD.res,GSDensity.res,roc,file = 'braincortex_ROC.RData')


olors = c(

   '#90C1F0', 
    '#FF9A92', 
    '#FFDDA0',  
    '#FFA85E',
    '#3882C9',
    '#C595D4'
)

library(ggpubr)
p=ggplot(data, aes(x = method, y = AUC)) +
  geom_boxplot(outlier.shape = NA,color='grey') +
  geom_jitter(aes(color = Layer),size=6) +
  guides(title = NULL) +
scale_color_manual(values=colors)+
  labs(x = NULL, y = 'AUC', title = 'brain cortex') +
  theme_minimal() +
  theme(
    axis.text.x = element_text( hjust = 1, vjust = 1, size = 10),
    axis.title = element_text(size = 12, face = "bold"),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    legend.position = "right",
    panel.background = element_blank(),
    plot.background = element_blank(),
    panel.border = element_rect(fill = NA, linewidth = 1),
    panel.grid.major.y = element_line(color = "grey90", linewidth = 0.3),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank()
  ) 
p


pdf('braincortex_boxplot.pdf',width = 6,height = 5)
p
dev.off()