library(tidyverse)
library(Seurat)
library(PaaSc)
library(reticulate)
use_condaenv(condaenv = "/sibcb1/bioinformatics/liaoxiqi/miniconda3/envs/R4.1/bin/python", required = TRUE)
ad <- import("anndata")


adata <- ad$read_h5ad("tabula-muris-senis-facs-official-raw-obj.h5ad") 
colnames(adata$X) <- adata$var$features

celltype = c('proerythroblast','granulocyte monocyte progenitor cell','classical monocyte','dendritic cell',"CD4-positive, alpha-beta T cell","CD8-positive, alpha-beta T cell",'regulatory T cell','medium spiny neuron','interneuron','neuron','oligodendrocyte','oligodendrocyte precursor cell','Bergmann glial cell','astrocyte','hepatocyte','chondrocyte','bladder cell','ventricular myocyte','pancreatic B cell')

cells = adata$obs %>% filter(cell_ontology_class %in% celltype) %>% rownames()

obj <- Seurat::CreateSeuratObject(counts = t(as.matrix(adata$X))[,cells], project = "TMS facs",
                                         meta.data =  adata$obs[cells,],
                                         min.cells = 0, min.features = 0)
obj = computeMCA(obj,nmcs = 20)
saveRDS(obj,file = 'TMS_facs_selected.RDS')

library(CelliD)
library(GSDensity)


mouse_gene_sets = readRDS('mouse_traits_geneset.RDS')
keggs = readRDS('mouse_kegg_geneset.RDS')


###PaaSc
gene_rate = getGeneRate(background.geneset = keggs, pathway.geneset = mouse_gene_sets,mode ='single')
regression_data = doRegression(obj, gene.rate = gene_rate)
score_data = computeScore(obj, regression.data = regression_data, pvalue = 0.05, weight = FALSE, normalize = "z-score")

write.table(score_data,'TMS_facs_selected_alltraits_PaaSc_score.tsv',sep='\t',col.names=T,row.names=T,quote=F)


###GSDensity
ce <- compute.mca(object = obj)
cells <- colnames(obj)
el <- compute.nn.edges(coembed = ce, nn.use = 300)
GSD <- run.rwr(el = el, gene_set = mouse_gene_sets, cells = cells)
GSD.res <- as.data.frame(GSD)
write.table(GSDensity.res,'TMS_facs_selected_alltraits_GSDensity_score.tsv',sep='\t',col.names=T,row.names=T,quote=F)


###CelliD
set.seed(123)
HGT <-RunCellHGT(X = obj,
            reduction = "mca",
            pathways = mouse_gene_sets,
            log.trans = T, minSize = 5,
            dims=1:20
        )
CelliD.res = HGT%>%t()%>%as.data.frame()
write.table(CelliD.res,'TMS_facs_selected_alltraits_CelliD_score.tsv',sep='\t',col.names=T,row.names=T,quote=F)


###Fig5A
Fig5A_df = PaaSc.res %>% mutate(celltype=ct_anno$abbreviation[match(meta$cell_ontology_class, ct_anno$CellType)])
Fig5A_df$celltype <- reorder(Fig5A_df$celltype, Fig5A_df[['UKB_460K.blood_LYMPHOCYTE_COUNT']], FUN = median, decreasing = TRUE)
Fig5A_df$tissue <- ct_anno$tissue[match(Fig5A_df$celltype, ct_anno$abbreviation)]

Fig5A <- ggplot(Fig5A_df, aes(x = celltype, y = UKB_460K.blood_LYMPHOCYTE_COUNT,color= tissue)) +
  geom_boxplot() +
  guides(title = NULL)+
  labs(x = NULL, y = 'PaaSc score', title = NULL) +
  scale_color_manual(values = c("blood/immune" = "#C17188", "brain" = "#9AB846", "others" = "#F5A433")) +
  theme_minimal() +
  theme(
    axis.text.x = element_text( angle = 60,hjust = 1, vjust = 1),
    panel.background = element_blank(),
    plot.background = element_blank(),
    panel.border = element_rect( fill = NA, linewidth = 1),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )  +coord_cartesian(ylim = c(-4, 3)) 
Fig5A



###### Fig5B
t='UKB_460K.blood_LYMPHOCYTE_COUNT'
threshold <- min(PaaSc.res[rownames(label_data[label_data[paste0(t,'.label')] == 'positive', ]), t])

Fig5B<-ggplot(PaaSc.res, aes(x = !!sym(t))) +
  geom_histogram(aes(y = after_stat(density)), 
                 binwidth = 0.2, 
                 fill = "slategray2", 
                 color = "black", 
                 alpha = 0.7, 
                 position = "identity") + 
  geom_density(color = "dodgerblue4", linewidth = 1) + 
  labs(title = NULL, x = NULL, y = "Density") + 
  annotate("text", x = 1, y = 0.7, label = paste0("threshold = ", signif(threshold, 4)), hjust = 0, vjust = 0, color = "black", size = 5) + 
  geom_vline(xintercept = threshold, color = "dodgerblue4", linetype = "dashed", linewidth = 1.5) + 
  theme_minimal() +  
  theme(
    panel.grid = element_blank(),  
    axis.title = element_text(size = 14),  
    axis.text = element_text(size = 12),   
    axis.line = element_line(linewidth = 1)    
  ) +
  scale_x_continuous(expand = c(0, 0), limits = c(-3, 3)) + 
  scale_y_continuous(expand = c(0, 0))

Fig5B


###Fig5C
Fig5C_df$Odds_ratio[Fig5C_df$Odds_ratio == 0] <- 1e-6
Fig5C_df$log_odds_ratio <- log(Fig5C_df$Odds_ratio)
Fig5C_df$padj = p.adjust(Fig5C_df$Pvalue,method = 'fdr')
Fig5C_df$padj[Fig5C_df$padj < 1e-5 ] <- 1e-5
Fig5C_df$logFDR = -log(Fig5C_df$padj)
Fig5C_df$tissue <- ct_anno$tissue[match(Fig5C_df$CellType, ct_anno$abbreviation)]
Fig5C_df$CellType <- factor(Fig5C_df$CellType, levels = Fig5C_df$CellType[order(Fig5C_df$log_odds_ratio)])


Fig5C <- ggplot(Fig5C_df, aes(x=log_odds_ratio, y=CellType))+geom_point( aes(color= logFDR,shape=tissue),fill=NA,size=4)+
  theme_minimal() + 
  theme(
       panel.border = element_rect( fill = NA, linewidth = 1),
    panel.grid = element_blank(),
    axis.text.x = element_text( hjust = 1),
    axis.title.y = element_text(angle = 0, hjust = 1), 
    axis.text.y = element_text(hjust = 1,size = 8)         
  ) +
geom_vline(aes(xintercept = 0),color="gray",linetype="dashed", size = 0.5)+
scale_color_distiller(palette ='Spectral',name ="-Log10(FDR)")+
  xlab("Log10(Odds Ratio)") +                                        
  ylab("") +                                  
  coord_fixed() 
Fig5C


###Fig5D
#PaaSc
gs_color=rev(c(rep('#C17188',9),rep('#9AB846',6),rep('#F5A433',7)))
ct_color=c(rep('#C17188',7),rep('#9AB846',7),rep('#F5A433',5))
Fig5D_PaaSc_df$ct_col = factor(ct_anno$tissue[match(Fig5D_PaaSc_df$CellType, ct_anno$abbreviation)],levels=c('blood/immune','brain','others'))
Fig5D_PaaSc_df$gs_col = factor(gs_anno$tissue[match(Fig5D_PaaSc_df$Traits, gs_anno$abbreviation)],levels=c('others','brain','blood/immune'))
Fig5D_PaaSc_df$Mean[Fig5D_PaaSc_df$Mean >3 ] =3
Fig5D_PaaSc_df$Mean[Fig5D_PaaSc_df$Mean < -3 ] = -3
Fig5D_PaaSc = ggplot(Fig5D_PaaSc_df, aes(x = interaction(CellType,ct_col), y = interaction(Traits,gs_col))) +
                        guides(y = "axis_nested",x = "axis_nested")+
                        geom_point(aes(fill=Mean),color='white',shape=22,size=5)+                               
                        scale_fill_distiller(palette ='RdBu')+
                        geom_text(aes(label = ifelse(FDR < 0.05, "*", "")), size = 5, color = "#696969") +
                        labs(title = 'PaaSc',x=NULL,y=NULL)+
                        theme_minimal() + 
                        theme(
                            panel.grid = element_blank(),                  
                            axis.text.x = element_text(angle = 45, hjust = 1, 
                                                   color = ct_color), 
                            axis.title.y = element_text(angle = 0, hjust = 1),
                            axis.text.y = element_text(hjust = 1, size = 7, 
                                                   color = gs_color),
                            ggh4x.axis.nestline.y = element_line(size = 0.5,colour = c("#F5A433","#9AB846","#C17188")), 
                            ggh4x.axis.nesttext.y = element_text(size = 9, colour = c("#F5A433","#9AB846","#C17188")), 
                            ggh4x.axis.nestline.x = element_line(size = 0.5,colour = c("#C17188","#9AB846","#F5A433")), 
                            ggh4x.axis.nesttext.x = element_text(size = 9, colour = c("#C17188","#9AB846","#F5A433"))
                      ) 
Fig5D_PaaSc


#GSDensity
Fig5D_GSDensity_df$ct_col = factor(ct_anno$tissue[match(Fig5D_GSDensity_df$CellType, ct_anno$abbreviation)],levels=c('blood/immune','brain','others'))
Fig5D_GSDensity_df$gs_col = factor(gs_anno$tissue[match(Fig5D_GSDensity_df$Traits, gs_anno$abbreviation)],levels=c('others','brain','blood/immune'))

Fig5D_GSDensity <- ggplot(Fig5D_GSDensity_df, aes(x = interaction(CellType,ct_col), y = interaction(Traits,gs_col))) +
        guides(x = "axis_nested",y="axis_nested")+
        geom_point(aes(fill=Mean),color='white',shape=22,size=5)+
        geom_text(aes(label = ifelse(FDR < 0.05, "*", "")), size = 5, color = "#696969") +
        scale_fill_gradient(high = "#d73027", low = "white")+
        geom_text(aes(label = ifelse(FDR < 0.05, "*", "")), size = 5, color = "#696969") +
                        labs(title = 'GSDensity',x=NULL,y=NULL)+
                        theme_minimal() + 
                        theme(
                            panel.grid = element_blank(),                  
                            axis.text.x = element_text(angle = 45, hjust = 1, 
                                                   color = ct_color), 
                            axis.title.y = element_text(angle = 0, hjust = 1),
                            axis.text.y = element_text(hjust = 1, size = 7, 
                                                   color = gs_color),
                            ggh4x.axis.nestline.y = element_line(size = 0.5,colour = c("#F5A433","#9AB846","#C17188")), 
                            ggh4x.axis.nesttext.y = element_text(size = 9, colour = c("#F5A433","#9AB846","#C17188")), 
                            ggh4x.axis.nestline.x = element_line(size = 0.5,colour = c("#C17188","#9AB846","#F5A433")), 
                            ggh4x.axis.nesttext.x = element_text(size = 9, colour = c("#C17188","#9AB846","#F5A433"))
                      ) 
Fig5D_GSDensity


#CelliD
Fig5D_CelliD_df$ct_col = factor(ct_anno$tissue[match(Fig5D_CelliD_df$CellType, ct_anno$abbreviation)],levels=c('blood/immune','brain','others'))
Fig5D_CelliD_df$gs_col = factor(gs_anno$tissue[match(Fig5D_CelliD_df$Traits, gs_anno$abbreviation)],levels=c('others','brain','blood/immune'))

Fig5D_CelliD <- ggplot(Fig5D_CelliD_df, aes(x = interaction(CellType,ct_col), y = interaction(Traits,gs_col))) +
        guides(x = "axis_nested",y= "axis_nested")+
        geom_point(aes(fill=Mean),color='white',shape=22,size=5)+
        scale_fill_gradient(high = "#d73027", low = "#FEFAFA")+
        geom_text(aes(label = ifelse(FDR < 0.05, "*", "")), size = 5, color = "#696969") +
                        labs(title = 'CelliD',x=NULL,y=NULL)+
                        theme_minimal() + 
                        theme(
                            panel.grid = element_blank(),                  
                            axis.text.x = element_text(angle = 45, hjust = 1, 
                                                   color = ct_color), 
                            axis.title.y = element_text(angle = 0, hjust = 1),
                            axis.text.y = element_text(hjust = 1, size = 7, 
                                                   color = gs_color),
                            ggh4x.axis.nestline.y = element_line(size = 0.5,colour = c("#F5A433","#9AB846","#C17188")), 
                            ggh4x.axis.nesttext.y = element_text(size = 9, colour = c("#F5A433","#9AB846","#C17188")), 
                            ggh4x.axis.nestline.x = element_line(size = 0.5,colour = c("#C17188","#9AB846","#F5A433")), 
                            ggh4x.axis.nesttext.x = element_text(size = 9, colour = c("#C17188","#9AB846","#F5A433"))
                      ) 
Fig5D_CelliD
