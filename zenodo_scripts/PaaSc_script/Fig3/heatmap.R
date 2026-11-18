library(tidyverse)

CD8T=read.table('/sibcb1/bioinformatics/fengyan/Project/data/CD8T_auc_values.txt')
rownames(CD8T)=NULL
colnames(CD8T)[2]='GSDensity'
colnames(CD8T)[5]='GSEA'
CD8T=CD8T%>%
        mutate(celltype='CD8 T',method=colnames(.)[apply(., 1, which.max)])
head(CD8T)
CD8T$method=factor(CD8T$method,levels = c('PaaSc','scGSEA','CelliD','GSDensity','AUCell','GSEA','VAM','GSVA'))


custom_colors <- setNames(
   c("#E89DA0", "#88CEE6", "#F6C8A8", "#B2D3A4", 
            "#9FBA95", "#E6CECF", "#B696B6", "#80C1C4"),  
    c('PaaSc','scGSEA','GSDensity','VAM','AUCell','CelliD','GSVA','GSEA') 
)
custom_colors
row_anno <- rowAnnotation(
    Method = CD8T[, 10],
    col = list(Method = custom_colors),
    annotation_name_gp = gpar(fontsize = 0),
    simple_anno_size = unit(1.3, "cm"),
    annotation_legend_param = list(
         Method= list(
            title = "Method",
            direction = "horizontal",
            title_position = "topcenter",
            legend_width = unit(15, "cm")  
        )
    )
    
)
p_CD8T=Heatmap(as.matrix(CD8T[, c('PaaSc','scGSEA','CelliD','GSDensity','AUCell','GSEA','VAM','GSVA')]), 
        row_split = CD8T[, 10],  
        col = colorRamp2(c(min(CD8T[, 1:8]), 0.8, max(CD8T[, 1:8])), 
                        c("#6395C7", "white", "firebrick3")),
        name = "AUC",
          cluster_columns = F,
        right_annotation = row_anno,
           #show_annotation_name = FALSE,
         
        show_row_names = FALSE, 
        cluster_row_slices = FALSE,  
        show_parent_dend_line = FALSE,
        show_row_dend = FALSE,
        show_column_dend = FALSE,
        heatmap_legend_param = list(
            direction = "horizontal",  
            title_position = "topcenter" 
        )
         ) 

pdf('CD8T.pdf',width = 5,height = 8)
draw(p_CD8T, heatmap_legend_side = "top")
dev.off()



Macro=read.table('/sibcb1/bioinformatics/fengyan/Project/data/Macro_auc_values.txt')
rownames(Macro)=NULL
colnames(Macro)[2]='GSDensity'
colnames(Macro)[5]='GSEA'

Macro=Macro%>%
        mutate(celltype='Macrophage',method=colnames(.)[apply(., 1, which.max)])
head(Macro)
Macro$method=factor(Macro$method,levels = c('PaaSc','scGSEA','GSDensity','VAM','AUCell','CelliD','GSEA','GSVA'))

row_anno <- rowAnnotation(
    Method = Macro[, 10],
    col = list(Method = custom_colors),
    annotation_name_gp = gpar(fontsize = 0),
    simple_anno_size = unit(1.3, "cm"),
    annotation_legend_param = list(
         Method= list(
            title = "Method",
            direction = "horizontal",
            title_position = "topcenter",
            legend_width = unit(15, "cm")  
        )
    )
    
)

p_Macro=Heatmap(as.matrix(Macro[, c('PaaSc','scGSEA','GSDensity','VAM','AUCell','CelliD','GSEA','GSVA')]), 
        row_split = Macro[, 10],  
        col = colorRamp2(c(min(CD8T[, 1:8]), 0.8, max(CD8T[, 1:8])), 
                        c("#6395C7", "white", "firebrick3")),
        name = "AUC",
          cluster_columns = F,
        right_annotation = row_anno,
           #show_annotation_name = FALSE,
         
        show_row_names = FALSE, 
        cluster_row_slices = FALSE,  
        show_parent_dend_line = FALSE,
        show_row_dend = FALSE,
        show_column_dend = FALSE,
        heatmap_legend_param = list(
            direction = "horizontal", 
            title_position = "topcenter"  
        )
         ) 
pdf('Macro.pdf',width = 5,height = 8)
draw(p_Macro, heatmap_legend_side = "top")
dev.off()


NK=read.table('/sibcb1/bioinformatics/fengyan/Project/data/NK_auc_values.txt')#%>%mutate(celltype='NK')
rownames(NK)=NULL
#colnames(NK)[2]='GSDensity'
NK=NK%>%mutate(celltype='NK',method=colnames(.)[apply(., 1, which.max)])
head(NK)
NK$method=factor(NK$method,levels = c('PaaSc','scGSEA','GSDensity','GSEA','AUCell','CelliD','VAM','GSVA'))
#head(NK)

row_anno <- rowAnnotation(
    Method = NK[, 10],
    col = list(Method = custom_colors),
    annotation_name_gp = gpar(fontsize = 0),
    simple_anno_size = unit(1.3, "cm"),
    annotation_legend_param = list(
         Method= list(
            title = "Method",
            direction = "horizontal",
            title_position = "topcenter",
            legend_width = unit(15, "cm")  
        )
    )
    
)

p_NK=Heatmap(as.matrix(NK[, c('PaaSc','scGSEA','GSDensity','GSEA','AUCell','CelliD','VAM','GSVA')]), 
        row_split = NK[, 10],  
        col = colorRamp2(c(min(CD8T[, 1:8]), 0.8, max(CD8T[, 1:8])), 
                        c("#6395C7", "white", "firebrick3")),
        name = "AUC",
          cluster_columns = F,
        right_annotation = row_anno,
           #show_annotation_name = FALSE,
         
        show_row_names = FALSE, 
        cluster_row_slices = FALSE,  
        show_parent_dend_line = FALSE,
        show_row_dend = FALSE,
        show_column_dend = FALSE,
        heatmap_legend_param = list(
            direction = "horizontal", 
            title_position = "topcenter" 
        )
         ) 

pdf('NK.pdf',width = 5,height = 8)
draw(p_NK, heatmap_legend_side = "top")
dev.off()


Treg=read.table('/sibcb1/bioinformatics/fengyan/Project/data/Treg_auc_values.txt')
rownames(Treg)=NULL
#colnames(NK)[2]='GSDensity'
Treg=Treg%>%mutate(celltype='Treg',method=colnames(.)[apply(., 1, which.max)])
head(Treg)
Treg$method=factor(Treg$method,levels = c('PaaSc','CelliD','VAM','GSVA','scGSEA','GSDensity','GSEA','AUCell'))
row_anno <- rowAnnotation(
    Method = Treg[, 10],
    col = list(Method = custom_colors),
    annotation_name_gp = gpar(fontsize = 0),
    simple_anno_size = unit(1.3, "cm"),
    annotation_legend_param = list(
         Method= list(
            title = "Method",
            direction = "horizontal",
            title_position = "topcenter",
            legend_width = unit(15, "cm") 
        )
    )
    
)

p_Treg=Heatmap(as.matrix(Treg[,  c('PaaSc','CelliD','VAM','GSVA','scGSEA','GSDensity','GSEA','AUCell')]), 
        row_split = Treg[, 10],  
        col = colorRamp2(c(min(CD8T[, 1:8]), 0.8, max(CD8T[, 1:8])), 
                        c("#6395C7", "white", "firebrick3")),
        name = "AUC",
          cluster_columns = F,
        right_annotation = row_anno,
           #show_annotation_name = FALSE,
         
        show_row_names = FALSE, 
        cluster_row_slices = FALSE,  
        show_parent_dend_line = FALSE,
        show_row_dend = FALSE,
        show_column_dend = FALSE,
        heatmap_legend_param = list(
            direction = "horizontal",  
            title_position = "topcenter"  
        )
         ) 

pdf('Treg.pdf',width = 5,height = 8)
draw(p_Treg, heatmap_legend_side = "top")
dev.off()