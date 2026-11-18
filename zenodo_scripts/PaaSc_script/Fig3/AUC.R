library(tidyverse)
library(data.table)
library(pROC)

AUCell_filelist <- list.files(pattern = "AUCell",path = "AUCell",full.names = T) 
#head(AUCell_filelist)
length(AUCell_filelist)
CelliD_filelist <- list.files(pattern = "CelliD",path = "CelliD",full.names = T) 
#CelliD_filelist
length(CelliD_filelist)
GDSdensity_filelist <- list.files(pattern = "GSDensity",path = "GSDensity",full.names = T) 
#GDSdensity_filelist
length(GDSdensity_filelist)
GSVA_filelist <- list.files(pattern = "GSVA",path = "GSVA",full.names = T) 
#GSVA_filelist
length(GSVA_filelist)
PaaSc_filelist <- list.files(pattern = "PaaSc",path = "PaaSc",full.names = T) 
#PaaSc_filelist
#head(PaaSc_filelist)
length(PaaSc_filelist)
ssGSEA_filelist <- list.files(pattern = "GSEA",path = "GSEA",full.names = T) 
#ssGSEA_filelist
length(ssGSEA_filelist)
VAM_filelist <- list.files(pattern = "VAM",path = "VAM/",full.names = T) 
length(VAM_filelist)
scGSEA_filelist <- list.files(pattern = "scGSEA",path = "scGSEA/",full.names = T) 
length(scGSEA_filelist)

CD8T_filename<-read.table("auc_name.txt")
head(CD8T_filename)
length(CD8T_filename)
CD8T_filename=as.character(CD8T_filename$V1)

auc_values<-data.frame()
#i=29
for (i in 1:length(CD8T_filename)){
    sample=CD8T_filename[i]
    print(sample)
###meta
    meta_filename<-paste0("/sibcb1/bioinformatics/hongyuyang/dataset/Tres/2.tisch_data/0.raw_data/",sample,"_CellMetainfo_table.tsv")
    meta<-fread(meta_filename)%>%
        rename(celltype='Celltype (major-lineage)')%>%
        mutate(celltype=ifelse(celltype== "DC","DC","others"))%>%select(Cell,celltype)
#head(meta)
    ###PaaSC
    PaaSC_filename<-paste0("PaaSc/",sample,"_PaaSc.txt")
    PaaSC_score<-fread(PaaSC_filename)
    ###GDS
    GDS_filename<-paste0("GSDensity/",sample,"_GSDensity.txt")
    GDS_score<-read.table(GDS_filename,sep="\t",header=T)
    ###CelliD
    CelliD_filename<-paste0("CelliD/",sample,"_CelliD.txt")
    CelliD_score<-fread(CelliD_filename)
    ###AUCell
    AUCell_filename<-paste0("AUCell/",sample,"_AUCell.txt")
    AUCell_score<-fread(AUCell_filename)
    ###ssGSEA
    ssGSEA_filename<-paste0("GSEA/",sample,"_GSEA.txt")
    ssGSEA_score<-fread(ssGSEA_filename)
    ###GSVA
    GSVA_filename<-paste0("GSVA/",sample,"_GSVA.txt")
    GSVA_score<-fread(GSVA_filename)
    ###VAM
    VAM_filename<-paste0("VAM/",sample,"_VAM.txt")
    VAM_score<-fread(VAM_filename)
    ###scGSEA
    scGSEA_filename<-paste0("scGSEA/",sample,"_scGSEA.txt")
    scGSEA_score<-fread(scGSEA_filename)
head(PaaSC_score);head(GDS_score);head(CelliD_score);head(AUCell_score);
head(ssGSEA_score);head(GSVA_score);head(VAM_score);head(scGSEA_score)
    ###score
    score<-inner_join(PaaSC_score, GDS_score, by = "Cell") %>%
         inner_join(CelliD_score, by = "Cell") %>%
         inner_join(AUCell_score, by = "Cell")%>%
         inner_join(ssGSEA_score, by = "Cell")%>%
         inner_join(GSVA_score, by = "Cell")%>%
         inner_join(VAM_score, by = "Cell")%>%
         inner_join(scGSEA_score, by = "Cell")%>%
         inner_join(meta,by="Cell")%>%as.data.frame()
#head(score)
    #write.table(score, paste0("~/project/PaaSc/out/CD8T_cell_marker/score/",sample,"_score.txt"), sep = '\t', row.names = FALSE, col.names = TRUE, quote = FALSE)
    methods <- c("PaaSc","GSDensity","CelliD" ,"AUCell",    "GSEA",      "GSVA"  , "VAM" , "scGSEA")
    for (n in seq_along(methods)) {
      method=methods[n]  
      roc_curve <- roc(score$celltype, score[[method]], 
                  levels = c("others", "DC"))
      auc_values[sample,method] <- auc(roc_curve)
    }
}
