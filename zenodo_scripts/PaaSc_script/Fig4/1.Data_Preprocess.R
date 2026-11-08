suppressPackageStartupMessages(library(PaaSc))
suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(reticulate))
use_condaenv(condaenv = "/sibcb1/bioinformatics/liaoxiqi/miniconda3/envs/R4.1/bin/python", required = TRUE)
ad <- import("anndata")

#-----[GSE175533]-----#
adata <- ad$read_h5ad("/sibcb1/bioinformatics/liaoxiqi/PaaSc_Data/Senescence/GSE175533_sceasy_hay.h5ad") 
counts=t(as.matrix(adata$X))
colnames(counts) <- rownames(adata$obs)
rownames(counts) <- rownames(adata$var)
GSE175533_object <- Seurat::CreateSeuratObject(counts = counts, project = "GSE175533",
                                         meta.data =  adata$obs,
                                         min.cells = 0, min.features = 0)
GSE175533_object = computeMCA(GSE175533_object, nmcs = 20)
saveRDS(GSE175533_object,file = '/sibcb1/bioinformatics/liaoxiqi/PaaSc_Data/Senescence/GSE175533_WI38.RDS')

#-----[GSE115301_10X]-----#
GSE115301_10X_counts = read.table('/sibcb1/bioinformatics/liaoxiqi/PaaSc_Data/Senescence/GSE115301_Growing_Sen_10x_count.txt.gz')
GSE115301_10X_meta = read.table('/sibcb1/bioinformatics/liaoxiqi/PaaSc_Data/Senescence/GSE115301_Growing_Sen_10x_metadata.txt.gz')
GSE115301_10X_object = CreateSeuratObject(counts = GSE115301_10X_counts, project = "GSE115301_10X", min.cells = 0, min.features = 0)
GSE115301_10X_object@meta.data$Group = ifelse(GSE115301_10X_meta$Condition1 == 'Senescence1','Sen','Growing')
GSE115301_10X_object = computeMCA(GSE115301_10X_object, nmcs = 20)
saveRDS(GSE115301_10X_object,file = '/sibcb1/bioinformatics/liaoxiqi/PaaSc_Data/Senescence/GSE115301_IMR90_10X.RDS')

#-----[GSE115301_timecourse]-----#
GSE115301_Timecourse_counts=read.table('/sibcb1/bioinformatics/liaoxiqi/PaaSc_Data/Senescence/GSE115301_Time_course_counts.txt',header = T,row.names = 1)
GSE115301_Timecourse_object <- CreateSeuratObject(counts = GSE115301_Timecourse_counts, project = "GSE115301_Timecourse", min.cells = 0, min.features = 0)
GSE115301_Timecourse_object@meta.data = GSE115301_Timecourse_object@meta.data %>% mutate(Time_course = case_when(
                    str_detect(colnames(GSE115301_Timecourse_object), "Sen") ~ "Senescent", 
                    str_detect(colnames(GSE115301_Timecourse_object), "Growing") ~ "Growing",
                    str_detect(colnames(GSE115301_Timecourse_object), "day2") ~ "Day2", 
                    str_detect(colnames(GSE115301_Timecourse_object), "day4") ~ "Day4",  
                    TRUE ~ NA_character_ )) %>% mutate(Group = ifelse(Time_course == 'Senescent','Sen','Growing'))
GSE115301_Timecourse_object = computeMCA(GSE115301_Timecourse_object, nmcs = 20)
saveRDS(GSE115301_Timecourse_object,file = '/sibcb1/bioinformatics/liaoxiqi/PaaSc_Data/Senescence/GSE115301_IMR90_Timecourse.RDS')

#-----[GSE119807]-----#
Growing_counts = read.table('/sibcb1/bioinformatics/liaoxiqi/PaaSc_Data/Senescence/GSM3384106_LowPDCtrl.dge.txt.gz',header = T)
Sen_counts = read.table('/sibcb1/bioinformatics/liaoxiqi/PaaSc_Data/Senescence/GSM3384107_LowPD50Gy.dge.txt.gz',header = T)
counts = inner_join(Growing_counts,Sen_counts,by = 'GENE') %>% column_to_rownames(var = 'GENE')
GSE119807_object = CreateSeuratObject(counts = counts, project = "GSE119807", min.cells = 0, min.features = 0)
GSE119807_object@meta.data$Group = rep(c('Growing','Sen'),each = 400)
GSE119807_object = computeMCA(GSE119807_object, nmcs = 20)
saveRDS(GSE119807_object,file = '/sibcb1/bioinformatics/liaoxiqi/PaaSc_Data/Senescence/GSE119807_HCA2.RDS')

#-----[GSE102090]-----#
suppressPackageStartupMessages(library(Matrix))
GENE = read.table("/sibcb1/bioinformatics/liaoxiqi/PaaSc_Data/Senescence/GSE102090_genes.tsv.gz")
Growing_barcodes = readLines('/sibcb1/bioinformatics/liaoxiqi/PaaSc_Data/Senescence/GSM2723761_D2P4_barcodes.tsv.gz')
Growing_counts = readMM('/sibcb1/bioinformatics/liaoxiqi/PaaSc_Data/Senescence/GSM2723761_D2P4_matrix.mtx.gz')
colnames(Growing_counts) = Growing_barcodes
Sen_counts = readMM('/sibcb1/bioinformatics/liaoxiqi/PaaSc_Data/Senescence/GSM2723762_D2P16_matrix.mtx.gz')
Sen_barcodes = readLines('/sibcb1/bioinformatics/liaoxiqi/PaaSc_Data/Senescence/GSM2723762_D2P16_barcodes.tsv.gz')
colnames(Sen_counts) = Sen_barcodes
counts = cbind(Growing_counts,Sen_counts)
rownames(counts) = GENE$V2
GSE102090_object = CreateSeuratObject(counts = counts, project = "GSE102090",min.cells = 0, min.features = 0)
GSE102090_object@meta.data$Group = c(rep("Growing", length(Growing_barcodes)), rep("Sen", length(Sen_barcodes)))
GSE102090_object = computeMCA(GSE102090_object, nmcs = 20)
saveRDS(GSE102090_object,file = '/sibcb1/bioinformatics/liaoxiqi/PaaSc_Data/Senescence/GSE102090_HUVEC.RDS')
