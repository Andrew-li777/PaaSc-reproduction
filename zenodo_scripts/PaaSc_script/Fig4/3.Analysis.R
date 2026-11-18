suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(clusterProfiler))
suppressPackageStartupMessages(library(pbapply))

#----- GSEA -----#
file_list <- list.files(pattern = "*.txt$",path = '/sibcb1/bioinformatics/liaoxiqi/project/senescence/SenCID/BTD/DESeq2',full.names= TRUE) 
sen_gmt<-clusterProfiler::read.gmt("/sibcb1/bioinformatics/liaoxiqi/project/senescence/SenCID/BTD/Senescence.gmt")

#calculate NES and Pvalue
gsea_res=lapply(file_list,function(x){
            cell_name <- gsub("_DESeq2\\.txt$", "", basename(x))  
            data <- read.table(x,header=T,sep='\t')%>%column_to_rownames(var="Gene")%>%select("log2FoldChange")%>% arrange(desc(log2FoldChange))%>% na.omit()
            id <- data$log2FoldChange
            names(id) <- rownames(data)
            gsea<-clusterProfiler::GSEA(geneList = id,
                                    TERM2GENE = sen_gmt,
                                    verbose = T,
                                    pvalueCutoff = 0.5)%>%as.data.frame()%>%mutate(sample=cell_name)
                 return(gsea)
})%>%bind_rows()


write.table(gsea_res,'/sibcb1/bioinformatics/liaoxiqi/PaaSc_Data/Senescence/bulk_SenMayo_GSEA.tsv',sep='\t',col.names=T,row.names=F,quote=F)

#----- monocle -----#
suppressPackageStartupMessages(library(monocle))

GSE115301_Timecourse_object = readRDS('/sibcb1/bioinformatics/liaoxiqi/PaaSc_Data/Senescence/GSE115301_IMR90_Timecourse.RDS')
DE = readLines('/sibcb1/bioinformatics/liaoxiqi/PaaSc_Data/Senescence/DE.txt')

HSMM_sample_sheet <- data.frame(sample=colnames(GSE115301_Timecourse_object),row.names=colnames(GSE115301_Timecourse_object))
HSMM_gene_annotation <- data.frame(gene_short_name = row.names(GSE115301_Timecourse_object), row.names = row.names(GSE115301_Timecourse_object))

pd <- new("AnnotatedDataFrame", data = HSMM_sample_sheet)
fd <- new("AnnotatedDataFrame", data = HSMM_gene_annotation)

HSMM <- newCellDataSet(as(as.matrix(GSE115301_Timecourse_object@assays$RNA@counts), "sparseMatrix"),
                phenoData = pd,
                featureData = fd,
                lowerDetectionLimit = 0.5,
                expressionFamily = negbinomial.size())
HSMM <- estimateSizeFactors(HSMM)
HSMM <- estimateDispersions(HSMM)
HSMM <- detectGenes(HSMM, min_expr = 0.1)
HSMM <- setOrderingFilter(HSMM, DE)
HSMM <- reduceDimension(HSMM,max_components =2,method = 'DDRTree')

PaaScscore = read.table('/sibcb1/bioinformatics/liaoxiqi/PaaSc_Data/Senescence/GSE115301_Timecourse_Score.tsv')%>%pull(PaaSc)

PaaScscore[PaaScscore > quantile(PaaScscore, 0.99)] = quantile(PaaScscore, 0.99)
PaaScscore[PaaScscore < quantile(PaaScscore, 0.01)] = quantile(PaaScscore, 0.01)

pData(HSMM)<-pData(HSMM)%>% mutate(time_course = GSE115301_Timecourse_object@meta.data$Time_course)%>%mutate(PaaSc = PaaScscore)
pData(HSMM)$time_course = factor(pData(HSMM)$time_course ,levels = c('Growing','Day2','Day4','Senescent'))

saveRDS(HSMM,file='/sibcb1/bioinformatics/liaoxiqi/PaaSc_Data/Senescence/GSE115301_Timecourse_monocle.RDS')

#----- correlation -----#
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(msigdbr))
suppressPackageStartupMessages(library(pbapply))

hallmarks <- msigdbr(species = "Homo sapiens", category = "H")%>%pull(gs_name)%>%unique()
scorefilelist <- list.files(pattern = "*_PaaSc.txt",path = "/sibcb1/bioinformatics/liaoxiqi/project/PaaSc/CD8T/Senecence_SenMayo",full.names = T) 

###Function of compute correlation by PID
corByPID <- function(data){
  pid <- unique(data$pid)
  res <- list()
  for(id in pid){
    df <- filter(data, pid == id)
    if(nrow(df) > 50){
      tmp <- cor.test(df$Sene, df$Hallmark)
      res[[id]] <- list(pid = id, cor = tmp$estimate, pvalue = tmp$p.value, N = nrow(df))
    }
  }
  out <- data.frame(pid = sapply(res, function(x) x$pid),
            N = sapply(res, function(x) x$N),
            cor = sapply(res, function(x) x$cor),
            pvalue = sapply(res, function(x) x$pvalue))
            rownames(out)=NULL
  return(out)
}

###computing by celltype and dataset
Cell = c('CD8T','CD4Tconv','B','NK','Mono/Macro')
correlation_metrics <-pblapply(scorefilelist,function(x){
    sample=sub(".*/(.*?)_PaaSc\\.txt$", "\\1", x)
    print(paste('Processing:',sample))
    scores=suppressWarnings(fread(x)) %>%column_to_rownames(var='V1')
    metafile <- paste0("/sibcb1/bioinformatics/hongyuyang/dataset/Tres/2.tisch_data/0.raw_data/",sample, "_CellMetainfo_table.tsv")
    meta <- fread(metafile)
    if (!("Patient" %in% colnames(meta))) {
            meta <- meta %>% rename(Patient = Sample)       
        }
    celldf=lapply(Cell,function(c){
        if (c %in% meta[['Celltype (major-lineage)']]){
            cells=meta%>%dplyr::filter(`Celltype (major-lineage)` == c)
            score_data=scores[rownames(scores) %in% cells$Cell, ]
            Sene = as.numeric(score_data$SenMayo) 
            hallmarkdf=lapply(hallmarks,function(h){
                Hallmark = as.numeric(score_data[[h]])
                pid <- as.character(cells$Patient)
                df <- data.frame(Sene, Hallmark, pid)
                final = corByPID(df)%>%mutate(hallmark=h,celltype=c,dataset=sample)
                return(final)
            })
            hallmarkdf = do.call(rbind, hallmarkdf)
            return(hallmarkdf)
        }else{
            return(NULL)
      } 
    })
    return(celldf)
})
correlation_df= bind_rows( correlation_metrics) 

write.table(correlation_df,'/sibcb1/bioinformatics/liaoxiqi/PaaSc_Data/senescence/all_pid_celltype_hallmarkVSsene_correlation.tsv',sep='\t',colnames=T,row.names=F,quote=F)