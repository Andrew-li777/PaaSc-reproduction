suppressPackageStartupMessages(library(PaaSc))
suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(CelliD))
suppressPackageStartupMessages(library(gsdensity))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(data.table))

#----- benchmark -----#
background_file = "/sibcb1/bioinformatics/shijiantao/rasc/Tres.kegg"
pathway_genes = readLines("/sibcb1/bioinformatics/shijiantao/rasc/data/Sen_signature_NC.txt")
gene_rate = getGeneRate(background.geneset = background_file, pathway.geneset = list(SenMayo = pathway_genes))

# load data
GSE115301_10X_object = readRDS('/sibcb1/bioinformatics/liaoxiqi/PaaSc_Data/Senescence/GSE115301_IMR90_10X.RDS')
GSE115301_Timecourse_object = readRDS('/sibcb1/bioinformatics/liaoxiqi/PaaSc_Data/Senescence/GSE115301_IMR90_Timecourse.RDS')
GSE119807_object = readRDS('/sibcb1/bioinformatics/liaoxiqi/PaaSc_Data/Senescence/GSE119807_HCA2.RDS')
GSE102090_object = readRDS('/sibcb1/bioinformatics/liaoxiqi/PaaSc_Data/Senescence/GSE102090_HUVEC.RDS')
GSE175533_object  = readRDS('/sibcb1/bioinformatics/liaoxiqi/PaaSc_Data/Senescence/GSE175533_WI38.RDS')

# Create a function to caculate activity score
get.Score <- function(glist,  object){
    ###PaaSc
    regression_data = doRegression(object, gene.rate = gene_rate)
    PaaSc.res = computeScore(object, regression.data = regression_data, pvalue = 0.05, weight = FALSE, normalize = "z-score")%>%setnames('PaaSc')
 
    ###GSDensity
    ce <- compute.mca(object = object)
    cells <- colnames(object)
    el <- compute.nn.edges(coembed = ce, nn.use = 300)
    GSD <- run.rwr(el = el, gene_set = glist, cells = cells)
    GSD.res <- as.data.frame(GSD)%>%setnames('GSDensity')
    
    ###CelliD
    set.seed(123)
    HGT <-
        RunCellHGT(
          X = object,
            reduction = "mca",
            pathways = list(SenMayo = pathway_genes),
            log.trans = T, minSize = 5,
            dims=1:20
        )
    CelliD.res = HGT%>%as.matrix()%>%as.data.frame()%>%setnames('CelliD')

    score = cbind(PaaSc.res,GSD.res)%>%cbind(CelliD.res)%>%mutate(Group = object@meta.data$Group)
      return(score)
}

# get score
GSE115301_10X = get.Score(pathway_genes,GSE115301_10X_object)
GSE115301_Timecourse = get.Score(pathway_genes,GSE115301_Timecourse_object)
GSE119807 = get.Score(pathway_genes,GSE119807_object)
GSE102090 = get.Score(pathway_genes,GSE102090_object)
GSE175533 = get.Score(pathway_genes,GSE175533_object)

# save
write.table(GSE115301_10X,'/sibcb1/bioinformatics/liaoxiqi/PaaSc_Data/Senescence/GSE115301_10X_Score.tsv',sep='\t',col.names=T,row.names=T,quote=F)
write.table(GSE115301_Timecourse,'/sibcb1/bioinformatics/liaoxiqi/PaaSc_Data/Senescence/GSE115301_Timecourse_Score.tsv',sep='\t',col.names=T,row.names=T,quote=F)
write.table(GSE119807,'/sibcb1/bioinformatics/liaoxiqi/PaaSc_Data/Senescence/GSE119807_Score.tsv',sep='\t',col.names=T,row.names=T,quote=F)
write.table(GSE102090,'/sibcb1/bioinformatics/liaoxiqi/PaaSc_Data/Senescence/GSE102090_Score.tsv',sep='\t',col.names=T,row.names=T,quote=F)
write.table(GSE175533,'/sibcb1/bioinformatics/liaoxiqi/PaaSc_Data/Senescence/GSE175533_Score.tsv',sep='\t',col.names=T,row.names=T,quote=F)



#----- Hallmark & SenMayo -----#
suppressPackageStartupMessages(library(msigdbr))
suppressPackageStartupMessages(library(pbapply))                   

hallmark_genes <- msigdbr(species = "Homo sapiens", category = "H")
hallmark_list <- split(hallmark_genes$gene_symbol, hallmark_genes$gs_name)
hallmark_list[["SenMayo"]] <- readLines('/sibcb1/bioinformatics/shijiantao/rasc/data/Sen_signature_NC.txt')

background_file = '/sibcb1/bioinformatics/shijiantao/rasc/Tres.kegg'
gene_rate = getGeneRate(background.geneset = background_file, pathway.geneset = hallmark_list,mode = "single")

filelist=list.files('/sibcb1/bioinformatics/shijiantao/Signaling/MCA','.RDS',full.names = T)
filelist=filelist[!grepl("mouse", filelist)]

pblapply(filelist,function(x){
    sample=sub(".*/(.*?)\\.RDS$", "\\1", x)
    object=readRDS(x)
    regression_data = doRegression(object, gene.rate = gene_rate)
    score_data = computeScore(object,regression.data = regression_data, pvalue = 0.05, weight = FALSE, normalize = "z-score")
    print(paste('Processed :',sample))
write.table(score_data,file=paste0('/sibcb1/bioinformatics/liaoxiqi/PaaSc_Data/Senescence/SenMayo_Hallmark_score/',sample,'_PaaSc.tsv'),sep='\t',row.names=T,col.names=T,quote=F)
})