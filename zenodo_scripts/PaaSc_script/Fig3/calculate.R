#bPaths('/sibcb1/bioinformatics/liaoxiqi/miniconda3/envs/R4.1/bin/')
library("optparse")


option_list = list(make_option("--R", type = "character", default = NULL, help  = "The Seurat object in RDS format."),
                                  # make_option("--S", type = "character", default = NULL, help = "A signature file."),
                                   make_option("--O", type = "character", default = NULL, help = "Background geneset."),
                                  # make_option("--B", type = "character", default = NULL, help = "The output folder."),
				  # make_option("--", type = "character", default = NULL, help = "Geneset name."),                                            
                                  make_option("--tag", type = "character", default = NULL, help = "Output tag"))
args <- parse_args(OptionParser(option_list=option_list))



### check parameters
#check <- any(sapply(c('R', 'O', 'B', 'N','tag'), function(z) is.null(args[[z]])))
#if(check)
#  stop('\nPlease set all required parameters [R, B ,O , tag].\n')
#check = all(sapply(c('S'), function(z) is.null(args[[z]])))
#if(check)
#  stop('\n Please specify a signature file or ID file.\n')


library(BiocParallel)
library(GSVA)
library(Seurat)
# library(SeuratData)
library(dplyr)
library(ggrepel)
library(ggplot2)
# to be bench marked
library(gsdensity)
library(CelliD)
# library(escape)
library(VAM)
library(DescTools)
library(GSEABase)
library(AUCell)
library(nnet)
library(reshape2)
library(RColorBrewer)
library(PaaSc)

library(gficf)
library(PaaSc)
library(Seurat)

### loading data 
tag <- paste0(args$tag)
object <- readRDS(args$R)


#preprocess_sr <- function(sr.fname){
#  sr <- readRDS(sr.fname)
  sr <- object %>%
    NormalizeData() %>% 
    FindVariableFeatures() %>%
    ScaleData() %>%
    RunPCA() %>%
    FindNeighbors() %>%
    RunUMAP(dims = 1:10)
#  return(sr)
#}


#sr <- preprocess_sr(args$tag)
# print(sr)
# pre-compute
# gsdensity
ce <- compute.mca(object = sr)
  # print(head(ce))
  cells <- colnames(sr)
  el <- compute.nn.edges(coembed = ce, nn.use = 300)
  # print(head(el))
  # cellid
  sr <- RunMCA(sr)
  # aucell
  cells_rankings <- AUCell_buildRankings(sr@assays$RNA@data, plotStats = F, splitByBlocks = TRUE)


object =sr
############
signature = read.table("DC_marker.txt",header=F)
glist = as.character(signature$V1)

  gsd <- run.rwr(el = el, gene_set = glist, cells = cells)
  gsd.res <- as.data.frame(gsd)
  
gsd.res$Cell=rownames(gsd.res)
colnames(gsd.res)[1]="GSDensity"
write.table(gsd.res, file = paste0("GSDensity/",args$tag,"_GSDensity.txt"), sep = '\t', row.names =F, col.names = TRUE, quote = FALSE)


  # gsea
  #gsea <- enrichIt(obj = object, 
  #                 gene.sets = list(glist = glist), groups = ncol(object))
  mt <- match(CharacterList(list(glist)), rownames(object))
  mapdgenesets <- as.list(mt[!is.na(mt)])
  
  a <- GSVA:::.gsva(expr = object@assays$RNA@data, gset.idx.list = mapdgenesets, abs.ranking = F,tau = 0.25, rnaseq = T,
                    method = "ssgsea", kcdf = "Poisson")
  gsea.res <- as.data.frame(t(a))
  rownames(gsea.res) <- colnames(object)
  colnames(gsea.res) <- 'gsea'
  
gsea.res$Cell=rownames(gsea.res)
colnames(gsea.res)[1]="GSEA"
write.table(gsea.res, file = paste0("GSEA/",args$tag,"_GSEA.txt"), sep = '\t', row.names =F, col.names = TRUE, quote = FALSE)
  


  b <- GSVA:::.gsva(expr = object@assays$RNA@data, gset.idx.list = mapdgenesets, method = "gsva", rnaseq = F, 
                    kernel = F, parallel.sz = 40, BPPARAM = SnowParam(), verbose = T)
  gsva.res <- as.data.frame(t(b))
  rownames(gsva.res) <- colnames(object)
  colnames(gsva.res) <- 'gsva'

 gsva.res$Cell=rownames(gsva.res)
colnames(gsva.res)[1]="GSVA"
write.table(gsva.res, file = paste0("GSVA/",args$tag,"_GSVA.txt"), sep = '\t', row.names =F, col.names = TRUE, quote = FALSE)


  # cellid
  
  HGT_glist <- RunCellHGT(object, pathways = list(glist = glist), dims = 1:30, minSize = 5)
  cellid.res <- as.data.frame(HGT_glist)
  colnames(cellid.res) <- 'cellid'
 cellid.res$Cell=rownames(cellid.res)
colnames(cellid.res)[1]="CelliD"
write.table(cellid.res, file = paste0("CelliD/",args$tag,"_CelliD.txt"), sep = '\t', row.names =F, col.names = TRUE, quote = FALSE)
  #
  
  geneSets <- list(geneSet1 = glist)
  cells_AUC <- AUCell_calcAUC(geneSets, rankings = cells_rankings)
  cells_AUC_df <- as.data.frame(t(cells_AUC@assays@data@listData$AUC))
  aucell.res <- as.data.frame(cells_AUC_df)
 aucell.res$Cell=rownames(aucell.res)
colnames(aucell.res)[1]="AUCell"
write.table(aucell.res, file = paste0("AUCell/",args$tag,"_AUCell.txt"), sep = '\t', row.names =F, col.names = TRUE, quote = FALSE)

  
  coll <- list()
  coll[[1]] <- intersect(rownames(object), glist)
  names(coll)[1] <- "VAM"
  gene.set.collection = createGeneSetCollection(gene.ids=rownames(object),
                                                gene.set.collection=coll)
  obj <- vamForSeurat(seurat.data=object,
                      gene.set.collection=gene.set.collection,
                      center=F, gamma=T, sample.cov=F, return.dist=F)
  vam.res <- as.data.frame(t(obj@assays$VAMcdf@data))
 vam.res$Cell=rownames(vam.res)
colnames(vam.res)[1]="VAM"
write.table(vam.res, file = paste0("VAM/",args$tag,"_VAM.txt"), sep = '\t', row.names =F, col.names = TRUE, quote = FALSE)




# Step1: compute MCA loading and embedding
#object <- CreateSeuratObject(counts =sr[['RNA']]@counts, project = tag, min.cells = 3, min.features = 200)
#object = computeMCA(object,nmcs=20)

# Step2: Get gene occurrence rate in background and signature genesets
#background_file = args$B
#signature = read.table("DC_marker.txt",header=T)
pathway_file = as.list(signature)
#pathway_name = sub(".*/(.*)\\.txt$", "\\1", pathway_file)

gene_rate2 = getGeneRate(background.geneset = "../../PaaSc/Tres.kegg", pathway.geneset = pathway_file)
regression_data2 = doRegression(object, gene.rate = gene_rate2)
paasc.res = computeScore(object, regression.data = regression_data2)
paasc.res$Cell=rownames(paasc.res)
colnames(paasc.res)[1]="PaaSc"
write.table(paasc.res, file = paste0("PaaSc/",args$tag,"_PaaSc.txt"), sep = '\t', row.names =F, col.names = TRUE, quote = FALSE)


#################################################
  library(Matrix)
  sparse_matrix <- object@assays$RNA@counts
  data <- gficf(M = sparse_matrix)

  data <- runPCA(data = data,dim = 10,use.odgenes = T)
  data <-runReduction(data = data,reduction = "umap",nt = 2,verbose = T,n_neighbors=150)
  data$embedded$ccl = sapply(
                            strsplit(x = rownames(data$embedded),
                                     split = "_",fixed = T)
                            ,function(x) x[1]
                          )

  gf_icf_sparse <- data$gficf
  # 执行 NMF 分解
  tmp <- RcppML::nmf(
    A = gf_icf_sparse,  # 传递稀疏矩阵
    k = 100,            # 设置分解的维度
    nonneg = TRUE,      # 启用非负约束
    maxit = 100,        # 最大迭代次数
    tol = 1e-04,        # 收敛容忍度
    verbose = TRUE      # 输出调试信息
  )



  data$scgsea = list()
  data$scgsea$nmf.w <- Matrix::Matrix(data = tmp$w, sparse = T)
  data$scgsea$nmf.h <- t(Matrix::Matrix(data = tmp$h, sparse = T))
        data$scgsea$es <- NULL
        data$scgsea$nes <- NULL
        data$scgsea$pval <- NULL
        data$scgsea$fdr <- NULL
        data$scgsea$pathways <- NULL
        data$scgsea$x <- NULL
        data$scgsea$stat <- NULL

    data$scgsea$pathways = as.list(signature)
    data$scgsea$es = Matrix::Matrix(data = 0, nrow = length(data$scgsea$pathways), ncol = ncol(data$scgsea$nmf.w))
    data$scgsea$nes = Matrix::Matrix(data = 0, nrow = length(data$scgsea$pathways), ncol = ncol(data$scgsea$nmf.w))
    data$scgsea$pval = Matrix::Matrix(data = 0, nrow = length(data$scgsea$pathways), ncol = ncol(data$scgsea$nmf.w))
    data$scgsea$fdr = Matrix::Matrix(data = 0, nrow = length(data$scgsea$pathways), ncol = ncol(data$scgsea$nmf.w))
    rownames(data$scgsea$es) = rownames(data$scgsea$nes) = rownames(data$scgsea$pval) = rownames(data$scgsea$fdr) = names(data$scgsea$pathways)

    nsim = 10000; nt = 0; minSize = 15; maxSize = Inf;verbose = TRUE; seed = 180582; nmf.k = 100; fdr.th = 0.05; gp = 0
    
    pb = utils::txtProgressBar(min = 0, max = ncol(data$scgsea$nmf.w), initial = 0, style = 3)
    nt_fgsea <- ceiling(length(data$scgsea$pathways)/100)
    nt_fgsea <- ifelse(nt_fgsea > nt, nt, nt_fgsea)
    bpparameters <- BiocParallel::SnowParam(nt_fgsea)
    
    rownames(data$scgsea$nmf.w)=rownames(data$rawCounts)
    for (i in 1:ncol(data$scgsea$nmf.w)) {
        df = as.data.frame(fgsea::fgseaMultilevel(pathways = data$scgsea$pathways, 
            stats = as.matrix(data$scgsea$nmf.w)[, i], nPermSimple = nsim, 
            gseaParam = gp, BPPARAM = bpparameters, minSize = minSize, 
            maxSize = maxSize))[, 1:7]
        data$scgsea$es[df$pathway, i] = df$ES
        data$scgsea$nes[df$pathway, i] = df$NES
        data$scgsea$pval[df$pathway, i] = df$pval
        data$scgsea$fdr[df$pathway, i] = df$padj
        utils::setTxtProgressBar(pb, i)
    }
    base::close(pb)

    ix = is.na(data$scgsea$nes)
    if (sum(ix) > 0) {
        data$scgsea$nes[ix] = 0
        data$scgsea$pval[ix] = 1
        data$scgsea$fdr[ix] = 1
    }
    data$scgsea$x = data$scgsea$nes
    data$scgsea$x[data$scgsea$x < 0 | data$scgsea$fdr >= fdr.th] = 0
    data$scgsea$x = Matrix::Matrix(data = data$scgsea$nmf.h %*% 
        t(data$scgsea$x), sparse = T)
    data$scgsea$stat = df[, c("pathway", "size")]
    data$scgsea$x = data$scgsea$x[, colSums(as.matrix(data$scgsea$x)) > 0]
    oo = as.data.frame(data$scgsea$x)
    colnames(oo)=c("scGSEA")
    oo$Cell <- rownames(sr@meta.data)

write.table(oo, file = paste0("scGSEA/",args$tag,"_scGSEA.txt"), sep = '\t', row.names = F, col.names = TRUE, quote = FALSE)
