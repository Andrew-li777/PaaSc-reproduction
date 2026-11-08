library(BiocParallel)
library(GSVA)
library(Seurat)
library(dplyr)
library(ggrepel)
library(ggplot2)
library(gsdensity)
library(CelliD)
library(VAM)
library(DescTools)
library(GSEABase)
library(AUCell)
library(nnet)
library(reshape2)
library(RColorBrewer)
library(PaaSc)


##https://github.com/qingnanl/GSDensity_manuscript_code/


##load_data

preprocess_sr <- function(sr.fname){
  sr <- readRDS(sr.fname)
  sr <- sr %>%
    NormalizeData() %>% 
    FindVariableFeatures() %>%
    ScaleData() %>%
    RunPCA() %>%
    FindNeighbors() %>%
    RunUMAP(dims = 1:10)
  return(sr)
}


#####################


  fname="pbmc3k"
  sr.fname <- paste0(fname, ".rds")
  sr.markers.fname <- paste0(fname, ".markers.rds")
  # obtain seurat object
  sr <- preprocess_sr(sr.fname)
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
  

marker = readRDS("pbmc3k.markers.rds")
clist = ls(marker)

igraph::igraph.options(sparsematrices = FALSE)
# get a function for calculate AUC
# get a function for calculate AUC
calc.auc <- function(input, 
                     pos.id, 
                     method.id,
                     sample.points = c(seq(from = 0, to = 0.2, by = 0.025), seq(from = 0.3, to = 1, by = 0.1))){
  total.pos <- sum(input[, 2] == pos.id)
  input.order <- input[order(input[, 1], decreasing = T), ]
  output.df <- as.data.frame(sample.points)
  output.df$positive.capture <- 0
  output.df$method <- method.id
  for (i in 2:length(sample.points)){
    pos.count <- sum(input.order[1:floor(nrow(input.order)*sample.points[i]), 2] == pos.id)
    output.df[i, "positive.capture"] <- pos.count/total.pos
  }
  return(output.df)
}


# get cell-level activity scores
# gsd
get.activity.list <- function(glist, pos.id, object, el, cells, cells_rankings){
  lst <- list()
  print(glist)
  gsd <- run.rwr(el = el, gene_set = glist, cells = cells)
  gsd.res <- as.data.frame(gsd)
  gsd.res$cell.type <- object$cell.type
  lst[['gsdensity']] <- gsd.res
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
  gsea.res$cell.type <- object$cell.type
  lst[['gsea']] <- gsea.res
  
  b <- GSVA:::.gsva(expr = object@assays$RNA@data, gset.idx.list = mapdgenesets, method = "gsva", rnaseq = F, 
                    kernel = F, parallel.sz = 40, BPPARAM = SnowParam(), verbose = T)
  gsva.res <- as.data.frame(t(b))
  rownames(gsva.res) <- colnames(object)
  colnames(gsva.res) <- 'gsva'
  gsva.res$cell.type <- object$cell.type
  lst[['gsva']] <- gsva.res
  # cellid
  
  HGT_glist <- RunCellHGT(object, pathways = list(glist = glist), dims = 1:30, minSize = 5)
  cellid.res <- as.data.frame(HGT_glist)
  colnames(cellid.res) <- 'cellid'
  cellid.res$cell.type <- object$cell.type
  lst[['cellid']] <- cellid.res
  #
  
  geneSets <- list(geneSet1 = glist)
  cells_AUC <- AUCell_calcAUC(geneSets, rankings = cells_rankings)
  cells_AUC_df <- as.data.frame(t(cells_AUC@assays@data@listData$AUC))
  aucell.res <- as.data.frame(cells_AUC_df)
  colnames(aucell.res) <- 'aucell'
  aucell.res$cell.type <- object$cell.type
  lst[['aucell']] <- aucell.res
  
  coll <- list()
  coll[[1]] <- intersect(rownames(object), glist)
  names(coll)[1] <- "VAM"
  gene.set.collection = createGeneSetCollection(gene.ids=rownames(object),
                                                gene.set.collection=coll)
  obj <- vamForSeurat(seurat.data=object,
                      gene.set.collection=gene.set.collection,
                      center=F, gamma=T, sample.cov=F, return.dist=F)
  vam.res <- as.data.frame(t(obj@assays$VAMcdf@data))
  vam.res$cell.type <- object$cell.type
  lst[['vam']] <- vam.res
  

  pathway_geneset = list(B_cell_activation = glist)
  gene_rate2 = getGeneRate(background.file ="Tres.kegg", pathway.geneset = pathway_geneset)
  regression_data2 = doRegression(object2, gene.rate = gene_rate2)
  paasc.res = computeScore(object2, regression.data = regression_data2)
  colnames(paasc.res) <- 'paasc'
  paasc.res$cell.type <- object$cell.type
  lst[['paasc']] <- paasc.res


  lst1 <- list()
  for (method in names(lst)){
    lst1[[method]] <- calc.auc(input = lst[[method]], method.id = method, pos.id = pos.id)
  }
  return(list(score.list = lst, auc.list = lst1))
}



#############



object=sr
benchmark.AUC <- function(object, clist, methods, markers, el, cells, cells_rankings){
  AUC.collection <- list()
  for (i in 1:length(clist)){
    # for cells

    glist <- marker[[clist[i]]]
    # print(glist)
    ga <- get.activity.list(glist = glist, pos.id = clist[i], object = object, 
                            el = el, cells = cells, cells_rankings = cells_rankings)
    AUC.list <- c()
    for (method in names(ga$auc.list)){
      AUC.list <- c(AUC.list, AUC(x = ga$auc.list[[method]][, 'sample.points'],
                                  y = ga$auc.list[[method]][, 'positive.capture']))
     }
    AUC.df <- data.frame(AUC.score = AUC.list)
    AUC.df$methods <- c("gsdensity","gsea", "gsva","cellid", "aucell", "vam","paasc")
    AUC.df$cell.type <- clist[i]
    AUC.collection[[clist[i]]] <- AUC.df
  }
  AUC.col.plot <- do.call(rbind, AUC.collection)
  return(AUC.col.plot)
}



AUC_results <- benchmark.AUC(object = sr, 
                             clist = clist, 
                             methods = methods, 
                             markers = marker, 
                             el = el, 
                             cells = cells, 
                             cells_rankings = cells_rankings)





  
  marker.list <- marker
  clist <- names(marker.list)
  print(clist)
  n.cell.types <- length(clist)
  # AUC.tar <- list()
  AUC.batch <- list()
  AUC.res <- benchmark.AUC(object = sr, clist = clist, 
                           methods = c("gsdensity", "gsea", "gsva","cellid", "aucell", "vam","paasc"), 
                           markers = marker.list, el = el, cells = cells, cells_rankings = cells_rankings)
  AUC.res$design <- "full.marker.set"
  # AUC.res$drop.out <- drop.out
  AUC.res$batch <- 1
  # AUC.res$target.type <- tar.type.name
  AUC.batch[["full.marker.set"]] <- AUC.res


ds.glist <- function(glist, all.genes, ratio, fixed, mode = "fixed", bottom = 5, s){
  set.seed(s)
  glist <- intersect(glist, all.genes)
  if (mode == "fixed"){
    sample.size <- fixed
  }
  if (mode == "ratio"){
    sample.size <- max(floor(length(glist) * ratio), bottom)
  }
  idx <- sample(x = 1:length(glist), size = sample.size, replace = F)
  out.glist <- glist[idx]
  return(out.glist)
}


  k.range = 5
  for (k in 1:k.range){
    ds.markers <- list()
    for (q in 1:n.cell.types){
      ds.markers[[q]] <- ds.glist(glist = marker.list[[q]], all.genes = rownames(sr), mode = "ratio", ratio = 0.9, s = k)
    }
    names(ds.markers) <- clist
    
    AUC.res <- benchmark.AUC(object = sr, clist = clist, 
                             methods = c("gsdensity", "gsea", "gsva","cellid", "aucell", "vam","paasc"), 
                             markers = ds.markers, el = el, cells = cells, cells_rankings = cells_rankings)
    AUC.res$design <- "90_percent.ds"
    # AUC.res$drop.out <- drop.out
    AUC.res$batch <- k
    # AUC.res$target.type <- tar.type.name
    AUC.batch[[paste0("90_percent.ds.", k)]] <- AUC.res
  }



output.auc <- function(info){
  auc <- info
  designs <- unique(auc$design)
  auc.ls <- list()
  for (des in designs){
    dat <- auc[auc$design == des, ] %>%
      group_by(methods) %>%
      summarise_at(vars(AUC.score), list(name = mean)) %>% 
      as.data.frame()
    dat <- as.data.frame(dat$name)
    colnames(dat) <- des
    auc.ls[[des]] <- dat
  }
  out <- do.call(cbind, auc.ls)
  rownames(out) <- c("aucell", "cellid", "gsdensity", "gsea", "gsva", "paasc", "vam")
  return(out)
}


  Out = output.auc(pbmc)
  Out$v1 = rownames(Out)
  Out = Out[,c(10,1:9)]
  plt.df=melt(Out)


plt.df$variable <- factor(plt.df$variable, 
                          levels = c("full.marker.set", "90_percent.ds", "80_percent.ds", 
                                     "70_percent.ds", "60_percent.ds","50_percent.ds", "40_percent.ds",
                                      "30_percent.ds", "20_percent.ds"))


p=ggplot(plt.df,aes(x=variable,y=value,color=v1))+
geom_smooth(se=F,size=1,span=0.25)+
scale_y_continuous(expand=c(0,0))+
theme_classic()+
scale_color_manual(values = c('steelblue','black','#56A902','#2e409a','#fcbe32','#e97f02','#f15c5c'))+
labs(x="Methods",y="Mean AUC Score")+
    theme_classic()+
      theme(axis.text.x = element_text(color = "black", 
                                     angle = 45, hjust = 0, size = 15),
          axis.text.y = element_text(color = "black", size = 15),
          legend.title = element_text(size = 15), 
          legend.text = element_text(size = 15),
          panel.background = element_blank(), 
          panel.border = element_rect(colour = "black", fill = NA, linewidth = 1.2))+
      coord_cartesian(ylim = c(0, 1))

print(p)




###### scGSEA
library(Matrix)
library(fgsea)
library(gficf)

#' Calculate scGSEA activity scores
#'
#' @param glist Gene set list
#' @param pos.id Pathway identifier
#' @param tmp NMF result object containing w and h matrices
#' @param data Main data object containing rawCounts
#' @param sr Seurat object containing cell type metadata
#' @param nsim Number of permutations (default=10000)
#' @param fdr.th FDR threshold (default=0.05)
#' @param minSize Minimum gene set size (default=5)
#' @param maxSize Maximum gene set size (default=Inf)
#' @param gp GSEA parameter (default=0)
#' 
#' @return A data.frame containing scGSEA scores and cell types

calculate_scgsea <- function(glist, pos.id, tmp, data, sr, 
                             nsim = 10000, fdr.th = 0.05,
                             minSize = 5, maxSize = Inf, gp = 0) {
  # 初始化scgsea结构
  data$scgsea <- list(
    nmf.w = Matrix::Matrix(tmp$w, sparse = TRUE),
    nmf.h = t(Matrix::Matrix(tmp$h, sparse = TRUE)),
    pathways = setNames(list(glist), pos.id),
    es = NULL, nes = NULL, pval = NULL, fdr = NULL,
    x = NULL, stat = NULL
  )
  
  # 初始化结果矩阵
  n_pathways <- length(data$scgsea$pathways)
  n_components <- ncol(data$scgsea$nmf.w)
  
  data$scgsea$es <- data$scgsea$nes <- 
  data$scgsea$pval <- data$scgsea$fdr <- 
  Matrix::Matrix(0, n_pathways, n_components)
  
  rownames(data$scgsea$es) <- rownames(data$scgsea$nes) <-
  rownames(data$scgsea$pval) <- rownames(data$scgsea$fdr) <- pos.id
  
  # 设置并行计算
  nt_fgsea <- ceiling(n_pathways / 100)
  bpparameters <- BiocParallel::SnowParam(nt_fgsea)
  rownames(data$scgsea$nmf.w) <- rownames(data$rawCounts)
  
  # 运行GSEA分析
  pb <- utils::txtProgressBar(0, n_components, style = 3)
  for(i in seq_len(n_components)) {
    stats <- as.matrix(data$scgsea$nmf.w)[, i]
    df <- fgsea::fgseaMultilevel(
      pathways = data$scgsea$pathways,
      stats = stats,
      nPermSimple = nsim,
      gseaParam = gp,
      BPPARAM = bpparameters,
      minSize = minSize,
      maxSize = maxSize
    )[, 1:7]
    
    data$scgsea$es[df$pathway, i] <- df$ES
    data$scgsea$nes[df$pathway, i] <- df$NES
    data$scgsea$pval[df$pathway, i] <- df$pval
    data$scgsea$fdr[df$pathway, i] <- df$padj
    
    utils::setTxtProgressBar(pb, i)
  }
  close(pb)
  
  # 后处理结果
  handle_nas <- function(x) replace(x, is.na(x), 0)
  data$scgsea$nes <- handle_nas(data$scgsea$nes)
  data$scgsea$pval <- handle_nas(data$scgsea$pval)
  data$scgsea$fdr <- handle_nas(data$scgsea$fdr)
  
  # 计算最终得分
  data$scgsea$x <- data$scgsea$nes
  data$scgsea$x[data$scgsea$x < 0 | data$scgsea$fdr >= fdr.th] <- 0
  data$scgsea$x <- Matrix::Matrix(
  data$scgsea$nmf.h %*% t(data$scgsea$x), 
    sparse = TRUE
  )
  
  # 生成输出数据框
  oo <- as.data.frame(data$scgsea$x)
  colnames(oo) <- "scGSEA"
  oo$cell.type <- sr@meta.data$cell.type
  
  return(oo)
}

get.activity.list <- function(glist, pos.id, object, el, cells, cells_rankings){
  lst <- list()
  print(glist)

  scgsea.res <- calculate_scgsea(
    glist = glist,
    pos.id = pos.id,
    tmp = tmp,         # 假设tmp已定义
    data = data,       # 假设data已定义
    sr = sr,           # 假设sr已定义
    nsim = 10000,
    fdr.th = 0.05
  )
  
  lst[['scGSEA']] <- scgsea.res


  lst1 <- list()
  for (method in names(lst)){
    lst1[[method]] <- calc.auc(input = lst[[method]], method.id = method, pos.id = pos.id)
  }
  return(list(score.list = lst, auc.list = lst1))
}

