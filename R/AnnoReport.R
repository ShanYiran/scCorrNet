#' Title
#'
#' @param exprMatrix
#' @param geneList
#' @param sampleNum
#' @param sampleSize
#' @param disMethod
#' @param pairNum
#' @param orderScore
#' @param IfSaveFile
#' @param saveFile
#'
#' @return
#' @export
#'
#' @examples
ConditionCorrNet <- function(exprMatrix,geneList,sampleNum = 10, sampleSize = 100,disMethod = 'spearman', pairNum = 'auto' ,orderScore = 'dis',IfSaveFile = F, saveFile = 'none'){
  p_value <- matrix(data = 0,nrow = dim(exprMatrix)[1],ncol = dim(exprMatrix)[1])
  p_value[lower.tri(p_value)] = NA
  cor_dis <- matrix(data = 0,nrow = dim(exprMatrix)[1],ncol = dim(exprMatrix)[1])
  cor_dis[lower.tri(cor_dis)] = NA
  diag(cor_dis) = NA
  diag(p_value) = NA
  for(i in 1:sampleNum){
    samplecell <- sample(colnames(exprMatrix),size =sampleSize)
    SubExprMatrix <- exprMatrix[,colnames(exprMatrix)%in%samplecell]
    SubCorr <- rcorr(t(SubExprMatrix))
    SubCorr$P[lower.tri(SubCorr$P)] = 0
    SubCorr$r[lower.tri(SubCorr$r)] = 0
    p_value <- p_value + SubCorr$P
    cor_dis <- cor_dis + SubCorr$r
  }
  p_value <- p_value/sampleNum
  cor_dis <- cor_dis/sampleNum
  if(orderScore=='p_value'){
    orderGenePair <- arrayInd(sort.list(p_value,decreasing = F)[1:min(dim(exprMatrix)[1]*(dim(exprMatrix)[1]-1)/2, pairNum)],dim(p_value))
  }
  else if(orderScore=='dis'){
    orderGenePair <- arrayInd(sort.list(cor_dis,decreasing = T)[1:min(dim(exprMatrix)[1]*(dim(exprMatrix)[1]-1)/2, pairNum)],dim(cor_dis))
  }
  gene1 <- rownames(exprMatrix)[orderGenePair[,1]]
  gene2 <- rownames(exprMatrix)[orderGenePair[,2]]
  p_sorted_value <- p_value[orderGenePair]
  dis_sorted_value <- cor_dis[orderGenePair]
  results <- cbind(gene1,gene2,p_sorted_value,dis_sorted_value)
  if(IfSaveFile==T){
    write.table(results,paste0(saveFile,'results.txt'),quote = F,sep = '\t',row.names = F,col.names = T)
  }
  return(results)
}
#' Title
#'
#' @param genePairResults
#'
#' @return
#' @export
#'
#' @examples
Enrichment <- function(genePairResults,
                       savePath){
  savePath <- normalizePath(savePath, "/")
  if(!dir.exists(file.path(savePath, 'report-figures/'))){
    dir.create(file.path(savePath, 'report-figures/'), recursive = T)
  }
  p.results <- list()
  gl <- c(results[,1],results[,2])
  gene <- bitr(gl, fromType = "SYMBOL",toType = c("ENTREZID","ENSEMBL"),OrgDb = "org.Hs.eg.db")
  gene <- gene$ENTREZID
  gene[duplicated(gene)]
  go <- enrichGO(gene, OrgDb = org.Hs.eg.db, ont = 'ALL',pAdjustMethod = 'BH',pvalueCutoff = 0.05,qvalueCutoff = 0.2,keyType = 'ENTREZID')
  p.results[["go.dot"]] <- dotplot(go,showCategory = 50)
  ggsave(filename = file.path(savePath, "report-figures/go_dotPlot.png"), p.results[["go.dot"]],
         width = 6, height = 6, dpi = 800)
  p.results[["go.bar"]] <- barplot(go,showCategory=20,drop=T)
  ggsave(filename = file.path(savePath, "report-figures/go_barPlot.png"), p.results[["go.bar"]],
         width = 6, height = 6, dpi = 800)
  kegg <- enrichKEGG(gene,organism = 'hsa', keyType = 'kegg', pvalueCutoff = 0.05,pAdjustMethod = 'BH', minGSSize = 10,maxGSSize = 500,qvalueCutoff = 0.2,use_internal_data = FALSE)
  p.results[["kegg.dot"]] <- barplot(kegg,showCategory=20,drop=T)
  ggsave(filename = file.path(savePath, "report-figures/kegg_barPlot.png"), p.results[["kegg.dot"]],
         width = 6, height = 6, dpi = 800)
  return(p.results)
}


#' Title
#'
#' @param exprMatrix
#' @param exprMatrix
#' @param geneList
#' @param sampleNum
#' @param sampleSize
#' @param disMethod
#' @param pairNum
#' @param orderScore
#' @param savePath
#'
#' @return
#' @export
#'
#' @examples
runScCorrNet <- function(exprMatrix,
                         geneList,
                         sampleNum = 10,
                         sampleSize = 100,
                         disMethod = 'spearman',
                         pairNum = 'auto' ,
                         orderScore = 'dis',
                         savePath){
  message("[", Sys.time(), "] START: RUN ScCorrNet")
  results <- as.list(environment())
  results[["savePath"]] <- savePath
  #------------Condition Network-----------
  genePair <- ConditionCorrNet(exprMatrix = exprMatrix,
                               geneList = geneList,
                               sampleNum = sampleNum,
                               sampleSize = sampleSize,
                               disMethod = disMethod,
                               pairNum = pairNum ,
                               orderScore = orderScore,
                               IfSaveFile = T,
                               savePath)
  results[["genePair"]] <- genePair
  #------------Enrichment------------
  results[["EnrichPlots"]] <- Enrichment(genePairResults = genePair,
                          savePath = savePath)

  #------------report-------------
  message("[", Sys.time(), "] -----: report generating")
  if(!dir.exists(file.path(savePath, 'report-figures/'))){
    dir.create(file.path(savePath, 'report-figures/'), recursive = T)
  }
  suppressWarnings(
    knit(system.file("rmd", "main-scCorrNet.Rmd", package = "scCorrNet"),
         file.path(savePath,'report-scCorrNet.md'), quiet = T)
  )
  markdownToHTML(file.path(savePath,'report-scCorrNet.md'),
                 file.path(savePath, 'report-scCorrNet.html'))
  return(results)
}

