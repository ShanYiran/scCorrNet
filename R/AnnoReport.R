#' Title
#'
#' @param exprMatrix
#' @param geneList
#' @param sampleNum
#' @param sampleSize
#' @param disMethod
#' @param pairNum
#' @param savePath
#' @param orderScore
#'
#' @return
#' @export
#'
#' @examples
ConditionCorrNet <- function(exprMatrix,
                             geneList,
                             sampleNum = 10,
                             sampleSize = 100,
                             disMethod = 'spearman',
                             pairNum = 'auto',
                             orderScore = 'dis',
                             savePath = 'none'){
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
  results <- as.data.frame(results)
   #------------plot p_value / dis_cor hist----------
  if(!dir.exists(file.path(savePath, 'report-figures/'))){
    dir.create(file.path(savePath, 'report-figures/'), recursive = T)
  }
  #hist plot
  p_thre <- 0.05
  p.results <- list()
  p.results[['hist_pvalue']] <- ggplot(data = results,aes(x = as.numeric(results$p_sorted_value))) +
    geom_histogram(bins = 200, fill = "#a788ab") +
    labs(x = "p_value", y = "Gene pair number")
  ggsave(filename = file.path(savePath, "report-figures/hist_pvalue.png"),p.results[['hist_pvalue']],
         width = 6, height = 6, dpi = 800)
  #hist plot
  p.results[['hist_disvalue']] <- ggplot(data = results,aes(x =as.numeric(results$dis_sorted_value) )) +
    geom_histogram(bins = 200, fill = "#a788ab") +
    labs(x = "dis_value", y = "Gene pair number")
  ggsave(filename = file.path(savePath, "report-figures/hist_disvalue.png"), p.results[['hist_disvalue']],
         width = 6, height = 6, dpi = 800)
  return(list(results = results,p.results = p.results))
}
#' Title
#'
#' @param genePairResults
#' @param savePath
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
  gl <- c(genePairResults$gene1,genePairResults$gene2)
  gene <- bitr(gl, fromType = "SYMBOL",toType = c("ENTREZID","ENSEMBL"),OrgDb = "org.Hs.eg.db")
  gene <- gene$ENTREZID
  gene[duplicated(gene)]
  go <- enrichGO(gene, OrgDb = org.Hs.eg.db, ont = 'ALL',pAdjustMethod = 'BH',pvalueCutoff = 0.05,qvalueCutoff = 0.2,keyType = 'ENTREZID')
  p.results[["go.dot"]] <- dotplot(go,showCategory = 50)
  ggsave(filename = file.path(savePath, "report-figures/go_dotPlot.png"), p.results[["go.dot"]],
         width = 8.5, height = 11, dpi = 800)
  p.results[["go.bar"]] <- barplot(go,showCategory=50,drop=T)
  ggsave(filename = file.path(savePath, "report-figures/go_barPlot.png"), p.results[["go.bar"]],
         width = 8.5, height = 11, dpi = 800)
  kegg <- enrichKEGG(gene,organism = 'hsa', keyType = 'kegg', pvalueCutoff = 0.05,pAdjustMethod = 'BH', minGSSize = 10,maxGSSize = 500,qvalueCutoff = 0.2,use_internal_data = FALSE)
  p.results[["kegg.dot"]] <- barplot(kegg,showCategory=50,drop=T)
  ggsave(filename = file.path(savePath, "report-figures/kegg_barPlot.png"), p.results[["kegg.dot"]],
         width = 8.5, height = 11, dpi = 800)
  return(p.results)
}

#' Title
#'
#' @param genePairResults
#' @param P_value_thre
#' @param Corr_dis_thre
#' @param savePath
#'
#' @return
#' @export
#'
#' @examples
ShowCorrScore <- function(genePairResults,P_value_thre = 0.05, Corr_dis_thre = 0.25, savePath){
  savePath <- normalizePath(savePath, "/")
  if(!dir.exists(file.path(savePath, 'report-figures/'))){
    dir.create(file.path(savePath, 'report-figures/'), recursive = T)
  }
  genePairResults$p_sorted_value <- as.numeric(genePairResults$p_sorted_value)
  genePairResults$dis_sorted_value <- as.numeric(genePairResults$dis_sorted_value)

  genePairResults$Is_P <- (genePairResults$p_sorted_value <= P_value_thre)
  genePairResults$Is_dis <- (genePairResults$dis_sorted_value >= Corr_dis_thre | genePairResults$dis_sorted_value <= -Corr_dis_thre)
  all_plot_cell <- (genePairResults$Is_P==TRUE & genePairResults$Is_dis==TRUE)
  p.results <- list()
  #print(head(genePairResults))
  all_plot_cell <- (genePairResults$Is_P==TRUE & genePairResults$Is_dis==TRUE)
  p.results[["CorrScore"]] <- ggplot(data = genePairResults,mapping = aes(x = dis_sorted_value,y = p_sorted_value,colour = all_plot_cell)) + geom_point()
  ggsave(plot =  p.results[["CorrScore"]],filename = file.path(savePath, "report-figures/ShowCorrScore.png"),
         width = 5, height = 7, dpi = 800)
  return(list(results = genePairResults,p.results = p.results))
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
#' @param P_value_thre
#' @param Corr_dis_thre
#'
#' @return
#' @export
#'
#' @examples
runScCorrNet <- function(exprMatrix,
                         geneList,
                         sampleNum = 10,
                         sampleSize = 100,
                         P_value_thre = 0.05,
                         Corr_dis_thre = 0.25,
                         disMethod = 'spearman',
                         pairNum = 'auto' ,
                         orderScore = 'dis',
                         savePath){
  message("[", Sys.time(), "] START: RUN ScCorrNet")
  results <- as.list(environment())
  results[["savePath"]] <- savePath
  #------------Condition Network-----------
  CCresults <- ConditionCorrNet(exprMatrix = exprMatrix,
                               geneList = geneList,
                               sampleNum = sampleNum,
                               sampleSize = sampleSize,
                               disMethod = disMethod,
                               pairNum = pairNum ,
                               orderScore = orderScore,
                               savePath = savePath)
  #print(head(CCresults))
  genePair <- CCresults$results
  results[["genePair"]] <- CCresults$p.results
  #------------Enrichment------------
  message("[", Sys.time(), "] START: RUN EnrichPlots")
  results[["EnrichPlots"]] <- Enrichment(genePairResults = genePair,
                          savePath = savePath)
  #------------plot Corr Score------------
  message("[", Sys.time(), "] START: RUN ShowCorrScore")
  tm <- ShowCorrScore(genePairResults = genePair,
                      P_value_thre = P_value_thre,
                      Corr_dis_thre = Corr_dis_thre,
                      savePath = savePath)

  genePair <- tm$results
  results[["CorrScorePlot"]] <- tm$p.results
  rm(tm)

  #------------save results------------
  write.table(genePair,paste0(savePath,'results.txt'),quote = F,sep = '\t',row.names = F,col.names = T)

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

