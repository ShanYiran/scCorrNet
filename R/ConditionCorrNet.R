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
ConditionCorrNet <- function(exprMatrix,geneList,sampleNum = 10, sampleSize = 100,disMethod = 'spearman', pairNum = 10000,orderScore = 'dis',IfSaveFile = F, saveFile = 'none'){
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



