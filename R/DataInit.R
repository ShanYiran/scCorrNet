#' Title
#'
#' @param exprMatrix
#' @param Genelistfile
#' @param IfSaveFile
#' @param SaveFileDir
#'
#' @return
#' @export
#'
#' @examples
#'
TenxDataInit <- function(exprMatrix, Genelistfile = 'none',
                         IfSaveFile = F, SaveFileDir = 'none'){
  if(Genelistfile=='none'){
    genelist <- read.table(system.file("txt", "Gene_list_all_with_init.txt", package = "scCorrNet"))
  }
  else genelist <- read.table(Genelistfile)
  genelist <- genelist[[1]]
  exprMatrix <- exprMatrix[rownames(exprMatrix) %in% genelist,]
  if(IfSaveFile == TRUE){
    write.csv(exprMatrix,file = SaveFileDir,quote = F,col.names = T,row.names = T)
  }
  return(list(exprMatrix = exprMatrix, genelist = genelist))
}

#' Title
#'
#' @param exprFile
#' @param sampleFile
#' @param CellType
#' @param Genelistfile
#' @param IfSaveFile
#' @param SaveFileDir
#'
#' @return
#' @export
#'
#' @examples
BulkDataInit <- function(exprFile, sampleFile, CellType = 'HCC', Genelistfile = 'none', IfSaveFile = F, SaveFileDir = 'none'){
  sample <- read.table(sampleFile,sep = '\t',stringsAsFactors = F)
  expr <- read.table(exprFile,sep = '\t',stringsAsFactors = F)
  message("[", Sys.time(), "] -----: " ,dim(expr)[1],' ',dim(expr)[2])
  sample <- sample[1:2,sample[2,]==CellType]
  if(Genelistfile=='none'){
    genelist <- read.table(system.file("txt", "Gene_list_all_with_init.txt", package = "scCorrNet"))
  }
  else genelist <- read.table(Genelistfile)
  genelist <- genelist[[1]]
  message("[", Sys.time(), "] -----: " ,genelist[3])
  tmcolname <- c(sample[1,],'Symbol')
  tmrowname <- c(genelist,'Symbol')
  expr <- expr[expr[,2]%in%tmrowname,expr[1,]%in%tmcolname]
  exprMatrix <- as.matrix(expr[2:dim(expr)[1],2:dim(expr)[2]])
  colnames(exprMatrix) <- expr[1,2:dim(expr)[2]]
  rownames(exprMatrix) <- expr[2:dim(expr)[1],1]
  message("[", Sys.time(), "] -----: " ,dim(exprMatrix)[1],' ',dim(exprMatrix)[2])
  message("[", Sys.time(), "] -----: " ,dim(expr)[1],' ',dim(expr)[2])
  if(IfSaveFile == TRUE){
    write.csv(exprMatrix,file = SaveFileDir,quote = F,col.names = T,row.names = T)
  }
  return(list(exprMatrix = exprMatrix, genelist = genelist))
}


