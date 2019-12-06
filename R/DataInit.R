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
  if(Genelistfile=='auto'){
    genelist <- read.table(system.file("txt", "Gene_list_all_with_init.txt", package = "scCorrNet"))
  }
  else if(Genelistfile=='none'){
    genelist <- c(list(rownames(exprMatrix)))
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
  if(Genelistfile=='auto'){
    genelist <- read.table(system.file("txt", "Gene_list_all_with_init.txt", package = "scCorrNet"))
  }
  else if(Genelistfile=='none'){
    genelist <- c(list(expr[,2]))
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
  #use seurat to scale var gene
  seu_expr <- CreateSeuratObject(exprMatrix)
  seu_expr <- NormalizeData(seu_expr)
  seu_expr <- FindVariableFeatures(seu_expr, selection.method = "vst", nfeatures = 2000)
  seu_expr <- ScaleData(seu_expr,features = VariableFeatures(seu_expr))
  expr_tenX = TenxDataInit(exprMatrix = seu_expr@assays$RNA@scale.data , IfSaveFile = T, SaveFileDir = paste0(outdir,'expr.csv'))
  exprMatrix<- expr_tenX$exprMatrix
  return(list(exprMatrix = exprMatrix, genelist = genelist))
}


