#' Title
#'
#' @param filelist
#' @param labelList
#' @param savePath
#'
#' @return
#' @export
#'
#' @examples
multiCorrNet <- function(filelist, labelList = NULL, savePath){
  labellist <- list()
  for(i in 1:length(filelist)){
    filelist[i] <- normalizePath(filelist[i], "/")
    filelist[i] <- paste0(filelist[i],"results.txt")
    tmlist <- strsplit(filelist[i],split = '/')[[1]]
    labellist <- c(labellist,list(tmlist[length(tmlist)-1]))
  }
  if(!is.null(labelList))
  labellist <- labelList
  GenePair <- data.frame()
  for(i in 1:length(filelist)){
    gpdata <- read.table(file = filelist[i],sep = '\t',header = T)
    tmlist1 <- paste0(gpdata$gene1,'-',gpdata$gene2)
    tmlist2 <- paste0(gpdata$gene2,'-',gpdata$gene1)
    tmlist <- c(tmlist1,tmlist2)
    tmlist <- as.data.frame(tmlist)
    tmlist <- rename(tmlist,c(tmlist = labellist[i]))
    GenePair <- plyr::rbind.fill(as.data.frame(t(GenePair)),as.data.frame(t(tmlist)))
    GenePair <- as.data.frame(t(GenePair))
  }
  colnames(GenePair) <- labellist
  #UniFeature <- Reduce(union, GenePair)
  IntFeature <- Reduce(intersect,GenePair)
  Uni_Gene_df <- separate(as.data.frame(IntFeature),IntFeature,c("Gene1","Gene2"),sep = '-')
  #------------save results------------
  write.table(Uni_Gene_df,paste0(savePath,'Uni_Gene_results.txt'),quote = F,sep = '\t',row.names = F,col.names = T)

  return(Uni_Gene_df)
}
