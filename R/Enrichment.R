#' Title
#'
#' @param genePairResults
#'
#' @return
#' @export
#'
#' @examples
Enrichment <- function(genePairResults){
  gl <- c(results[,1],results[,2])
  gene <- bitr(gl, fromType = "SYMBOL",toType = c("ENTREZID","ENSEMBL"),OrgDb = "org.Hs.eg.db")
  gene <- gene$ENTREZID
  gene[duplicated(gene)]
  go <- enrichGO(gene, OrgDb = org.Hs.eg.db, ont = 'ALL',pAdjustMethod = 'BH',pvalueCutoff = 0.05,qvalueCutoff = 0.2,keyType = 'ENTREZID')
  p <- dotplot(go,showCategory = 50)
  return(p)
}
