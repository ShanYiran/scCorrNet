#' Title
#'
#' @param genePairResults
#'
#' @return
#' @export
#'
#' @examples
Enrichment <- function(genePairResults,savePath){
  savePath <- normalizePath(savePath, "/")
  if(!dir.exists(file.path(savePath, 'report-figures/'))){
    dir.create(file.path(savePath, 'report-figures/'), recursive = T)
  }
  p.results <- list()
  gl <- c(genePairResults[,1],genePairResults[,2])
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


