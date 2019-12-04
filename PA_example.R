#PA_data_Pipline_without_immune_cells
options(stringsAsFactors = F)

tm <- load("D:/Gu_lab/PA/data/immune_cells.rda")
immune_cells <- tm

indir <- 'D:/Gu_lab/PA/data/P1T_1122/'
outdir <- 'D:/Gu_lab/PA/data/P1T_1122/scCorrNet/'
P1T_1122 <- load(paste0(indir,'P1T_1122_standard.rda'))
res_P1 <- res
data_P1 <- res$obj
data_P1 <- RenameCells(data_P1,add.cell.id = 'P1T')
cell_list <- setdiff(colnames(data_P1),immune_cells)
data_P1 <- SubsetData(data_P1,cells = cell_list)
data_P1 <- FindVariableFeatures(data_P1, selection.method = "vst", nfeatures = 2000)
data_P1 <- ScaleData(data_P1, features = VariableFeatures(data_P1))
em_P1 <- data_P1@assays$RNA@scale.data
expr_P1 = TenxDataInit(exprMatrix = em_P1 , IfSaveFile = T, SaveFileDir = paste0(outdir,'expr.csv'))
results <-runScCorrNet(exprMatrix = expr_P1$exprMatrix,
                       geneList = expr_P1$genelist,
                       sampleNum = 10,
                       sampleSize = 100,
                       P_value_thre = 0.05,
                       Corr_dis_thre = 0.25,
                       disMethod = 'spearman',
                       pairNum = 50000 ,
                       orderScore = 'dis',
                       savePath = outdir)


indir <- 'D:/Gu_lab/PA/data/P2T_1122/'
outdir <- 'D:/Gu_lab/PA/data/P2T_1122/scCorrNet/'
P2T_1122 <- load(paste0(indir,'P2T_1122_standard.rda'))
res_P2 <- res
data_P2 <- res$obj
data_P2 <- RenameCells(data_P2,add.cell.id = 'P2T')
cell_list <- setdiff(colnames(data_P2),immune_cells)
data_P2 <- SubsetData(data_P2,cells = cell_list)
data_P2 <- FindVariableFeatures(data_P2, selection.method = "vst", nfeatures = 2000)
data_P2 <- ScaleData(data_P2, features = VariableFeatures(data_P2))
em_P2 <- data_P2@assays$RNA@scale.data
expr_P2 = TenxDataInit(exprMatrix = em_P2 , IfSaveFile = T, SaveFileDir = paste0(outdir,'expr.csv'))
results <-runScCorrNet(exprMatrix = expr_P2$exprMatrix,
                       geneList = expr_P2$genelist,
                       sampleNum = 10,
                       sampleSize = 100,
                       P_value_thre = 0.05,
                       Corr_dis_thre = 0.25,
                       disMethod = 'spearman',
                       pairNum = 50000 ,
                       orderScore = 'dis',
                       savePath = outdir)



indir <- 'D:/Gu_lab/PA/data/P3T_1122/'
outdir <- 'D:/Gu_lab/PA/data/P3T_1122/scCorrNet/'
P3T_1122 <- load(paste0(indir,'P3T_1122_standard.rda'))
res_P3 <- res
data_P3 <- res$obj
data_P3 <- RenameCells(data_P3,add.cell.id = 'P3T')
cell_list <- setdiff(colnames(data_P3),immune_cells)
data_P3 <- SubsetData(data_P3,cells = cell_list)
data_P3 <- FindVariableFeatures(data_P3, selection.method = "vst", nfeatures = 2000)
data_P3 <- ScaleData(data_P3, features = VariableFeatures(data_P3))
em_P3 <- data_P3@assays$RNA@scale.data
expr_P3 = TenxDataInit(exprMatrix = em_P3 , IfSaveFile = T, SaveFileDir = paste0(outdir,'expr.csv'))
results <-runScCorrNet(exprMatrix = expr_P3$exprMatrix,
                       geneList = expr_P3$genelist,
                       sampleNum = 10,
                       sampleSize = 100,
                       P_value_thre = 0.05,
                       Corr_dis_thre = 0.25,
                       disMethod = 'spearman',
                       pairNum = 50000 ,
                       orderScore = 'dis',
                       savePath = outdir)









#多样本




tm <- multiCorrNet(filelist = c("D:/Gu_lab/PA/data/P1T_1122/scCorrNet/",
                                "D:/Gu_lab/PA/data/P2T_1122/scCorrNet/",
                                "D:/Gu_lab/PA/data/P3T_1122/scCorrNet/",
                                "D:/Gu_lab/PA/data/P4T_1122/scCorrNet/",
                                "D:/Gu_lab/PA/data/P5T_1122/scCorrNet/",
                                "D:/Gu_lab/PA/data/P6T_1122/scCorrNet/"),savePath = "D:/Gu_lab/PA/data/scCorrNet/")

length(tm[,1])
tm2 <- multiCorrNet(filelist = c("D:/Gu_lab/PA/data/P1T_1122/scCorrNet/","D:/Gu_lab/PA/data/P2T_1122/scCorrNet/"),savePath = "D:/Gu_lab/PA/data/scCorrNet/")
tm3 <- multiCorrNet(filelist = c("D:/Gu_lab/PA/data/P3T_1122/scCorrNet/","D:/Gu_lab/PA/data/P2T_1122/scCorrNet/"),savePath = "D:/Gu_lab/PA/data/scCorrNet/")
tm3 <- multiCorrNet(filelist = c("D:/Gu_lab/PA/data/P3T_1122/scCorrNet/","D:/Gu_lab/PA/data/P1T_1122/scCorrNet/"),savePath = "D:/Gu_lab/PA/data/scCorrNet/")
tm4 <- multiCorrNet(filelist = c("D:/Gu_lab/PA/data/P3T_1122/scCorrNet/","D:/Gu_lab/PA/data/P1T_1122/scCorrNet/","D:/Gu_lab/PA/data/P2T_1122/scCorrNet/"),savePath = "D:/Gu_lab/PA/data/scCorrNet/")


genePairResults <- c(tm$Gene1,tm$Gene2)
gene <- bitr(genePairResults, fromType = "SYMBOL",toType = c("ENTREZID","ENSEMBL"),OrgDb = "org.Hs.eg.db")
gene <- gene$ENTREZID
gene[duplicated(gene)]
go <- enrichGO(gene, OrgDb = org.Hs.eg.db, ont = 'ALL',pAdjustMethod = 'BH',pvalueCutoff = 0.05,qvalueCutoff = 0.2,keyType = 'ENTREZID')
dotplot(go,showCategory = 50)
go@ontology
kegg <- enrichKEGG(gene,organism = 'hsa', keyType = 'kegg', pvalueCutoff = 0.05,pAdjustMethod = 'BH', minGSSize = 10,maxGSSize = 500,qvalueCutoff = 0.2,use_internal_data = FALSE)
barplot(kegg,showCategory=50,drop=T)
kegg@result$Description


#overlap var Gene
OVG <- Reduce(intersect,c(list(rownames(expr_P1$exprMatrix)),list(rownames(expr_P2$exprMatrix)),list(rownames(expr_P3$exprMatrix)),list(rownames(expr_P4$exprMatrix)),list(rownames(expr_P5$exprMatrix)),list(rownames(expr_P6$exprMatrix))))
gene <- bitr(OVG, fromType = "SYMBOL",toType = c("ENTREZID","ENSEMBL"),OrgDb = "org.Hs.eg.db")
gene <- gene$ENTREZID
gene[duplicated(gene)]
go <- enrichGO(gene, OrgDb = org.Hs.eg.db, ont = 'ALL',pAdjustMethod = 'BH',pvalueCutoff = 0.05,qvalueCutoff = 0.2,keyType = 'ENTREZID')
dotplot(go,showCategory = 50)
go@ontology
kegg <- enrichKEGG(gene,organism = 'hsa', keyType = 'kegg', pvalueCutoff = 0.05,pAdjustMethod = 'BH', minGSSize = 10,maxGSSize = 500,qvalueCutoff = 0.2,use_internal_data = FALSE)
barplot(kegg,showCategory=50,drop=T)
kegg@result$Description

