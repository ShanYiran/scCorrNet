options(stringsAsFactors = F )
rm(list = ls())
library(Seurat)
library(dplyr)
library(rlist)
library(ggplot2)
library(Hmisc)
library(clusterProfiler)
library(scCorrNet)
library(knitr)
library(cowplot)
library(diptest)
library(markdown)
library(reshape)
library(tidyr)
library(plyr)
#PA_data_Pipline_without_immune_cells

tm <- load("D:/Gu_lab/PA/data/immune_cells.rda")

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


indir <- 'D:/Gu_lab/PA/data/P4T_1122/'
outdir <- 'D:/Gu_lab/PA/data/P4T_1122/scCorrNet/'
P4T_1122 <- load(paste0(indir,'P4T_1122_standard.rda'))
res_P4 <- res
data_P4 <- res$obj
data_P4 <- RenameCells(data_P4,add.cell.id = 'P4T')
cell_list <- setdiff(colnames(data_P4),immune_cells)
data_P4 <- SubsetData(data_P4,cells = cell_list)
data_P4 <- FindVariableFeatures(data_P4, selection.method = "vst", nfeatures = 2000)
data_P4 <- ScaleData(data_P4, features = VariableFeatures(data_P4))
em_P4 <- data_P4@assays$RNA@scale.data
expr_P4 = TenxDataInit(exprMatrix = em_P4 , IfSaveFile = T, SaveFileDir = paste0(outdir,'expr.csv'))
results <-runScCorrNet(exprMatrix = expr_P4$exprMatrix,
                       geneList = expr_P4$genelist,
                       sampleNum = 10,
                       sampleSize = 100,
                       P_value_thre = 0.05,
                       Corr_dis_thre = 0.25,
                       disMethod = 'spearman',
                       pairNum = 50000 ,
                       orderScore = 'dis',
                       savePath = outdir)


indir <- 'D:/Gu_lab/PA/data/P5T_1122/'
outdir <- 'D:/Gu_lab/PA/data/P5T_1122/scCorrNet/'
P5T_1122 <- load(paste0(indir,'P5T_1122_standard.rda'))
res_P5 <- res
data_P5 <- res$obj
data_P5 <- RenameCells(data_P5,add.cell.id = 'P5T')
cell_list <- setdiff(colnames(data_P5),immune_cells)
data_P5 <- SubsetData(data_P5,cells = cell_list)
data_P5 <- FindVariableFeatures(data_P5, selection.method = "vst", nfeatures = 2000)
data_P5 <- ScaleData(data_P5, features = VariableFeatures(data_P5))
em_P5 <- data_P5@assays$RNA@scale.data
expr_P5 = TenxDataInit(exprMatrix = em_P5 , IfSaveFile = T, SaveFileDir = paste0(outdir,'expr.csv'))
results <-runScCorrNet(exprMatrix = expr_P5$exprMatrix,
                       geneList = expr_P5$genelist,
                       sampleNum = 10,
                       sampleSize = 100,
                       P_value_thre = 0.05,
                       Corr_dis_thre = 0.25,
                       disMethod = 'spearman',
                       pairNum = 50000 ,
                       orderScore = 'dis',
                       savePath = outdir)


indir <- 'D:/Gu_lab/PA/data/P6T_1122/'
outdir <- 'D:/Gu_lab/PA/data/P6T_1122/scCorrNet/'
P6T_1122 <- load(paste0(indir,'P6T_1122_standard.rda'))
res_P6 <- res
data_P6 <- res$obj
data_P6 <- RenameCells(data_P6,add.cell.id = 'P6T')
cell_list <- setdiff(colnames(data_P6),immune_cells)
data_P6 <- SubsetData(data_P6,cells = cell_list)
data_P6 <- FindVariableFeatures(data_P6, selection.method = "vst", nfeatures = 2000)
data_P6 <- ScaleData(data_P6, features = VariableFeatures(data_P6))
em_P6 <- data_P6@assays$RNA@scale.data
expr_P6 = TenxDataInit(exprMatrix = em_P6 , IfSaveFile = T, SaveFileDir = paste0(outdir,'expr.csv'))
results <-runScCorrNet(exprMatrix = expr_P6$exprMatrix,
                       geneList = expr_P6$genelist,
                       sampleNum = 10,
                       sampleSize = 100,
                       P_value_thre = 0.05,
                       Corr_dis_thre = 0.25,
                       disMethod = 'spearman',
                       pairNum = 50000 ,
                       orderScore = 'dis',
                       savePath = outdir)

indir <- 'D:/Gu_lab/PA/data/P6T2_1122/'
outdir <- 'D:/Gu_lab/PA/data/P6T2_1122/scCorrNet/'
P6T2_1122 <- load(paste0(indir,'P6T2_1122_standard.rda'))
res_P6_2 <- res
data_P6_2 <- res$obj
data_P6_2 <- RenameCells(data_P6_2,add.cell.id = 'P6T2')
cell_list <- setdiff(colnames(data_P6_2),immune_cells)
data_P6_2 <- SubsetData(data_P6_2,cells = cell_list)
data_P6_2 <- FindVariableFeatures(data_P6_2, selection.method = "vst", nfeatures = 2000)
data_P6_2 <- ScaleData(data_P6_2, features = VariableFeatures(data_P6_2))
em_P6_2 <- data_P6_2@assays$RNA@scale.data
expr_P6 = TenxDataInit(exprMatrix = em_P6_2 , IfSaveFile = T, SaveFileDir = paste0(outdir,'expr.csv'))
results <-runScCorrNet(exprMatrix = expr_P6$exprMatrix,
                       geneList = expr_P6$genelist,
                       sampleNum = 10,
                       sampleSize = 100,
                       P_value_thre = 0.05,
                       Corr_dis_thre = 0.25,
                       disMethod = 'spearman',
                       pairNum = 50000 ,
                       orderScore = 'dis',
                       savePath = outdir)


indir <- 'D:/Gu_lab/PA/data/P7T_1122/'
outdir <- 'D:/Gu_lab/PA/data/P7T_1122/scCorrNet/'
P7T_1122 <- load(paste0(indir,'P7T_1122_standard.rda'))
res_P7 <- res
data_P7 <- res$obj
data_P7 <- RenameCells(data_P7,add.cell.id = 'P7T')
cell_list <- setdiff(colnames(data_P7),immune_cells)
data_P7 <- SubsetData(data_P7,cells = cell_list)
data_P7 <- FindVariableFeatures(data_P7, selection.method = "vst", nfeatures = 2000)
data_P7 <- ScaleData(data_P7, features = VariableFeatures(data_P7))
em_P7 <- data_P7@assays$RNA@scale.data
expr_P7 = TenxDataInit(exprMatrix = em_P7 , IfSaveFile = T, SaveFileDir = paste0(outdir,'expr.csv'))
results <-runScCorrNet(exprMatrix = expr_P7$exprMatrix,
                       geneList = expr_P7$genelist,
                       sampleNum = 10,
                       sampleSize = 100,
                       P_value_thre = 0.05,
                       Corr_dis_thre = 0.25,
                       disMethod = 'spearman',
                       pairNum = 50000 ,
                       orderScore = 'dis',
                       savePath = outdir)

indir <- 'D:/Gu_lab/PA/data/P8T_1122/'
outdir <- 'D:/Gu_lab/PA/data/P8T_1122/scCorrNet/'
P8T_1122 <- load(paste0(indir,'P8T_1122_standard.rda'))
res_P8 <- res
data_P8 <- res$obj
data_P8 <- RenameCells(data_P8,add.cell.id = 'P8T')
cell_list <- setdiff(colnames(data_P8),immune_cells)
data_P8 <- SubsetData(data_P8,cells = cell_list)
data_P8 <- FindVariableFeatures(data_P8, selection.method = "vst", nfeatures = 2000)
data_P8 <- ScaleData(data_P8, features = VariableFeatures(data_P8))
em_P8 <- data_P8@assays$RNA@scale.data
expr_P8 = TenxDataInit(exprMatrix = em_P8 , IfSaveFile = T, SaveFileDir = paste0(outdir,'expr.csv'))
results <-runScCorrNet(exprMatrix = expr_P8$exprMatrix,
                       geneList = expr_P8$genelist,
                       sampleNum = 10,
                       sampleSize = 100,
                       P_value_thre = 0.05,
                       Corr_dis_thre = 0.25,
                       disMethod = 'spearman',
                       pairNum = 50000 ,
                       orderScore = 'dis',
                       savePath = outdir)

indir <- 'D:/Gu_lab/PA/data/P12T_1122/'
outdir <- 'D:/Gu_lab/PA/data/P12T_1122/scCorrNet/'
P12T_1122 <- load(paste0(indir,'P12T_1122_standard.rda'))
res_P12 <- res
data_P12 <- res$obj
data_P12 <- RenameCells(data_P12,add.cell.id = 'P12T')
cell_list <- setdiff(colnames(data_P12),immune_cells)
data_P12 <- SubsetData(data_P12,cells = cell_list)
data_P12 <- FindVariableFeatures(data_P12, selection.method = "vst", nfeatures = 2000)
data_P12 <- ScaleData(data_P12, features = VariableFeatures(data_P12))
em_P12 <- data_P12@assays$RNA@scale.data
expr_P12 = TenxDataInit(exprMatrix = em_P12 , IfSaveFile = T, SaveFileDir = paste0(outdir,'expr.csv'))
results <-runScCorrNet(exprMatrix = expr_P12$exprMatrix,
                       geneList = expr_P12$genelist,
                       sampleNum = 10,
                       sampleSize = 100,
                       P_value_thre = 0.05,
                       Corr_dis_thre = 0.25,
                       disMethod = 'spearman',
                       pairNum = 50000 ,
                       orderScore = 'dis',
                       savePath = outdir)

indir <- 'D:/Gu_lab/PA/data/P21T_1122/'
outdir <- 'D:/Gu_lab/PA/data/P21T_1122/scCorrNet/'
P21T_1122 <- load(paste0(indir,'P21T_1122_standard.rda'))
res_P21 <- res
data_P21 <- res$obj
data_P21 <- RenameCells(data_P21,add.cell.id = 'P21T')
cell_list <- setdiff(colnames(data_P21),immune_cells)
data_P21 <- SubsetData(data_P21,cells = cell_list)
data_P21 <- FindVariableFeatures(data_P21, selection.method = "vst", nfeatures = 2000)
data_P21 <- ScaleData(data_P21, features = VariableFeatures(data_P21))
em_P21 <- data_P21@assays$RNA@scale.data
expr_P21 = TenxDataInit(exprMatrix = em_P21 , IfSaveFile = T, SaveFileDir = paste0(outdir,'expr.csv'))
results <-runScCorrNet(exprMatrix = expr_P21$exprMatrix,
                       geneList = expr_P21$genelist,
                       sampleNum = 10,
                       sampleSize = 100,
                       P_value_thre = 0.05,
                       Corr_dis_thre = 0.25,
                       disMethod = 'spearman',
                       pairNum = 50000 ,
                       orderScore = 'dis',
                       savePath = outdir)

indir <- 'D:/Gu_lab/PA/data/P22T_1122/'
outdir <- 'D:/Gu_lab/PA/data/P22T_1122/scCorrNet/'
P22T_1122 <- load(paste0(indir,'P22T_1122_standard.rda'))
res_P22 <- res
data_P22 <- res$obj
data_P22 <- RenameCells(data_P22,add.cell.id = 'P22T')
cell_list <- setdiff(colnames(data_P22),immune_cells)
data_P22 <- SubsetData(data_P22,cells = cell_list)
data_P22 <- FindVariableFeatures(data_P22, selection.method = "vst", nfeatures = 2000)
data_P22 <- ScaleData(data_P22, features = VariableFeatures(data_P22))
em_P22 <- data_P22@assays$RNA@scale.data
expr_P22 = TenxDataInit(exprMatrix = em_P22 , IfSaveFile = T, SaveFileDir = paste0(outdir,'expr.csv'))
results <-runScCorrNet(exprMatrix = expr_P22$exprMatrix,
                       geneList = expr_P22$genelist,
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

