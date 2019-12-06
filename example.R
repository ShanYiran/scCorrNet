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

indir <- 'D:/Gu_lab/HCCsc/HCCsc_new/data/Hcc_bulk_DEGene/Multi_data/HCCDB-1/'
outdir <- 'D:/Gu_lab/HCCsc/HCCsc_new/result/scCorrNet/HCCDB1/'
expr = BulkDataInit(paste0(indir,'GSE22058-GPL6793.gene.txt'), paste0(indir,'GSE22058.sample.txt'), CellType = 'HCC', Genelistfile = 'none', IfSaveFile = T, SaveFileDir = paste0(outdir,'expr.csv'))
#test small data
exprMatrix<- expr_HCCDB1$exprMatrix
results <-runScCorrNet(exprMatrix = exprMatrix,
                       geneList = expr_HCCDB1$genelist,
                       sampleNum = 10,
                       sampleSize = 100,
                       P_value_thre = 0.05,
                       Corr_dis_thre = 0.25,
                       disMethod = 'spearman',
                       pairNum = 50000 ,
                       orderScore = 'dis',
                       savePath = outdir)


indir <- 'D:/Gu_lab/HCCsc/HCCsc_new/data/Hcc_bulk_DEGene/Multi_data/HCCDB-3/'
outdir <- 'D:/Gu_lab/HCCsc/HCCsc_new/result/scCorrNet/HCCDB3/'
expr = BulkDataInit(paste0(indir,'GSE25097.gene.txt'), paste0(indir,'GSE25097.sample.txt'), CellType = 'HCC', Genelistfile = 'none', IfSaveFile = T, SaveFileDir = paste0(outdir,'expr.csv'))

#test small data
small_exprMatrix<- expr$exprMatrix[1:1000,]
results <-runScCorrNet(exprMatrix = small_exprMatrix,
                       geneList = expr$genelist,
                       sampleNum = 10,
                       sampleSize = 100,
                       P_value_thre = 0.05,
                       Corr_dis_thre = 0.25,
                       disMethod = 'spearman',
                       pairNum = 50000 ,
                       orderScore = 'dis',
                       savePath = outdir)

