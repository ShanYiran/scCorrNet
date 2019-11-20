options(stringsAsFactors = F )
rm(list = ls())
library(Seurat)
library(dplyr)
library(rlist)
library(ggplot2)
library(Hmisc)
library(clusterProfiler)
library(scCorrNet)
indir <- 'D:/Gu_lab/HCCsc/HCCsc_new/data/Hcc_bulk_DEGene/Multi_data/HCCDB-1/'
outdir <- 'D:/Gu_lab/HCCsc/HCCsc_new/result/scCorrNet/HCCDB1/'
expr = BulkDataInit(paste0(indir,'GSE22058-GPL6793.gene.txt'), paste0(indir,'GSE22058.sample.txt'), CellType = 'HCC', Genelistfile = 'none', IfSaveFile = T, SaveFileDir = paste0(outdir,'expr.csv'))
results <- ConditionCorrNet(expr$exprMatrix,expr$genelist,sampleNum = 10, sampleSize = 100,orderScore = 'dis',IfSaveFile = T, saveFile = paste0(outdir,'results.txt'))
