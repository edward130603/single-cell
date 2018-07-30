#MAST Transcript set enrichment analysis

##load libraries
library("MultiAssayExperiment")
library("SummarizedExperiment")
library(MAST)
library(parallel)
options(mc.cores = detectCores() - 1)
library(data.table)

##load data
data1 = readRDS("data/GSE63818-GPL16791.rds")
#download from "http://imlspenticton.uzh.ch/robinson_lab/conquer/data-mae/GSE63818-GPL16791.rds"

##extract 
data1_tx = experiments(data1)[["tx"]] #transcripts
cData = colData(data1) #column info
fData = rowData(data1@ExperimentList$tx) #transcript info
mData = metadata(data1) #metadata

##prepare data
scaRaw = FromMatrix(log2(assays(data1_tx)[["TPM"]]+1), cData, fData)

filterCrit = mData$salmon_summary$percent_mapped > 60
sca = subset(scaRaw,filterCrit) #quality
sca = sca[which(freq(sca)> 0.2), ]

cdr2 = colSums(assay(sca)>0)
colData(sca)$cngeneson = scale(cdr2)

sca = subset(sca, sca$characteristics_ch1 == "developmental stage: 7 week gestation")
sca_name = paste0("tmp/sca_", Sys.Date(), ".RData")
save(sca, file = sca_name)
colData(sca)$source_name_ch1 = sample(colData(sca)$source_name_ch1) 

##run model
zlmCond = zlm(~source_name_ch1 + cngeneson, sca)
zlm_name = paste0("tmp/zlm_", Sys.Date(), ".RData")
save(zlmCond, file = zlm_name)
load(zlm_name)

##bootstrap
boots = bootVcov1(zlmCond, R = 50)
boots_name = paste0("tmp/boots_", Sys.Date(), ".RData")
save(boots, file = boots_name)
load(boots_name)

##create gene sets
genes = unique(mcols(sca)$gene)
sets = vector("list", length(genes))
names(sets) = genes

for (i in genes){
  sets[[i]] = which(mcols(sca)$gene == i)
}
sets_name = paste0("tmp/sets_", Sys.Date(), ".RData")
save(sets, file = sets_name)
load(sets_name)

##gsea
gsea = gseaAfterBoot(zlmCond, boots, sets, CoefficientHypothesis('source_name_ch1Somatic Cells'))
gsea_name = paste0("tmp/gsea_", Sys.Date(), ".RData")
save(gsea, file = gsea_name)
load(gsea_name)

z_stat_comb <- summary(gsea) #all genes
sigModules <- z_stat_comb[combined_adj<.05] #significant genes
hist(z_stat_comb$combined_P, main = "Distribution of MAST-TSEA p-values",
     xlab = "p")

